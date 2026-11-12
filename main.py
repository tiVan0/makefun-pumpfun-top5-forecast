#!/usr/bin/env python3
"""Lifecycle forecaster for make.fun pumpfun-top5-forecast-01.

Approach: survival + graduation, NOT a GBM Monte Carlo on log-cap.

Every coin in the history is modelled as a life: it is born on its bonding
curve, it may graduate to its canonical PumpSwap pool (a regime change with a
repricing jump), it grows fast while young and then its growth decays with
age, and it may die (its cap collapses to the curve floor).  Coins are
simulated forward from the cut to each of the seven answer reads (read from
task.json, never hardcoded) with a fixed seed, so two runs are byte-identical.

On top of the known coins, each read also draws a field of newborn coins
(coins launched after the history ends, which no forecast can name): their
number is Poisson and their caps are lognormal around the current top of the
chart.  The newborns displace known coins, so the place probabilities that
come out of the simulation are honest about how much of each place is
unnameable -- that mass is left unnamed, which is what scores when a
newcomer takes a place.

Two more edges over a plain persistence forecast:
  * coins that are almost certainly dead for top-5 purposes but very likely
    still on the answer chart are named at p=0 at a place, which takes them
    out of that place's unnamed share and concentrates it on the coins that
    could actually take the place (newborns included);
  * size (bin) probabilities come from the simulated distribution of the
    cap at each place, pooled over paths and reads: wide, and shifted up
    with the horizon, instead of a narrow normal around the last mark's cap.

Standard library only. No network. No model calls. Deterministic.
Reads $MAKEFUN_INPUT_DIR (else /input), writes $MAKEFUN_OUTPUT_DIR (else
/output).
"""

import json
import math
import os
import random
from datetime import datetime, timezone

SEED = 987654321
N_PATHS = 2500

# log10 of a market cap in SOL at which a dead curve coin sits.
FLOOR_LOG = 1.45
# A coin's cap can never be simulated below this (a few SOL).
MIN_LOG = 1.0
# A coin's cap is capped here (10^6.5 SOL) to keep bins sane.
MAX_LOG = 6.5

# --- lifecycle parameters, estimated from the practice history itself ---
# Pool-coin growth: young pool coins grow ~+0.45 log10/hour, decaying with
# age (tau 1.5 h).  Volatility also decays with age (tau 3 h): a coin that
# is hours old no longer random-walks with a young coin's volatility, or
# the lognormal tail explodes over a 7-13 h horizon and the simulated caps
# run away to 10^6 SOL.  The integrated variance over a horizon is closed
# form for an exponentially decaying instantaneous vol.
POOL_DRIFT0 = 0.45
POOL_DRIFT_TAU = 1.5
POOL_VOL0 = 0.45
CURVE_VOL0 = 0.25
VOL_TAU = 3.0
# Death hazards per hour.  Curve coins die fast.  Pool coins die slowly at
# first and faster as they age, but the hazard saturates: a coin that has
# survived for hours is a strong coin, not a certain corpse.
CURVE_DEATH_HAZARD = 0.10
POOL_DEATH_BASE = 0.04
POOL_DEATH_AGE = 0.03
POOL_DEATH_AGE_CAP = 4.0
# Graduation hazard per hour for a curve coin, by log-cap.  Graduation is
# concentrated in the highest curve coins (the practice saw 6 graduations,
# all from coins at log-cap 2.2-2.5); a coin flagged complete graduates
# imminently.
GRAD_HAZARD = ((2.2, 1.2), (2.0, 0.4), (1.8, 0.08), (float("-inf"), 0.01))
GRAD_HAZARD_COMPLETE = 2.0
# On graduation the cap reprices: practice graduations moved log-cap by
# -0.86..+0.54, mean about -0.30.
GRAD_REPRICE_MEAN = -0.30
GRAD_REPRICE_SD = 0.60
# Newborn field: coins launched after the history ends, which no forecast
# can name.  Their number per read grows with the horizon (the post-cut
# cohort accumulates) and so does the cap a newborn needs to crack the
# top 5.  In the practice, 2 of the top 5 at the last read were launched
# after the history ended and a third was never in the top 200.
NEWBORN_LAMBDA0 = 2.0
NEWBORN_LAMBDA_PER_HOUR = 0.2
NEWBORN_MU_PER_HOUR = 0.05
NEWBORN_SD = 0.45
# Only coins with at least this simulated probability are named with p > 0.
MIN_P = 0.003
# Size bins are named when their simulated probability is at least this.
MIN_SIZE_P = 0.002


def parse_time(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def logcap_sol(lamports):
    return math.log10(int(lamports)) - 9.0


def floor6(x):
    """Round a probability down to 6 decimals (keeps place sums <= 1)."""
    return math.floor(x * 1000000.0 + 1e-9) / 1000000.0


def grad_hazard(logcap, complete):
    if complete:
        return GRAD_HAZARD_COMPLETE
    for threshold, hazard in GRAD_HAZARD:
        if logcap >= threshold:
            return hazard
    return 0.01


def pool_drift(age_h, horizon_h, slope):
    """Integrated log-cap drift of a pool coin over horizon_h hours.

    Drift starts at slope (shrunk toward the young-pool prior) and decays
    exponentially with age; the integral over the horizon is closed form.
    """
    d0 = 0.5 * slope + 0.5 * (POOL_DRIFT0 * math.exp(-age_h / 1.0) + 0.05)
    return d0 * POOL_DRIFT_TAU * (1.0 - math.exp(-horizon_h / POOL_DRIFT_TAU))


def vol_sd(v0, age_h, horizon_h):
    """Std dev of a coin's log-cap move over horizon_h hours.

    The instantaneous volatility decays exponentially with age (tau 3 h),
    so the integrated variance saturates instead of growing like sqrt(h)
    forever: an old coin barely moves, a young one swings.
    """
    a0 = math.exp(-2.0 * age_h / VOL_TAU)
    a1 = math.exp(-2.0 * (age_h + horizon_h) / VOL_TAU)
    var = v0 * v0 * (VOL_TAU / 2.0) * (a0 - a1)
    return math.sqrt(max(0.0, var))


def pool_death_prob(age_h, horizon_h):
    """P(a pool coin dies within horizon_h hours).

    The hazard rises with age but saturates (a coin that survived for
    hours is a strong coin); the average hazard over the horizon is the
    trapezoid of the two ends.
    """
    def lam(a):
        return POOL_DEATH_BASE + POOL_DEATH_AGE * min(
            max(0.0, a - 1.0), POOL_DEATH_AGE_CAP)
    avg = 0.5 * (lam(age_h) + lam(age_h + horizon_h))
    return 1.0 - math.exp(-avg * horizon_h)


def exp_time(rng, hazard):
    """Draw from Exp(hazard) using a uniform; returns inf for hazard 0."""
    if hazard <= 0.0:
        return float("inf")
    u = rng.random()
    if u <= 0.0:
        u = 1e-12
    return -math.log(1.0 - u) / hazard


def poisson(rng, lam):
    """Knuth's algorithm; lam is small (a few)."""
    limit = math.exp(-lam)
    k = 0
    p = 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


class Coin(object):
    __slots__ = ("mint", "logcap", "venue", "age_h", "slope", "complete",
                 "last_rank")

    def __init__(self, mint, logcap, venue, age_h, slope, complete,
                 last_rank):
        self.mint = mint
        self.logcap = logcap
        self.venue = venue
        self.age_h = age_h
        self.slope = slope
        self.complete = complete
        self.last_rank = last_rank


def build_coins(marks, cut):
    """Build one Coin per mint seen in the history, with a shrunk slope."""
    traj = {}
    for m in marks:
        t = parse_time(m["at"])
        for r in m["rows"]:
            mint = r["mint"]
            entry = traj.get(mint)
            if entry is None:
                entry = {"obs": [], "created": r["created_at"],
                         "last": None}
                traj[mint] = entry
            entry["obs"].append((t, logcap_sol(r["market_cap_lamports"]),
                                 r["venue"]))
            entry["last"] = r
    coins = []
    for mint, entry in traj.items():
        obs = entry["obs"]
        last = entry["last"]
        created = parse_time(entry["created"])
        age_h = max(0.0, (cut - created).total_seconds() / 3600.0)
        # Slope from the last up to 4 observations (per hour), shrunk
        # toward 0: short trajectories carry little signal.
        slope = 0.0
        pts = obs[-4:]
        if len(pts) >= 2:
            t0 = pts[0][0]
            xs = [(p[0] - t0).total_seconds() / 3600.0 for p in pts]
            ys = [p[1] for p in pts]
            n = len(pts)
            mx = sum(xs) / n
            my = sum(ys) / n
            den = sum((x - mx) ** 2 for x in xs)
            if den > 1e-9:
                raw = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
                # Clip wild slopes and shrink toward 0.
                raw = max(-0.6, min(0.6, raw))
                slope = raw * (0.5 if n >= 3 else 0.3)
        venue = last["venue"]
        if venue == "migrating":
            venue = "pool"
        coins.append(Coin(mint, logcap_sol(last["market_cap_lamports"]),
                          venue, age_h, slope, bool(last.get("complete")),
                          int(last["rank"])))
    coins.sort(key=lambda c: (-c.logcap, c.mint))
    return coins


def simulate_coin(rng, coin, horizon_h):
    """Simulate one coin's log-cap horizon_h hours after the cut."""
    if coin.venue == "pool":
        if rng.random() < pool_death_prob(coin.age_h, horizon_h):
            return FLOOR_LOG + rng.gauss(0.0, 0.05)
        drift = pool_drift(coin.age_h, horizon_h, coin.slope)
        cap = coin.logcap + drift \
            + vol_sd(POOL_VOL0, coin.age_h, horizon_h) * rng.gauss(0.0, 1.0)
        return max(MIN_LOG, min(MAX_LOG, cap))
    # Curve coin: race between death and graduation.
    lam_d = CURVE_DEATH_HAZARD
    lam_g = grad_hazard(coin.logcap, coin.complete)
    t_die = exp_time(rng, lam_d)
    t_grad = exp_time(rng, lam_g)
    if t_die < t_grad and t_die < horizon_h:
        return FLOOR_LOG + rng.gauss(0.0, 0.05)
    if t_grad < horizon_h:
        pre = vol_sd(CURVE_VOL0, coin.age_h, t_grad) * rng.gauss(0.0, 1.0)
        reprice = GRAD_REPRICE_MEAN + GRAD_REPRICE_SD * rng.gauss(0.0, 1.0)
        rest = horizon_h - t_grad
        age_at_grad = coin.age_h + t_grad
        if rng.random() < pool_death_prob(age_at_grad, rest):
            return FLOOR_LOG + rng.gauss(0.0, 0.05)
        drift = pool_drift(age_at_grad, rest, 0.0)
        cap = coin.logcap + pre + reprice + drift \
            + vol_sd(POOL_VOL0, age_at_grad, rest) * rng.gauss(0.0, 1.0)
        return max(MIN_LOG, min(MAX_LOG, cap))
    cap = coin.logcap \
        + vol_sd(CURVE_VOL0, coin.age_h, horizon_h) * rng.gauss(0.0, 1.0)
    return max(MIN_LOG, min(MAX_LOG, cap))


def bin_of(logcap, edges):
    """Size bin 0..81 for a log10 cap in SOL (edges are 2.00..6.00)."""
    if logcap < edges[0]:
        return 0
    if logcap >= edges[-1]:
        return len(edges)
    # Binary search over the 81 edges.
    lo, hi = 0, len(edges)
    while lo < hi:
        mid = (lo + hi) // 2
        if edges[mid] <= logcap:
            lo = mid + 1
        else:
            hi = mid
    return lo


def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")
    with open(os.path.join(in_dir, "task.json"), "r", encoding="utf-8") as f:
        task = json.load(f)
    places = int(task["places"])
    edges = [float(e) for e in task["size"]["edges"]]
    reads = [parse_time(s) for s in task["reads"]]
    marks = []
    with open(os.path.join(in_dir, "history.ndjson"), "r",
              encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = json.loads(line)
            if m.get("ok") and isinstance(m.get("rows"), list) and m["rows"]:
                marks.append(m)
    if not marks:
        raise SystemExit("the history holds no ok mark with rows")
    marks.sort(key=lambda m: m["at"])
    cut = parse_time(marks[-1]["at"])
    horizons = []
    for r in reads:
        h = (r - cut).total_seconds() / 3600.0
        if h < 0.0:
            h = 0.0
        horizons.append(h)

    coins = build_coins(marks, cut)
    last_chart = marks[-1]["rows"]
    last_chart.sort(key=lambda r: r["rank"])
    last_caps = [logcap_sol(r["market_cap_lamports"])
                 for r in last_chart[:3]]
    newborn_mu = sorted(last_caps)[len(last_caps) // 2]

    # Candidates: the top 100 by last cap, every pool/migrating coin, and
    # every fast grower; capped so the simulation stays fast.
    cand = {}
    for c in coins[:100]:
        cand[c.mint] = c
    for c in coins:
        if c.venue == "pool":
            cand[c.mint] = c
        elif c.slope > 0.2:
            cand[c.mint] = c
    cand_list = sorted(cand.values(), key=lambda c: (-c.logcap, c.mint))
    cand_list = cand_list[:150]
    cand_index = {c.mint: i for i, c in enumerate(cand_list)}

    rng = random.Random(SEED)
    n_reads = len(reads)
    total = N_PATHS * n_reads
    # place_counts[place][coin_idx] and size_counts[place][bin].
    place_counts = [[0] * len(cand_list) for _ in range(places)]
    size_counts = [[0] * (len(edges) + 1) for _ in range(places)]

    for _ in range(N_PATHS):
        for h in horizons:
            # Simulate the known coins.
            sims = []
            for i, c in enumerate(cand_list):
                sims.append((simulate_coin(rng, c, h), i))
            # Newborn field: coins launched after the cut, which no
            # forecast can name.  The field grows with the horizon (the
            # post-cut cohort accumulates) and so does the cap a newborn
            # needs to crack the top 5.
            lam_new = NEWBORN_LAMBDA0 + NEWBORN_LAMBDA_PER_HOUR * h
            mu_new = newborn_mu + NEWBORN_MU_PER_HOUR * h
            k_new = poisson(rng, lam_new)
            for _ in range(k_new):
                sims.append((mu_new + NEWBORN_SD * rng.gauss(0.0, 1.0), -1))
            sims.sort(key=lambda t: (-t[0], t[1]))
            for k in range(places):
                cap, idx = sims[k]
                size_counts[k][bin_of(cap, edges)] += 1
                if idx >= 0:
                    place_counts[k][idx] += 1

    # Place probabilities: average over paths and reads.
    probs = []  # (mint, place, p)
    named = [set() for _ in range(places)]
    for k in range(places):
        for i, c in enumerate(cand_list):
            p = place_counts[k][i] / total
            if p >= MIN_P:
                probs.append((c.mint, k + 1, floor6(p)))
                named[k].add(c.mint)
    # p=0 rows: coins very likely on the answer chart but very unlikely to
    # take a top-5 place.  They leave the place's unnamed share to the
    # coins that could actually take it.  Pool coins persist on the chart;
    # flat, low curve coins sit at the floor and cannot reach the top 5.
    # Never name a fast grower or a near-graduation coin at p=0.
    # (MAKEFUN_NO_P0=1 disables this for ablation testing.)
    if os.environ.get("MAKEFUN_NO_P0") != "1":
        p0_pool = [c for c in coins if c.venue == "pool"
                   and c.mint not in cand_index]
        p0_pool.sort(key=lambda c: (c.last_rank, c.mint))
        p0_curve = [c for c in coins
                    if c.venue == "curve" and c.slope <= 0.0
                    and c.logcap < 2.0 and 25 <= c.last_rank <= 90]
        p0_curve.sort(key=lambda c: (c.last_rank, c.mint))
        p0 = p0_pool + p0_curve
        for k in range(places):
            rows_here = sum(1 for (m, pl, p) in probs if pl == k + 1)
            for c in p0:
                if rows_here >= 100:
                    break
                if c.mint in named[k]:
                    continue
                probs.append((c.mint, k + 1, 0.0))
                named[k].add(c.mint)
                rows_here += 1
    probs.sort(key=lambda t: (t[1], -t[2], t[0]))
    probabilities = [{"mint": m, "place": pl, "p": p}
                     for (m, pl, p) in probs]

    # top5: rank by expected place score.
    score = {}
    for k in range(places):
        for i, c in enumerate(cand_list):
            p = place_counts[k][i] / total
            if p > 0.0:
                score[c.mint] = score.get(c.mint, 0.0) + p * (places + 1 - (k + 1))
    ranked = sorted(
        score.items(),
        key=lambda kv: (-kv[1],
                        -place_counts[0][cand_index[kv[0]]] / total,
                        kv[0]))
    top5 = [m for m, _ in ranked[:places]]
    # Fall back to the last mark's top 5 if the simulation named nothing.
    if len(top5) < places:
        seen = set(top5)
        for r in last_chart:
            if len(top5) >= places:
                break
            if r["mint"] not in seen:
                top5.append(r["mint"])
                seen.add(r["mint"])

    # Sizes: pooled simulated distribution of the cap at each place.
    sizes = []
    for k in range(places):
        tot = sum(size_counts[k])
        for b in range(len(edges) + 1):
            p = size_counts[k][b] / tot if tot else 0.0
            if p >= MIN_SIZE_P:
                sizes.append({"place": k + 1, "bin": b, "p": floor6(p)})
    sizes.sort(key=lambda s: (s["place"], -s["p"], s["bin"]))

    forecast = {"probabilities": probabilities, "top5": top5, "sizes": sizes}
    out_path = os.path.join(out_dir, task["output"])
    os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(json.dumps(forecast) + "\n")


if __name__ == "__main__":
    main()
