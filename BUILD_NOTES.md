# BUILD NOTES — pumpfun-top5-forecast-01 entry

## What it is

`main.py` is a **lifecycle (survival + graduation) forecaster**, not a GBM /
jump-diffusion Monte Carlo on log-cap (which is what most of the 1,440
entries do, judging by their titles). Each coin in the history is simulated
as a life:

1. **Curve phase** — mostly flat at the floor (log-cap ≈ 1.45); it can die
   (hazard 0.10/h) or graduate (hazard by log-cap: 1.2/h above 2.2, 0.4/h
   at 2.0–2.2, 0.08/h at 1.8–2.0, 0.01/h below; 2.0/h when flagged
   `complete`).
2. **Graduation** — a regime change with a repricing jump (mean −0.30,
   sd 0.60 in log10; practice graduations moved −0.86…+0.54).
3. **Pool phase** — fast growth while young (drift ≈ +0.45 log10/h,
   decaying with age, tau 1.5 h), then stagnation; death hazard rises
   with age but saturates (0.04 + 0.03·min(age−1, 4) per hour).
4. **Volatility decays with age** (tau 3 h) for both phases. This matters:
   with constant vol the lognormal tail explodes over a 7–13 h horizon and
   the simulated place-1 cap ran away to ≥10^6 SOL (bin 81 held p≈0.29 in
   an early version). With decaying vol the place-1 size peaks at bins
   46–50 (≈20k–32k SOL) at long horizons.
5. **Newborn field** — coins launched after the history ends, which no
   forecast can name. Per read, K ~ Poisson(λ(h)) coins with log-caps ~
   Normal(μ(h), 0.45), λ(h) = 2.0 + 0.2·h, μ(h) = (last mark's top-3 median
   log-cap) + 0.05·h. They displace known coins, so the place
   probabilities coming out of the simulation are honest about how much
   of each place is unnameable — that mass is left unnamed.

Two more edges:

- **p=0 rows.** Coins very likely still on the answer chart but very
  unlikely to take a top-5 place (flat low curve coins ranked 25–90, and
  low pool coins outside the candidate set) are named at p=0 at each
  place. That takes them out of the place's unnamed share and concentrates
  it on the coins that could actually take the place. Never a fast grower
  or a near-graduation coin.
- **Sizes from the simulation.** The bin distribution at each place is the
  simulated distribution of that place's cap, pooled over paths and reads:
  wide, and shifted up with the horizon — not a narrow normal around the
  last mark's cap like the house baseline (sd 0.35).

## Key parameters (all in `main.py`, top of file)

- Seed `987654321`, 2,500 paths, 7 reads from `task.json` (never hardcoded).
- Candidates: top 100 by last-mark cap + every pool coin + every fast
  grower (slope > 0.2/h), capped at 150.
- Pool drift 0.45/h decaying (tau 1.5 h); pool vol 0.45, curve vol 0.25,
  both decaying (tau 3 h). Death: curve 0.10/h; pool 0.04 + 0.03·age.
- Newborns: λ = 2.0 + 0.2·h, μ = top-3 median + 0.05·h, σ = 0.45.
- Naming thresholds: p ≥ 0.003 for probabilities, p ≥ 0.002 for size bins;
  probabilities floored to 6 decimals (keeps place sums ≤ 1).

## Tests (all on the kit's practice data)

- **Format + determinism:** forecast valid (exactly 3 keys, place sums
  ≤ 1, ≤100 rows/place, top5 = 5 distinct history mints, bins 0–81,
  44,366 bytes ≤ 1 MB). Two runs byte-identical. Runtime ≈ 24 s on the
  practice input (limit 300 s).
- **Backtest** (history = marks 1–3, i.e. 20:10–20:30; answer reads =
  marks 4–5, 20:40/20:50; scorer implements the rules' log score with
  the 0.0001 floor):

  | forecast | order | size | total | agreement |
  |---|---|---|---|---|
  | **mine (lifecycle + p=0)** | −14.95 | −10.59 | **−25.53** | 25 |
  | mine, no p=0 | −15.28 | −10.96 | −26.24 | 21 |
  | house baseline | −18.59 | −14.89 | −33.48 | 25 |
  | even spread | | | −45.06 | |

  The model beats the baseline by **+7.94** on the backtest, on both
  reads (read 1: −22.97 vs −29.65; read 2: −28.10 vs −37.31). The p=0
  rows add +0.71 and no p=0 coin took a place in the backtest.
- **Long-horizon sanity** (same practice history, reads faked to 7–13 h
  after the cut): named mass falls to 0.40–0.55 per place (unnamed
  0.45–0.60 — honest about newborns dominating at long horizons),
  place-1 size peaks at bins 46–50, bin 81 mass negligible. Before the
  decaying-vol fix, bin 81 held p≈0.29 and place-1 named mass was 0.98.

## Known weaknesses (read before trusting it)

- **Long-horizon calibration is untested.** Only the 10–20 min backtest
  is empirical; the 7–13 h behaviour relies on the structure (decaying
  vol/drift, saturating death hazard, growing newborn field), not on data.
  If real 7–13 h survival is higher than modelled, the forecast is too
  humble (too much unnamed); if lower, too confident on known coins.
- **Newborn field is calibrated from one data point** (2 newborns in the
  practice's last-read top 5, a third coin never in the top 200). λ and μ
  are guesses extrapolated with horizon, not fitted.
- **Hazards come from 50 minutes of practice data** (6 graduations, 19
  pool→pool transitions). The final's 11-hour history will give far better
  estimates, but the code does not re-estimate hazards from the input —
  the numbers are fixed constants. (The slopes and ages ARE read from the
  input; only the hazard/drift constants are fixed.)
- **Only coins in the history's rows can be named.** A coin below the
  top 200 that explodes — one did exactly that in the practice — cannot
  be named; it is covered only through the unnamed share.
- **p=0 rows are a calculated risk.** Safe in the backtest, but if a p=0
  coin takes a place at a read, that place scores the floor (−9.21)
  there. The set deliberately excludes fast growers and near-graduation
  coins to keep that risk small.
- **Backtest has only 2 reads.** The +7.94 gap is large and consistent
  across both reads, but the sample is small; do not over-read it.
- **Sizes are pooled over reads.** The same size distribution is scored
  at every read, so it must cover both the 7 h and the 13 h cap levels;
  it is wide by design, which costs a little at each individual read
  against a perfectly centred narrow distribution.
