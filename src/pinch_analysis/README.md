# Pinch with lowered T_HENE_in

PFHX-5 pinch + area pipeline at the baseline target $\Delta T_\mathrm{min} =
2.48$ K (rounded to 2.5 K in chapter prose), with the He-Ne inlet temperature
lowered slightly for each ejector case to make the baseline target feasible.

## What changed vs the v2 pipeline

Only the He-Ne inlet temperature into PFHX-5, set per case:

| case      | $f$  | $\eta$ | $T_\mathrm{HeNe,in}$ [K] | $\Delta T_\mathrm{min}$ target |
|---        |---   |---     |---                       |---                              |
| baseline  | --   | --     | 27.45 (Moro reference)   | 2.48 K (natural, no solver)    |
| low_eta   | 0.20 | 0.20   | 26.20                    | 2.48 K (target met)            |
| high_eta  | 0.10 | 0.30   | 26.70                    | 2.48 K (target met)            |

The $T_\mathrm{HeNe,in}$ values were chosen from a separate sweep that
mapped $\Delta T_\mathrm{min}$ vs $T_\mathrm{HeNe,in}$ at each ejector
operating point. Each value sits just below the threshold where the
2.48 K target becomes feasible.

The dT_min fallback ladder (1.90 K primary, 1.0 K fallback) used in the
previous pipeline is removed because the lowered $T_\mathrm{HeNe,in}$
makes the single 2.48 K target feasible at both ejector cases.

The integral-method area calculation is computed internally as a sanity
check but is no longer printed in the summary; the side-by-side table
reports only single LMTD and piecewise LMTD ($N=100$).

## Important caveat

$T_\mathrm{HeNe,in}$ is not a free design variable in the Moro
architecture. It is set by the He-Ne Brayton expander outlet
temperature, which depends on the Brayton expander pressure ratio and
the upstream cooling. Lowering $T_\mathrm{HeNe,in}$ in practice requires
either a higher Brayton expander pressure ratio or additional upstream
cooling, both of which increase Brayton compressor work. The PFHX-5-side
improvement reported here must be traded against this Brayton penalty
before any conclusion that "lowering $T_\mathrm{HeNe,in}$ helps" the
overall plant SEC. That trade-off is not modelled here.

## Directory layout

```
Pinch_with_lowered_THENEin/
  ejector/
    ejector_solve.m              unchanged from v13 pipeline
    run_pinch_inputs.m           patched: writes T_HENE_in_K column per case
    refpropm.m, rp_proto*.m
  pinch/
    pinch_PFHX5.m                patched: reads T_HENE_in_K per case,
                                  no dT_min fallback (single 2.48 K target)
    area_two_cases.m             patched: integral method dropped from
                                  printed summary; constant-dT_min note
    plot_pinch_PFHX5_final.py    patched: dT_min labels rounded to 1 decimal
                                  (so 2.48 K renders as 2.5 K)
    PLOT_SETTINGS.py             unchanged
    refpropm.m, rp_proto*.m
  README.md
```

## Workflow

```matlab
cd ejector
run_pinch_inputs                 % writes pinch_inputs.csv (2 rows + new col)
```

```matlab
cd ../pinch
pinch_PFHX5                      % auto-finds inputs, writes 3 pinch_curves CSVs
area_two_cases                   % UA + length per case (single LMTD + pw100)
```

```bash
python3 plot_pinch_PFHX5_final.py    % 3-panel composite curves figure
```

## Expected output

Console summary from `pinch_PFHX5.m` should report dT_min ~ 2.48 K for
both ejector cases, mdot_HENE/baseline ~ 1.27 (low_eta) and ~ 1.07
(high_eta), with the pinch in the interior of PFHX-5 (Q_pinch in the
700-1200 kW range, well above the 25 kW cold-end pinch seen at
$T_\mathrm{HeNe,in} = 27.45$ K).

Console summary from `area_two_cases.m` should report L_pw100 values
that are now meaningfully closer to the baseline 6.0 m than the
$T_\mathrm{HeNe,in} = 27.45$ K case produced (which gave 18.6 m and
11.8 m at 1.0 K and 1.90 K respectively). With both cases at
constant 2.48 K and only $\sim$1 K lower $T_\mathrm{HeNe,in}$, the
area ratios should drop substantially.

The 3-panel pinch figure will show all three cases at the same
operating margin ($\Delta T_\mathrm{min} = 2.5$ K rounded), making the
curve-geometry differences (Widom curvature) immediately visible
without dT_min confounding.

## What was unchanged

- ejector_solve.m -- the ejector solution does not depend on $T_\mathrm{HeNe,in}$
- run_baseline -- baseline keeps $T_\mathrm{HeNe,in} = 27.45$ K (the Moro reference)
- All helper functions (build_HENE, pinch_residual, pinch_metrics, ...)
- PLOT_SETTINGS.py
