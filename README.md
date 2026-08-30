# RobustGaSP-MATLAB

A from-scratch MATLAB implementation of

> M. Gu, X. Wang and J. O. Berger (2018).
> **Robust Gaussian Stochastic Process Emulation.**
> *The Annals of Statistics* 46(6A), 3038–3066. ([arXiv:1708.04738](https://arxiv.org/abs/1708.04738))

with the API and default settings of the R package
[**RobustGaSP**](https://cran.r-project.org/package=RobustGaSP) (Gu, Palomo and Berger, *The R Journal*, 2019).

The **Optimization Toolbox** drives the posterior mode search (`fmincon`, or `fminunc`
when unconstrained) and the **Statistics and Machine Learning Toolbox** supplies the
quantiles, designs and random draws. An **optional** C++ (MEX) layer, compiled once,
removes the two dominant costs of a fit. Where a toolbox is absent the code still runs
on built-in fallbacks — `rgasp_toolbox_report` says which path is live, and the test
suite runs itself both ways.

---

## 1. What is implemented

| Capability | Function |
|---|---|
| Scalar-output emulator (marginal posterior mode) | `rgasp` |
| Vector/field-output emulator (PP GaSP) | `ppgasp` |
| Predictive mean, sd, 95% interval | `rgasp_predict`, `ppgasp_predict` |
| Joint sample paths from the predictive process | `rgasp_simulate` |
| Leave-one-out cross validation | `rgasp_loo` |
| Inert-input detection | `rgasp_inert_inputs` |
| RMSE / coverage / interval-length criteria | `rgasp_validate` |
| Printed model summary | `rgasp_summary` |
| Latin hypercube designs | `rgasp_lhs` (`lhsdesign`), `rgasp_lhs_det` |
| Path setup | `setup_robustgasp` |
| Which optional components are live | `rgasp_toolbox_report`, `rgasp_use_toolboxes` |
| Optional C++ acceleration | `rgasp_compile_mex`, `rgasp_have_mex` |

Estimation methods: `post_mode` (default — marginal posterior mode with the jointly
robust prior), `mmle` (marginal likelihood, no prior) and `mle` (profile likelihood).
Correlation functions: `matern_5_2` (default), `matern_3_2`, `pow_exp`, isotropic or
anisotropic, with or without a nugget, with a zero mean or an arbitrary trend basis,
and mixed kernels across input dimensions.

## 2. Quick start

```matlab
run('/path/to/RobustGaSP-MATLAB/setup_robustgasp.m')

x = linspace(0, 10, 15)';
y = higdon_1_data(x);

model = rgasp(x, y);            % marginal posterior mode, Matern 5/2
rgasp_summary(model);

xt   = linspace(0, 10, 400)';
pred = rgasp_predict(model, xt);

plot(xt, pred.mean, 'b-', xt, pred.lower95, 'b:', xt, pred.upper95, 'b:', x, y, 'ro');
rgasp_validate(pred, higdon_1_data(xt));
```

Vector output (a field per run):

```matlab
Y  = environ_4_data(D);         % 40-by-200
mp = ppgasp(D, Y);
pp = ppgasp_predict(mp, Dtest); % nt-by-200 mean / sd / lower95 / upper95
```

## 3. The model and the estimator

For design `x^D = (x_1,…,x_n)` and outputs `y^D`,

```
y(x) = h(x) θ + z(x),      z ~ GaSP(0, σ² c(·,·)),      c(x_i,x_j) = ∏_l c_l(|x_il − x_jl|)
```

With the objective prior `π(θ, σ²) ∝ 1/σ²` both are integrated out analytically, giving
the marginal likelihood used for estimation (eq. 2.5–2.7 of the paper):

```
L(γ | y^D) ∝ |R|^(−1/2) · |Hᵀ R⁻¹ H|^(−1/2) · (S²)^(−(n−q)/2)
S²  = y^Dᵀ Q y^D,   Q = R⁻¹ − R⁻¹H (Hᵀ R⁻¹ H)⁻¹ Hᵀ R⁻¹
```

`R = R(β) + η I` with the **inverse range** parameterization `β_l = 1/γ_l`, and the
optimization variable is `ξ_l = log β_l` — the parameterization the paper shows to give
robust estimates.

The default estimator is the mode of the marginal posterior formed with the
**jointly robust prior** (`prior_choice = "ref_approx"` in the R package):

```
π_JR(β, η) ∝ ( Σ_l C_l β_l + η )^a · exp{ −b ( Σ_l C_l β_l + η ) }
a = 0.2,   b = (a + p)/n^(1/p),   C_l = (max_i x_il − min_i x_il)/n^(1/p)
```

Prediction is a Student-t (eq. 2.10–2.13):

```
y(x*) | y^D, γ̂  ~  T( ŷ(x*), σ̂² c**(x*), n − q )
ŷ(x*)  = h(x*) θ̂ + r(x*)ᵀ R⁻¹ (y^D − H θ̂)
c**    = 1 + η − rᵀR⁻¹r + (h(x*) − HᵀR⁻¹r)ᵀ (HᵀR⁻¹H)⁻¹ (h(x*) − HᵀR⁻¹r)
```

**PP GaSP.** For an `n × k` output matrix the correlation over the input space is shared
by all `k` coordinates while each coordinate keeps its own mean vector and variance:

```
L ∝ |R|^(−k/2) |HᵀR⁻¹H|^(−k/2) ∏_{i=1..k} (S_i²)^(−(n−q)/2)
```

so one fit serves an entire output field, and `sigma2_hat` is a `k`-vector.

### Gradients

Analytic gradients with respect to `ξ` are supplied to the optimizer. With
`Ṙ_l = ∂R/∂β_l` and `W_l = Ṙ_l Q`,

```
∂ℓ/∂β_l = −(k/2)·tr(W_l) + ((n−q)/2)·Σ_i (Qy_i)ᵀ Ṙ_l (Qy_i) / S_i²
```

and `∂/∂ξ_l = β_l · ∂/∂β_l`. Because `R` is a Hadamard product, `Ṙ_l = R ∘ (d log c_l/dβ_l)`,
which is evaluated from closed forms that never divide by a (possibly underflowed)
correlation entry — e.g. for Matérn 5/2, `d log c/dβ = −√5 d · u(1+u)/(3+3u+u²)`, `u = √5βd`.

All gradients are verified against central differences in the test suite
(108 configurations, max relative error `4.4e-7`).

## 4. Numerical safeguards

* **Lower bound on `log β`** (`lowerBound`, default `true`). The likelihood is flat as
  `β → 0`, where `R` degenerates to a matrix of ones. `rgasp_search_lb` finds the common
  correlation floor `q` such that `β_l = −log(q)/max_ij|x_il − x_jl|` gives
  `cond(R) = 10^16`, by a golden-section search on the logit of `q` — the same
  construction as `search_LB_prob` in the R package.
* **Multi-start optimization** (`numInitialValues`, default 2), with the same two default
  starting points as the R package (`50·exp(LB)` and `(a+p)/(p·C_l·b)/2`).
* **Conditioning guard.** Parameter values that make `R` numerically singular return a
  large penalty so the line search backtracks instead of inverting a singular matrix.
* **Optimizer.** `optimizer` defaults to `'auto'`, which selects `fmincon`
  (interior-point, analytic gradient, limited-memory BFGS Hessian) from the
  Optimization Toolbox, or `fminunc` when `lowerBound` is `false`. `'lbfgs'` selects the
  built-in projected L-BFGS — two-loop recursion, active-set projection, Armijo
  backtracking along the projected arc — which is also the fallback where the toolbox is
  absent; `'neldermead'` is a derivative-free option. On the borehole function the two
  gradient-based solvers agree on the mode to the digit printed:

  | n | `fmincon` log posterior | `lbfgs` log posterior | `fmincon` nRMSE | `lbfgs` nRMSE |
  |---|---|---|---|---|
  | 100 | −238.3668 | −238.3686 | 0.00723 | 0.00722 |
  | 250 | −383.8629 | −383.8638 | 0.00217 | 0.00217 |
  | 500 | −429.6746 | −429.6754 | 0.00054 | 0.00054 |
  | 1000 | 18.7473 | 18.7294 | 0.00023 | 0.00023 |

## 5. Toolbox use

| Toolbox function | Used for | Fallback when absent |
|---|---|---|
| `fmincon` | default optimizer for the posterior mode | `rgasp_lbfgs` |
| `fminunc` | unconstrained fits (`lowerBound = false`) | `rgasp_lbfgs` |
| `optimoptions` | solver option objects | `optimset` |
| `tinv` | Student-t predictive quantiles | bracketed inversion of the t CDF via `betainc` |
| `norminv` | normal quantiles (`method = 'mle'`) | `erfinv` |
| `lhsdesign` | Latin hypercube designs (`rgasp_lhs`) | `rgasp_lhs_det` |
| `mvnrnd` | joint Gaussian draws in `rgasp_simulate` | Cholesky + `randn` |
| `chi2rnd` | the t scaling in `rgasp_simulate` | sum of squared `randn` |
| `pdist2` | isotropic distance matrices | explicit `bsxfun` loop |
| `quantile`, `boxplot` | demo summaries and figures | linear interpolation / hand-drawn box plot |

`rgasp_toolbox_report()` prints this table with the live status of each entry.

Two notes on why the fallbacks are still there rather than deleted:

* `rgasp_use_toolboxes(false)` forces every fallback, so you can check that a result does
  not depend on a toolbox being installed. `run_all_tests('fallback')` runs the whole
  suite that way, and `test_toolbox_paths` compares the two paths directly.
* Writing the fallback found a real bug worth recording. The obvious one-line t-quantile,
  `betaincinv` on the Beta(ν/2, ½) representation, is **not** reliable in the far tail:
  for ν = 38 and p = 0.001 it returns a root whose actual tail probability is 0.015, not
  0.001. That would have corrupted any interval requested at a level far from 95%. The
  fallback now brackets and bisects the CDF and is checked against it in the tests; the
  toolbox `tinv` was correct throughout.

The test suite calls `rgasp_lhs_det` rather than `rgasp_lhs`, so recorded test and demo
results are byte-identical across installations; `rgasp_lhs` prefers `lhsdesign`.

## 6. Optional C++ acceleration

```matlab
rgasp_compile_mex();     % once; needs a compiler configured for MEX
```

Nothing requires this. Every MEX function has a pure-MATLAB fallback, `useMex` is a
per-call option, and `test_mex_equivalence` checks the two paths against each other.
Two functions are compiled:

* `rgasp_corr_mex` — the separable correlation matrix, evaluated **directly from the
  design matrices**. The product over the `p` dimensions is formed element by element,
  so the `n`-by-`n` result is written once instead of `p` times and the `p`
  per-dimension distance matrices are never allocated at all.
* `rgasp_gradterms_mex` — both per-dimension gradient terms in one fused pass, with
  `Rdot_l = R ∘ (d log c_l/dβ_l)` never materialized and the loop over `l` innermost, so
  the `n`-by-`n` arrays are streamed once in total rather than once per dimension.

The MATLAB side was reorganized to match: the trend correction enters only through the
low-rank factor `B` with `BBᵀ = R⁻¹H(HᵀR⁻¹H)⁻¹HᵀR⁻¹`, so `Q` is never formed; and the
nugget-free `R` is recovered from `R̃` on the diagonal instead of being kept separately.
The working set is **three** `n`-by-`n` arrays regardless of `p`.

Measured on this build (Octave 8.4, 2 cores, `p = 8`, Matérn 5/2, 2 restarts):

| n | full fit, pure `.m` | full fit, MEX | speed-up | one obj+grad `.m` | MEX | speed-up | peak `n×n` RAM `.m` → MEX |
|---|---|---|---|---|---|---|---|
| 250 | 5.13 s | 2.47 s | 2.1× | 0.0355 s | 0.0124 s | 2.9× | 5 MB → 1 MB |
| 500 | 20.1 s | 7.20 s | 2.8× | 0.148 s | 0.0393 s | 3.8× | 21 MB → 6 MB |
| 1000 | 70.5 s | 20.4 s | 3.5× | 0.801 s | 0.166 s | 4.8× | 84 MB → 23 MB |
| 2000 | 203 s | 37.7 s | 5.4× | 4.95 s | 0.790 s | 6.3× | 336 MB → 92 MB |

`bench_rgasp([250 500 1000 2000], 8)` reproduces the table. MATLAB's own elementwise
maths is faster than Octave's, so expect smaller ratios there — the memory column, and
the removal of the `p·n²` distance-matrix allocation, carry over unchanged.

The objective and the gradient agree between the two paths to machine precision
(`max |Δf| = 0`, `max |Δg| = 3.2e-14` over 216 configurations). Full *fits* can still
stop at slightly different points once `n` is large, because inert directions leave the
posterior almost flat and the line search then compounds a 1e-14 difference: at
`n = 500, p = 8` the two fits differ by 1% in `β̂` for the inert coordinates, `5e-3` in
log posterior, `1.6e-6` (in units of `sd(y)`) in predictive mean, and not at all in
RMSE or coverage to four figures. Both are valid modes on the same ridge.

### Why not the C++ from the official MATLAB package?

[RobustGaSP-in-Matlab](https://github.com/MengyangGu/RobustGaSP-in-Matlab) is **GPL-3**
(see its `LICENSE`), ships its C++ inside `RobustGaSP_Matlab.zip`, and requires the
Optimization Toolbox. Two reasons this port does not use it:

1. **Licensing.** GPL-3 is strongly copyleft. Linking that code in would make this whole
   toolbox GPL-3 on redistribution. The MEX here is written from the papers and carries
   no code from either the R or the MATLAB package.
2. **It is the slower algorithm.** Its gradient forms
   `W_l = (R⁻¹Ṙ_l)ᵀ − Ṙ_l·R⁻¹H(HᵀR⁻¹H)⁻¹HᵀR⁻¹` as a full `n`-by-`n` product for every
   input dimension — `O(p·n³)` per gradient evaluation. Precomputing `R⁻¹` once and
   reducing each dimension to traces is `O(n³) + O(p·n²)`. Transcribed into MATLAB so
   the languages cancel out, that algorithm runs 2.4× (`n = 100`) to 3.4× (`n = 1000`)
   slower than the one used here, before any C++ is involved.

## 7. Reproducing the paper

`demo_robustness` reproduces the paper's central point directly. Along the ray
`ξ = log(1/γ)` for the Lim et al. (2002) function embedded in four inputs
(two of which are inert), `n = 16`:

| ξ | 5 | 8 | 11 | 14 |
|---|---|---|---|---|
| profile likelihood (`mle`) | 32.8904 | 32.8904 | 32.8904 | 32.8904 |
| marginal likelihood (`mmle`) | 32.2210 | 32.2210 | 32.2210 | 32.2210 |
| marginal posterior (`post_mode`) | 615 | 1.18e4 | 2.36e5 | 4.74e6 |

Both likelihood criteria are **constant to machine precision** on the `γ → 0` plateau, so
their maximizer can sit anywhere on it — where `R → I` and the emulator collapses to its
fitted mean. The jointly robust prior grows without bound there, so the posterior mode is
always interior.

The consequence, over 20 repeated Latin hypercube designs (out-of-sample normalized RMSE,
`lowerBound` switched off so the estimators are unconstrained). Whether a given solver
actually walks onto the plateau is a property of the solver, so both are shown:

*Built-in projected L-BFGS (identical in every installation):*

| | `mle` | `mmle` | `post_mode` |
|---|---|---|---|
| median | 0.1219 | 0.0798 | **0.0717** |
| mean | 0.2862 | 0.0861 | **0.0766** |
| worst case | 1.0095 | 0.1554 | **0.1466** |
| # designs with nRMSE > 0.3 | **4** | 0 | 0 |

*Optimization Toolbox solver:*

| | `mle` | `mmle` | `post_mode` |
|---|---|---|---|
| median | 0.0987 | 0.0798 | **0.0717** |
| mean | 0.1033 | 0.0862 | **0.0767** |
| worst case | 0.1717 | 0.1554 | **0.1466** |
| # designs with nRMSE > 0.3 | 0 | 0 | 0 |

`fmincon`'s convergence tests stop it before it reaches the flat region, so it never
produces the catastrophic fits — but it is still worse than the posterior mode on average
and in the worst case, and nothing about the objective has changed: the maximizer of the
likelihood remains undefined on the plateau. `post_mode` is unaffected by the choice of
solver, which is the point.

## 8. Demos

| Script | What it shows |
|---|---|
| `demo_higdon_1d` | 1-d emulation, predictive band, leave-one-out diagnostics |
| `demo_borehole_8d` | 8-d emulation (n = 80): RMSE 0.34, nRMSE 0.008, coverage 0.94 on 500 held-out runs; the diagnostic flags `r`, `Tu`, `Tl` as inert |
| `demo_noisy_nugget` | noisy data, nugget estimated jointly (recovers noise variance 0.0085 vs true 0.01) |
| `demo_ppgasp_environ` | 200-output environmental spill field from 40 runs, nRMSE 0.002, coverage 0.97 |
| `demo_robustness` | the two panels described in §7 |

Each takes an optional output directory and writes a PNG there:
`demo_higdon_1d('figures')`.

The exact figures in the last digit depend on whether the MEX layer is compiled — see
the note at the end of §6 — but the emulator quality does not.

## 9. Tests

```matlab
addpath tests; ok = run_all_tests();
```

| Test | Coverage |
|---|---|
| `test_kernels` | closed forms, symmetry, unit diagonal, `d log c/dβ` vs finite differences, product structure, isotropic distances |
| `test_gradients` | analytic vs numerical gradient over 3 kernels × 3 methods × nugget on/off × zero-mean on/off × k ∈ {1,3} × 3 parameter vectors |
| `test_optimizer` | quadratic, Rosenbrock, an active lower bound, Nelder–Mead fallback |
| `test_toolbox_paths` | `tinv` vs the fallback at 6 degrees of freedom × 9 probabilities, with both checked against the t CDF; `norminv`; `pdist2` vs an explicit loop; the quantile definition; every available optimizer on one fit; `rgasp_simulate` centring and spread |
| `test_rgasp_fit_predict` | interpolation at design points; predictive mean/sd/intervals against an independent brute-force implementation; borehole accuracy and coverage; nugget recovery; `mmle`/`mle`; zero mean, linear trend, isotropic, all three kernels; `rgasp_simulate` centring |
| `test_ppgasp` | k = 1 reduces exactly to `rgasp`; field accuracy; interpolation; column-wise agreement with per-column fits at the shared range estimate |
| `test_loo_and_utils` | LOO against a refitted model, LOO calibration, inert-input detection on planted dummy inputs and on borehole, validation-metric definitions |
| `test_paper_behaviour` | flat likelihood plateau, interior posterior mode with vanishing gradient, robustness over repeated designs |

Result on Octave 8.4: **7 passed, 0 failed** (about 14 s).

## 10. Option reference

`rgasp(design, response, 'Name', Value, ...)` — R names in parentheses.

| Name | Default | Notes |
|---|---|---|
| `trend` | `ones(n,1)` | `h(x^D)`, n-by-q |
| `zeroMean` (`zero.mean`) | `false` | |
| `nugget` | `0` | fixed nugget-variance ratio η |
| `nuggetEst` (`nugget.est`) | `false` | estimate η jointly |
| `rangePar` (`range.par`) | `[]` | fix γ, skip optimization |
| `method` | `'post_mode'` | `'post_mode'` \| `'mmle'` \| `'mle'` |
| `a`, `b` | `0.2`, `(a+p)/n^(1/p)` | jointly robust prior |
| `kernelType` (`kernel_type`) | `'matern_5_2'` | string or 1-by-p cellstr |
| `alpha` | `1.9` | `pow_exp` only |
| `isotropic` | `false` | |
| `lowerBound` (`lower_bound`) | `true` | |
| `numInitialValues` (`num_initial_values`) | `2` | |
| `maxEval` (`max_eval`) | `max(30, 20+5p)` | |
| `initialValues` (`initial_values`) | `[]` | rows of starting `log β` (and `log η`) |
| `optimizer` (`optimization`) | `'auto'` | `'auto'` picks `fmincon`/`fminunc`; also `'lbfgs'`, `'neldermead'` |
| `optimTolFun` | `1e-8` | solver optimality tolerance |
| `optimTolX` | `1e-10` | solver step tolerance |
| `condNumUB` | `1e16` | conditioning target for the lower bound |
| `useMex` | `true` | set `false` to force the pure MATLAB code paths |

`rgasp_predict(model, testing_input, ...)`: `'trend'`, `'intervalData'` (default `true`,
adds η so the interval covers a noisy observation), `'level'` (default `0.95`).

## 11. Differences from the R package

The statistical model, defaults and reported quantities match RobustGaSP. Four
implementation choices differ deliberately:

1. **Matérn 3/2 derivative.** The R/C++ `matern_3_2_deriv` computes
   `−√3 d·R + √3 d·exp(−√3 β d)`. The second term is `∂c_l/∂β_l` rather than
   `R·(∂ log c_l/∂β_l)`, so it is correct only when `p = 1`; for a product correlation it
   omits the other dimensions' factors. This port uses `R ∘ (d log c_l/dβ_l)` throughout,
   which matches finite differences in every dimension (see `test_gradients`).
2. **Leave-one-out.** `leave_one_out_rgasp` in R forms `R_tilde = L Lᵀ + nugget`, which
   adds the nugget to *every* entry of the matrix (it is already on the diagonal of
   `L Lᵀ`). This port uses `L Lᵀ`.
3. **Lower-bound search.** R minimizes `(cond(R) − 10^16)²`, which is numerically flat
   over most of the search interval. This port minimizes
   `(log cond(R) − log 10^16)²` — the same root, a far better conditioned search.
4. **Reproducible pseudo-random starts.** Starting values beyond the second use a small
   built-in LCG (`rgasp_lcg`) instead of R's `set.seed`/`runif`, so results are identical
   across MATLAB and Octave.
5. **Lower-bound search on large designs.** `cond()` is an SVD, and the golden-section
   search calls it ~40 times; at `n` in the thousands that alone would dominate a fit.
   The search therefore runs on a subsample of at most 400 design rows. Identical to the
   full computation for `n ≤ 400`, and the bound is a safeguard rather than part of the
   model.

One further deviation is worth stating plainly. Whether the profile likelihood actually
*lands* on the degenerate plateau of §7 depends on the optimizer: the built-in projected
L-BFGS walks onto it on 4 of 20 designs, `fmincon`/`fminunc` on none of them. The
plateau — the property the paper proves — is a fact about the objective and is identical
either way; the failure count is not. The demo and the test therefore report both and
make their assertions against the built-in solver, which behaves identically everywhere.

Not ported: the `periodic_gauss` / `periodic_exp` kernels, and the exact reference priors
(`ref_xi`, `ref_gamma`) — the jointly robust prior (`ref_approx`) is the package default
and the one the paper recommends. `mmle` and `mle` are available for comparison.

## 12. Files

```
setup_robustgasp.m          add the toolbox to the path
rgasp.m  ppgasp.m           model fitting
rgasp_predict.m  ppgasp_predict.m
rgasp_simulate.m  rgasp_loo.m  rgasp_inert_inputs.m
rgasp_validate.m  rgasp_summary.m
rgasp_lhs.m                 Latin hypercube designs (prefers lhsdesign)
rgasp_compile_mex.m         build the optional C++ acceleration
rgasp_toolbox_report.m      what is available and what is being used
rgasp_use_toolboxes.m       force the built-in fallbacks on or off
internal/                   estimation engine
  gasp_fit.m                shared fitting driver
  gasp_predict.m            shared predictive distribution
  rgasp_objective.m         negative log marginal posterior + gradient
  rgasp_corr.m  rgasp_corr_1d.m  rgasp_dlogcorr_1d.m  rgasp_R0.m
  rgasp_lbfgs.m  rgasp_neldermead.m  rgasp_search_lb.m
  rgasp_options.m  rgasp_tinv.m  rgasp_norminv.m
  rgasp_corr_inputs.m       MEX / MATLAB dispatcher for the correlation
  rgasp_have_mex.m  rgasp_kernel_code.m  rgasp_toolbox.m
  rgasp_optimoptions.m      fmincon/fminunc options (optimoptions or optimset)
  rgasp_quantile.m
  rgasp_lhs_det.m  rgasp_lcg.m  rgasp_mute_warnings.m
  mex/                      rgasp_corr_mex.cpp, rgasp_gradterms_mex.cpp
testfun/                    higdon_1_data, limetal_2_data, dettepepel_3_data,
                            friedman_5_data, borehole (+ borehole_ranges),
                            environ_4_data
demos/                      five demo scripts
tests/                      run_all_tests + nine test files + bench_rgasp
figures/                    PNGs produced by the demos
```

## 13. References

- Gu, M., Wang, X., Berger, J. O. (2018). Robust Gaussian stochastic process emulation.
  *Annals of Statistics* 46(6A), 3038–3066. arXiv:1708.04738.
- Gu, M., Palomo, J., Berger, J. O. (2019). RobustGaSP: Robust Gaussian stochastic process
  emulation in R. *The R Journal* 11(1), 112–136.
- Gu, M. (2019). Jointly robust prior for Gaussian stochastic process in emulation,
  calibration and validation. *Bayesian Analysis* 14(3), 877–905.
- Gu, M., Berger, J. O. (2016). Parallel partial Gaussian process emulation for computer
  models with massive output. *Annals of Applied Statistics* 10(3), 1317–1347.

The R package RobustGaSP is GPL-2 licensed; this is an independent MATLAB implementation
written from the papers and the published algorithm descriptions.
