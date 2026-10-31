# Derivation ledger: `xtpqardl()` (xtpqardl 1.0.5)

Copied from the audit ledger of 1.0.4 (`/home/claude/rel/audit/xtpqardl/ledger_xtpqardl.md`) and updated to the corrected code. Line numbers refer to `R/xtpqardl.R` of 1.0.5. Papers: Cho, Kim and Shin (2015) JoE 188 (CKS); Koenker and Bassett (1982) Econometrica 50(1), 43-61 (KB); Pesaran, Shin and Smith (1999) JASA 94 (PSS); Hausman (1978) Econometrica 46(6). "Verification" names the testthat test (tests/testthat/test-fixes-1-0-5.R unless stated) that pins the computation to a hand calculation or oracle, with its tolerance.

## A. Model and estimating regression

| Paper equation / page | Mathematical expression | Code file:line | Convention adopted | Source of convention | Verification status |
|---|---|---|---|---|---|
| CKS eq. (6) (ECM form): dY_t = alpha + zeta (Y_{t-1} - beta' X_{t-1}) + sum phi_j dY_{t-j} + sum theta_j' dX_{t-j} + U_t(tau) | regression of d_y on (lr_y, lr_x, L1..L(p-1) d_y, d_x lags 0..q-1, constant); rho = coef(lr_y), theta = coef(lr_x) | xtpqardl.R:587-626 (`.build_panel_regression`) | unrestricted ARDL in ECM form (PSS (2), CKS (4)) | paper | unchanged from 1.0.4 (matched `rq.fit` at 0 in the audit); re-pinned: "DFE joint covariance ... hand-built" compares `dfe_coef` to `rq.fit` on the hand design at 1e-10; "lags are built by the time index" compares rho and short run to `rq.fit` at 1e-8 (DFE and MG) |
| CKS eq. (3): beta(tau) = gamma(tau)(1 - sum phi)^{-1} | beta = -theta / rho | xtpqardl.R:1064 (DFE), 839 (MG/PMG) | identical to CKS plug-in | paper | "DFE joint covariance ... hand-built" at 1e-10 |
| Lag counting | lags of d_y: 1..p-1; of each d_x: 0..q_j-1; `lr` unlagged | xtpqardl.R:603-615 | PSS counting; user supplies the lagged levels | package choice, documented | by construction; lags now by time (Section E) |
| Constant | MG/PMG: column of ones; DFE: intercept + N-1 panel contrasts (first panel reference) | xtpqardl.R:617-620, 1024-1028 | DFE always has panel intercepts | package choice | same column space as one dummy per panel (audit: slopes equal at 5e-16); the contrast parameterisation is now the one estimated, so `all = TRUE` in `wald_test` restricts exactly quantreg's set |
| Effective sample | complete rows per panel after lagging; `n_dropped` counts the removed rows | xtpqardl.R:622-626 | listwise deletion per panel, documented (the leading rows lost to the user's differencing are the usual case) | package choice | "lags are built by the time index": n_obs equals the hand complete-case count |
| Formula terms | `model.matrix` on the sorted data with `na.pass`; expanded columns attached as numeric variables | xtpqardl.R:539-547 | log(), I(), a:b allowed; `lr` plain column names | FIX_RULES | "collinear regressors, transformed formulas ...": log(xpos), I(d_x^2) fit equals `rq.fit` on the hand-transformed design at 1e-8 |

## B. Dynamic fixed effects (DFE): covariances

| Paper equation / page | Mathematical expression | Code file:line | Convention adopted | Source of convention | Verification status |
|---|---|---|---|---|---|
| KB (1982) Theorem 3.1, p. 46: covariance Omega(u, F) (x) D^{-1}, omega_ij = [min(u_i,u_j) - u_i u_j]/(f(Q(u_i)) f(Q(u_j))), D = lim n^{-1} X'X, iid errors (A1-A3) | Cov(b(tau_a), b(tau_b)) = (min(tau_a,tau_b) - tau_a tau_b) s_a s_b (X'X)^{-1}, s_tau the sparsity of `summary.rq(se = "iid")` (1/`$scale`); written as w H_a^{-1} (X'X) H_b^{-1} with H_a^{-1} = s_a (X'X)^{-1} | xtpqardl.R:1115-1161 (`.dfe_joint_cov`), default `dfe_cov = "iid"`; blocks assembled at 1155 | the paper's iid form; one estimator for every block (diagonal and cross); audit size 4.3 / 5.0 / 7.0 percent | paper (KB Theorem 3.1); the sparsity estimator is quantreg's | "DFE joint covariance ... hand-built": full pm x pm matrix equals the hand formula with the sparsity recipe written out at 1e-10; diagonal block equals `summary.rq(se = "iid")$cov` at 1e-10; rho_V and beta_V equal the hand blocks and delta transform at 1e-10 |
| "nid" and "ker" options | w H_a^{-1} (X'X) H_b^{-1} with H_a^{-1} = `summary.rq(se = type)$Hinv` | xtpqardl.R:1141, 1155 | quantreg's density estimates and bandwidth rule (h halved while tau +/- h leaves (0,1)); documented warning on the audit's under-statement of variance with panel dummies | package choice following quantreg | "nid and ker DFE covariances are quantreg's ...": every block equals w Hinv_a J Hinv_b from `summary.rq` at 1e-10; diagonal equals `$cov` |
| Delta method for beta = -theta/rho: G = [theta/rho^2, -(1/rho) I_k] | V_beta(a,b) = G_a C_ab G_b' for every pair (a, b) | xtpqardl.R:1067-1082 | CKS (10) linearised; cross blocks included | paper | hand-built test at 1e-10 |
| Mixed matrix of 1.0.4 | removed | | never mix estimators | FIX_RULES | the hand test checks the cross block uses the same sparsities as the diagonal |
| Failure of the estimator at a quantile | blocks NA, `panel_cov_status$ok = FALSE`, one warning naming tau and reason | xtpqardl.R:1046-1057, 1130-1149 | NA with warning, no substitute | FIX_RULES | "collinear regressors, transformed formulas, and NA covariance blocks": exact linear y gives NA rho_V with the warning, Wald and KB tests NA |

## C. Mean group (MG) and pooled mean group (PMG)

| Paper equation / page | Mathematical expression | Code file:line | Convention adopted | Source of convention | Verification status |
|---|---|---|---|---|---|
| PSS (1999) MG: average of panel estimates; variance sum_i (b_i - b_bar)(b_i - b_bar)'/(N(N-1)) | `colMeans`; `crossprod(D)/(n(n-1))` | xtpqardl.R:859-867 | denominator N(N-1); N < 3 warns | paper | unchanged (audit: 0); test-xtpqardl.R "MG equals the average"; "MG with fewer than three panels warns" |
| Panel covariance for the PMG weights | `summary.rq(se = "iid")`, then "nid", "ker", then `.iid_cov_fallback`; a switch is reported once | xtpqardl.R:633-660, 849-857 | fallback chain kept, no longer silent | package choice | warning text pinned in test-wald-pmg-inputs.R "PMG on very short panels" |
| PMG pooling beta_P = (sum W_i)^{-1} sum W_i beta_i, V_P = (sum W_i)^{-1} | same; panels without a usable covariance are named in a warning | xtpqardl.R:884-929 | minimum-distance analogue of PSS, documented | package choice | unchanged (audit: 1.1e-16) |
| PMG stage 2 | ECT = lr_y - beta_P' lr_x imposed, rho mean group | xtpqardl.R:931-954 | as PSS | paper | unchanged (audit: 1.1e-16) |
| Hausman (1978): h = d'(V_MG - V_PMG)^{-1} d | per quantile; strict: NA with warning unless positive definite (eigenvalues > 1e-8 max); ginv: Moore-Penrose over the kept eigenvalues, df = rank, warning | xtpqardl.R:693-738 (`.hausman_one`), 956-966 | same rule as ardlverse 2.1.1; the 1.0.4 statistic set the unknown PMG cross-quantile blocks to zero and dropped negative eigenvalues silently | FIX_RULES | "Hausman test is NA with a warning ...": zero difference (equal weights) gives NA with the warning; hand indefinite matrix diag(2, -1), d = (1, 1): strict NA, ginv statistic 1/2 with df 1; positive definite case equals d^2/(V_MG - V_PMG) at 1e-10 |

## D. Half-life, IRF, lag selection

| Paper equation / page | Mathematical expression | Code file:line | Convention adopted | Source of convention | Verification status |
|---|---|---|---|---|---|
| Half-life ln(0.5)/ln(1 + rho), -1 < rho < 0 | same | xtpqardl.R:472-478 | NA outside (-1, 0) | package choice | test-xtpqardl.R |
| IRF (1 + rho)^t | same | wald_test.R:537 | computed for -2 < rho < 0; the oscillating range is now documented | package choice | by reading |
| Lag selection by BIC/AIC | least squares per panel, pooled RSS and n; penalty counts all coefficients: K = sum over panels of ncol(X) | xtpqardl.R:1164-1220, line 1199 | documented as least squares on the mean regression | package choice | "lag selection penalises all estimated coefficients": K = 4N for p = q = 1 and 6N for p = q = 2; IC equals n log(RSS/n) + K log n from `lm` per panel at 1e-10 |

## E. Missing values, unbalanced panels, inputs

| Item | Code file:line | Behaviour | Verification |
|---|---|---|---|
| Time index | `time` numeric; unit step = smallest positive within-panel difference; integer index `.tidx` | xtpqardl.R:488-516 | non-numeric time is an error ("lags are built by the time index" test) |
| Duplicated (id, time) | error naming the pairs | xtpqardl.R:500-506 | same test |
| Gap inside a panel | error unless `allow_gaps = TRUE`; lags matched on `.tidx - j` | xtpqardl.R:517-535, 557-560, 605-612 | same test: with a gap at t = 15 and q = 2 the fit equals `rq.fit` on the hand time-aligned design at 1e-8, n equal; the row-position lag would differ |
| Interior NA value | row and the rows whose lags need it dropped; count in `n_dropped` | xtpqardl.R:622-626 | by construction (match returns NA) |
| Exactly collinear regressors | error naming the columns (pooled DFE; panel for MG/PMG) | xtpqardl.R:566-577, 775, 1032 | "collinear regressors ..." test |
| Formula with log(), I(), a:b | expanded with model.matrix | xtpqardl.R:539-547 | same test, hand `rq.fit` at 1e-8 |
| Global side effects | none (no set.seed, no options) | grep | verified |
