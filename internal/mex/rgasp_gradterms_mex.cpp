/* rgasp_gradterms_mex.cpp
 *
 * Fused evaluation of the two per-dimension gradient terms of the log
 * marginal posterior.  With Rdot_l = R o (d log c_l / d beta_l) and
 * Q = R~^{-1} - B B'  (B B' = R~^{-1}X (X'R~^{-1}X)^{-1} X'R~^{-1}):
 *
 *     tr(l)   = tr(Rdot_l Q)   = sum_ij Ri_ij Rdot_ij
 *                              - sum_c sum_ij B_ic Rdot_ij B_jc
 *     quad(l) = sum_c (QY_c' Rdot_l QY_c) / S2_c
 *
 *   [tr, quad] = rgasp_gradterms_mex(input, Rt, nugget, Ri, B, QY, S2, ...
 *                                    beta, ktcode, alpha, isotropic, useTrend)
 *
 * Rdot_l is never materialized and the n-by-n arrays Rt and Ri are streamed
 * once in total rather than once per input dimension: the loop over l is the
 * innermost one.  Neither the low-rank trend correction B B' nor the
 * projection matrix Q is ever formed as an n-by-n matrix.
 *
 * Rt carries the nugget on its diagonal; the nugget-free correlation needed
 * for Rdot is recovered by subtracting it there.
 *
 * Part of RobustGaSP-MATLAB.  Written from the published algorithms; contains
 * no code from the R or MATLAB RobustGaSP packages.
 */

#include "mex.h"
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

static const double SQ5 = 2.2360679774997896964;
static const double SQ3 = 1.7320508075688772935;

/* d log c / d beta, in forms that stay finite for large u */
static inline double dlogcorr1d(double d, double b, int code, double a)
{
    if (code == 3) { const double u = SQ5*b*d; return -SQ5*d*(u*(1.0+u))/(3.0 + 3.0*u + u*u); }
    if (code == 2) { const double u = SQ3*b*d; return -SQ3*d*u/(1.0+u); }
    if (b <= 0.0) return 0.0;
    return -a*std::pow(b, a-1.0)*std::pow(d, a);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 12) mexErrMsgIdAndTxt("rgasp:gradmex",
        "Usage: [tr, quad] = rgasp_gradterms_mex(input, Rt, nugget, Ri, B, QY, S2, "
        "beta, ktcode, alpha, isotropic, useTrend)");

    const double *X   = mxGetPr(prhs[0]);
    const double *Rt  = mxGetPr(prhs[1]);
    const double  nu  = mxGetScalar(prhs[2]);
    const double *Ri  = mxGetPr(prhs[3]);
    const double *B   = mxGetPr(prhs[4]);
    const double *QY  = mxGetPr(prhs[5]);
    const double *S2  = mxGetPr(prhs[6]);
    const double *beta  = mxGetPr(prhs[7]);
    const double *kt    = mxGetPr(prhs[8]);
    const double *alpha = mxGetPr(prhs[9]);
    const bool iso      = (mxGetScalar(prhs[10]) != 0.0);
    const bool useTrend = (mxGetScalar(prhs[11]) != 0.0);

    const mwSize n  = mxGetM(prhs[1]);
    const mwSize pd = mxGetN(prhs[0]);                 /* design dimension  */
    const mwSize pe = iso ? 1 : pd;                    /* effective ranges  */
    const mwSize k  = mxGetN(prhs[5]);                 /* output columns    */
    const mwSize q  = (useTrend && !mxIsEmpty(prhs[4])) ? mxGetN(prhs[4]) : 0;

    plhs[0] = mxCreateDoubleMatrix(pe, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(pe, 1, mxREAL);
    double *trOut = mxGetPr(plhs[0]);
    double *qdOut = mxGetPr(plhs[1]);

    std::vector<double> trAcc(pe, 0.0);
    std::vector<double> qdAcc((size_t)pe*k, 0.0);      /* per (l, column)   */
    std::vector<double> tvAcc((size_t)pe*(q ? q : 1), 0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<double> trL(pe, 0.0);
        std::vector<double> qdL((size_t)pe*k, 0.0);
        std::vector<double> tvL((size_t)pe*(q ? q : 1), 0.0);
        std::vector<double> dl(pe);

#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
        for (long j = 0; j < (long)n; ++j) {
            const double *Rc = Rt + (size_t)j*n;
            const double *Ic = Ri + (size_t)j*n;
            for (mwSize i = 0; i < n; ++i) {

                double Rij = Rc[i];
                if (i == (mwSize)j) Rij -= nu;          /* nugget-free value */
                if (Rij == 0.0) continue;

                if (iso) {
                    double s = 0.0;
                    for (mwSize l = 0; l < pd; ++l) {
                        const double t = X[(size_t)l*n + i] - X[(size_t)l*n + j];
                        s += t*t;
                    }
                    dl[0] = dlogcorr1d(std::sqrt(s), beta[0], (int)kt[0], alpha[0]);
                } else {
                    for (mwSize l = 0; l < pe; ++l) {
                        const double d = std::fabs(X[(size_t)l*n + i] - X[(size_t)l*n + j]);
                        dl[l] = dlogcorr1d(d, beta[l], (int)kt[l], alpha[l]);
                    }
                }

                const double Iij = Ic[i];
                for (mwSize l = 0; l < pe; ++l) {
                    const double dR = Rij * dl[l];
                    trL[l] += Iij * dR;
                    for (mwSize c = 0; c < k; ++c)
                        qdL[(size_t)c*pe + l] += QY[(size_t)c*n + i] * dR * QY[(size_t)c*n + j];
                    for (mwSize c = 0; c < q; ++c)
                        tvL[(size_t)c*pe + l] += B[(size_t)c*n + i] * dR * B[(size_t)c*n + j];
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (mwSize l = 0; l < pe; ++l) trAcc[l] += trL[l];
            for (size_t t = 0; t < qdL.size(); ++t) qdAcc[t] += qdL[t];
            for (size_t t = 0; t < tvL.size(); ++t) tvAcc[t] += tvL[t];
        }
    }

    for (mwSize l = 0; l < pe; ++l) {
        double tr = trAcc[l];
        for (mwSize c = 0; c < q; ++c) tr -= tvAcc[(size_t)c*pe + l];
        trOut[l] = tr;
        double s = 0.0;
        for (mwSize c = 0; c < k; ++c) s += qdAcc[(size_t)c*pe + l] / S2[c];
        qdOut[l] = s;
    }
}
