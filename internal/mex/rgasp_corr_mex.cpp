/* rgasp_corr_mex.cpp
 *
 * Separable (product) correlation matrix, computed directly from the design
 * matrices so that the p per-dimension distance matrices are never
 * materialized.
 *
 *   R = rgasp_corr_mex(input1, input2, beta, ktcode, alpha, isotropic)
 *
 *   input1     n1-by-p design (rows of the output)
 *   input2     n2-by-p design (columns of the output)
 *   beta       p-vector (or scalar when isotropic) of inverse ranges
 *   ktcode     p-vector: 1 = pow_exp, 2 = matern_3_2, 3 = matern_5_2
 *   alpha      p-vector of roughness parameters (pow_exp only)
 *   isotropic  logical scalar; when true a single beta acts on the Euclidean
 *              distance across all p coordinates
 *
 * The whole product over dimensions is evaluated element by element, so the
 * n1-by-n2 result is written exactly once instead of p times.
 *
 * Part of RobustGaSP-MATLAB.  Written from the published algorithms; contains
 * no code from the R or MATLAB RobustGaSP packages.
 */

#include "mex.h"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

static const double SQ5 = 2.2360679774997896964;
static const double SQ3 = 1.7320508075688772935;

static inline double corr1d(double d, double b, int code, double a)
{
    if (code == 3) { const double u = SQ5*b*d; return (1.0 + u + u*u/3.0)*std::exp(-u); }
    if (code == 2) { const double u = SQ3*b*d; return (1.0 + u)*std::exp(-u); }
    return std::exp(-std::pow(b*d, a));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6) mexErrMsgIdAndTxt("rgasp:corrmex",
        "Usage: R = rgasp_corr_mex(input1, input2, beta, ktcode, alpha, isotropic)");

    const double *x1 = mxGetPr(prhs[0]);
    const double *x2 = mxGetPr(prhs[1]);
    const double *beta  = mxGetPr(prhs[2]);
    const double *kt    = mxGetPr(prhs[3]);
    const double *alpha = mxGetPr(prhs[4]);
    const bool iso = (mxGetScalar(prhs[5]) != 0.0);

    const mwSize n1 = mxGetM(prhs[0]);
    const mwSize n2 = mxGetM(prhs[1]);
    const mwSize p  = mxGetN(prhs[0]);
    if (mxGetN(prhs[1]) != p)
        mexErrMsgIdAndTxt("rgasp:corrmex", "input1 and input2 need the same number of columns.");

    plhs[0] = mxCreateDoubleMatrix(n1, n2, mxREAL);
    double *R = mxGetPr(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long j = 0; j < (long)n2; ++j) {
        double *Rc = R + (size_t)j*n1;
        for (mwSize i = 0; i < n1; ++i) {
            if (iso) {
                double s = 0.0;
                for (mwSize l = 0; l < p; ++l) {
                    const double t = x1[(size_t)l*n1 + i] - x2[(size_t)l*n2 + j];
                    s += t*t;
                }
                Rc[i] = corr1d(std::sqrt(s), beta[0], (int)kt[0], alpha[0]);
            } else {
                double v = 1.0;
                for (mwSize l = 0; l < p; ++l) {
                    const double d = std::fabs(x1[(size_t)l*n1 + i] - x2[(size_t)l*n2 + j]);
                    v *= corr1d(d, beta[l], (int)kt[l], alpha[l]);
                }
                Rc[i] = v;
            }
        }
    }
}
