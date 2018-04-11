#include "mex.h"
#include <stdint.h>

/* JE_ADD_HESSIAN_LAYER
 * Adds a layer of triangular components of the Hessian matrix.
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /*
   * mat_out: trihess
   * mat_in1 = AtQ2{j} or R2_j * R2_j'
   * mat_in2 = B(j, :)' * B(j, :) or AtR2inv{j}
   */
  double    *mat_out, *mat_in1, *mat_in2;
  double    A_ij;
  mwSize    mr, p_j, r, mrr, hidx[3], kidx;
  uint32_t *rng;
  
  if (nlhs > 0) {
    mexErrMsgTxt("This function does not create any output.");
  }
  
  if (nrhs != 4) {
    mexErrMsgTxt("This function requires 4 inputs.");
  }
  
  // Pointers to the input variables.
  mat_out = mxGetPr(prhs[0]);
  mat_in1 = mxGetPr(prhs[1]);
  mat_in2 = mxGetPr(prhs[2]);
  rng = (unsigned int*) mxGetData(prhs[3]);
  
  // Get rank and nnz.
  mr = mxGetM(prhs[0]);
  r = mxGetM(prhs[2]);
  mrr = mr * r;
  p_j = mxGetM(prhs[3]);
  
  // Compute the triangular elements.
  for (mwSize j = 0; j < p_j; j++) {
    hidx[0] = mrr * rng[j];
    for (mwSize i = 0; i <= j; i++) {
      A_ij = mat_in1[i * p_j + j];
      hidx[1] = hidx[0] + r * rng[i];
      for (mwSize l = 0; l < r; l++) {
        hidx[2] = hidx[1] + l * mr;
        kidx = l * r;
        for (mwSize k =0; k <= l; k++) {
          mat_out[hidx[2] + k] += A_ij * mat_in2[kidx + k];
        }
      }
    }
  }
}