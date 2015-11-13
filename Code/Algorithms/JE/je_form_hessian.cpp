#include "mex.h"

/* JE_FORM_HESSIAN
 * Adds a layer of triangular components of the Hessian matrix.
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /*
   * Input 1: trihess
   * Input 2: r
   */
  double    *hess, *rankfill;
  mwSize    m, r, mr, mrr, hidx[6];
  mwSize    i, j, k, l;
  mwSize    *rng;
  
  if (nlhs > 0) {
    mexErrMsgTxt("This function does not create any output.");
  }
  
  if (nrhs < 2) {
    mexErrMsgTxt("This function requires at least 2 inputs.");
  } else if (nrhs > 3) {
    mexErrMsgTxt("This function requires at most 3 inputs.");
  }
  
  // Pointers to the input variables.
  hess = mxGetPr(prhs[0]);
  
  // Get rank and nnz.
  mr = mxGetM(prhs[0]);
  r = mxGetScalar(prhs[1]);
  m = mr / r;
  mrr = mr * r;
  
  // If Okatani's constraint is added, add it to current Hessian.
  if (nrhs == 3) {
    rankfill = mxGetPr(prhs[2]);
    for (j = 0; j < m; j++) {
      hidx[0] = j * mrr;
      hidx[1] = j * m;
      for (i = 0; i <= j; i++) {
        hidx[2] = hidx[0] + i * r;
        hidx[3] = hidx[1] + i;
        for (k = 0; k < r; k++) {
          hess[hidx[2] + k * mr + k] += rankfill[hidx[3]];
        }
      }
    }
  }
  
  // Fill the lower triangular components of the upper triangluar Hessian blocks.
  for (i = 0; i < m; i++) {
    hidx[0] = i * mrr;
    for (j = 0; j <= i; j++) {
      hidx[1] = hidx[0] + j * r;
      if (hess[hidx[1] + mr] == 0) {
        continue;
      } else {
        for (k = 1; k < r; k++) {
          hidx[2] = hidx[1] + k;
          hidx[3] = hidx[1] + k * mr;
          for (l = 0; l < k; l++) {
            hess[hidx[2] + l * mr] = hess[hidx[3] + l];
          }
        }
      }
    }
  }
  
  // Fill the lower tringular Hessian blocks.
  for (i = 1; i < m; i++) {
    hidx[0] = i * r;
    hidx[1] = i * mrr;
    for (j = 0; j < i; j++) {
      hidx[2] = hidx[0] + j * mrr;
      hidx[3] = hidx[1] + j * r;
      if (hess[hidx[3] + mr] == 0) {
        if (nrhs == 2) {
          continue;
        } else {
          for (k = 0; k < r; k++) {
            hess[hidx[2] + k * mr + k] = hess[hidx[3] + k * mr + k];
          }
        }
      } else {
        for (k = 0; k < r; k++) {
          hidx[4] = hidx[2] + k;
          hidx[5] = hidx[3] + k * mr;
          for (l = 0; l < r; l++) {
            // hess[i * r + j * mrr + l * mr + k] = hess[i * mrr + j * r + k * mr + l];
            hess[hidx[4] + l * mr] = hess[hidx[5] + l];
          }
        }
      }
    }
  }
}