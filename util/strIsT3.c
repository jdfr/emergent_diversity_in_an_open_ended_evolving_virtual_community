#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  short int *cad, *init, c;
  size_t nin;
  int i, anid, nm;
  double ress;
  mxArray *cell;
  
  init = cad = (short int *) mxGetPr(prhs[0]);
  nm = mxGetN(prhs[0])*mxGetM(prhs[0]);
  anid = 0;
  for(i=0; i<nm; i++) {
    c = (*(cad++));
    anid += (c=='[') - (c==']');
    if (anid<0) {
      mexPrintf("ERROR: BAD GENOME: %s\n", init);
      mexErrMsgTxt("BAD GENOME!!!");
    }
    if (anid>1) {
      plhs[0] = mxCreateLogicalScalar(false);
      return;
    }
  }
  plhs[0] = mxCreateLogicalScalar(true);
}
