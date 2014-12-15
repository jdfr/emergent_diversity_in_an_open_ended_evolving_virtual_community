#include "mex.h"
#include "matrix.h"

#define MIN3(a,b,c) (a < b ? (a < c ? a : c) : (b < c ? b : c))
 
int levenshtein(short int *s1, size_t l1,
    short int *s2, size_t l2) {
  int i, j;
  size_t len = (l1 + 1) * (l2 + 1);
  /*char*/short int *p1, *p2;
  unsigned int d1, d2, d3, *d, *dp, res;
 
  if (l1 == 0) {
    return l2;
  } else if (l2 == 0) {
    return l1;
  }
 
  d = (unsigned int*)mxMalloc(len * sizeof(unsigned int));
 
  *d = 0;
  for(i = 1, dp = d + l2 + 1;
      i < l1 + 1;
      ++i, dp += l2 + 1) {
    *dp = (unsigned) i;
  }
  for(j = 1, dp = d + 1;
      j < l2 + 1;
      ++j, ++dp) {
    *dp = (unsigned) j;
  }
 
  for(i = 1, p1 = s1, dp = d + l2 + 2;
      i < l1 + 1;
      ++i, p1++, ++dp) {
    for(j = 1, p2 = s2;
        j < l2 + 1;
        ++j, p2++, ++dp) {
       /*mexPrintf("Comp (%d,%d)-(%d,%d)... ", i, *p1, j, *p2);*/
      if((*p1)==(*p2)) {/*!comp(p1, p2)) {*/
        *dp = *(dp - l2 - 2);
        /*mexPrintf("!=: %d\n", *dp);*/
      } else {
        d1 = *(dp - 1) + 1;
        d2 = *(dp - l2 - 1) + 1;
        d3 = *(dp - l2 - 2) + 1;
        *dp = MIN3(d1, d2, d3);
        /*mexPrintf("==: %d,%d,%d -> %d\n", d1, d2, d3, *dp);*/
      }
    }
  }
  res = *(dp - 2);
 
  dp = NULL;
  mxFree(d);
  return res;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  short int *v1, *v2;
  size_t nin;
  int res, nm1, nm2, i,z;
  double ress;
  mxArray *cell;
  
  /*mexPrintf("Mira char: %d\n", sizeof(char));
  mexPrintf("Mira shortint: %d\n", sizeof(short int));*/
  v1 = (short int *) mxGetPr(prhs[0]);
  v2 = (short int *) mxGetPr(prhs[1]);
  nm1 = mxGetN(prhs[0])*mxGetM(prhs[0]);
  nm2 = mxGetN(prhs[1])*mxGetM(prhs[1]);
  res = levenshtein(v1, nm1, v2, nm2);
  ress = res;
  plhs[0] = mxCreateDoubleScalar(ress);
}
