#include "mex.h"
#include "matrix.h"

void calculateChars(const mxArray *chararray, int results[]) {
    mxChar *chars;
    int i, nchars;
    
    for (i=0;i<5;i++) {
      results[i]=0;
    }
    
    nchars = (int)(mxGetN(chararray)*mxGetM(chararray));
    chars  = mxGetChars(chararray);
    
    for (i=0;i<nchars;i++) {
      switch ((char)(*chars)) {
        case 'G':
          results[0]++;
          break;
        case '[':
          results[1]++;
          break;
        case '+':
          results[2]++;
          break;
        case '-':
          results[3]++;
      }
      chars++;
    }
    results[4] = nchars;

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int nchars, i, j, jj, *out, *res, ress[5], *allres;/*, NG, NB, NP, NM;*/
    mxArray *thiscell;

    if (mxIsChar(prhs[0])) {
      plhs[0] = mxCreateNumericMatrix(1, 5, mxINT32_CLASS, mxREAL);
      out = (int*)mxGetData(plhs[0]);
      calculateChars(prhs[0], out);
    } else if (mxIsCell(prhs[0])) {
      nchars = (int)(mxGetN(prhs[0])*mxGetM(prhs[0]));
      plhs[0] = mxCreateNumericMatrix(nchars, 5, mxINT32_CLASS, mxREAL);
      out = (int*)mxGetData(plhs[0]);
      for (i=0;i<nchars;i++) {
        thiscell = mxGetCell(prhs[0], i);
        calculateChars(thiscell, ress);
        jj = i;
        for (j=0;j<5;j++) {
          out[jj] = ress[j];
          jj += nchars;
        }
      }
    } else {
      mexErrMsgTxt("Input of incorrect type!!!!");
    }
}


