#include <math.h>
#include "mex.h"

//pack each byte from a logical image into a bit in a reduced image
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *logic_img;
    mxLogical *input_ptr;
    mxArray *packed_img;
    uint8_T *output_ptr, *output_buffer;
    mwSize rows, row, cols, col, packed_rows, outsz[2];
    int cuentabits;

    if (nrhs != 1) {
        return;
    }
    logic_img      = prhs[0];
    rows           = mxGetM(logic_img);
    cols           = mxGetN(logic_img);
    packed_rows    = (rows+7)/8;
    outsz[0]       = packed_rows;
    outsz[1]       = cols;
    packed_img     = mxCreateNumericArray(2, outsz, mxUINT8_CLASS, mxREAL);
    if (packed_img == NULL) {
        return;
    }
    plhs[0]        = packed_img;

    input_ptr      = (mxLogical *) mxGetData(logic_img);        
    output_buffer  = (uint8_T *) mxGetData(packed_img);

    for (col = 0; col < cols; col++) {
        cuentabits = 0;
        output_ptr = output_buffer + packed_rows*col;
        for (row = 0; row < rows; row++) {
            if (*input_ptr != 0) {
                *output_ptr |= 1 << cuentabits;
            }
            cuentabits = (cuentabits+1)%8;
            if (cuentabits == 0) {
                output_ptr++;
            }
            input_ptr++;
        }
    }


}
