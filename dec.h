#include "hls_stream.h"
#include "hls_linear_algebra.h"
#include "hls_math.h"
//#include "hls_matrix_multiply.h"
#include "math.h"
#include "stdio.h"
#include "stdint.h"
#define matrix_size 10
#define row 10
#define column 10
#define Inp_matrix 9
using namespace hls;

typedef float matrix_t;
void decompose(matrix_t input[row][column]);


