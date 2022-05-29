#include "dec.h"
using namespace hls;

int main ()
{

	matrix_t array[Inp_matrix] = { 3,2,1,4};
	matrix_t tem_input_2d_array[matrix_size][matrix_size] = {0};
	tem_input_2d_array[0][0]=array[0];
	tem_input_2d_array[0][1]=array[1];
	tem_input_2d_array[1][0]=array[2];
	tem_input_2d_array[1][1]=array[3];

	printf("                                           OUR INPUT MATRIX\n");
	for(int i=0;i<matrix_size;i++)
	{
		for(int j=0;j<matrix_size;j++)
		{
			printf("%1.4f ",tem_input_2d_array[i][j]);
		}
		printf("\n");
	}
	decomposition(array);
	return 0;
}
