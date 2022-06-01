#include "dec.h"
using namespace hls;

int main ()
{

	matrix_t array[row][column] = {{1,1,1,1,1,1,1,1,1,0},{1,1,1,1,1,1,1,1,0,1},{1,1,1,1,1,1,1,0,1,1},{1,1,1,1,1,1,0,1,1,1},{1,1,1,1,1,0,1,1,1,1},
									{1,1,1,1,1,0,1,1,1,1},{1,1,1,0,1,1,1,1,1,1},{1,1,0,1,1,1,1,1,1,1},{1,0,1,1,1,1,1,1,1,1},{0,1,1,1,1,1,1,1,1,1}};

	printf("                                           OUR INPUT MATRIX\n");
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
		{
			printf("%1.4f ",array[i][j]);
		}
		printf("\n");
	}
	decompose(array);
	return 0;
}
