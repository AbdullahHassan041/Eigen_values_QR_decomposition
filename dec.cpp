#include <iostream>
#include "dec.h"
#include "hls/linear_algebra/utils/x_hls_matrix_utils.h"
//#include "hls/linear_algebra/utils/x_hls_matrix_tb_utils.h"
#include "hls_stream.h"

using namespace hls;

static matrix_t A_cal[matrix_size][matrix_size]={0};
static matrix_t A_cal_temp2[Inp_matrix]={0};

matrix_t h[8]={0};

matrix_t * eigen_calculations(matrix_t *inp_a)
{
	float c=0;
	float sin_theta = 0;
	float cos_theta = 0;
	float check = -(*(inp_a+2)) /  (*(inp_a+0));

	if(check < 0)
	{
		 c = -1;
	}
	else
	{
		 c = 1;
	}

	sin_theta=c * (*(inp_a+2)) / x_sqrt((pow(*(inp_a+0),2)+pow(*(inp_a+2),2)));
	cos_theta=    *(inp_a+0)   / x_sqrt((pow(*(inp_a+0),2)+pow(*(inp_a+2),2)));

	h[0] =cos_theta;
	h[1] =sin_theta;
	h[2] =-1 * (sin_theta);
	h[3] =cos_theta;
	h[4] =( *(inp_a+0)) * cos_theta - ( *(inp_a+2)) * sin_theta;
	h[5] =( *(inp_a+1)) * cos_theta - ( *(inp_a+3)) * sin_theta;
	h[6] =( *(inp_a+0)) * sin_theta + ( *(inp_a+2)) * cos_theta;
	h[7] =( *(inp_a+1)) * sin_theta + ( *(inp_a+3)) * cos_theta;

	return h;
}


float * dummy_function(matrix_t *A)
{
	float diff1=0;
	float diff2=0;
    uint8_t count = 0;

	do
	{
		count++;

		matrix_t A_cal_temp[matrix_size][matrix_size];
		matrix_t Q[matrix_size][matrix_size];
		matrix_t R[matrix_size][matrix_size];


		matrix_t *array_ptr;

		array_ptr = eigen_calculations(A);

		Q[0][0]=array_ptr[0];
		Q[0][1]=array_ptr[1];
		Q[1][0]=array_ptr[2];
		Q[1][1]=array_ptr[3];
		R[0][0]=array_ptr[4];
		R[0][1]=array_ptr[5];
		R[1][0]=array_ptr[6];
		R[1][1]=array_ptr[7];

		hls::matrix_multiply<hls::NoTranspose,hls::NoTranspose,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_t,matrix_t>(Q, R,A_cal);
		hls::matrix_multiply<hls::NoTranspose,hls::NoTranspose,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_t,matrix_t>(R, Q,A_cal_temp);

		diff1 =A_cal_temp[0][0]-A_cal[0][0];
		diff2 =A_cal_temp[1][1]-A_cal[1][1];

		for(int i=0;i < matrix_size;i++)
		{
			for(int j=0;j < matrix_size;j++)
			{
				A_cal[i][j]=A_cal_temp[i][j];
			}
		}

        A[0]= A_cal_temp[0][0];
		A[1]= A_cal_temp[0][1];
		A[2]= A_cal_temp[1][0];
		A[3]= A_cal_temp[1][1];

	}while(diff1 > 0.05);
	printf("\ncount%d \n",count);
	printf("\nNow difference between first diagonal element is: %1.3f\n",diff1);
	printf("Now difference between second diagonal element is: %1.3f\n",diff2);
	printf("First  Eigen value(Lamda_1)  is: %1.3f\n",A[0]);
	printf("second Eigen value(Lamda_2)  is: %1.3f\n",A[3]);
	printf("We got our desired Eigen values(Lamda) After %d iterations",count);
	return A;
}

int decomposition(matrix_t *A,matrix_t *A_recon)
{
#pragma HLS INTERFACE m_axi bundle=Master_bus port=A_recon offset=slave
#pragma HLS INTERFACE m_axi bundle=Master_bus port=A offset=slave
	A_recon = dummy_function(A);
	printf("\nMAIN \n");
	for(int j=0;j<4;j++)
	{
		printf("%1.4f \n",A_recon[j]);
	}
	return 0;
}
