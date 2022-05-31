#include "dec.h"
#include <iostream>
#include "stdio.h"
#include "stdint.h"
#include "hls/linear_algebra/utils/x_hls_matrix_utils.h"
#include "hls_stream.h"
using namespace hls;

float e1[row];
float e2[row];
float e3[row];
float r[9]={0};
///////////////////////////////////functions/////////////////////////////////////////////////////
float * caculate_ee1(float *ptr_ee1)
{
	float temp1_pow=0;
	float det_e1 =0;
	for(int i = 0; i <3; i++)
		{
		temp1_pow=temp1_pow+(pow(ptr_ee1[i],2));
		}
	det_e1 =x_sqrt(temp1_pow);
	for(int i = 0; i <3; i++)
		{
		e1[i]=(*(ptr_ee1+i))/ det_e1;
		}
	return e1;
}

float * caculate_ee2(float *ptr_ee2,float *ptr_ee1)
{
	float temp2_pow=0;
	float temp2_multiply1=0;
	float temp2_multiply2[row];
	float det_e2 =0;
	for(int i = 0; i <3; i++)
		{
		temp2_multiply1=temp2_multiply1+((*(ptr_ee2+i))*ptr_ee1[i]);
		temp2_multiply2[i]=((*(ptr_ee2+i))-(temp2_multiply1* ptr_ee1[i]));
		temp2_pow=temp2_pow+(pow(temp2_multiply2[i],2));
		}
	det_e2 =x_sqrt(temp2_pow);
	for(int i = 0; i <3; i++)
		{
		e2[i]=(temp2_multiply2[i])/ det_e2;
		}
	return e2;
}

float * caculate_ee3(float *ptr_ee3,float *ptr_ee1,float *ptr_ee2)
{
	float temp3_pow=0;
	float det_e3 =0;
	float temp3_multiply1=0;
	float temp3_multiply2=0;
	float temp3_multiply3[row];
	float temp3_multiply4[row];
	float temp3_multiply5[row];
	for(int i = 0; i <3; i++)
	{
		temp3_multiply1=temp3_multiply1+((*(ptr_ee3+i))*ptr_ee1[i]);
		temp3_multiply2=temp3_multiply2+((*(ptr_ee3+i))*ptr_ee2[i]);
	}
	for(int i = 0; i <3; i++)
	{
		temp3_multiply3[i]=((ptr_ee1[i]) * temp3_multiply1);
		temp3_multiply4[i]=((ptr_ee2[i]) * temp3_multiply2);
		temp3_multiply5[i]=( (*(ptr_ee3+i)) - (temp3_multiply3[i]) - (temp3_multiply4[i]) );
		temp3_pow=temp3_pow+(pow(temp3_multiply5[i],2));
	}
	det_e3 =x_sqrt(temp3_pow);
	for(int i = 0; i <3; i++)
	{
		e3[i]=(temp3_multiply5[i])/ det_e3;
	}
	return e3;
}

float * upper_triagular_matrix(float *ptr_col1,float *ptr_col2,float *ptr_col3,
		float *ptr_ee1,float *ptr_ee2,float *ptr_ee3)
{
	float temp_r_00 = 0;
	float temp_r_01 = 0;
	float temp_r_02 = 0;
	float temp_r_11 = 0;
	float temp_r_12 =0;
	float temp_r_22 =0;
	for(int i=0;i<3;i++)
	{
	temp_r_00=temp_r_00+((*(ptr_col1+i))*ptr_ee1[i]);
	}
	for(int i=0;i<3;i++)
	{
	temp_r_01=temp_r_01+((*(ptr_col2+i))*ptr_ee1[i]);
	}
	for(int i=0;i<3;i++)
	{
	temp_r_02=temp_r_02+((*(ptr_col3+i))*ptr_ee1[i]);
	}
	for(int i=0;i<3;i++)
	{
	temp_r_11=temp_r_11+((*(ptr_col2+i))*ptr_ee2[i]);
	}
	for(int i=0;i<3;i++)
	{
	temp_r_12=temp_r_12+((*(ptr_col3+i))*ptr_ee2[i]);
	}
	for(int i=0;i<3;i++)
	{
	temp_r_22=temp_r_22+((*(ptr_col3+i))*ptr_ee3[i]);
	}

	r[0]=temp_r_00;
	r[1]=temp_r_01;
	r[2]=temp_r_02;
	r[4]=temp_r_11;
	r[5]=temp_r_12;
	r[8]=temp_r_22;
	return r;
}
////////////////////////////////////////////////////////////////////////////////////


void Status_check_call(matrix_t input[row][column])
{
	float col1[row];
	float col2[row];
	float col3[row];
	matrix_t Q[row][column];
	matrix_t R[row][column];
	matrix_t QR[row][column]={0};
	matrix_t A_temp[row][column]={0};
	float diff1=0;
	float diff2=0;
	float iterations=0;
	do
	{
		iterations++;
		for(int i = 0; i <3; i++)
		{
			col1[i]=input[i][0];
		}

		for(int i = 0; i <3; i++)
		{
			col2[i]=input[i][1];
		}
		for(int i = 0; i <3; i++)
		{
			col3[i]=input[i][2];
		}
		////////////////for col1 calculations/////////////////////////
		float *ptr_e1;
		ptr_e1=caculate_ee1(col1);
		////////////////for col1 calculations/////////////////////////
		float *ptr_e2;
		ptr_e2=caculate_ee2(col2,ptr_e1);
		////////////////for col1 calculations/////////////////////////
		float *ptr_e3;
		ptr_e3=caculate_ee3(col3,ptr_e1,ptr_e2);
		////////////////////////////////////////////////////////////////

		///////////for r calculation/////////////////////////////////
		float *ptr_r;
		ptr_r = upper_triagular_matrix(col1,col2,col3,ptr_e1,ptr_e2,ptr_e3);

		for(int i=0;i<3;i++)
		{
			R[0][i]=*(ptr_r+i);
		}
		 int k=3;
		for(int i=0;i<3;i++)
		{
			if(i<1)
			{
				R[1][i]=0;
			}
			else
			{
				R[1][i]=*(ptr_r+k);
			}
		   k++;
		}
		int l=6;
		for(int i=0;i<3;i++)
		{
			if(i<2)
			{
				R[2][i]=0;
			}
			else
			{
				R[2][i]=*(ptr_r+l);
			}
		   l++;
		}
		/////////////////////////////////filling Q/////////////////

		Q[0][0]=ptr_e1[0];
		Q[0][1]=ptr_e2[0];
		Q[0][2]=ptr_e3[0];

		Q[1][0]=ptr_e1[1];
		Q[1][1]=ptr_e2[1];
		Q[1][2]=ptr_e3[1];

		Q[2][0]=ptr_e1[2];
		Q[2][1]=ptr_e2[2];
		Q[2][2]=ptr_e3[2];

	///////////////////////////////////////////////////////////////////////////
		printf("           Before Decomposition \n");

		printf("Actual Input\n");
		for(int j = 0; j <3; j++)
			{
				for(int i=0; i<3;i++)
				{
					printf("%1.4f  ",input[i][j]);
				}
				printf("\n");
			}
		printf("           After Decomposition \n");
		printf("Q=\n");

		for(int i = 0; i <3; i++)
			{
				for(int j=0; j<3;j++)
				{
					printf("%1.4f  ",Q[i][j]);
				}
				printf("\n");
			}

		printf("R=\n");

		for(int i = 0; i <3; i++)
			{
				for(int j=0; j<3;j++)
				{
					printf("%1.4f  ",R[i][j]);
				}
				printf("\n");
			}

		/////////////////////////////////////////////////////////////////////////////////////

		/////////////////////////////////////for retrieveing/////////////////////////////////////
		for(int j = 0; j <3; j++)
			{
				for(int i=0; i<3;i++)
				{
					//QR[0]=QR[0] + (Q[0][j] * R[i][0]);
					QR[0][0]=((QR[0][0]) + ( (Q[0][j]) * R[i][0]) );
				}
			}

		printf("           After multiplying Q & R \n");
		printf("QR=\n");
		hls::matrix_multiply<hls::NoTranspose,hls::NoTranspose,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_t,matrix_t>(Q, R, QR);
		hls::matrix_multiply<hls::NoTranspose,hls::NoTranspose,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_t,matrix_t>(R, Q, A_temp);
		diff1 =A_temp[0][0]-QR[0][0];
		diff2 =A_temp[1][1]-QR[1][1];
		printf("first diif is:%1.4f  \n",diff1);
		printf("second diif is:%1.4f  \n",diff2);
		for(int i=0;i < matrix_size;i++)
		{
			for(int j=0;j < matrix_size;j++)
			{
				QR[i][j]=A_temp[i][j];
			}
		}
		for(int i=0;i < matrix_size;i++)
		{
			for(int j=0;j < matrix_size;j++)
			{
				input[i][j]=A_temp[i][j];
			}
		}

	}while(diff1 > 0.05);
	printf("Total iterations are :%1.4f  \n",iterations);
	printf("First  Eigen value(Lamda_1)  is: %1.3f\n",QR[0][0]);
	printf("second Eigen value(Lamda_2)  is: %1.3f\n",QR[1][1]);
	printf("third Eigen value(Lamda_2)  is: %1.3f\n",QR[2][2]);
	for(int j = 0; j <3; j++)
		{
			for(int i=0; i<3;i++)
			{
				printf("%1.4f  ",QR[i][j]);
			}
			printf("\n");
		}
}
////////////////////////////////////////main_functions//////////////////////////////////////////
void decompose(matrix_t input[row][column])
{
	Status_check_call(input);
}
///////////////////////////////////////////////////////////////////////////////////////
