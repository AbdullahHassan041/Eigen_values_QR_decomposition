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
float e4[row];
float e5[row];
float e6[row];
float e7[row];
float e8[row];
float e9[row];
float e10[row];
float r[55]={0};
///////////////////////////////////functions/////////////////////////////////////////////////////
float * caculate_ee1(float *ptr_ee1)
{
	float temp1_pow=0;
	float det_e1 =0;
	for(int i = 0; i <column; i++)
		{
		temp1_pow=temp1_pow+(pow(ptr_ee1[i],2));
		}
	det_e1 =x_sqrt(temp1_pow);
	for(int i = 0; i <column; i++)
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
	for(int i = 0; i <column; i++)
		{
		temp2_multiply1=temp2_multiply1+((*(ptr_ee2+i))*ptr_ee1[i]);
		temp2_multiply2[i]=((*(ptr_ee2+i))-(temp2_multiply1* ptr_ee1[i]));
		temp2_pow=temp2_pow+(pow(temp2_multiply2[i],2));
		}
	det_e2 =x_sqrt(temp2_pow);
	for(int i = 0; i <column; i++)
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
	for(int i = 0; i <column; i++)
	{
		temp3_multiply1=temp3_multiply1+((*(ptr_ee3+i))*ptr_ee1[i]);
		temp3_multiply2=temp3_multiply2+((*(ptr_ee3+i))*ptr_ee2[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp3_multiply3[i]=((ptr_ee1[i]) * temp3_multiply1);
		temp3_multiply4[i]=((ptr_ee2[i]) * temp3_multiply2);
		temp3_multiply5[i]=( (*(ptr_ee3+i)) - (temp3_multiply3[i]) - (temp3_multiply4[i]) );
		temp3_pow=temp3_pow+(pow(temp3_multiply5[i],2));
	}
	det_e3 =x_sqrt(temp3_pow);
	for(int i = 0; i <column; i++)
	{
		e3[i]=(temp3_multiply5[i])/ det_e3;
	}
	return e3;
}

float * caculate_ee4(float *ptr_ee4,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3)
{
	float temp4_pow=0;
	float det_e4 =0;
	float temp4_multiply1=0;
	float temp4_multiply2=0;
	float temp4_multiply3=0;
	float temp4_multiply4[row];
	float temp4_multiply5[row];
	float temp4_multiply6[row];
	float temp4_multiply7[row];
	for(int i = 0; i <column; i++)
	{
		temp4_multiply1=temp4_multiply1+((*(ptr_ee4+i))*ptr_ee1[i]);
		temp4_multiply2=temp4_multiply2+((*(ptr_ee4+i))*ptr_ee2[i]);
		temp4_multiply3=temp4_multiply3+((*(ptr_ee4+i))*ptr_ee3[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp4_multiply4[i]=((ptr_ee1[i]) * temp4_multiply1);
		temp4_multiply5[i]=((ptr_ee2[i]) * temp4_multiply2);
		temp4_multiply6[i]=((ptr_ee3[i]) * temp4_multiply3);
		temp4_multiply7[i]=( (*(ptr_ee4+i)) - (temp4_multiply4[i]) - (temp4_multiply5[i]) - (temp4_multiply6[i]) );
		temp4_pow=temp4_pow+(pow(temp4_multiply7[i],2));
	}
	det_e4 =x_sqrt(temp4_pow);
	for(int i = 0; i <column; i++)
	{
		e4[i]=(temp4_multiply7[i])/ det_e4;
	}
	return e4;
}

float * caculate_ee5(float *ptr_ee5,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4)
{
	float temp5_pow=0;
	float det_e5 =0;
	float temp5_multiply1=0;
	float temp5_multiply2=0;
	float temp5_multiply3=0;
	float temp5_multiply4=0;
	float temp5_multiply5[row];
	float temp5_multiply6[row];
	float temp5_multiply7[row];
	float temp5_multiply8[row];
	float temp5_multiply9[row];
	for(int i = 0; i <column; i++)
	{
		temp5_multiply1=temp5_multiply1+((*(ptr_ee5+i))*ptr_ee1[i]);
		temp5_multiply2=temp5_multiply2+((*(ptr_ee5+i))*ptr_ee2[i]);
		temp5_multiply3=temp5_multiply3+((*(ptr_ee5+i))*ptr_ee3[i]);
		temp5_multiply4=temp5_multiply3+((*(ptr_ee5+i))*ptr_ee4[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp5_multiply5[i]=((ptr_ee1[i]) * temp5_multiply1);
		temp5_multiply6[i]=((ptr_ee2[i]) * temp5_multiply2);
		temp5_multiply7[i]=((ptr_ee3[i]) * temp5_multiply3);
		temp5_multiply8[i]=((ptr_ee4[i]) * temp5_multiply4);
		temp5_multiply9[i]=( (*(ptr_ee5+i)) - (temp5_multiply5[i]) - (temp5_multiply6[i]) - (temp5_multiply7[i]) - (temp5_multiply8[i]) );
		temp5_pow=temp5_pow+(pow(temp5_multiply9[i],2));
	}
	det_e5 =x_sqrt(temp5_pow);
	for(int i = 0; i <column; i++)
	{
		e5[i]=(temp5_multiply9[i])/ det_e5;
	}
	return e5;
}

float * caculate_ee6(float *ptr_ee6,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4,float *ptr_ee5)
{
	float temp6_pow=0;
	float det_e6 =0;
	float temp6_multiply1=0;
	float temp6_multiply2=0;
	float temp6_multiply3=0;
	float temp6_multiply4=0;
	float temp6_multiply5=0;
	float temp6_multiply6[row];
	float temp6_multiply7[row];
	float temp6_multiply8[row];
	float temp6_multiply9[row];
	float temp6_multiply10[row];
	float temp6_multiply11[row];
	for(int i = 0; i <column; i++)
	{
		temp6_multiply1=temp6_multiply1+((*(ptr_ee6+i))*ptr_ee1[i]);
		temp6_multiply2=temp6_multiply2+((*(ptr_ee6+i))*ptr_ee2[i]);
		temp6_multiply3=temp6_multiply3+((*(ptr_ee6+i))*ptr_ee3[i]);
		temp6_multiply4=temp6_multiply4+((*(ptr_ee6+i))*ptr_ee4[i]);
		temp6_multiply5=temp6_multiply5+((*(ptr_ee6+i))*ptr_ee5[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp6_multiply6[i]=((ptr_ee1[i]) * temp6_multiply1);
		temp6_multiply7[i]=((ptr_ee2[i]) * temp6_multiply2);
		temp6_multiply8[i]=((ptr_ee3[i]) * temp6_multiply3);
		temp6_multiply9[i]=((ptr_ee4[i]) * temp6_multiply4);
		temp6_multiply10[i]=((ptr_ee5[i]) * temp6_multiply5);
		temp6_multiply11[i]=( (*(ptr_ee6+i)) - (temp6_multiply6[i]) - (temp6_multiply7[i]) - (temp6_multiply8[i]) - (temp6_multiply9[i])- (temp6_multiply10[i]) );
		temp6_pow=temp6_pow+(pow(temp6_multiply11[i],2));
	}
	det_e6 =x_sqrt(temp6_pow);
	for(int i = 0; i <column; i++)
	{
		e6[i]=(temp6_multiply11[i])/ det_e6;
	}
	return e6;
}

float * caculate_ee7(float *ptr_ee7,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4,float *ptr_ee5,float *ptr_ee6)
{
	float temp7_pow=0;
	float det_e7 =0;
	float temp7_multiply1=0;
	float temp7_multiply2=0;
	float temp7_multiply3=0;
	float temp7_multiply4=0;
	float temp7_multiply5=0;
	float temp7_multiply6=0;
	float temp7_multiply7[row];
	float temp7_multiply8[row];
	float temp7_multiply9[row];
	float temp7_multiply10[row];
	float temp7_multiply11[row];
	float temp7_multiply12[row];
	float temp7_multiply13[row];
	for(int i = 0; i <column; i++)
	{
		temp7_multiply1=temp7_multiply1+((*(ptr_ee7+i))*ptr_ee1[i]);
		temp7_multiply2=temp7_multiply2+((*(ptr_ee7+i))*ptr_ee2[i]);
		temp7_multiply3=temp7_multiply3+((*(ptr_ee7+i))*ptr_ee3[i]);
		temp7_multiply4=temp7_multiply4+((*(ptr_ee7+i))*ptr_ee4[i]);
		temp7_multiply5=temp7_multiply5+((*(ptr_ee7+i))*ptr_ee5[i]);
		temp7_multiply6=temp7_multiply6+((*(ptr_ee7+i))*ptr_ee6[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp7_multiply7[i]=((ptr_ee1[i]) * temp7_multiply1);
		temp7_multiply8[i]=((ptr_ee2[i]) * temp7_multiply2);
		temp7_multiply9[i]=((ptr_ee3[i]) * temp7_multiply3);
		temp7_multiply10[i]=((ptr_ee4[i]) * temp7_multiply4);
		temp7_multiply11[i]=((ptr_ee5[i]) * temp7_multiply5);
		temp7_multiply12[i]=((ptr_ee6[i]) * temp7_multiply6);
		temp7_multiply13[i]=( (*(ptr_ee7+i)) - (temp7_multiply7[i]) - (temp7_multiply8[i]) - (temp7_multiply9[i]) - (temp7_multiply10[i]) - (temp7_multiply11[i]) - (temp7_multiply12[i]) );
		temp7_pow=temp7_pow+(pow(temp7_multiply13[i],2));
	}
	det_e7 =x_sqrt(temp7_pow);
	for(int i = 0; i <column; i++)
	{
		e7[i]=(temp7_multiply13[i])/ det_e7;
	}
	return e7;
}


float * caculate_ee8(float *ptr_ee8,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4,float *ptr_ee5,float *ptr_ee6,float *ptr_ee7)
{
	float temp8_pow=0;
	float det_e8 =0;
	float temp8_multiply1=0;
	float temp8_multiply2=0;
	float temp8_multiply3=0;
	float temp8_multiply4=0;
	float temp8_multiply5=0;
	float temp8_multiply6=0;
	float temp8_multiply7=0;
	float temp8_multiply8[row];
	float temp8_multiply9[row];
	float temp8_multiply10[row];
	float temp8_multiply11[row];
	float temp8_multiply12[row];
	float temp8_multiply13[row];
	float temp8_multiply14[row];
	float temp8_multiply15[row];
	for(int i = 0; i <column; i++)
	{
		temp8_multiply1=temp8_multiply1+((*(ptr_ee8+i))*ptr_ee1[i]);
		temp8_multiply2=temp8_multiply2+((*(ptr_ee8+i))*ptr_ee2[i]);
		temp8_multiply3=temp8_multiply3+((*(ptr_ee8+i))*ptr_ee3[i]);
		temp8_multiply4=temp8_multiply4+((*(ptr_ee8+i))*ptr_ee4[i]);
		temp8_multiply5=temp8_multiply5+((*(ptr_ee8+i))*ptr_ee5[i]);
		temp8_multiply6=temp8_multiply6+((*(ptr_ee8+i))*ptr_ee6[i]);
		temp8_multiply7=temp8_multiply7+((*(ptr_ee8+i))*ptr_ee7[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp8_multiply8[i]=((ptr_ee1[i]) * temp8_multiply1);
		temp8_multiply9[i]=((ptr_ee2[i]) * temp8_multiply2);
		temp8_multiply10[i]=((ptr_ee3[i]) * temp8_multiply3);
		temp8_multiply11[i]=((ptr_ee4[i]) * temp8_multiply4);
		temp8_multiply12[i]=((ptr_ee5[i]) * temp8_multiply5);
		temp8_multiply13[i]=((ptr_ee6[i]) * temp8_multiply6);
		temp8_multiply14[i]=((ptr_ee7[i]) * temp8_multiply7);
		temp8_multiply15[i]=( (*(ptr_ee8+i)) - (temp8_multiply8[i]) - (temp8_multiply9[i]) - (temp8_multiply10[i]) - (temp8_multiply11[i]) - (temp8_multiply12[i]) - (temp8_multiply13[i])- (temp8_multiply14[i]) );
		temp8_pow=temp8_pow+(pow(temp8_multiply15[i],2));
	}
	det_e8 =x_sqrt(temp8_pow);
	for(int i = 0; i <column; i++)
	{
		e8[i]=(temp8_multiply15[i])/ det_e8;
	}
	return e8;
}

float * caculate_ee9(float *ptr_ee9,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4,float *ptr_ee5,float *ptr_ee6,float *ptr_ee7,float *ptr_ee8)
{
	float temp9_pow=0;
	float det_e9 =0;
	float temp9_multiply1=0;
	float temp9_multiply2=0;
	float temp9_multiply3=0;
	float temp9_multiply4=0;
	float temp9_multiply5=0;
	float temp9_multiply6=0;
	float temp9_multiply7=0;
	float temp9_multiply8=0;
	float temp9_multiply9[row];
	float temp9_multiply10[row];
	float temp9_multiply11[row];
	float temp9_multiply12[row];
	float temp9_multiply13[row];
	float temp9_multiply14[row];
	float temp9_multiply15[row];
	float temp9_multiply16[row];
	float temp9_multiply17[row];
	for(int i = 0; i <column; i++)
	{
		temp9_multiply1=temp9_multiply1+((*(ptr_ee9+i))*ptr_ee1[i]);
		temp9_multiply2=temp9_multiply2+((*(ptr_ee9+i))*ptr_ee2[i]);
		temp9_multiply3=temp9_multiply3+((*(ptr_ee9+i))*ptr_ee3[i]);
		temp9_multiply4=temp9_multiply4+((*(ptr_ee9+i))*ptr_ee4[i]);
		temp9_multiply5=temp9_multiply5+((*(ptr_ee9+i))*ptr_ee5[i]);
		temp9_multiply6=temp9_multiply6+((*(ptr_ee9+i))*ptr_ee6[i]);
		temp9_multiply7=temp9_multiply7+((*(ptr_ee9+i))*ptr_ee7[i]);
		temp9_multiply8=temp9_multiply8+((*(ptr_ee9+i))*ptr_ee8[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp9_multiply9[i]=((ptr_ee1[i]) * temp9_multiply1);
		temp9_multiply10[i]=((ptr_ee2[i]) * temp9_multiply2);
		temp9_multiply11[i]=((ptr_ee3[i]) * temp9_multiply3);
		temp9_multiply12[i]=((ptr_ee4[i]) * temp9_multiply4);
		temp9_multiply13[i]=((ptr_ee5[i]) * temp9_multiply5);
		temp9_multiply14[i]=((ptr_ee6[i]) * temp9_multiply6);
		temp9_multiply15[i]=((ptr_ee7[i]) * temp9_multiply7);
		temp9_multiply16[i]=((ptr_ee8[i]) * temp9_multiply8);

		temp9_multiply17[i]=( (*(ptr_ee9+i)) - (temp9_multiply9[i]) - (temp9_multiply10[i]) - (temp9_multiply11[i]) - (temp9_multiply12[i]) - (temp9_multiply13[i]) - (temp9_multiply14[i]) - (temp9_multiply15[i])  - (temp9_multiply16[i]) );
		temp9_pow=temp9_pow+(pow(temp9_multiply17[i],2));
	}
	det_e9 =x_sqrt(temp9_pow);
	for(int i = 0; i <column; i++)
	{
		e9[i]=(temp9_multiply17[i])/ det_e9;
	}
	return e9;
}


float * caculate_ee10(float *ptr_ee10,float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4,float *ptr_ee5,float *ptr_ee6,float *ptr_ee7,float *ptr_ee8,float *ptr_ee9)
{
	float temp10_pow=0;
	float det_e10 =0;
	float temp10_multiply1=0;
	float temp10_multiply2=0;
	float temp10_multiply3=0;
	float temp10_multiply4=0;
	float temp10_multiply5=0;
	float temp10_multiply6=0;
	float temp10_multiply7=0;
	float temp10_multiply8=0;
	float temp10_multiply9=0;
	float temp10_multiply10[row];
	float temp10_multiply11[row];
	float temp10_multiply12[row];
	float temp10_multiply13[row];
	float temp10_multiply14[row];
	float temp10_multiply15[row];
	float temp10_multiply16[row];
	float temp10_multiply17[row];
	float temp10_multiply18[row];
	float temp10_multiply19[row];
	for(int i = 0; i <column; i++)
	{
		temp10_multiply1=temp10_multiply1+((*(ptr_ee10+i))*ptr_ee1[i]);
		temp10_multiply2=temp10_multiply2+((*(ptr_ee10+i))*ptr_ee2[i]);
		temp10_multiply3=temp10_multiply3+((*(ptr_ee10+i))*ptr_ee3[i]);
		temp10_multiply4=temp10_multiply4+((*(ptr_ee10+i))*ptr_ee4[i]);
		temp10_multiply5=temp10_multiply5+((*(ptr_ee10+i))*ptr_ee5[i]);
		temp10_multiply6=temp10_multiply6+((*(ptr_ee10+i))*ptr_ee6[i]);
		temp10_multiply7=temp10_multiply7+((*(ptr_ee10+i))*ptr_ee7[i]);
		temp10_multiply8=temp10_multiply8+((*(ptr_ee10+i))*ptr_ee8[i]);
		temp10_multiply9=temp10_multiply9+((*(ptr_ee10+i))*ptr_ee9[i]);
	}
	for(int i = 0; i <column; i++)
	{
		temp10_multiply10[i]=((ptr_ee1[i]) * temp10_multiply1);
		temp10_multiply11[i]=((ptr_ee2[i]) * temp10_multiply2);
		temp10_multiply12[i]=((ptr_ee3[i]) * temp10_multiply3);
		temp10_multiply13[i]=((ptr_ee4[i]) * temp10_multiply4);
		temp10_multiply14[i]=((ptr_ee5[i]) * temp10_multiply5);
		temp10_multiply15[i]=((ptr_ee6[i]) * temp10_multiply6);
		temp10_multiply16[i]=((ptr_ee7[i]) * temp10_multiply7);
		temp10_multiply17[i]=((ptr_ee8[i]) * temp10_multiply8);
		temp10_multiply18[i]=((ptr_ee9[i]) * temp10_multiply9);

		temp10_multiply19[i]=( (*(ptr_ee10+i)) - (temp10_multiply10[i]) - (temp10_multiply11[i]) - (temp10_multiply12[i]) - (temp10_multiply13[i]) - (temp10_multiply14[i]) - (temp10_multiply15[i]) - (temp10_multiply16[i])  - (temp10_multiply17[i])  - (temp10_multiply18[i]));
		temp10_pow=temp10_pow+(pow(temp10_multiply19[i],2));
	}
	det_e10 =x_sqrt(temp10_pow);
	for(int i = 0; i <column; i++)
	{
		e10[i]=(temp10_multiply19[i])/ det_e10;
	}
	return e10;
}


float * upper_triagular_matrix(float *ptr_col1,float *ptr_col2,float *ptr_col3,float *ptr_col4,float *ptr_col5,float *ptr_col6,float *ptr_col7,float *ptr_col8,
		float *ptr_col9,float *ptr_col10,
		float *ptr_ee1,float *ptr_ee2,float *ptr_ee3,float *ptr_ee4,float *ptr_ee5,float *ptr_ee6,float *ptr_ee7,float *ptr_ee8,float *ptr_ee9,float *ptr_ee10)
{
	float temp_r_00 = 0;
	float temp_r_01 = 0;
	float temp_r_02 = 0;
	float temp_r_11 = 0;
	float temp_r_12 =0;
	float temp_r_22 =0;
	///////////////////////////////////////row1//////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[0]=r[0]+((*(ptr_col1+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[1]=r[1]+((*(ptr_col2+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[2]=r[2]+((*(ptr_col3+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[3]=r[3]+((*(ptr_col4+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[4]=r[4]+((*(ptr_col5+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[5]=r[5]+((*(ptr_col6+i))*ptr_ee1[i]);
	}
	//
	for(int i=0;i<row;i++)
	{
		r[6]=r[6]+((*(ptr_col7+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[7]=r[7]+((*(ptr_col8+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[8]=r[8]+((*(ptr_col9+i))*ptr_ee1[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[9]=r[9]+((*(ptr_col10+i))*ptr_ee1[i]);
	}
	/////////////////////////////////////////row2/////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[11]=r[11]+((*(ptr_col2+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[12]=r[12]+((*(ptr_col3+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[13]=r[13]+((*(ptr_col4+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[14]=r[14]+((*(ptr_col5+i))*ptr_ee2[i]);
	}
	//
	for(int i=0;i<row;i++)
	{
		r[15]=r[15]+((*(ptr_col6+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[16]=r[16]+((*(ptr_col7+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[17]=r[17]+((*(ptr_col8+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[18]=r[18]+((*(ptr_col9+i))*ptr_ee2[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[19]=r[19]+((*(ptr_col10+i))*ptr_ee2[i]);
	}
	///////////////////////////////////////row3//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[22]=r[22]+((*(ptr_col3+i))*ptr_ee3[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[23]=r[23]+((*(ptr_col4+i))*ptr_ee3[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[24]=r[24]+((*(ptr_col5+i))*ptr_ee3[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[25]=r[25]+((*(ptr_col6+i))*ptr_ee3[i]);
	}
	//
	for(int i=0;i<row;i++)
	{
		r[26]=r[26]+((*(ptr_col7+i))*ptr_ee3[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[27]=r[27]+((*(ptr_col8+i))*ptr_ee3[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[28]=r[28]+((*(ptr_col9+i))*ptr_ee3[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[29]=r[29]+((*(ptr_col10+i))*ptr_ee3[i]);
	}

	///////////////////////////////////////row4//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[33]=r[33]+((*(ptr_col4+i))*ptr_ee4[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[34]=r[34]+((*(ptr_col5+i))*ptr_ee4[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[35]=r[35]+((*(ptr_col6+i))*ptr_ee4[i]);
	}
	//
	for(int i=0;i<row;i++)
	{
		r[36]=r[36]+((*(ptr_col7+i))*ptr_ee4[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[37]=r[37]+((*(ptr_col8+i))*ptr_ee4[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[38]=r[38]+((*(ptr_col9+i))*ptr_ee4[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[39]=r[39]+((*(ptr_col10+i))*ptr_ee4[i]);
	}
	///////////////////////////////////////row5//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[44]=r[44]+((*(ptr_col5+i))*ptr_ee5[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[45]=r[45]+((*(ptr_col6+i))*ptr_ee5[i]);
	}
	//
	for(int i=0;i<row;i++)
	{
		r[46]=r[46]+((*(ptr_col7+i))*ptr_ee5[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[47]=r[47]+((*(ptr_col8+i))*ptr_ee5[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[48]=r[48]+((*(ptr_col9+i))*ptr_ee5[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[49]=r[49]+((*(ptr_col10+i))*ptr_ee5[i]);
	}
	///////////////////////////////////////row6//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[55]=r[55]+((*(ptr_col6+i))*ptr_ee6[i]);
	}
	//
	for(int i=0;i<row;i++)
	{
		r[56]=r[56]+((*(ptr_col7+i))*ptr_ee6[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[57]=r[57]+((*(ptr_col8+i))*ptr_ee6[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[58]=r[58]+((*(ptr_col9+i))*ptr_ee6[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[59]=r[59]+((*(ptr_col10+i))*ptr_ee6[i]);
	}
	///////////////////////////////////////row7//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[66]=r[66]+((*(ptr_col7+i))*ptr_ee7[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[67]=r[67]+((*(ptr_col8+i))*ptr_ee7[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[68]=r[68]+((*(ptr_col9+i))*ptr_ee7[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[69]=r[69]+((*(ptr_col10+i))*ptr_ee7[i]);
	}
	///////////////////////////////////////row8//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[77]=r[77]+((*(ptr_col8+i))*ptr_ee8[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[78]=r[78]+((*(ptr_col9+i))*ptr_ee8[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[79]=r[79]+((*(ptr_col10+i))*ptr_ee8[i]);
	}
	///////////////////////////////////////row9//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[88]=r[88]+((*(ptr_col9+i))*ptr_ee9[i]);
	}
	for(int i=0;i<row;i++)
	{
		r[89]=r[89]+((*(ptr_col10+i))*ptr_ee9[i]);
	}
	///////////////////////////////////////row10//////////////////////////////////////////////////////////////
	for(int i=0;i<row;i++)
	{
		r[99]=r[99]+((*(ptr_col10+i))*ptr_ee10[i]);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return r;
}
////////////////////////////////////////////////////////////////////////////////////


void Status_check_call(matrix_t input[row][column])
{
	float col1[column];
	float col2[column];
	float col3[column];
	float col4[column];
	float col5[column];
	float col6[column];
	float col7[column];
	float col8[column];
	float col9[column];
	float col10[column];

	matrix_t Q[row][column];
	matrix_t R[row][column];
	matrix_t QR[row][column];
	matrix_t A_temp[row][column];
	float diff1=0;
	float diff2=0;
	float iterations=0;
	float first_check=1;
	do
	{
		iterations++;
		/////////////////////////for filling columns///////////////////////////////////////////
		for(int i = 0; i <column; i++)
		{
			col1[i]=input[i][0];
		}

		for(int i = 0; i <column; i++)
		{
			col2[i]=input[i][1];
		}
		for(int i = 0; i <column; i++)
		{
			col3[i]=input[i][2];
		}
		for(int i = 0; i <column; i++)
		{
			col4[i]=input[i][3];
		}

		for(int i = 0; i <column; i++)
		{
			col5[i]=input[i][4];
		}
		for(int i = 0; i <column; i++)
		{
			col6[i]=input[i][5];
		}
		for(int i = 0; i <column; i++)
		{
			col7[i]=input[i][6];
		}

		for(int i = 0; i <column; i++)
		{
			col8[i]=input[i][7];
		}
		for(int i = 0; i <column; i++)
		{
			col9[i]=input[i][8];
		}
		for(int i = 0; i <column; i++)
		{
			col10[i]=input[i][9];
		}
		/////////////////////////////////////////////////////////////////////////////////////////////
		////////////////for col1 calculations/////////////////////////
		float *ptr_e1;
		ptr_e1=caculate_ee1(col1);

		////////////////for col2 calculations/////////////////////////
		float *ptr_e2;
		ptr_e2=caculate_ee2(col2,ptr_e1);

		////////////////for col3 calculations/////////////////////////
		float *ptr_e3;
		ptr_e3=caculate_ee3(col3,ptr_e1,ptr_e2);

		////////////////for col4 calculations/////////////////////////
		float *ptr_e4;
		ptr_e4=caculate_ee4(col4,ptr_e1,ptr_e2,ptr_e3);

		////////////////for col5 calculations/////////////////////////
		float *ptr_e5;
		ptr_e5=caculate_ee5(col5,ptr_e1,ptr_e2,ptr_e3,ptr_e4);

		////////////////for col6 calculations/////////////////////////
		float *ptr_e6;
		ptr_e6=caculate_ee6(col6,ptr_e1,ptr_e2,ptr_e3,ptr_e4,ptr_e5);

		////////////////for col7 calculations/////////////////////////
		float *ptr_e7;
		ptr_e7=caculate_ee7(col7,ptr_e1,ptr_e2,ptr_e3,ptr_e4,ptr_e5,ptr_e6);

		////////////////for col8 calculations/////////////////////////
		float *ptr_e8;
		ptr_e8=caculate_ee8(col8,ptr_e1,ptr_e2,ptr_e3,ptr_e4,ptr_e5,ptr_e6,ptr_e7);

		////////////////for col9 calculations/////////////////////////
		float *ptr_e9;
		ptr_e9=caculate_ee9(col9,ptr_e1,ptr_e2,ptr_e3,ptr_e4,ptr_e5,ptr_e6,ptr_e7,ptr_e8);

		////////////////for col10 calculations/////////////////////////
		float *ptr_e10;
		ptr_e10=caculate_ee10(col10,ptr_e1,ptr_e2,ptr_e3,ptr_e4,ptr_e5,ptr_e6,ptr_e7,ptr_e8,ptr_e9);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////

		///////////for r calculation//////////////////////////////////////////////////////////////////////////////
		float *ptr_r;
		ptr_r = upper_triagular_matrix(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,ptr_e1,ptr_e2,ptr_e3,ptr_e4,ptr_e5,ptr_e6,ptr_e7,ptr_e8,ptr_e9,ptr_e10);
		////////////////////////////////////Filling r1////////////////////////////////////////////////////////////
		for(int i=0;i<row;i++)
		{
			R[0][i]=*(ptr_r+i);
		}
		////////////////////////////////////Filling r2////////////////////////////////////////////////////////////
		 int a2=10;
		for(int i=0;i<row;i++)
		{
			if(i<1)
			{
				R[1][i]=0;
			}
			else
			{
				R[1][i]=*(ptr_r+a2);
			}
		   a2++;
		}
		////////////////////////////////////Filling r3////////////////////////////////////////////////////////////
		int a3=20;
		for(int i=0;i<row;i++)
		{
			if(i<2)
			{
				R[2][i]=0;
			}
			else
			{
				R[2][i]=*(ptr_r+a3);
			}
		   a3++;
		}
		////////////////////////////////////Filling r4////////////////////////////////////////////////////////////
		int a4=30;
		for(int i=0;i<row;i++)
		{
			if(i<3)
			{
				R[3][i]=0;
			}
			else
			{
				R[3][i]=*(ptr_r+a4);
			}
		   a4++;
		}
		////////////////////////////////////Filling r5////////////////////////////////////////////////////////////
		int a5=40;
		for(int i=0;i<row;i++)
		{
			if(i<4)
			{
				R[4][i]=0;
			}
			else
			{
				R[4][i]=*(ptr_r+a5);
			}
		   a5++;
		}
		////////////////////////////////////Filling r6////////////////////////////////////////////////////////////
		int a6=50;
		for(int i=0;i<row;i++)
		{
			if(i<5)
			{
				R[5][i]=0;
			}
			else
			{
				R[5][i]=*(ptr_r+a6);
			}
		   a6++;
		}
		////////////////////////////////////Filling r7////////////////////////////////////////////////////////////
		int a7=60;
		for(int i=0;i<row;i++)
		{
			if(i<6)
			{
				R[6][i]=0;
			}
			else
			{
				R[6][i]=*(ptr_r+a7);
			}
		   a7++;
		}
		////////////////////////////////////Filling r8////////////////////////////////////////////////////////////
		int a8=70;
		for(int i=0;i<row;i++)
		{
			if(i<7)
			{
				R[7][i]=0;
			}
			else
			{
				R[7][i]=*(ptr_r+a8);
			}
		   a8++;
		}
		////////////////////////////////////Filling r9////////////////////////////////////////////////////////////
		int a9=80;
		for(int i=0;i<row;i++)
		{
			if(i<8)
			{
				R[8][i]=0;
			}
			else
			{
				R[8][i]=*(ptr_r+a9);
			}
		   a9++;
		}
		////////////////////////////////////Filling r10////////////////////////////////////////////////////////////
		int a10=90;
		for(int i=0;i<row;i++)
		{
			if(i<9)
			{
				R[9][i]=0;
			}
			else
			{
				R[9][i]=*(ptr_r+a10);
			}
		   a10++;
		}
		//////////////////////////////////////////end of Filling///////////////////////////////////////////////////////////////////

		/////////////////////////////////filling Q/////////////////

		for(int i=0;i < row;i++)
		{
			for(int j=0;j < 1;j++)
			{
				Q[i][0]=ptr_e1[i];
				Q[i][1]=ptr_e2[i];
				Q[i][2]=ptr_e3[i];
				Q[i][3]=ptr_e4[i];
				Q[i][4]=ptr_e5[i];
				Q[i][5]=ptr_e6[i];
				Q[i][6]=ptr_e7[i];
				Q[i][7]=ptr_e8[i];
				Q[i][8]=ptr_e9[i];
				Q[i][9]=ptr_e10[i];
			}
		}
	////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////for retrieveing/////////////////////////////////////
		if(first_check==1)
		{
			hls::matrix_multiply<hls::NoTranspose,hls::NoTranspose,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_t,matrix_t>(Q, R, QR);
			++first_check;
		}
		else
		{
			for(int i=0;i < row;i++)
			{
				for(int j=0;j < column;j++)
				{
					QR[i][j]=A_temp[i][j];
				}
			}
		}
		hls::matrix_multiply<hls::NoTranspose,hls::NoTranspose,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_size,matrix_t,matrix_t>(R, Q, A_temp);
		diff1 =A_temp[0][0]-QR[0][0];
		diff2 =A_temp[1][1]-QR[1][1];
		printf("           After Decomposition \n");
		printf("A_temp = \n");
		for(int j = 0; j <row; j++)
			{
				for(int i=0; i<column;i++)
				{
					printf("%1.4f  ",A_temp[i][j]);
				}
				printf("\n");
			}
		for(int i=0;i < row;i++)
		{
			for(int j=0;j < column;j++)
			{
				input[i][j]=A_temp[i][j];
			}
		}

	}while((diff1 > 0.05) || (diff2 > 0.05) );
	printf("           OUT FINAL MATRIX IS \n");
	printf("QR=\n");
	for(int j = 0; j <row; j++)
		{
			for(int i=0; i<column;i++)
			{
				printf("%1.4f  ",QR[i][j]);
			}
			printf("\n");
		}
	for(int i=0;i < row;i++)
	{
		printf("%d Eigen value(Lamda)  is: %1.3f\n",i,QR[i][i]);
	}
	printf("Total iterations are :%1.4f  \n",iterations);
	printf("first diif is:%1.4f  \n",diff1);
	printf("second diif is:%1.4f  \n",diff2);
}
////////////////////////////////////////main_functions//////////////////////////////////////////
void decompose(matrix_t input[row][column])
{
	Status_check_call(input);
}
///////////////////////////////////////////////////////////////////////////////////////
