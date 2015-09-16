
//#include <string.h>
#include <math.h>
#include <float.h>
#include "Globals.h"
//#include "ExtGlobals.h"
#include "rsr_utils.h"
#define square_gc(val) val*val

void matmul(int m1rows, int m1cols, double  *M1, \
			int m2rows, int m2cols, double  *M2, \
			double *M3)
{
	int i,j,k;

	for(i=0;i<m1rows;i++)
    {
       for(j=0;j<m2cols;j++)
	   {
			*(M3+i*m2cols+j) = (double)0.0;
		    
			for(k=0;k<m1cols;k++)
	 		{
	           (*(M3+i*m2cols+j))  +=   \
				   (double)(*(M1+i*m1cols+k)) * \
				   (double)(*(M2+k*m2cols+j));
	 		}
	   }
	}
}

void mat_trans( double mk[3][3], double gk[3][3])
{
	int i,j;
	for(i= 0; i <3; i++)
	{
		for (j =0; j<3; j++)
		{
			gk[i][j]=mk[j][i];
		}
	}
}



void swap(double *m1, double *m2)
{

	*m1 = *m1 + *m2;
	*m2 = *m1 - *m2;
	*m1 = *m1 - *m2;
}


void transpose(int rows, int cols, double *M1, double *M2)
{
	int i,j;
	for(i = 0; i< rows; i++)
	{
		for(j =0; j< cols; j++)
		{
			if(M1 != M2)
				(*(M2+j*rows+i)) = (*(M1+i*cols+j));
			else
			{
				// inplace transpose
				if( j > i )
					swap((M1+i*cols+j), (M2+j*rows+i));
			}

		}
	}
}

void matadd(int rows, int cols,double *M1,double *M2,double *M3)
{
	int i,j;

	for(i=0; i<rows; i++)
	{
		for(j=0;j<cols; j++)
		{
			(*(M3+i*cols+j)) = (*(M1+i*cols+j)) + (*(M2+i*cols+j));
		}
	}
}

void matsub(int rows, int cols,double *M1,double *M2,double *M3)
{
	int i,j;

	for(i=0; i<rows; i++)
	{
		for(j=0;j<cols; j++)
		{
			(*(M3+i*cols+j)) = (*(M1+i*cols+j)) - (*(M2+i*cols+j));
		}
	}
}

/*
void matclone(int rows,int cols,float *M1,float *M2)
{
	int i,j;

	for(i=0; i<rows; i++)
	{
		for(j=0;j<cols; j++)
		{
			(*(M2+i*cols+j)) = (*(M1+i*cols+j));
		}
	}
} 
*/

void matmulint(int m1rows, int m1cols, double *M1,double scalar, double *M3)
{
	int i,j;

	for(i=0;i<m1rows;i++)
	{
	   for(j=0;j<m1cols;j++)
	   {
			   (*(M3+i*m1cols+j)) = (*(M1+i*m1cols+j)) * scalar;
	   }
	}
}

void cross(double *v4,double *v5,double *v6)
{
	v6[0] = (double)v4[1] * (double)v5[2] - (double)v4[2] * (double)v5[1];
	v6[1] = (double)v4[2] * (double)v5[0] - (double)v4[0] * (double)v5[2];
	v6[2] = (double)v4[0] * (double)v5[1] - (double)v4[1] * (double)v5[0];
}

void dcm2euler_spt(double *mm, float *si, float *phi, float *theta)
{
	double temp;

	temp = (*(mm+2*3+1));
	if (fabs((*(mm+2*3+1))) > 1.0)
		(*(mm+2*3+1))=temp/fabs(temp);

	/* The Sequence is Si-Phi-Theta */
	*si			= atan2(-(*(mm+0*3+1)),(*(mm+1*3+1)) ) ;  /* Si  */
	*phi		= asin((*(mm+2*3+1)));                    /* Phi */
	*theta   = atan2(-(*(mm+2*3+0)),(*(mm+2*3+2)) ) ;  /* Theta */

}

void dcm2euler_stp(double *mm, float *si, float *theta, float *phi)
{
	double temp;

	temp = (*(mm+2*3+0));
	if (fabs((*(mm+2*3+0))) > 1.0)
		(*(mm+2*3+0))=temp/fabs(temp);

	/* The Sequence is Si-Phi-Theta */
	*si				= atan2((*(mm+1*3+0)),(*(mm+0*3+0)) ) ;	/* Si  */
	*theta		= asin(-(*(mm+2*3+0)));								/* theta */
	*phi			= atan2((*(mm+2*3+1)),(*(mm+2*3+2)) ) ;  /* phi */

}

void euler2dcm_spt(float si,float phi,float theta,double *mm)  /* si,phi,theta(y,r,p) */ /* rci sequence */
{
		double sz, cz;
		double sy, cy;
		double sx, cx;
	
		sz = sin(si); cz = cos(si);
		sy = sin(theta); cy = cos(theta);
		sx = sin(phi); cx = cos(phi);

		/* 1st row                              2nd row                               3rd row         */
	     *(mm+0) =  cy*cz-sy*sx*sz;	*(mm+3) = cy*sz+sy*sx*cz;		*(mm+6) =  -sy*cx;
         *(mm+1) =              -sz*cx;	*(mm+4) =              cz*cx;		*(mm+7) =       sx;
         *(mm+2) =  sy*cz+cy*sx*sz;	*(mm+5) = sy*sz-cy*sx*cz;		*(mm+8) =   cy*cx;

}

void euler2dcm_stp(float si,float theta,float phi,double *mm)  /* si,theta, phi(y,p,r) */ /* aerospace sequence */
{
		double sz, cz;
		double sy, cy;
		double sx, cx;
	
		sz = sin(si); cz = cos(si);
		sy = sin(theta); cy = cos(theta);
		sx = sin(phi); cx = cos(phi);

		/* 1st row                              2nd row                               3rd row         */
	     *(mm+0) =				cy*cz;	*(mm+3) =				  cy*sz;		*(mm+6) =    -sy;
         *(mm+1) =  sy*sx*cz-sz*cx;	*(mm+4) = sy*sx*sz+cz*cx;		*(mm+7) = cy*sx;
         *(mm+2) =  sy*cx*cz+sz*sx;	*(mm+5) = sy*cx*sz-cz*sx;		*(mm+8) = cy*cx;

}


void quat2dcm(double *quat, double *tdcm)
{
	(*(tdcm+0*3+0)) = square(*(quat+0))+square(*(quat+1))-square(*(quat+2))-square(*(quat+3));
	(*(tdcm+0*3+1)) = 2*((*(quat+1))*(*(quat+2)) - (*(quat+0))*(*(quat+3)));
	(*(tdcm+0*3+2)) = 2*((*(quat+1))*(*(quat+3)) + (*(quat+0))*(*(quat+2)));
	(*(tdcm+1*3+0)) = 2*((*(quat+1))*(*(quat+2)) + (*(quat+0))*(*(quat+3)));
	(*(tdcm+1*3+1)) = square(*(quat+0))-square(*(quat+1))+square(*(quat+2))-square(*(quat+3));
	(*(tdcm+1*3+2)) = 2*((*(quat+2))*(*(quat+3)) - (*(quat+0))*(*(quat+1)));
	(*(tdcm+2*3+0)) = 2*((*(quat+1))*(*(quat+3)) - (*(quat+0))*(*(quat+2)));
	(*(tdcm+2*3+1)) = 2*((*(quat+2))*(*(quat+3)) + (*(quat+0))*(*(quat+1)));
	(*(tdcm+2*3+2)) = square(*(quat+0))-square(*(quat+1))-square(*(quat+2))+square(*(quat+3));
}


void dcm2quat(double *tdcm, double *quat)
{
    
    
	double tr;
	double p1, p2, p3, p4;
	double mx;

	int i;

	tr = (*(tdcm+0*3+0)) + (*(tdcm+1*3+1)) + (*(tdcm+2*3+2));

	p1 = 1 + tr;
	p2 = 1 + 2 * (*(tdcm+0*3+0)) -tr;
	p3 = 1 + 2 * (*(tdcm+1*3+1)) -tr;
	p4 = 1 + 2 * (*(tdcm+2*3+2)) -tr;
	
	mx = p1;
	if(p2 > mx) mx = p2;
	if(p3 > mx) mx = p3;
	if(p4 > mx) mx = p4;
	
 	if(mx == p1)
	{
	   	(*(quat+0)) = 0.5 * sqrt(p1);
		(*(quat+1)) = ((*(tdcm+2*3+1)) - (*(tdcm+1*3+2)))/(4*(*(quat+0)));
    	(*(quat+2)) = ((*(tdcm+0*3+2)) - (*(tdcm+2*3+0)))/(4*(*(quat+0)));
    	(*(quat+3)) = ((*(tdcm+1*3+0)) - (*(tdcm+0*3+1)))/(4*(*(quat+0)));
	}
	else if(mx == p2)
	{
	  	(*(quat+1)) = 0.5 *sqrt(p2);
    	(*(quat+2)) = ((*(tdcm+1*3+0)) + (*(tdcm+0*3+1)))/(4*(*(quat+1)));
    	(*(quat+3)) = ((*(tdcm+0*3+2)) + (*(tdcm+2*3+0)))/(4*(*(quat+1)));
		(*(quat+0)) = ((*(tdcm+2*3+1)) - (*(tdcm+1*3+2)))/(4*(*(quat+1)));
	}
	else if(mx == p3)
	{
		(*(quat+2)) = 0.5 *sqrt(p3);
    	(*(quat+3)) = ((*(tdcm+2*3+1)) + (*(tdcm+1*3+2)))/(4*(*(quat+2)));
    	(*(quat+0)) = ((*(tdcm+0*3+2)) - (*(tdcm+2*3+0)))/(4*(*(quat+2)));
    	(*(quat+1)) = ((*(tdcm+1*3+0)) + (*(tdcm+0*3+1)))/(4*(*(quat+2)));
	}
	else if(mx == p4)
	{
		(*(quat+3)) = 0.5 *sqrt(p4);
    	(*(quat+0)) = ((*(tdcm+1*3+0)) - (*(tdcm+0*3+1)))/(4*(*(quat+3)));
    	(*(quat+1)) = ((*(tdcm+0*3+2)) + (*(tdcm+2*3+0)))/(4*(*(quat+3)));
    	(*(quat+2)) = ((*(tdcm+2*3+1)) + (*(tdcm+1*3+2)))/(4*(*(quat+3)));
	}
	
	if((*(quat+0)) < 0)
	{
		for(i=0;i<4;i++)
		{
			(*(quat+i)) = - (*(quat+i));
		}
	}

/*	
	(*(quat+0)) = sqrt( 1.0 + (*(tdcm+0*3+0)) + (*(tdcm+1*3+1)) + (*(tdcm+2*3+2)) )/2.0;
	(*(quat+1)) = ( (*(tdcm+2*3+1)) - (*(tdcm+1*3+2)) )/( 4.0 * (*(quat+0)) );
	(*(quat+2)) = ( (*(tdcm+0*3+2)) - (*(tdcm+2*3+0)) )/( 4.0 * (*(quat+0)) );
	(*(quat+3)) = ( (*(tdcm+1*3+0)) - (*(tdcm+0*3+1)) )/( 4.0 * (*(quat+0)) );
	
*/
}

void dcm2quat_double(double *tdcm, double *quat)
{

	double tr;
	double p1, p2, p3, p4;
	double mx;
	int i;

	tr = (*(tdcm+0*3+0)) + (*(tdcm+1*3+1)) + (*(tdcm+2*3+2));

	p1 = 1 + tr;
	p2 = 1 + 2 * (*(tdcm+0*3+0)) -tr;
	p3 = 1 + 2 * (*(tdcm+1*3+1)) -tr;
	p4 = 1 + 2 * (*(tdcm+2*3+2)) -tr;

	mx = p1;
	if(p2 > mx) mx = p2;
	if(p3 > mx) mx = p3;
	if(p4 > mx) mx = p4;

	if(mx == p1)
	{
		(*(quat+0)) = 0.5 *sqrt(p1);
		(*(quat+1)) = ((*(tdcm+2*3+1)) - (*(tdcm+1*3+2)))/(4*(*(quat+0)));
		(*(quat+2)) = ((*(tdcm+0*3+2)) - (*(tdcm+2*3+0)))/(4*(*(quat+0)));
		(*(quat+3)) = ((*(tdcm+1*3+0)) - (*(tdcm+0*3+1)))/(4*(*(quat+0)));
	}
	else if(mx == p2)
	{
		(*(quat+1)) = 0.5 *sqrt(p2);
		(*(quat+2)) = ((*(tdcm+1*3+0)) + (*(tdcm+0*3+1)))/(4*(*(quat+1)));
		(*(quat+3)) = ((*(tdcm+0*3+2)) + (*(tdcm+2*3+0)))/(4*(*(quat+1)));
		(*(quat+0)) = ((*(tdcm+2*3+1)) - (*(tdcm+1*3+2)))/(4*(*(quat+1)));
	}
	else if(mx == p3)
	{
		(*(quat+2)) = 0.5 *sqrt(p3);
		(*(quat+3)) = ((*(tdcm+2*3+1)) + (*(tdcm+1*3+2)))/(4*(*(quat+2)));
		(*(quat+0)) = ((*(tdcm+0*3+2)) - (*(tdcm+2*3+0)))/(4*(*(quat+2)));
		(*(quat+1)) = ((*(tdcm+1*3+0)) + (*(tdcm+0*3+1)))/(4*(*(quat+2)));
	}
	else if(mx == p4)
	{
		(*(quat+3)) = 0.5 *sqrt(p4);
		(*(quat+0)) = ((*(tdcm+1*3+0)) - (*(tdcm+0*3+1)))/(4*(*(quat+3)));
		(*(quat+1)) = ((*(tdcm+0*3+2)) + (*(tdcm+2*3+0)))/(4*(*(quat+3)));
		(*(quat+2)) = ((*(tdcm+2*3+1)) + (*(tdcm+1*3+2)))/(4*(*(quat+3)));
	}

	if((*(quat+0)) < 0)
	{
		for(i=0;i<4;i++)
		{
			(*(quat+i)) = - (*(quat+i));
		}
	}

	/*
	(*(quat+0)) = sqrt( 1.0 + (*(tdcm+0*3+0)) + (*(tdcm+1*3+1)) + (*(tdcm+2*3+2)) )/2.0;
	(*(quat+1)) = ( (*(tdcm+2*3+1)) - (*(tdcm+1*3+2)) )/( 4.0 * (*(quat+0)) );
	(*(quat+2)) = ( (*(tdcm+0*3+2)) - (*(tdcm+2*3+0)) )/( 4.0 * (*(quat+0)) );
	(*(quat+3)) = ( (*(tdcm+1*3+0)) - (*(tdcm+0*3+1)) )/( 4.0 * (*(quat+0)) );
	*/
}


void quat2euler_spt(double *quat,float *si, float *phi, float *theta)
{
   *si      = atan2(2*((*(quat+0))*(*(quat+3))-(*(quat+1))*(*(quat+2))), square(*(quat+0))-square(*(quat+1))+square(*(quat+2))-square(*(quat+3)));
   *phi    = asin(2*((*(quat+2))*(*(quat+3))+(*(quat+0))*(*(quat+1))));
   *theta = atan2(2*(*(quat+0)* *(quat+2)- *(quat+1)* *(quat+3)), square(*(quat+0))-square(*(quat+1))-square(*(quat+2))+square(*(quat+3)));
}


void quat2euler_stp(double *quat,float *si, float *theta, float *phi)
{
   *si			= atan2(2.0*((*(quat+1))*(*(quat+2))+(*(quat+0))*(*(quat+3))), square(*(quat+0))+square(*(quat+1))-square(*(quat+2))-square(*(quat+3)));
   *theta    = asin(-2.0*((*(quat+1))*(*(quat+3))-(*(quat+0))*(*(quat+2))));
   *phi		= atan2(2.0*(*(quat+2)* *(quat+3)+ *(quat+0)* *(quat+1)), square(*(quat+0))-square(*(quat+1))-square(*(quat+2))+square(*(quat+3)));
}


/* yaw roll pitch convention */
void euler2quat_spt(float si,float phi,float theta,double *quat)
{
		double cy,cr,cp,sy,sr,sp;

        si		  = si/2;
        phi	  = phi/2;
        theta = theta/2;

        cy = cos(si);
        sy = sin(si);
        cr = cos(phi);
        sr = sin(phi);
        cp = cos(theta);
        sp = sin(theta);

        (*(quat+0))=cy * cr * cp - sy * sr * sp;
        (*(quat+1))=cy * sr * cp - sy * cr * sp;
        (*(quat+2))=cy * cr * sp + sy * sr * cp;
        (*(quat+3))=cy * sr * sp + sy * cr * cp;
}

void euler2quat_stp(float si, float theta, float phi, double *quat)
{
	double cy,cr,cp,sy,sr,sp;

        si		  = si/2;
        phi	  = phi/2;
        theta = theta/2;

        cy = cos(si);
        sy = sin(si);
        cr = cos(phi);
        sr = sin(phi);
        cp = cos(theta);
        sp = sin(theta);

		quat[0] =  cy*cp*cr + sy*sp*sr;
        quat[1] =  cy*cp*sr - sy*sp*cr;
        quat[2] =  cy*sp*cr + sy*cp*sr;
        quat[3] =  sy*cp*cr - cy*sp*sr;

}


/*
float quat_norm(float *quat_in)//??looks like bug
{
	float norm;
	norm = square(*(quat_in+0));
	norm += square(*(quat_in+1));
	norm += square(*(quat_in+2));
	norm += square(*(quat_in+3));
	
	return norm;
}*/

void quat_norm(double *q1)
{
	int i;
	double sum;

	sum = ((q1[0] * q1[0]) + (q1[1] * q1[1]) + (q1[2] * q1[2]) + (q1[3] * q1[3]));
	sum = sqrt(sum);

	for (i = 0; i < 4; i++)
		q1[i] = q1[i] / sum;

}
/*
double quat_mod(double *quat_in)
{
	float mod;
	mod = sqrt(quat_norm(quat_in));
	return mod;
}
*/

void quat_normalize(double *quat_in)
{
	double err_coeff;
	err_coeff = square(*(quat_in+0));
	err_coeff += square(*(quat_in+1));
	err_coeff += square(*(quat_in+2));
	err_coeff += square(*(quat_in+3));
	err_coeff = err_coeff -1;
	err_coeff = 0.5*err_coeff;

	(*(quat_in+0)) = (1-err_coeff) * (*(quat_in+0));
	(*(quat_in+1)) = (1-err_coeff) * (*(quat_in+1));
	(*(quat_in+2)) = (1-err_coeff) * (*(quat_in+2));
	(*(quat_in+3)) = (1-err_coeff) * (*(quat_in+3));
}

/*
void rot2quat(float *rot, double *quat)
{
	float rot_abs;
	float sin_term;
	

	rot_abs = sqrt( (*(rot+0))*(*(rot+0)) + (*(rot+1))*(*(rot+1)) + (*(rot+2))*(*(rot+2)) );
	
	if(rot_abs == 0.0)
	{
	    (*(quat+0)) = 1.0;
		(*(quat+1)) = 0.0;
		(*(quat+2)) = 0.0;
		(*(quat+3)) = 0.0;
		
		return;
	}

	sin_term = sin(rot_abs/2.0)/rot_abs;

	(*(quat+0)) = cos(rot_abs/2.0);
	(*(quat+1)) = sin_term * (*(rot+0));
	(*(quat+2)) = sin_term * (*(rot+1));
	(*(quat+3)) = sin_term * (*(rot+2));
}

void quat2rot(double *quat, float *rot)
{
    float norm_qv;
	float tterm;

	if(*(quat+0) < 0.999)
	{
		norm_qv = sqrt(1.0 - square(*(quat+0)));
		tterm	  = 2.0 * acos(*(quat+0)) / norm_qv;

		*(rot+0)   = tterm * *(quat+1);
		*(rot+1)   = tterm * *(quat+2);
		*(rot+2)   = tterm * *(quat+3);
	}
	else
	{
		*(rot+0)   = 2.0  * (*(quat+1));
		*(rot+1)   = 2.0  * (*(quat+2));
		*(rot+2)   = 2.0  * (*(quat+3));
	}
    
}
*/
void quat_mult(double *q1, double *q2, double *q3)
{
	(*(q3+0)) = (*(q1+0)) * (*(q2+0)) - (*(q1+1)) * (*(q2+1)) \
				- (*(q1+2)) * (*(q2+2)) - (*(q1+3)) * (*(q2+3));
	(*(q3+1)) = (*(q1+0)) * (*(q2+1)) + (*(q1+1)) * (*(q2+0)) \
				+ (*(q1+2)) * (*(q2+3)) - (*(q1+3)) * (*(q2+2));
	(*(q3+2)) = (*(q1+0)) * (*(q2+2)) + (*(q1+2)) * (*(q2+0)) \
				+ (*(q1+3)) * (*(q2+1)) - (*(q1+1)) * (*(q2+3));
	(*(q3+3)) = (*(q1+0)) * (*(q2+3)) + (*(q1+3)) * (*(q2+0)) \
				+ (*(q1+1)) * (*(q2+2)) - (*(q1+2)) * (*(q2+1));
}
/*
void quat_transform(double *quat, float *in_vec, float *out_vec)
{

	float inv_quat[4];
	float temp_in_vec[4];
	float temp_q[4];
	float temp_out_vec[4];

	int i;

	inv_quat[0] = (*(quat+0));
	temp_in_vec[0] = 0.0;
	temp_out_vec[0] = 0.0;

	for(i=1;i<4;i++)
	{
		inv_quat[i] = -(*(quat+i));
		temp_in_vec[i] = (*(in_vec+i-1));
		temp_out_vec[i] = 0.0;
	}

	quat_mult((float *)temp_in_vec, (float *)inv_quat, (float *)temp_q);
	quat_mult((float *)quat, (float *)temp_q, (float *)temp_out_vec);

	for(i=0;i<3;i++)
		(*(out_vec+i)) = temp_out_vec[i+1];
}

*/

void quat_inv(double *q_in, double *q_inv)
{
	(*(q_inv+0)) =  (*(q_in + 0));
	(*(q_inv+1)) = -(*(q_in + 1));
	(*(q_inv+2)) = -(*(q_in + 2));
	(*(q_inv+3)) = -(*(q_in + 3));
}






void ecef2ned_dcm(double latm, double longm, double *Cecef2ned)
{
	double sin_lat, cos_lat;
	double sin_long, cos_long;

	sin_lat = sin(latm);
	sin_long = sin(longm);
	cos_lat = cos(latm);
	cos_long = cos(longm);

	*(Cecef2ned+0*3+0) = -sin_lat * cos_long;
	*(Cecef2ned+0*3+1) = -sin_lat * sin_long;
	*(Cecef2ned+0*3+2) = cos_lat;
	*(Cecef2ned+1*3+0) = -sin_long;
	*(Cecef2ned+1*3+1) = cos_long;
	*(Cecef2ned+1*3+2) = 0.0;
	*(Cecef2ned+2*3+0) = -cos_lat*cos_long;
	*(Cecef2ned+2*3+1) = -cos_lat*sin_long;
	*(Cecef2ned+2*3+2) = -sin_lat;
}

void ecef2ned(double latm, double longm, double *vin, double *vout)
{
	double Cecef2ned[3][3];

	ecef2ned_dcm(latm, longm, (double *)Cecef2ned);

	matmul(3,3,(double *)Cecef2ned, 3,1,vin, vout);

}

void ned2ecef_dcm(double latm, double longm, double *Cned2ecef)
{
	double sin_lat, cos_lat;
	double sin_long, cos_long;

	sin_lat = sin(latm);
	sin_long = sin(longm);
	cos_lat = cos(latm);
	cos_long = cos(longm);

	*(Cned2ecef+0*3+0) = -sin_lat * cos_long;
	*(Cned2ecef+0*3+1) = -sin_long;
	*(Cned2ecef+0*3+2) = -cos_lat * cos_long;
	*(Cned2ecef+1*3+0) = -sin_lat * sin_long;
	*(Cned2ecef+1*3+1) = cos_long;
	*(Cned2ecef+1*3+2) = -cos_lat * sin_long;
	*(Cned2ecef+2*3+0) =  cos_lat;
	*(Cned2ecef+2*3+1) = 0.0;
	*(Cned2ecef+2*3+2) = -sin_lat;
}

void ned2ecef(double latm, double longm, double *vin, double *vout)
{
	double Cned2ecef[3][3];

	ned2ecef_dcm(latm, longm, (double *)Cned2ecef);

	matmul(3,3,(double *)Cned2ecef, 3,1,vin, vout);

}

void ned2ecef_q(double latm, double longm, double *qned2ecef)
{
	double sin_lat, cos_lat;
	double sin_long, cos_long;

	sin_lat = sin(-pi/4.0 - latm/2.0);
	sin_long = sin(longm/2.0);
	cos_lat = cos(-pi/4.0 - latm/2.0);
	cos_long = cos(longm/2.0);

	*(qned2ecef+0) = cos_lat * cos_long;
	*(qned2ecef+1) = -sin_lat * sin_long;
	*(qned2ecef+2) = sin_lat * cos_long;
	*(qned2ecef+3) = cos_lat * sin_long;
}

void ecef2ned_q(double latm, double longm, double *qecef2ned)
{
	double sin_lat, cos_lat;
	double sin_long, cos_long;

	sin_lat = sin(-pi/4.0 - latm/2.0);
	sin_long = sin(longm/2.0);
	cos_lat = cos(-pi/4.0 - latm/2.0);
	cos_long = cos(longm/2.0);

	*(qecef2ned+0) = cos_lat * cos_long;
	*(qecef2ned+1) = sin_lat * sin_long;
	*(qecef2ned+2) = -sin_lat * cos_long;
	*(qecef2ned+3) = -cos_lat * sin_long;
}




void nwv2ecef_q(double lat1, double long1, double *q4)
{
	double sin_long,sin_lat,cos_long,cos_lat,sqrt_2;

	sqrt_2   = sqrt(2.0);
	sin_long = sin(long1/2.0);
	cos_long = cos(long1/2.0);
	sin_lat  = sin(lat1 /2.0);
	cos_lat  = cos(lat1 /2.0);

	q4[0]=   (sin_long/sqrt_2)*(cos_lat + sin_lat);
	q4[1]=  -(cos_long/sqrt_2)*(cos_lat - sin_lat);
	q4[2]=  -(sin_long/sqrt_2)*(cos_lat - sin_lat);
	q4[3]=  -(cos_long/sqrt_2)*(cos_lat + sin_lat);
}


void ecef2nwv_q(double t_lat,double t_long,double *q)
{

	double sqtwo;
	double sin_lat_bytwo1,cos_lat_bytwo1, sin_long_bytwo1,cos_long_bytwo1;

	sqtwo = sqrt(2.0);

	sin_lat_bytwo1  = sin(t_lat/2);
	cos_lat_bytwo1  = cos(t_lat/2);
	cos_long_bytwo1 = cos(t_long/2);
	sin_long_bytwo1 = sin(t_long/2);

	q[0] =  (sin_long_bytwo1/sqtwo)*(cos_lat_bytwo1+sin_lat_bytwo1);
	q[1] =  (cos_long_bytwo1/sqtwo)*(cos_lat_bytwo1-sin_lat_bytwo1);
	q[2] =  (sin_long_bytwo1/sqtwo)*(cos_lat_bytwo1-sin_lat_bytwo1);
	q[3] =  (cos_long_bytwo1/sqtwo)*(cos_lat_bytwo1+sin_lat_bytwo1);
}

void ecef2nwv_dcm(double latm, double longm, double *Cecef2nwv)
{
	double sin_lat, cos_lat;
	double sin_long, cos_long;

	sin_lat = sin(latm);
	sin_long = sin(longm);
	cos_lat = cos(latm);
	cos_long = cos(longm);

	*(Cecef2nwv+0*3+0) = -sin_lat * cos_long;
	*(Cecef2nwv+0*3+1) = -sin_lat * sin_long;
	*(Cecef2nwv+0*3+2) = cos_lat;
	*(Cecef2nwv+1*3+0) = sin_long;
	*(Cecef2nwv+1*3+1) = -cos_long;
	*(Cecef2nwv+1*3+2) = 0.0;
	*(Cecef2nwv+2*3+0) = cos_lat * cos_long;
	*(Cecef2nwv+2*3+1) = cos_lat * sin_long;
	*(Cecef2nwv+2*3+2) = sin_lat;
}

void ecef2nwv(double latm, double longm, double *vin, double *vout)
{
	double Cecef2nwv[3][3];

	ecef2nwv_dcm(latm, longm, (double *)Cecef2nwv);

	matmul(3,3,(double *)Cecef2nwv, 3,1,vin, vout);

}


void nwv2ecef_dcm(double latm, double longm, double *Cnwv2ecef)
{
	double sin_lat, cos_lat;
	double sin_long, cos_long;

	sin_lat = sin(latm);
	sin_long = sin(longm);
	cos_lat = cos(latm);
	cos_long = cos(longm);

	*(Cnwv2ecef+0*3+0) = -sin_lat * cos_long;
	*(Cnwv2ecef+0*3+1) = sin_long;
	*(Cnwv2ecef+0*3+2) = cos_lat * cos_long;
	*(Cnwv2ecef+1*3+0) = -sin_lat * sin_long;
	*(Cnwv2ecef+1*3+1) = -cos_long;
	*(Cnwv2ecef+1*3+2) = cos_lat * sin_long;
	*(Cnwv2ecef+2*3+0) = cos_lat;
	*(Cnwv2ecef+2*3+1) = 0.0;
	*(Cnwv2ecef+2*3+2) = sin_lat;

}

void nwv2ecef(double latm, double longm, double *vin, double *vout)
{
	double Cnwv2ecef[3][3];

	nwv2ecef_dcm(latm, longm, (double *)Cnwv2ecef);

	matmul(3,3,(double *)Cnwv2ecef, 3,1,vin, vout);

}
/*
void print_matrix(int rows, int cols, double *mat, char *name)
{

	int i,j;

	printf("\n%s\n", name);
	for(i = 0;i<rows; i++)
	{
		printf("\n\t");
		for(j=0;j<cols; j++)
			printf("%e ", *(mat+i*cols+j));
	}
}
*/

/****************************************************/


/*
  function to compute lu decomposition of a given matrix
  input : temp_M - matrix to be lu decomposed, size of input
  output : indx - matrix containing lu indices
  result stored in temp_M
*/
/*
void ludcmp(int size, double *temp_M,int *indx)
{
	int i, imax, j, k, d;
	int rs;
	double big, dum, sum, temp;
	//double vv[6];
	double *vv;
	double TINY = 1.0e-20;


	vv = (double *)calloc(size,sizeof(double));


	rs = size;
	d = 1;

		
	for(i=0; i< rs; i++)
	{
		big = 0.0;
		for(j = 0; j < rs; j++)
			   if((temp=fabs((*(temp_M+i*rs+j)))) > big) big = temp;
		 		
		if(big == 0.0)
        	return;
        // printf("matrix is singular to working precision\n"); 
		vv[i] = 1.0/big;
	}
	
	

	for(j = 0; j< rs; j++)
	{
		for(i=0; i< j; i++)
		{
			sum = (*(temp_M+i*rs+j)); 
			for(k=0; k<i; k++) sum -= (*(temp_M+i*rs+k)) * (*(temp_M+k*rs+j));
			(*(temp_M+i*rs+j)) = sum;
		}
		big = 0.0;
		for(i=j; i<rs; i++)
		{
			sum = (*(temp_M+i*rs+j));
			for(k=0; k<j; k++)
				sum -= (*(temp_M+i*rs+k)) * (*(temp_M+k*rs+j));
			(*(temp_M+i*rs+j)) = sum;
			if((dum = vv[i] * fabs(sum)) >= big)
			{
				big = dum;
				imax = i;
			}
		}
		if(j!=imax)
		{
			for(k=0; k<rs; k++)
			{
				dum = (*(temp_M+imax*rs+k));
				(*(temp_M+imax*rs+k)) = (*(temp_M+j*rs+k));
				(*(temp_M+j*rs+k)) =  dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if((*(temp_M+j*rs+j)) == 0.0)
		{
		 	(*(temp_M+j*rs+j)) = TINY;
		    //printf("\n Warning : Matrix is singular to working precision \n");
		}
		if( j!= rs) 
		{
			dum = 1.0/(*(temp_M+j*rs+j));
			for(i=j+1; i<rs; i++) (*(temp_M+i*rs+j)) *= dum;
		}	
	}


	free(vv);
}
*/

/* 
  lu back substitution
  A x = b
  input : temp_M - matrix A 
          indx - matrix of indices
          b - the b vector
*/
/*
void lubksb(int size, double *temp_M,int *indx, double *b)
{
	int i;
	int ip, j;
	double sum;
	int rs;
	rs = size;
	

	for(i=0; i<rs; i++)
	{
		ip = (int)(*(indx+i));
		sum = (*(b+ip));
		(*(b+ip)) = (*(b+i));
		for(j=0; j < i; j++)  sum -= (*(temp_M+i*rs+j))* (*(b+j));
			   
			
		
		(*(b+i)) = sum;
		
	}

	for(i=rs-1; i>=0; i--)
	{
		sum = (*(b+i));
		for(j=i+1; j<rs; j++) 
			sum -= (*(temp_M+i*rs+j)) * (*(b+j));
		(*(b+i)) = sum/(*(temp_M+i*rs+i));
	}
}
*/
/* 
  input : the sqaure matrix, size of the matrix (rows or cols)
  output : the inverse matrix
*/
/*
void invert_splu(int size, double *temp_M,double *M_inv)
{
	int i,j;

	
	double  *M_lud;
	int	   *indx;
	double  *col;
	

	int rs = size;
	
	
	//double M_lud[36];
	//int    indx[6];
	//double col[6];
	
		M_lud = (double *)calloc(size*size, sizeof(double));
		indx  = (int *)calloc(size, sizeof(int));
		col   = (double *)calloc(size, sizeof(double));
	

	for(i=0;i<size*size;i++)
			M_lud[i] = *(temp_M+i);
		
	ludcmp(size,(double *)M_lud, (int *)indx);
	
	for(j=0;j<rs; j++)
	{
		for(i=0; i<rs; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(size,(double *)M_lud,(int *)indx, (double *)col);
		for(i=0; i<rs; i++)
			*(M_inv+i*rs+j) = (double)col[i];
		     
	}


	free(indx);
	free(col);
	free(M_lud);

}

*/
/**********************************************************
 * Function Name: to_cross_form
 * Input: Vector(3x1)
 * Output: Matrix(3x3)
 ***********************************************************/
void to_cross_form(double *matrix, double *vector)
{
	(*(matrix + 0 * 3 + 0)) = 0.0;
	(*(matrix + 0 * 3 + 1)) = -(*(vector + 2));
	(*(matrix + 0 * 3 + 2)) = (*(vector + 1));

	(*(matrix + 1 * 3 + 0)) = (*(vector + 2));
	(*(matrix + 1 * 3 + 1)) = 0.0;
	(*(matrix + 1 * 3 + 2)) = -(*(vector + 0));

	(*(matrix + 2 * 3 + 0)) = -(*(vector + 1));
	(*(matrix + 2 * 3 + 1)) = (*(vector + 0));
	(*(matrix + 2 * 3 + 2)) = 0.0;
}
/*
void init_m(int size, float *mat, ...)
{
	int i,j;
	float val;
	va_list ap;

	va_start(ap, mat);

	for(i=0;i<size; i++)
	{
		for(j=0;j<size;j++)
		{
			*(mat+i*size+j)  = va_arg(ap,float);
		}
	}

	va_end(ap);
}
*/
float trace(int size, float *mat)
{
	int i;
	float tr;

	tr = 0.0;

	for(i=0;i<size;i++)
	{
		tr += mat[i*size+i];
	}

	return tr;
}

float max(int size, float *vec)
{
	int i;
	float maxi;

	maxi = *(vec+0);

	for(i=1;i<size;i++)
		if(*(vec+i)>maxi)
			maxi = *(vec+i);

	return maxi;

}

double norm(int size, double *vec)
{
	int i;
	double nr;

	nr = 0.0;
	for(i=0;i<size;i++)
		nr += square(vec[i]);

	nr = sqrt(nr);

	return nr;
}

void init(double a, double b, double c, double *v3)
{
	v3[0] = a;
	v3[1] = b;
	v3[2] = c;
}



double gravity_ecef(double *p_ecef_pos, double *p_ecef_gravity)
{
	double rad1,rad2;
	double vv[3],vv1[3],r_vec[3];
	double omega[3];
	double p_g_ecef_mag;

	double p_R;

	p_R = norm(3,p_ecef_pos);

	rad1 = (R0/p_R) * (R0/p_R);
	rad2 = p_ecef_pos[2]/p_R;

	p_ecef_gravity[0] = -(mue/(R0*R0))*(p_ecef_pos[0]/p_R)*((rad1*(1.0+epsilon))+(jconst2*(rad1*rad1))*(1.0-5.0*(rad2*rad2)));
	p_ecef_gravity[1] = -(mue/(R0*R0))*(p_ecef_pos[1]/p_R)*((rad1*(1.0+epsilon))+(jconst2*(rad1*rad1))*(1.0-5.0*(rad2*rad2)));
	p_ecef_gravity[2] = -(mue/(R0*R0))*(p_ecef_pos[2]/p_R)*((rad1*(1.0+epsilon))+(jconst2*(rad1*rad1))*(3.0-5.0*(rad2*rad2)));

  
	init(p_ecef_pos[0],p_ecef_pos[1],p_ecef_pos[2],r_vec);
	init(0.0,0.0,0.0,vv);
	init(0.0,0.0,w_e, omega);
	cross(omega,r_vec,vv);
	cross(omega,vv,vv1);
	matsub(3,1,(double *)p_ecef_gravity, (double *)vv1, (double *)p_ecef_gravity);
	p_g_ecef_mag = norm(3,p_ecef_gravity);
	return p_g_ecef_mag;
}
/*
void gr_britting_ned(float latd, float alt, float *Ned_gravity)
{
   int i;
   float  co_lat,r0,j3,j4;
   float  g_colat,g_rad;
   float  vv[3],vv1[3],rad[3],grav[3];
   float  latc, dev_ang;
   float  omega_ned[3][1];
   float	radi, radialt;
   
   j3 = -2.3e-6;
   j4 = -1.8e-6;
   
   r0		= R;
   radi		= R * (1.0-flat * square(sin(latd)));
   radialt	= radi + alt;

   latc		= (latd - flat * sin(2.0*latd));

   co_lat	= (pi/2.0)  -  latc;
   
   // colatitude component
   g_colat = 3.0*(mue/(radialt*radialt))*(pow((r0/radialt),2))*sin(co_lat)*cos(co_lat)*( jconst + ( (j3/2.0)*(r0/radialt)*(1/cos(co_lat))*(5.0*(pow((cos(co_lat)),2))-1) ) + ( (5.0/6.0)*j4*(pow((r0/r),2))*(7.0*(pow((cos(co_lat)),2))-3.0) )  );  
   // radial component
   g_rad   = -(mue/(radialt*radialt))*( 1 - ( (3.0/2.0)*jconst*(pow((r0/radialt),2))*(3.0*(pow((cos(co_lat)),2))-1) ) - ( 2.0*j3*(pow((r0/radialt),3))*(cos(co_lat))*(5.0*(pow((cos(co_lat)),2))-3.0) ) - ( (5.0/8.0)*j4*(pow((r0/radialt),4))*(35.0*(pow((cos(co_lat)),4))-30.0*(pow((cos(co_lat)),2))+3.0) )  );  

   dev_ang  = (latd - latc);
 
   grav[0] = -g_colat * cos(dev_ang) - g_rad * sin(dev_ang);
   grav[1] =  0.0;
   grav[2] =  g_colat * sin(dev_ang) - g_rad * cos(dev_ang);
   
   rad[0] = -radi * sin(dev_ang);
   rad[1] =  0.0;
   rad[2] = -(radi * cos(dev_ang) + alt);

   omega_ned[0][0] = w_e*cos(latd);
   omega_ned[1][0] = 0.0;
   omega_ned[2][0] = -w_e*sin(latd);
   
   cross((float *)omega_ned,rad,vv);
   cross((float *)omega_ned,vv,vv1);

   for(i=0; i<3;i++)
		Ned_gravity[i] = grav[i] - vv1[i];  // gravity : g - wxwxr    
}
*/

float GeoCentric_to_GeoDetic( float Lc)
{
	int dcount = 0;
	float error, Ld_temp,Ld;
	
	error = 1;
	Ld = 0;//FIXME Ld default value to be put
	Ld_temp = Lc;
    while( (error>1.0e-10) && (dcount <4) )
	{
	   dcount++;
	   Ld = Lc + flat*sin(2.0*Ld_temp);
	   error = fabs( Ld_temp - Ld );
	   Ld_temp = Ld;
	}

	return Ld;
}

/*
float GeoDetic_to_GeoCentric(float Ld, float alt)
{
	float ls; // geodetic latitude at surface 
	float sgd, cgd, sls, cls;

	float Lc;

	float sinlam, geoc_radius;

	ls = atan2( (1.0-square(eccen))*sin(Ld), cos(Ld) );

	sgd = sin(Ld);
	cgd = cos(Ld);
	sls = sin(ls);
	cls = cos(ls);

	sinlam = sls;
	geoc_radius  = sqrt(( square(R) )/( 1 + (1/(square( 1 - f )) - 1)*square(sinlam) ));

	Lc = atan2(geoc_radius*sls + alt*sgd, geoc_radius*cls + alt*cgd);

	return Lc;
}

*/
/* lat long in radians, alt in meters */
void lla2ecef(double lat, double longm, double alt, double ecef_pos[3])
{
	double e2;
	double sinphi, cosphi;
	//double temp;

	e2 = square(e);
	sinphi = sin(lat);
	cosphi = cos(lat);

	//temp = sqrt(1.0);
	//temp = N(lat);
	*(ecef_pos+0) = (N(lat) + alt) * cosphi * cos(longm);
	*(ecef_pos+1) = (N(lat) + alt) * cosphi * sin(longm);
	*(ecef_pos+2) = (N(lat) * (1 - e2) + alt) * sinphi;
}


#if 0
/* lat long in radians, alt in meters */
void ecef2lla(double *ecef_pos, double *latm, double *longm, double *altm)
{

	double a, b, e2, ep2;
	double x, y, z;

	double rho, beta, phi;
	double sin_beta, cos_beta;
	double betaNew;
	double sinphi;
	unsigned int count;
	
	a  = R;					/* Semimajor axis */
	e2 = square(eccen);		/* Square of first eccentricity */
	ep2 = e2 / (1 - e2);    /* Square of second eccentricity */
	b = r;			        /* Semiminor axis */

	x = *(ecef_pos+0);
	y = *(ecef_pos+1);
	z = *(ecef_pos+2);

	/* extract  Longitude */
	*longm = atan2(y,x);

	/* Distance from Z-axis */
	rho = sqrt(square(x)+square(y));

	/* Bowring's formula for initial parametric (beta) and geodetic (phi) latitudes */
	beta = atan2(z, (1 - f) * rho);
	sin_beta = sin(beta);
	cos_beta = cos(beta);
	phi = atan2(z + b * ep2 * square(sin_beta)*sin_beta,rho - a * e2  * square(cos_beta) * cos_beta);

	/* 
	 * Fixed-point iteration with Bowring's formula
	 * (typically converges within two or three iterations)
	 */
	betaNew = atan2((1 - f)*sin(phi), cos(phi));
	count = 0;
	while(beta != betaNew && count < 5)
	{
		beta = betaNew;
		phi = atan2(z + b * ep2 * square(sin_beta)*sin_beta,rho - a * e2  * square(cos_beta) * cos_beta);
		betaNew = atan2((1 - f)*sin(phi), cos(phi));
		count = count + 1;
	}

	*latm = phi;

	/* Calculate ellipsoidal height from the final value for latitude */
	sinphi = sin(phi);
	*altm = rho * cos(phi) + (z + e2 * N(phi) * sinphi) * sinphi - N(phi);
}



void rot2dcm(float *rot, float *dcm)
{
	int i,j;
	float mag_rot;
	float cos_mag_rot;
	float sin_mag_rot;
	float unit_rot[3];

	mag_rot = norm(3,rot);


	if(mag_rot < 1e-15)
	{
		cos_mag_rot = 1.0;
		sin_mag_rot = 0.0;
	}
	else
	{
		cos_mag_rot = cos(mag_rot);
		sin_mag_rot = sin(mag_rot);
	}

	/* derive the unit vector from the rotation vector */
	unit_rot[0] = rot[0]/mag_rot;
	unit_rot[1] = rot[1]/mag_rot;
	unit_rot[2] = rot[2]/mag_rot;

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			if(i == j)
			{
				*(dcm+i*3+j) = 1.0 * cos_mag_rot;
			}
			else
			{
				*(dcm+i*3+j) = 0.0;
			}

			if(cos_mag_rot != 1.0)
			{
				*(dcm+i*3+j) += (1.0 - cos_mag_rot) * unit_rot[i] * unit_rot[j];
				/* cross product form 
				 * of sin(mag_rot) * to_cross_form(unit_rot)
				 */
				if(i==0 && j==1)
					*(dcm+i*3+j) += sin_mag_rot * (-unit_rot[2]);
				else if(i==0 && j==2)
					*(dcm+i*3+j) += sin_mag_rot * (unit_rot[1]);
				else if(i==1 && j==0)
					*(dcm+i*3+j) += sin_mag_rot * (unit_rot[2]);
				else if(i==1 && j==2)
					*(dcm+i*3+j) += sin_mag_rot * (-unit_rot[0]);
				else if(i==2 && j==0)
					*(dcm+i*3+j) += sin_mag_rot * (-unit_rot[1]);
				else if(i==2 && j==1)
					*(dcm+i*3+j) += sin_mag_rot * (unit_rot[0]);

			}
		}
	}

}


void dcm2rot(float *dcm, float *rot)
{
	int i,j;
	float A[3][3],S[3][3];
	float magn_rot;
	float trm;
	float k;

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			A[i][j] = 0.5 * ( dcm[i*3+j] - dcm[j*3+i]);
			S[i][j] = 0.5 * ( dcm[i*3+j] + dcm[j*3+i]);
		}
	}

	magn_rot = atan2(sqrt(square(A[2][1])+square(A[0][2])+square(A[1][0])),((trace(3, (float *)dcm)-1.0)/2));


	if(sin(magn_rot) > 1e-05)
	{
    		trm = magn_rot/sin(magn_rot);

    		*(rot+0) = trm * A[2][1];
    		*(rot+1) = trm * A[0][2];
    		*(rot+2) = trm * A[1][0];
		return;
	}
	else
	{
		k = S[0][0];
		if( S[1][1] > k)
			k = S[1][1];
		if( S[2][2] > k)
			k = S[2][2];
    
    	if(k == S[0][0])
		{
			*(rot+0) = sqrt((S[0][0]+1.0)/2.0);
          	*(rot+1) = S[0][1]/(2.0* *(rot+0));
          	*(rot+2) = S[0][2]/(2.0* *(rot+0));
		}
		else if(k == S[1][1])
		{
        	*(rot+1) = sqrt((S[1][1]+1.0)/2.0);   
        	*(rot+0) = S[0][1]/(2.0* *(rot+1));
        	*(rot+2) = S[1][2]/(2.0* *(rot+1));
		}
		else if(k == S[2][2])
		{
        	*(rot+2) = sqrt((S[2][2]+1.0)/2.0);   
        	*(rot+0) = S[0][2]/(2.0* *(rot+2));
        	*(rot+1) = S[1][2]/(2.0* *(rot+2));
		}


		for(i=0;i<3;i++)
	    		*(rot+i) = magn_rot * *(rot+i);

	}
	return;
}
#endif
void ned2nwv_dcm(double *dcm)
{
	int i, j;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			if(i==j)
				*(dcm+i*3+j) = 1.0;
			else
				*(dcm+i*3+j) = 0.0;
		}

	*(dcm+4) = -1.0;
	*(dcm+8) = -1.0;
}

void nwv2ned_dcm(double *dcm)
{
	int i, j;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			if(i==j) *(dcm+i*3+j) = 1.0;
			else *(dcm+i*3+j) = 0.0;
		}

	*(dcm+4) = -1.0;
	*(dcm+8) = -1.0;
}

void ned2nwv_q(double *q)
{
	*(q+0) = 0.0;
	*(q+1) = 1.0;
	*(q+2) = 0.0;
	*(q+3) = 0.0;
}

void nwv2ned_q(double *q)
{
	*(q+0) = 0.0;
	*(q+1) = 1.0;
	*(q+2) = 0.0;
	*(q+3) = 0.0;
}

void ned2nwv(double *vned, double *vnwv)
{
	
	double dcm[3][3] = {{1.0, 0.0, 0.0}, {0.0, -1.0,  0.0}, {0.0, 0.0, -1.0}};

	matmul(3,3,(double *)dcm, 3, 1, vned, vnwv);

}

void nwv2ned(double *vnwv, double *vned)
{
	
	double dcm[3][3] = {{1.0, 0.0, 0.0}, {0.0, -1.0,  0.0}, {0.0, 0.0, -1.0}};

	matmul(3,3,(double *)dcm, 3, 1, vnwv, vned);

}


unsigned int  ComputeChecksum( unsigned int *data, unsigned int cnt)
{
	unsigned int Checksum=0x00,i;
	for(i=0;i<cnt;i++)
		Checksum = (Checksum) ^ (data[i]);
	
	return Checksum;
}
/*
void Delay(unsigned int count)
{
    unsigned int i;    
	for(i = 0; i< count; i++);
	
}*/
void invert3(int size, double *temp_M, double *M_inv)
{
	double a,b,c,d,es,fs,i,j,k;
	double det;

	a = (*(temp_M+0*3+0));
	b = (*(temp_M+0*3+1));
	c = (*(temp_M+0*3+2));
	d = (*(temp_M+1*3+0));
	es = (*(temp_M+1*3+1));
	fs = (*(temp_M+1*3+2));
	i = (*(temp_M+2*3+0));
	j = (*(temp_M+2*3+1));
	k = (*(temp_M+2*3+2));

	det = a*(es*k-j*fs) - b*(d*k-i*fs) + c*(d*j - i*es);

	if(det >= 1.0e-10)
	{
		(*(M_inv+0*3+0)) =   (es*k - j*fs)/det;
		(*(M_inv+0*3+1)) =  -(b*k - c*j)/det;
		(*(M_inv+0*3+2)) =   (b*fs - es*c)/det;
		(*(M_inv+1*3+0)) =  -(d*k - i*fs)/det;
		(*(M_inv+1*3+1)) =   (a*k - i*c)/det;
		(*(M_inv+1*3+2)) =  -(a*fs - d*c)/det;
		(*(M_inv+2*3+0)) =   (d*j - i*es)/det;
		(*(M_inv+2*3+1)) =  -(a*j - b*i)/det;
		(*(M_inv+2*3+2)) =   (a*es - d*b)/det;
	}
}

/*************************************************************************
FUNCTION NAME: matmul_double
DESCRIPTION  : matrix multiplication of double
INPUT		 : no.of rows, columns of matrix 1 and 2 and two matrix
OUTPUT		 : void
**************************************************************************/
void matmul_double(int m1rows, int m1cols, double *M1,
			int m2rows, int m2cols, double *M2,
			double *M3)
{
	int i,j,k;

	for(i=0;i<m1rows;i++)
	{
	   for(j=0;j<m2cols;j++)
	   {
			*(M3+i*m2cols+j) = (double)0.0;

			for(k=0;k<m1cols;k++)
			{

			   (*(M3+i*m2cols+j))  +=   \
	               (double)(*(M1+i*m1cols+k)) * \
	               (double)(*(M2+k*m2cols+j));
	 		}
	   }
	}
}

/*************************************************************************
FUNCTION NAME: transpose_double
DESCRIPTION  : matrix transpose_double
INPUT		 : no.of rows, columns of matrix 1 and 2 and two matrix
OUTPUT		 : void
**************************************************************************/
void transpose_double(int rows, int cols, double *M1, double *M2)
{
	int i,j;
	for(i = 0; i< rows; i++)
	{
		for(j =0; j< cols; j++)
		{
			(*(M2+j*rows+i)) = (*(M1+i*cols+j));
		}
	}
}

/*************************************************************************
FUNCTION NAME: matadd_double
DESCRIPTION  : matrix addition
INPUT		 : no.of rows, columns of matrix 1 and 2 and two matrix
OUTPUT		 : void
**************************************************************************/
void matadd_double(int rows, int cols,double *M1,double *M2,double *M3)
{
	int i,j;
	for(i=0; i<rows; i++)
	{
		for(j=0;j<cols; j++)
		{
			(*(M3+i*cols+j)) = (*(M1+i*cols+j)) + (*(M2+i*cols+j));
		}
	}
}


void matmulint_double(int m1rows, int m1cols, double *M1, double scalar, double *M3)
{
	int i,j;
	for(i=0;i<m1rows;i++)
	{
	   for(j=0;j<m1cols;j++)
	   {
			   (*(M3+i*m1cols+j)) = (*(M1+i*m1cols+j)) * scalar;
	   }
	}
}

void cross_gc( double *v4,double *v5,double *v6)
{
	v6[0] = (double)v4[1] * (double)v5[2] - (double)v4[2] * (double)v5[1];
	v6[1] = (double)v4[2] * (double)v5[0] - (double)v4[0] * (double)v5[2];
	v6[2] = (double)v4[0] * (double)v5[1] - (double)v4[1] * (double)v5[0];
}

void get_euler_spt(double *mm, double *angles)
{
	double temp;

	temp = (*(mm+2*3+1));
	if (fabs((*(mm+2*3+1))) > 1.0)
		(*(mm+2*3+1))=temp/fabs(temp);

	/* The Sequence is Si-Phi-Theta*/
	(*(angles+0))   = atan2(-(*(mm+0*3+1)),(*(mm+1*3+1)) ) ;  /* Si */
	(*(angles+1))   = asin((*(mm+2*3+1)));                    /* Phi*/
	(*(angles+2))   = atan2(-(*(mm+2*3+0)),(*(mm+2*3+2)) ) ;  /* Theta */
}

void get_euler_stp(double *mm, double *angles)
{
	double temp;

	temp = (*(mm+2*3+0));
	if (fabs((*(mm+2*3+0))) > 1.0)
		(*(mm+2*3+0))=temp/fabs(temp);

	/* The Sequence is Si-Theta-Phi*/
	(*(angles+0))   = atan2((*(mm+1*3+0)),(*(mm+0*3+0)) ) ;  /* Si*/
	(*(angles+1))   = asin(-(*(mm+2*3+0)));                  /* Theta*/
	(*(angles+2))   = atan2((*(mm+2*3+1)),(*(mm+2*3+2)) ) ;  /* Phi*/
}



void matmulvec_double(double *M1, double *M2, double *M3 )
{
	M3[0] = (*(M1+0)) * M2[0] + (*(M1+1)) * M2[1] + (*(M1+2)) * M2[2];
	M3[1] = (*(M1+3)) * M2[0] + (*(M1+4)) * M2[1] + (*(M1+5)) * M2[2];
	M3[2] = (*(M1+6)) * M2[0] + (*(M1+7)) * M2[1] + (*(M1+8)) * M2[2];
}

void rot2quat(double *rot, double *quat)
{
	double rot_abs;
	double sin_term;

	rot_abs = sqrt( (*(rot+0))*(*(rot+0)) + (*(rot+1))*(*(rot+1)) + (*(rot+2))*(*(rot+2)) );

/*	if(rot_abs <=  1e-04)         old value - 0.0 */
/*	if(rot_abs ==  0.0) */
	if(rot_abs <=  1e-6)
	{
		(*(quat+0)) = 1.0;
		(*(quat+1)) = 0.0;
		(*(quat+2)) = 0.0;
		(*(quat+3)) = 0.0;

		return;
	}

	sin_term = sin(rot_abs/2.0)/rot_abs;

	(*(quat+0)) = cos(rot_abs/2.0);
	(*(quat+1)) = sin_term * (*(rot+0));
	(*(quat+2)) = sin_term * (*(rot+1));
	(*(quat+3)) = sin_term * (*(rot+2));
}

void quatmult(double *q1, double *q2, double *q3)
{
	(*(q3+0)) = (*(q1+0)) * (*(q2+0)) - (*(q1+1)) * (*(q2+1)) \
				- (*(q1+2)) * (*(q2+2)) - (*(q1+3)) * (*(q2+3));
	(*(q3+1)) = (*(q1+0)) * (*(q2+1)) + (*(q1+1)) * (*(q2+0)) \
				+ (*(q1+2)) * (*(q2+3)) - (*(q1+3)) * (*(q2+2));
	(*(q3+2)) = (*(q1+0)) * (*(q2+2)) + (*(q1+2)) * (*(q2+0)) \
				+ (*(q1+3)) * (*(q2+1)) - (*(q1+1)) * (*(q2+3));
	(*(q3+3)) = (*(q1+0)) * (*(q2+3)) + (*(q1+3)) * (*(q2+0)) \
				+ (*(q1+1)) * (*(q2+2)) - (*(q1+2)) * (*(q2+1));
}

void quat_transform(double *quat, double *in_vec, double *out_vec)
{

	double inv_quat[4];
	double temp_in_vec[4];
	double temp_q[4];
	double temp_out_vec[4];

	int i;

	inv_quat[0] = (*(quat+0));
	temp_in_vec[0] = 0.0;
	temp_out_vec[0] = 0.0;

	for(i=1;i<4;i++)
	{
		inv_quat[i] = -(*(quat+i));
		temp_in_vec[i] = (*(in_vec+i-1));
		temp_out_vec[i] = 0.0;
	}

	quatmult((double *)temp_in_vec, (double *)inv_quat, (double *)temp_q);
	quatmult((double *)quat, (double *)temp_q, (double *)temp_out_vec);

	for(i=0;i<3;i++)
	{
		(*(out_vec+i)) = temp_out_vec[i+1];
	}
}


void invert2(int size, double *temp_M, double *M_inv)
{
	double a,b,c,d;
	double det;

	a = (*(temp_M+0*2+0));
	b = (*(temp_M+0*2+1));
	c = (*(temp_M+1*2+0));
	d = (*(temp_M+1*2+1));


	det = a * d - b * c;

	if(det >= 1.0e-10)
	{
		(*(M_inv+0*2+0)) =    d/det;
		(*(M_inv+0*2+1)) =   -b/det;
		(*(M_inv+1*2+0)) =   -c/det;
		(*(M_inv+1*2+1)) =    a/det;
	}
}
/*
 * GC UTILS
 */
void dcm_to_quat(double *tdcm, double *quat)
{

	double tr;
	double p1, p2, p3, p4;
	double mx;
	int i;

	tr = (*(tdcm+0*3+0)) + (*(tdcm+1*3+1)) + (*(tdcm+2*3+2));

	p1 = 1 + tr;
	p2 = 1 + 2 * (*(tdcm+0*3+0)) -tr;
	p3 = 1 + 2 * (*(tdcm+1*3+1)) -tr;
	p4 = 1 + 2 * (*(tdcm+2*3+2)) -tr;

	mx = p1;
	if(p2 > mx) mx = p2;
	if(p3 > mx) mx = p3;
	if(p4 > mx) mx = p4;

	if(mx == p1)
	{
		(*(quat+0)) = 0.5 *sqrt(p1);
		(*(quat+1)) = ((*(tdcm+2*3+1)) - (*(tdcm+1*3+2)))/(4*(*(quat+0)));
		(*(quat+2)) = ((*(tdcm+0*3+2)) - (*(tdcm+2*3+0)))/(4*(*(quat+0)));
		(*(quat+3)) = ((*(tdcm+1*3+0)) - (*(tdcm+0*3+1)))/(4*(*(quat+0)));
	}
	else if(mx == p2)
	{
		(*(quat+1)) = 0.5 *sqrt(p2);
		(*(quat+2)) = ((*(tdcm+1*3+0)) + (*(tdcm+0*3+1)))/(4*(*(quat+1)));
		(*(quat+3)) = ((*(tdcm+0*3+2)) + (*(tdcm+2*3+0)))/(4*(*(quat+1)));
		(*(quat+0)) = ((*(tdcm+2*3+1)) - (*(tdcm+1*3+2)))/(4*(*(quat+1)));
	}
	else if(mx == p3)
	{
		(*(quat+2)) = 0.5 *sqrt(p3);
		(*(quat+3)) = ((*(tdcm+2*3+1)) + (*(tdcm+1*3+2)))/(4*(*(quat+2)));
		(*(quat+0)) = ((*(tdcm+0*3+2)) - (*(tdcm+2*3+0)))/(4*(*(quat+2)));
		(*(quat+1)) = ((*(tdcm+1*3+0)) + (*(tdcm+0*3+1)))/(4*(*(quat+2)));
	}
	else if(mx == p4)
	{
		(*(quat+3)) = 0.5 *sqrt(p4);
		(*(quat+0)) = ((*(tdcm+1*3+0)) - (*(tdcm+0*3+1)))/(4*(*(quat+3)));
		(*(quat+1)) = ((*(tdcm+0*3+2)) + (*(tdcm+2*3+0)))/(4*(*(quat+3)));
		(*(quat+2)) = ((*(tdcm+2*3+1)) + (*(tdcm+1*3+2)))/(4*(*(quat+3)));
	}

	if((*(quat+0)) < 0)
	{
		for(i=0;i<4;i++)
		{
			(*(quat+i)) = - (*(quat+i));
		}
	}

	/*
	(*(quat+0)) = sqrt( 1.0 + (*(tdcm+0*3+0)) + (*(tdcm+1*3+1)) + (*(tdcm+2*3+2)) )/2.0;
	(*(quat+1)) = ( (*(tdcm+2*3+1)) - (*(tdcm+1*3+2)) )/( 4.0 * (*(quat+0)) );
	(*(quat+2)) = ( (*(tdcm+0*3+2)) - (*(tdcm+2*3+0)) )/( 4.0 * (*(quat+0)) );
	(*(quat+3)) = ( (*(tdcm+1*3+0)) - (*(tdcm+0*3+1)) )/( 4.0 * (*(quat+0)) );
	*/
}
void euler2dcm_spt_double(float si,float phi,float theta,double *mm)  /* si,phi,theta(y,r,p) */ /* rci sequence */
{
		float sz, cz;
		float sy, cy;
		float sx, cx;

		sz = sin(si); cz = cos(si);
		sy = sin(theta); cy = cos(theta);
		sx = sin(phi); cx = cos(phi);

		/* 1st row                              2nd row                               3rd row         */
	     *(mm+0) =  cy*cz-sy*sx*sz;	*(mm+3) = cy*sz+sy*sx*cz;		*(mm+6) =  -sy*cx;
         *(mm+1) =              -sz*cx;	*(mm+4) =              cz*cx;		*(mm+7) =       sx;
         *(mm+2) =  sy*cz+cy*sx*sz;	*(mm+5) = sy*sz-cy*sx*cz;		*(mm+8) =   cy*cx;

}
void quat_to_dcm(double *quat, double *tdcm)
{
	(*(tdcm+0*3+0)) = square_gc(*(quat+0))+square_gc(*(quat+1))-square_gc(*(quat+2))-square_gc(*(quat+3));
	(*(tdcm+0*3+1)) = 2*((*(quat+1))*(*(quat+2)) - (*(quat+0))*(*(quat+3)));
	(*(tdcm+0*3+2)) = 2*((*(quat+1))*(*(quat+3)) + (*(quat+0))*(*(quat+2)));
	(*(tdcm+1*3+0)) = 2*((*(quat+1))*(*(quat+2)) + (*(quat+0))*(*(quat+3)));
	(*(tdcm+1*3+1)) = square_gc(*(quat+0))-square_gc(*(quat+1))+square_gc(*(quat+2))-square_gc(*(quat+3));
	(*(tdcm+1*3+2)) = 2*((*(quat+2))*(*(quat+3)) - (*(quat+0))*(*(quat+1)));
	(*(tdcm+2*3+0)) = 2*((*(quat+1))*(*(quat+3)) - (*(quat+0))*(*(quat+2)));
	(*(tdcm+2*3+1)) = 2*((*(quat+2))*(*(quat+3)) + (*(quat+0))*(*(quat+1)));
	(*(tdcm+2*3+2)) = square_gc(*(quat+0))-square_gc(*(quat+1))-square_gc(*(quat+2))+square_gc(*(quat+3));
}

