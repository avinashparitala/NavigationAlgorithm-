// ECEFNRT.cpp : hybrid Navigation & Hybrid NAvigation in ECEF; Leveling in NED
// Date : 19-July-2012
// Prog : Saumya


#include<stdio.h>
#include <math.h>

#include "Globals.h"
#include "rsr_utils.h"
#include "Func_Decl.h"


#define GYRX_BIAS 1.0//10.0
#define GYRY_BIAS -0.0//10.0
#define GYRZ_BIAS -0.0//10.0
#define ACCX_BIAS 0.0//500.0
#define ACCY_BIAS 0.0//500.0
#define ACCZ_BIAS 0.0//500.0


/**************************************************************************
FUNCTION NAME: acc_model
DESCRIPTION  : Accelerometer model
INPUT		 : void
OUTPUT		 : void
 **************************************************************************/
void acc_model(void)
{

	double temp1[3],temp2[3],acc_mam[3][3];

	acc_mam[0][1] = uCalibData.CALIBDATA.acc_mam01;
	acc_mam[0][2] = uCalibData.CALIBDATA.acc_mam02;

	acc_mam[1][0] = uCalibData.CALIBDATA.acc_mam10;
	acc_mam[1][2] = uCalibData.CALIBDATA.acc_mam12;

	acc_mam[2][0] = uCalibData.CALIBDATA.acc_mam20;
	acc_mam[2][1] = uCalibData.CALIBDATA.acc_mam21;

	if(Delta_Vel_cal[0] < 0.0)
	{
		acc_mam[0][0] = uCalibData.CALIBDATA.acc_sfxn;
	}
	else
	{
		acc_mam[0][0] = uCalibData.CALIBDATA.acc_sfxp;
	}

	if(Delta_Vel_cal[1] < 0.0)
	{
		acc_mam[1][1] = uCalibData.CALIBDATA.acc_sfyn;
	}
	else
	{
		acc_mam[1][1] = uCalibData.CALIBDATA.acc_sfyp;
	}

	if(Delta_Vel_cal[2] < 0.0)
	{
		acc_mam[2][2] = uCalibData.CALIBDATA.acc_sfzn;
	}
	else
	{
		acc_mam[2][2] = uCalibData.CALIBDATA.acc_sfzp;
	}

	matmul(3, 3, (double *)acc_mam, 3, 1, (double *)Delta_Vel_cal, (double *)temp1);

	temp2[0] = uCalibData.CALIBDATA.acc_nl[0] * Delta_Vel_cal[0]*Delta_Vel_cal[0];
	temp2[1] = uCalibData.CALIBDATA.acc_nl[1] * Delta_Vel_cal[1]*Delta_Vel_cal[1];
	temp2[2] = uCalibData.CALIBDATA.acc_nl[2] * Delta_Vel_cal[2]*Delta_Vel_cal[2];

	Delta_Vel_cal[0] = uCalibData.CALIBDATA.p_acc_bias[0] + temp1[0] + temp2[0];
	Delta_Vel_cal[1] = uCalibData.CALIBDATA.p_acc_bias[1] + temp1[1] + temp2[1];
	Delta_Vel_cal[2] = uCalibData.CALIBDATA.p_acc_bias[2] + temp1[2] + temp2[2];
}

/**************************************************************************
FUNCTION NAME: gyro_model
DESCRIPTION  : GYRO Model
INPUT		 : void
OUTPUT		 : void
 **************************************************************************/
void gyro_model(void)
{

	double temp1[3],gyro_mam[3][3];

	gyro_mam[0][0] = uCalibData.CALIBDATA.gyro_sf[0];
	gyro_mam[0][1] = uCalibData.CALIBDATA.gyro_mam01;
	gyro_mam[0][2] = uCalibData.CALIBDATA.gyro_mam02;

	gyro_mam[1][0] = uCalibData.CALIBDATA.gyro_mam10;
	gyro_mam[1][1] = uCalibData.CALIBDATA.gyro_sf[1];
	gyro_mam[1][2] = uCalibData.CALIBDATA.gyro_mam12;

	gyro_mam[2][0] = uCalibData.CALIBDATA.gyro_mam20;
	gyro_mam[2][1] = uCalibData.CALIBDATA.gyro_mam21;
	gyro_mam[2][2] = uCalibData.CALIBDATA.gyro_sf[2];

	temp1[0] = uCalibData.CALIBDATA.p_gyro_bias[0] + uCalibData.CALIBDATA.gyro_sf[0]*Delta_Angle_cal[0];
	temp1[1] = uCalibData.CALIBDATA.p_gyro_bias[1] + uCalibData.CALIBDATA.gyro_sf[1]*Delta_Angle_cal[1];
	temp1[2] = uCalibData.CALIBDATA.p_gyro_bias[2] + uCalibData.CALIBDATA.gyro_sf[2]*Delta_Angle_cal[2];


	matmul(3, 3, (double *)gyro_mam, 3, 1, (double *)temp1, (double *)Delta_Angle_cal);
}


void put_model(void)
{
	int i;
	/*
	 * sensor to calibration frame
	 * Xc = -Zs
	 * Yc = Ys
	 * Zc = Xs
	 */
	Delta_Angle_cal[0] = -Delta_Angle[2];
	Delta_Angle_cal[1] =  Delta_Angle[1];
	Delta_Angle_cal[2] =  Delta_Angle[0];

	Delta_Vel_cal[0] = -Delta_Vel[2] ;
	Delta_Vel_cal[1] =  Delta_Vel[1] ;
	Delta_Vel_cal[2] =  Delta_Vel[0] ;


	// For Calibration
	if (calib_flag == 1 )
	{
		for(i=0;i<3;i++)
		{
			Calib_Angle[i] += Delta_Angle_cal[i];
			Calib_Vel[i] += Delta_Vel_cal[i];
			//put_s("\tM",2);
		}
	}

	acc_model();
	gyro_model();


	/*
	 * calibration frame to Sensor Frame
	 * Xs = Zc
	 * Ys = Yc
	 * Zs = -Xc
	 */


	/*TODO : DO CALIBRATION FRAME CONVERSION TO INS FRAME*/
	/*
	 * Xins = -Xc
	 * Yins = -Yc
	 * Zins =  Zc
	 */
#if 1
	pure_Delta_Angle[0] =  -Delta_Angle_cal[0];
	pure_Delta_Angle[1] =  -Delta_Angle_cal[1];
	pure_Delta_Angle[2] =   Delta_Angle_cal[2];

	pure_Delta_Vel[0] =  -Delta_Vel_cal[0] ;
	pure_Delta_Vel[1] =  -Delta_Vel_cal[1] ;
	pure_Delta_Vel[2] =   Delta_Vel_cal[2] ;
#endif

}


void SensorDataRead()
{


	int j;

	//ReadPulses_RNAV(); //Reading sensor pulses
	Delta_Angle[0]  = t_delta_theta_2p5ms[0]; //Gyros
	Delta_Angle[1]  = t_delta_theta_2p5ms[1];
	Delta_Angle[2]  = t_delta_theta_2p5ms[2];

	Delta_Vel[0] = t_delta_velocity_2p5ms[0]; //Accelerometers
	Delta_Vel[1] = t_delta_velocity_2p5ms[1];
	Delta_Vel[2] = t_delta_velocity_2p5ms[2];

	put_model(); //Sensor Data Modeling

	/********* Ideal Sensor Equations***********/
	pure_Delta_Angle[0]  = earth_rate * cos(latm)/400 + (GYRX_BIAS * cdr / 3600 / 400);
	pure_Delta_Angle[1]  = 0.0 + (GYRY_BIAS * cdr / 3600 / 400);
	pure_Delta_Angle[2]  = -earth_rate * sin(latm)/400+ (GYRZ_BIAS * cdr / 3600 / 400);

	pure_Delta_Vel[0] = (ACCX_BIAS * 1e-6 * 10/400);
	pure_Delta_Vel[1] = (ACCY_BIAS * 1e-6 * 10/400);
	pure_Delta_Vel[2] = -(pure_g_ecef_mag/400 )+ (ACCZ_BIAS * 1e-6 * 10/400);


	// END OF GYRO MODEL

	if (qcnt == 0)				 //sample 1
	{

		p_alp1[0] = pure_Delta_Angle[0];
		p_alp1[1] = pure_Delta_Angle[1];
		p_alp1[2] = pure_Delta_Angle[2];


	}							 //sample 2
	else if (qcnt == 1)
	{

		p_alp2[0] = pure_Delta_Angle[0];
		p_alp2[1] = pure_Delta_Angle[1];
		p_alp2[2] = pure_Delta_Angle[2];



	}							 //sample 3
	else if (qcnt == 2)
	{
		p_alp3[0] = pure_Delta_Angle[0];
		p_alp3[1] = pure_Delta_Angle[1];
		p_alp3[2] = pure_Delta_Angle[2];

	}							 //sample 4
	else if (qcnt == 3)
	{
		p_alp4[0] = pure_Delta_Angle[0];
		p_alp4[1] = pure_Delta_Angle[1];
		p_alp4[2] = pure_Delta_Angle[2];

	}

	for (j = 0; j < 3; j++)
	{
		p_velo[j] = p_velo[j] + pure_Delta_Vel[j];
		p_Ang[j] = p_Ang[j] + pure_Delta_Angle[j];

	}



	velcnt++;

	if (velcnt == 8)
	{
		velcnt = 0;
		for (j = 0; j < 3; j++)
		{
			p_velo_20ms[j] = p_velo[j];
		}
		accumulator();
	}


}

void accumulator()
{
	int i;

	for (i = 0; i < 3; i++)
	{
		accum1[i] += p_Ang[i];
		accum1[3 + i] += p_velo[i];

		p_velo[i] = 0.0;
		p_Ang[i] = 0.0;

	}
}


void navig(void)
{

	if ((nav_flag == 1) && (level_flag == 0 ))
	{

		SensorDataRead(); //Reading data from Sensors

		qcnt++;					 // no. of samples
		rtime = rcnt * 0.0025;	 // real time
		count++;				 // total no. of cycles
		rcnt++;

		if (qcnt == 4)
		{
			qcnt = 0;

			pure_nav_attup(); //Attitude update routine


			cnt_10ms++;

			if (cnt_10ms == 1)
			{
				quat2dcm((double *)p_q_body2ecef,(double*)p_dcm);
				pure_g_ecef(); //Computation of Gravity vector based on lat and long

			}

			if (cnt_10ms == 2)
			{
				cnt_10ms = 0;

				pure_vel_update(); //Velocity update
				pure_pos_update(); //Position update

			}
		}

	}


}
void NavConstantInit(void)
{
	pi	=	4.0*atan(1.0);
	cdr =	0.017453292519943295769236907684886;
	r0	=	6378135.0;
	earth_rate=	(15.0409*cdr/3600.0);
	eccen=(1.0/298.257);
	g_local=9.78335;
	jconst=0.0010823;
	jconst2=0.0016370;
	mue=3.986005e+14;
	del_t=0.0025;
}

void varinit(void)
{
	int i;
	/*
	 * Resetting all flags
	 */
	Intr1_Cnt=0;
	Intr2_Cnt=0;
	IRQ1Flag = 1;
	IRQ2Flag = 1;

	WSZ = 34;
	TA_cnt =0;


	count = 0;
	qcnt = 0;

	velcnt = 0;
	rtime = 0.0;
	rcnt = 0;
	cnt_10ms = 0;




	latm = MasterLat;
	longm = MasterLon;

	epsilon = 0.0;
	four_delt = 4.0 * del_t;
	eight_delt = 8.0 * del_t;
	cdr_delt = cdr * del_t;
	cdr_delt_ms = cdr_delt / 3600;

	for(i=0;i<32;i++){
		Array_SA[i] = 0;
	}

	for(i=0;i<3;i++)
	{
		velo_ref_y[i] = 0.0;
		velo_ref_yold[i] = 0.0;;
		velo_ref_x[i] = 0.0;
		velo_ref_xold[i] = 0.0;

		pure_vel[i] = 0.0;

		p_velo_20ms[i] = 0.0;
		p_velo[i] = 0.0;

		pure_v_old[i] = 0.0;
		p_Ang[i] = 0.0;

		pure_gyro_drift[i] = 0.0;
		pure_acc_residu[i] = 0.0;

	}

#if 0

	/* these are known misalignment angles between M and S -
	 * Measured w.r.t Master to give DCM from slave to Master.
	 * Beware they are not between slave to NED */
	known_si    =  0.0 * cdr;
	known_theta =  0.0 * cdr;
	known_phi   =  0.0 * cdr;

	euler2dcm_stp(0, 0, 0, (double*)CSkew_est);
	transpose(3, 3, (double*)CSkew_est, (double*)CSkew_est_T);

	euler2dcm_stp(known_si, known_theta, known_phi, (double*)CS2M_K);
	transpose(3, 3, (double*)CS2M_K, (double*)CM2S_K);

	euler2dcm_stp(THDG, PITCH, ROLL, (double*)Cb2ned_M);
	matmul(3, 3, (double*)Cb2ned_M, 3, 3, (double*)CS2M_K, (double*)Cb2ned_S);

	if(ta_flag==1 && nav_flag==1)

	{
		dcm2quat((double*)Cb2ned_S, (double *)p_q_body2ned);

	}

	else if(ta_flag ==0 && level_flag==1)
#endif
	{

		euler2quat_spt(mdl_si,mdl_phi,mdl_theta,(double *)p_q_body2ned);


		p_si = mdl_si;
		p_phi = mdl_phi;
		p_theta = mdl_theta;



	}

	ned2ecef_q(latm, longm,(double*) q_ned2ecef);
	quat_mult((double*)q_ned2ecef,(double*)p_q_body2ned, (double*)p_q_body2ecef);


	/*
	 * Modification after Manjit discussion
	 */
	quat2dcm((double *)p_q_body2ecef,(double*)p_dcm);


	quat2dcm((double *)q_ned2ecef,(double*)p_dcm_n);
	matmul(3,3, (double*)p_dcm_n,3,1,(double*)MasterVel,(double*)pure_vel);




	pure_v_old[0] = pure_vel[0];
	pure_v_old[1] = pure_vel[1];
	pure_v_old[2] = pure_vel[2];

	init(0.0, 0.0, 0.0, p_velo_20ms);

	init(0.0, 0.0, 0.0, p_velo);


	init(0.0,0.0,0.0,pure_gyro_drift);
	init(0.0,0.0,0.0,pure_acc_residu);



	for (i = 0; i < 3; i++)
	{
		p_alp1[i] = 0.0;    p_alp2[i] = 0.0;    p_alp3[i] = 0.0;    p_alp4[i] = 0.0;

	}

	for (i = 0; i < 3; i++)
		Delta_Angle[i] = 0.0;

	for (i = 0; i < 6; i++)
		accum1[i] = 0.0;

	init(0.0, 0.0, earth_rate, omega);	 //earth rate vector ECEF

	//used in levelling
	Ned_omega[0] = earth_rate * cos(latm);
	Ned_omega[1] = 0.0;
	Ned_omega[2] = -earth_rate *sin(latm);

	for (i = 0; i < 3; i++)
		omg_dub[i] = 2.0 * omega[i];

	r_init = r0 * (1.0 - eccen * (sin(latm) * sin(latm)));


	pure_R = r_init + MasterAlt; // altitude;


	lla2ecef(latm,longm,MasterAlt,(double *)pure_ecef_pos); //input is geodetic



	pure_g_ecef();

	/****  for epsilon estimation   ****/

	init(0.0, 0.0, -pure_g_ecef_mag, Ned_gravity_detic);

}								 //end of varinit()




void pure_g_ecef(void)
{
	int i;
	double rad1, rad2;
	double vv[3], vv1[3], r_vec[3];

	rad1 = (r0 * r0) / (pure_R * pure_R);
	rad2 = pure_ecef_pos[2] / pure_R;

	pure_ecef_gravity[0] = -(mue / (r0 * r0)) * (pure_ecef_pos[0] / pure_R) * ((rad1 * (1 + epsilon)) + (jconst2 * (rad1 * rad1)) * (1 - 5 * (rad2 * rad2)));
	pure_ecef_gravity[1] = -(mue / (r0 * r0)) * (pure_ecef_pos[1] / pure_R) * ((rad1 * (1 + epsilon)) + (jconst2 * (rad1 * rad1)) * (1 - 5 * (rad2 * rad2)));
	pure_ecef_gravity[2] = -(mue / (r0 * r0)) * (pure_ecef_pos[2] / pure_R) * ((rad1 * (1 + epsilon)) + (jconst2 * (rad1 * rad1)) * (3 - 5 * (rad2 * rad2)));


	for (i = 0; i < 3; i++)
	{
		r_vec[i] = pure_ecef_pos[i];
		vv[i] = 0.0;
	}

	cross(omega, r_vec, vv);
	cross(omega, vv, vv1);
	matsub(3,1,(double*)pure_ecef_gravity, (double*)vv1, (double*)pure_ecef_gravity);

	pure_g_ecef_mag = sqrt((pure_ecef_gravity[0] * pure_ecef_gravity[0]) + (pure_ecef_gravity[1] * pure_ecef_gravity[1]) + (pure_ecef_gravity[2] * pure_ecef_gravity[2]));

	matmulint(3,1,(double*)pure_ecef_gravity,eight_delt,(double*)pure_ecef_gravity);
}	//end of g_ecef()



void pure_nav_attup(void)
{
	int i;
	double a1 = 1 / 3, a2 = 1 / 3, a3 = 1 / 3, b1 = 0.5, b2 = 0.5;
	double ap1[3], ap2[3], ap3[3], ap4[3], ap5[3], ap6[3];

	double k1=54.0/105.0,k2=92.0/105.0,k3=214.0/105.0;
	double temp_q[4],qt[4],qtinv[4];
	double q_tmp[4];
	double s, beta2, sum, beta[3];

	for (i = 0; i < 3; i++)		 // reason for floating point earth_rateor
	{
		ap1[i] = 0.0; ap2[i] = 0.0; ap3[i] = 0.0;
		ap4[i] = 0.0; ap5[i] = 0.0; ap6[i] = 0.0;
		beta[i] = 0.0;
	}
	for (i = 0; i < 4; i++)
	{
		temp_q[i] = 0.0;
		q_tmp[i] = 0.0;
	}

	cross(p_alp1,p_alp4,ap1);
	cross(p_alp1,p_alp3,ap2);
	cross(p_alp2,p_alp3,ap3);
	cross(p_alp3,p_alp4,ap4);
	cross(p_alp2,p_alp4,ap5);
	cross(p_alp1,p_alp2,ap6);

	for(i=0;i<3;i++)
	{
		beta[i] = p_alp1[i]+p_alp2[i]+p_alp3[i]+p_alp4[i]+k1*ap1[i]+k2*(b1*ap2[i]+b2*ap5[i])+k3*(a1*ap6[i]+a2*ap3[i]+a3*ap4[i]);
	}
	beta2 = (pow(beta[0],2)) + (pow(beta[1],2)) + (pow(beta[2],2));
	sum   = 0.5 - beta2/48.0;

	temp_q[0] = 1.0-beta2/8.0+pow(beta2,2)/480.0;
	temp_q[1] = sum*beta[0];
	temp_q[2] = sum*beta[1];
	temp_q[3] = sum*beta[2];

	quat_mult((double *)p_q_body2ecef,(double *)temp_q,(double *)q_tmp);

	for(i=0;i< 4;i++)
	{
		p_q_body2ecef[i] = q_tmp[i];
	}

	//earth rate corr
	for (i = 0; i < 3; i++)
		beta[i] = 4 * del_t * omega[i];

	beta2 = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
	s = 0.5 - beta2 / 48.0;

	qt[0] = 1.0 - beta2 / 8.0 + (beta2 * beta2) / 480.0;
	qt[1] = s * beta[0];
	qt[2] = s * beta[1];
	qt[3] = s * beta[2];

	quat_inv(qt,qtinv);
	quat_mult(qtinv, p_q_body2ecef, q_tmp);

	for (i = 0; i < 4; i++)
		p_q_body2ecef[i] = q_tmp[i];


	quat_norm((double *)p_q_body2ecef);

}								 //end of nav_attup


// POSTION UPDATE
void pure_pos_update(void)
{
	int  i;
	double temp1[3];

	for (i = 0; i < 3; i++)
	{
		temp1[i] = pure_vav[i] * eight_delt;

		pure_ecef_pos[i] = pure_ecef_pos[i]+temp1[i];

	}

	pure_R = sqrt((pure_ecef_pos[0] * pure_ecef_pos[0]) + (pure_ecef_pos[1] * pure_ecef_pos[1]) + (pure_ecef_pos[2] * pure_ecef_pos[2]));
	get_Pure_Euler_angles();


	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",rtime,pure_ecef_pos[0],pure_ecef_pos[1],pure_ecef_pos[2],p_si,p_theta,p_phi);
}



void get_Pure_Euler_angles(void)
{
	double tdcm[3][3];
	double Cbody2ned[3][3],Cecef2ned[3][3];

	// free willy
	FinsLat = asin(pure_ecef_pos[2]/pure_R);

	FinsLon = atan2(pure_ecef_pos[1],pure_ecef_pos[0]);

	ecef2ned_dcm(FinsLat, FinsLon, (double*)Cecef2ned);

	quat2dcm((double *)p_q_body2ecef, (double*)tdcm);

	matmul(3, 3, (double*)Cecef2ned, 3, 3, (double*)tdcm, (double*)Cbody2ned);


	dcm2euler_stp((double*)Cbody2ned, &p_si,&p_theta,&p_phi); //soumi changed on 1st May..
	//dcm2euler_spt((double*)Cbody2ned, &p_si,&p_phi,&p_theta);

	p_si = p_si / cdr;
	p_theta = p_theta / cdr;
	p_phi = p_phi / cdr;

	//kgdb_init();
//	breakpoint();

}

void pure_vel_update(void)
{
	int i;
	double temp[3];

	matmul(3, 3, (double *)p_dcm, 3, 1, (double *)p_velo_20ms, (double *)pure_vel);

	//earth rate correction
	cross(omg_dub, pure_v_old, temp);
	matmulint(3,1,(double *)temp, eight_delt, (double *)temp);
	matsub(3,1,(double *)pure_vel, (double *)temp, (double *)pure_vel);
	matadd(3,1,(double *)pure_vel, (double *)pure_ecef_gravity, (double *)pure_vel);

	matadd(3,1,(double *)pure_vel, (double *)pure_v_old, (double *)pure_vel);

	matadd(3,1,(double *)pure_vel, (double *)pure_v_old, (double *)pure_vav);
	matmulint(3,1,(double *)pure_vav, 0.5, (double *)pure_vav);

	for (i = 0; i < 3; i++)
		pure_v_old[i] = pure_vel[i];


}


int main(void)
{
	NavConstantInit();  //Navigation constants initialisation

	/* Data Required to be loaded from GUI*/
	MasterLat = 17.0*cdr;
	MasterLon=78.0*cdr;
	MasterAlt= 500.0;
	NavTime = 300.0;
	mdl_si=0.0;
	mdl_phi=0.0;
	mdl_theta = 0.0;

/* Variables initialisation*/
	varinit();

// These flags required for executing navigation run
	nav_flag = 1;
	level_flag = 0;

//Initialised to computed euler angles
	FinsLat = MasterLat;
	FinsLon = MasterLon;

	while (rcnt < NavTime * 400)
	{
		navig(); //Navigation Module

	}

return 0;

}






