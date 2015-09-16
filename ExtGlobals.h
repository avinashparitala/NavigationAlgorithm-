#ifndef _EXTGLOBALS_H_
#define _EXTGLOBALS_H_

#define WORD unsigned short int
extern int Array_SA[32];//For Double buffering

extern int Sprintf(char *outbuf, const char *fmt, ...);

extern double pi, cdr, r0, earth_rate,eccen,g_local,jconst,jconst2, mue,del_t,epsilon;

/*********************************************************************************
		Variables in Leveling_ned.c ECEFNRT.c rsr_utils.c
**********************************************************************************/
extern unsigned int Mission_Data_flag;
extern unsigned int qcnt,velcnt;

extern unsigned int cnt_10ms ,count;
extern volatile unsigned int rcnt;
extern volatile int nav_flag,level_flag,calib_flag,ta_flag;


extern float NavTime,LevTime;
extern float rtime;
extern float four_delt,eight_delt,cdr_delt,cdr_delt_ms;
extern float mdl_si,mdl_phi,mdl_theta;
extern float EulerAngles[3];
extern float p_theta, p_phi, p_si,h_theta, h_phi, h_si;
extern float  period;



extern double p_alp1[3],p_alp2[3],p_alp3[3],p_alp4[3];
extern double FinsLat, FinsLon, FinsAlt;

extern double p_velo[3], h_velo[3];
extern double pure_vel[3],pure_vav[3];

extern double pure_ecef_gravity[3];
extern double omg_dub[3],omega[3];

extern double p_dcm[3][3];

extern double pure_acc_residu[3],pure_gyro_drift[3];
extern double latm,longm;

extern double pure_g_ecef_mag;

extern double r_init,pure_R;

extern double pure_v_old[3],p_velo_20ms[3];
extern double Calib_Angle[3], Calib_Vel[3];

extern double init_dcm_body2ned[3][3],init_dcm_body_ned[3][3];
extern double velo_ref_x[3] , velo_ref_xold[3] , velo_ref_y[3], velo_ref_yold[3];
extern double Ned_omega[3];
extern double p_Ang[3],accum1[6];

extern double Ned_gravity_detic[3];
extern double MasterLat, MasterLon, MasterAlt;
extern double MasterVel[3];

extern double Delta_Angle[3],pure_Delta_Angle[3],Delta_Angle_cal[3];
extern double Delta_Vel[3],pure_Delta_Vel[3],Delta_Vel_cal[3];

extern double delta_velocity_2p5ms_float[3];


extern double p_q_body2ecef[4];
extern double q_ned2body[4];
extern double q_ecef2ned[4];
extern double p_q_body2ned[4];
extern double q_ned2ecef[4];
extern double pure_ecef_pos[3];

extern unsigned int INSDataPre[6];




/****************************************************************************
		Variable in SensorData.c
*****************************************************************************/

extern unsigned char SensDataFIFO[64],SensDataFIFO_LE[64];
extern unsigned int SensorDataValid,SensStatus,SDataInvalidCnt,SensCheck,RNAV_DataError;
extern unsigned int HeadCnt,ChkCnt,RnavOk;
extern double t_delta_theta_2p5ms[3],t_delta_velocity_2p5ms[3];


/****************************************************************************
		Variable in ICD.c
*****************************************************************************/
extern unsigned short int subadd1buf[32],subadd2buf[32],subadd3buf[32],subadd4buf[32],subadd5buf[32],subadd6buf[32],
	 subadd7buf[32],subadd8buf[32],subadd9buf[32],subadd10buf[32],subadd11buf[32],subadd12buf[32],
	 subadd13buf[32],subadd14buf[32],subadd15buf[32],subadd16buf[32],subadd17buf[32],subadd18buf[32],
	 subadd19buf[32],subadd20buf[32],subadd21buf[32],subadd22buf[32],subadd23buf[32],subadd24buf[32],
	 subadd25buf[32],subadd26buf[32],subadd27buf[32],subadd28buf[32],subadd29buf[32],subadd30buf[32],
	 subadd31buf[32];

extern unsigned int pos_test_flag,TestWord,command_flag;
extern volatile union _PPCR16_STATUS_
{
	WORD PPCR16STATUS;
	struct
	{
		unsigned short int reserved : 1;
		unsigned short int GPSValid : 1;
		unsigned short int abort_cond : 1;
		unsigned short int Pos_Valid: 1;
		unsigned short int Pos_Corr : 1;
		unsigned short int GPS_Fix  : 1;
		unsigned short int Maneuver : 1;
		unsigned short int TA_Status: 2;
		unsigned short int MData_Status : 2;
		unsigned short int INS_Hlth : 1;
		unsigned short int Cmd_Num  : 4;
	}INSstatus;

}uINSstatus;

extern volatile union _INS_HEALTH_
{
	WORD HEALTH;
	struct
	{
		unsigned short int reserved :12;
		unsigned short int CalDataChkSum : 1;
		unsigned short int NavAppChkSum : 1;
		unsigned short int BusMemChk : 1;
		unsigned short int ALUChk : 1;
	}INShealth;

}uINShealth;

extern union _CALIB_DATA_
{
	unsigned char CalData[272];
	struct
	{
		double del_t;
		double p_acc_bias[3];
		double acc_sfxp;
		double acc_sfxn;
		double acc_sfyp;
		double acc_sfyn;
		double acc_sfzp;
		double acc_sfzn;
		double acc_mam01;
		double acc_mam02;
		double acc_mam10;
		double acc_mam12;
		double acc_mam20;
		double acc_mam21;
		double acc_nl[3];
		double p_gyro_bias[3];
		double gyro_mam01;
		double gyro_mam02;
		double gyro_mam10;
		double gyro_mam12;
		double gyro_mam20;
		double gyro_mam21;
		double gyro_sf[3];
		double pp_corr;
		double caldate;
		double uno;
	}CALIBDATA;
}uCalibData;


extern WORD APPSW_Checksum[2],CALDATA_Checksum[2];


/****************************************************************************
		Variable in main.c for hardware initialisation
*****************************************************************************/
extern volatile unsigned char *ppcsioaddr_charptr;
extern volatile unsigned short int *ppcsioaddr_shortintptr,*SioStatus,*Sio3Status;
extern volatile unsigned int *ppcsioaddr1;
extern volatile unsigned long  int *ppcsioaddr_longintptr;
extern volatile unsigned short int *ppcsio2addr_shortintptr,*ppcsio1addr_shortintptr;
extern volatile unsigned char *ppcsio1addr_charptr,*ppcsio2addr_charptr,*ppc3sioaddr_charptr;

extern volatile unsigned long int *int_sourceptr,*int_destptr;
extern volatile unsigned short int *addr;
extern volatile unsigned char *Sio2ptr,*Sio3ptr,*Sio4ptr;
extern volatile unsigned char *InhibitReset;
extern volatile unsigned short int *AceMem,*AceReg,*Bus1Ptr,*Bus2Ptr;
extern volatile unsigned short int *DopPtr1,*DopPtr2,*DipPtr1,*DipPtr2,*AdcRdptr;

extern volatile unsigned long int *sp_store,*sp_intret,decrementer_ptr;

extern volatile unsigned int IRQ1Flag,IRQ2Flag;
extern volatile unsigned long int Intr1_Cnt,Intr2_Cnt;


/****************************************************************************
		Variable in TA_Attitude.c
*****************************************************************************/


/*****1553B related********/

//WORD data_rec1[32],data_rec2[32], data_rec3[32],data_rec4[1];

extern WORD data_rec1[32],data_rec2[32],data_rec4[1];
//extern short data_rec3[32];
WORD data_rec3[32];

extern float data_rec11[32],data_rec29[32],file_time[4];
extern unsigned int stkptrA_old,stkptrA,stkoffset;
extern unsigned int rt_bsw,cmd_word;
extern unsigned int timer_cnt;




// Master Params
extern double w_e_ie[3];

extern double sum2MasterVel;

extern float THDG, PITCH, ROLL;
extern float yaw_prev;

extern float MasterTime, prevMasterTime;

// end Master Params

extern unsigned int Rcvchksum;

// KF Related Parameters
extern double ta_states[12][1];
extern double P_ta[12][12];
extern double eye12_ta[12][12];
extern double Q_ta[12][12], R_ta[3][3];
extern double kf_P[12],kf_Q[12],kf_R[3];
extern float t_del_t;

// TA count iterator
extern unsigned int TA_cnt;

// transformations
extern double q_b2ned_M[4];
extern double Cned2ecef[3][3];
extern double Cb2ned_M[3][3], Cb2ned_S[3][3];

extern double Cb2ecef_M[3][3], Cecef2b_M[3][3];
extern double Cb2ecef_S[3][3];

// estimated bias states
extern double ta_biases[3];
extern unsigned int ta_biases_set;


unsigned int Rcvchksum_flag;
unsigned short Time_Elapsed;


// estimated misalignment matrices
extern float known_si, known_theta, known_phi;
extern double CSkew_est[3][3], CSkew_est_T[3][3];
extern double CS2M_K[3][3], CM2S_K[3][3],p_dcm_n[3][3];

extern float mis_angles[3];


// maneuver detectors
extern int          tic ;
extern float        rtime_diff, rtime_prev;
extern float        yaw_diff, yaw_prev;
extern float        min_ta_time ;

// Right Angled Triangles for keeping the pulses intact
extern volatile double  slave_RATriangle_Angles[34][3];
extern volatile double slave_RATriangle_Vels[34][3];
extern int WSZ ;


// WatchTimers
extern volatile float timerWatch;
extern float priorTimerWatch;

// Latch Counter (prior to the now Recv Master Data,
// so many latch_cnts before Master has committed his update)
extern int  latch_cnts;
extern int  no_timer_rsts;
extern int  no_Mdata_pkts;

// Latched Slave Quaternion
extern double LatchSlaveQuat[4];

// thrift angles and velocities;
extern double *sum_thrift_ang_b;
extern double *sum_thrift_vels_b;
extern double thrift_vels_e[3];
extern unsigned int   rcvMdata_rcnt;

extern double thrift_ang_b[3];
extern unsigned int master_data_health;
extern float telm_tdel;


/*
 * TA Variables
 */

extern unsigned int ta_run_flag,inu03_flag,sync_flag;

#endif
