#ifndef _GLOBALS_H_
#define _GLOBALS_H_	

#define WORD unsigned short int
int Array_SA[32];//For Double buffering

extern int Sprintf(char *outbuf, const char *fmt, ...);

double pi, cdr, r0, earth_rate,eccen,g_local,jconst,jconst2, mue,del_t,epsilon;

/*********************************************************************************
		Variables in Leveling_ned.c ECEFNRT.c rsr_utils.c
**********************************************************************************/
unsigned int Mission_Data_flag;
unsigned int qcnt,velcnt;

unsigned int cnt_10ms ,count;
volatile unsigned int rcnt;
volatile int nav_flag,level_flag,calib_flag,ta_flag;


float NavTime,LevTime;
float rtime;
float four_delt,eight_delt,cdr_delt,cdr_delt_ms;
float mdl_si,mdl_phi,mdl_theta;
float EulerAngles[3];
float p_theta, p_phi, p_si,h_theta, h_phi, h_si;
float  period;



double p_alp1[3],p_alp2[3],p_alp3[3],p_alp4[3];
double FinsLat, FinsLon, FinsAlt;

double p_velo[3], h_velo[3];
double pure_vel[3],pure_vav[3];

double pure_ecef_gravity[3];
double omg_dub[3],omega[3];

double p_dcm[3][3];

double pure_acc_residu[3],pure_gyro_drift[3];
double latm,longm;

double pure_g_ecef_mag;

double r_init,pure_R;

double pure_v_old[3],p_velo_20ms[3];
double Calib_Angle[3], Calib_Vel[3];

double init_dcm_body2ned[3][3],init_dcm_body_ned[3][3];
double velo_ref_x[3] , velo_ref_xold[3] , velo_ref_y[3], velo_ref_yold[3];
double Ned_omega[3];
double p_Ang[3],accum1[6];

double Ned_gravity_detic[3];
double MasterLat, MasterLon, MasterAlt;
double MasterVel[3];

double Delta_Angle[3],pure_Delta_Angle[3],Delta_Angle_cal[3];
double Delta_Vel[3],pure_Delta_Vel[3],Delta_Vel_cal[3];

double delta_velocity_2p5ms_float[3];


double p_q_body2ecef[4];
double q_ned2body[4];
double q_ecef2ned[4];
double p_q_body2ned[4];
double q_ned2ecef[4];
double pure_ecef_pos[3];

unsigned int INSDataPre[6];




/****************************************************************************
		Variable in SensorData.c
*****************************************************************************/

unsigned char SensDataFIFO[64],SensDataFIFO_LE[64];
unsigned int SensorDataValid,SensStatus,SDataInvalidCnt,SensCheck,RNAV_DataError;
unsigned int HeadCnt,ChkCnt,RnavOk;
double t_delta_theta_2p5ms[3],t_delta_velocity_2p5ms[3];


/****************************************************************************
		Variable in ICD.c
*****************************************************************************/
WORD subadd1buf[32],subadd2buf[32],subadd3buf[32],subadd4buf[32],subadd5buf[32],subadd6buf[32],
	 subadd7buf[32],subadd8buf[32],subadd9buf[32],subadd10buf[32],subadd11buf[32],subadd12buf[32],
	 subadd13buf[32],subadd14buf[32],subadd15buf[32],subadd16buf[32],subadd17buf[32],subadd18buf[32],
	 subadd19buf[32],subadd20buf[32],subadd21buf[32],subadd22buf[32],subadd23buf[32],subadd24buf[32],
	 subadd25buf[32],subadd26buf[32],subadd27buf[32],subadd28buf[32],subadd29buf[32],subadd30buf[32],
	 subadd31buf[32];
	
unsigned int pos_test_flag,TestWord,command_flag;
volatile union _PPCR16_STATUS_
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

volatile union _INS_HEALTH_
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

union _CALIB_DATA_
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


WORD APPSW_Checksum[2],CALDATA_Checksum[2];


/****************************************************************************
		Variable in main.c for hardware initialisation
*****************************************************************************/
volatile unsigned char *ppcsioaddr_charptr;
volatile unsigned short int *ppcsioaddr_shortintptr,*SioStatus,*Sio3Status;
volatile unsigned int *ppcsioaddr1;
volatile unsigned long  int *ppcsioaddr_longintptr;
volatile unsigned short int *ppcsio2addr_shortintptr,*ppcsio1addr_shortintptr;
volatile unsigned char *ppcsio1addr_charptr,*ppcsio2addr_charptr,*ppc3sioaddr_charptr;
	
volatile unsigned long int *int_sourceptr,*int_destptr;
volatile unsigned short int *addr;
volatile unsigned char *Sio2ptr,*Sio3ptr,*Sio4ptr;
volatile unsigned char *InhibitReset;
volatile unsigned short int *AceMem,*AceReg,*Bus1Ptr,*Bus2Ptr;
volatile unsigned short int *DopPtr1,*DopPtr2,*DipPtr1,*DipPtr2,*AdcRdptr;

volatile unsigned long int *sp_store,*sp_intret,decrementer_ptr;

volatile unsigned int IRQ1Flag,IRQ2Flag;
volatile unsigned long int Intr1_Cnt,Intr2_Cnt;


/****************************************************************************
		Variable in TA_Attitude.c
*****************************************************************************/


/*****1553B related********/

//WORD data_rec1[32],data_rec2[32], data_rec3[32],data_rec4[1];

WORD data_rec1[32],data_rec2[32],data_rec4[1];
//short data_rec3[32];
WORD data_rec3[32];
float data_rec11[32],data_rec29[32],file_time[4];
unsigned int stkptrA_old,stkptrA,stkoffset;
unsigned int rt_bsw,cmd_word;
unsigned int timer_cnt;




// Master Params
double w_e_ie[3];

double sum2MasterVel;

float THDG, PITCH, ROLL;
float yaw_prev;

float MasterTime, prevMasterTime;

// end Master Params

unsigned int Rcvchksum;

// KF Related Parameters
double ta_states[12][1];
double P_ta[12][12];
double eye12_ta[12][12];
double Q_ta[12][12], R_ta[3][3];
double kf_P[12],kf_Q[12],kf_R[3];
float t_del_t;

// TA count iterator
unsigned int TA_cnt;

// transformations
double q_b2ned_M[4];
double Cned2ecef[3][3];
double Cb2ned_M[3][3], Cb2ned_S[3][3];

double Cb2ecef_M[3][3], Cecef2b_M[3][3];
double Cb2ecef_S[3][3];

// estimated bias states
double ta_biases[3];
unsigned int ta_biases_set;


unsigned int Rcvchksum_flag;
unsigned short Time_Elapsed;


// estimated misalignment matrices
float known_si, known_theta, known_phi;
double CSkew_est[3][3], CSkew_est_T[3][3];
double CS2M_K[3][3], CM2S_K[3][3],p_dcm_n[3][3];

float mis_angles[3];


// maneuver detectors
int          tic ;
float        rtime_diff, rtime_prev;
float        yaw_diff, yaw_prev;
float        min_ta_time ;

// Right Angled Triangles for keeping the pulses intact
volatile double  slave_RATriangle_Angles[34][3];
volatile double slave_RATriangle_Vels[34][3];
int WSZ ;


// WatchTimers
volatile float timerWatch;
float priorTimerWatch;

// Latch Counter (prior to the now Recv Master Data,
// so many latch_cnts before Master has committed his update)
int  latch_cnts;
int  no_timer_rsts;
int  no_Mdata_pkts;

// Latched Slave Quaternion
double LatchSlaveQuat[4];

// thrift angles and velocities;
double *sum_thrift_ang_b;
double *sum_thrift_vels_b;
double thrift_vels_e[3];
unsigned int   rcvMdata_rcnt;

double thrift_ang_b[3];
unsigned int master_data_health;
float telm_tdel;


/*
 * TA Variables
 */

unsigned int ta_run_flag,inu03_flag,sync_flag;



#define w_e              7.292115e-05
#define Rer              6378137.0
#define ecc              0.0818191908426



#endif
