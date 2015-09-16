void ECEFNRT(void);
void NavConstantInit(void);
void levelling(void);
void Mission_DataLoad(void);
void Mission_DataDump(void);
void timer1ISR(void);
void SensorDataRead(void);
void Write_Inc_Ang_Vel(void);
void varinit(void);
void put_model(void);
//function declaration in Pnav

void accumulator(void);

void pure_nav_attup(void);
void pure_pos_update(void);
void pure_vel_update(void);
void pure_g_ecef(void);
void get_Pure_Euler_angles(void);

void hybrid_nav_attup(void);
void hybrid_pos_update(void);
void hybrid_vel_update(void);
void hybrid_g_ecef(void);
void get_Hybrid_Euler_angles(void);

void insitu_drift(void);
void IRQ0_ISR(int);
void Navigation_post_lev(void);
void navigation_post_ta(void);
void NAV_cleardata(void);


//leveling
void LevelQuaternionUpdate(void);
void LevelQuatSlew(void);
void LevelNedVelocityResolve(void);
void compute_dcm(double jj[4], double m[3][3]);
void compute_euler_angles(double mm[3][3], double angles[3]);
void vect_cross_mult( double v4[3],double v5[3],double *v6);
void init_vect(double a,double b,double c,double *v3);
void vectintmul(double k,double v7[3], double *v8);
void mat_trans( double mk[3][3], double gk[3][3]);
void vectadd(double v8[3],double v9[3],double *vp1);



/********************************************************************************
	Functions in main.cpp
********************************************************************************/
#define WORD unsigned short int
void clear_data_block(void);
void Register_Allocation(void);
void delay(unsigned int );
void INT1_disable(void);
void INT1_unmask(void);
void INT2_unmask(void);
void Sio2IpcoreInit(void);
void inittimer0(void);
void siu_int_controller_init(void);
void InitRt(void);
void WriteData(WORD ,WORD data_buff[],WORD);
void ReadData(WORD ,WORD data_buff[],WORD );
void Iwait(void);
void timer1Handler(void);
void Iwait2(void);
void timer2Handler(void);
void Delay(unsigned int );
void put_char(char ABYTE);   
void put_s(char *str,WORD length);
void FlashEpromTest(void);

extern WORD MemRead(WORD Offset);
extern void SelectModule(unsigned char);

/********************************************************************************
	Functions in SensorData.cpp
********************************************************************************/
void BigEndiantoLittleEndian(void);
void ReadPulses_RNAV(void);

/********************************************************************************
	Functions in ICD.cpp
********************************************************************************/
void write_debug_data(void);
void write_debug_data2(void);
void write_nav_angles(void);
void Write_lev_data(void);
void write_status(void);
void write_pos_vel(void);
void Write_KF_Corrections(void);
void Send_TA_data_1(void);

void Write_GPS_data_SA12(void);
void Write_GPS_data_SA13(void);
void Write_GPS_data_SA14(void);
void Write_GPS_data_SA15(void);
void Write_GPS_data_SA16(void);
void Write_GPS_data_SA17(void);
void Write_GPS_data_SA18(void);
void Write_GPS_data_SA19(void);
void Write_GPS_data_SA20(void);
void Write_GPS_data_SA21(void);
void Write_GPS_data_SA22(void);
void Write_GPS_data_SA23(void);
void Write_GPS_data_SA24(void);
void Write_GPS_data_SA25(void);
void Write_GPS_data_SA26(void);

void position_test(void);
void PBIT(void);
void ResetFIFO(void);

/********************************************************************************
	Functions in Utils1553.cpp
********************************************************************************/

void WriteDouble_PPC(int noe, double *a, WORD *array1);
void ReadDouble_PPC(int noe, WORD *data, double *op);
void WriteFloat_PPC(int noe, float  *a,WORD *Temp_Buff);
void ReadFloat_PPC(int noe, WORD *data, float *op);
void ReadInt_PPC(int noe, WORD *data, long int *op);
void WriteInt_PPC(int noe, long int *a,WORD *Temp_Buff);

/********************************************************************************
	Functions in GPSAcq.c
********************************************************************************/
void GPSDATAAcq(void);
void read_gps_data(void);
void read_gps_dataSA1(void);
void message(int ,int ,int ,int ,unsigned short int *);
void FrameGPSMessage(void);
void SchedlGPSData(int mess_no);
void SchedlAllGPSMsgs(void);
void Set_GPS_DATA_AVBL(void);
//checking last three sample
void setGpsValidity(void);
void ByteXor(WORD x);
void gDataCheckSum1(void);
void gDataCheckSum2(void);
void gDataCheckSum3(void);
void AssignFix_Dop(void);
void CmptWeaponSysTime(void);
void Update3DFixLoss(void);

void GPSReadData(void);

/********************************************************************************
	Functions in HybridNav.c
********************************************************************************/
void Hybrid_navigation(void);
void att_corr(void);
void pos_corr(void);
void vel_corr(void);
void Assign_GPS_Data(void);
void AssignG3OMData(void);
void Init_KF_Var(void);
void LimitKFCorr(void);
void Init_first_KF_var(void);

void nkprpc(int *nkrini,int nkuc,int nkipn,float nkcps);
void nkfilm( double nktlos[28][3],double nktplos[28][3]);
void nkprpm( double nkprdt,double nktlos[28][3],double nkpevr[28],int * nkmedt,long int* nkdrtm,long int* nkprtm);
void nkcntl_output(void);


/********************************************************************************
	Functions in GC.c
********************************************************************************/


/* function definitions */
void rtd(double *rrr, double *ddd);
void one_shot_L_vectors(int iter, double norm_chi, double lat_detic);
void coarse_gc(void);
void fine_gc(void);

/* utility functions */
void cross_gc( double *v4,double *v5,double *v6);
void get_euler_stp(double *mm, double *angles);
void get_euler_spt(double *, double *);
void dcm_to_quat(double *tdcm, double *quat);
void matmul_double(int m1rows, int m1cols, double *M1,int m2rows, int m2cols, double *M2, double *M3);
void transpose_double(int rows, int cols, double *M1, double *M2);
void matadd_double(int rows, int cols, double *M1,double *M2,double *M3);
void matmulint_double(int m1rows, int m1cols, double *M1, double scalar, double *M3);
void matmulvec(double *M1, double *M2, double *M3);
void rot2quat(double *rot, double *quat);
void quatmult(double *q1, double *q2, double *q3);
void quat_transform(double *quat, double *in_vec, double *out_vec);
void invert3(int size, double *temp_M, double *M_inv);

void quat_to_dcm(double *quat, double *tdcm);

void invert2(int size, double *temp_M, double *M_inv);

void init_gc(void);

void euler2dcm_spt_double(float si,float phi,float theta,double *mm);



/********************************************************************************
	Functions in TA_Attitude.c
********************************************************************************/
// KF Functions
void form_stm_TA(float, double *);
void form_H_TA(double *);
void running_kf(void);
void Receive_Data(void);
void Reset_Timer_Now(void);
void Take_MasterData(void);
void TA_Attitude(void);
float get_filter_update_interval(void);
void form_RATriangle(void);
void Transfer_Alignment(void);
int get_latch_cnts(void);

// Utilities
/*
void quatinv(double *, double *);
void dcm2quat(double *, double *);
void quatrotate(float *, float *, float *);
void to_cross_form(double *, double *);
void matmul(int, int, double *, int, int, double *, double *);
void matmulint(int, int, double *, double, double *);
void matadd(int, int, double *, double *, double *);
void transpose(int, int, double *, double *);
void invert3(int, double *, double *);
void ned2ecef_dcm(double, double, double *);
void euler2dcm_stp(float si, float theta, float phi, double *mm);
void quat2euler(double *, float *, float *, float *);
void lla2ecef(double lat, double longm, double alt, double *ecef_pos);
*/
void init_KF_TA(int, double *, int, double *, int, double *, int bias_reset);
float BNR1toFloat(short, float, int);
float BC1toFloat(short, float, int);
float BC2toFloat(long int, float, int);
void Filter_Para_load(void);
void poll_bc_rt_status(void);
void Receive_Data(void);
void Reset_Timer_Now(void);
void formsie_LatchSlaveQuat(volatile int cntr);
void is_nav(void);
void Take_MasterData(void);
void record_biases(void);
void reset_att_states(void);
void reset_bias_states(void);
void record_CSkew(void);
void TA_Attitude(void);
void form_measurement(double *Z1);
float get_filter_update_interval(void);
void running_kf(void);
void form_stm_TA(float t_del_t, double *TA_STM);
void form_RATriangle(void);
void Transfer_Alignment(void);
void form_H_TA(double *H);

void ReadFullSim(void);

