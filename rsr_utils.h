#ifndef _RSR_UTILS_H_
#define _RSR_UTILS_H_

// standard constants and constant expressions
#define R0							6378137.0
#define r							6356752.3142
//#define f							3.352810671830990e-003
#define flat					    3.352810671830990e-003//f							/*** (1.0/298.257) ***/
#define	e							0.0818191908426
//#define eccen					e
#define w_e						7.292115e-05
//#define jconst					0.0010823				/*** oblateness constant  ***/
//#define	mue						3.986005e+14			/*** univ gravity const * mass of earth ***/
//#define jconst2					1.637e-03				/*** for eci gravity model ***/
#define gravity_hayes		9.809658054
//#define epsilon					0
//#define pi							3.141592653589793
//#define cdr                         1.745329251994330e-002

#define square(val) val*val

#define M(lat)			R0 * (1.0-square(e)) / ( ( 1.0 - square(e)*square(sin(lat)) ) * sqrt( ( 1.0 - square(e)*square(sin(lat)) ) ) )
#define N(lat)			R0 /  sqrt( ( 1.0 - square(e)*square(sin(lat)) ) ) 

// geographic utilities

// position utilities
void lla2ecef(double lat, double longm, double alt, double ecef_pos[3]);
void ecef2lla(double *ecef_pos, double *latm, double *longm, double *altm);

// ned related routines
void ecef2ned_dcm(double latm, double longm, double *Cecef2ned);
void ned2ecef_dcm(double latm, double longm, double *Cned2ecef);
void ned2ecef_q(double latm, double longm, double *qned2ecef);
void ecef2ned_q(double latm, double longm, double *qecef2ned);
void ecef2ned(double latm, double longm, double *vin, double *vout);
void ned2ecef(double latm, double longm, double *vin, double *vout);
void getNEDOmega(float latd, float *nedOmega);
void getNEDGravity(float latd, float alt, float *nedGravity);

// nwv related routines
void ecef2nwv_q(double latm, double longm, double *qecef2nwv);
void nwv2ecef_q(double latm, double longm, double *qnwv2ecef);
void ecef2nwv_dcm(double latm, double longm, double *Cecef2nwv);
void nwv2ecef_dcm(double latm, double longm, double *Cnwv2ecef);
void ecef2nwv(double latm, double longm, double *vin, double *vout);
void getNWVOmega(float latd, float *nwvOmega);
void getNWVGravity(float latd, float alt, float *nwvGravity);

// between ned and nwv
void ned2nwv_dcm(double *dcm);
void nwv2ned_dcm(double *dcm);
void ned2nwv_q(double *q);
void nwv2ned_q(double *q);
void ned2nwv(double *vned, double *vnwv);
void nwv2ned(double *vnwv, double *vned);

// just pure geographic ones
float GeoDetic_to_GeoCentric(float Ld, float alt);
float GeoCentric_to_GeoDetic(float Lc);

// gravity models
void		gr_britting_ned(float lat, float alt, float *Ned_gravity);
double	    gravity_ecef(double *p_ecef_pos, double *p_ecef_gravity);

// Quaternion utilities
void		quat_inv(double *q_in, double *q_inv);
void		quat_transform(double *quat, double *in_vec, double *out_vec);
void		quat_mult(double *q1, double *q2, double *q3);
void		quat_normalize(double *quat_in);
void	    quat_norm(double *quat);
void	    quat_mod(double *quat_in);

// Rotation transformation utilities
void dcm2quat(double *tdcm, double *quat);
void quat2dcm(double *quat, double *tdcm);
void euler2dcm_spt(float si,float phi,float theta,double *mm);  /* si,phi,theta(y,r,p) */ /* rci sequence */
void dcm2euler_spt(double *mm, float *si, float *phi, float *theta);
void euler2dcm_stp(float si,float theta,float phi,double *mm);  /* si,theta,phi(y,p,r) */ /* aerospace sequence */
void dcm2euler_stp(double *mm, float *si, float *theta, float *phi);
void quat2rot(double *quat, float *rot);

void euler2quat_spt(float si,float phi,float theta,double *quat);  
void euler2quat_stp(float si, float theta, float phi, double *quat);
void quat2euler_spt(double *quat,float *si, float *phi, float *theta);   
void quat2euler_stp(double *quat,float *si, float *theta, float *phi);
void rot2dcm(float *rot, float *dcm);
void dcm2rot(float *dcm, float *rot);

// Vector utilities
void		cross(double *v4,double *v5,double *v6);
double	    norm(int size, double *vec);
float	    max(int size, float *vec);
void		init(double a, double b, double c, double *v3);



// Matrix utilities
void		init_m(int size, double *mat, ...);
void		matmul(int m1rows, int m1cols, double  *M1, \
		int m2rows, int m2cols, double  *M2, \
		double *M3);
//void		invert2(int size, float *temp_M, float *M_inv);
//void		invert3(int size, float *temp_M, float *M_inv);
void		invert_splu(int size, float *temp_M,float *M_inv);
void		matmulint(int m1rows, int m1cols, double *M1,double scalar, double *M3);
void		matclone(int rows,int cols,float *M1,float *M2);
void		matadd(int rows, int cols,double *M1,double *M2,double *M3);
void		matsub(int rows, int cols,double *M1,double *M2,double *M3);
void		transpose(int rows, int cols, double *M1, double *M2);
void		to_cross_form(double *, double *);
float	    trace(int size, float *mat);
void		print_matrix(int rows, int cols, float *mat, char *name);



// miscellaneous
unsigned int  ComputeChecksum( unsigned int *data, unsigned int cnt);
void				Delay(unsigned int count);

#endif
