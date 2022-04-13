#ifndef MAGcal_h_
#define MAGcal_h_

#include <math.h>
#include <string.h>
#include <stdint.h>
#include "esp32-hal-log.h"


#define TIMEOUT_MSEC 14

#define MAGBUFFSIZE 400 // Freescale's lib needs at least 392

// compile time constants that are private to this file
#define CORRUPTMATRIX 0.001F      // column vector modulus limit for rotation matrix

// vector components



#ifdef __cplusplus
extern "C"{
#endif

typedef struct {
	float x;
	float y;
	float z;
	//int valid;
} Point_t;

typedef struct {
	float q0; // w
	float q1; // x
	float q2; // y
	float q3; // z
} Quaternion_t;
extern Quaternion_t current_orientation;



void apply_calibration(int16_t rawx, int16_t rawy, int16_t rawz, Point_t *out);
int MagCal_Run(void);



// magnetic calibration & buffer structure
typedef struct {
    float V[3];                  // current hard iron offset x, y, z, (uT)
    float invW[3][3];            // current inverse soft iron matrix
    float B;                     // current geomagnetic field magnitude (uT)
    float FourBsq;               // current 4*B*B (uT^2)
    float FitError;              // current fit error %
    float FitErrorAge;           // current fit error % (grows automatically with age)
    float trV[3];                // trial value of hard iron offset z, y, z (uT)
    float trinvW[3][3];          // trial inverse soft iron matrix size
    float trB;                   // trial value of geomagnetic field magnitude in uT
    float trFitErrorpc;          // trial value of fit error %
    float A[3][3];               // ellipsoid matrix A
    float invA[3][3];            // inverse of ellipsoid matrix A
    float matA[10][10];          // scratch 10x10 matrix used by calibration algorithms
    float matB[10][10];          // scratch 10x10 matrix used by calibration algorithms
    float vecA[10];              // scratch 10x1 vector used by calibration algorithms
    float vecB[4];               // scratch 4x1 vector used by calibration algorithms
    int8_t ValidMagCal;          // integer value 0, 4, 7, 10 denoting both valid calibration and solver used
    float BpFast[3][MAGBUFFSIZE];   // uncalibrated magnetometer readings
    int8_t  valid[MAGBUFFSIZE];        // 1=has data, 0=empty slot
    int16_t MagBufferCount;           // number of magnetometer readings
} MagCalibration_t;


MagCalibration_t magcal;


void f3x3matrixAeqI(float A[][3]);
void fmatrixAeqI(float *A[], int16_t rc);
void f3x3matrixAeqScalar(float A[][3], float Scalar);
void f3x3matrixAeqInvSymB(float A[][3], float B[][3]);
void f3x3matrixAeqAxScalar(float A[][3], float Scalar);
void f3x3matrixAeqMinusA(float A[][3]);
float f3x3matrixDetA(float A[][3]);
void eigencompute(float A[][10], float eigval[], float eigvec[][10], int8_t n);
void fmatrixAeqInvA(float *A[], int8_t iColInd[], int8_t iRowInd[], int8_t iPivot[], int8_t isize);
void fmatrixAeqRenormRotA(float A[][3]);

void fUpdateCalibration4INV(MagCalibration_t *MagCal);
void fUpdateCalibration7EIG(MagCalibration_t *MagCal);
void fUpdateCalibration10EIG(MagCalibration_t *MagCal);


#define SENSORFS 100
#define OVERSAMPLE_RATIO 4


// accelerometer sensor structure definition
#define G_PER_COUNT 0.0001220703125F  // = 1/8192
typedef struct
{
	float Gp[3];           // slow (typically 25Hz) averaged readings (g)
	float GpFast[3];       // fast (typically 200Hz) readings (g)
} AccelSensor_t;

// magnetometer sensor structure definition
#define UT_PER_COUNT 0.1F
typedef struct
{
	float Bc[3];           // slow (typically 25Hz) averaged calibrated readings (uT)
	float BcFast[3];       // fast (typically 200Hz) calibrated readings (uT)
} MagSensor_t;

// gyro sensor structure definition
#define DEG_PER_SEC_PER_COUNT 0.0625F  // = 1/16
typedef struct
{
	float Yp[3];                           // raw gyro sensor output (deg/s)
	float YpFast[OVERSAMPLE_RATIO][3];     // fast (typically 200Hz) readings
} GyroSensor_t;


#ifdef __cplusplus
} // extern "C"
#endif

#endif
