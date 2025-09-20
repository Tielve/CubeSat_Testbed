#ifndef __KALMAN_H__
#define __KALMAN_H__
#include <math.h>
#include <string.h>

typedef struct{
    float q[4];   // attitude quaternion (w,x,y,z), body←inertial
    float b[3];   // gyro bias [rad/s]
    float P[6][6]; // covariance
    float noise_g;   // gyro noise [rad/s/√Hz]
    float bias_g;  // gyro bias RW [rad/s^2/√Hz]
    float noise_a; // accel meas dir noise (normalized) ~0.03..0.1
    float noise_m; // mag meas dir noise (normalized) ~0.03..0.1
    float dt; // time
} Kalman;

void init_kalman(Kalman *k);
void vec_cross(float v1[3], float v2[3], float res[3]);
void skew(float v1[3], float res[3][3]);
void transpose3(float A[3][3], float res[3][3]);
void transpose6(float A[6][6], float res[6][6]);
void inverse6(float A[6][6], float res[6][6]);
void kalman_predict(Kalman *k, float g[3], float P_est[6][6]);
void kalman_correction(Kalman *k, float a[3], float m[3], float P_est[6][6]);

#endif
