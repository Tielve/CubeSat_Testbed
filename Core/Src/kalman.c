#include "kalman.h"

void init_kalman(Kalman *k) {
  // Start with identity
  memset(k->q, 0, sizeof(k->q));
  k->q[0] = 1;
  k->q[1] = 0;
  k->q[2] = 0;
  k->q[3] = 0;

  // Gyro bias
  memset(k->b, 0, sizeof(k->b));
  k->b[0] = .005f;
  k->b[1] = .005f;
  k->b[2] = .005f;

  // Populate P with covariances
  memset(k->P, 0, sizeof(k->P));
  for (int i = 0; i < 3; i++) {
    k->P[i][i] = .1745f * .1745f;
    k->P[i + 3][i + 3] = .00873f * .00873f;
  }

  // Set noise and dt
  k->noise_g = 0.02f;
  k->bias_g = 0.002f;
  k->noise_a = 0.02f;
  k->noise_m = 0.05f;
  k->dt = 0.002f;
}
void vec_cross(float v1[3], float v2[3], float res[3]) {
  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

float vec_dot(float v1[3], float v2[3]) {
  return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

float magnitude(float v1[3]) { return sqrt(vec_dot(v1, v1)); }

float magnitudeq(float v1[4]) {
  return sqrt((v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3] * v1[3]));
}

void skew(float v1[3], float A[3][3]) {
  float temp[3][3] = {
      {0, -v1[2], v1[1]}, {v1[2], 0, -v1[0]}, {-v1[1], v1[0], 0}};
  A = temp;
}

void transpose3(float A[3][3], float C[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      C[i][j] = A[j][i];
}

void transpose6(float A[6][6], float C[6][6]) {
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      C[i][j] = A[j][i];
}

void mat_mul3(float A[3][3], float B[3][3], float C[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      float sum = 0;
      for (int k = 0; k < 3; k++) {
        sum += A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}

void mat_mul6(float A[6][6], float B[6][6], float C[6][6]) {
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      float sum = 0;
      for (int k = 0; k < 6; k++) {
        sum += A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}

void mat_add6(float A[6][6], float B[6][6], float C[6][6]) {
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
}

void mat_sub6(float A[6][6], float B[6][6], float C[6][6]) {
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      C[i][j] = A[i][j] - B[i][j];
    }
  }
}

void cholesky(float A[6][6], float B[6][6]) {
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j <= i; j++) {
      float sum = 0;
      if (j == i) {
        for (int k = 0; k < j; k++)
          sum += B[j][k] * B[j][k];
        B[j][j] = sqrtf(A[j][j] - sum);
      } else {
        for (int k = 0; k < j; k++)
          sum += (B[i][k] * B[j][k]);
        B[i][j] = (A[i][j] - sum) / B[j][j];
      }
    }
  }
}

void kalman_predict(Kalman *k, float g[3], float P_est[6][6]) {
  float theta = 0;
  float ug[3];
  // Remove bias from gyro readings
  ug[0] = g[0] - k->b[0];
  ug[1] = g[1] - k->b[1];
  ug[2] = g[2] - k->b[2];

  // Normalize gyro vector
  float mag_g = magnitude(ug);
  if (mag_g == 0) {
    theta = 0;
    ug[0] = 1;
    ug[1] = 1;
    ug[2] = 1;
  } else {
    theta = mag_g * k->dt;
    float norm_mag_g = 1 / mag_g;
    ug[0] *= norm_mag_g;
    ug[1] *= norm_mag_g;
    ug[2] *= norm_mag_g;
  }

  // Compute deltaq
  float dq[4];
  dq[0] = cos(theta / 2);
  dq[1] = ug[0] * sin(theta / 2);
  dq[2] = ug[1] * sin(theta / 2);
  dq[3] = ug[2] * sin(theta / 2);

  // Compute and normalize qk+1
  float qk[4];
  qk[0] = (k->q[0] * dq[0]) - (dq[1] * k->q[1]) - (dq[2] * k->q[2]) -
          (dq[3] * k->q[3]);
  qk[1] = (dq[0] * k->q[1]) + (dq[1] * k->q[0]) + (dq[2] * k->q[3]) -
          (dq[3] * k->q[2]);
  qk[2] = (dq[0] * k->q[2]) - (dq[1] * k->q[3]) + (dq[2] * k->q[0]) +
          (dq[3] * k->q[1]);
  qk[3] = (dq[0] * k->q[3]) + (dq[1] * k->q[2]) - (dq[2] * k->q[1]) +
          (dq[3] * k->q[0]);

  float mag_q = magnitudeq(qk);
  mag_q = 1 / mag_q;
  k->q[0] = qk[0] * mag_q;
  k->q[1] = qk[1] * mag_q;
  k->q[2] = qk[2] * mag_q;
  k->q[3] = qk[3] * mag_q;

  // Find state transition matrix
  float phi[6][6] = {{1, ug[2] * k->dt, -ug[1] * k->dt, -k->dt, 0, 0},
                     {-ug[2] * k->dt, 1, ug[0] * k->dt, 0, -k->dt, 0},
                     {ug[1] * k->dt, -ug[0] * k->dt, 1, 0, 0, -k->dt},
                     {0, 0, 0, 1, 0, 0},
                     {0, 0, 0, 0, 1, 0},
                     {0, 0, 0, 0, 0, 1}};

  // Find noise covariance
  float Qd[6][6] = {{k->noise_g * k->noise_g * k->dt, 0, 0, 0, 0, 0},
                    {0, k->noise_g * k->noise_g * k->dt, 0, 0, 0, 0},
                    {0, 0, k->noise_g * k->noise_g * k->dt, 0, 0, 0},
                    {0, 0, 0, k->bias_g * k->bias_g * k->dt, 0, 0},
                    {0, 0, 0, 0, k->bias_g * k->bias_g * k->dt, 0},
                    {0, 0, 0, 0, 0, k->bias_g * k->bias_g * k->dt}};

  // Solve for P_est
  float phi_t[6][6];
  transpose6(phi, phi_t);
  mat_mul6(phi, k->P, P_est);
  mat_mul6(P_est, phi_t, phi);
  mat_add6(phi, Qd, P_est);
}

void kalman_correction(Kalman *k, float a[3], float m[3], float P_est[6][6]) {
  float Rq[3][3] = {{k->q[0] * k->q[0] + k->q[1] * k->q[1] - k->q[2] * k->q[2] -
                         k->q[3] * k->q[3],
                     2 * (k->q[1] * k->q[2] - k->q[0] * k->q[3]),
                     2 * (k->q[1] * k->q[3] + k->q[0] * k->q[2])},
                    {2 * (k->q[1] * k->q[2] + k->q[0] * k->q[3]),
                     k->q[0] * k->q[0] - k->q[1] * k->q[1] + k->q[2] * k->q[2] -
                         k->q[3] * k->q[3],
                     2 * (k->q[2] * k->q[3] - k->q[0] * k->q[1])},
                    {2 * (k->q[1] * k->q[3] - k->q[0] * k->q[2]),
                     2 * (k->q[2] * k->q[3] + k->q[0] * k->q[1]),
                     k->q[0] * k->q[0] - k->q[1] * k->q[1] - k->q[2] * k->q[2] +
                         k->q[3] * k->q[3]}};
  float gw[3] = {0, 0, -1};
  float mw[3] = {1, 0, 0};

  // Calculate residual
  float r[6] = {a[0] - Rq[0][2] * gw[0], a[1] - Rq[1][2] * gw[1],
                a[2] - Rq[2][2] * gw[2], m[0] - Rq[0][0],
                m[1] - Rq[0][1],         m[2] - Rq[0][2]};
  float s_gw[3][3] = {0};
  float s_mw[3][3] = {0};
  skew(gw, s_gw);
  skew(mw, s_mw);

  float Rgw[3][3] = {0}, Rmw[3][3] = {0};
  mat_mul3(Rq, s_gw, Rgw);
  mat_mul3(Rq, s_mw, Rmw);

  // Find measurement Jacobian
  float H[6][6] = {0}; // Set H to zero
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      H[i][j] = -Rgw[i][j];
      H[i + 3][j] = -Rmw[i][j];
    }
  }

  float noise_R[6][6] = {0};
  for (int i = 0; i < 3; i++) {
    noise_R[i][i] = k->noise_a * k->noise_a;
    noise_R[i + 3][i + 3] = k->noise_m * k->noise_m;
  }
  float Ht[6][6] = {0};
  transpose6(H, Ht);

  // Find intermediate/covariance projected onto measurement Jacobian
  float U[6][6] = {0};
  mat_mul6(P_est, Ht, U);

  // Start solving for K with S, and Cholesky decomp
  float S[6][6] = {0};
  mat_mul6(H, U, S);
  mat_add6(S, noise_R, S);

  float L[6][6] = {0};
  float Lt[6][6] = {0};
  cholesky(S, L);
  transpose6(L, Lt);

  // Forward Solve L * Y = U
  float Y[6][6] = {0};
  for (int j = 0; j < 6; ++j) {
    for (int i = 0; i < 6; ++i) {
      float s = U[i][j];
      for (int k = 0; k < i; ++k) {
        s -= L[i][k] * U[k][j];
      }
      Y[i][j] = s / L[i][i];
    }
  }

  // Backward Solve L^T * K = Y
  float K[6][6] = {0};
  for (int j = 0; j < 6; ++j) {
    for (int i = 5; i >= 0; --i) {
      float s = Y[i][j];
      for (int k = i + 1; k < 6; ++k) {
        s -= L[k][i] * K[k][j];
      }
      K[i][j] = s / L[i][i];
    }
  }

  // dx = K * r
  float dx[6] = {0};
  for (int i = 0; i < 6; ++i) {
    float s = 0.0f;
    for (int j = 0; j < 6; ++j) {
      s += K[i][j] * r[j];
    }
    dx[i] = s;
  }

  // delta theta and delta bias
  float dth[3] = {dx[0], dx[1], dx[2]};
  float db[3] = {dx[3], dx[4], dx[5]};

  // Find delta quaternion using delta theta
  float dq[4] = {1, 0.5f * dth[0], 0.5f * dth[1], 0.5f * dth[2]};

  // Find new quaternion at time step k
  float qk[4];
  qk[0] = (k->q[0] * dq[0]) - (dq[1] * k->q[1]) - (dq[2] * k->q[2]) -
          (dq[3] * k->q[3]);
  qk[1] = (dq[0] * k->q[1]) + (dq[1] * k->q[0]) + (dq[2] * k->q[3]) -
          (dq[3] * k->q[2]);
  qk[2] = (dq[0] * k->q[2]) - (dq[1] * k->q[3]) + (dq[2] * k->q[0]) +
          (dq[3] * k->q[1]);
  qk[3] = (dq[0] * k->q[3]) + (dq[1] * k->q[2]) - (dq[2] * k->q[1]) +
          (dq[3] * k->q[0]);

  // Normalize qk
  float mag_q = magnitudeq(qk);
  mag_q = 1 / mag_q;
  k->q[0] = qk[0] * mag_q;
  k->q[1] = qk[1] * mag_q;
  k->q[2] = qk[2] * mag_q;
  k->q[3] = qk[3] * mag_q;

  // Add/Subtract Bias
  k->b[0] += db[0];
  k->b[1] += db[1];
  k->b[2] += db[2];

  // Find and assign new covariance matrix P
  float Kt[6][6], UKt[6][6], P_new[6][6];
  transpose6(K, Kt);
  mat_mul6(U, Kt, UKt);
  mat_sub6(P_est, UKt, P_new);
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      k->P[i][j] = P_new[i][j];
    }
  }
}
