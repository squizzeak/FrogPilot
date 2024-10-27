#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6740860442353781028) {
   out_6740860442353781028[0] = delta_x[0] + nom_x[0];
   out_6740860442353781028[1] = delta_x[1] + nom_x[1];
   out_6740860442353781028[2] = delta_x[2] + nom_x[2];
   out_6740860442353781028[3] = delta_x[3] + nom_x[3];
   out_6740860442353781028[4] = delta_x[4] + nom_x[4];
   out_6740860442353781028[5] = delta_x[5] + nom_x[5];
   out_6740860442353781028[6] = delta_x[6] + nom_x[6];
   out_6740860442353781028[7] = delta_x[7] + nom_x[7];
   out_6740860442353781028[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_551001988110241683) {
   out_551001988110241683[0] = -nom_x[0] + true_x[0];
   out_551001988110241683[1] = -nom_x[1] + true_x[1];
   out_551001988110241683[2] = -nom_x[2] + true_x[2];
   out_551001988110241683[3] = -nom_x[3] + true_x[3];
   out_551001988110241683[4] = -nom_x[4] + true_x[4];
   out_551001988110241683[5] = -nom_x[5] + true_x[5];
   out_551001988110241683[6] = -nom_x[6] + true_x[6];
   out_551001988110241683[7] = -nom_x[7] + true_x[7];
   out_551001988110241683[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7585604638564912887) {
   out_7585604638564912887[0] = 1.0;
   out_7585604638564912887[1] = 0;
   out_7585604638564912887[2] = 0;
   out_7585604638564912887[3] = 0;
   out_7585604638564912887[4] = 0;
   out_7585604638564912887[5] = 0;
   out_7585604638564912887[6] = 0;
   out_7585604638564912887[7] = 0;
   out_7585604638564912887[8] = 0;
   out_7585604638564912887[9] = 0;
   out_7585604638564912887[10] = 1.0;
   out_7585604638564912887[11] = 0;
   out_7585604638564912887[12] = 0;
   out_7585604638564912887[13] = 0;
   out_7585604638564912887[14] = 0;
   out_7585604638564912887[15] = 0;
   out_7585604638564912887[16] = 0;
   out_7585604638564912887[17] = 0;
   out_7585604638564912887[18] = 0;
   out_7585604638564912887[19] = 0;
   out_7585604638564912887[20] = 1.0;
   out_7585604638564912887[21] = 0;
   out_7585604638564912887[22] = 0;
   out_7585604638564912887[23] = 0;
   out_7585604638564912887[24] = 0;
   out_7585604638564912887[25] = 0;
   out_7585604638564912887[26] = 0;
   out_7585604638564912887[27] = 0;
   out_7585604638564912887[28] = 0;
   out_7585604638564912887[29] = 0;
   out_7585604638564912887[30] = 1.0;
   out_7585604638564912887[31] = 0;
   out_7585604638564912887[32] = 0;
   out_7585604638564912887[33] = 0;
   out_7585604638564912887[34] = 0;
   out_7585604638564912887[35] = 0;
   out_7585604638564912887[36] = 0;
   out_7585604638564912887[37] = 0;
   out_7585604638564912887[38] = 0;
   out_7585604638564912887[39] = 0;
   out_7585604638564912887[40] = 1.0;
   out_7585604638564912887[41] = 0;
   out_7585604638564912887[42] = 0;
   out_7585604638564912887[43] = 0;
   out_7585604638564912887[44] = 0;
   out_7585604638564912887[45] = 0;
   out_7585604638564912887[46] = 0;
   out_7585604638564912887[47] = 0;
   out_7585604638564912887[48] = 0;
   out_7585604638564912887[49] = 0;
   out_7585604638564912887[50] = 1.0;
   out_7585604638564912887[51] = 0;
   out_7585604638564912887[52] = 0;
   out_7585604638564912887[53] = 0;
   out_7585604638564912887[54] = 0;
   out_7585604638564912887[55] = 0;
   out_7585604638564912887[56] = 0;
   out_7585604638564912887[57] = 0;
   out_7585604638564912887[58] = 0;
   out_7585604638564912887[59] = 0;
   out_7585604638564912887[60] = 1.0;
   out_7585604638564912887[61] = 0;
   out_7585604638564912887[62] = 0;
   out_7585604638564912887[63] = 0;
   out_7585604638564912887[64] = 0;
   out_7585604638564912887[65] = 0;
   out_7585604638564912887[66] = 0;
   out_7585604638564912887[67] = 0;
   out_7585604638564912887[68] = 0;
   out_7585604638564912887[69] = 0;
   out_7585604638564912887[70] = 1.0;
   out_7585604638564912887[71] = 0;
   out_7585604638564912887[72] = 0;
   out_7585604638564912887[73] = 0;
   out_7585604638564912887[74] = 0;
   out_7585604638564912887[75] = 0;
   out_7585604638564912887[76] = 0;
   out_7585604638564912887[77] = 0;
   out_7585604638564912887[78] = 0;
   out_7585604638564912887[79] = 0;
   out_7585604638564912887[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6456710635463231286) {
   out_6456710635463231286[0] = state[0];
   out_6456710635463231286[1] = state[1];
   out_6456710635463231286[2] = state[2];
   out_6456710635463231286[3] = state[3];
   out_6456710635463231286[4] = state[4];
   out_6456710635463231286[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6456710635463231286[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6456710635463231286[7] = state[7];
   out_6456710635463231286[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1468259014266406367) {
   out_1468259014266406367[0] = 1;
   out_1468259014266406367[1] = 0;
   out_1468259014266406367[2] = 0;
   out_1468259014266406367[3] = 0;
   out_1468259014266406367[4] = 0;
   out_1468259014266406367[5] = 0;
   out_1468259014266406367[6] = 0;
   out_1468259014266406367[7] = 0;
   out_1468259014266406367[8] = 0;
   out_1468259014266406367[9] = 0;
   out_1468259014266406367[10] = 1;
   out_1468259014266406367[11] = 0;
   out_1468259014266406367[12] = 0;
   out_1468259014266406367[13] = 0;
   out_1468259014266406367[14] = 0;
   out_1468259014266406367[15] = 0;
   out_1468259014266406367[16] = 0;
   out_1468259014266406367[17] = 0;
   out_1468259014266406367[18] = 0;
   out_1468259014266406367[19] = 0;
   out_1468259014266406367[20] = 1;
   out_1468259014266406367[21] = 0;
   out_1468259014266406367[22] = 0;
   out_1468259014266406367[23] = 0;
   out_1468259014266406367[24] = 0;
   out_1468259014266406367[25] = 0;
   out_1468259014266406367[26] = 0;
   out_1468259014266406367[27] = 0;
   out_1468259014266406367[28] = 0;
   out_1468259014266406367[29] = 0;
   out_1468259014266406367[30] = 1;
   out_1468259014266406367[31] = 0;
   out_1468259014266406367[32] = 0;
   out_1468259014266406367[33] = 0;
   out_1468259014266406367[34] = 0;
   out_1468259014266406367[35] = 0;
   out_1468259014266406367[36] = 0;
   out_1468259014266406367[37] = 0;
   out_1468259014266406367[38] = 0;
   out_1468259014266406367[39] = 0;
   out_1468259014266406367[40] = 1;
   out_1468259014266406367[41] = 0;
   out_1468259014266406367[42] = 0;
   out_1468259014266406367[43] = 0;
   out_1468259014266406367[44] = 0;
   out_1468259014266406367[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1468259014266406367[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1468259014266406367[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1468259014266406367[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1468259014266406367[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1468259014266406367[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1468259014266406367[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1468259014266406367[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1468259014266406367[53] = -9.8000000000000007*dt;
   out_1468259014266406367[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1468259014266406367[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1468259014266406367[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1468259014266406367[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1468259014266406367[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1468259014266406367[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1468259014266406367[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1468259014266406367[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1468259014266406367[62] = 0;
   out_1468259014266406367[63] = 0;
   out_1468259014266406367[64] = 0;
   out_1468259014266406367[65] = 0;
   out_1468259014266406367[66] = 0;
   out_1468259014266406367[67] = 0;
   out_1468259014266406367[68] = 0;
   out_1468259014266406367[69] = 0;
   out_1468259014266406367[70] = 1;
   out_1468259014266406367[71] = 0;
   out_1468259014266406367[72] = 0;
   out_1468259014266406367[73] = 0;
   out_1468259014266406367[74] = 0;
   out_1468259014266406367[75] = 0;
   out_1468259014266406367[76] = 0;
   out_1468259014266406367[77] = 0;
   out_1468259014266406367[78] = 0;
   out_1468259014266406367[79] = 0;
   out_1468259014266406367[80] = 1;
}
void h_25(double *state, double *unused, double *out_5279605555161969149) {
   out_5279605555161969149[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4190957892875844299) {
   out_4190957892875844299[0] = 0;
   out_4190957892875844299[1] = 0;
   out_4190957892875844299[2] = 0;
   out_4190957892875844299[3] = 0;
   out_4190957892875844299[4] = 0;
   out_4190957892875844299[5] = 0;
   out_4190957892875844299[6] = 1;
   out_4190957892875844299[7] = 0;
   out_4190957892875844299[8] = 0;
}
void h_24(double *state, double *unused, double *out_6297754289708247101) {
   out_6297754289708247101[0] = state[4];
   out_6297754289708247101[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4587603197229703794) {
   out_4587603197229703794[0] = 0;
   out_4587603197229703794[1] = 0;
   out_4587603197229703794[2] = 0;
   out_4587603197229703794[3] = 0;
   out_4587603197229703794[4] = 1;
   out_4587603197229703794[5] = 0;
   out_4587603197229703794[6] = 0;
   out_4587603197229703794[7] = 0;
   out_4587603197229703794[8] = 0;
   out_4587603197229703794[9] = 0;
   out_4587603197229703794[10] = 0;
   out_4587603197229703794[11] = 0;
   out_4587603197229703794[12] = 0;
   out_4587603197229703794[13] = 0;
   out_4587603197229703794[14] = 1;
   out_4587603197229703794[15] = 0;
   out_4587603197229703794[16] = 0;
   out_4587603197229703794[17] = 0;
}
void h_30(double *state, double *unused, double *out_8227331545517454110) {
   out_8227331545517454110[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2725732448615772456) {
   out_2725732448615772456[0] = 0;
   out_2725732448615772456[1] = 0;
   out_2725732448615772456[2] = 0;
   out_2725732448615772456[3] = 0;
   out_2725732448615772456[4] = 1;
   out_2725732448615772456[5] = 0;
   out_2725732448615772456[6] = 0;
   out_2725732448615772456[7] = 0;
   out_2725732448615772456[8] = 0;
}
void h_26(double *state, double *unused, double *out_3225658277856608465) {
   out_3225658277856608465[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7932461211749900523) {
   out_7932461211749900523[0] = 0;
   out_7932461211749900523[1] = 0;
   out_7932461211749900523[2] = 0;
   out_7932461211749900523[3] = 0;
   out_7932461211749900523[4] = 0;
   out_7932461211749900523[5] = 0;
   out_7932461211749900523[6] = 0;
   out_7932461211749900523[7] = 1;
   out_7932461211749900523[8] = 0;
}
void h_27(double *state, double *unused, double *out_2242819503717636464) {
   out_2242819503717636464[0] = state[3];
}
void H_27(double *state, double *unused, double *out_550969136815347545) {
   out_550969136815347545[0] = 0;
   out_550969136815347545[1] = 0;
   out_550969136815347545[2] = 0;
   out_550969136815347545[3] = 1;
   out_550969136815347545[4] = 0;
   out_550969136815347545[5] = 0;
   out_550969136815347545[6] = 0;
   out_550969136815347545[7] = 0;
   out_550969136815347545[8] = 0;
}
void h_29(double *state, double *unused, double *out_8830647430193931627) {
   out_8830647430193931627[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1162393590054203488) {
   out_1162393590054203488[0] = 0;
   out_1162393590054203488[1] = 1;
   out_1162393590054203488[2] = 0;
   out_1162393590054203488[3] = 0;
   out_1162393590054203488[4] = 0;
   out_1162393590054203488[5] = 0;
   out_1162393590054203488[6] = 0;
   out_1162393590054203488[7] = 0;
   out_1162393590054203488[8] = 0;
}
void h_28(double *state, double *unused, double *out_1299309650815974029) {
   out_1299309650815974029[0] = state[0];
}
void H_28(double *state, double *unused, double *out_6244792607123734062) {
   out_6244792607123734062[0] = 1;
   out_6244792607123734062[1] = 0;
   out_6244792607123734062[2] = 0;
   out_6244792607123734062[3] = 0;
   out_6244792607123734062[4] = 0;
   out_6244792607123734062[5] = 0;
   out_6244792607123734062[6] = 0;
   out_6244792607123734062[7] = 0;
   out_6244792607123734062[8] = 0;
}
void h_31(double *state, double *unused, double *out_6396303292197231531) {
   out_6396303292197231531[0] = state[8];
}
void H_31(double *state, double *unused, double *out_4160311930998883871) {
   out_4160311930998883871[0] = 0;
   out_4160311930998883871[1] = 0;
   out_4160311930998883871[2] = 0;
   out_4160311930998883871[3] = 0;
   out_4160311930998883871[4] = 0;
   out_4160311930998883871[5] = 0;
   out_4160311930998883871[6] = 0;
   out_4160311930998883871[7] = 0;
   out_4160311930998883871[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6740860442353781028) {
  err_fun(nom_x, delta_x, out_6740860442353781028);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_551001988110241683) {
  inv_err_fun(nom_x, true_x, out_551001988110241683);
}
void car_H_mod_fun(double *state, double *out_7585604638564912887) {
  H_mod_fun(state, out_7585604638564912887);
}
void car_f_fun(double *state, double dt, double *out_6456710635463231286) {
  f_fun(state,  dt, out_6456710635463231286);
}
void car_F_fun(double *state, double dt, double *out_1468259014266406367) {
  F_fun(state,  dt, out_1468259014266406367);
}
void car_h_25(double *state, double *unused, double *out_5279605555161969149) {
  h_25(state, unused, out_5279605555161969149);
}
void car_H_25(double *state, double *unused, double *out_4190957892875844299) {
  H_25(state, unused, out_4190957892875844299);
}
void car_h_24(double *state, double *unused, double *out_6297754289708247101) {
  h_24(state, unused, out_6297754289708247101);
}
void car_H_24(double *state, double *unused, double *out_4587603197229703794) {
  H_24(state, unused, out_4587603197229703794);
}
void car_h_30(double *state, double *unused, double *out_8227331545517454110) {
  h_30(state, unused, out_8227331545517454110);
}
void car_H_30(double *state, double *unused, double *out_2725732448615772456) {
  H_30(state, unused, out_2725732448615772456);
}
void car_h_26(double *state, double *unused, double *out_3225658277856608465) {
  h_26(state, unused, out_3225658277856608465);
}
void car_H_26(double *state, double *unused, double *out_7932461211749900523) {
  H_26(state, unused, out_7932461211749900523);
}
void car_h_27(double *state, double *unused, double *out_2242819503717636464) {
  h_27(state, unused, out_2242819503717636464);
}
void car_H_27(double *state, double *unused, double *out_550969136815347545) {
  H_27(state, unused, out_550969136815347545);
}
void car_h_29(double *state, double *unused, double *out_8830647430193931627) {
  h_29(state, unused, out_8830647430193931627);
}
void car_H_29(double *state, double *unused, double *out_1162393590054203488) {
  H_29(state, unused, out_1162393590054203488);
}
void car_h_28(double *state, double *unused, double *out_1299309650815974029) {
  h_28(state, unused, out_1299309650815974029);
}
void car_H_28(double *state, double *unused, double *out_6244792607123734062) {
  H_28(state, unused, out_6244792607123734062);
}
void car_h_31(double *state, double *unused, double *out_6396303292197231531) {
  h_31(state, unused, out_6396303292197231531);
}
void car_H_31(double *state, double *unused, double *out_4160311930998883871) {
  H_31(state, unused, out_4160311930998883871);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
