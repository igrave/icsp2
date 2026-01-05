#ifndef UTILITIES_H
#define UTILITIES_H

#include <RcppEigen.h>
#include <vector>

// Maths utilities
double max(double a, double b);
double max(int a, int b);
double min(double a, double b);
double max(Eigen::VectorXd v);
double signVal(double x);

// Vector operations
void mult_vec(double a, std::vector<double> &vec);
void add_vec(double a, std::vector<double> &vec);
void add_vec(std::vector<double> &a, std::vector<double> &vec);
void add_vec(double lambda, std::vector<double> &a, std::vector<double> &vec);

// PAVA algorithm
void weighted_pool(double *y, double *w, int start, int stop);
void weighted_pava(double *y, double *w, int *numberParameters);
void pava(double *y, double *w, int *np);

void pavaForOptim(std::vector<double> &d1, std::vector<double> &d2, std::vector<double> &x, std::vector<double> &prop_delta);
void pavaForOptim(Eigen::VectorXd &d1, Eigen::VectorXd &d2, Eigen::VectorXd &x, Eigen::VectorXd &prop_delta);
void addIfNeeded(std::vector<int> &points, int l, int r, int max);

// Matrix/vector conversions
void copyRmatrix_intoEigen(SEXP r_mat, Eigen::MatrixXd &e_mat);
Rcpp::NumericMatrix eigen2RMat(Eigen::MatrixXd &e_mat);
void Rvec2eigen(SEXP r_vec, Eigen::VectorXd &e_vec);
Rcpp::NumericVector eigen2RVec(Eigen::VectorXd &e_vec);

// Distribution functions
double ic_dloglogistic(double x, double a, double b);
double ic_ploglogistic(double x, double a, double b);
double ic_qloglogistic(double p, double a, double b);
double ic_dlnorm(double x, double mu, double s);
double ic_plnorm(double x, double mu, double s);
double ic_dgeneralgamma(double x, double mu, double s, double Q);
double ic_pgeneralgamma(double q, double mu, double s, double Q);
double ic_qgeneralgamma(double p, double mu, double s, double Q);

// SEXP utilities
void SEXP2doubleVec(SEXP R_vec, std::vector<double> &c_vec);
void doubleVec2SEXP(std::vector<double> &c_vec, SEXP R_vec);
void SEXPIndex2intIndex(SEXP R_Inds, std::vector<int> &c_inds);
void indexVec2SEXP(std::vector<int> &c_vec, SEXP R_vec);

// Index operations
void findIndexDiffs(std::vector<int> &in1, std::vector<int> &in2,
                    std::vector<int> &in1not2, std::vector<int> &in2not1);
void drop_index(int d_ind, std::vector<int> &indVec);
void add_index(int a_ind, std::vector<int> &indVec);
std::vector<int> getSEXP_MatDims(SEXP R_mat);

// Statistical utilities
void getPosNegIndices(std::vector<double> &vals, std::vector<int> &isPos, std::vector<int> &isNeg);
void getRelValIndices(double relVal, std::vector<double> &vals, std::vector<int> &subIndex,
                      std::vector<int> &above, std::vector<int> &below, int *max, int *min);
double directional_derv(std::vector<double> &derv, std::vector<double> &delta);
void makeUnitVector(std::vector<double> &v);
void getUniqInts(int i1, int i2, std::vector<int> &uniqInts, 
                 std::vector<std::vector<int>> &vec_vec, std::vector<bool> &usedVec);

// Interval utilities
int isValueInInterval(double val, double l, double r);
int isValueInInterval(double val, int ind, std::vector<double>& lvec, std::vector<double>& rvec);
int findSurroundingVals(double val, std::vector<double>& leftVec, 
                        std::vector<double>& rightVec, bool isLeft);

// Matrix utilities
Eigen::MatrixXd xtx(Eigen::MatrixXd &x);
Eigen::MatrixXd xtx(Eigen::MatrixXd &x, int row_start, int row_end);
Eigen::MatrixXd copyRows(Eigen::MatrixXd &x, int row_start, int row_end);

#endif // UTILITIES_H
