#ifndef UTILITIES_H
#define UTILITIES_H

#include <RcppEigen.h>
#include <vector>

// Maths utilities
double max(double a, double b);
double max(int a, int b);
double min(double a, double b);
double signVal(double x);

// Vector operations
void mult_vec(double a, std::vector<double> &vec);
void add_vec(double lambda, std::vector<double> &a, std::vector<double> &vec);

// PAVA algorithm
void weighted_pool(double *y, double *w, int start, int stop);
void weighted_pava(double *y, double *w, int *numberParameters);
void pava(double *y, double *w, int *np);
void pava_monotone(double* y, double* w, int* np);

void pavaForOptim(std::vector<double> &d1, std::vector<double> &d2, std::vector<double> &x, std::vector<double> &prop_delta);

// Matrix/vector conversions
void copyRmatrix_intoEigen(SEXP r_mat, Eigen::MatrixXd &e_mat);

// Statistical utilities
double directional_derv(const std::vector<double> &derv, const std::vector<double> &delta);
void makeUnitVector(std::vector<double> &v);


// Interval utilities
int isValueInInterval(double val, double l, double r);
int isValueInInterval(double val, int ind, std::vector<double>& lvec, std::vector<double>& rvec);
int findSurroundingVals(double val, std::vector<double>& leftVec, 
                        std::vector<double>& rightVec, bool isLeft);

#endif // UTILITIES_H
