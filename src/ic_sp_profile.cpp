//
//  ic_sp_gradDescent.cpp
//  
//
//  Created by Cliff Anderson Bergman on 10/9/15.
//
//
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <RcppEigen.h>
#include "ic_sp_ch.h"
#include "utilities.h"

using namespace std;
using namespace Rcpp;
using namespace Eigen;

/*   Do profile likelihood search for one covariate in one direction */

void icm_Abst::profile_llk_search(double target_llk, int cov_i, int direction){
// First let's get a limit based on the hessian diagonal
// search_width = 2.3 / sqrt(-diag(object$hessian))

double step = 2.3 / pow(-reg_d2(cov_i, cov_i), 1/2) * (direction == 0 ? -1 : 1);
reg_par[cov_i] += step;
updateCovars = FALSE;
update_etas();
run(50, 10e-10, true, 5);
int j = 1;
// R_zeroin2()
// set new betas
// update profile likelihood
// refit
}
