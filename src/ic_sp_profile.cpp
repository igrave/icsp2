#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <RcppEigen.h>
#include "ic_sp_ch.h"
#include "utilities.h"

double zeroin2(double ax, double bx, double fa, double fb,
               double (*f)(double, void*), void *info,
               double *Tol, int *Maxit);

using namespace std;
using namespace Rcpp;
using namespace Eigen;

struct ProfileInfo {
    icm_Abst* obj;
    int cov_i;
    double target_llk;
};

static double profile_target_wrapper(double try_cov, void* info) {
    ProfileInfo* pi = static_cast<ProfileInfo*>(info);
    pi->obj->reg_par[pi->cov_i] = try_cov;
    pi->obj->profile_cov_idx = pi->cov_i;
    pi->obj->updateCovars = true;
    pi->obj->update_etas();
    double new_llk = pi->obj->run(50, 1e-10, true, 5);
    return new_llk - pi->target_llk;
}

static double profile_target_wrapper(double try_cov, void* info) { 
    ProfileInfo* pi = static_cast<ProfileInfo*>(info); 
    pi->obj->reg_par[pi->cov_i] = try_cov; 
    pi->obj->profile_cov_idx = pi->cov_i; 
    pi->obj->updateCovars = true; 
    double new_llk = pi->obj->run(50, 1e-10, true, 5); 
    return new_llk - pi->target_llk; 
} 
 
void icm_Abst::profile_llk_search(double target_llk, double ref_llk, int cov_i, int direction){ 
    double diag = -reg_d2(cov_i, cov_i); 
    if (diag <= 0.0) { 
        Rprintf("warning: non-negative Hessian diagonal for covariate %d\n", cov_i); 
        return; 
    } 
    double fa = ref_llk - target_llk; 
    double ax = reg_par[cov_i]; 
    double step = 1.1 * 2 * fa / sqrt(diag) * (direction == 0 ? -1.0 : 1.0); 
    double bx = ax + step; 
 
    ProfileInfo info = { this, cov_i, target_llk }; 
 
    double fb = 1; 
    while (fb >= 0) { 
      fb = profile_target_wrapper(bx, &info);  // should be < 0 (past boundary) 
      if (fb >= 0) bx += 0.5/sqrt(diag); 
    } 
    double tol = 1e-6; 
    int maxit = 50; 
 
    double root = zeroin2(ax, bx, fa, fb, 
                          profile_target_wrapper, &info, 
                          &tol, &maxit); 
 
    reg_par[cov_i] = root; 
    profile_cov_idx = -1;  // unlock 
} 

// Adapted from R_zeroin2() under GPL
// src/library/stats/src/zeroin.c
double zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit)				/* Max # of iterations */
{
    double a, b, c, fc;
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = *Tol;

    if(fa == 0.0) { *Tol = 0.0; *Maxit = 0; return a; }
    if(fb == 0.0) { *Tol = 0.0; *Maxit = 0; return b; }

    while(maxit--) {
        double prev_step = b - a;
        double tol_act;
        double p, q;
        double new_step;

        if(std::fabs(fc) < std::fabs(fb)) {
            a = b;  b = c;  c = a;
            fa = fb; fb = fc; fc = fa;
        }
        tol_act = 2.0 * DBL_EPSILON * std::fabs(b) + tol / 2.0;
        new_step = (c - b) / 2.0;

        if(std::fabs(new_step) <= tol_act || fb == 0.0) {
            *Maxit -= maxit;
            *Tol = std::fabs(c - b);
            return b;
        }

        if(std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
            double t1, cb, t2;
            cb = c - b;
            if(a == c) {
                t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                q = fa / fc; t1 = fb / fc; t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if(p > 0.0) q = -q; else p = -p;
            if(p < (0.75 * cb * q - std::fabs(tol_act * q) / 2.0)
               && p < std::fabs(prev_step * q / 2.0))
                new_step = p / q;
        }

        if(std::fabs(new_step) < tol_act) {
            new_step = (new_step > 0.0) ? tol_act : -tol_act;
        }
        a = b; fa = fb;
        b += new_step;
        fb = f(b, info);
        if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa;
        }
    }
    *Tol = std::fabs(c - b);
    *Maxit = -1;
    return b;
}