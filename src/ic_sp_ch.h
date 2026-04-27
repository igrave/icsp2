//
//  ic_sp_ch.h
//  
//
//  Created by Cliff Anderson Bergman on 5/25/15.
//
//

#ifndef ____ic_sp_ch__
#define ____ic_sp_ch__
/*#include "../Eigen_local/Dense"
#include <stdio.h>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>  */

//using namespace std;
//#include "../icenReg_files/basicUtilities.cpp"

#include <stdio.h>
#include <vector>
#include <cmath>
#include <omp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Compute log(1 - exp(-x)) for x > 0, numerically stable.
// Uses the formulation from Mächler (2012).
inline double log1mexp(double x) {
    if (x <= M_LN2) {
        return std::log(-std::expm1(-x));
    } else {
        return std::log1p(-std::exp(-x));
    }
}

class node_info{
public:
    std::vector<int> l;      //vector that indicates the observations for which this node is the left side
    std::vector<int> r;      //vector that indicated the observations for which this node is the right side
//    double par;         //log cumulative hazard
};

class obInf{
public:
    int l,r;
    double pob;
    double log_pob;
};


class icm_Abst{
public:
    virtual void update_p_ob(int s, int i);
    
    
    double sum_llk_all(); //done, not checked
    // calculates the entire likelihood function.
    // Does not update eta or hazards!
    double sum_llk(int s);    // calculates likelihood for a single stratum


    double par_llk(int s, int ind);     //done, not checked
    // only calculates partial likelihood based on an active index
    
    std::vector<std::vector<obInf>> obs_inf;
    std::vector<std::vector<node_info>> node_inf;
    
    void numericBaseDervsAllRaw(int s, std::vector<double> &d1, std::vector<double> &d2);
    
    void icm_addPar(int s, std::vector<double> &delta);

    void numericBaseDervsOne(int s, int raw_ind, std::vector<double> &d);
    
    void analytical_dobs_dch(int s, std::vector<double> &d1, std::vector<double> &d2);

    void update_etas();
	virtual void stablizeBCH() = 0;
    void recenterBCH();
	
    void icm_step();
    void icm_step_s(int s);
    
    void numericRegDervs();
    void covar_nr_step();
    
    virtual double basHaz2CondS(double ch, double eta) = 0;     //done
    virtual double baseS2CondS(double s, double eta) = 0;
    virtual double base_d1_contr(double h, double pob, double eta) = 0; //done, not checked
    virtual double reg_d1_lnk(double ch, double xb, double log_p) = 0;
    virtual double reg_d2_lnk(double ch, double xb, double log_p) = 0;
    virtual double reg_d3_lnk(double ch, double xb, double log_p) = 0;

    // contributions of single observations to the likelihood derivatives wrt baseP and baseCH
    virtual double dllk_dp_i(double s_l, double s_r, double eta, double pob,  bool left, bool right) = 0;
    virtual std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double log_pob, bool left) = 0;
    
    void calcAnalyticRegDervs(Eigen::MatrixXd &hess, Eigen::VectorXd &d1);
    void calcFinalRegContr(Eigen::MatrixXd &hess, Eigen::VectorXd &d1, Eigen::VectorXd &d3);
    void rawDervs2ActDervs();
    
    std::vector<Eigen::VectorXd>     baseCH;     //Vector of baseline log cumulative hazards.
                                    //baseH[0] fixed to -Inf, baseH[k-1] = Inf
	std::vector<double> intercept;				//used for numerical stabilization
	
    std::vector<Eigen::VectorXd>     backupCH;   //used to save values in optimization steps
    Eigen::VectorXd     propVec;    //used for proposition step during NR update on regression parameters
 /*   Eigen::VectorXd     H_d1;       //Vector of derivatives for CH's
    Eigen::MatrixXd     H_d2;       //Hessian for CH's          */
    std::vector<Eigen::VectorXd>     base_p_obs; //Baseline probability of each observation  //initialized
    std::vector<Eigen::VectorXd>     etas;       //linear combination of regression parameters   //initialized
    std::vector<Eigen::VectorXd>     expEtas;    //exp(etas) //initialized
    Eigen::VectorXd     reg_par;    //regression parameters //initialized
    std::vector<Eigen::MatrixXd>     covars;     //covariates        //initialized
    Eigen::VectorXd     reg_d1;     //first derivatives of regression parameters        //initialized
    Eigen::MatrixXd     reg_d2;     //Hessian for derivatives       //initialized
    Eigen::VectorXd     reg_d3;     //third derivatives of regression parameters
//    Eigen::VectorXd     reg_d2;     //second derivatives: ignoring off diagonals!

    std::vector<std::vector<double>> w;
    
    int n_strata;             //number of strata
    
    double maxBaseChg;      //Max change in baseline parameters during icm step
    double h;
    bool hasCovars;
    bool updateCovars;

    int derivMethod; // 1 = numeric, 2 = auto, 3 = auto vectorized
    
    bool startGD;
    std::vector<std::vector<double>> baseS;
    std::vector<std::vector<double>> baseP;
    std::vector<std::vector<double>> baseP_backup;
    //std::vector<double> d_cond_S_left;  // IG not used?
    //std::vector<double> d_cond_S_right; //IG not used?
    std::vector<std::vector<double>> base_p_derv;
    std::vector<std::vector<double>> base_p_derv2;			// For computing 2nd derivative
    std::vector<std::vector<double>> base_p_2ndDerv;
    std::vector<std::vector<double>> prop_p;
    double llk_from_p(int s);
    double numeric_p_der(int i);
    
    double dervConS_fromBaseS(double s, double eta);
    void baseCH_2_baseS(int s);
    void baseS_2_baseP(int s);
    void baseP_2_baseS(int s);
    void baseS_2_baseCH(int s);
    void calc_cond_S_derv();
    void calc_base_p_derv();
    double getMaxScaleSize(const std::vector<double> &p, const std::vector<double> &prop_p);
    void gradientDescent_step();
    void auto_base_p_derv(int s);
    void analytical_dobs_dp(int s);
    // void experimental_step();
    // void EM_step();
    
    std::vector<std::vector<double>> dob_dp_both;
    std::vector<std::vector<double>> dob_dp_rightOnly;

    std::vector<std::vector<int>> isActive; // for gradientDescent_step()

	double run(int maxIter, double tol, bool useGA, int baselineUpdates);
    
    void numeric_dobs_dp(int s, bool forGA);
    //void numeric_dobs2_d2p();
    
    virtual double cal_log_obs(double s1, double s2, double eta);
    
    
    std::vector<std::vector<bool>> usedVec;
    
    double almost_inf;
    std::vector<int>  failedGA_counts;
    int iter;
    int numBaselineIts;
    bool useFullHess;
    
    double exchangeAndUpdate(double delta, int i1, int i2);
    // REQUIRES baseP BEING UP TO DATE!!!
    
    std::vector<int> exchangeIndices;
    
    void checkCH(int s);

    virtual icm_Abst* clone() const = 0;
    
    void last_p_update();
    void vem();
    void exchange_p_opt(int i1, int i2);
    void vem_sweep();
    void vem_sweep2();

    void profile_llk_search(double target_llk, int cov_i, int direction);
};

void setup_icm(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, SEXP R_strata, icm_Abst* icm_obj);
//function for setting up a actSet_Abst class

void cumhaz2p_hat(Eigen::VectorXd &ch, std::vector<double> &p);

class icm_ph : public icm_Abst{
public:
    double basHaz2CondS(double ch, double eta){
        if(ch == R_NegInf)  return(1);
        if(ch == R_PosInf)  return(0);
        return(exp(-exp(ch + eta) )) ;}
    
    double baseS2CondS(double s, double eta){
        if(s >= 1.0) return(1.0);
        if(s <= 0.0) return(0.0);
/*        double expEta = exp(eta);
        double ans = pow(s, expEta);    */
        double logCH = log( -log(s) );
        double ans = exp(-exp(logCH + eta));
        return(ans);
    }
    
    double base_d1_contr(double ch, double pob, double eta){
        double expVal = -exp(eta + ch);
        double logAns = eta + ch + expVal - log(pob);
        return (-exp(logAns));
    }
    
    double reg_d1_lnk(double ch, double xb, double log_p){
        double term1 = -exp(ch + xb);
        return(-exp(term1 + ch + xb - log_p));
    }
    double reg_d2_lnk(double ch, double xb, double log_p){
        double term1 = -exp(ch + xb);
        double term2 = exp(term1 - log_p);
        return(term1 * term2 + term1 * term1 *term2);
    }
	double reg_d3_lnk(double ch, double xb, double log_p){
        double term1 = -exp(ch + xb);
        double term2 = exp(term1 - log_p);
        double term1_sq = term1 * term1;
        return(term1 * term2 + 3 * term1_sq * term2 + term1 * term1_sq * term2);
    }

	void stablizeBCH(){
        for(int s = 0; s < n_strata; s++){
            int k = baseCH[s].size();
		    double thisChange = baseCH[s][k-2] - 2.0;
		    intercept[s] += thisChange;
		    for(int i = 1; i < (k-1); i++){
    			baseCH[s][i] -= thisChange;
	    	} 
        }
        update_etas();
	}
	
    double dllk_dp_i(double s_l, double s_r, double eta, double pob, bool left, bool right){
        // no derivative terms
        if (!left && !right) return(0.0);
        
        /*
        Possibilities:
        sl ==1 && sr == llk_0 := 0
        sl < 1 && sr == 0 && p in sl
        sl < 1 && sr == 0 && p not in sl := 0

        sl == 1 && sr > 0 && p in sr
        sl == 1 && sr > 0 && p not in sr := 0
      
        sl < 1 && sr > 0 && p in sr, in sl
        sl < 1 && sr > 0 && p in sr, not in sl
        sl < 1 && sr > 0 && p not in sr, not in sl := 0
        
        
        */
        if (eta == 0.0) {
            if (left && right) {
                return(0.0);
            } else if (left) {
                return(-1.0 / s_l);
            } else if (right) {
                return(1.0 / (s_l - s_r));
            } else {
                Rcpp::Rcout << "Error in dllk_dp_i: both left and right are false!" << std::endl;
                return(0.0);
            }
        }

        double r_term, l_term;
        if (right) {
            r_term = exp(log(s_r) * (exp(eta) - 1));
        }
        double ans = 0;
        if (left && right){
            // exp(eta) * (r_term - l_term) / exp(pob)
            l_term = exp(log(s_l) * (exp(eta) - 1));
            //ans = exp(eta + log(r_term - l_term) - pob);
            ans = exp(eta) * (r_term - l_term) / exp(pob);
        } else if (left) { // but not right, can simplify
            // exp(eta) * l_term / exp(pob)
            // = - exp(eta) / s_l
            ans = - exp(eta) / s_l;
        } else if (right) { // but not left
            // exp(eta) * r_term / exp(pob)
            ans = exp(eta + log(r_term) - pob);
        } else {
            Rcpp::Rcout << "Error in dllk_dp_i: both left and right are false!" << std::endl;
        }
       
        if (ISNAN(ans)) {
            Rcpp::Rcout << "Warning: dllk_dp_i returned NaN!" << std::endl;
            ans = 0.0;
        }
        return(ans);
    }

    std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double log_pob, bool left){
        std::vector<double> ans(2);
        double d1, d2;
        // Derivative of S w.r.t. H: dS/dH = -exp(H+eta) * exp(-exp(H+eta))
        //   = -exp(a - exp(a)) where a = H + eta
        // d1 = dS/dH / P = -exp(a - exp(a) - log_pob)
        if (left) {
            double a = ch_l + eta;
            if (ch_r == R_PosInf) {
               d1 = -exp(a);
               d2 = d1;
            } else {
                double log_numer = a - exp(a);  // log(|dS_l/dH_l|)
                double log_abs_d1 = log_numer - log_pob;
                if (log_abs_d1 < -36.0 || log_numer < -700.0) {
                    d1 = 0.0; d2 = 0.0;
                } else {
                    d1 = -exp(log_abs_d1);
                    d2 = d1 * (1.0 - exp(a)) - d1 * d1;
                }
            }
        } else {
            double a = ch_r + eta;
            double log_numer = a - exp(a);
            double log_abs_d1 = log_numer - log_pob;
            if (log_abs_d1 < -36.0 || log_numer < -700.0) {
                d1 = 0.0; d2 = 0.0;
            } else {
                d1 = exp(log_abs_d1);
                d2 = d1 * (1.0 - exp(a)) - d1 * d1;
            }
        }

        ans[0] = d1;
        ans[1] = d2;
        return(ans);
    };

    // Numerically stable log(S(s1|eta) - S(s2|eta)) for PH model.
    // Computes log(s1^nu - s2^nu) where nu = exp(eta) in log-space,
    // avoiding underflow when s^nu is near zero.
    double cal_log_obs(double s1, double s2, double eta) {
        double nu = exp(eta);
        if (s1 >= 1.0 && s2 <= 0.0) {
            return 0.0;
        }
        if (s2 <= 0.0) {
            return nu * log(s1);
        }
        if (s1 >= 1.0) {
            // log(1 - s2^nu) = log(1 - exp(nu*log(s2)))
            return log1mexp(-nu * log(s2));
        }
        // General: log(s1^nu - s2^nu)
        //   = nu*log(s1) + log(1 - exp(nu*log(s2/s1)))
        //   = nu*log(s1) + log1mexp(nu*log(s1/s2))
        double log_ratio = log(s1) - log(s2);
        return nu * log(s1) + log1mexp(nu * log_ratio);
    }

    // Numerically stable update_p_ob for PH model.
    // Uses direct computation when possible (bit-exact with old code),
    // falls back to log-space when pob underflows to 0.
    void update_p_ob(int s, int i) {
        double chl = baseCH[s][ obs_inf[s][i].l ];
        double chr = baseCH[s][ obs_inf[s][i].r + 1 ];
        double eta = etas[s][i];
        double sl = basHaz2CondS(chl, eta);
        double sr = basHaz2CondS(chr, eta);
        obs_inf[s][i].pob = sl - sr;
        if (obs_inf[s][i].pob > 0.0) {
            obs_inf[s][i].log_pob = log(obs_inf[s][i].pob);
        } else {
            // pob underflowed to 0 or negative (monotonicity violation).
            // Use log-space: log(pob) = -exp(a) + log1mexp(exp(b) - exp(a))
            double ea = exp(chl + eta);
            double eb = exp(chr + eta);
            double diff = eb - ea;
            if (diff > 0.0) {
                obs_inf[s][i].log_pob = -ea + log1mexp(diff);
                obs_inf[s][i].pob = exp(obs_inf[s][i].log_pob);
            } else {
                obs_inf[s][i].log_pob = R_NegInf;
                obs_inf[s][i].pob = 0.0;
            }
        }
    }

    icm_Abst* clone() const override { return new icm_ph(*this); }
    virtual ~icm_ph(){};
};


class icm_po : public icm_Abst{
public:
    double basHaz2CondS(double ch, double eta){
        if(ch == R_NegInf)  return(1);
        if(ch == R_PosInf)  return(0);
        double mu = exp(ch);
        double s = exp(-mu);
        double s_nu = exp(eta - mu);
        return( (s_nu) / (s_nu - s + 1)) ;}
    
    double baseS2CondS(double s, double eta){
        double nu = exp(eta);
        double s_nu = s * nu;
        return((s_nu)/ (s_nu - s + 1) );
    }
    
    double base_d1_contr(double ch, double pob, double eta){
        double s = exp(-exp(ch));
        double s_nu = exp(eta - exp(ch));
        double logAns = -log(pob) - 2 * log(s_nu - s + 1) + ch - exp(h);
        return (-exp(logAns));
    }
    
    double reg_d1_lnk(double ch, double xb, double log_p){
        double s = exp(-exp(ch));
        double a = exp(xb-exp(ch));
        double ans = exp( log(a *(1-s)) - 2 * log( a - s + 1) - log_p);
        return(ans);
//        return( a * (1-s) / pow(a - s + 1, 2.0) );
    }
    double reg_d2_lnk(double ch, double xb, double log_p){
        double s = exp(-exp(ch));
        double a = exp(xb-exp(ch));
        double b = a - s + 1;
        double top =  (a * (1 - s) * b - 2 * a * a *(1-s)) ;
        double bottom = pow(b, 3.0);
        double ans = top/(bottom * exp(log_p));
        return(ans);
    }
    double reg_d3_lnk(double ch, double xb, double log_p){
        double s = exp(-exp(ch));
        double a = exp(xb - exp(ch));
        double b = a - s + 1;
        double top = a * (1.0 - s) * (b * b - 6.0 * a * b + 6.0 * a * a);
        double bottom = pow(b, 4.0) * exp(log_p);
        return(top / bottom);
    }

    
	void stablizeBCH(){}
	
    double dllk_dp_i(double s_l, double s_r, double eta, double pob, bool left, bool right){
        // no derivative terms
        if (!left && !right) return(0.0);
        if (eta == 0.0) {
            if (left && right) {
                return(0.0);
            } else if (left) {
                return(-1.0 / s_l);
            } else if (right) {
                return(1.0 / (s_l - s_r));
            } else {
                Rcpp::Rcout << "Error in dllk_dp_i: both left and right are false!" << std::endl;
                return(0.0);
            }
        }
        double l_term, r_term;
        if (left) {
            double denom = s_l * (exp(eta) - 1) + 1;
            l_term = s_l * (1-exp(eta)) / (denom * denom) + 1 / denom;

        } else {
            l_term = 0;
        }

        if (right) {
            double denom = s_r * (exp(eta) - 1) + 1;
            r_term = s_r  * (1-exp(eta)) / (denom * denom) + 1 / denom;
        } else {
            r_term = 0;
        }

       // double ret = exp(eta) * (r_term - l_term) / exp(pob)

        //logged version
        //double ans = exp(eta + log(r_term - l_term) - pob);
        
        double ans = exp(eta) * (r_term - l_term) / exp(pob);
        return(ans);
    }


/* std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double pob, bool left){
        std::vector<double> ans(9);
        double d1, d2;
        
          double eHl = exp(ch_l);
          double eHr = exp(ch_r);
          double eeHr = exp(eHr);
          double eeHl = exp(eHl);
          double eEta = exp(eta);
          double diff_ee = eeHl - eeHr;
          double inv_diff_ee = 1/diff_ee;
          double term_r = eeHr + eEta - 1;
          double term_l = eeHl + eEta - 1;
        
        // Special cases:
        if (left && ch_l == R_NegInf) {
          d1 = 0;
          d2 = 0;
        } else if (!left && ch_r == R_PosInf) {
          d1 = 0;
          d2 = 0;
        } else {
          // anything common for all calculations?
          // put calcs back here later
        
          if (left) {
            double inv_term_l = 1/ term_l;
            if (ch_r == R_PosInf) {
              d1 = eeHl * eHl * inv_term_l;
              d2 = eeHl * eHl * (
                  eeHl * eHl * inv_term_l * inv_term_l -
                  (1 + eHl) * inv_term_l
              );
            } else {
            
              d1 = eeHl * eHl * term_r * inv_diff_ee * inv_term_l;

              d2 = eeHl * eHl * (
                  (diff_ee - eeHr * eHl) * inv_diff_ee * inv_diff_ee +
                  eeHl * eHl * inv_term_l * inv_term_l -
                  (1 + eHl) * inv_term_l
                );
                }// left special case
                //end left
            } else { // right normal case
              double inv_term_r = 1 / term_r;
            
            // do we need special case for ch_l == R_NegInf?
            d1 = eeHr * eHr * term_l * -inv_diff_ee * inv_term_r;
            
            d2 = eeHr * eHr * (
                -(diff_ee + eeHl * eHr) * inv_diff_ee * inv_diff_ee +
                eeHr * eHr * inv_term_r*inv_term_r -
                (1 + eHr) * inv_term_r
            );
         }
        }
        //ans[0] = eHl;
        //ans[1] = eeHl;
        //ans[2] = eHr;
        //ans[3] = eeHr;
        //ans[4] = eEta;
        //ans[5] = 1/term_l;
        //ans[6] = 1/term_r;
        ans[0] = d1;
        ans[1] = d2;
        return(ans);
    };
 */

  std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double log_pob, bool left) {
        std::vector<double> ans(2);
        double d1 = 0.0, d2 = 0.0;
        
        // Define common terms
        double eHl = exp(ch_l);
        double eHr = exp(ch_r);
        double eeHr = exp(eHr);
        double eeHl = exp(eHl);
        double eEta = exp(eta);
        
        // diff_ee = E_l - E_r (This is practically always negative)
        double diff_ee = eeHl - eeHr; 
        double inv_diff_ee = 1.0 / diff_ee;
        
        // D terms: D = exp(eta) - 1 + E
        double term_r = eeHr + eEta - 1.0;
        double term_l = eeHl + eEta - 1.0;
        double inv_term_l = 1.0 / term_l;
        double inv_term_r = 1.0 / term_r;

        // Derivative of E w.r.t H: d(eeH)/dH = eeH * eH
        double X_l_prime = eeHl * eHl;
        double X_r_prime = eeHr * eHr;

        if (left) {
            if (ch_l == R_NegInf) {
                d1 = 0.0; d2 = 0.0;
            } else if (ch_r == R_PosInf) {
                // Case: Right boundary is infinity -> S(Right) = 0
                // Proportional Odds Survival S(t) = exp(eta) / ( exp(eta) - 1 + exp(exp(H)) )
                // P = S(Left); log L = eta - log(term_l)
                // d1 = - exp(eta)/term_l^2 * (eeHl * eHl) * (term_l/exp(eta)) ?? No, simpler:
                // d/dH log(S(Left)) = - (1/term_l) * (eeHl * eHl)
                d1 = - X_l_prime * inv_term_l;
                
                // d2 = d1 * d/dH[ log(d1) ]
                // log(|d1|) = H + eH - log(term_l)
                // deriv = 1 + eH - (eeHl * eHl)/term_l
                double d_log_d1 = (1.0 + eHl) - (X_l_prime * inv_term_l);
                d2 = d1 * d_log_d1;
                
            } else {
                // Normal Left Case
                // d1 = (eeHl * eHl * term_r) / ( (eeHl - eeHr) * term_l )
                d1 = X_l_prime * term_r * inv_diff_ee * inv_term_l;

                // d2 = d1 * d/dHl [ log |d1| ]
                // log |d1| = H + eH + log(term_r) - log|diff_ee| - log(term_l)
                // Note: term_r is constant w.r.t H_l
                // d/dHl log|d1| = 1 + eHl - d/dHl(log|E_l - E_r|) - d/dHl(log term_l)
                //               = 1 + eHl - (X_l_prime / diff_ee) - (X_l_prime / term_l)
                
                double d_log_d1 = (1.0 + eHl) - (X_l_prime * inv_diff_ee) - (X_l_prime * inv_term_l);
                d2 = d1 * d_log_d1;
            }
        } else { // Right Boundary
             if (ch_r == R_PosInf) {
                 d1 = 0.0; d2 = 0.0;
             } else {
                 // Normal Right Case
                 // d1 = - (eeHr * eHr * term_l) / ( (eeHl - eeHr) * term_r )
                 //    = (eeHr * eHr * term_l) * ( - inv_diff_ee) * inv_term_r
                 d1 = - X_r_prime * term_l * inv_diff_ee * inv_term_r;

                 // d2 = d1 * d/dHr [ log |d1| ]
                 // log |d1| = H + eH + const - log|diff_ee| - log(term_r)
                 // deriv = 1 + eHr - ( - X_r_prime / diff_ee ) - ( X_r_prime / term_r )
                 // Note: d/dHr(diff_ee) = - X_r_prime. 
                 // So d/dHr log|diff_ee| = (1/diff_ee) * (-X_r_prime)
                 
                 double d_log_d1 = (1.0 + eHr) + (X_r_prime * inv_diff_ee) - (X_r_prime * inv_term_r);
                 d2 = d1 * d_log_d1;
             }
        }
        ans[0] = d1;
        ans[1] = d2;
        return(ans);
    };

    // Numerically stable log(S(s1|eta) - S(s2|eta)) for PO model.
    // Uses the identity: S_l - S_r = nu*(s_l - s_r) / (D_l * D_r)
    // where D = s*(nu-1) + 1, avoiding cancellation.
    double cal_log_obs(double s1, double s2, double eta) {
        if (s1 >= 1.0 && s2 <= 0.0) return 0.0;
        double nu = exp(eta);
        double D1 = s1 * (nu - 1.0) + 1.0;
        double D2 = s2 * (nu - 1.0) + 1.0;
        return eta + log(s1 - s2) - log(D1) - log(D2);
    }

    icm_Abst* clone() const override { return new icm_po(*this); }
    virtual ~icm_po(){};
};

#endif /* defined(____ic_sp_cm__) */
