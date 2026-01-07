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
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

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
};


class icm_Abst{
public:
    void update_p_ob(int s, int i);    //done, not checked
    
    
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
    void numericBaseDervsAllAct(int s, std::vector<double> &d1, std::vector<double> &d2);
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

    // contributions of single observations to the likelihood derivatives wrt baseP and baseCH
    virtual double dllk_dp_i(double s_l, double s_r, double eta, double pob,  bool left, bool right) = 0;
    virtual std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double pob, bool left) = 0;
    
    void calcAnalyticRegDervs(Eigen::MatrixXd &hess, Eigen::VectorXd &d1);
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
    double getMaxScaleSize( std::vector<double> &p, std::vector<double> &prop_p);
    void gradientDescent_step();
    void auto_base_p_derv(int s);
    void analytical_dobs_dp(int s);
    // void experimental_step();
    // void EM_step();
    
    std::vector<std::vector<double>> dob_dp_both;
    std::vector<std::vector<double>> dob_dp_rightOnly;

	double run(int maxIter, double tol, bool useGA, int baselineUpdates);
    
    void numeric_dobs_dp(int s, bool forGA);
    //void numeric_dobs2_d2p();
    
    double cal_log_obs(double s1, double s2, double eta);
    
    
    std::vector<std::vector<bool>> usedVec;
    
    double almost_inf;
    int failedGA_counts;
    int iter;
    int numBaselineIts;
    bool useFullHess;
    
    double exchangeAndUpdate(double delta, int i1, int i2);
    // REQUIRES baseP BEING UP TO DATE!!!
    
    std::vector<int> exchangeIndices;
    
    void checkCH(int s);
    
    void last_p_update();
    void vem();
    void exchange_p_opt(int i1, int i2);
    void vem_sweep();
    void vem_sweep2();
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

    std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double pob, bool left){
        std::vector<double> ans(2);
        double d1, d2;
        double ech;

        if (left) {
            ech = exp(ch_l + eta);
            if (ch_r == R_PosInf) {
               d1 = -ech;
               d2 = d1;
            } else {
                d1 = -(ech * exp(-ech)) / pob;
                d2 = d1 * (1 - ech) - d1 * d1;
            }
        } else {
            ech = exp(ch_r + eta);
            d1 = (ech * exp(-ech)) / pob;
            d2 = d1 * (1 - ech) - d1 * d1;
        }

        ans[0] = d1;
        ans[1] = d2;
        return(ans);
    };

	
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


    std::vector<double> dllk_dch_i(double ch_l, double ch_r, double eta, double pob, bool left){
        std::vector<double> ans(2);
        ans[0] = 0.0;
        ans[1] = 0.0;
        return(ans);
    };


    virtual ~icm_po(){};
};

extern "C" {
SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType,
 			  SEXP R_w, SEXP R_strata, SEXP R_use_GD, SEXP R_maxiter,
 			  SEXP R_baselineUpdates, SEXP R_useFullHess, SEXP R_updateCovars,
 			  SEXP R_initialRegVals, SEXP R_derivMethod);
    SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals);
}
#endif /* defined(____ic_sp_cm__) */
