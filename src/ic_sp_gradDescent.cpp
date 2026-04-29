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

/*   Gradient Ascent Step */

void icm_Abst::baseCH_2_baseS(int s){
    int k = baseCH[s].size();
    baseS[s].resize(k);
    baseS[s][0] = 1.0;
    baseS[s][k-1] = 0.0;
    for(int i = 1; i < (k-1); i++){
        baseS[s][i] = exp(-exp(baseCH[s][i]));
    }
}

void icm_Abst::baseS_2_baseP(int s){
    int k = baseS[s].size() - 1;
    baseP[s].resize(k);
    for(int i = 0; i < k; i++){
        baseP[s][i] = baseS[s][i] - baseS[s][i+1];
    }
}

void icm_Abst::baseP_2_baseS(int s){
    int k = baseP[s].size();
    baseS[s].resize(k + 1);
    baseS[s][0] = 1.0;
    for(int i = 1; i < k; i++){
        baseS[s][i] = baseS[s][i-1] - baseP[s][i-1];
    }
    baseS[s][k] = 0.0;
}

void icm_Abst::baseS_2_baseCH(int s){
    int k = baseS[s].size();
    baseCH[s].resize(k);
    baseCH[s][0] = R_NegInf;
    baseCH[s][k-1] = R_PosInf;
    for(int i = 1; i < (k-1); i++){
        baseCH[s][i] = log(-log(baseS[s][i]));
    }
}

double icm_Abst::llk_from_p(int s){
    baseP_2_baseS(s);
    return sum_llk_from_S(s);
}

double icm_Abst::sum_llk_from_S(int s) { 
    int n = obs_inf[s].size(); 
    double ans = 0; 
    for (int i = 0; i < n; i++) { 
        double sl = baseS[s][obs_inf[s][i].l]; 
        double sr = baseS[s][obs_inf[s][i].r + 1]; 
        double eta = etas[s][i]; 
        double log_pob = cal_log_obs(sl, sr, eta); 
        ans += log_pob * w[s][i]; 
    } 
    if (ISNAN(ans)) ans = R_NegInf; 
    return ans; 
} 
 
 
double icm_Abst::getMaxScaleSize(const std::vector<double> &p, const std::vector<double> &prop_p){
    double max_scale = 2.0;
    int k = p.size();
    int k2 = prop_p.size();
    if(k != k2){
        Rprintf("warning: k != k2 in getMaxScaleSize k = %d, k2 = %d\n", k, k2);
        return(0.0);
    }
    double thisMax = 1.0;
    for(int i = 0; i < k; i++){
        if(prop_p[i] != 0 && p[i] > 0){
            thisMax = max(-p[i]/prop_p[i], (1.0 - p[i]) / prop_p[i]);
            if(ISNAN(thisMax)){
                thisMax = 1.0;
            }
        }
        max_scale = min(max_scale, thisMax);
    }
    
    return(max_scale);
}



void icm_Abst::gradientDescent_step(){
	if(std::accumulate(failedGA_counts.begin(), failedGA_counts.end(), 0) > 500){return;}
	
#ifdef _OPENMP
//#pragma omp parallel for schedule(dynamic)
#endif
    for(int s = 0; s < n_strata; s++){
        backupCH[s] = baseCH[s];
        baseCH_2_baseS(s);
        baseS_2_baseP(s);
        
        //double org_llk = sum_llk(s);
        double org_llk = sum_llk_from_S(s); 
        

        if (derivMethod >= 10) {
            analytical_dobs_dp(s);
        } else {
            numeric_dobs_dp(s, true);
        }

        int k = base_p_derv[s].size();

        prop_p[s].resize(k);
        double prop_mean = 0;
        int act_sum = 0;
        double new_llk;

        isActive[s].assign(k, 0);
        
        for(int i = 0; i < k; i++){
            if(baseP[s][i] > 0 && !ISNAN(base_p_derv[s][i]) ){
                isActive[s][i] = 1;
                act_sum++;
            }
            else { isActive[s][i] = 0; }
        }

        for(int i = 0; i < k; i++){
            if(isActive[s][i]){ prop_mean += base_p_derv[s][i]; }
        }

        if(act_sum == 0){
            failedGA_counts[s]++;
            baseCH[s] = backupCH[s];
            continue;
        }
        prop_mean = prop_mean / act_sum;

        for(int i = 0; i < k; i++){
           if(isActive[s][i]){ prop_p[s][i] = base_p_derv[s][i] - prop_mean;}
            else {prop_p[s][i] = 0.0;}
        }
        
        makeUnitVector(prop_p[s]);
        double scale_max = getMaxScaleSize(baseP[s], prop_p[s]);

        for(int i = 0; i < k; i++){ prop_p[s][i] *= -1.0; }
        scale_max = min(scale_max, getMaxScaleSize(baseP[s], prop_p[s]));
        for(int i = 0; i < k; i++){ prop_p[s][i] *= -1.0; }

        double delta_val = scale_max/2.0;
    
        delta_val = min(delta_val, h);
        delta_val = delta_val/10.0;
    
        double analytic_dd = directional_derv(base_p_derv[s], prop_p[s]);

        if(delta_val == 0){
            failedGA_counts[s]++;
            baseCH[s] = backupCH[s];
            continue; // continue with the next stratum
        }

        add_vec(delta_val, prop_p[s], baseP[s]);
        double llk_h = llk_from_p(s);
        add_vec(-2.0 * delta_val, prop_p[s], baseP[s]);
        double llk_l = llk_from_p(s);
        add_vec(delta_val, prop_p[s], baseP[s]);
        double llk_0 = llk_from_p(s);

        double d1 = ( llk_h - llk_l ) / ( 2 * delta_val );
        double d2 = (llk_h + llk_l - 2.0 * llk_0 ) / (delta_val * delta_val);

        if(iter % 2 ==0){ d1 = analytic_dd; }
        
        if(!R_FINITE(d2) || std::abs(d2) <= 1e-8) {
            failedGA_counts[s]++;
            baseCH[s] = backupCH[s];
            continue;
        }

        delta_val = -d1/d2;
	
        if(!(R_FINITE(delta_val) && delta_val > 0.0)){
            failedGA_counts[s]++;
            baseCH[s] = backupCH[s];
            continue;
        }

        scale_max = getMaxScaleSize(baseP[s], prop_p[s]);
        delta_val = min(delta_val, scale_max);
        add_vec(delta_val, prop_p[s], baseP[s]);

        new_llk = llk_from_p(s);
        mult_vec(-1.0, prop_p[s]);
        int tries = 0;
    
        double this_delta = delta_val;
    
        while(tries < 5 && new_llk < llk_0){
            tries++;
            this_delta = this_delta/2;
            add_vec(this_delta, prop_p[s], baseP[s]);
            new_llk = llk_from_p(s);
        }
        if(new_llk < llk_0){
            failedGA_counts[s]++;
            baseCH[s] = backupCH[s];
		    continue;
        }

        if(org_llk > new_llk){
		    failedGA_counts[s]++;
		    baseCH[s] = backupCH[s];
	    } else {
            baseP_2_baseS(s); 
            baseS_2_baseCH(s); 
        }
    } // end of loop over strata, s
}




double icm_Abst::cal_log_obs(double s1, double s2, double eta){
    double l = baseS2CondS(s1, eta);
    double r = baseS2CondS(s2, eta);
    if (l - r <= 0){
       return(R_NegInf);
    }
    return(log(l - r) );
}



void icm_Abst::numeric_dobs_dp(int s, bool forGA){    
    int p_k = baseS[s].size();
    int k = p_k - 1;
    int n = etas[s].size();
    dob_dp_both[s].resize(n);
    dob_dp_rightOnly[s].resize(n);
    int lind, rind;
	double h_mult = 0.0001;
   	h *= h_mult;

	if(forGA){
	    double sl, sr, llk_h,llk_l, this_eta, this_h;    
   	 
   		for(int i = 0; i < n; i++){
    	    sl = baseS[s][ obs_inf[s][i].l];
    	    sr = baseS[s][ obs_inf[s][i].r + 1];
    	    this_eta = etas[s][i];
    	    if(sl == 1.0 && sr == 0.0){
    	        dob_dp_rightOnly[s][i] = 0;
    	        dob_dp_both[s][i] = 0;
    	    }
    	    else if(sr == 0){
    	        dob_dp_rightOnly[s][i] = 0;
    	        this_h = min(sl/2.0, h);
    	        sl -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
    	        sl += this_h * 2.0;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        dob_dp_both[s][i] = (llk_h - llk_l) / (2 * this_h);
    	    }
    	    else if( sl == 1.0 ){
    	        this_h = min(sr / 2.0, h);
    	        sr -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
    	        sr += 2.0 * this_h;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        dob_dp_both[s][i] = (llk_h - llk_l)/(2*this_h);
    	        dob_dp_rightOnly[s][i] = dob_dp_both[s][i];
    	    }
    	    else{
    	        this_h = min(sr /2.0, h);
    	        sr -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
    	        sr += 2.0 * this_h;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        sr -= this_h;
    	        dob_dp_rightOnly[s][i] = (llk_h - llk_l)/(2*this_h);
    	        sr -= this_h;
    	        sl -= this_h;
    	        llk_h = cal_log_obs(sl, sr, this_eta);
            
    	        sr += 2.0 * this_h;
    	        sl += 2.0 * this_h;
    	        llk_l = cal_log_obs(sl, sr, this_eta);
    	        dob_dp_both[s][i] = (llk_h - llk_l)/(2*this_h);
    	        
    	    }
    	}
    }
    else{
    	for(int i = 0; i < p_k; i++){
    		dob_dp_both[s][i] = 0;
    		dob_dp_rightOnly[s][i] = 0;
    	}
    	double thisProb;
    	double num_n = n;
    	for(int i = 0; i < n; i++){
    		lind = obs_inf[s][i].l;
    		rind = obs_inf[s][i].r + 1;
    		thisProb = baseS[s][lind] - baseS[s][rind];
    		
            if(thisProb > 0.0){
               dob_dp_rightOnly[s][i] = 1.0 / (num_n * thisProb);
            } else {
                dob_dp_rightOnly[s][i] = 0.0;
            }
    	}
    }

    base_p_derv[s].resize(k);
	
    int k_l, k_r;
    node_info* nd;
    for(int j = k-1; j >=0; j--){
        nd = &node_inf[s][j+1];
        k_r = nd->r.size();
        k_l = nd->l.size();
        if(j != k-1){
            base_p_derv[s][j] = base_p_derv[s][j+1];
        }
        else{
            base_p_derv[s][j] = 0;
        }
        for(int i = 0; i < k_r; i++){
            rind = nd->r[i];
            base_p_derv[s][j] += dob_dp_rightOnly[s][rind] * w[s][rind];
        }
        for(int i = 0; i < k_l; i++){
            lind = nd->l[i];
            base_p_derv[s][j] -= dob_dp_rightOnly[s][lind] * w[s][lind];
            base_p_derv[s][j] += dob_dp_both[s][lind] * w[s][lind];
        }
    }

    h = h/h_mult;

}

void icm_Abst::analytical_dobs_dp(int s){    
    int p_k = baseS[s].size();
    int k = p_k - 1;
    int n = etas[s].size();
    base_p_derv[s].resize(k);
    int lind, rind;

    dob_dp_both[s].assign(n, 0.0);
    dob_dp_rightOnly[s].assign(n, 0.0);

    for (int j = 0; j < k; j++) {
        base_p_derv[s][j] = 0;
    }

    for (int i = 0; i < n; i++) {
        double sl = baseS[s][ obs_inf[s][i].l];
        double sr = baseS[s][ obs_inf[s][i].r + 1];
        double eta = etas[s][i];
        double pob = cal_log_obs(sl, sr, eta);
        int lind =  obs_inf[s][i].l;
        int rind =  obs_inf[s][i].r + 1;

        // Skip observations with degenerate probabilities
        if (pob == R_NegInf) continue;

        if (lind > 0 && rind < k) {
            // both sides could contribute
            dob_dp_both[s][i] = dllk_dp_i(sl, sr, eta, pob, true, true);
        }
        if (lind > 0 && rind == k) {
            // left only contribution
            dob_dp_both[s][i] = dllk_dp_i(sl, sr, eta, pob, true, false);
        }
        if (rind < k) {
            // right side only contribution
            dob_dp_rightOnly[s][i] = dllk_dp_i(sl, sr, eta, pob, false, true);
        }
    }

    base_p_derv[s].resize(k);

    int k_l, k_r;
    node_info* nd;
    for(int j = k-1; j >=0; j--){
        nd = &node_inf[s][j+1];
        k_r = nd->r.size();
        k_l = nd->l.size();
        if(j != k-1){
            base_p_derv[s][j] = base_p_derv[s][j+1];
        }
        else{
            base_p_derv[s][j] = 0;
        }
        for(int i = 0; i < k_r; i++){
            rind = nd->r[i];
            base_p_derv[s][j] += dob_dp_rightOnly[s][rind] * w[s][rind];
        }
        for(int i = 0; i < k_l; i++){
            lind = nd->l[i];
            base_p_derv[s][j] -= dob_dp_rightOnly[s][lind] * w[s][lind];
            base_p_derv[s][j] += dob_dp_both[s][lind] * w[s][lind];
        }
    }
}
