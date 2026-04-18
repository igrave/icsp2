//
//  ic_sp_cm.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/25/15.
//
//

#include "ic_sp_ch.h"
#include "utilities.h"


/*      LIKELIHOOD TOOLS        */
void icm_Abst::update_p_ob(int s, int i){
    double chl = baseCH[s][ obs_inf[s][i].l ];
    double chr = baseCH[s][ obs_inf[s][i].r +1 ];
    double eta = etas[s][i];
    obs_inf[s][i].pob = basHaz2CondS(chl, eta) - basHaz2CondS(chr, eta);
    obs_inf[s][i].log_pob = log(obs_inf[s][i].pob);
}

double icm_Abst::sum_llk(int s){
    int n = obs_inf[s].size();
    double ans = 0;
    for(int i = 0; i < n; i++){
        update_p_ob(s, i);
        ans += obs_inf[s][i].log_pob * w[s][i];
    }
    if(ISNAN(ans)) {ans = R_NegInf;}
    return(ans);
}

double icm_Abst::sum_llk_all(){
    double ans = 0;
    for(int s = 0; s < n_strata; s++){
        ans += sum_llk(s);
    }
    return(ans);
}

double icm_Abst::par_llk(int s, int ind){
    int num_l = node_inf[s][ind].l.size();
    int num_r = node_inf[s][ind].r.size();
    double ans = 0;
    int thisInd;
    for(int i = 0; i < num_l; i++){
        thisInd = node_inf[s][ind].l[i];
        update_p_ob(s, thisInd);
        ans+= obs_inf[s][thisInd].log_pob * w[s][thisInd];
    }
    for(int i = 0; i < num_r; i++){
        thisInd = node_inf[s][ind].r[i];
        update_p_ob(s, thisInd);
        ans+= obs_inf[s][thisInd].log_pob * w[s][thisInd];
    }
    if(ISNAN(ans)) ans = R_NegInf;
    return(ans);
}

void icm_Abst::update_etas(){
    for(int s = 0; s < n_strata; s++){
        etas[s] = covars[s] * reg_par;
        for(int i = 0; i < etas[s].size(); i++){
            etas[s][i] += intercept[s];
            expEtas[s][i] = exp(etas[s][i] );
        }
    }
}


/* Recenter the baseline cumulative hazard after the optimisation is finished.
   Executed once in ic_sp_ch() after run returns.
   This could loop internally over strata */
void icm_Abst::recenterBCH(){
    for(int s = 0; s < n_strata; s++){
        int k = baseCH[s].size();
        for(int i = 1; i < (k-1); i++){
            baseCH[s][i] += intercept[s];
        } 
    }
}

void cumhaz2p_hat(Eigen::VectorXd &ch, std::vector<double> &p){
    int k = ch.size();
    std::vector<double> S(k);
    p.resize(k-1);
    for(int i = 0; i < k; i++){ S[i] = exp(-exp(ch[i])); }
    for(int i = 0; i < (k-1); i++){ p[i] = S[i+1] - S[i]; }
}


void icm_Abst::icm_addPar(int s, std::vector<double> &delta){
    int p_k = delta.size();
    int a_k = baseCH[s].size();
    if( (p_k+2) != a_k){Rprintf("in icm_addPar, delta is not the same length as actIndex!\n");return;}
    for(int i = 0; i < p_k; i++){ baseCH[s][i+1] += delta[i]; }
}




/*      INITIALIZATION TOOLS    */
void setup_icm(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, SEXP R_strata,
                SEXP R_RegPars, icm_Abst* icm_obj){
    icm_obj->h = 0.0001;
    icm_obj->almost_inf = 1.0/icm_obj->h;

    Rcpp::IntegerVector strata(R_strata);
    int nS = strata[0];
    icm_obj->n_strata = nS;
    // check inputs consistent for # strata
    if(Rf_length(Rlind) != Rf_length(Rrind)){
        Rprintf("length of Rlind and Rrind not equal\n");
        return;
    }
    if(Rf_length(Rlind) != Rf_length(R_w)){
        Rprintf("length of Rlind and R_w not equal\n");
        return;
    }
    if(Rf_length(Rlind) != Rf_length(RCovars)){
        Rprintf("length of Rlind and RCovars not equal\n");
        return;
    }

    icm_obj->base_p_obs.resize(nS);
    icm_obj->etas.resize(nS);
    icm_obj->expEtas.resize(nS);
    icm_obj->w.resize(nS);
    icm_obj->intercept.resize(nS);
    icm_obj->covars.resize(nS);
    
    int reg_k;

    for(int s = 0; s < icm_obj->n_strata; s++){
        int n = Rf_length(VECTOR_ELT(Rlind, s));
        if(n != Rf_length(VECTOR_ELT(Rrind, s))){Rprintf("length of Rlind and Rrind not equal\n"); return;}

        icm_obj->base_p_obs[s].resize(n);
        icm_obj->etas[s].resize(n);
        icm_obj->expEtas[s].resize(n);
        icm_obj->w[s].resize(n);

        icm_obj->intercept[s] = 0.0;

        for(int i = 0; i < n; i++){
            icm_obj->etas[s][i]       = 0;
            icm_obj->expEtas[s][i]    = 1;
            icm_obj->base_p_obs[s][i] = 0;
            icm_obj->w[s][i]          = REAL(VECTOR_ELT(R_w,s))[i];
        }

        copyRmatrix_intoEigen(VECTOR_ELT(RCovars, s), icm_obj->covars[s]);
    
        int reg_k_s = icm_obj->covars[s].cols();
        if (s == 0) {
            reg_k = reg_k_s;
            if(reg_k == 0) icm_obj->hasCovars = false; else icm_obj->hasCovars = true;
        } else {
            if (reg_k != reg_k_s) {
                Rprintf("Covariates have different number of columns across strata!\n");
                return;
            }
        }
        if(reg_k_s > 0){
            if(n != icm_obj->covars[s].rows()) {Rprintf("covar rows not equal to n!\n"); return;}
        }
    }
       
    icm_obj->reg_d1.resize(reg_k);
    icm_obj->reg_d2.resize(reg_k, reg_k);
    icm_obj->reg_d3.resize(reg_k);
    icm_obj->reg_par.resize(reg_k);
    double* regParPtr = REAL(R_RegPars);
    for(int i = 0; i < reg_k; i++){ icm_obj->reg_par[i] = regParPtr[i]; }
    

    icm_obj->baseCH.resize(nS);
    icm_obj->backupCH.resize(nS);
    icm_obj->baseS.resize(nS);
    icm_obj->baseP.resize(nS);
    icm_obj->baseP_backup.resize(nS);
    
    icm_obj->obs_inf.resize(nS);
    icm_obj->node_inf.resize(nS);
    icm_obj->usedVec.resize(nS);
    icm_obj->dob_dp_both.resize(nS);
    icm_obj->dob_dp_rightOnly.resize(nS);

    icm_obj->base_p_2ndDerv.resize(nS);
    icm_obj->base_p_derv.resize(nS);
    icm_obj->base_p_derv2.resize(nS);
    icm_obj->prop_p.resize(nS);

    for(int s = 0; s < nS; s++){
        int n = Rf_length(VECTOR_ELT(Rrind, s));
        int maxInd = 0;
        for(int i = 0; i < n; i++){
            maxInd = max(maxInd, INTEGER(VECTOR_ELT(Rrind, s))[i]);
        }

        icm_obj->baseCH[s].resize(maxInd + 2);
        

        for(int i = 0; i <= maxInd; i++){ 
            icm_obj->baseCH[s][i] = R_NegInf; 
        }
        icm_obj->baseCH[s][maxInd+1] = R_PosInf;
        icm_obj->baseS[s].resize(maxInd + 2);
        icm_obj->baseS[s][0] = 1.0;
        icm_obj->baseS[s][maxInd+1] = 0;
        int this_l, this_r;

        icm_obj->obs_inf[s].resize(n);
        icm_obj->node_inf[s].resize(maxInd + 2);

        for(int i = 0; i < n; i++){
            this_l = INTEGER(VECTOR_ELT(Rlind, s))[i];
            this_r = INTEGER(VECTOR_ELT(Rrind, s))[i];
            icm_obj->obs_inf[s][i].l = this_l;
            icm_obj->obs_inf[s][i].r = this_r;
            icm_obj->node_inf[s][this_l].l.push_back(i);
            icm_obj->node_inf[s][this_r + 1].r.push_back(i);
        }

        double stepSize = -1.0/(1.0 + icm_obj->baseS[s].size() );
        double curVal = 1.0;

        for(int i = 1; i < (maxInd+1); i++){
            curVal += stepSize;
            icm_obj->baseS[s][i] = curVal;
        }
        std::vector<double> this_S = icm_obj->baseS[s];

        icm_obj->baseS_2_baseCH(s);     //TURN OFF IF WANT TO SWICH TO CH START


        icm_obj->usedVec[s].resize(n);
        for(int i = 0; i < n; i++){icm_obj->usedVec[s][i] = false;}
    }
                
    icm_obj->startGD = false;
    icm_obj->failedGA_counts = 0;
    icm_obj->iter = 0;
    icm_obj->numBaselineIts = 5;
    
}



/*      OPTIMIZATION TOOLS      */
void icm_Abst::numericBaseDervsOne(int s, int raw_ind, std::vector<double> &dvec){
    dvec.assign(3, 0.0);
    
    if(raw_ind <= 0 || raw_ind >= (baseCH[s].size()- 1)){Rprintf("warning: inappropriate choice of ind for numericBaseDervs ind = %d\n", raw_ind); return;}
    
    h = h / 25.0;
    
    baseCH[s][raw_ind] += h;
    double llk_h = par_llk(s, raw_ind);
    baseCH[s][raw_ind] -= 2*h;
    double llk_l = par_llk(s, raw_ind);
    baseCH[s][raw_ind] += h;
    double llk_st = par_llk(s, raw_ind);
    
    if(llk_l == R_NegInf){
        llk_l = llk_st;
        baseCH[s][raw_ind] += h/2.0;
        llk_st = par_llk(s, raw_ind);
        baseCH[s][raw_ind] -= h/2.0;
    }
    
    if(llk_h == R_NegInf){
        llk_h = llk_st;
        baseCH[s][raw_ind] -= h/2.0;
        llk_st = par_llk(s, raw_ind);
        baseCH[s][raw_ind] += h/2.0;
    }
    
    dvec[0] = (llk_h - llk_l)/(2.0*h);
    dvec[1] = (llk_h + llk_l - 2.0 * llk_st) / (h * h);
    dvec[2] = llk_st;
        
    if(dvec[1] == R_NegInf || ISNAN(dvec[1]) ){
        h = h/100.0;
        
        baseCH[s][raw_ind] += h;
        double llk_h = par_llk(s, raw_ind);
        baseCH[s][raw_ind] -= 2*h;
        double llk_l = par_llk(s, raw_ind);
        baseCH[s][raw_ind] += h;
        double llk_st = par_llk(s, raw_ind);
        
        dvec[0] = (llk_h - llk_l)/(2.0*h);
        dvec[1] = (llk_h + llk_l - 2.0 * llk_st) / (h * h);
        dvec[2] = llk_st;
        h *=100.0;
    }
    
    h = h * 25.0;
}

void icm_Abst::numericBaseDervsAllRaw(int s, std::vector<double> &d1, std::vector<double> &d2){
    int k = baseCH[s].size() - 2;
    d1.resize(k);
    d2.resize(k);
    
    std::vector<double> ind_dervs(3);
    for(int i = 0; i < k; i++){
        numericBaseDervsOne(s, i + 1, ind_dervs);
        d1[i] = ind_dervs[0];
        d2[i] = ind_dervs[1];
        
    }
}



void icm_Abst::analytical_dobs_dch(int s, std::vector<double> &d1, std::vector<double> &d2){
    int k = baseCH[s].size() - 2;
    d1.resize(k);
    d2.resize(k);

    for(int param = 0; param < k; ++param){
        d1[param] = 0.0;
        d2[param] = 0.0;
        
        int baseCH_idx = param + 1;

        // Left boundary
        int num_l = node_inf[s][baseCH_idx].l.size();
        for(int i = 0; i < num_l; ++i){
            int obs = node_inf[s][baseCH_idx].l[i];
            double chl = baseCH[s][obs_inf[s][obs].l];
            double chr = baseCH[s][obs_inf[s][obs].r + 1];
            double eta = etas[s][obs];

            std::vector<double> derivs(2);
            derivs = dllk_dch_i(chl, chr, eta, obs_inf[s][obs].log_pob, true);
            if(std::isnan(derivs[0]) || std::isnan(derivs[1])){
                Rcpp::Rcout << "NaN detected in analytical_dobs_dch at stratum " << s << ", observation " << obs << std::endl;
            }
            d1[param] += derivs[0];
            d2[param] += derivs[1];
        }

        // Right boundary
        int num_r = node_inf[s][baseCH_idx].r.size();
        for(int i = 0; i < num_r; ++i){
            int obs = node_inf[s][baseCH_idx].r[i];
            double chl = baseCH[s][obs_inf[s][obs].l];
            double chr = baseCH[s][obs_inf[s][obs].r + 1];
            double eta = etas[s][obs];
            std::vector<double> derivs(2);
            derivs = dllk_dch_i(chl, chr, eta, obs_inf[s][obs].log_pob, false);
            
            if(std::isnan(derivs[0]) || std::isnan(derivs[1])){
                Rcpp::Rcout << "NaN detected in analytical_dobs_dch at stratum " << s << ", observation " << obs << std::endl;
            }
            d1[param] += derivs[0];
            d2[param] += derivs[1];
        }
    }
}




void icm_Abst::icm_step(){
    for(int s = 0; s < n_strata; s++){
        icm_step_s(s);
    }
}

void icm_Abst::icm_step_s(int s){
        backupCH[s] = baseCH[s];
        double llk_st = sum_llk(s);
        
        std::vector<double> d1;
        std::vector<double> d2;
        

        if (derivMethod == 1) {
            // Use raw numeric derivatives
            numericBaseDervsAllRaw(s, d1, d2);
        } else if (derivMethod == 2) {
         // Use raw numeric derivatives
            analytical_dobs_dch(s, d1, d2);
        } else if (derivMethod == 13) {
         // Use raw numeric derivatives
         if (iter == 1) {
             numericBaseDervsAllRaw(s, d1, d2);
         } else {
            analytical_dobs_dch(s, d1, d2);
         }
         } else if (derivMethod == 14) {
         // Use raw numeric derivatives
             analytical_dobs_dch(s, d1, d2);
        
        } else if (derivMethod == 11) {
            // Use raw numeric derivatives
            numericBaseDervsAllRaw(s, d1, d2);
        } else if (derivMethod == 12) {
            // Use analytical differentiation
            analytical_dobs_dch(s, d1, d2);
        } else {

            Rcpp::Rcout  << derivMethod << "Invalid derivation method selected.\n";
            return;
        }
     
        int thisSize = d1.size();

        for(int i = 0; i < thisSize; i ++){
            if(d2[i] == R_NegInf){d2[i] = -almost_inf;}
            if(ISNAN(d2[i]))    {
            //	Rprintf("warning: d2 isnan! \n");
                baseCH[s] = backupCH[s];
            return;
        }
        if(d2[i] >= 0) {
            int sum_neg = 0;
            double sum_neg_d2s = 0.0;
            for(int j = 0; j < thisSize; j++){
                if(d2[j] < 0){
                    sum_neg++;
                    sum_neg_d2s += d2[j];
                }
            }
            double mean_neg_d2s = sum_neg_d2s / sum_neg;
            if(ISNAN(mean_neg_d2s) ){mean_neg_d2s = -1.0;}
            for(int j = 0; j < thisSize; j++){
                if(d2[j] >= 0){d2[j] = mean_neg_d2s;}
            }
        }
    }
    std::vector<double> x(d1.size());
    int x_k = x.size();
    int baseCH_k = baseCH[s].size();
    if(x_k != baseCH_k - 2){Rprintf("warning: x.size()! = actIndex.size()\n"); return;}
    thisSize = baseCH[s].size() - 2;
    for(int i = 0; i < thisSize; i++){x[i] = baseCH[s][i + 1];}
    std::vector<double> prop(d1.size());

    pavaForOptim(d1, d2, x, prop);

    icm_addPar(s, prop);
    checkCH(s);

    double llk_new = sum_llk(s);
    mult_vec(-1.0, prop);
    int tries = 0;
    while(llk_st > llk_new && tries < 5){
        tries++;
        mult_vec(0.5, prop);
        icm_addPar(s, prop);
        checkCH(s);        

        llk_new = sum_llk(s);
    }

    if(llk_new < llk_st){
        baseCH[s] = backupCH[s];
        llk_new = sum_llk(s);
        
        int numNAs = 0;
        double sumAbsProp = 0;
        for(int i = 0; i < thisSize; i++){
            if(ISNAN(prop[i])){
                numNAs++;
            }
            else{
                sumAbsProp += abs(prop[i]);
            }
        }
        mult_vec(0, prop);
    }
    maxBaseChg = 0;
    for(int i = 0; i < thisSize; i++){
        maxBaseChg = max(maxBaseChg, abs(prop[i]) );
    }
}


void icm_Abst::calcAnalyticRegDervs(Eigen::MatrixXd &hess, Eigen::VectorXd &d1){
    int k = reg_par.size();

    hess.resize(k, k);
    d1.resize(k);
    for(int i = 0; i < k; i++){
        d1[i] = 0;
        hess(i,i) = 0;
        if(useFullHess){
            for(int j = 0; j < i; j++){hess(i,j) = 0.0; hess(j,i) = 0.0;}
        }
    }

    for(int s = 0; s < n_strata; s++){
        int n = etas[s].size();
    
        Eigen::VectorXd l_cont(n);
        Eigen::VectorXd r_cont(n);
        Eigen::VectorXd totCont(n);

        Eigen::VectorXd l_cont2(n);
        Eigen::VectorXd r_cont2(n);
        Eigen::VectorXd totCont2(n);


        int lind, rind;
        double l_ch, r_ch, eta, pob, log_p;
        for(int i = 0; i < n; i++){
            l_cont[i]  = 0;
            r_cont[i]  = 0;
            l_cont2[i] = 0;
            r_cont2[i] = 0;

            lind = obs_inf[s][i].l;
            rind = obs_inf[s][i].r;
            pob  = obs_inf[s][i].pob;
            log_p = obs_inf[s][i].log_pob;
            l_ch = baseCH[s][lind];
            r_ch = baseCH[s][rind + 1];
            eta  = etas[s][i];
            if(l_ch > R_NegInf){
                l_cont[i]  = reg_d1_lnk(l_ch, eta, log_p);
                l_cont2[i] = reg_d2_lnk(l_ch, eta, log_p);
            }
            if(r_ch < R_PosInf){
                r_cont[i]  = -reg_d1_lnk(r_ch, eta, log_p);
                r_cont2[i] = -reg_d2_lnk(r_ch, eta, log_p);
            }
            totCont[i] = l_cont[i] + r_cont[i];
            totCont2[i] = l_cont2[i] + r_cont2[i] - totCont[i] * totCont[i];
        }

        double this_covar;
        double this_w;
        double this_w_covar;
        double this_totCont;
        double this_totCont2;
        for(int i = 0; i < n; i++){
            this_w = w[s][i];
            this_totCont = totCont[i];
            this_totCont2 = totCont2[i];
            for(int a = 0; a < k; a++){
                this_covar = covars[s](i,a);
                this_w_covar = this_w * this_covar;
                d1[a] += this_w_covar * this_totCont;
                if(useFullHess){
                    for(int b = 0; b < a; b++){
                        hess(a,b) += this_w_covar * covars[s](i,b) * this_totCont2;
                        hess(b,a) = hess(a,b);
                    }
                }
                hess(a,a) += this_w_covar * this_covar * this_totCont2;
            }
        }
    }
}

void icm_Abst::calcFinalRegContr(Eigen::MatrixXd &hess, Eigen::VectorXd &d1, Eigen::VectorXd &d3){
    int k = reg_par.size();

    hess.resize(k, k);
    d1.resize(k);
    d3.resize(k);
    for(int i = 0; i < k; i++){
        d1[i] = 0;
        d3[i] = 0;
        hess(i,i) = 0;
        if(useFullHess){
            for(int j = 0; j < i; j++){hess(i,j) = 0.0; hess(j,i) = 0.0;}
        }
    }

    for(int s = 0; s < n_strata; s++){
        int n = etas[s].size();
    
        Eigen::VectorXd l_cont(n);
        Eigen::VectorXd r_cont(n);
        Eigen::VectorXd totCont(n);

        Eigen::VectorXd l_cont2(n);
        Eigen::VectorXd r_cont2(n);
        Eigen::VectorXd totCont2(n);

        Eigen::VectorXd l_cont3(n);
        Eigen::VectorXd r_cont3(n);
        Eigen::VectorXd totCont3(n);

        int lind, rind;
        double l_ch, r_ch, eta, pob, log_p;
        for(int i = 0; i < n; i++){
            l_cont[i]  = 0;
            r_cont[i]  = 0;
            l_cont2[i] = 0;
            r_cont2[i] = 0;
            l_cont3[i] = 0;
            r_cont3[i] = 0;

            lind = obs_inf[s][i].l;
            rind = obs_inf[s][i].r;
            pob  = obs_inf[s][i].pob;
            log_p = obs_inf[s][i].log_pob;
            l_ch = baseCH[s][lind];
            r_ch = baseCH[s][rind + 1];
            eta  = etas[s][i];
            if(l_ch > R_NegInf){
                l_cont[i]  = reg_d1_lnk(l_ch, eta, log_p);
                l_cont2[i] = reg_d2_lnk(l_ch, eta, log_p);
                l_cont3[i] = reg_d3_lnk(l_ch, eta, log_p);
            }
            if(r_ch < R_PosInf){
                r_cont[i]  = -reg_d1_lnk(r_ch, eta, log_p);
                r_cont2[i] = -reg_d2_lnk(r_ch, eta, log_p);
                r_cont3[i] = -reg_d3_lnk(r_ch, eta, log_p);
            }
            totCont[i] = l_cont[i] + r_cont[i];
            totCont2[i] = l_cont2[i] + r_cont2[i] - totCont[i] * totCont[i];
            totCont3[i] = l_cont3[i] + r_cont3[i]
                          - 3.0 * totCont[i] * totCont2[i]
                          - totCont[i] * totCont[i] * totCont[i];
        }

        double this_covar;
        double this_w;
        double this_w_covar;
        double this_totCont;
        double this_totCont2;
        double this_totCont3;
        for(int i = 0; i < n; i++){
            this_w = w[s][i];
            this_totCont = totCont[i];
            this_totCont2 = totCont2[i];
            this_totCont3 = totCont3[i];
            // Compute ||z_i||^2 for cross-term d3 accumulation
            double z_sq_sum = 0.0;
            for(int b = 0; b < k; b++){
                z_sq_sum += covars[s](i,b) * covars[s](i,b);
            }
            for(int a = 0; a < k; a++){
                this_covar = covars[s](i,a);
                this_w_covar = this_w * this_covar;
                d1[a] += this_w_covar * this_totCont;
                d3[a] += this_w_covar * z_sq_sum * this_totCont3;
                if(useFullHess){
                    for(int b = 0; b < a; b++){
                        hess(a,b) += this_w_covar * covars[s](i,b) * this_totCont2;
                        hess(b,a) = hess(a,b);
                    }
                }
                hess(a,a) += this_w_covar * this_covar * this_totCont2;
            }
        }
    }
}
 
 
void icm_Abst::covar_nr_step(){
    int k = reg_par.size();
    calcAnalyticRegDervs(reg_d2, reg_d1);
    double lk_0 = sum_llk_all();

/*    for(int i = 0; i < k; i++){
        if(reg_d2[i] >= -0.0000001 || ISNAN(reg_d2[i])){
            reg_d2[i] = -100.00;
        }
        if(ISNAN(reg_d1[i]) ){reg_d1[i] = 0;}
    }       */
    
    propVec.resize(k);
    if(useFullHess){
      propVec = -reg_d2.fullPivLu().solve(reg_d1);
      
      double err = (reg_d2*propVec + reg_d1).norm() / reg_d1.norm();
      if(err > .001){
        for(int i = 0; i < k; i++){
          propVec[i] = 0;
          if(reg_d2(i,i) < 0)   propVec[i] = -reg_d1[i] / reg_d2(i,i);
          else propVec[i] = signVal(reg_d1[i]) * 0.01;
          if(ISNAN(propVec[i])) propVec[i] = 0;
        }
      }
      
    }
    else{for(int i = 0; i < k; i++){propVec[i] = -reg_d1[i]/reg_d2(i,i);}}
    int tries = 0;
    reg_par += propVec;
    propVec *= -1;
    update_etas();
    double lk_new = sum_llk_all();
    while(lk_new < lk_0 && tries < 10){
        tries++;
        propVec *= 0.5;
        reg_par += propVec;
        update_etas();
        lk_new = sum_llk_all();
    }
}


/*      CALLING ALGORITHM FROM R     */
// [[Rcpp::export]]
SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType,
              SEXP R_w, SEXP R_strata, SEXP R_use_GD, SEXP R_maxiter,
              SEXP R_baselineUpdates, SEXP R_useFullHess, SEXP R_updateCovars,
              SEXP R_initialRegVals, SEXP R_derivMethod,
              SEXP R_baselineStart) {
    icm_Abst* optObj;
    bool useGD = LOGICAL(R_use_GD)[0] == TRUE;
    
    if(INTEGER(fitType)[0] == 1){
        optObj = new icm_ph;
    }
    else if(INTEGER(fitType)[0] == 2){
        optObj = new icm_po;
    }
    else { Rprintf("fit type not supported\n");return(R_NilValue);}
    optObj->updateCovars = LOGICAL(R_updateCovars)[0] == TRUE;
    double llk_new = R_NegInf;
    
    Rcpp::IntegerVector derivOptions(R_derivMethod);
    int restart = 0;
    for (int derivMethod : derivOptions) {
        if (restart > 0) {
            Rprintf("Error found (see warnings). Restarting optimization with derivative method %d\n", derivMethod);
        }
        try {
            setup_icm(Rlind, Rrind, Rcovars, R_w, R_strata, R_initialRegVals, optObj);

            // Warm-start baseline from previous fit if provided
            if (R_baselineStart != R_NilValue) {
                for (int s = 0; s < optObj->n_strata; s++) {
                    SEXP s_vec = VECTOR_ELT(R_baselineStart, s);
                    int k = Rf_length(s_vec);
                    double* s_ptr = REAL(s_vec);
                    optObj->baseS[s].resize(k);
                    for (int i = 0; i < k; i++) {
                        optObj->baseS[s][i] = s_ptr[i];
                    }
                    optObj->baseS_2_baseCH(s);
                }
            }

            optObj->useFullHess = LOGICAL(R_useFullHess)[0] == TRUE;
            optObj->derivMethod = derivMethod;//INTEGER(R_derivMethod)[0];
    
            double tol = pow(10.0, -10.0);
            int maxIter = INTEGER(R_maxiter)[0];
            int baselineUpdates = INTEGER(R_baselineUpdates)[0];
    
            llk_new = optObj->run(maxIter, tol, useGD, baselineUpdates);
            if (llk_new == R_NegInf) {
                throw std::runtime_error("Log-likelihood is -Inf after optimization.");
            }
            if (restart > 0) {
                Rprintf("Optimization successfully completed with derivative method %d.", derivMethod);
            }
            break;
        } catch (...) {
           warning("Error encountered with derivative method %d. \n", derivMethod);
           restart++;
        }
    }
    if (llk_new == R_NegInf) {
        throw Rcpp::exception("Final log-likelihood is -Inf.");
    }
    
    std::vector<std::vector<double>> p_hat; 
    p_hat.resize(optObj->n_strata);
    optObj->recenterBCH();
    
    for(int s = 0; s < optObj->n_strata; s++){
        cumhaz2p_hat(optObj->baseCH[s], p_hat[s]);
    }
    
    
    
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 8));
    SEXP R_pans = PROTECT(Rf_allocVector(VECSXP,p_hat.size()));
    SEXP R_coef = PROTECT(Rf_allocVector(REALSXP, optObj->reg_par.size()));
    SEXP R_fnl_llk = PROTECT(Rf_allocVector(REALSXP, 1));
    SEXP R_its = PROTECT(Rf_allocVector(REALSXP, 1));
    SEXP R_score = PROTECT(Rf_allocVector(REALSXP, optObj->reg_par.size()));
    SEXP R_hessian = PROTECT(Rf_allocMatrix(REALSXP, optObj->reg_par.size(), optObj->reg_par.size()));
    SEXP R_d3 = PROTECT(Rf_allocVector(REALSXP, optObj->reg_par.size()));
    SEXP R_subj_llk = PROTECT(Rf_allocVector(VECSXP, optObj->n_strata));

    for (size_t i = 0; i < p_hat.size(); ++i) {
        const std::vector<double>& inner = p_hat[i];
        SEXP inner_vec = PROTECT(Rf_allocVector(REALSXP, inner.size()));
        std::copy(inner.begin(), inner.end(), REAL(inner_vec));
        SET_VECTOR_ELT(R_pans, i, inner_vec);
        UNPROTECT(1); // unprotect inner_vec after assigning to list
    }

    for(int i = 0; i < optObj->reg_par.size(); i++){
        REAL(R_coef)[i] = optObj->reg_par[i];
        REAL(R_score)[i] = optObj->reg_d1[i];
        REAL(R_d3)[i] = optObj->reg_d3[i];
    }

    for (int s = 0; s < optObj->n_strata; s++) {
        int n = optObj->obs_inf[s].size();
        SEXP R_slk_s = PROTECT(Rf_allocVector(REALSXP, n));
        for (int i = 0; i < n; i++) {
            REAL(R_slk_s)[i] = optObj->obs_inf[s][i].log_pob;
        }
        SET_VECTOR_ELT(R_subj_llk, s, R_slk_s);
        UNPROTECT(1);
    }
    
    
    // Copy Hessian matrix (column-major order for R)
    for(int i = 0; i < optObj->reg_par.size(); i++){
        for(int j = 0; j < optObj->reg_par.size(); j++){
            REAL(R_hessian)[i + j * optObj->reg_par.size()] = optObj->reg_d2(i, j);
        }
    }
    
    REAL(R_fnl_llk)[0] = llk_new;
    REAL(R_its)[0] = optObj->iter;
    
    SET_VECTOR_ELT(ans, 0, R_pans);
    SET_VECTOR_ELT(ans, 1, R_coef);
    SET_VECTOR_ELT(ans, 2, R_fnl_llk);
    SET_VECTOR_ELT(ans, 3, R_its);
    SET_VECTOR_ELT(ans, 4, R_score);
    SET_VECTOR_ELT(ans, 5, R_hessian);
    SET_VECTOR_ELT(ans, 6, R_d3);
    SET_VECTOR_ELT(ans, 7, R_subj_llk);
    
    UNPROTECT(9);

    
    if(INTEGER(fitType)[0] == 1){
        icm_ph* deleteObj = static_cast<icm_ph*>(optObj);
        delete deleteObj;
    }
    else if(INTEGER(fitType)[0] == 2){
        icm_po* deleteObj = static_cast<icm_po*>(optObj);
        delete deleteObj;
    }
    
    
    return(ans);

}

void icm_Abst::checkCH(int s){
    int k = baseCH[s].size();
    for(int i = 1; i < k; i++){
        if(baseCH[s][i] < baseCH[s][i-1]){
            baseCH[s][i] = baseCH[s][i-1]; 
        }
    }
}

double icm_Abst::run(int maxIter, double tol, bool useGD, int baselineUpdates){
    iter = 0;
    bool metOnce = false;
    double llk_old = R_NegInf;
 
    bool regNon0 = false;
    int reg_k = reg_par.size();
    for(int i = 0; i < reg_k; i++){
        if(reg_par[i] != 0 ){ regNon0 = true; } 
    }
    
    if(regNon0){update_etas();}
    double llk_new = sum_llk_all(); // global log-likelihood

    if(regNon0){
        if(hasCovars){stablizeBCH();}
        if(useGD){ gradientDescent_step();}
        icm_step();
        if(useGD){ gradientDescent_step();}		
        icm_step();
        llk_new = sum_llk_all(); // recompute after warm-up steps
    }
 
    while(iter < maxIter && (llk_new - llk_old) > tol){
        iter++;
        llk_old = llk_new;
        
        //Rprintf("%.7f\n", llk_old);
        if(hasCovars && updateCovars){ covar_nr_step(); }

        for(int i = 0; i < baselineUpdates; i++)  {
            if(hasCovars){stablizeBCH();}
            icm_step();
            if(useGD){ gradientDescent_step(); }
        }
            
        llk_new = sum_llk_all();
        if(llk_new - llk_old > tol){metOnce = false;}
        if(metOnce == false){
            if(llk_new - llk_old <= tol){
                metOnce = true;
                llk_old = llk_old - 2 * tol;
            }
        }
 
       if((llk_new - llk_old) < -0.001 ){
           Rprintf("warning: likelihood decreased! difference = %f\n", llk_new - llk_old);
       }
    }
    
    // Update final derivatives for return to R
    calcFinalRegContr(reg_d2, reg_d1, reg_d3);
    
    return(llk_new);
}


/*      GETTING MAXIMAL INTERSECTIONS       */
// [[Rcpp::export]]
SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals){
    //NOTE: R_AllVals MUST be sorted!!
    int k = LENGTH(R_AllVals);
    std::vector<double> mi_l;
    std::vector<double> mi_r;
    
    mi_l.reserve(k);
    mi_r.reserve(k);
    
    bool foundLeft = false;
    double last_left = R_NegInf;
    
    double* c_AllVals = REAL(R_AllVals);
    
    for(int i = 0; i < k; i++){
        if(!foundLeft)                      foundLeft = LOGICAL(isL)[i] == TRUE;
        if(LOGICAL(isL)[i] == TRUE)         last_left = c_AllVals[i];
        if(foundLeft){
            if(LOGICAL(isR)[i] == TRUE){
                mi_l.push_back(last_left);
                mi_r.push_back(c_AllVals[i]);
                foundLeft = false;
            }
        }
    }
    int tbulls = mi_l.size();
    
    int n = LENGTH(lVals);
    SEXP l_ind = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_ind = PROTECT(Rf_allocVector(INTSXP, n));
    
    int* cl_ind = INTEGER(l_ind);
    int* cr_ind = INTEGER(r_ind);
    double* clVals = REAL(lVals);
    double* crVals = REAL(rVals);
    
    double this_Lval, this_Rval;
    
    for(int i = 0; i < n; i++){
        this_Lval = clVals[i];
        cl_ind[i] = findSurroundingVals(this_Lval, mi_l, mi_r, true);
        this_Rval = crVals[i];
        cr_ind[i] = findSurroundingVals(this_Rval, mi_l, mi_r, false);
     }
    
    
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP Rl_mi = PROTECT(Rf_allocVector(REALSXP, tbulls));
    SEXP Rr_mi = PROTECT(Rf_allocVector(REALSXP, tbulls));
    for(int i = 0; i < tbulls; i++){
        REAL(Rl_mi)[i] = mi_l[i];
        REAL(Rr_mi)[i] = mi_r[i];
    }
    SET_VECTOR_ELT(ans, 0, l_ind);
    SET_VECTOR_ELT(ans, 1, r_ind);
    SET_VECTOR_ELT(ans, 2, Rl_mi);
    SET_VECTOR_ELT(ans, 3, Rr_mi);
    UNPROTECT(5);
    return(ans);
}
