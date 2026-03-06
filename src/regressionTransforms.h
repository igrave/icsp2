//
//  regressionTransforms.h
//  
//
//  Created by Cliff Anderson Bergman on 11/17/15.
//
//

#ifndef ____regressionTransforms__
#define ____regressionTransforms__

#include <Rcpp.h>

class condProbCal{
public:
    SEXP* baselineInfo;
    double (*getBaseSurv)(double q, SEXP bli);
    double (*transformSurv)(double s, double nu);
    double (*getBaseQ)(double s, SEXP bli);
    double (*transform_p)(double q, double nu);
    std::vector<double> preppedParams;
    condProbCal(SEXP regType, SEXP baseType, SEXP bli);
    bool isOK;
};


// SURVIVAL TRANSFORMS
inline double propHazTrans(double s, double nu){
	if(s == 0 || s == 1) return(s);
    double ans = pow(s, nu);
    return( ans );
}

inline double propOddsTrans(double s, double nu){
	if(s == 0 || s == 1) return(s);
    double ans;
    double prod = s * nu;
    ans = prod/(prod - s + 1);
    return(ans);
}

inline double noTrans(double s, double nu){ return(s); }

inline double transform_p_ph(double p, double nu){
	if(p == 0 || p == 1) return(p);
    double log_s = log(1.0 - p);
    double log_trans = log_s / nu;
    return(1.0 - exp(log_trans));
}
inline double transform_p_po(double p, double nu){
	if(p == 1 || p == 0) return(p);
    double s = 1.0 - p;
    return(1.0 - s * (1/nu) / (s * 1/nu - s +1));
}

inline double transform_p_none(double p, double nu){return(p);}


// BASELINE MODELS
double getNonParSurv(double q, SEXP bli);
double getNonParQ(double q, SEXP bli);




SEXP s_regTrans(SEXP times, SEXP etas,
                    SEXP bli, SEXP regType, SEXP baseType);
SEXP q_regTrans(SEXP q, SEXP etas,
                    SEXP bli, SEXP regType, SEXP baseType);

#endif /* defined(____regressionTransforms__) */