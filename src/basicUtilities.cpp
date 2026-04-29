//
//  basicUtilities.cpp
//  
//
//  Created by Cliff Anderson Bergman on 5/16/15.
//
//
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <RcppEigen.h>

using namespace std;
using namespace Rcpp;
using namespace Eigen;

#include "utilities.h"

#define SMALLNUMBER 0.00000000001;

double max(double a, double b){
    if(a > b) return(a);
    return(b);
}

double max(int a, int b){
    if(a > b) return(a);
    return(b);
}

template <typename T>
int getMaxIndex(T &v){
    double mVal = R_NegInf;
    int mInd = 0;
    int k = v.size();
    for(int i = 0; i < k; i++){
        if(v[i] == max(mVal, v[i])) {
            mVal = v[i];
            mInd = i;
        }
    }
    return(mInd);
}

double min(double a, double b){
    if(a < b) return (a);
    return(b);
}

void mult_vec(double a, std::vector<double> &vec){
    int thisSize = vec.size();
    for(int i = 0; i < thisSize; i++)
        vec[i] *= a;
}

void add_vec(double lambda, std::vector<double> &a, std::vector<double> &vec){
    int thisSize = vec.size();
    int thisSize2 = a.size();
    if(thisSize != thisSize2){
         Rprintf("warning: sizes do not match in add_vec\n");
             return;
    }
    for(int i = 0; i < thisSize; i++)
        vec[i] += a[i] * lambda;
}

double signVal(double x){
    if(x > 0) return 1.0;
    return -1.0;
}

void copyRmatrix_intoEigen(SEXP r_mat, Eigen::MatrixXd &e_mat){
    SEXP Rdims = Rf_getAttrib(r_mat, R_DimSymbol);
    PROTECT(Rdims);
    int nRows = INTEGER(Rdims)[0];
    int nCols = INTEGER(Rdims)[1];
    
    e_mat.resize(nRows, nCols);
    for(int i = 0; i < nRows; i++){
        for(int j = 0; j < nCols; j++)
            e_mat(i, j) = REAL(r_mat)[i + j*nRows];
    }
    UNPROTECT(1);
}

void pavaForOptim(std::vector<double> &d1, std::vector<double> &d2, std::vector<double> &x, std::vector<double> &prop_delta){
    int k = d1.size();
    int d2_size = d2.size();
    int x_size = x.size();
    if(k != d2_size || k!= x_size){ Rprintf("incorrect sizes provided to pavaForOptim\n"); return;}
    prop_delta.resize(k);
    std::vector<double> y(k);
    std::vector<double> w(k);
    
    for(int i = 0; i < k; i++){
        y[i] = -d1[i]/d2[i] + x[i];
        w[i] = d2[i]/2;
    }
    int k_sign = k;
    pava( &y[0], &w[0], &k_sign );
    for(int i = 0; i < k; i++){
        prop_delta[i] = y[i] - x[i];
    }
}



// Statistical utilities

double directional_derv(const std::vector<double> &derv, const std::vector<double> &delta){
    int k = derv.size();
    int k2 = delta.size();
    if(k != k2){
        Rprintf("warning: sizes don't match in directional_derv\n");
        return(0.0);
    }
     double ans = 0.0;
    for(int i = 0; i < k; i++){
        ans += derv[i] * delta[i];
    }
    return(ans);
}


void makeUnitVector(std::vector<double> &v){
    double sum = 0;
    int k = v.size();
    for(int i = 0; i < k; i++){
        sum += fabs(v[i]);
    }
    for(int i = 0; i < k; i++){
        v[i] = v[i]/sum;
    }
}


int isValueInInterval(double val, double l, double r){
	if(val < l) return(-1);
	if(val > r) return(1);
	return(0);
}

int isValueInInterval(double val, int ind, 
					  std::vector<double>& lvec, std::vector<double>& rvec){
	return(isValueInInterval(val, lvec[ind], rvec[ind]));					  
}

int findSurroundingVals(double val, std::vector<double>& leftVec,
						std::vector<double>& rightVec, bool isLeft){
	
	int a = 0;
	int b = leftVec.size()-1;
	if(b == 0){return(0);}
	if(isValueInInterval(val, R_NegInf, rightVec[0]) == 0) return(0);
	if(isValueInInterval(val, leftVec[b], R_PosInf) == 0) return(b);
	
/*	a++;
	b--;	*/
	
	int maxTries = b;
	
	int propInd = (a + b)/2;
	int tries = 0;
	int testVal;
	while( b - a > 1 && tries < maxTries){
		tries++;
		propInd = (a + b)/2;
		testVal = isValueInInterval(val, propInd, leftVec, rightVec);
		if(testVal == 0){ return(propInd);}
		if(testVal == -1){ b = propInd; }
		else{ a = propInd; }
	}
	if(a == b){
		Rprintf("this is very surprising... a = %d, size = %zu\n", a, leftVec.size());
		return(a);
	}
	if( isLeft ) return(b); 
	return(a);
}
