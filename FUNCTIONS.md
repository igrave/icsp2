# C++ Functions Reference

| Source File | Function Signature | Called By Other Functions |
|------------|-------------------|---------------------------|
| basicUtilities.cpp | `double max(double a, double b)` | Yes (getMaxIndex, max overload, getMaxScaleSize, addIfNeeded, pavaForOptim) |
| basicUtilities.cpp | `double max(int a, int b)` | Yes (used in ic_sp_ch.cpp) |
| basicUtilities.cpp | `template <typename T> int getMaxIndex(T &v)` | No |
| basicUtilities.cpp | `double min(double a, double b)` | Yes (getMaxScaleSize, gradientDescent_step, numeric_dobs_dp) |
| basicUtilities.cpp | `double max(Eigen::VectorXd v)` | No |
| basicUtilities.cpp | `void mult_vec(double a, std::vector<double> &vec)` | Yes (icm_step_s, gradientDescent_step) |
| basicUtilities.cpp | `void add_vec(double a, std::vector<double> &vec)` | No |
| basicUtilities.cpp | `void add_vec(std::vector<double> &a, std::vector<double> &vec)` | No |
| basicUtilities.cpp | `void add_vec(double lambda, std::vector<double> &a, std::vector<double> &vec)` | Yes (gradientDescent_step) |
| basicUtilities.cpp | `double signVal(double x)` | No |
| basicUtilities.cpp | `void copyRmatrix_intoEigen(SEXP r_mat, Eigen::MatrixXd &e_mat)` | Yes (setup_icm) |
| basicUtilities.cpp | `Rcpp::NumericMatrix eigen2RMat(Eigen::MatrixXd &e_mat)` | No |
| basicUtilities.cpp | `void Rvec2eigen(SEXP r_vec, Eigen::VectorXd &e_vec)` | No |
| basicUtilities.cpp | `Rcpp::NumericVector eigen2RVec(Eigen::VectorXd &e_vec)` | No |
| basicUtilities.cpp | `double ic_dloglogistic(double x, double a, double b)` | No |
| basicUtilities.cpp | `double ic_ploglogistic(double x, double a, double b)` | No |
| basicUtilities.cpp | `double ic_qloglogistic(double p, double a, double b)` | No |
| basicUtilities.cpp | `double ic_dlnorm(double x, double mu, double s)` | Yes (ic_dgeneralgamma) |
| basicUtilities.cpp | `double ic_plnorm(double x, double mu, double s)` | Yes (ic_pgeneralgamma) |
| basicUtilities.cpp | `double ic_dgeneralgamma(double x, double mu, double s, double Q)` | No |
| basicUtilities.cpp | `double ic_pgeneralgamma(double q, double mu, double s, double Q)` | No |
| basicUtilities.cpp | `double ic_qgeneralgamma(double p, double mu, double s, double Q)` | No |
| basicUtilities.cpp | `void pavaForOptim(std::vector<double> &d1, std::vector<double> &d2, std::vector<double> &x, std::vector<double> &prop_delta)` | Yes (icm_step_s) |
| basicUtilities.cpp | `void pavaForOptim(Eigen::VectorXd &d1, Eigen::VectorXd &d2, Eigen::VectorXd &x, Eigen::VectorXd &prop_delta)` | No |
| basicUtilities.cpp | `void addIfNeeded(std::vector<int> &points, int l, int r, int max)` | No |
| basicUtilities.cpp | `SEXP pava(SEXP R_d1, SEXP R_d2, SEXP R_x)` | No (R callable) |
| basicUtilities.cpp | `SEXP fastNumericInsert(SEXP newVals, SEXP target, SEXP indices)` | No (R callable) |
| basicUtilities.cpp | `void SEXP2doubleVec(SEXP R_vec, std::vector<double> &c_vec)` | No |
| basicUtilities.cpp | `void doubleVec2SEXP(std::vector<double> &c_vec, SEXP R_vec)` | No |
| basicUtilities.cpp | `void SEXPIndex2intIndex(SEXP R_Inds, std::vector<int> &c_inds)` | No |
| basicUtilities.cpp | `void indexVec2SEXP(std::vector<int> &c_vec, SEXP R_vec)` | No |
| basicUtilities.cpp | `void findIndexDiffs(std::vector<int> &in1, std::vector<int> &in2, std::vector<int> &in1not2, std::vector<int> &in2not1)` | No |
| basicUtilities.cpp | `void drop_index(int d_ind, std::vector<int> &indVec)` | No |
| basicUtilities.cpp | `void add_index(int a_ind, std::vector<int> &indVec)` | No |
| basicUtilities.cpp | `std::vector<int> getSEXP_MatDims(SEXP R_mat)` | No |
| basicUtilities.cpp | `void getPosNegIndices(std::vector<double> &vals, std::vector<int> &isPos, std::vector<int> &isNeg)` | No |
| basicUtilities.cpp | `void getRelValIndices(double relVal, std::vector<double> &vals, std::vector<int> &subIndex, std::vector<int> &above, std::vector<int> &below, int *max, int *min)` | No |
| basicUtilities.cpp | `double directional_derv(std::vector<double> &derv, std::vector<double> &delta)` | Yes (gradientDescent_step) |
| basicUtilities.cpp | `void makeUnitVector(std::vector<double> &v)` | Yes (gradientDescent_step) |
| basicUtilities.cpp | `void getUniqInts(int i1, int i2, std::vector<int> &uniqInts, std::vector<std::vector<int>> &vec_vec, std::vector<bool> &usedVec)` | No |
| basicUtilities.cpp | `int isValueInInterval(double val, double l, double r)` | Yes (isValueInInterval overload, findSurroundingVals) |
| basicUtilities.cpp | `int isValueInInterval(double val, int ind, std::vector<double>& lvec, std::vector<double>& rvec)` | Yes (findSurroundingVals) |
| basicUtilities.cpp | `int findSurroundingVals(double val, std::vector<double>& leftVec, std::vector<double>& rightVec, bool isLeft)` | Yes (findMI) |
| basicUtilities.cpp | `Eigen::MatrixXd xtx(Eigen::MatrixXd &x)` | No |
| basicUtilities.cpp | `Eigen::MatrixXd xtx(Eigen::MatrixXd &x, int row_start, int row_end)` | No |
| basicUtilities.cpp | `Eigen::MatrixXd copyRows(Eigen::MatrixXd &x, int row_start, int row_end)` | No |
| ic_sp_ch.cpp | `void icm_Abst::update_p_ob(int s, int i)` | Yes (sum_llk, par_llk) |
| ic_sp_ch.cpp | `double icm_Abst::sum_llk(int s)` | Yes (sum_llk_all, icm_step_s, run, llk_from_p, gradientDescent_step) |
| ic_sp_ch.cpp | `double icm_Abst::sum_llk_all()` | Yes (covar_nr_step, run) |
| ic_sp_ch.cpp | `double icm_Abst::par_llk(int s, int ind)` | Yes (numericBaseDervsOne) |
| ic_sp_ch.cpp | `void icm_Abst::update_etas()` | Yes (covar_nr_step, stablizeBCH, run) |
| ic_sp_ch.cpp | `void icm_Abst::recenterBCH()` | Yes (ic_sp_ch) |
| ic_sp_ch.cpp | `void cumhaz2p_hat(Eigen::VectorXd &ch, std::vector<double> &p)` | Yes (ic_sp_ch) |
| ic_sp_ch.cpp | `void icm_Abst::icm_addPar(int s, std::vector<double> &delta)` | Yes (icm_step_s) |
| ic_sp_ch.cpp | `void setup_icm(SEXP Rlind, SEXP Rrind, SEXP RCovars, SEXP R_w, SEXP R_strata, SEXP R_RegPars, icm_Abst* icm_obj)` | Yes (ic_sp_ch) |
| ic_sp_ch.cpp | `void icm_Abst::numericBaseDervsOne(int s, int raw_ind, std::vector<double> &dvec)` | Yes (numericBaseDervsAllAct, numericBaseDervsAllRaw) |
| ic_sp_ch.cpp | `void icm_Abst::numericBaseDervsAllAct(int s, std::vector<double> &d1, std::vector<double> &d2)` | No |
| ic_sp_ch.cpp | `void icm_Abst::numericBaseDervsAllRaw(int s, std::vector<double> &d1, std::vector<double> &d2)` | Yes (icm_step_s) |
| ic_sp_ch.cpp | `void icm_Abst::analytical_dobs_dch(int s, std::vector<double> &d1, std::vector<double> &d2)` | Yes (icm_step_s) |
| ic_sp_ch.cpp | `void icm_Abst::icm_step()` | Yes (run) |
| ic_sp_ch.cpp | `void icm_Abst::icm_step_s(int s)` | Yes (icm_step) |
| ic_sp_ch.cpp | `void icm_Abst::calcAnalyticRegDervs(Eigen::MatrixXd &hess, Eigen::VectorXd &d1)` | Yes (covar_nr_step) |
| ic_sp_ch.cpp | `void icm_Abst::covar_nr_step()` | Yes (run) |
| ic_sp_ch.cpp | `SEXP ic_sp_ch(SEXP Rlind, SEXP Rrind, SEXP Rcovars, SEXP fitType, SEXP R_w, SEXP R_strata, SEXP R_use_GD, SEXP R_maxiter, SEXP R_baselineUpdates, SEXP R_useFullHess, SEXP R_updateCovars, SEXP R_initialRegVals, SEXP R_derivMethod)` | No (R callable) |
| ic_sp_ch.cpp | `void icm_Abst::checkCH(int s)` | Yes (icm_step_s) |
| ic_sp_ch.cpp | `double icm_Abst::run(int maxIter, double tol, bool useGD, int baselineUpdates)` | Yes (ic_sp_ch) |
| ic_sp_ch.cpp | `SEXP findMI(SEXP R_AllVals, SEXP isL, SEXP isR, SEXP lVals, SEXP rVals)` | No (R callable) |
| ic_sp_ch.h | `double icm_ph::basHaz2CondS(double ch, double eta)` | Yes (update_p_ob) |
| ic_sp_ch.h | `double icm_ph::baseS2CondS(double s, double eta)` | Yes (cal_log_obs) |
| ic_sp_ch.h | `double icm_ph::base_d1_contr(double ch, double pob, double eta)` | No |
| ic_sp_ch.h | `double icm_ph::reg_d1_lnk(double ch, double xb, double log_p)` | Yes (calcAnalyticRegDervs) |
| ic_sp_ch.h | `double icm_ph::reg_d2_lnk(double ch, double xb, double log_p)` | Yes (calcAnalyticRegDervs) |
| ic_sp_ch.h | `void icm_ph::stablizeBCH()` | Yes (run) |
| ic_sp_ch.h | `double icm_ph::dllk_dp_i(double s_l, double s_r, double eta, double pob, bool left, bool right)` | Yes (analytical_dobs_dp) |
| ic_sp_ch.h | `std::vector<double> icm_ph::dllk_dch_i(double ch_l, double ch_r, double eta, double pob, bool left)` | Yes (analytical_dobs_dch) |
| ic_sp_ch.h | `double icm_po::basHaz2CondS(double ch, double eta)` | Yes (update_p_ob) |
| ic_sp_ch.h | `double icm_po::baseS2CondS(double s, double eta)` | Yes (cal_log_obs) |
| ic_sp_ch.h | `double icm_po::base_d1_contr(double ch, double pob, double eta)` | No |
| ic_sp_ch.h | `double icm_po::reg_d1_lnk(double ch, double xb, double log_p)` | Yes (calcAnalyticRegDervs) |
| ic_sp_ch.h | `double icm_po::reg_d2_lnk(double ch, double xb, double log_p)` | Yes (calcAnalyticRegDervs) |
| ic_sp_ch.h | `void icm_po::stablizeBCH()` | Yes (run) |
| ic_sp_ch.h | `double icm_po::dllk_dp_i(double s_l, double s_r, double eta, double pob, bool left, bool right)` | Yes (analytical_dobs_dp) |
| ic_sp_ch.h | `std::vector<double> icm_po::dllk_dch_i(double ch_l, double ch_r, double eta, double pob, bool left)` | Yes (analytical_dobs_dch) |
| ic_sp_gradDescent.cpp | `void icm_Abst::baseCH_2_baseS(int s)` | Yes (gradientDescent_step, llk_from_p) |
| ic_sp_gradDescent.cpp | `void icm_Abst::baseS_2_baseP(int s)` | Yes (gradientDescent_step, llk_from_p) |
| ic_sp_gradDescent.cpp | `void icm_Abst::baseP_2_baseS(int s)` | Yes (llk_from_p, numeric_dobs_dp) |
| ic_sp_gradDescent.cpp | `void icm_Abst::baseS_2_baseCH(int s)` | Yes (llk_from_p, setup_icm) |
| ic_sp_gradDescent.cpp | `double icm_Abst::llk_from_p(int s)` | Yes (gradientDescent_step) |
| ic_sp_gradDescent.cpp | `double icm_Abst::getMaxScaleSize(std::vector<double> &p, std::vector<double> &prop_p)` | Yes (gradientDescent_step) |
| ic_sp_gradDescent.cpp | `void icm_Abst::gradientDescent_step()` | Yes (run) |
| ic_sp_gradDescent.cpp | `double icm_Abst::cal_log_obs(double s1, double s2, double eta)` | Yes (numeric_dobs_dp, analytical_dobs_dp) |
| ic_sp_gradDescent.cpp | `void icm_Abst::numeric_dobs_dp(int s, bool forGA)` | Yes (gradientDescent_step) |
| ic_sp_gradDescent.cpp | `void icm_Abst::analytical_dobs_dp(int s)` | Yes (gradientDescent_step) |
| myPAVAalgorithm.cpp | `void weighted_pool(double *y, double *w, int start, int stop)` | Yes (weighted_pava) |
| myPAVAalgorithm.cpp | `void weighted_pava(double *y, double *w, int *numberParameters)` | Yes (pava) |
| myPAVAalgorithm.cpp | `void pava(double *y, double *w, int *np)` | Yes (pavaForOptim) |
| icsp2-package.cpp | None (empty file, just includes) | N/A |

---

## R Functions

| Source File | Function Signature | Called By Other Functions |
|------------|-------------------|---------------------------|
| data.R | None (data documentation only) | N/A |
| ic_sp.R | None (commented out legacy code) | N/A |
| ic_sp2.R | `ic_sp_settings(useGA = TRUE, maxIter = 10000, baseUpdates = 5, regStart = NULL, derivMethod = c(12, 1), updateCovars = TRUE)` | Yes (ic_sp2) |
| ic_sp2.R | `ic_sp2(formula, data, weights, subset, na.action, B = c(0, 1), settings = ic_sp_settings(), model = c("ph", "po"))` | No (user-facing) |
| ic_sp2.R | `ic_sp_ph` | No (alias for ic_sp2) |
| ic_sp2.R | `ic_sp_po` | No (alias for ic_sp2) |
| ic_sp2.R | `.fit_ic_sp(x, y, weights, strata, model_type, other_info)` | Yes (ic_sp2, vcov.ic_po, getBS_coef) |
| ic_sp2.R | `vcov.ic_po(object, constant = 1, ...)` | No (S3 method) |
| icsp2-package.R | None (package documentation only) | N/A |
| internal_utilities.R | `adjust_intervals(B = c(0, 1), surv_matrix, eps = 10^-10)` | Yes (ic_sp2) |
| internal_utilities.R | `findMaximalIntersections(lower, upper)` | Yes (.fit_ic_sp) |
| internal_utilities.R | `bs_sampleData(rawDataEnv, weights)` | No |
| internal_utilities.R | `getBS_coef(sampDataEnv, callText = 'ic_ph', other_info)` | No |
| internal_utilities.R | `expandX(formula, data, fit)` | Yes (get_etas, diag_covar) |
| internal_utilities.R | `removeSurvFromFormula(formula, ind = 2)` | Yes (expandY) |
| internal_utilities.R | `expandY(formula, data, fit)` | Yes (getResponse, imputeCens) |
| internal_utilities.R | `getResponse(fit, newdata = NULL)` | No |
| internal_utilities.R | `make_par_fitList(y_mat, x_mat, parFam = "gamma", link = "po", leftCen = 0, rightCen = Inf, uncenTol = 10^-6, regnames, weights, callText)` | No |
| internal_utilities.R | `imputeCensoredData_exp(l, u, impInfo, dist = 'web', maxVal)` | No |
| internal_utilities.R | `fullParamFit(formula, data, param_y, dist = 'weibull')` | No |
| internal_utilities.R | `fullParamFit_exp(formula, data, param_y, rightCenVal, dist = 'weibull')` | No |
| internal_utilities.R | `simPars_fromFit(fit, web = TRUE)` | No |
| internal_utilities.R | `simRawTimes(b1 = 0.5, b2 = -0.5, n = 100, shape1 = 2, shape2 = 2)` | No |
| internal_utilities.R | `simRawExpTimes(b1 = 0.5, b2 = -0.5, n = 100, rate = 1)` | No |
| internal_utilities.R | `subSampleData(data, max_n_use, weights)` | Yes (diag_baseline, diag_covar) |
| internal_utilities.R | `removeVarFromCall(call, varName)` | Yes (diag_covar) |
| internal_utilities.R | `splitData(data, varName, splits, splitFun, weights)` | Yes (splitAndFit) |
| internal_utilities.R | `makeFactorSplitInfo(vals, levels)` | Yes (diag_covar) |
| internal_utilities.R | `makeNumericSplitInfo(vals, cuts)` | Yes (diag_covar) |
| internal_utilities.R | `splitAndFit(newcall, data, varName, splitInfo, fitFunction, model, weights)` | Yes (diag_covar) |
| internal_utilities.R | `s_exp(x, par)` | Yes (get_s_fun) |
| internal_utilities.R | `s_weib(x, par)` | Yes (get_s_fun) |
| internal_utilities.R | `s_gamma(x, par)` | Yes (get_s_fun) |
| internal_utilities.R | `s_lnorm(x, par)` | Yes (get_s_fun) |
| internal_utilities.R | `s_loglgst(x, par)` | Yes (get_s_fun) |
| internal_utilities.R | `get_etas(fit, newdata = NULL, reg_pars = NULL)` | Yes (getSCurves, getFitEsts, sample_etas_and_base, lines.icenReg_fit, predict.icenReg_fit) |
| internal_utilities.R | `get_s_fun(fit)` | No |
| internal_utilities.R | `po_link(s, nu)` | Yes (get_link_fun) |
| internal_utilities.R | `ph_link(s, nu)` | Yes (get_link_fun) |
| internal_utilities.R | `no_link(s, nu)` | Yes (get_link_fun) |
| internal_utilities.R | `get_link_fun(fit)` | Yes (getSCurves) |
| internal_utilities.R | `findUpperBound(val = 1, x, s_fun, link_fun, fit, eta, baseline = NULL)` | No |
| internal_utilities.R | `subsetData_ifNeeded(i, data)` | No |
| internal_utilities.R | `getVarNames_fromFormula(formula)` | Yes (diag_covar) |
| internal_utilities.R | `getVar_fromRHS(rhs, names)` | Yes (getVarNames_fromFormula) |
| internal_utilities.R | `getFormula(object)` | Yes (diag_baseline, diag_covar) |
| internal_utilities.R | `getData(fit)` | Yes (diag_baseline, diag_covar) |
| internal_utilities.R | `get_tbull_mid_q(p, s_t, tbulls)` | No |
| internal_utilities.R | `get_tbull_mid_p(q, s_t, tbulls)` | No |
| internal_utilities.R | `PCAFit2OrgParFit(PCA_info, PCA_Hessian, PCA_parEsts, numIdPars)` | No |
| internal_utilities.R | `getNumCovars(object)` | No |
| internal_utilities.R | `getSurvProbs(times, etas, baselineInfo, regMod, baseMod)` | Yes (getFitEsts) |
| internal_utilities.R | `getSurvTimes(p, etas, baselineInfo, regMod, baseMod)` | Yes (getFitEsts) |
| internal_utilities.R | `getSamplablePars(fit)` | Yes (sampleSurv, ic_sample, imputeCens) |
| internal_utilities.R | `getSamplableVar(fit)` | Yes (sampleSurv, ic_sample) |
| internal_utilities.R | `sampBayesPar(fit)` | Yes (sampleSurv, ic_sample, imputeCens) |
| internal_utilities.R | `sampPars(mean, var)` | Yes (sampleSurv, ic_sample, imputeCens) |
| internal_utilities.R | `getBSParSample(fit)` | Yes (sampleSurv, ic_sample, imputeCens) |
| internal_utilities.R | `setSamplablePars(fit, coefs)` | Yes (sampleSurv, ic_sample, imputeCens) |
| internal_utilities.R | `fastNumericInsert(newVals, target, indices)` | Yes (fastMatrixInsert) |
| internal_utilities.R | `fastMatrixInsert(newVals, targMat, rowNum = NULL, colNum = NULL)` | Yes (imputeCens) |
| internal_utilities.R | `updateDistPars(vals, max_n)` | Yes (dGeneralGamma, qGeneralGamma, pGeneralGamma) |
| internal_utilities.R | `getMaxLength(thisList)` | Yes (dGeneralGamma, qGeneralGamma, pGeneralGamma) |
| internal_utilities.R | `addIfMissing(val, name, list)` | Yes (addListIfMissing, plot.icenReg_fit, plot.sp_curves, plot.ic_npList) |
| internal_utilities.R | `addListIfMissing(listFrom, listInto)` | Yes (plot.icenReg_fit) |
| internal_utilities.R | `readingCall(mf)` | Yes (ic_sp2) |
| internal_utilities.R | `makeIntervals(y, mf)` | Yes (make_xy) |
| internal_utilities.R | `check_weights(model.frame)` | Yes (ic_sp2) |
| internal_utilities.R | `checkStrata(strata, yMat)` | No |
| internal_utilities.R | `check_matrix(x)` | Yes (ic_sp2) |
| internal_utilities.R | `icr_nrow(x)` | Yes (plot.icenReg_fit, lines.icenReg_fit) |
| internal_utilities.R | `icr_ncol(x)` | No |
| internal_utilities.R | `icr_colMeans(x)` | No |
| internal_utilities.R | `get_dataframe_row(df, row)` | Yes (survCIs reference class) |
| internal_utilities.R | `default_baseline(fit)` | No |
| internal_utilities.R | `default_reg_pars(fit)` | No |
| internal_utilities.R | `sample_in_interval(fit, newdata, lower_time, upper_time)` | No |
| internal_utilities.R | `sample_pars(fit, samples = 100)` | No |
| internal_utilities.R | `sample_etas_and_base(fit, samples, newdata)` | Yes (sampleSurv) |
| internal_utilities.R | `subtractOffset(new_x, offset)` | Yes (get_etas) |
| internal_utilities.R | `checkFor_cluster(form)` | No |
| internal_utilities.R | `make_xy(frml, df)` | No |
| internal_utilities.R | `icColMeans(x)` | No |
| plottingFunctions.R | `plot.icenReg_fit(x, y, newdata = NULL, fun = 'surv', plot_legend = T, cis = T, ci_level = 0.9, survRange = c(0.025, 1), evalPoints = 200, lgdLocation = lgd_default(fun), xlab = "time", ...)` | No (S3 method) |
| plottingFunctions.R | `lines.surv_cis(x, col, include_cis = T, fun = "surv", ...)` | No (S3 method) |
| plottingFunctions.R | `lgd_default(fun_type = "surv")` | Yes (plot.icenReg_fit, plot.ic_npList, diag_covar) |
| plottingFunctions.R | `lines.icenReg_fit(x, y, newdata = NULL, fun = 'surv', cis = F, ci_level = 0.9, survRange = c(0.025, 1), evalPoints = 20, ...)` | Yes (plot.icenReg_fit) |
| plottingFunctions.R | `plot.sp_curves(x, sname = 'baseline', xRange = NULL, ...)` | No (S3 method) |
| plottingFunctions.R | `lines.ic_npList(x, fitNames = NULL, fun = "surv", ...)` | Yes (plot.ic_npList) |
| plottingFunctions.R | `plot.ic_npList(x, fitNames = NULL, fun = "surv", lgdLocation = lgd_default(fun), plot_legend = T, ... )` | No (S3 method) |
| plottingFunctions.R | `lines.sp_curves(x, sname = 'baseline', fun = "surv", ...)` | Yes (plot.sp_curves) |
| referenceClasses.R | `icenReg_fit` (reference class) | N/A (class definition) |
| referenceClasses.R | `sp_fit_class` (reference class) | N/A (class definition) |
| referenceClasses.R | `ic_np_class` (reference class) | N/A (class definition) |
| referenceClasses.R | `ic_ph_class` (reference class) | N/A (class definition) |
| referenceClasses.R | `ic_po_class` (reference class) | N/A (class definition) |
| referenceClasses.R | `par_class` (reference class) | N/A (class definition) |
| referenceClasses.R | `bayes_fit` (reference class) | N/A (class definition) |
| referenceClasses.R | `icenRegSummary` (reference class) | N/A (class definition) |
| referenceClasses.R | `ic_npList` (reference class) | N/A (class definition) |
| referenceClasses.R | `surv_cis` (reference class) | Yes (survCIs) |
| simulationFunctions.R | `simEventTime(linPred = 0, model = 'ph', dist = qweibull, paramList = list(shape = 1, scale = 1))` | Yes (simIC_weib, simCS_weib, simDC_weib) |
| simulationFunctions.R | `simIC_weib(n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, scale = 2, inspections = 2, inspectLength = 2.5, rndDigits = NULL, prob_cen = 1)` | No (user-facing) |
| simulationFunctions.R | `simCS_weib(n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, scale = 2)` | No (user-facing) |
| simulationFunctions.R | `simDC_weib(n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2, scale = 2, lowerLimit = 0.75, upperLimit = 2)` | No (user-facing) |
| user_utilities.R | `names.icenReg_fit(x)` | No (S3 method) |
| user_utilities.R | `getSCurves(fit, newdata = NULL)` | Yes (plot.icenReg_fit, lines.icenReg_fit, diag_covar, diag_baseline, getFitEsts) |
| user_utilities.R | `getSCurves.default(fit, newdata = NULL)` | No (S3 method) |
| user_utilities.R | `range.sp_curves_list(..., na.rm = FALSE, finite = FALSE)` | No (S3 method) |
| user_utilities.R | `getSCurves.ic_np(fit, newdata = NULL)` | No (S3 method) |
| user_utilities.R | `summary.icenReg_fit(object, ...)` | No (S3 method) |
| user_utilities.R | `summary.ic_npList(object, ...)` | No (S3 method) |
| user_utilities.R | `simIC_weib(n = 100, b1 = 0.5, b2 = -0.5, model = 'ph', shape = 2, scale = 2, inspections = 2, inspectLength = 2.5, rndDigits = NULL, prob_cen = 1)` | No (duplicate/user-facing) |
| user_utilities.R | `simICPO_beta(n = 100, b1 = 1, b2 = -1, inspections = 1, shape1 = 2, shape2 = 2, rndDigits = NULL)` | No |
| user_utilities.R | `diag_covar(object, varName, data, model, weights = NULL, yType = 'meanRemovedTransform', factorSplit = TRUE, numericCuts, col, xlab, ylab, main, lgdLocation = NULL)` | No (user-facing) |
| user_utilities.R | `getFitEsts(fit, newdata = NULL, p, q)` | Yes (plot.icenReg_fit, lines.icenReg_fit, diag_covar, predict.icenReg_fit, imputeCens, ic_sample, icqqplot) |
| user_utilities.R | `diag_baseline(object, data, model = 'ph', weights = NULL, dists = c('exponential', 'weibull', 'gamma', 'lnorm', 'loglogistic', 'generalgamma'), cols = NULL, lgdLocation = 'bottomleft', useMidCovars = T)` | No (user-facing) |
| user_utilities.R | `predict.icenReg_fit(object, type = 'response', newdata = NULL, ...)` | No (S3 method) |
| user_utilities.R | `imputeCens(fit, newdata = NULL, imputeType = 'fullSample', samples = 5)` | No (user-facing) |
| user_utilities.R | `ic_sample(fit, newdata = NULL, sampleType = 'fullSample', samples = 5)` | No (user-facing) |
| user_utilities.R | `sampleSurv_slow(fit, newdata, p = NULL, q = NULL, samples = 100)` | No |
| user_utilities.R | `sampleSurv(fit, newdata = NULL, p = NULL, q = NULL, samples = 100)` | No (user-facing) |
| user_utilities.R | `dGeneralGamma(x, mu, s, Q)` | No |
| user_utilities.R | `qGeneralGamma(p, mu, s, Q)` | No |
| user_utilities.R | `pGeneralGamma(q, mu, s, Q)` | No |
| user_utilities.R | `icqqplot(par_fit)` | No |
| user_utilities.R | `cs2ic(time, eventOccurred)` | No (user-facing) |
| user_utilities.R | `survCIs(fit, newdata = NULL, p = NULL, q = NULL, ci_level = 0.95, MC_samps = 4000)` | Yes (plot.icenReg_fit, lines.icenReg_fit) |
