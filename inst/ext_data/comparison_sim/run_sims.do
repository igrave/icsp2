/* Stata do-file for running stintcox on all replications */
cd "/home/isaac/stata/share/icsp2/comparison_sim"
program define run_stintcox
  args scenario n_reps has_strata has_age

  local n_coefs = 1 + ("`has_age'" == "1")

  forvalues i = 1/`n_reps' {
    local fname "sim_data/`scenario'/rep_`=string(`i', "%03.0f")'.csv"
    import delimited "`fname'", clear

    * Handle right-censored (Inf) values
    replace u = "." if u == "Inf"
    destring u, replace

    if "`has_strata'" == "1" {
      capture confirm numeric variable strata_num
      if _rc != 0 {
        egen strata_num = group(strata)
      }
    }
    
    stintcox trt `=cond("`has_age'"=="1", "age", "")' ///
      , interval(l u) `=cond("`has_strata'"=="1", "strata(strata_num)", "")' nohr emtolerance(1e-8)  emiterate(5000) vce(opg)

    * Store coefficients and SEs
    matrix b = e(b)
    matrix V = e(V)
    matrix c = e(converged)
    matrix it = e(ic)
    matrix se = J(1, `n_coefs', .)
    forvalues c = 1/`n_coefs' {
      matrix se[1, `c'] = sqrt(V[`c', `c'])
    }

    if `i' == 1 {
      matrix all_b = b
      matrix all_se = se
      matrix all_c = c
      matrix all_it = it
    }
    else {
      matrix all_b = all_b \ b
      matrix all_se = all_se \ se
      matrix all_c = all_c \ c
      matrix all_c = all_it \ it
    }

    * Save baseline survival from first replication
    *if `i' == 1 {
    *  predict double surv0, basesurv
    *  predict double haz0, basehaz
    *  gen double time_grid = _t
    *  keep time_grid surv0 haz0
    *  duplicates drop
    *  sort time_grid
    *  export delimited "sim_data/`scenario'/stata_baseline.csv", replace
    *}
  }

  * Export coefficient and SE matrices to CSV
  clear
  svmat all_b, names(b)
  svmat all_se, names(se)
  svmat all_c, names(c)
  svmat all_it, names(it)
  gen rep = _n
  export delimited "output_stata/`scenario'/stata_coefs.csv", replace
end

*run_stintcox S1_null 100 0 0
timer on 1

run_stintcox S2_moderate 100 0 0
timer off 1
display "--- Process Finished at: " c(current_time) " ---"
timer list 1
*
*run_stintcox S3_strong 100 0 0
*run_stintcox S4_weak 100 0 0
*run_stintcox S5_moderate_age 100 0 1
*run_stintcox S6_strat_homo 100 1 0
*run_stintcox S7_strat_diffbase 100 1 0
*run_stintcox S8_strat_hetero 20 1 0
