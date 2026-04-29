# ic test with sas type works with strata

    Code
      print(result)
    Output
      
      Stratified Sun's log-rank test for interval-censored data 
      ========================================================= 
      
      Call:
      ic_logrank(formula = Surv(l, u, type = "interval2") ~ arm + strata(group), 
          data = stratified_data, n_samples = 1000, type = "sas")
      
      Sample sizes by group:
            Strata
      Group   A  B
        trtA 52 50
        trtB 48 50
      
      Log-rank Statistics:
           trtA      trtB 
       24.03557 -24.03557 
      
        Chi-squared statistic: Q = 28.24
        Degrees of freedom:    df = 1
        P-value:               p = 1.073e-07
      
      Variance-covariance matrix:
               trtA     trtB
      trtA  20.4593 -20.4593
      trtB -20.4593  20.4593
      
      Calculated with  1000  samples

# ic test with hly type works with strata

    Code
      print(result)
    Output
      
      Stratified Sun's log-rank test for interval-censored data 
      ========================================================= 
      
      Call:
      ic_logrank(formula = Surv(l, u, type = "interval2") ~ arm + strata(group), 
          data = stratified_data, n_samples = 1000, type = "hly")
      
      Sample sizes by group:
            Strata
      Group   A  B
        trtA 52 50
        trtB 48 50
      
      Log-rank Statistics:
          trtA     trtB 
       24.0572 -24.0572 
      
        Chi-squared statistic: Q = 28.29
        Degrees of freedom:    df = 1
        P-value:               p = 1.046e-07
      
      Variance-covariance matrix:
               trtA     trtB
      trtA  20.4593 -20.4593
      trtB -20.4593  20.4593
      
      Calculated with  1000  samples

