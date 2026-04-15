# print works for ic_sp2 / PH model works for sim data

    Code
      print(result)
    Output
      Call:
      ic_sp_ph(formula = Surv(l, u, type = "interval2") ~ x1 + x2, 
          data = sim_data, control = ic_sp_control(derivMethod = 1))
      
      Coefficients:
              x1         x2 
       0.8030607 -0.8795480 
      
      Log-likelihood: -80.9277
      Number of iterations: 17

# print works for stratified fit

    Code
      print(result)
    Output
      Call:
      ic_sp_ph(formula = Surv(l, u, type = "interval2") ~ grp + strata(strata), 
          data = md)
      
      Coefficients:
          grpge 
      0.6133488 
      
      Log-likelihood: -74.768
      Number of iterations: 23

