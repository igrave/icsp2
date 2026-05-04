# cpp_logrank matches survdiff for 2 groups

    Code
      print(res)
    Output
      
      Log-rank test
      Sample sizes by group:
        1   2 
      100 100 
      
                N Observed Expected     O-E
      Group 1 100       73  74.9074 -1.9074
      Group 2 100       72  70.0926  1.9074
      
      Variance-covariance matrix:
               1        2
      1  35.7382 -35.7382
      2 -35.7382  35.7382

# cpp_logrank matches survdiff for 3 groups

    Code
      print(res)
    Output
      
      Log-rank test
      Sample sizes by group:
        1   2   3 
      100 100 100 
      
                N Observed Expected      O-E
      Group 1 100       78  94.8389 -16.8389
      Group 2 100       70  81.6248 -11.6248
      Group 3 100       81  52.5364  28.4636
      
      Variance-covariance matrix:
               1        2        3
      1  53.9661 -33.3797 -20.5864
      2 -33.3797  51.9912 -18.6115
      3 -20.5864 -18.6115  39.1980

