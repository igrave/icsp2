%let root = /home/gravesti/icsp2_compare;

/* SAS macro for running PROC ICLIFETEST on all replications */
%macro run_iclifetest(scenario, n_reps=100, has_strata=0);
  %do i = 1 %to &n_reps;
    proc import datafile="&root/sim_data/&scenario/rep_%sysfunc(putn(&i, z3.)).csv"
     
 	out=work.rep_&i dbms=csv replace;
	GUESSINGROWS=200;
	run;
   
	/* Convert u from character to numeric (Inf -> .) */
    data work.rep_&i;
      set work.rep_&i(rename=(u=u_char));
      if upcase(u_char) = "INF" then u = .;
      else u = input(u_char, best.);
      drop u_char;
    run;

  ods output HomTests = _ictests_i (keep = Test ChiSq ProbChiSq
                                      rename = (ChiSq = ChiSquare));
  
  proc iclifetest data=work.rep_&i impute(seed=1000 NIMTEST=10000)
	%if &i = 1  %then %do;
        plots=(survival) outsurv=work.surv_1
    %end;
    conftype=linear;
   time (l, u);
   %if &has_strata = 1 %then %do;
        strata strata;
   %end;
   test trt / weight=sun;
   run;
    ods output close;

/*
    proc iclifetest data=work.rep_&i
      %if &i = 1 %then %do;
        plots=(survival) outsurv=work.surv_1
      %end;
      ;
      time (l, u);
      %if &has_strata = 1 %then %do;
        strata strata;
      %end;
      test trt;
      ods output LogRankTest=work.lr_&i;
    run;
*/
    /* Add replication index */
   data work.lr_&i;
      set _ictests_i;
      rep = &i;
    run;
  %end;

  /* Stack all log-rank results and export */
  data work.&scenario._lr;
    set work.lr_1 - work.lr_&n_reps;
  run;

  proc export data=work.&scenario._lr
    outfile="&root/output/&scenario/sas_logrank.csv"
    dbms=csv replace;
  run;

  /* Export NPMLE survival from first replication */
  proc export data=work.surv_1
    outfile="&root/output/&scenario/sas_npmle_surv.csv"
    dbms=csv replace;
  run;
%mend;



%run_iclifetest(S1_null, has_strata=0);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S2_moderate, has_strata=0);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S3_strong, has_strata=0);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S4_weak, has_strata=0);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S5_moderate_age, has_strata=0);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S6_strat_homo, has_strata=1);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S7_strat_diffbase, has_strata=1);
proc datasets lib=work kill nolist;
run;
quit;
%run_iclifetest(S8_strat_hetero, has_strata=1);
proc datasets lib=work kill nolist;
run;
quit;