 // Random walk model for survey averaging using log-scale process and observation errors
DATA_SECTION
  init_int styr
  init_int endyr
  ivector yrs(styr,endyr);
  !! yrs.fill_seqadd(styr,1);
  init_int nobs
  init_ivector yrs_srv(1,nobs)
  init_vector srv_est(1,nobs)
  init_vector srv_cv(1,nobs)
  vector srv_sd(1,nobs)
  number meany
  vector yvar(1,nobs)
  vector yconst(1,nobs)
 
 LOC_CALCS
    logdat( styr);
    logdat( endyr);
    logdat(nobs);
    logdat( yrs_srv);
    logdat( srv_est);
    logdat( srv_cv);
    if (mean(srv_cv)>5) srv_cv = elem_div(srv_cv,srv_est+0.0001);
    srv_sd = elem_prod(srv_cv,srv_cv) + 1.;
    srv_sd = sqrt(log(srv_sd));
    yvar = elem_prod(srv_sd,srv_sd);
    yconst = log(2.0*M_PI*yvar);
 END_CALCS
 
PARAMETER_SECTION
  init_number logSdLam(1)
  sdreport_vector biomsd(styr,endyr);
  sdreport_vector biomA(styr,endyr);
  random_effects_vector biom(styr,endyr);
  // init_vector biom(styr,endyr);
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0.0;
  for(int i=styr+1; i<=endyr; ++i)
  {
    step(biom(i-1),biom(i),logSdLam);
  }

  for(int i=1; i<=nobs; ++i)
  {
    obs(biom(yrs_srv(i)),i);
  }
  if (sd_phase()) 
  {
    biomA = exp(biom);
    biomsd = biom;
  }
	if (mceval_phase())
  {
    biomA = exp(biom);
    biomsd = biom;
	  cout <<biomA<<endl;
	}
	jnll += square(logSdLam+1.5);// Modest penalty to keep process error from getting too large...

SEPARABLE_FUNCTION void step(const dvariable& biom1, const dvariable& biom2, const dvariable& logSdLam)
  dvariable var=exp(2.0*logSdLam);
  jnll+=0.5*(log(2.0*M_PI*var)+square(biom2-biom1)/var);

SEPARABLE_FUNCTION void obs(const dvariable& biom, int i)
  jnll+=0.5*(yconst(i) + square(biom-log(srv_est(i)))/yvar(i));

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(3000);

REPORT_SECTION
  biomsd = biom;
  report << biomsd <<endl;
  // biomsd = meany*biom;
  // report << "yr"  <<endl;
  // report <<  yr   <<endl;
  // report << "est" <<endl;
GLOBALS_SECTION
  #include <admodel.h>
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("rwout.rep");
  ofstream checkin("log_read.rep");
  #undef logdat
  #define logdat(object) checkin << #object "\n" << object << endl;
  adstring sppname;

FINAL_SECTION
  /*
  dvar_vector UCI = elem_prod(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector LCI = elem_div(biomsd,exp(1.96*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector upp90th =  elem_prod(biomsd,exp(1.645*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  dvar_vector low90th =  elem_div(biomsd,exp(1.645*sqrt(log(1.+elem_div(elem_prod(biomsd.sd,biomsd.sd),
                         elem_prod(biomsd,biomsd))))));
  */
  dvar_vector UCI = exp(biomsd+1.96*biomsd.sd);
  dvar_vector LCI = exp(biomsd-1.96*biomsd.sd);
  dvar_vector upp90th = exp(biomsd+1.645*biomsd.sd);
  dvar_vector low90th = exp(biomsd-1.645*biomsd.sd);

  write_R(yrs_srv);
  write_R(srv_est);
  write_R(srv_sd);
  write_R(yrs);
  write_R(LCI);
  write_R(biomA);
  write_R(UCI);
  write_R(low90th);
  write_R(upp90th);
  write_R(biomsd);
  write_R(biomsd.sd);

  mysum.close();
