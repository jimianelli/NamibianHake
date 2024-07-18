#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include <admodel.h>
  ofstream csvrep("nh_out.csv");
  #undef csv_rep
  #define csv_rep(object) csvrep << #object "\n" << object << endl;
  ofstream Rreport("nh_R.rep");
  #undef R_report
  #define R_report(object) Rreport << #object "\n" << object << endl;
  ofstream log_input("log_input.rep");
  #undef log_input
  #define log_input(object) log_input << "#"<<#object "\n" << object << endl;
  adstring model_name;
  adstring datafilename;
  adstring model_number;
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <nh.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
    /*
  if (ad_comm::argc > 1) {
    int on=0;
		if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-steepness"))>-1) {
      if (on>ad_comm::argc-2 | ad_comm::argv[on+1][0] == '-') {
        cerr << "Invalid number of arguments for command line option -steepness; option ignored" << endl;
      }
      else {
        h_in = std::stod(ad_comm::argv[on+1]);
				// atoi(ad_comm::argv[on+1]);
        cout<<  "Currently using "<<adstring(ad_comm::argv[on+1])<<" as steepness value"<<endl;
      }
    }
  }
  */
 *(ad_comm::global_datafile) >>  model_number; 
 *(ad_comm::global_datafile) >>  model_name;      //e.g. q08
 *(ad_comm::global_datafile) >>  datafilename;    //e.g. namhakdata.dat 
  NProj.allocate("NProj");
  first_yr.allocate("first_yr");
  last_yr.allocate("last_yr");
  plus_grp.allocate("plus_grp");
  Minfage.allocate("Minfage");
  est_M.allocate("est_M");
  est_Minf.allocate("est_Minf");
  use_multinomial.allocate("use_multinomial");
  est_add_seal_sigma.allocate("est_add_seal_sigma");
  est_prop.allocate("est_prop");
  est_Steep.allocate("est_Steep");
  steep_in.allocate("steep_in");
  steep_cv.allocate("steep_cv");
  Ksp_fraction.allocate("Ksp_fraction");
  NSelPeriods.allocate("NSelPeriods");
  SelYrs.allocate(1,NSelPeriods,1,2,"SelYrs");
  est_Sel.allocate("est_Sel");
  est_cSel.allocate("est_cSel");
  est_Sel_slope.allocate("est_Sel_slope");
  est_Sel_slopeS.allocate("est_Sel_slopeS");
  Age_Sel_slopeS.allocate("Age_Sel_slopeS");
  Age_Sel_slope.allocate("Age_Sel_slope");
  Rec_est.allocate("Rec_est");
  SigmaRec.allocate("SigmaRec");
  RecYr1.allocate("RecYr1");
  RecYr2.allocate("RecYr2");
  CAAMinus.allocate("CAAMinus");
  CAAPlus.allocate("CAAPlus");
  NCpueSeries.allocate("NCpueSeries");
  UseCPUE.allocate(1,NCpueSeries,"UseCPUE");
  CPUE_Sigma.allocate(1,NCpueSeries,"CPUE_Sigma");
  CpueIndx.allocate(1,NCpueSeries,"CpueIndx");
  CpueWght.allocate(1,NCpueSeries,"CpueWght");
  NSurveySeries.allocate("NSurveySeries");
  UseSurvey.allocate(1,NSurveySeries,"UseSurvey");
  SurveyIndx.allocate(1,NSurveySeries,"SurveyIndx");
  q_input.allocate(1,NSurveySeries,"q_input");
  q_chngYr.allocate(1,NSurveySeries,"q_chngYr");
  q_CV.allocate(1,NSurveySeries,"q_CV");
  q_mean.allocate(1,NSurveySeries,"q_mean");
  est_addvar.allocate("est_addvar");
  CAASMinus.allocate(1,NSurveySeries,"CAASMinus");
  CAASPlus.allocate(1,NSurveySeries,"CAASPlus");
  FutSurCV.allocate("FutSurCV");
  FutBiasCV.allocate("FutBiasCV");
  ExtraCV.allocate("ExtraCV");
  FutSurBias.allocate("FutSurBias");
  NKspPeriods.allocate("NKspPeriods");
  KspYrs.allocate(1,NKspPeriods,1,2,"KspYrs");
  if (KspYrs(1,1) != first_yr) cout << "Warning: first Ksp year must be " << first_yr << endl;  
  if (KspYrs(NKspPeriods,2) != last_yr) cout << "Warning: last Ksp year must be " << last_yr << endl;  
  weight.allocate("weight");
  select.allocate("select");
  propF.allocate("propF");
  propW.allocate("propW");
  TrwdayF.allocate("TrwdayF");
  TrwdayW.allocate("TrwdayW");
  fishTimF.allocate("fishTimF");
  fishTimW.allocate("fishTimW");
  EmpF.allocate("EmpF");
  EmpW.allocate("EmpW");
  EmpFFF.allocate("EmpFFF");
  EmpFFW.allocate("EmpFFW");
  ProfF.allocate("ProfF");
  ProfW.allocate("ProfW");
  dr.allocate("dr");
  Pri.allocate("Pri");
  AocF.allocate("AocF");
  AocW.allocate("AocW");
  AocFFF.allocate("AocFFF");
  AocFFW.allocate("AocFFW");
  ApF.allocate("ApF");
  ApW.allocate("ApW");
  eof.allocate("eof");
  log_input(model_number);
  log_input(model_name);
  log_input(datafilename);
  log_input(NProj);
  log_input(first_yr);
  log_input(last_yr);
  log_input(plus_grp);
  log_input(Minfage);
  log_input(est_M);
  log_input(est_Minf);
  log_input(use_multinomial);
  log_input(est_add_seal_sigma);
  log_input(est_prop);
  log_input(est_Steep);
  log_input(steep_in);
  log_input(steep_cv);
  log_input(Ksp_fraction);
  log_input(NSelPeriods);
  log_input(SelYrs);
  log_input(est_Sel);
  log_input(est_cSel);
  log_input(est_Sel_slope);
  log_input(est_Sel_slopeS);
  log_input(Age_Sel_slopeS);
  log_input(Age_Sel_slope);
  log_input(Rec_est);
  log_input(SigmaRec);
  log_input(RecYr1);
  log_input(RecYr2);
  
  log_input(CAAMinus);
  log_input(CAAPlus);
  log_input(NCpueSeries);
  log_input( UseCPUE);
  log_input( CPUE_Sigma);
  log_input( CpueIndx);
  log_input(CpueWght);
   
  log_input(NSurveySeries);
  log_input( UseSurvey);
  log_input( SurveyIndx);
  log_input(q_input);
  log_input( q_chngYr);
  log_input(q_CV);
  log_input(q_mean);
  log_input(est_addvar);
  log_input( CAASMinus);
  log_input( CAASPlus);
  
  // ** For projections **
  // ---------------------
  log_input(FutSurCV);
  log_input(FutBiasCV);
  log_input(ExtraCV);
  log_input(FutSurBias);
  
  log_input(NKspPeriods);
  log_input(KspYrs);
  log_input(weight);
  log_input(select);
  log_input(propF);
  log_input(propW);
  log_input(TrwdayF);
  log_input(TrwdayW);
  log_input(fishTimF);
  log_input(fishTimW);
  log_input(EmpF);
  log_input(EmpW);
  log_input(EmpFFF);
  log_input(EmpFFW);
  log_input(ProfF);
  log_input(ProfW);
  log_input(dr);
  log_input(Pri);
  log_input(AocF);
  log_input(AocW);
  log_input(AocFFF);
  log_input(AocFFW);
  log_input(ApF);
  log_input(ApW);
  log_input(eof);
 if (eof!=54321) {cout<<"Oh shit..."<<endl;exit(1);}
 ad_comm::change_datafile_name(datafilename);
  Wstrt.allocate(first_yr,last_yr,-1,plus_grp,"Wstrt");
  Wmid.allocate(first_yr,last_yr,-1,plus_grp,"Wmid");
  Catch.allocate(first_yr,last_yr,"Catch");
  CAA.allocate(first_yr,last_yr,-2,plus_grp,"CAA");
  F_N.allocate(first_yr,last_yr);
  CPUE.allocate(first_yr,last_yr,0,NCpueSeries,"CPUE");
  Survey.allocate(first_yr,last_yr,0,NSurveySeries*2,"Survey");
  SurvCAA.allocate(1,NSurveySeries,first_yr,last_yr,-2,plus_grp,"SurvCAA");
  S_N.allocate(1,NSurveySeries,first_yr,last_yr);
  seali.allocate(first_yr,last_yr,"seali");
  sealCV.allocate(first_yr,last_yr,"sealCV");
  mat.allocate(0,plus_grp,"mat");
  SSB.allocate(1,20);
  eof1.allocate("eof1");
 if (eof1!=12345) {cout<<"Oh crap..."<<endl;exit(1);}
Nsim = -1;
  rndCPUE.allocate(last_yr,last_yr+NProj);
  rndSurv.allocate(last_yr+1,last_yr+NProj+1);
  rndRec.allocate(last_yr+1,last_yr+NProj+1);
}

void model_parameters::initializationfunction(void)
{
  Minf.set_initial_value(0.35606928);
  sig_seal_added.set_initial_value(0.001);
  Addvar.set_initial_value(0.00001);
  par_B0.set_initial_value(8.);
  Prop.set_initial_value(1.);
  Steep.set_initial_value(steep_in);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  int NselPar;
  NselPar = 2*NSelPeriods;
  par_B0.allocate(-1.0,13.0,"par_B0");
  M1.allocate(0.1,2.0,est_M,"M1");
  Minf.allocate(0.1,0.5,est_Minf,"Minf");
  Sage.allocate(1,NselPar,0.5,10.0,est_cSel,"Sage");
  SelSlope.allocate(1,NSelPeriods+1,0.0,1.0,est_Sel_slope,"SelSlope");
  SurvA50.allocate(0.0,plus_grp,est_Sel,"SurvA50");
  SurvAd.allocate(0.001,100.0,est_Sel,"SurvAd");
  Addvar.allocate(0.0,1.0,est_addvar,"Addvar");
  sig_seal_added.allocate(0.0,1.0,est_add_seal_sigma,"sig_seal_added");
  RecPar.allocate(RecYr1,RecYr2,Rec_est,"RecPar");
  Prop.allocate(0.2,2.0,est_prop,"Prop");
  Steep.allocate(0.21,0.99,est_Steep,"Steep");
  SelSlopeS.allocate(0.0,1.0,est_Sel_slopeS,"SelSlopeS");
  N.allocate(first_yr,last_yr+NProj+1,0,plus_grp,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  NM.allocate(first_yr,last_yr+NProj+1,0,plus_grp,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  S.allocate(first_yr,last_yr+NProj,0,plus_grp,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  Sfut.allocate(0,plus_grp,"Sfut");
  #ifndef NO_AD_INITIALIZE
    Sfut.initialize();
  #endif
  SurvS.allocate(0,plus_grp,"SurvS");
  #ifndef NO_AD_INITIALIZE
    SurvS.initialize();
  #endif
  F.allocate(first_yr,last_yr+NProj,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  M.allocate(0,plus_grp,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  Nvirg.allocate(0,plus_grp,"Nvirg");
  #ifndef NO_AD_INITIALIZE
    Nvirg.initialize();
  #endif
  Ntemp.allocate(0,plus_grp,"Ntemp");
  #ifndef NO_AD_INITIALIZE
    Ntemp.initialize();
  #endif
  Mtemp.allocate(0,plus_grp,"Mtemp");
  #ifndef NO_AD_INITIALIZE
    Mtemp.initialize();
  #endif
  Ntemp2.allocate(0,plus_grp,"Ntemp2");
  #ifndef NO_AD_INITIALIZE
    Ntemp2.initialize();
  #endif
  SAPred.allocate(1,NSurveySeries,first_yr,last_yr+NProj,0,plus_grp,"SAPred");
  #ifndef NO_AD_INITIALIZE
    SAPred.initialize();
  #endif
  CAPred.allocate(first_yr,last_yr+NProj,0,plus_grp,"CAPred");
  #ifndef NO_AD_INITIALIZE
    CAPred.initialize();
  #endif
  CAPred2.allocate(first_yr,last_yr+NProj,0,plus_grp,"CAPred2");
  #ifndef NO_AD_INITIALIZE
    CAPred2.initialize();
  #endif
  Recruit_Res.allocate(first_yr,last_yr+NProj+1,"Recruit_Res");
  #ifndef NO_AD_INITIALIZE
    Recruit_Res.initialize();
  #endif
  SigR_out.allocate("SigR_out");
  #ifndef NO_AD_INITIALIZE
  SigR_out.initialize();
  #endif
  Alpha.allocate(first_yr,last_yr,"Alpha");
  #ifndef NO_AD_INITIALIZE
    Alpha.initialize();
  #endif
  Beta.allocate(first_yr,last_yr,"Beta");
  #ifndef NO_AD_INITIALIZE
    Beta.initialize();
  #endif
  Spawn.allocate(first_yr,last_yr+NProj+1,"Spawn");
  #ifndef NO_AD_INITIALIZE
    Spawn.initialize();
  #endif
  Bexp.allocate(first_yr,last_yr+NProj,"Bexp");
  #ifndef NO_AD_INITIALIZE
    Bexp.initialize();
  #endif
  BBeg.allocate(first_yr,last_yr+NProj,"BBeg");
  #ifndef NO_AD_INITIALIZE
    BBeg.initialize();
  #endif
  SurvMid.allocate(first_yr,last_yr+NProj,"SurvMid");
  #ifndef NO_AD_INITIALIZE
    SurvMid.initialize();
  #endif
  SurvBeg.allocate(first_yr,last_yr+NProj,"SurvBeg");
  #ifndef NO_AD_INITIALIZE
    SurvBeg.initialize();
  #endif
  TotalB.allocate(first_yr,last_yr+NProj+1,"TotalB");
  #ifndef NO_AD_INITIALIZE
    TotalB.initialize();
  #endif
  BMid.allocate(first_yr,last_yr+NProj+1,"BMid");
  #ifndef NO_AD_INITIALIZE
    BMid.initialize();
  #endif
  seali_pred.allocate(first_yr,last_yr+NProj,"seali_pred");
  #ifndef NO_AD_INITIALIZE
    seali_pred.initialize();
  #endif
  BIO.allocate(first_yr,last_yr+NProj,"BIO");
  #ifndef NO_AD_INITIALIZE
    BIO.initialize();
  #endif
  Mbio.allocate(first_yr,last_yr+NProj+1,0,plus_grp,"Mbio");
  #ifndef NO_AD_INITIALIZE
    Mbio.initialize();
  #endif
  Recruits.allocate(first_yr,last_yr+NProj+1,"Recruits");
  #ifndef NO_AD_INITIALIZE
    Recruits.initialize();
  #endif
  futCatch.allocate(last_yr-1,last_yr+NProj+1,"futCatch");
  #ifndef NO_AD_INITIALIZE
    futCatch.initialize();
  #endif
  TAC.allocate(first_yr,last_yr+NProj+1,"TAC");
  #ifndef NO_AD_INITIALIZE
    TAC.initialize();
  #endif
  CPUEfut.allocate(last_yr,last_yr+NProj+1,"CPUEfut");
  #ifndef NO_AD_INITIALIZE
    CPUEfut.initialize();
  #endif
  TACF.allocate(last_yr,last_yr+NProj+1,"TACF");
  #ifndef NO_AD_INITIALIZE
    TACF.initialize();
  #endif
  TACW.allocate(last_yr,last_yr+NProj+1,"TACW");
  #ifndef NO_AD_INITIALIZE
    TACW.initialize();
  #endif
  VesselF.allocate(last_yr,last_yr+NProj+1,"VesselF");
  #ifndef NO_AD_INITIALIZE
    VesselF.initialize();
  #endif
  VesselW.allocate(last_yr,last_yr+NProj+1,"VesselW");
  #ifndef NO_AD_INITIALIZE
    VesselW.initialize();
  #endif
  CostF.allocate(last_yr,last_yr+NProj+1,"CostF");
  #ifndef NO_AD_INITIALIZE
    CostF.initialize();
  #endif
  CostW.allocate(last_yr,last_yr+NProj+1,"CostW");
  #ifndef NO_AD_INITIALIZE
    CostW.initialize();
  #endif
  CostT.allocate(last_yr,last_yr+NProj+1,"CostT");
  #ifndef NO_AD_INITIALIZE
    CostT.initialize();
  #endif
  RevF.allocate(last_yr,last_yr+NProj+1,"RevF");
  #ifndef NO_AD_INITIALIZE
    RevF.initialize();
  #endif
  RevW.allocate(last_yr,last_yr+NProj+1,"RevW");
  #ifndef NO_AD_INITIALIZE
    RevW.initialize();
  #endif
  RevT.allocate(last_yr,last_yr+NProj+1,"RevT");
  #ifndef NO_AD_INITIALIZE
    RevT.initialize();
  #endif
  ProfitF.allocate(last_yr,last_yr+NProj+1,"ProfitF");
  #ifndef NO_AD_INITIALIZE
    ProfitF.initialize();
  #endif
  ProfitW.allocate(last_yr,last_yr+NProj+1,"ProfitW");
  #ifndef NO_AD_INITIALIZE
    ProfitW.initialize();
  #endif
  ProfitT.allocate(last_yr,last_yr+NProj+1,"ProfitT");
  #ifndef NO_AD_INITIALIZE
    ProfitT.initialize();
  #endif
  EmployF.allocate(last_yr,last_yr+NProj+1,"EmployF");
  #ifndef NO_AD_INITIALIZE
    EmployF.initialize();
  #endif
  EmployW.allocate(last_yr,last_yr+NProj+1,"EmployW");
  #ifndef NO_AD_INITIALIZE
    EmployW.initialize();
  #endif
  EmployFFF.allocate(last_yr,last_yr+NProj+1,"EmployFFF");
  #ifndef NO_AD_INITIALIZE
    EmployFFF.initialize();
  #endif
  EmployFFW.allocate(last_yr,last_yr+NProj+1,"EmployFFW");
  #ifndef NO_AD_INITIALIZE
    EmployFFW.initialize();
  #endif
  EmployT.allocate(last_yr,last_yr+NProj+1,"EmployT");
  #ifndef NO_AD_INITIALIZE
    EmployT.initialize();
  #endif
  PV.allocate(last_yr,last_yr+NProj+1,"PV");
  #ifndef NO_AD_INITIALIZE
    PV.initialize();
  #endif
  TempX.allocate(last_yr,last_yr+NProj+1,"TempX");
  #ifndef NO_AD_INITIALIZE
    TempX.initialize();
  #endif
  TempF.allocate(last_yr,last_yr+NProj+1,"TempF");
  #ifndef NO_AD_INITIALIZE
    TempF.initialize();
  #endif
  TempW.allocate(last_yr,last_yr+NProj+1,"TempW");
  #ifndef NO_AD_INITIALIZE
    TempW.initialize();
  #endif
  PresValue.allocate("PresValue");
  #ifndef NO_AD_INITIALIZE
  PresValue.initialize();
  #endif
  CPUE_Like.allocate("CPUE_Like");
  #ifndef NO_AD_INITIALIZE
  CPUE_Like.initialize();
  #endif
  Survey_Like.allocate("Survey_Like");
  #ifndef NO_AD_INITIALIZE
  Survey_Like.initialize();
  #endif
  CAA_Likelihood.allocate("CAA_Likelihood");
  #ifndef NO_AD_INITIALIZE
  CAA_Likelihood.initialize();
  #endif
  CAAS_Likelihood.allocate("CAAS_Likelihood");
  #ifndef NO_AD_INITIALIZE
  CAAS_Likelihood.initialize();
  #endif
  RecRes_Likelihood.allocate("RecRes_Likelihood");
  #ifndef NO_AD_INITIALIZE
  RecRes_Likelihood.initialize();
  #endif
  Oneyearold_Likelihood.allocate("Oneyearold_Likelihood");
  #ifndef NO_AD_INITIALIZE
  Oneyearold_Likelihood.initialize();
  #endif
  count1.allocate("count1");
  #ifndef NO_AD_INITIALIZE
  count1.initialize();
  #endif
  count2.allocate("count2");
  #ifndef NO_AD_INITIALIZE
  count2.initialize();
  #endif
  count4.allocate("count4");
  #ifndef NO_AD_INITIALIZE
  count4.initialize();
  #endif
  count3.allocate("count3");
  #ifndef NO_AD_INITIALIZE
  count3.initialize();
  #endif
  count5.allocate("count5");
  #ifndef NO_AD_INITIALIZE
  count5.initialize();
  #endif
  countAll.allocate("countAll");
  #ifndef NO_AD_INITIALIZE
  countAll.initialize();
  #endif
  Akaike.allocate("Akaike");
  #ifndef NO_AD_INITIALIZE
  Akaike.initialize();
  #endif
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  KspSTD.allocate("KspSTD");
  Pred_Rec.allocate(1,20,"Pred_Rec");
  hSTD.allocate("hSTD");
  #ifndef NO_AD_INITIALIZE
  hSTD.initialize();
  #endif
  Cur_90.allocate("Cur_90");
  Cur_B0.allocate("Cur_B0");
  Cur_Bmsy.allocate("Cur_Bmsy");
  aveRY_90.allocate("aveRY_90");
  aveRY_last5.allocate("aveRY_last5");
  KexpSTD.allocate("KexpSTD");
  #ifndef NO_AD_INITIALIZE
  KexpSTD.initialize();
  #endif
  Bspcur.allocate("Bspcur");
  #ifndef NO_AD_INITIALIZE
  Bspcur.initialize();
  #endif
  Bexpcur.allocate("Bexpcur");
  #ifndef NO_AD_INITIALIZE
  Bexpcur.initialize();
  #endif
  Spmsy.allocate("Spmsy");
  #ifndef NO_AD_INITIALIZE
  Spmsy.initialize();
  #endif
  Bexpmsy.allocate("Bexpmsy");
  #ifndef NO_AD_INITIALIZE
  Bexpmsy.initialize();
  #endif
  MSY.allocate("MSY");
  #ifndef NO_AD_INITIALIZE
  MSY.initialize();
  #endif
  DepspSTD.allocate("DepspSTD");
  #ifndef NO_AD_INITIALIZE
  DepspSTD.initialize();
  #endif
  DepexSTD.allocate("DepexSTD");
  #ifndef NO_AD_INITIALIZE
  DepexSTD.initialize();
  #endif
  Bcur_Bspmsy.allocate("Bcur_Bspmsy");
  #ifndef NO_AD_INITIALIZE
  Bcur_Bspmsy.initialize();
  #endif
  Bcur_Bexpmsy.allocate("Bcur_Bexpmsy");
  #ifndef NO_AD_INITIALIZE
  Bcur_Bexpmsy.initialize();
  #endif
  MSYL_Ksp.allocate("MSYL_Ksp");
  #ifndef NO_AD_INITIALIZE
  MSYL_Ksp.initialize();
  #endif
  MSYL_Kexp.allocate("MSYL_Kexp");
  #ifndef NO_AD_INITIALIZE
  MSYL_Kexp.initialize();
  #endif
  MSTD.allocate(0,plus_grp,"MSTD");
  #ifndef NO_AD_INITIALIZE
    MSTD.initialize();
  #endif
  Spre1.allocate(0,plus_grp,"Spre1");
  #ifndef NO_AD_INITIALIZE
    Spre1.initialize();
  #endif
  Spre2.allocate(0,plus_grp,"Spre2");
  #ifndef NO_AD_INITIALIZE
    Spre2.initialize();
  #endif
  Spost.allocate(0,plus_grp,"Spost");
  #ifndef NO_AD_INITIALIZE
    Spost.initialize();
  #endif
  Ss.allocate(0,plus_grp,"Ss");
  #ifndef NO_AD_INITIALIZE
    Ss.initialize();
  #endif
  qCPU.allocate(1,NCpueSeries,"qCPU");
  #ifndef NO_AD_INITIALIZE
    qCPU.initialize();
  #endif
  qSurvPre.allocate(1,NSurveySeries,"qSurvPre");
  #ifndef NO_AD_INITIALIZE
    qSurvPre.initialize();
  #endif
  qSurvPost.allocate(1,NSurveySeries,"qSurvPost");
  #ifndef NO_AD_INITIALIZE
    qSurvPost.initialize();
  #endif
  qSeal.allocate("qSeal");
  #ifndef NO_AD_INITIALIZE
  qSeal.initialize();
  #endif
  SigCPU.allocate(1,NCpueSeries,"SigCPU");
  #ifndef NO_AD_INITIALIZE
    SigCPU.initialize();
  #endif
  SigCAA_com.allocate("SigCAA_com");
  #ifndef NO_AD_INITIALIZE
  SigCAA_com.initialize();
  #endif
  SigCAA_surv.allocate(1,NSurveySeries,"SigCAA_surv");
  #ifndef NO_AD_INITIALIZE
    SigCAA_surv.initialize();
  #endif
  SigSeal_Index.allocate("SigSeal_Index");
  #ifndef NO_AD_INITIALIZE
  SigSeal_Index.initialize();
  #endif
  AddvarSTD.allocate("AddvarSTD");
  #ifndef NO_AD_INITIALIZE
  AddvarSTD.initialize();
  #endif
  Fmsy.allocate("Fmsy");
  #ifndef NO_AD_INITIALIZE
  Fmsy.initialize();
  #endif
  Rmsy.allocate("Rmsy");
  #ifndef NO_AD_INITIALIZE
  Rmsy.initialize();
  #endif
  RYSTD.allocate("RYSTD");
  #ifndef NO_AD_INITIALIZE
  RYSTD.initialize();
  #endif
  Bmsy.allocate("Bmsy");
  #ifndef NO_AD_INITIALIZE
  Bmsy.initialize();
  #endif
  MSYL_K.allocate("MSYL_K");
  #ifndef NO_AD_INITIALIZE
  MSYL_K.initialize();
  #endif
  RanNum.allocate("RanNum");
  #ifndef NO_AD_INITIALIZE
  RanNum.initialize();
  #endif
  r.allocate("r");
  #ifndef NO_AD_INITIALIZE
  r.initialize();
  #endif
  DepSTD.allocate(first_yr,last_yr,"DepSTD");
  #ifndef NO_AD_INITIALIZE
    DepSTD.initialize();
  #endif
  negpen.allocate("negpen");
  #ifndef NO_AD_INITIALIZE
  negpen.initialize();
  #endif
  prior.allocate("prior");
  #ifndef NO_AD_INITIALIZE
  prior.initialize();
  #endif
  sel_pen.allocate("sel_pen");
  #ifndef NO_AD_INITIALIZE
  sel_pen.initialize();
  #endif
  neg2.allocate("neg2");
  #ifndef NO_AD_INITIALIZE
  neg2.initialize();
  #endif
  Bstd.allocate(first_yr,last_yr,"Bstd");
  Rstd.allocate(first_yr,last_yr,"Rstd");
  Fstd.allocate(first_yr,last_yr,"Fstd");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  int Age,Year,Iser,AMinus,APlus;
  double Temp;
  
  //**Clean up commercial CAA data**
  //---------------------------------
  for (Year=first_yr;Year<=last_yr;Year++)
  {
     //** Total the catch-at-age and store the total**      
     //---------------------------------------------
     Temp = 0;
     for (Age=0;Age<=plus_grp;Age++)
      Temp += CAA(Year,Age);
     CAA(Year,-1) = Temp; 
     if (Temp > 0)
      {
              
       // Implement the minus group
       //-------------------------------
       for (Age=0;Age<CAAMinus;Age++)
        {
          CAA(Year,CAAMinus) += CAA(Year,Age);
          CAA(Year,Age) = 0;
        }
      
      //** Implement the plus group**
      //-------------------------------
       for (Age=plus_grp;Age>CAAPlus;Age--)
        {
          CAA(Year,CAAPlus) += CAA(Year,Age);
          CAA(Year,Age) = 0;
        }
       //** Renormalise the CAA data**
       //----------------------------
        for (Age=0;Age<=plus_grp;Age++)
         CAA(Year,Age) = CAA(Year,Age) / Temp;
      }
    }
   
   //-----------------------------------------
   for (Iser=1; Iser<=NSurveySeries;Iser++)
    for (Year=first_yr;Year<=last_yr;Year++)
     {
       //** Total the catch-at-age and store the total**
      //--------------------------------------------
      
      Temp = 0;
      for (Age=0;Age<=plus_grp;Age++)
        Temp += SurvCAA(Iser,Year,Age);
      SurvCAA(Iser,Year,-1) = Temp;
      //cout<<Year<<" "<<SurvCAA(Iser,Year,-1)<<endl; 
      if (Temp>0)
       {
        //** Implement the minus group**
        //--------------------------
        AMinus = CAASMinus(Iser);
        for (Age=0;Age<AMinus;Age++)
         {
           SurvCAA(Iser,Year,AMinus) +=SurvCAA(Iser,Year,Age);
           SurvCAA(Iser,Year,Age) = 0;
         }
        
        // **Implement the minus group**
        //------------------------------
        APlus = CAASPlus(Iser);
        for (Age=plus_grp;Age>APlus;Age--)
         {
           SurvCAA(Iser,Year,APlus) += SurvCAA(Iser,Year,Age);
           SurvCAA(Iser,Year,Age) = 0;
         }
     
        //** Renormalise the CAA data**
        //-------------------------
        if (Temp > 0)
        for (Age=0;Age<=plus_grp;Age++)
          SurvCAA(Iser,Year,Age) = SurvCAA(Iser,Year,Age) / Temp;
       }
    }
     
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  SAPred.initialize();
  CAPred.initialize();
  CAPred2.initialize();
  int Age,Year;
  negpen.initialize();
  Specify_Select();
     //cout<<"Specify Selectivity ok"<<endl;
  Specify_M();
     //cout<<"Specify M ok"<<endl;
  Set_Recruitment_Residuals();
     //cout<<"Set Recruitment Residuals ok"<<endl;
  Calculate_Initial_Age_Structure();
     //cout<<"Calculate initial age-structure ok"<<endl;
  Project_forward();
     //cout<<"Project forward ok"<<endl;
  Calc_CPUE_Like();
     //cout<<"CPUE likelihood ok"<<endl;
  Calc_Survey_Like();
     //cout<<"Survey likelihood ok"<<endl;
  Calc_CAA_Likelihood();
     //cout<<"CAA likelihood ok"<<endl;
  Calc_CAAS_Likelihood();  
     //cout<<"Survey CAA likelihood ok"<<endl;
  Calc_RecRes_Likelihood();
     //cout<<"Recruitment residual likelihood"<<endl;
  Calc_Oneyearold_Likelihood();
    //cout<<"Oneyearold_Likelihood"<<endl;
  Get_SRCurve();
  Get_stdev();
  if (sd_phase()) 
  {
    Get_MSY();
  }
    //cout<<"Get MSY ok"<<endl;
     //cout<<"Get sd report numbers ok"<<endl;
  // Slight penalty on changes of selex over time
  sel_pen.initialize();
  for (int SelPeriod=1;SelPeriod<=NSelPeriods;SelPeriod++) {
    int Year = SelYrs(SelPeriod,1);
  sel_pen = 12.5*norm2(log(S(Year)-log(S(Year+1))));
  }
  obj_fun = CPUE_Like + Survey_Like + CAA_Likelihood + CAAS_Likelihood + RecRes_Likelihood + 
  weight*Oneyearold_Likelihood + negpen + sel_pen  + prior;
  //f = CPUE_Like + Survey_Like + CAA_Likelihood + CAAS_Likelihood + RecRes_Likelihood + negpen;
   Get_MSY();
  if (mceval_phase())
  {
     Do_Projections();
    ofstream out6("CurDep1990.out",ios::app); 
    for (Year=1990; Year<=last_yr; Year++)
      out6<<TotalB(last_yr)/TotalB(1990)<< " ";
    out6<<" "<<endl;
  }
  // Nsim = Nsim +1;
  // r = N(1985,0);
  // random_number_generator rndRec(3+Nsim);
  // RanNum = randu(rndRec);
#ifdef DEBUG
  std::cout << "DEBUG: Total gradient stack used is " << gradient_structure::get()->GRAD_STACK1->total() << " out of " << gradient_structure::get_GRADSTACK_BUFFER_SIZE() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->GRAD_LIST->total_addresses() << " out of " << gradient_structure::get_MAX_DLINKS() << std::endl;;
  std::cout << "DEBUG: Total dvariable address used is " << gradient_structure::get()->ARR_LIST1->get_max_last_offset() << " out of " << gradient_structure::get_ARRAY_MEMBLOCK_SIZE() << std::endl;;
#endif
}

void model_parameters::Specify_Select(void)
{
  int Year,Age,Yr1,Yr2, SelPeriod;
  dvariable Temp, Rage, slope, intc, Smax, RageS, SmaxS, TempS;
  // ==+==+==+== COMMERCIAL SELECTIVITY ==+==+==+==
  // Work through each selectivity period, first set up the selectivity
  // vector for the first year of the period, then copy it to all years
  for (SelPeriod=1;SelPeriod<=NSelPeriods;SelPeriod++)
  {
    // **Set up logistic selectivity for the first year of the period**
    //-------------------------------------------------------------
    Year = SelYrs(SelPeriod,1);
    Temp = log(19.0)/Sage(2*SelPeriod);
    Rage = 0;
    Smax = 0;
    for (Age=0; Age<=plus_grp; Age++)
    {
      S(Year,Age) = 1.0/(1+mfexp(-1*Temp*(Rage-Sage(2*SelPeriod-1))));
      if (Age >= Age_Sel_slope) S(Year,Age) = S(Year,Age) * mfexp(-SelSlope(SelPeriod)*(Age-Age_Sel_slope));
      if (S(Year,Age) > Smax) Smax = S(Year,Age);
      Rage = Rage + 1;
    }
    for (Age=0;Age<=plus_grp;Age++)
      S(Year,Age) = S(Year,Age) / Smax;
    // Copy to the rest of the period
    for (Year=SelYrs(SelPeriod,1)+1; Year<=SelYrs(SelPeriod,2); Year++)
     for (Age=0;Age<=plus_grp;Age++)
      S(Year,Age) = S(SelYrs(SelPeriod,1),Age);
  }
  for (SelPeriod=1;SelPeriod<=NSelPeriods-1;SelPeriod++)
  {
    for (Age=0;Age<=plus_grp;Age++)
     {
      Yr1 = SelYrs(SelPeriod,2);
      Yr2 = SelYrs(SelPeriod+1,1);
      slope = (S(Yr1,Age)-S(Yr2,Age))/(Yr1-Yr2);
      intc = S(Yr1,Age) - slope * Yr1;
      for (Year=SelYrs(SelPeriod,2)+1;Year<=SelYrs(SelPeriod+1,1)-1;Year++)
       S(Year,Age) = intc + slope * Year;
    }
    for (Year=SelYrs(SelPeriod,2)+1;Year<=SelYrs(SelPeriod+1,1)-1;Year++)
    {
      Smax = 0;
      for (Age=0;Age<=plus_grp;Age++)
       if (S(Year,Age) > Smax) Smax = S(Year,Age);
      for (Age=0;Age<=plus_grp;Age++)
       S(Year,Age) = S(Year,Age) / Smax;
    }
  } 
  //Set up the selectivity for future projections
  for (Age=0; Age<= plus_grp; Age++)
   Sfut(Age) = S(last_yr,Age);
  // ==+==+==+== SURVEY SELECTIVITY ==+==+==+==
  TempS = log(19.0)/SurvAd;
  RageS = 0;
  SmaxS = 0;
  if (select == 1)
  {
    for (Age=0; Age<=plus_grp; Age++)
    {
      SurvS(Age) = 1.0/(1+mfexp(-1*TempS*(RageS-SurvA50)));
      if (SurvS(Age) > SmaxS) SmaxS = SurvS(Age);
      RageS = RageS + 1;
    }
  }
  if (select == 0)
  {
    for (Age=0; Age<=plus_grp; Age++)
    {
      SurvS(Age) = 1.0/(1+mfexp(-1*TempS*(RageS-SurvA50)));
      if (Age >= Age_Sel_slopeS) 
        SurvS(Age) = SurvS(Age) * mfexp(-SelSlopeS*(Age-Age_Sel_slopeS));
      if (SurvS(Age) > SmaxS) 
        SmaxS = SurvS(Age);
      RageS = RageS + 1;
   }
  }
  for (Age=0;Age<=plus_grp;Age++)
    SurvS(Age) = SurvS(Age) / SmaxS; 
}

void model_parameters::Specify_M(void)
{
  int Age;
  if (est_M >0)
  {
   for (Age=0; Age<=plus_grp; Age++)
    M(Age) = M1;
  }
  else
  {
  for (Age=1; Age<=plus_grp; Age++)
    M(Age) = Minfage*Minf/pow(Age, 0.321921041); // Designed to give M=0.5 at age 3 and 0.4 by age 6
    // M(Age) = Minf * Minfage/Age;
  M(0) = M(1) * 2;
  }
}

void model_parameters::Set_Recruitment_Residuals(void)
{
  int Year;
  for (Year=first_yr; Year<RecYr1; Year++)
    Recruit_Res(Year).initialize();
  for (Year=RecYr1; Year<=RecYr2; Year++)
  {
    Recruit_Res(Year) = RecPar(Year);
  }
  for (Year=RecYr2+1; Year<=last_yr+1; Year++)
    Recruit_Res(Year).initialize() ;
}

void model_parameters::Calculate_Initial_Age_Structure(void)
{
  int Age, year, strt1,end1,strt2,end2,strt3,end3;
  dvariable SPR0,R0;
  dvariable slope, intc,Ni;
  // **First get the VIRGIN age-structure**
  //-----------------------------------
  Nvirg(0) = 1;
  for (Age=1; Age<=plus_grp; Age++)
    Nvirg(Age) = Nvirg(Age-1)*mfexp(-M(Age-1));
  Nvirg(plus_grp) = Nvirg(plus_grp) / (1.0 - mfexp(-M(plus_grp)));
  // **Find the virgin spawner biomass-per-recruit and hence R0**
  //-----------------------------------------------------------
    //model two Ksp -periods
  //---------------------------------------------------------------
  strt1 = KspYrs(NKspPeriods-1,1);
  end1  = KspYrs(NKspPeriods-1,2);
  strt2 = KspYrs(NKspPeriods,1);
  end2  = KspYrs(NKspPeriods,2);
  strt3 = (KspYrs((NKspPeriods-1),2))+1;
  end3  = (KspYrs(NKspPeriods,1))-1;
  for (year=strt1;year<=end1;year++)
  {
    SPR0 = 0;
    for (Age=1; Age<=plus_grp; Age++)
      SPR0 += mat(Age)*Wstrt(KspYrs(NKspPeriods-1,1),Age)*Nvirg(Age);
    R0          = mfexp(par_B0)/SPR0;
    Alpha(year) = (4.0 * Steep * R0)/(5 * Steep - 1.0);
    Beta(year)  = (mfexp(par_B0))*(1.0 - Steep)/ (5 * Steep - 1.0);
  }
  for (year=strt2;year<=end2;year++)
  {
    SPR0 = 0;
    for (Age=1; Age<=plus_grp; Age++)
     SPR0      += mat(Age)*Wstrt(KspYrs(NKspPeriods,1),Age)*Nvirg(Age);
    R0          = mfexp(par_B0)*Prop/SPR0;
    Alpha(year) = (4.0 * Steep * R0)/(5 * Steep - 1.0);
    Beta(year)  = (mfexp(par_B0))*Prop*(1.0 - Steep)/ (5 * Steep - 1.0);
  }
  // ** Interpolate linearly between different Ksp Periods **
  // i.e. for years in between set Ksp periods, set the alpha and beta by linear interpolation
  //------------------------------------------------------------- 
  for (year=strt3;year<=end3;year++)
  {
    slope       = 0;
    intc        = 0;
    slope       = (Alpha(KspYrs(NKspPeriods-1,2))-Alpha(KspYrs(NKspPeriods,1)))/(KspYrs(NKspPeriods-1,2)-KspYrs(NKspPeriods,1));
    intc        = Alpha(KspYrs(NKspPeriods-1,2)) - slope * KspYrs(NKspPeriods-1,2);
    Alpha(year) = intc + slope * year;
    slope       = 0;
    intc        = 0;
    slope       = (Beta(KspYrs(NKspPeriods-1,2))-Beta(KspYrs(NKspPeriods,1)))/(KspYrs(NKspPeriods-1,2)-KspYrs(NKspPeriods,1));
    intc        = Beta(KspYrs(NKspPeriods-1,2)) - slope * KspYrs(NKspPeriods-1,2);
    Beta(year)  = intc + slope * year;
  }
  //-----------------------------------------------------------  
  Spawn(first_yr)  = 0.; 
  TotalB(first_yr) = 0.;
  BMid(first_yr)   = 0.;
  for (Age=0;Age<=plus_grp; Age++)
  {
    N(first_yr,Age)   = mfexp(par_B0)/SPR0 * Nvirg(Age);
    Ni               += N(first_yr,Age);                                                  //Add numbers of fish for first year
    Spawn(first_yr)  += N(first_yr,Age)*Wstrt(first_yr,Age)*mat(Age);
    TotalB(first_yr) += (N(first_yr,Age)*Wstrt(first_yr,Age));
    BMid(first_yr)   += N(first_yr,Age)*Wmid(first_yr,Age)*mfexp(-M(Age)/2);
  }
  // ** Set the spawner biomass at the start of year 1 **
  //-------------------------------------------------------------
  Spawn(first_yr)    = Ksp_fraction*mfexp(par_B0);
  Recruits(first_yr) = N(first_yr,0)*Wstrt(first_yr,0)*(1-mat(0));
}

void model_parameters::Project_forward(void)
{
  int Year, Age, II;
  dvariable Fish,Fish2,Temp;
  // **DO for each year of the projection period**
  //---------------------------------------------------
  for (Year=first_yr; Year<=last_yr; Year++)
  {
    //** Begin-year exploitable and survey biomasses**
    //----------------------------------------------
    BBeg(Year) = 0; SurvBeg(Year) = 0; TotalB(Year) = 0;
    for (Age=1; Age<=plus_grp;Age++)
    {
      Temp          = Wstrt(Year,Age)*N(Year,Age);
      BBeg(Year)    = BBeg(Year) + S(Year,Age)*Temp;
      SurvBeg(Year) = SurvBeg(Year) + SurvS(Age)*Temp;
      TotalB(Year)  = TotalB(Year)+Temp;
    }
    // ** Make sure pop do not go negative **
    // --------------------------------------
    BBeg(Year)    = posfun(BBeg(Year),0.10,negpen);
    SurvBeg(Year) = posfun(SurvBeg(Year),0.10,negpen);
    // Mild penalty to tame Sage
		prior  = 2*norm2(Sage - 2.);
		// prior += 12.5*norm2(SelSlope - 1.);
		// prior += square(SurvA50 - 4.);
    // **First remove half of natural mortality**
    //--------------------------------------------
    for (Age=0; Age<=plus_grp; Age++)
    {
      Mtemp(Age) = N(Year,Age)*(1-mfexp(-M(Age)/2.0));   //Number of fish that died due to Natural mortality by midyear
      Ntemp(Age) = N(Year,Age)*mfexp(-M(Age)/2.0);
    }
    //** Exploitable biomass and fishing mortality**
    //----------------------------------------------
    Bexp(Year) = 0;
    for (Age=0; Age<=plus_grp; Age++)
    {
       Bexp(Year) = Bexp(Year) + Wmid(Year,Age)*S(Year,Age)*Ntemp(Age);
    }
    Fish = Catch(Year)/Bexp(Year);
    if (Fish>0.95)
    {
       //cout<<"!!! Fishing mortality greater than 0.95 in year "<<Year<<endl;
       negpen += square(Fish - 0.95)*100.;
       Fish    = 0.95;
    }
    F(Year) = Fish;
    //** Remove the catch and compute survey selected biomass**
    //---------------------------------------------------------
    SurvMid(Year) = 0;
    for (Age=0; Age <= plus_grp; Age++)
     {
       CAPred(Year,Age) = Ntemp(Age)*S(Year,Age)*Fish;
       if (Year<=1988)
         CAPred2(Year,Age) = Ntemp(Age)*S(Year,Age)*Fish2;
       else
         CAPred2(Year,Age) = Ntemp(Age)*S(1988,Age)*Fish2;
       SurvMid(Year) = SurvMid(Year) + SurvS(Age)*Wmid(Year,Age)*Ntemp(Age)*(1.0-S(Year,Age)*Fish/2.0);
       Ntemp(Age)    = Ntemp(Age) - CAPred(Year,Age);
       Ntemp2(Age)   = Ntemp(Age) - CAPred2(Year,Age);
     }
    //** Adjust Bexp to the middle of the year**
    //------------------------------------------
    Bexp(Year)   = Bexp(Year) * (1.0 - Fish/2.0);
    Bexp(Year)   = posfun(Bexp(Year),0.10,negpen);
    TotalB(Year) = posfun(TotalB(Year),0.10,negpen);
    // **Remove the rest of natural mortality**
    //------------------------------------------
    for (Age=0; Age<=plus_grp; Age++)
     {
      Ntemp(Age) = Ntemp(Age)*mfexp(-M(Age)/2.0);
      Mtemp(Age) = Mtemp(Age) + Ntemp(Age)*(1-mfexp(-M(Age)/2.0));
     }
    // **Update the age-structure**
    //--------------------------------
    N(Year+1,plus_grp)  = Ntemp(plus_grp)+Ntemp(plus_grp-1);
    NM(Year+1,plus_grp) = Mtemp(plus_grp)+Mtemp(plus_grp-1);
    for (Age=1;Age<plus_grp;Age++)
     {
       N(Year+1,Age)  = Ntemp(Age-1);
       NM(Year+1,Age) = Mtemp(Age-1);
     }
    // **Find the biomass by age that died due to M for J-P**
    //------------------------------------------------
   for (Age=0; Age<=plus_grp; Age++)
    {
     Mbio(Year,Age) = Wstrt(Year,Age)*NM(Year,Age);
     //cout<<Year<<"  "<<Age<<" "<<Mbio(Year+1,Age)<<endl;
    }
    // **Find the spawner biomass**
    //-------------------------------
    Spawn(Year+1) = 0;
    for (Age=1; Age<=plus_grp; Age++)
     Spawn(Year+1) = Spawn(Year+1) + mat(Age)*Wstrt(Year,Age)*N(Year+1,Age);
     Spawn(Year+1) = posfun(Spawn(Year+1),0.10,negpen);
    //** Generate the recruitment**
    //---------------------------------------
     N(Year+1,0) = ((Alpha(Year)*Spawn(Year+1))/(Beta(Year) + Spawn(Year+1)))*mfexp(Recruit_Res(Year));
  }
   //cout<<" here1  "<<endl;
}

void model_parameters::Calc_Oneyearold_Likelihood(void)
{
  int Year,Iser;
  dvariable qtot1;
  int nobs;
  dvariable ss;
  dvariable SSQ;
  dvariable Error;
  dvariable CVSeal;
  // This estimates recruitment from the sealscat data
     Oneyearold_Likelihood = 0;
     qtot1 = 0; nobs = 0; 
     for (Year=first_yr; Year<=last_yr-2; Year++)
     {
        seali_pred(Year) = N(Year,0) ;
        if (seali(Year) > 0)
        {
           nobs++;
           qtot1 += log(seali(Year)/seali_pred(Year));
        }
      }
      qtot1 = mfexp(qtot1/nobs);
     // Now find the likelihood
     //-------------------------
     ss.initialize();
     SSQ.initialize();
     for (Year=first_yr; Year<=last_yr; Year++)
     {
       if (seali(Year) > 0)
       {
         CVSeal = sealCV(Year) + sig_seal_added;
         Error  = log(seali(Year)) - log(qtot1* seali_pred(Year));
         Error *= Error;
         SSQ   += Error;
         ss    += Error/(2.0*CVSeal*CVSeal) + log(CVSeal);
       }
     }
      qSeal = qtot1;
      SigSeal_Index = sqrt(SSQ/nobs);
      Oneyearold_Likelihood = ss;
}

void model_parameters::Calc_CPUE_Like(void)
{
  int Year,Iser;
  dvariable qtot;
  int nobs;
  dvariable ss;
  dvariable Error;
  dvariable Sigma;
  dvariable Temp;
  // Initialise CPUE likelihood
  CPUE_Like.initialize() ;
  count1=0;
  for (Iser=1; Iser<=NCpueSeries; Iser++)
  {
    // Find ML for q
    qtot.initialize();
    double sigma_temp   = CPUE_Sigma(Iser);
    double sigmasq_temp = sigma_temp*sigma_temp;
    nobs = 0;
    for (Year=first_yr; Year<=last_yr; Year++)
    {
      if (CPUE(Year,Iser) > 0)
      {
        nobs++; 
        if (CpueIndx(Iser) == 1) BIO(Year) = Bexp(Year);
        if (CpueIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
        if (CpueIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
        qtot += log(CPUE(Year,Iser)/BIO(Year));
      }
    }
    qtot = mfexp(qtot/nobs);
      // Now find the sum of squares bit
    ss.initialize();
    for (Year=first_yr; Year<=last_yr; Year++)
    {
      if (CPUE(Year,Iser) > 0)
      {
        Error = log(CPUE(Year,Iser)) - log(qtot*BIO(Year));
        ss   += Error*Error;
      }
    } 
    // Find the likelihood component
    Sigma = sqrt(ss/nobs);
    if (UseCPUE(Iser) > 0) 
    {
      if (sigma_temp>0)
        Temp =  0.5 * ss /sigmasq_temp + nobs*log(sigma_temp); // Use assumed sigmas input
      else
        Temp = nobs*log(Sigma) + nobs/2.0;                    // Estimate sigmas internally
      // Account for CPUE weighting
      CPUE_Like += CpueWght(Iser)*Temp;
    } 
    // Store variables for output
    qCPU(Iser) = qtot;
    SigCPU(Iser) = Sigma;
    count1=nobs;
  }  
}

void model_parameters::Calc_Survey_Like(void)
{
  int Year,Iser;
  dvariable qtot1;
  dvariable qtot2;
  dvariable den1;
  dvariable den2;
  dvariable dez1;
  dvariable dez2;
  dvariable temp;
  dvariable temp1;
  dvariable temp2;
  dvariable ss;
  dvariable Error;
  dvariable CVSurv2;
  // Initialise Survey likelihood
  Survey_Like = 0;
  count2=0;
  for (Iser=1; Iser<=NSurveySeries; Iser++)
   {
    // First find the survey-q
    if (q_input(Iser) > 0)
    {
      qtot1 = q_input(Iser);
      qtot2 = q_input(Iser);
    }
    else
    {
      qtot1 = 0; qtot2 = 0; den1 = 0; den2 = 0;
      for (Year=first_yr; Year<=last_yr; Year++)
       if (Survey(Year,Iser*2-1) > 0)
       {
         if (SurveyIndx(Iser) == 1) BIO(Year) = Bexp(Year);
         if (SurveyIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
         if (SurveyIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
         CVSurv2 = Addvar+pow(Survey(Year,Iser*2),2.0);
         if (Year <= q_chngYr(Iser))
         {
           qtot1 = qtot1 + log(Survey(Year,Iser*2-1)/BIO(Year))/CVSurv2;
           den1  = den1 + 1.0/CVSurv2;
         }
         else  
         {
           qtot2 = qtot2 + log(Survey(Year,Iser*2-1)/BIO(Year))/CVSurv2;
           den2  = den2 + 1.0/CVSurv2;
         }
       }
       dez1  = -1.0/pow(q_CV(Iser),2.0);
       dez2  = -1.0/pow(q_CV(Iser),2.0);
       den1  = den1 + 1.0/pow(q_CV(Iser),2.0);
       den2  = den2 + 1.0/pow(q_CV(Iser),2.0);
       qtot1 = qtot1 + q_mean(Iser)/pow(q_CV(Iser),2.0);
       qtot2 = qtot2 - q_mean(Iser)/pow(q_CV(Iser),2.0);
       temp  = den1*den2 - dez1*dez2;
       temp1 = den2*qtot1 - dez1*qtot2;
       temp2 = den1*qtot2 - dez2*qtot1;
       qtot1 = mfexp(temp1/temp);
       qtot2 = mfexp(temp2/temp);
     } 
     // Now find the sum of squares bit
     ss = 0;
     for (Year=first_yr; Year<=last_yr; Year++)
      if (Survey(Year,Iser*2-1) > 0)
      {
        count2 += count2 +1;
        if (SurveyIndx(Iser) == 1) BIO(Year) = Bexp(Year);
        if (SurveyIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
        if (SurveyIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
        CVSurv2 = Addvar+pow(Survey(Year,Iser*2),2.0);
        if (Year <= q_chngYr(Iser)) 
          Error = log(Survey(Year,Iser*2-1)) - log(qtot1*BIO(Year));
        else
          Error = log(Survey(Year,Iser*2-1)) - log(qtot2*BIO(Year));
        ss = ss + Error*Error/(2.0*CVSurv2) + log(sqrt(CVSurv2));
     }
     qSurvPre(Iser)  = qtot1;
     qSurvPost(Iser) = qtot2;
     if (UseSurvey(Iser) > 0)
     {
       Survey_Like = Survey_Like + ss;
       Error = log(qtot1) - q_mean(Iser) - log(qtot2);
       Survey_Like += 0.5*Error*Error/(q_CV*q_CV);
  //   cout << Survey_Like << endl;
     }
   }
}

void model_parameters::Calc_CAA_Likelihood(void)
{
  int Year,Age;
  int Sum1;
  dvariable Sum2;
  dvariable Sum3;
  dvariable Sum4;
  dvariable Sigma;
  dvariable Total;
  dvariable Residual;
  // Initial CAA Likelihood
  CAA_Likelihood = 0;
  Sum1 = 0; Sum2 = 0; Sum3 = 0;Sum4.initialize();
  count3=0;
  for (Year=first_yr;Year<=last_yr;Year++)
   if (CAA(Year,-1) > 0)
    {
     count3+=count3+1;
     //**Create a predicted minus group**
     //---------------------------------
     for (Age=0;Age<CAAMinus;Age++)
      {
       CAPred(Year,CAAMinus) += CAPred(Year,Age);
       // CAPred(Year,Age) = 0; // Wrong
      }
     //**Create a predicted plus group**
     //---------------------------------
     for (Age=plus_grp;Age>CAAPlus;Age--)
      {
       CAPred(Year,CAAPlus) += CAPred(Year,Age);
       // CAPred(Year,Age) = 0; // Wrong
      }
      Total = 0;
      // for (Age=CAAMinus;Age<=CAAPlus;Age++)
      // Total += CAPred(Year,Age);
     //**Renormalize the CAA data**
      CAPred(Year)(CAAMinus,CAAPlus) /= sum(CAPred(Year)(CAAMinus,CAAPlus)) ;
     //-----------------------------
     // for (Age=CAAMinus;Age<=CAAPlus;Age++)
      // CAPred(Year,Age) /= Total;
     // **Now accumulate the total variables (Sum1 is a counter)**
     //----------------------------------------------------------
     for (Age=CAAMinus;Age<=CAAPlus;Age++)
      {
        if (CAPred(Year,Age) >0) 
         if (CAA(Year,Age) > 0)
          {
            Sum1++;
            if (use_multinomial>0)
            {
              Sum4 -= CAA(Year,-2)*CAA(Year,Age)*log(CAPred(Year,Age) +1e-4);
            }
            Sum2 += log(CAA(Year,Age));
            Residual = log(CAA(Year,Age) / CAPred(Year,Age));
            Sum3 += CAA(Year,Age)*Residual*Residual;
          }
       }
     }
   if (Sum3/Sum1 <= 0)
     Sigma=0;
   else
     Sigma = sqrt(Sum3/Sum1);
   if (use_multinomial>0)
     CAA_Likelihood = Sum4;
   else
     CAA_Likelihood = Sum1*log(Sigma) - 0.5*Sum2 + Sum1/2.0; 
   SigCAA_com = Sigma;
}

void model_parameters::Calc_CAAS_Likelihood(void)
{
  int Year,Age,Iser,AMinus,APlus;
  int Sum1;
  dvariable Sum2;
  dvariable Sum3;
  dvariable Sum4;
  dvariable Sigma;
  dvariable Total;
  dvariable Residual;
  dvariable ZZ;
  // **Initial Survey CAA Likelihood**
  //----------------------------------
  CAAS_Likelihood = 0;
  for (Iser=1;Iser<=NSurveySeries;Iser++)
   {
    Sum1 = 0; Sum2 = 0; Sum3 = 0;Sum4.initialize();
    count4=0;
    for (Year=first_yr;Year<=last_yr;Year++)
     if (SurvCAA(Iser,Year,-1) > 0)
      {
       count4+=count4+1;
       Total = 0;
       for (Age=0;Age<=plus_grp;Age++)
        {
         //** Mid-year correction**
         //------------------------
         ZZ = mfexp(-M(Age)/2.0)*(1.0-S(Year,Age)*F(Year)/2.0);
         if (SurveyIndx(Iser) == 1) SAPred(Iser,Year,Age) = S(Year,Age)*N(Year,Age)*ZZ;
         if (SurveyIndx(Iser) == 2) SAPred(Iser,Year,Age) = SurvS(Age)*N(Year,Age)*ZZ;
         if (SurveyIndx(Iser) == 3) SAPred(Iser,Year,Age) = SurvS(Age)*N(Year,Age);
        }
       //** Now incorporate the survey "minus group" **
       //---------------------------------------
       AMinus = CAASMinus(Iser);
       for (Age=0;Age<AMinus;Age++)
       {
         SAPred(Iser,Year,AMinus) += SAPred(Iser,Year,Age);
         // SAPred(Iser,Year,Age) = 0; // Strictly speaking, this isn't allowed
       }
       //** Now incorporate the survey "plus group" **
       //---------------------------------------
       APlus = CAASPlus(Iser);
       for (Age=plus_grp;Age<APlus;Age--)
       {
         SAPred(Iser,Year,APlus) += SAPred(Iser,Year,Age);
         // SAPred(Iser,Year,Age) = 0; // Strictly speaking, this isn't allowed
       }
      //**Rescale**
      //----------------
       // for (Age=AMinus;Age<=APlus;Age++)
        // SAPred(Iser,Year,Age) /= Total;  
      SAPred(Iser,Year)(AMinus,APlus) /= sum(SAPred(Iser,Year)(AMinus,APlus)) ;
       //** Now accumulate the total variables**
       //---------------------------------------
       for (Age=AMinus;Age<=APlus;Age++)
        if (SurvCAA(Iser,Year,Age) > 0)
         if (SAPred(Iser,Year,Age)> 0)
         {
          Sum1++;
          Sum2 += log(SurvCAA(Iser,Year,Age));
          Residual = log(SurvCAA(Iser,Year,Age) / SAPred(Iser,Year,Age));
          Sum3 += SurvCAA(Iser,Year,Age)*Residual*Residual;
          if (use_multinomial>0)
          {
            Sum4 -= SurvCAA(Iser,Year,-2)*SurvCAA(Iser,Year,Age)*log(SAPred(Iser,Year,Age) +1e-4);
          }
         }
      } 
    if (Sum3/Sum1 <=0)
     Sigma = 0;
    else
    Sigma = sqrt(Sum3/Sum1);
    Total = Sum1*log(Sigma) - 0.5*Sum2 + Sum1/2.0; 
    if (use_multinomial>0)
      CAAS_Likelihood = Sum4;
    else
      CAAS_Likelihood = CAAS_Likelihood + Total;
      //cout<<" CAAS_Likelihood "<<CAAS_Likelihood<<endl;
    SigCAA_surv(Iser) = Sigma;
   }   
}

void model_parameters::Calc_RecRes_Likelihood(void)
{
  int Year;
  int nyear_RecRes;
  RecRes_Likelihood = 0.;
  count5 = 0;
  for (Year=first_yr; Year<=last_yr; Year++)
   {
    //RecRes_Likelihood += log(SigmaRec) + square(Recruit_Res(Year))/(2*square(SigmaRec));
    RecRes_Likelihood += square(Recruit_Res(Year))/(2*square(SigmaRec));
    //cout<<" Rec_Res "<<Recruit_Res(Year)<<endl;
   }
  // Calculate the sigmaR output
  nyear_RecRes = 0;
  SigR_out = 0.0;
  for (Year=RecYr1; Year<= RecYr2; Year++)
   {
    count5+=count5+1;
    nyear_RecRes += 1;
    SigR_out += square(Recruit_Res(Year));
   }
  if (SigR_out == 0)
   SigR_out = 99999;
  else
   SigR_out = sqrt(SigR_out/nyear_RecRes);
}

void model_parameters::Get_SRCurve(void)
{
 /* just fill in a curve for plotting purposes */
  // For plotting from SR curve
  SSB(20) = value( mfexp(par_B0));
  // cout<<SSB(20)<<endl;exit(1);
  for (int i=1;i<=20;i++)
  {
    SSB(i) = SSB(20)*i/20;
    Pred_Rec(i) = (SSB(i)*Alpha(last_yr))/(Beta(last_yr) + SSB(i)) *mfexp(-(SigmaRec*SigmaRec)/2.0);
  }
}

void model_parameters::Get_MSY(void)
{
 /*Function calculates used in calculating MSY and MSYL for a designated component of the
  population, given values for stock recruitment and selectivity...  Fmsy is the trial 
  value of MSY example of the use of "funnel" to reduce the amount of storage for 
  derivative calculations */
  dvariable Fdmsy;
  dvariable Stmp;
  dvariable Rtmp;
  dvariable Bexptmp;
  dvariable Btmp;
  double df=1.e-7;
  dvariable F1=.05;
  dvariable F2;
  dvariable F3;
  dvariable yld1;
  dvariable yld2;
  dvariable yld3;
  dvariable dyld;
  dvariable dyldp;
  // Newton Raphson stuff to go here
  for (int ii=1;ii<=4;ii++)
  {
    //cout<<"-----"<<ii<<"-----"<<endl;
    F2     = F1 + df*.5;
    F3     = F2 - df;
    //F1    = double(ii)/400;
    yld1   = yield(F1, Stmp,Bexptmp,Btmp);
    //cout <<ii<<" " <<F1<<" "<< Stmp <<" "<<yld1<<" "<<dyldp<<" "<< endl; }
    yld2   = yield(F2,Stmp,Bexptmp,Btmp);
    yld3   = yield(F3,Stmp,Bexptmp,Btmp);
    dyld   = (yld2 - yld3)/df;                          // First derivative (to find the root of this)
    dyldp  = (yld2 + yld3 - 2.*yld1)/(.25*df*df) + 1.e-10;   // Second derivative (for Newton Raphson)
    F1    -= dyld/dyldp;
  }
 // Reset funnel variable
  Fdmsy    = F1;
  Fmsy     = Fdmsy;
  MSY      = yield(Fmsy,Stmp,Bexptmp,Btmp);
  // cout<<      MSY<<" "<<Fmsy<<" "<<Stmp<<" "<<Bexptmp<<" "<<Btmp<<endl;
  Spmsy    = Stmp;
  Bexpmsy  = Bexptmp;
  Bmsy     = Btmp;
  Bcur_Bspmsy = Spawn(last_yr)/Spmsy;
  Bcur_Bexpmsy = Bexp(last_yr)/Bexpmsy;
  MSYL_Ksp    = Spmsy/mfexp(par_B0);
  MSYL_Kexp   = Bexpmsy/Bexp(first_yr);
  MSYL_K      = Bmsy/TotalB(first_yr);
  Rmsy     = Rtmp;
   //cout<<" MSy "<<MSY<<endl;
}

dvariable model_parameters::yield(dvariable& Ftmp, dvariable& Stmp,dvariable& Bexptmp,dvariable& Btmp)
{
  int Age, Year;
  dvariable Rtmp;
  dvariable ZZtmp;
  dvariable yield;
  dvar_vector Natmp(0,plus_grp);
  //Compute the equilibrium N
  Natmp(0) = 1.;
  for (Age=1; Age<=plus_grp; Age++)
   {
    ZZtmp = mfexp(-M(Age-1))*(1.0 - Ftmp*Sfut(Age-1));
    Natmp(Age) = Natmp(Age-1)*ZZtmp;
    //cout<<" S "<<Sfut<<endl;
   }
  ZZtmp = mfexp(-M(plus_grp))*(1.0-Ftmp*Sfut(plus_grp));
  Natmp(plus_grp) = Natmp(plus_grp)/(1.0 - ZZtmp);
  //cout<<"OK 2"<<endl;
  //Compute the model outputs
  yield = 0.0;
  Stmp = 0.0;
  Bexptmp = 0.0;
  Btmp = 0.0;
  for (Age=0; Age<=plus_grp; Age++)
   {
    yield += Wmid(last_yr,Age)*Sfut(Age)*Ftmp*Natmp(Age)*mfexp(-M(Age)/2.0);
    Stmp += mat(Age)*Wstrt(last_yr,Age)*Natmp(Age);
    ZZtmp = mfexp(-M(Age)/2.0)*(1.0-Ftmp*Sfut(Age)/2.0);
    Bexptmp += Wmid(last_yr,Age)*Sfut(Age)*Natmp(Age)*ZZtmp;
    Btmp += Wmid(last_yr,Age)*Natmp(Age);
   }
  //Compute the recruitement
  Rtmp = (Alpha(last_yr)*Stmp-Beta(last_yr))/Stmp;
  //cout<<"Yield "<<yield<<" Rtmp "<<Rtmp<<" Stmp "<<Stmp<<" Alpha "<<Alpha(last_yr)<<" Beta "<<Beta(last_yr)<<endl;
  //Adjust by recruitment 
  yield *=  Rtmp;
  Stmp *=  Rtmp;
  Bexptmp *=  Rtmp;
  Btmp *= Rtmp;
  return yield;
}

dvariable model_parameters::RY(int Year)
{
  int II,Age;
  dvariable Fmin;
  dvariable Fmax;
  dvariable FF;
  dvar_vector Ntemp(0,plus_grp);
  dvar_vector NT(0,plus_grp);
  dvariable Cpred;
  dvariable Spawno;
  dvariable test;
  dvariable RYtest;
  dvariable RYtmp;
  Fmin = 0.;
  Fmax = 10.;
  for (II = 1; II<=30; II++)
   {
    FF = (Fmin + Fmax)/2.0;
    // Project one year ahead and calculate the spawner biomass
    Cpred = 0.;
    for (Age = 0; Age<=plus_grp; Age++)
     {
      Cpred += Wmid(last_yr,Age)*Sfut(Age)*FF*N(Year,Age)*mfexp(-M(Age)/2.0);
      Ntemp(Age) = N(Year,Age)*mfexp(-M(Age))*(1.0-Sfut(Age)*FF);
     }
    // Update the age-structure
    NT(plus_grp) = Ntemp(plus_grp) + Ntemp(plus_grp-1);
    for (Age = 1; Age<=plus_grp-1; Age++)
     NT(Age) = Ntemp(Age-1);
    // Find the spawner biomass
    Spawno =0.0;
    for (Age = 1; Age<=plus_grp; Age++)
     Spawno += mat(Age)*Wstrt(last_yr,Age)*NT(Age);
    // Check for convergence and update F
    if (fabs(Spawno - Spawn(Year))<0.001)
     {
      test = 1;
      RYtest = Cpred;
     }
    if (Spawno>Spawn(Year))
     {
      Fmin = FF;
     } 
    else
     {
      Fmax = FF;
     }
  }
  if (test == 1) 
   RYtmp = RYtest;
  else
   RYtmp = Cpred;
  return RYtmp;
}

void model_parameters::Get_stdev(void)
{
  int Age, Year;
  // *** STD report ***
  KspSTD = mfexp(par_B0);
  KexpSTD = Bexp(first_yr);
  hSTD = Steep;
  Bspcur = Spawn(last_yr);      
  Bexpcur = Bexp(last_yr);      
  DepspSTD = Spawn(last_yr)/Spawn(first_yr);
  DepexSTD = Bexp(last_yr)/Bexp(first_yr);
  for (Age=0; Age<=plus_grp; Age++)
  {
    MSTD(Age) = M(Age);
    Spre1(Age) = S(1964,Age);    
    Spre2(Age) = S(1989,Age);    
    Spost(Age) = S(last_yr,Age);    
    Ss(Age)    = SurvS(Age);  
  }
	// Terminal year total biomass over 1990 value 
  Cur_90   = TotalB(last_yr)/TotalB(1990);
	// Terminal year spawning biomass over Bzero value 
  Cur_B0   = Spawn(last_yr)/mfexp(par_B0);
	// Terminal year spawning biomass over Bmsy value 
  Cur_Bmsy = Spawn(last_yr)/Spmsy;
  int no = 0; 
	aveRY_90.initialize();
  for (Year=1990;Year<=last_yr;Year++) {
    aveRY_90 += RY(Year);
    no = no + 1;
  }
  aveRY_90 /= no;
	no=0;
	aveRY_last5.initialize();
  for (Year=last_yr-5;Year<=last_yr-1;Year++) {
    aveRY_last5 += RY(Year);
    no = no + 1;
  }
	aveRY_last5 /= no;
  for (Year=first_yr; Year<=last_yr; Year++)
  {
    DepSTD(Year) = Spawn(Year)/Spawn(first_yr);
    Bstd(Year)=Spawn(Year);
    Rstd(Year) = N(Year,0);
    Fstd(Year) = F(Year); 
  }
  AddvarSTD = sqrt(Addvar);
}

void model_parameters::Do_Projections(void)
{
  int Year, Age, II, ii, noT;
  dvariable Fish,Temp,Tmp1,error, aveRYT;
  dvariable FutBias, TheCV2, FCHF, FCHW;
    ofstream out0("main.out");//,ios::app);
    ofstream out1("Depletion.out");//,ios::app);
    ofstream out2("Dep1990.out");//,ios::app);
    ofstream out4("fingerplots.out");//,ios::app);
    ofstream out5("CPUE.out");//,ios::app);
    ofstream out6("Economic.out");//, ios::app);
    ofstream out7("B_Bmsy.out");//,ios::app);
    ofstream out8("Sp_SPmsy.out");//ios::app);
    ofstream out9("Sp_SPmsy(1990).out");//,ios::app);
  int u=0;
  PresValue=0;
  // *** DO for each year of the projection period ***
  // -------------------------------------------------
  for (Year=last_yr; Year<=last_yr+NProj; Year++)
  {
    aveRYT = 0;  
    {
      aveRYT = aveRYT + RY(Year-1) + RY(Year-2)+ RY(Year-3)+ RY(Year-4)+RY(Year);
      aveRYT = aveRYT/5;
    }
    futCatch(last_yr-1) = 148;
    futCatch(Year) = 0.8*aveRYT;
    futCatch(last_yr) = 150;
    if (futCatch(Year)>1.1*futCatch(Year-1))
    {
      futCatch(Year)=1.1*futCatch(Year-1);
    }     
    if (futCatch(Year)<0.9*futCatch(Year-1))
    {
      futCatch(Year)=0.9*futCatch(Year-1);
    }
    // *** Future selectivity same as last year's ***
    // ----------------------------------------------
    for (Age=0; Age<=plus_grp; Age++)
    {
      S(Year,Age) = S(last_yr,Age);
    }
        // *** Exploitable biomass and fishing mortality ***
    // -------------------------------------------------
    Bexp(Year) = 0;
    for (Age=0; Age<=plus_grp; Age++)
    {
      Bexp(Year) += Wmid(last_yr,Age)*S(Year,Age)*N(Year,Age)*mfexp(-M(Age)/2.0);
    }
    Fish = futCatch(Year)/Bexp(Year);
    // ** Make sure fishing mortality is less than 1 **
    // ------------------------------------------------
    if (Fish>0.95)
    {
       //cout<<"!!! Fishing mortality greater than 0.95 in year "<<Year<<endl;
       Fish = 0.95;
    }
    F(Year) = Fish;
    // *** Remove the catch and compute survey selected biomass ***
    // ------------------------------------------------------------
    SurvMid(Year) = 0;
    for (Age=0; Age<=plus_grp; Age++)
    {
      CAPred(Year,Age) = N(Year,Age)*mfexp(-M(Age)/2.0)*S(Year,Age)*Fish;
      SurvBeg(Year) += SurvS(Age)*Wstrt(last_yr,Age)*Ntemp(Age)*(1.0-S(Year,Age)*Fish/2.0);    
      Ntemp(Age) = N(Year,Age)*mfexp(-M(Age)/2.0) - CAPred(Year,Age);
    }
    // *** Adjust Bexp to the middle of the year ***
    // ---------------------------------------------
    Bexp(Year) *= (1.0 - Fish/2.0);
    // *** Remove the rest of natural mortality ***
    // --------------------------------------------
    for (Age=0; Age<=plus_grp; Age++)
      Ntemp(Age) = Ntemp(Age)*mfexp(-M(Age)/2.0);
    // *** Update the age-structure ***
    // --------------------------------
    N(Year+1,plus_grp) = Ntemp(plus_grp)+Ntemp(plus_grp-1);
    for (Age=1; Age<plus_grp; Age++)
      N(Year+1,Age) = Ntemp(Age-1);
    // *** Find the spawner biomass ***
    // --------------------------------
    Spawn(Year+1) = 0;
    TotalB(Year+1) =0;
    for (Age=1; Age<=plus_grp; Age++)
    {
      TotalB(Year+1) += Wstrt(last_yr,Age)*N(Year+1,Age);
      Spawn(Year+1) += mat(Age)* Wstrt(last_yr,Age)*N(Year+1,Age);  
    }
    // *** Generate the recruitment ***
    // --------------------------------
    N(Year+1,0) = ((Spawn(Year+1)*Alpha(last_yr))/(Beta(last_yr) + Spawn(Year+1)))*mfexp(Recruit_Res(Year+1));
    //*****Calculate Economic stuff*****
    //-------------------------------------------------------------------------------------------------
    TAC(Year)=futCatch(Year)*1000;                   //projected TAC
    CPUEfut(Year)=(Bexp(Year)*qCPU(3))/1000;   //projected CPUE
    //cout<<Year<<" "<<CPUEfut(Year)<< "   "<<TAC(Year)<<endl;
    TACF(Year)=TAC(Year)*propF;           //TAC allocated to Freezer
    TACW(Year)=TAC(Year)*propW;           //TAC allocated to wetfish
    TempF(Year)=TACF(Year)/(CPUEfut(Year)*1.4);                     //number of hours allocated to freezers (1.29*CPUE)
    TempW(Year)=TACW(Year)/(CPUEfut(Year)*1);                      //number of hours allocated to wetfish  ((0.71 *CPUE)
    //cout<<Year<<" "<<TempW(Year)<<" "<<TempF(Year)<<" "<<TACW(Year)<<endl;
   //****number of vessels needed ****
    VesselF(Year) =TempF(Year)/fishTimF/TrwdayF; 
    VesselW(Year) =TempW(Year)/fishTimW/TrwdayW;
     //cout<<Year<<" "<<VesselF(Year)<<" "<<VesselW(Year)<<" "<<fishTimF<<" "<<fishTimW<<" "<<TrwdayW<<endl;
   //****Cost, revenue and profit***********
    FCHF=AocF/fishTimF/TrwdayF*1000000;   ///Fishing cost per hour for freezer
    FCHW=AocW/fishTimW/TrwdayF*1000000;   //fihing cost per hour for wetfish
    CostF(Year)=((TempF(Year)*FCHF)+(TACF(Year)*AocFFF))/1000000;
    CostW(Year)= ((TempW(Year)*FCHW)+(TACW(Year)*AocFFW))/1000000;
    CostT(Year)=CostF(Year) + CostW(Year);
     //cout<<AocF<<" "<<fishTimF<< " "<<TrwdayF<<" "<<FCHF<<endl;
    RevF(Year)=TACF(Year)*ApF/1000000;
    RevW(Year)=TACW(Year)*ApW/1000000;
    RevT(Year)=RevF(Year)+RevW(Year);
    ProfitF(Year) =RevF(Year)-CostF(Year);
    ProfitW(Year) =RevW(Year)-CostW(Year);
    ProfitT(Year) =ProfitF(Year)+ ProfitW(Year);
   //*****employment figures***********
   //cout<<Year<<" "<<RevF(Year)<<"  "<<RevW(Year)<<" "<<RevT(Year)<<endl;
    EmployF(Year) =VesselF(Year)*EmpF;
    EmployW(Year) =VesselW(Year)*EmpW;
    EmployFFF(Year)= TACF(Year)/1000*EmpFFF;
    EmployFFW(Year)=TACW(Year)/1000*EmpFFW;
    EmployT(Year) = EmployF(Year)+EmployW(Year)+EmployFFF(Year)+EmployFFW(Year);
    //cout<<Year<<" "<<EmployF(Year)<<" "<<EmployW(Year)<<" "<<EmployFF(Year)<<" "<<EmpF<<"  "<<EmpW<<endl;
    //*****calculate Present value******
    //cout<< Year<<" "<<ProfitT(Year)<<endl;
    ProfitT(Year)=ProfitT(Year)*pow((1 + Pri),u);
    //cout<<Year<<" "<<ProfitT(Year)<<endl;
    PV(Year)=PV(Year)+ProfitT(Year)/pow((1+dr),u);
    u=u+1;
     //cout<<Year<<" "<<RevT(Year)<<" "<<CostT(Year)<<" "<<ProfitT(Year)<< " "<<PV(Year)<<endl;
     PresValue = PresValue + PV(Year);
     // *** Save the catch and spawning biomass ***
    // -------------------------------------------
    out0<<
		Year<<" ,Catch ,"         <<futCatch(Year)<<endl<<
		Year<<" ,Depletion,"      <<Spawn(Year)/mfexp(par_B0)<<endl<<
    Year<<" ,Depletion_1990," <<TotalB(Year)/TotalB(1990)<<endl<<
    Year<<" ,TAC, "          <<TAC(Year)<<endl<<
		Year<<" ,PresVal ,"       <<PV(Year)<<endl<<
		Year<<" ,EmployT,"        <<EmployT(Year)<<endl<<
		Year<<" ,ProfitT,"        <<ProfitT(Year)<<endl<<
		Year<<" ,VesselF,"        <<VesselF(Year)<<endl<<
		Year<<" ,VesselW,"        <<VesselW(Year)<<" "<<endl;
    out1<<Year<<" "<<futCatch(Year)<<" "<<Spawn(Year)/mfexp(par_B0)<<endl;
    out2<<Year<<" "<<futCatch(Year)<<" "<<TotalB(Year)/TotalB(1990)<<endl;
    out4<<Year<<"  "<<futCatch(Year)<<" "<<TotalB(Year)/TotalB(1990)<<" "<<Bexp(Year)<<endl;
    out5<<Year<<" "<<futCatch(Year)<<" "<<(Bexp(Year)*qCPU(3))<< endl;
    out6<<Year<<" "<<TAC(Year)<<"  "<<PV(Year)<<" "<<EmployT(Year)<<" " <<ProfitT(Year)<<" "<<VesselF(Year)<<" "<<VesselW(Year)<<" "<<endl;
    out7<<Year<<" "<<futCatch(Year)<<" "<<TotalB(Year)<<" "<<Bmsy<<" "<<TotalB(Year)/Bmsy<<endl;
    out8<<Year<<" "<<futCatch(Year)<<" "<<Spawn(Year)<<" "<<Spmsy<<" "<<Spawn(Year)/Spmsy<<endl;
    out9<<Year<<" "<<futCatch(Year)<<" "<<Spawn(1990)<<" "<<Spmsy<<" "<<Spawn(Year)/Spawn(1990)<<endl;
  }
  out1<<" "<<endl;
  out2<<" "<<endl;
  out4<<" "<<endl;
  out5<<" "<<endl;
  out6<<" "<<endl;
  out7<<" "<<endl;
  out8<<" "<<endl;
  out9<<" "<<endl;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  save_gradients(gradients);
  int Iser, Year,Age,i,no, p;
  dvariable aveRY, n;
  countAll =count1 + count2 +count3 + count4 + count5;
  aveRY = 0;
  no = 0; 
  for (Year=1990;Year<=last_yr;Year++)
  {
    aveRY = aveRY + RY(Year);
    no = no + 1;
  }
  aveRY = aveRY/no;
  //report <<"aveRy "<<aveRY<<endl;
  report <<"Model: "<<model_number<<" "<<model_name<<" Weight= "<<weight<<" select= "<<select<<endl;
  report <<"Datafile_used: "<<datafilename<<endl;
  report <<"SigmaR_input SigmaR_output MSY Bsp/K"<<endl;
  report <<SigmaRec<<" "<<SigR_out<<" "<<MSY<<" "<<Spawn(last_yr)/mfexp(par_B0)<<endl;
  report <<" "<<endl; 
  report <<"Alpha_first "<<Alpha(first_yr)<<" Beta_first "<<Beta(first_yr)<<" Alpha_last "<<Alpha(last_yr)<<" Beta_last "<<Beta(last_yr)<<" h "<<Steep<<endl;
  // ==+==+==+== Summary statistics ==+==+==+==
  report <<"'ln(Likelihoods) "<<endl;
  report << "overall      "<<obj_fun<< endl;
  //cout << "overall      "<<CPUE_Like+Survey_Like+CAA_Likelihood+CAAS_Likelihood+RecRes_Likelihood<< endl;
  report << "CPUE         "<<CPUE_Like << endl;
  report << "Survey       "<<Survey_Like << endl;
  report << "CAA          "<<CAA_Likelihood << endl;
  report << "CAAsurv      "<<CAAS_Likelihood << endl;
  report << "SRresidual   "<<RecRes_Likelihood << endl;
  report << "negpen       "<<negpen            << endl;
  report << "selpen       "<<sel_pen            << endl;
  report << "prior        "<<prior             << endl;
  report << "One_year_old's_biomass "<<Oneyearold_Likelihood <<endl;
  report << "Num_parameters_Estimated "<<initial_params::nvarcalc()<<endl;
  p=initial_params::nvarcalc();
  n=countAll;
  Akaike=(2*obj_fun+2*p*(n/(n-p-1)));
  report<<"Akaike_info_crit "<<Akaike<<endl;
  report << " " << endl;
  report << "Ksp                "<<mfexp(par_B0)<<endl;
  report << "Kexp               "<<Bexp(first_yr)<<endl;
  report << "TotalB("<<last_yr<<")            "<<TotalB(last_yr)<<endl;
  report << "Bsp("<<last_yr<<")            "<<Spawn(last_yr)<<endl;
  //cout << "Bsp("<<last_yr<<")            "<<Spawn(last_yr)<<endl;
  report << "Bexp("<<last_yr<<")           "<<Bexp(last_yr)<<endl;
  report << "h                  "<<Steep<<endl;
  //cout << "h                  "<<Steep<<endl;
  report << "TotalB(MSY)           "<<Bmsy<<endl;
  report << "Bsp(MSY)           "<<Spmsy<<endl;
  report << "Bexp(MSY)          "<<Bexpmsy<<endl;
  report << "MSY                "<<MSY<<endl;
  report << "Bsp("<<last_yr<<")/Ksp        "<<Spawn(last_yr)/mfexp(par_B0)<<endl;
  report << "Bexp("<<last_yr<<")/Kexp      "<<Bexp(last_yr)/Bexp(first_yr)<<endl;
  report << "         "<<endl;
  report << "TotalB(last_yr)/TotalB(1990)      "<<TotalB(last_yr)/TotalB(1990)<<endl;
  // cout<<"Total B "<<TotalB(last_yr)<<endl;
  report << "Bsp("<<last_yr<<")/Bsp(msy)   "<<Spawn(last_yr)/Spmsy<<endl;
  report<<"aveRY                "<<aveRY<<endl;
  report<<"Catchlastyr-1/aveRY   "<<Catch(last_yr-1)/aveRY<<endl;
  report << "         "<<endl;
  report << "Bexp("<<last_yr<<")/Bexp(msy) "<<Bexp(last_yr)/Bexpmsy<<endl;
  report << "MSYL/Ksp           "<<Spmsy/mfexp(par_B0)<<endl;
  //cout << "MSYL/Ksp           "<<Spmsy/mfexp(par_B0)<<endl;
  report << "MSYL/Kexp          "<<Bexpmsy/Bexp(first_yr)<<endl;
  // cout << "M                  "<<M(5)<<endl;
  report << "M                  "<<M(5)<<endl;
  // ---------- Selectivity at age ----------
  report<<"Age S(1968) S(1989) S("<<last_yr<<") survS"<<endl;
  for (Age=0; Age<=plus_grp; Age++)
   report <<Age<<" "<<S(1968,Age)<<" "<<S(1989,Age)<<" "<<S(last_yr,Age)<<" "<<SurvS(Age)<<endl;
  report << " "<<endl;
  //
  // ---------- q's and sigma's ----------
  report << "sigmaCPUE1      " << SigCPU(1) << endl;
  report << "sigmaCPUE2      " << SigCPU(2) << endl;
  report << "sigmaCPUE_GLM   " << SigCPU(3) << endl;
  report << "sigmaCPUE_7Vessel "<<SigCPU(6) <<endl;
  report << " " << endl;
  report << "qCPUE1          " << qCPU(1) << endl;
  report << "qCPUE2          " << qCPU(2) << endl;
  report << "qCPUE_GLM       " << qCPU(3) << endl;
  report << "qCPUE_7Vessel   " << qCPU(6) << endl;
  report << "qSeal_index      "<<qSeal<< endl;
  report << " " << endl;  
  report << "sigmaSpanWint   " << SigCPU(4) << endl;
  report << "sigmaSpanSum    " << SigCPU(5) << endl;
  report << " " << endl; 
  report << "SpanSurveyqWin  " << qCPU(4) << endl;
  report << "SpanSurveyqSum  " << qCPU(5) << endl;
  report << "Addvariance    "<<Addvar<<endl;
  report << " "<<endl;
  report << "NanSurveyqSum   " << qSurvPre(1) <<" "<< qSurvPost(1) << " "<<mean(qSurvPre(1)*SurvS(3,6))<<" "<<mean(qSurvPost(1)*SurvS(3,6))<<endl;
  report << "NanSurveyqWin   " << qSurvPre(2) <<" "<< qSurvPost(2) << " "<<mean(qSurvPre(2)*SurvS(3,6))<<" "<<mean(qSurvPost(2)*SurvS(3,6))<<endl;
  // cout << "NanSurveyqSum   " << qSurvPre(1) <<" "<< qSurvPost(1) << endl;
  // cout << "NanSurveyqWin   " << qSurvPre(2) <<" "<< qSurvPost(2) <<endl;
  report << " " << endl;
  report << "sigCAA_com      " << SigCAA_com << endl;
  report << "sigCAA_summer   " << SigCAA_surv(1) << endl;
  report << "sigCAA_winter   " << SigCAA_surv(2) << endl;
  report << "sigSeal_Index   " << SigSeal_Index << endl;
  report << "Added_Sigma_Seal    "<<sig_seal_added<<" prop "<< Prop <<endl;
  report << " "<<endl;
  report << " " << endl;
  //
  // ---------- Natural Mortality ----------
  report<<"Natural mortality by age "<<endl;
  for (Age=0; Age<=plus_grp; Age++)
   report << "M("<<Age<<")       "<<M(Age)<<endl;
  report <<"Minf  "<<Minf<<"  Minfage "<< Minfage <<endl;
  //
  //
  // ==+==+==+== Fit to abundance data ==+==+==+==
  report<<" *** Fit to abundance data *** "<<endl;
  for (Iser=1; Iser<=NCpueSeries; Iser++)
    report<<"CPUEseries"<<Iser<<" q "<<qCPU(Iser)<<" Sigma "<<SigCPU(Iser)<<endl;
  for (Iser=1; Iser<=NSurveySeries; Iser++)
    report<<"Survseries"<<Iser<<" qPred "<<qSurvPre(Iser)<<" qPost "<<qSurvPost(Iser)<<endl;
  report<<" "<<endl;
  report<<"Year obs1 pred1 obs2 pred2 obs3 pred3 obs4 pred4 obs5 pred5 obs6 pred6 NanObs1 NanPred1 NanObs2 NanPred2 Sealobs Sealpred"<<endl;
  for (Year=first_yr; Year<=last_yr; Year++)
  {
    report<<Year<<" ";
    // *** CPUE fit ***
    // ----------------
    for (Iser=1; Iser<=NCpueSeries; Iser++)
    {
      if (CpueIndx(Iser) == 1) BIO(Year) = Bexp(Year);
      if (CpueIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
      if (CpueIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
        report<<CPUE(Year,Iser)<<" "<<qCPU(Iser)*BIO(Year)<<" ";
    }
    // *** Survey fit ***
    // ------------------
    for (Iser=1; Iser<=NSurveySeries; Iser++)
    {
      if (SurveyIndx(Iser) == 1) BIO(Year) = Bexp(Year);
      if (SurveyIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
      if (SurveyIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
      if (Year<=q_chngYr(Iser))
        report<<Survey(Year,Iser*2-1)<<" "<<qSurvPre(Iser)*BIO(Year)<<" ";
      else
        report<<Survey(Year,Iser*2-1)<<" "<<qSurvPost(Iser)*BIO(Year)<<" ";
     }
     report<<seali(Year)/qSeal<<" "<<seali_pred(Year)<<" "<<N(Year,0)<<" ";
     report<<" "<<endl;
   }
  report << " " << endl;
  //
  // 
  // ==+==+==+== Some more values ==+==+==+==
  report << "Virgin Age Structure" << endl << Nvirg << endl;
  report << " " << endl;
  report << "Numbers at Age" <<endl;
  for (Year=first_yr;Year<=last_yr;Year++) 
  {
    report<<Year<<" ";
    for (Age=0; Age<=plus_grp; Age++)
      report<<N(Year,Age)<<" ";
    report<<" "<<endl;
  }
  report << " " << endl;
  report << "Year Bsp Bexp F RecRes Catch Depletion TotalB"<<endl;
  for (Year=first_yr; Year<=last_yr; Year++)
  {
    report <<Year<<" "<<Spawn(Year)<<" "<<Bexp(Year)<<" "<<F(Year) <<" "<<Recruit_Res(Year)<<" "<<Catch(Year)<<"  "<<Spawn(Year)/Spawn(first_yr)<<"  "<<TotalB(Year)<<endl;
  }
  report<<" "<<endl;
  //
  // 
  report<<"Commercial selectivity"<<endl;
  for (Year=first_yr;Year<=last_yr;Year++) 
  {
    report<<Year<<" ";
    for (Age=0; Age<=plus_grp; Age++)
      report<<S(Year,Age)<<" ";
    report<<" "<<endl;
  }
  report<<" "<<endl;
  for (Year=first_yr;Year<=last_yr;Year++)
  { 
     if (sum(CAA(Year)(0,8))>0)
     {
        report<<Year<<"  "<<CAA(Year)(0,8) <<" ";
        for (i=0;i<CAAMinus;i++)
          report<< " 0  ";
        report<<CAPred(Year)(CAAMinus,CAAPlus)<<" ";
        for (i=CAAPlus+1;i<=plus_grp;i++)
          report<< " 0  ";
        report<<endl;
     }
  }
  //report<<" "<<endl;
  for (Iser=1;Iser<=NSurveySeries;Iser++)
  {
    report<<"'----------Survey"<<Iser<<"'----------"<<endl;
    for (Year=first_yr;Year<=last_yr;Year++)
    {
      if (sum(SurvCAA(Iser, Year)(0,CAASPlus(Iser)))>0)
      {
        report<<Year<<"  "<<SurvCAA(Iser,Year)(0,8)<<" " ;
        for (i=0;i<CAASMinus(Iser);i++)
          report<< " 0  ";
        report<<SAPred(Iser,Year)(CAASMinus(Iser),CAASPlus(Iser))<<" ";
        for (i=CAASPlus(Iser)+1;i<=plus_grp;i++)
          report<< " 0  ";
        report<<endl;
          // report<<Year<<" "<<SurvCAA(Iser,Year)(0,8)<<" "<<SAPred(Iser,Year)(0,8)             <<endl;
      }     
    }
    //report<<" "<<endl;
  }  
  ofstream RYout("ReplacYield.dat");
   RYout<<"Replacement Yield by year"<<endl;
   RYout<<"MSY="<<"  "<<MSY<<endl;
   RYout<<"Year ReplacYield Bsp Catch"<<endl;
   for (Year=first_yr;Year<=last_yr;Year++)
     RYout<<Year<<" "<<RY(Year)<<" "<<Spawn(Year)<<endl;
   RYout.close();
  report<<"Replacement Yield by year"<<endl;
  report<<"Year ReplacYield Bsp Catch"<<endl; 
  for (Year=first_yr;Year<=last_yr;Year++)
     report<<Year<<" "<<RY(Year)<<" "<<Spawn(Year)<<" "<<Catch(Year)<<endl;
  report<<"'----------Stock-recruitment_curve----"<<endl;
  report << "0 "<<SSB<<endl;
  report << "0 "<<Pred_Rec<<endl;
  report<<"'----------Catch-at-Age----------"<<endl;
  report<<"'----------CommercialCAA----------"<<endl;
  for (Year=first_yr;Year<=last_yr;Year++)
  { 
    for (Age=CAAMinus; Age<=CAAPlus; Age++)
      if (CAA(Year,Age)>0)
        if (CAPred(Year,Age)>0)
       {
        report<<Year<<"  "<<Age<<" "<<CAA(Year,Age)<<" "<<CAPred(Year,Age)<<" "<<(log(CAA(Year,Age))-log(CAPred(Year,Age)))/(SigCAA_com/sqrt(CAPred(Year,Age)))<<" "<<endl;
        //cout<<Year<<"  "<<Age<<" "<<CAA(Year,Age)<<" "<<CAPred(Year,Age)<<" "<<(log(CAA(Year,Age))-log(CAPred(Year,Age)))/(SigCAA_com/sqrt(CAPred(Year,Age)))<<" "<<endl;
      }
  }
  //report<<" "<<endl;
  for (Iser=1;Iser<=NSurveySeries;Iser++)
  {
    report<<"'----------Survey"<<Iser<<"'----------"<<endl;
    for (Year=first_yr;Year<=last_yr;Year++)
    {
      int s = CAASMinus(Iser);
      int t = CAASPlus(Iser);
      for (Age=s; Age<=t; Age++) 
      if (SurvCAA(Iser,Year,Age)>0)
        if (SAPred(Iser,Year,Age)>0)
        {
          report<<Year<<" "<<Age<<" "<<SurvCAA(Iser,Year,Age)<<" "<<SAPred(Iser,Year,Age)<<" "<<(log(SurvCAA(Iser,Year,Age))-log(SAPred(Iser,Year,Age)))/(SigCAA_surv(Iser)/sqrt(SAPred(Iser,Year,Age)))<<" "<<endl;
          //cout<<Year<<" "<<Age<<" "<<SurvCAA(Iser,Year,Age)<<endl;
        }     
     }
    //report<<" "<<endl;
  }  
  report<<" "<<endl;
  report<<" "<<endl;
  report<<" biomass of fish by age "<<endl;
   for (Year=first_yr+1;Year<=last_yr;Year++)
   {
    report<<Year<<" ";
    for (Age=0; Age<=plus_grp; Age++)
      report<<N(Year,Age)*Wstrt(Year,Age)<<" ";
    report<<" "<<endl;
   }
  report<<" biomass of fish that died due to natural mortality by age for J-P Roux "<<endl;
  for (Year=first_yr+1;Year<=last_yr;Year++)
   {
    report<<Year<<" "<<TotalB(Year)<<" ";
    for (Age=1; Age<=plus_grp; Age++)
      report<<Mbio(Year,Age)<<" ";
    report<<" "<<endl;
   }
  report<<"   " <<endl;
  report<<"----- The End -----"<<endl;
}

void model_parameters::final_calcs()
{
  ofstream out4("stats.out",ios::app);
  Do_Projections();
  cout<<"Present Value "<<PresValue<<endl;
  //cout<<" f "<<f<<endl;
  out4<<model_number<<" "<<model_name<<" "<<TotalB(last_yr)/TotalB(first_yr)<<" "<<TotalB(last_yr)/TotalB(1990)<<"  " <<Akaike<<"   "<<Spawn(last_yr)/Spmsy<<"   "<<PresValue<<endl;
  //"model_number","model_name","TotalB_depl", "TotalB_1990",  " Akaike",   "SSB/Bmsy",   "PresValue"
  Rreport<<model_number<<" "<<model_name<<" "<<TotalB(last_yr)/TotalB(first_yr)<<" "<<TotalB(last_yr)/TotalB(1990)<<"  " <<Akaike<<"   "<<Spawn(last_yr)/Spmsy<<"   "<<PresValue<<endl;
  Rreport << "ObjFun"<<endl<<obj_fun<< endl;
  R_report(CPUE_Like);
  R_report(Survey_Like);
  R_report(CAA_Likelihood);
  R_report(CAAS_Likelihood);
  R_report(RecRes_Likelihood);
  R_report(Oneyearold_Likelihood);
  Rreport << "Npars"<<endl<<initial_params::nvarcalc()<<endl;
  R_report(Akaike);
  R_report(KspSTD);
  R_report(KexpSTD);
  R_report(TotalB);
  Rreport << "Bsp"<<endl<<Spawn<<endl;
  R_report(Bexp);
  R_report(Steep);
  R_report(Bmsy);
  R_report(Spmsy);
  R_report(Bexpmsy);
  R_report(MSY);
  R_report(DepspSTD);
  R_report(DepexSTD);
	Rreport <<"Year"<<endl;
  for (int Year=first_yr; Year<=last_yr; Year++)
	  Rreport <<Year<< " ";
	Rreport<<endl;
 // ==+==+==+== Fit to abundance data ==+==+==+==
  for (int Iser=1; Iser<=NCpueSeries; Iser++){
    Rreport<<"Obs_CPUE_"<<Iser<<endl;
    for (int Year=first_yr; Year<=last_yr; Year++)
      Rreport<<CPUE(Year,Iser)<<" ";
    Rreport<<endl<<"e_CPUE_"<<Iser<<endl;
    for (int Year=first_yr; Year<=last_yr; Year++){
      if (CpueIndx(Iser) == 1) BIO(Year) = Bexp(Year);
      if (CpueIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
      if (CpueIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
			Rreport<<qCPU(Iser)*BIO(Year)<<" ";
		}
		Rreport<<endl;
	  }
  for (int Iser=1; Iser<=NSurveySeries; Iser++)
  {
    Rreport<<"Obs_Survey_"<<Iser<<endl;
    for (int Year=first_yr; Year<=last_yr; Year++)
      Rreport<<Survey(Year,Iser*2-1)<<" ";
    Rreport<<endl<<"Pre_Survey_"<<Iser<<endl;
    for (int Year=first_yr; Year<=last_yr; Year++)
		{
      if (SurveyIndx(Iser) == 1) BIO(Year) = Bexp(Year);
      if (SurveyIndx(Iser) == 2) BIO(Year) = SurvMid(Year);
      if (SurveyIndx(Iser) == 3) BIO(Year) = SurvBeg(Year);
      if (Year<=q_chngYr(Iser))
			  Rreport<<qSurvPre(Iser)*BIO(Year)<<" ";
      else
			  Rreport<<qSurvPost(Iser)*BIO(Year)<<" ";
     }
     Rreport<<" "<<endl;
   }
   R_report(N);
   R_report(S);
   R_report(CAA);
   R_report(CAPred);
  // For plotting from SR curve
 	R_report(SSB);
 	R_report(Pred_Rec);
 	R_report(Pred_Rec.sd);
  R_report(Bstd);
  R_report(Bstd.sd);
 //  R_report(Rstd.sd);
  Rreport << "Rstd_sd"<<endl<<Rstd.sd<< endl;
  R_report(Rstd);
  R_report(Steep);
	csvrep << "Variable, Year, value, ymin, ymax" <<endl;
  csvrep << "Cur_B0"     <<",NA,"<< Cur_B0      <<","<< Cur_B0-1.96*Cur_B0.sd <<","<< Cur_B0+1.96*Cur_B0.sd <<endl;
  csvrep << "Cur_Bmsy"   <<",NA,"<< Cur_Bmsy    <<","<< Cur_Bmsy-1.96*Cur_Bmsy.sd <<","<< Cur_Bmsy+1.96*Cur_Bmsy.sd <<endl;
  csvrep << "Cur_90"     <<",NA,"<< Cur_90      <<","<< Cur_90-1.96*Cur_90.sd <<","<< Cur_90+1.96*Cur_90.sd <<endl;
  csvrep << "aveRY_90"   <<",NA,"<< aveRY_90    <<","<< aveRY_90-1.96*aveRY_90.sd <<","<< aveRY_90+1.96*aveRY_90.sd <<endl;
  csvrep << "aveRY_last5"<<",NA,"<< aveRY_last5 <<","<< aveRY_last5-1.96*aveRY_last5.sd <<","<< aveRY_last5+1.96*aveRY_last5.sd <<endl;
  double lbfy=value(Bstd(first_yr)/exp(2.*sqrt(log(1+square(Bstd.sd(first_yr))/square(Bstd(first_yr))))));
  double ubfy=value(Bstd(first_yr)*exp(2.*sqrt(log(1+square(Bstd.sd(first_yr))/square(Bstd(first_yr))))));
  for (int Year=first_yr;Year<=last_yr;Year++){
    double lb=value(Rstd(Year)/exp(2.*sqrt(log(1+square(Rstd.sd(Year))/square(Rstd(Year))))));
    double ub=value(Rstd(Year)*exp(2.*sqrt(log(1+square(Rstd.sd(Year))/square(Rstd(Year))))));
    csvrep << "R"  <<","<<Year<<","<< Rstd(Year) <<","<< lb <<","<< ub <<endl;
    lb=value(Bstd(Year)/exp(2.*sqrt(log(1+square(Bstd.sd(Year))/square(Bstd(Year))))));
    ub=value(Bstd(Year)*exp(2.*sqrt(log(1+square(Bstd.sd(Year))/square(Bstd(Year))))));
    csvrep << "SSB"<<","<<Year<<","<< Bstd(Year) <<","<< lb <<","<< ub <<endl;
    lb=value(Bstd(Year)/exp(2.*sqrt(log(1+square(Bstd.sd(Year))/square(Bstd(Year))))));
    ub=value(Bstd(Year)*exp(2.*sqrt(log(1+square(Bstd.sd(Year))/square(Bstd(Year))))));
    csvrep << "Depletion"<<","<<Year<<","<< Bstd(Year)/Bstd(first_yr) <<","<< lb/lbfy <<","<< ub/ubfy <<endl;
    csvrep << "B_Bmsy"<<","<<Year<<","<< Spawn(Year)/Spmsy  <<","<< "NA" <<","<< "NA" <<endl;
    csvrep << "Catch"<<","<<Year<<","<< Catch(Year) <<","<< "NA"<<","<< "NA"<<endl;
    csvrep << "RY"<<","<<Year<<","<< RY(Year) <<","<< "NA"<<","<< "NA"<<endl;
    csvrep << "Catch_RY"<<","<<Year<<","<< Catch(Year)/RY(Year) <<","<< "NA"<<","<< "NA"<<endl;
	}
  R_report(Cur_B0);
  R_report(Cur_B0.sd);
  R_report(Cur_Bmsy);
  R_report(Cur_Bmsy.sd);
  R_report(Cur_90);
  R_report(Cur_90.sd);
  R_report(aveRY_90);
  R_report(aveRY_90.sd);
  R_report(aveRY_last5);
  R_report(aveRY_last5.sd);
  R_report(Steep);
  R_report(prior);
  R_report(negpen);
	for (int Iser=1;Iser<=NSurveySeries;Iser++) {
	  Rreport <<"survey"<<Iser<<"_Phat"<<endl;
    for (int Year=first_yr;Year<=last_yr;Year++) {
      int s = CAASMinus(Iser);
      int t = CAASPlus(Iser);
      Rreport<< Year <<" ";
      for (int Age=s; Age<=t; Age++) 
      {
        Rreport<< SAPred(Iser,Year,Age)<<" ";
      }     
			Rreport<<endl;
    }     
	  Rreport <<"survey"<<Iser<<"_Pobs"<<endl;
    for (int Year=first_yr;Year<=last_yr;Year++) {
      int s = CAASMinus(Iser);
      int t = CAASPlus(Iser);
      Rreport<< Year <<" ";
      for (int Age=s; Age<=t; Age++) 
      {
        Rreport<< SurvCAA(Iser,Year,Age)<<" ";
      }     
			Rreport<<endl;
    }     
	  Rreport <<"fishery_Pobs"<<endl;
    for (int Year=first_yr;Year<=last_yr;Year++) {
      int s = CAAMinus;
      int t = CAAPlus;
      Rreport<< Year <<" ";
      for (int Age=s; Age<=t; Age++) 
      {
        Rreport<< CAA(Year,Age)<<" ";
      }     
			Rreport<<endl;
    }     
	  Rreport <<"fishery_Phat"<<endl;
    for (int Year=first_yr;Year<=last_yr;Year++) {
      int s = CAAMinus;
      int t = CAAPlus;
      Rreport<< Year <<" ";
      for (int Age=s; Age<=t; Age++) 
      {
        Rreport<< CAPred(Year,Age)<<" ";
      }     
			Rreport<<endl;
    }     
  }     
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(200000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  gradient_structure::set_MAX_NVAR_OFFSET(500);
  arrmblsize = 300000;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint = defaults::iprint;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
