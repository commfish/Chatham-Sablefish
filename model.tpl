// ------------------------------------------------------------------------------ //
//               Age-structured model for Chatham Strait Sablefish                //
//                                                                                //
//                                  VERSION 0.2                                   //
//                                   Dec  2015                                    //
//                                                                                //
//                                   AUTHORS                                      //
//                                Kray Van Kirk                                   //
//                            kray.vankirk@alaska.gov                             //
//                                                                                //
//                   Built on a Pacific Ocean Perch template                      //
//                                Dana Hanselman                                  //
//                            dana.hanselman@noaa.gov                             //
//                                                                                //
//                    Scripts in FINAL and GLOBAL sections                        //
//                                Steven Martell                                  //
//                           martell.steve@gmail.com                              //
//                                                                                //
 
// ------------------------------------------------------------------------------ //
//-- CHANGE LOG:                                                                --//
//--  Dec 2016 -                                                                --//
//--               :Hessian-derived estimates of variance from binomial         --//
//--                 likelihoods for both survey and fishery MR                 --//
//--               :Survey MR included as distinct recapture event              --//
//--               :Estimated MR abundance conditioned on selectivity vectors   --//
//--  Dec 2015 -                                                                --//
//--               :penalized log-likelihoods (p.l.l.) and penalties            --//
//--               :reduction in penalties and constraints due to p.l.l.        --//
//--               :normal LL for cpue and MR                                   --//
//--               :input M and selectivities (NMFS)                            --//
//--                                                                            --//
//--  Feb 2015 -                                                                --//
//--               :sex-specific                                                --//
//--               :separate fits to 2 survey CPUE series                       --//
//--                 (1988-1996, 1-hr soak; 1997 - present, 3-hr soak)          --//
//--               :sex-specific natural mortality estimated                    --//
//--               :modification for potential code distribution                --//
//--               :sex-specific gear selectivity   (survey and fishery)        --//
//--               :sex-specific full-recruitment fishing mortality F           --//
//--               :year 1 abundance = Year 1 Age 1 decremented by -M           --//
//--                                   for each subsequent age                  --//   
//--                                                                            --//
// ------------------------------------------------------------------------------ //


DATA_SECTION  

  // |--------------------------------------------------------------------------|
  // |MCMC OUTPUT FILE
  // |--------------------------------------------------------------------------|

     !!CLASS ofstream evalout("evalout.prj");
     !!time(&start);
     

  // |--------------------------------------------------------------------------|
  // | STRINGS FOR INPUT FILES                                                  |
  // |--------------------------------------------------------------------------|
  // |
  // | DataFile               : data to condition the assessment model    
  // | ControlFile            : controls for years, phases, and block options 

     init_adstring DataFile;      
     init_adstring ControlFile;    

  // | BaseFileName           : file prefix used for all  model output
  // | ReportFileName         : file name to which report file is printed

     !! BaseFileName = stripExtension(DataFile);  
     !! ReportFileName = BaseFileName + adstring(".rep");

     !! cout<<"You are modeling the "<<BaseFileName<<" stock of sablefish"<<endl;
     !! cout<<""<<endl;


  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM CONTROL FILE                                             |
  // |--------------------------------------------------------------------------|
  // | DEBUG_FLAG             : Boolean Flag used for manual debugging
  // | nboot                  : number of parametric bootstraps to run

     !! ad_comm::change_datafile_name(ControlFile);

     init_int DEBUG_FLAG;
     init_int nboot;


 // |-------------------------------------------------------------------------|
 // | DESIGN MATRIX FOR PARAMETER CONTROLS                                    |
 // |-------------------------------------------------------------------------|
 // | - theta_DM -> theta is a vector of estimated parameters.

    init_int n_theta;
    init_matrix theta_DM(1,n_theta,1,7);
    vector    theta_ival(1,n_theta);
    vector      theta_lb(1,n_theta);
    vector      theta_ub(1,n_theta);
    ivector    theta_phz(1,n_theta);
    ivector theta_iprior(1,n_theta);
    vector      theta_p1(1,n_theta);
    vector      theta_p2(1,n_theta);
    !! theta_ival = column(theta_DM,1);
    !! theta_lb  = column(theta_DM,2);
    !! theta_ub  = column(theta_DM,3);
    !! theta_phz = ivector(column(theta_DM,4));
    !! theta_iprior = ivector(column(theta_DM,5));
    !! theta_p1 = column(theta_DM,6);
    !! theta_p2 = column(theta_DM,7);



 // |---------------------------------------------------------------------------|
// | Miscellaneous Controls
// |---------------------------------------------------------------------------|
// | nMiscCont  » Number of controls to read in.
// | dMiscCont  » Vector of miscellaneous controls,


   init_int nMiscCont;
   init_vector dMiscCont(1,nMiscCont);

   number sigma_catch
   number sigr
   number sig1
   number Mm;
   number Mf;
   
   int     ph_rec;
   int     ph_init;
   int     ph_F;
   int     ph_spr;


    //EOF Marker
    init_int eof_ctl;
    
 LOCAL_CALCS

          sigma_catch = dMiscCont(1);
          sigr        = dMiscCont(2);
          sig1        = dMiscCont(3);
          Mm          = dMiscCont(4);
          Mf          = dMiscCont(5);
          ph_rec      = dMiscCont(6);
          ph_init     = dMiscCont(7);
          ph_F        = dMiscCont(8);
          ph_spr      = dMiscCont(9);


    if(eof_ctl==42) cout << BaseFileName<<".ctl has been read correctly!"<<endl;
    else 
    {    
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|      Red alert! Captain to bridge! The .ctl file is compromised!     |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof_ctl<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .ctl file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }

 END_CALCS


  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM DATA FILE                                                |
  // |--------------------------------------------------------------------------|
  // | This calls the data file

     !! ad_comm::change_datafile_name(DataFile);


  // |--------------------------------------------------------------------------|
  // | MODEL STRUCTURAL DIMENSIONS                                              |
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nages       -> number of ages
  // | 2) styr        -> data start year
  // | 3) endyr       -> data end year
  // | 4) recage      -> age at recruitment into the model
  // |
  // | 7) spawn_fract ->spawning month (set to February [2])
  // | 8) srv_fract   ->survey month (set to August [8])
  // | 9) fsh_fract   ->month when fishery takes place (set to October [10])
  // |    
               
     init_int      nages
     init_int      styr
     init_int      endyr
     init_int      recage

     init_number   spawn_fract
     init_number   srv1_fract
     init_number   fshy_fract

     vector        yy(styr,endyr)
     vector        aa(1,nages)


  // |--------------------------------------------------------------------------|
  // | Gear selectivity for fishery and survey (using NOAA values)
  // |--------------------------------------------------------------------------|
  // | Selectivity curves from Dana Hanselman @ NOAA
  
     init_vector	fish_sel_f(1,nages)
     init_vector	fish_sel_m(1,nages)
     init_vector	srv1_sel_f(1,nages)
     init_vector	srv1_sel_m(1,nages)

  
  // |--------------------------------------------------------------------------|
  // | PHYSIOLOGY                                                               
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) p_mature       -> proportion of females mature-at-age
  // | 2) wt_avg_fshy_f  -> mean weight-at-age in commercial fishery, females
  // | 3) wt_avg_fshy_m  -> mean weight-at-age in commercial fishery, males
  // | 4) wt_avg_srv1_f  -> mean weight-at-age in longline survey, females
  // | 5) wt_avg_srv1_m  -> mean weight-at-age in longline survey, males

     init_vector   p_mature(1,nages)
     init_vector   wt_avg_fshy_f(1,nages)
     init_vector   wt_avg_fshy_m(1,nages)
     init_vector   wt_avg_srv1_f(1,nages)
     init_vector   wt_avg_srv1_m(1,nages)


  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE COMMERCIAL CATCH                                   
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) obs_catchf     -> total annual catch, females
  // | 2) obs_catchm     -> total annual catch, males

     init_vector   obs_catchf(styr,endyr)
     init_vector   obs_catchm(styr,endyr)


  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE SURVEY CPUE (1988 - 1996, 1 hr soak)               
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_cpue1      -> number of years of CPUE estimates
  // | 2) yrs_cpue1       -> specific years CPUE estimated
  // | 3) obs_cpue1_biom  -> mean CPUE by year
  // | 4) obs_cpue1_var    -> variance of CPUE

     init_int      nyrs_cpue1
     init_ivector  yrs_cpue1(1,nyrs_cpue1)
     init_vector   obs_cpue1_biom(1,nyrs_cpue1)
     init_vector   obs_cpue1_var(1,nyrs_cpue1)



  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE SURVEY CPUE (1997 +, 3 hr soak)                    
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_cpue2      -> number of years of CPUE estimates
  // | 2) yrs_cpue2       -> specific years CPUE estimated
  // | 3) obs_cpue2_biom  -> mean CPUE by year
  // | 4) obs_cpue2_var    -> variance of CPUE

     init_int      nyrs_cpue2
     init_ivector  yrs_cpue2(1,nyrs_cpue2)    

     init_vector   obs_cpue2_biom(1,nyrs_cpue2)
     init_vector   obs_cpue2_var(1,nyrs_cpue2) 
 

  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE COMMERCIAL CPUE                                    
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_fshy_cpue     -> number of years of CPUE estimates
  // | 2) yrs_fshy_cpue      -> specific years CPUE estimated
  // | 3) obs_fshy_cpue      -> mean CPUE by year
  // | 4) obs_fshy_cpue_var   -> variance of CPUE
  
     init_int      nyrs_fshy_cpue 
     init_ivector  yrs_fshy_cpue(1,nyrs_fshy_cpue)
     init_vector   obs_fshy_cpue(1,nyrs_fshy_cpue) 
     init_vector   obs_fshy_cpue_var(1,nyrs_fshy_cpue)


  // |--------------------------------------------------------------------------|
  // | MARK-RECAPTURE ESTIMATES OF EXPLOITABLE ABUNDANCE                        
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_MR       -> number of years of MR estimates
  // | 2) yrs_MR        -> specific years MR estimated
  // | 3) obs_MR_abund  -> mark-recapture estimates of abundance
  // | 4) obs_MR_var    -> variance of MR estimates

     init_int      nyrs_MR
     init_ivector  yrs_MR(1,nyrs_MR)
     init_vector   obs_MR_abund(1,nyrs_MR) 
     init_vector   obs_MR_var(1,nyrs_MR)
     init_vector   obs_MR_abund_SRV(1,nyrs_MR)
     init_vector   obs_MR_var_SRV(1,nyrs_MR) 


  // |--------------------------------------------------------------------------|
  // | COMMERCIAL LONGLINE FISHERY AGE COMPOSITION                              
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_fish_age       -> number of years of fishery age comps
  // | 2) yrs_fish_age        -> specific years age comps
  // | 3) nsamples_fish_agef  -> some measure of relative sample size FEMALES
  // | 4) oac_fishf           -> female age comp data
  // | 5) nsamples_fish_agem  -> some measure of relative sample size MALES
  // | 6) oac_fishm           -> male age comp data
  // |
  // | Relative sample size = square root of number of fish aged to scale data impact

     init_int      nyrs_fish_age 
     init_vector   yrs_fish_age(1,nyrs_fish_age) 
     init_vector   nsamples_fish_agef(1,nyrs_fish_age)   
     init_matrix   oac_fishf(1,nyrs_fish_age,1,nages) 
     init_vector   nsamples_fish_agem(1,nyrs_fish_age)
     init_matrix   oac_fishm(1,nyrs_fish_age,1,nages)


  // |--------------------------------------------------------------------------|
  // | LONGLINE SURVEY AGE COMPOSITION                                          
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_srv1_age       -> number of years of survey age comps
  // | 2) yrs_srv1_age        -> specific years age comps
  // | 3) nsamples_srv1_agef  -> some measure of relative sample size FEMALES
  // | 4) oac_srv1f           -> female age comp data
  // | 5) nsamples_srv1_agem  -> some measure of relative sample size MALES
  // | 6) oac_srv1m           -> male age comp data
  // |
  // | Relative sample size = square root of number of fish aged to scale data impact

     init_int      nyrs_srv1_age 	
     init_ivector  yrs_srv1_age(1,nyrs_srv1_age)   
     init_vector   nsamples_srv1_agef(1,nyrs_srv1_age) 
     init_matrix   oac_srv1f(1,nyrs_srv1_age,1,nages)  
     init_vector   nsamples_srv1_agem(1,nyrs_srv1_age) 
     init_matrix   oac_srv1m(1,nyrs_srv1_age,1,nages)  
  
  // |--------------------------------------------------------------------------|
  // | AGEING ERROR MATRIX                                                      
  // |--------------------------------------------------------------------------|
  // |
  // | proportion at reader age (columns) given true age (rows)
  // |
  // | This was re-calculated using ADU data and a script from Pete Hulson
  // | kvk - 12/10/2015

     init_matrix   ageage(1,nages,1,nages)
   
  // |  Multinomial "offset", length must be at least as long as number of multinomial likelihoods
  // |  Calculate "offset" for multinomials - fishery age, survey age
  // |  "Offset" value lets the multinomial likelihood equal zero when the observed and predicted are equal
  // |  First step is to ensure that the data are expressed as proportions

     vector offset(1,4);       
  
     number oo;
     !!oo = 0.000001;
 
     int i
     int j

     init_number eof_dat
  

 LOCAL_CALCS

    if(eof_dat==42) cout << BaseFileName<<".dat has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   ** Red alert! Captain to bridge! The .dat file is compromised! **  |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof_dat<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .dat file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }

     yy.fill_seqadd(styr,1);		        //Year vector
     aa.fill_seqadd(recage,1);                  //Age vector
  
     spawn_fract = (spawn_fract - 1) / 12;	// Fraction of year at which spawning occurs
     srv1_fract = (srv1_fract - 1) / 12;	// Fraction of year at which survey occurs
     fshy_fract = (fshy_fract - 1) / 12;	// Fraction of year at which fishery occurs


    offset.initialize();


    for (i=1; i<=nyrs_fish_age; i++)
      {
        oac_fishm(i)/=sum(oac_fishm(i));
        oac_fishf(i)/=sum(oac_fishf(i));
        offset(1) -= nsamples_fish_agem(i)*((oac_fishm(i) + 0.00001)*log(oac_fishm(i) + 0.00001));
        offset(2) -= nsamples_fish_agef(i)*((oac_fishf(i) + 0.00001)*log(oac_fishf(i) + 0.00001));
      }

     for (i=1; i<=nyrs_srv1_age; i++)
       {
         oac_srv1m(i)/=sum(oac_srv1m(i));
         oac_srv1f(i)/=sum(oac_srv1f(i));
         offset(3) -= nsamples_srv1_agem(i)*((oac_srv1m(i) + 0.00001)*log(oac_srv1m(i) + 0.00001));
         offset(4) -= nsamples_srv1_agef(i)*((oac_srv1f(i) + 0.00001)*log(oac_srv1f(i) + 0.00001));
       }

 END_CALCS



INITIALIZATION_SECTION

 // Nothing to see here... move along...

PARAMETER_SECTION

  // |---------------------------------------------------------------------------|
  // | POPULATION PARAMETERS
  // |---------------------------------------------------------------------------|
  // | - theta(1) -> log_mean_rec
  // | - theta(2) -> log_mean_y1
  // | - theta(3) -> log_avg_F
  // | - theta(4) -> log_comm_q
  // | - theta(5) -> log_surv_q1
  // | - theta(6) -> log_surv_q2

     init_bounded_number_vector theta(1,n_theta,theta_lb,theta_ub,theta_phz);
     number log_mean_rec;
     number log_mean_y1;
     number log_avg_F;
     number log_comm_q;
     number log_surv_q1;
     number log_surv_q2;

     number comm_q;
     number surv_q1;
     number surv_q2;



  // |---------------------------------------------------------------------------------|
  // | RECRUITMENT AND MORTALITY PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | 
  // | - log_mean_rec  <- mean recruitment
  // | - log_rec_dev   <- recruitment deviation vector   

     init_bounded_vector   init_pop(1,nages-1,-5,5,ph_init);
     init_bounded_vector   log_rec_dev(styr,endyr,-5,5,ph_rec);

     //vector        l_var_cpue1(1,nyrs_cpue1)
     //vector        l_var_cpue2(1,nyrs_cpue2)
     //vector        l_var_fshy(1,nyrs_fshy_cpue)
     vector        l_var(1,nyrs_MR)
     vector        l_var_SRV(1,nyrs_MR)


  // |---------------------------------------------------------------------------------|
  // | FISHING MORTALITY PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | 
  // | - log_F_devsf        <- annual fishing mortality deviations females
  // | - log_F_devsm        <- annual fishing mortality deviations males
  
     init_bounded_vector  log_F_devsf(styr,endyr,-5.,5.,ph_F);
     init_bounded_vector  log_F_devsm(styr,endyr,-5.,5.,ph_F);
    
     vector               Fmortf(styr,endyr);
     vector               Fmortm(styr,endyr);


  // |---------------------------------------------------------------------------------|
  // | POPULATION VECTORS
  // |---------------------------------------------------------------------------------|
  // | - tot_n          -> total abundance
  // | - tot_nf         -> total abundance females
  // | - tot_nm         -> total abundance males
  // | - tot_biomass    -> total biomass
  // | - pred_catch     -> total annual catch
  // | - pred_catchf    -> total annual catch females
  // | - pred_catchm    -> total annual catch males
  // |
  // | - natage_proj    -> numbers at age for projected year
  // | - natage_projf   -> numbers at age for projected year females
  // | - natage_projm   -> numbers at age for projected year males
	
     sdreport_vector  tot_n(styr,endyr)
     sdreport_vector  tot_nf(styr,endyr)
     sdreport_vector  tot_nm(styr,endyr)
     sdreport_vector  tot_biomass(styr,endyr)
     vector           pred_catch(styr,endyr)		
     vector           pred_catchf(styr,endyr)
     vector           pred_catchm(styr,endyr)

     sdreport_vector  natage_proj(1,nages);
     sdreport_vector  natage_projf(1,nages);
     sdreport_vector  natage_projm(1,nages);


  // |---------------------------------------------------------------------------------|
  // | POPULATION MATRICES
  // |---------------------------------------------------------------------------------|
  // | - natage         -> abundance-at-age
  // | - natagef        -> abundance-at-age females
  // | - natagem        -> abundance-at-age males
  // | - catege         -> catch-at-age
  // | - catagef        -> catch-at-age females
  // | - catagem        -> catch-at-age males
  // |
  // | - Zf             -> total mortality-at-age females
  // | - Ff             -> fishing mortality-at-age females
  // | - Sf             -> survival-at-age females
  // | - Zm             -> total mortality-at-age males
  // | - Fm             -> fishing mortality-at-ageemales
  // | - Sm             -> survival-at-age males
  // | 
  // | - eac_fishf      -> proportion-at-age fishery females
  // | - eac_fishm      -> proportion-at-age fishery males
  // | - eac_srv1f      -> proportion-at-age survey females
  // | - eac_srv1m      -> proportion-at-age survey males

     matrix  natage(styr,endyr,1,nages);
     matrix  natagef(styr,endyr,1,nages);
     matrix  natagem(styr,endyr,1,nages);
     
     matrix  catage(styr,endyr,1,nages);
     matrix  catagef(styr,endyr,1,nages);	
     matrix  catagem(styr,endyr,1,nages);
     
     matrix  Zf(styr,endyr,1,nages); 
     matrix  Ff(styr,endyr,1,nages); 
     matrix  Sf(styr,endyr,1,nages); 
     matrix  Zm(styr,endyr,1,nages);
     
     matrix  Fm(styr,endyr,1,nages); 
     matrix  Sm(styr,endyr,1,nages); 

     matrix  eac_fishf(1,nyrs_fish_age,1,nages); 
     matrix  eac_fishm(1,nyrs_fish_age,1,nages);
     matrix  eac_srv1f(1,nyrs_srv1_age,1,nages); 
     matrix  eac_srv1m(1,nyrs_srv1_age,1,nages);


  // |---------------------------------------------------------------------------------|
  // | STATISTICAL VECTORS
  // |---------------------------------------------------------------------------------|
  // | - effn_fish_agef       -> effective numbers-at-age fishery females
  // | - effn_fish_agem       -> effective numbers-at-age fishery males
  // | - sdnr_fish_agef       -> SDNR numbers-at-age fishery females
  // | - sdnr_fish_agem       -> SDNR numbers-at-age fishery males

     vector effn_fish_agef(1,nyrs_fish_age)
     vector effn_fish_agem(1,nyrs_fish_age)
     vector sdnr_fish_agef(1,nyrs_fish_age)
     vector sdnr_fish_agem(1,nyrs_fish_age)

  // |---------------------------------------------------------------------------------|
  // | MARK-RECAPTURE AND CPUE
  // |---------------------------------------------------------------------------------|
  // | - pred_cpue1      -> Survey CPUE (both sexes) 1988 - 1996
  // | - pred_cpue2      -> Survey CPUE (both sexes) 1997 - present
  // |
  // | - pred_fshy_pue   -> Fishery CPUE (both sexes)
  // |
  // | - pred_MR_abund   -> MR abundance

     vector  pred_cpue1(1,nyrs_cpue1); 
     vector  pred_cpue2(1,nyrs_cpue2);

     vector  pred_fshy_cpue(1,nyrs_fshy_cpue);

     sdreport_vector  pred_MR_abund(1,nyrs_MR);
     sdreport_vector  pred_MR_abund_SRV(1,nyrs_MR);

     //matrix  mr_f(1,nyrs_MR,1,nages-2);
     //matrix  mr_m(1,nyrs_MR,1,nages-2);

     matrix  mr_f_cpue(1,nyrs_fshy_cpue,1,nages-2);
     matrix  mr_m_cpue(1,nyrs_fshy_cpue,1,nages-2);

     vector  mr_f_proj(1,nages-2);
     vector  mr_m_proj(1,nages-2);
     vector  wt_f_proj(1,nages-2);
     vector  wt_m_proj(1,nages-2);


     number  s_srv1f;   // fraction females surviving from beginning of year to time of survey
     number  s_fshyf;   // fraction females surviving from beginning of year to time of fishery
     number  s_srv1m;   // fraction males surviving from beginning of year to time of survey
     number  s_fshym;   // fraction males surviving from beginning of year to time of fishery

  
  // |---------------------------------------------------------------------------------|
  // | SDREPORT QUANTITIES - legacy code: revisit these and clean up 2015
  // |---------------------------------------------------------------------------------|
  // |  - Age4_biom          <- total Age 4+ biomass
  // |  - expl_biom          <- exploitable biomass
  // |  - q_cpue1            <- q for survey 1 catchability 
  // |  - q_cpue2            <- q for survey 2 catchability 
  // |  - q_fshy            <- q for fishery catchability both sexes
  // |  - pred_rec           <- predicted recruitments
  // |  - avg_rec            <- mean recruitment 
  // |  - spbiom_trend       <- trend in biomass over last 2 years (B(t)/B(t-1); t=endyr)
  // |  - expl_biom_proj     <- next year's projected exploitable biomass 
  // |  - Age4_biom_proj     <- next year's projected Age 4+ biomass 
  // |  - spawn_biom_proj    <- next year's projected spawning biomass
  // |  - spawn_biom	   <- spawning biomass vector
  // |  - F50		   <- F50
  // |  - F40		   <- F40
  // |  - F35		   <- F35

     sdreport_vector  Age4_biom(styr,endyr);
     sdreport_vector  expl_biom(styr,endyr);
     sdreport_number  q_cpue1;
     sdreport_number  q_cpue2;
     sdreport_number  q_fshy;
     sdreport_vector  recruitment(styr,endyr);
     sdreport_number  avg_rec;
   //sdreport_number  spbiom_trend;
     vector           expl_abund_proj_comp(1,nages-2);
     sdreport_number  expl_abund_proj; 
     sdreport_number  Age4_biom_proj;
     vector           expl_biom_proj_comp(1,nages);
     sdreport_number  expl_biom_proj; 
     matrix        spawn_biom(styr,endyr,1,nages);
     sdreport_number  F50;
     sdreport_number  F40;
     sdreport_number  F35;


  // |---------------------------------------------------------------------------------|
  // | SPR / ABC
  // |---------------------------------------------------------------------------------|
  // |  - mF50            <- F50
  // |  - mF40            <- F40
  // |  - mF35            <- F35
  // |  - SB0             <- Spawning biomass per recruit at F = 0
  // |  - SBF50           <- Spawning biomass per recruit at F = 50
  // |  - SBF40           <- Spawning biomass per recruit at F = 40
  // |  - SBF35           <- Spawning biomass per recruit at F = 35
  // |  - SSB_0           <- total spawning biomass at F = 0
  // |  - SSB_40          <- total spawning biomass at F = 40
  // |
  // |  - Nspr            <- # spawners per recruit @ varying F
  // |  - wt_mature_s     <- weight of mature fish survey
  // |  - wt_mature_f     <- weight of mature fish fishery
  // |  - wt_fsh_sel_m	<- 
  // |  - wt_fsh_sel_f    <- 
  // |
  // |  - Fself		<- Full recruitment fishing mortality at F40 females
  // |  - Fselm		<- Full recruitment fishing mortality at F40 males
  // |  - ABC             <- Estimate of projected ABC
  // |  - Z1              <- total mortality at F = 40 females
  // |  - Z2              <- total mortality at F = 40 males
  // |  - S1              <- survival at F = 40 females
  // |  - S2              <- survival at F = 40 males
           //number mF50;
           //number mF40;
           //number mF35;
           //number              SB0;
           //number              SBF50;
           //number              SBF40;
           //number              SBF35;
           //number              SSB_0;
           //number              SSB_40;
     number    B40;
         
           //matrix              Nspr(1,4,1,nages);	
     vector              wt_mature_s(1,nages);
     vector              wt_mature_f(1,nages);
     vector              wt_fsh_sel_m(1,nages);
     vector              wt_fsh_sel_f(1,nages);
  
     vector              Fself(1,nages-2);
     vector              Fselm(1,nages-2);
     vector              ABC_comp(1,nages);
     sdreport_number     ABC;
     vector              ABC_all(styr,endyr);
     vector              Z1(1,nages-2);	
     vector              S1(1,nages-2);		
     vector              Z2(1,nages-2);	
     vector              S2(1,nages-2);

  // |-------------------------------------------------------------------------|
  // | PROJECTED POPULATIONS, SPR AND F LEVELS                                 |
  // |-------------------------------------------------------------------------|
  // | mF                 -> vector of potential F levels
  // | SBx                -> spawning biomass from mF            
  // | Bx                 -> biomass from mF
  // |
  // | Nspr               -> Number of projected spawners
  // | N0                 -> Number projected with F = 0

     init_bounded_vector mF(1,7,0.01,1)
     vector SBx(1,8)
     vector Bx(1,8)
					       
     matrix Nspr(1,8,1,nages)			
     matrix N0(1,8,1,nages)


  // |---------------------------------------------------------------------------------|
  // | LIKELIHOODS
  // |---------------------------------------------------------------------------------|
  // | - sprpen             <- penalty for SPR function to estimate F values
  // | - priors             <- vector of priors (q, M, etc.)
  // | - surv_like          <- survey likelihoods (CPUE1, CPUE2, fishery CPUE, MR)
  // | - age_like           <- age composition likelihoods (fishery and survey-at-age)
  // | - rec_like           <- recruitment likelihood
  // | - ssqcatch           <- total annual catch likelihood
  // | - F_mort_regularity  <- constraint on F level variability
  // | 
  // | - Like               <- add it all up
  // | - obj_fun            <- ADMB objective function call

     number  sprpen;
     vector  penalties(1,4);
     vector  surv_like(1,5);
     vector  age_like(1,4);
     number  rec_like;
     vector  catch_like(1,2);
     number  F_mort_regularity;

     number  Like;
     objective_function_value obj_fun;	

						
 LOCAL_CALCS

   wt_mature_f = elem_prod(p_mature,wt_avg_fshy_f);
   wt_mature_s = elem_prod(p_mature,wt_avg_srv1_f);

   wt_fsh_sel_f = elem_prod(fish_sel_f,wt_avg_fshy_f);
   wt_fsh_sel_m = elem_prod(fish_sel_m,wt_avg_fshy_m);


  // |-------------------------------------------------------------------------|
  // | Survival fractions to decrement abundance relative to annual timing
  // |-------------------------------------------------------------------------|
  // |
  // | This will make a difference if distinct natural mortalities are used
  // | between males and females. If both are input at 0.1, then one decrement
  // | for the fishery and another for the survey is all that is required

     s_srv1f  = mfexp(-1*srv1_fract*Mf);
     s_fshyf  = mfexp(-1*fshy_fract*Mf);
     s_srv1m  = mfexp(-1*srv1_fract*Mm);
     s_fshym  = mfexp(-1*fshy_fract*Mm);

 END_CALCS


PROCEDURE_SECTION
     initializeModelParameters();

     Mortality();
            if(DEBUG_FLAG == 1) cout<<"**Mortality**"<<endl;
     Abundance();	
            if(DEBUG_FLAG == 1) cout<<"**Abundance**"<<endl;
     Catch();	
            if(DEBUG_FLAG == 1) cout<<"**Catch**"<<endl;
     Predicted();
            if(DEBUG_FLAG == 1) cout<<"**Predicted**"<<endl;
		
     if (last_phase())
       {
         spr();
     	    if(DEBUG_FLAG == 1) cout<<"**spr**"<<endl;
                
         Projection();
            if(DEBUG_FLAG == 1) cout<<"**Projection**"<<endl;

         Get_ABC();
            if(DEBUG_FLAG == 1) cout<<"**ABC**"<<endl;
       }

     Objective_Function();


     if(mceval_phase())	
       {
         writePosteriorSamples();

         evalout<<tot_n<<" "<<
              tot_nf<<" "<<
              tot_nm<<" "<<
              pred_MR_abund<<""<<
              pred_MR_abund_SRV<<""<<
              tot_biomass<<" "<<
              recruitment<<" "<<
              expl_abund_proj<<" "<<
              expl_biom_proj<<" "<<
              endl;
       }

FUNCTION void writePosteriorSamples()
	/**
	- This function is only envoked when the -mceval
		command line option is implemented.
	*/
	//ofstream ofs0("bmT.ps",ios::app);
	//ofs0<<tot_biomassT<<endl;

FUNCTION void initializeModelParameters()
	//fpen = 0;

	log_mean_rec          = theta(1);
	log_mean_y1           = theta(2);
	log_avg_F             = theta(3);
        log_comm_q            = theta(4);
        log_surv_q1           = theta(5);
        log_surv_q2           = theta(6);


FUNCTION Mortality
  Fmortf.initialize();
  Fmortm.initialize();
  Ff.initialize();
  Fm.initialize();
  Zf.initialize();
  Zm.initialize();
  Sf.initialize();
  Sm.initialize();

  // |----------------------------------------------------------------------|
  // | ANNUAL DEVIATIONS FOR COMMERCIAL FISHERY
  // |----------------------------------------------------------------------|
  
     Fmortf        = mfexp(log_avg_F +  log_F_devsf);
     Fmortm        = mfexp(log_avg_F +  log_F_devsm);


  // |----------------------------------------------|
  // | Fishing Mortality matrices
  // |----------------------------------------------|
  // | Note here different mortalities by gender
  // | due to selectivity - input natural mortality
  // | is the same for both as per Johnson and Quinn
  // | (1988) and Hanselman et al. (2014)
  
     for (i=styr; i<=endyr; i++)
       {
         for (j=1 ; j<=nages; j++)
           {
             Ff(i,j) = Fmortf(i) * fish_sel_f(j);
             Fm(i,j) = Fmortm(i) * fish_sel_m(j);
           }
       }		


  // |----------------------------------------------|
  // | Total Mortality & Survival
  // |----------------------------------------------|
     
     Zf  = Ff+Mf;	
     Sf  = mfexp(-1.0*Zf);
  
     Zm  = Fm+Mm;	
     Sm  = mfexp(-1.0*Zm);


FUNCTION Abundance
  natagef.initialize();
  natagem.initialize();
  natage.initialize();
  tot_nm.initialize();
  tot_nf.initialize();
  tot_n.initialize();

  // |-------------------------------------------------------------------------|
  // | YEAR 1
  // |-------------------------------------------------------------------------|
  // | Equal proportions at recruitment are assumed
  // |
  // | Recruitment variability set to 1.2 as per Sigler et al. (2002) and
  // | Hanselman et al. (2014)

     natagef(styr,1)= 0.5 * (mfexp(log_rec_dev(styr) + log_mean_rec  //FEMALES
                            + (sigr*sigr)/2)/2);

     natagem(styr,1)= 0.5 * (mfexp(log_rec_dev(styr) + log_mean_rec  //MALES
                            + (sigr*sigr)/2)/2);

     for(j=1;j<=nages-1;j++)
       {
         natagef(styr,j+1)= 0.5 * (mfexp(log_mean_y1 + init_pop(j)   //FEMALES
                                  + (sig1*sig1)/2)/2);
                                  
         natagem(styr,j+1)= 0.5 * (mfexp(log_mean_y1 + init_pop(j)   //MALES
                                  + (sig1*sig1)/2)/2);
       }
        

  // |-------------------------------------------------------------------------|
  // | Subsequent years
  // |-------------------------------------------------------------------------|
  // | The ++ notation increments index range by 1 so that computations
  // |  for (1,nages-1) correspond to index range (2,nages)

     for (i=styr+1; i<=endyr; i++)
       {
         natagef(i,1) = 0.5 * mfexp(log_rec_dev(i) + log_mean_rec + (sigr*sigr)/2)/2;  
         natagem(i,1) = 0.5 * mfexp(log_rec_dev(i) + log_mean_rec + (sigr*sigr)/2)/2; 

       for(j=2; j<=nages; j++)
         {
           natagef(i,j)  = natagef(i-1,j-1)*Sf(i-1,j-1);
           natagem(i,j)  = natagem(i-1,j-1)*Sm(i-1,j-1); 
         }                                                                         
										  
       natagef(i,nages)    += natagef(i,nages)*Sf(i,nages);
       natagem(i,nages)    += natagem(i,nages)*Sm(i,nages);   
       }


    // |-----------------------------------------------------------------------|
    // | Summations
    // |-----------------------------------------------------------------------|

       for(i=styr; i<=endyr; i++)
         {
           for(j=1; j<=nages; j++)
             {
               natage(i,j) = natagef(i,j)+natagem(i,j);
               tot_nm(i)   = sum(natagem(i));
               tot_nf(i)   = sum(natagef(i));
               tot_n(i)    = sum(natage(i));
             }
         }

       
FUNCTION Catch
  pred_catchf.initialize();
  pred_catchm.initialize();
  pred_catch.initialize();
  catagef.initialize();
  catagem.initialize();
  catage.initialize();

  // |----------------------------------------------------------------------|
  // | COMMERCIAL CATCH-AT-AGE (note that age comps for OF are farther below)
  // |----------------------------------------------------------------------|
  // | Total annual catch from summed catch-at-age

     for (i=styr; i<=endyr; i++)
       {
         catagef(i)     = elem_div(elem_prod(elem_prod(natagef(i),Ff(i)),(1.-Sf(i))),Zf(i));
         pred_catchf(i) =  (catagef(i)*wt_avg_fshy_f);
    
         catagem(i)     = elem_div(elem_prod(elem_prod(natagem(i),Fm(i)),(1.-Sm(i))),Zm(i));
         pred_catchm(i) =  (catagem(i)*wt_avg_fshy_m);
    
         catage(i)      = catagef(i) + catagem(i);
         pred_catch(i)  = pred_catchf(i) + pred_catchm(i);
       }

 
FUNCTION Predicted



  // |-------------------------------------------------------------------------|
  // | Catchabilities
  // |-------------------------------------------------------------------------|
  // | Note that catchabilities are assumed *constant* between sexes
  // | Differences in retention are accounted as selectivity

     surv_q1 = mfexp(log_surv_q1);
     surv_q2 = mfexp(log_surv_q2); 
     comm_q  = mfexp(log_comm_q);


  // |-------------------------------------------------------------------------|
  // | MARK-RECAPTURE ABUNDANCE - SURVEY YEARS
  // |-------------------------------------------------------------------------|
  // | Note that the time decrement IS implemented here. The mark-recapture survey
  // | Chapman estimate of abundance is a snap-shot of abundance at the
  // | beginning of the commercial fishery.
  // |
  // | Note that NO catchability coefficient is used.
  // |
  // | The '1000' multiplier is for scale. Weights are kilograms, so abundance
  // | scale is in 1000s. The Chapman point estimates and variances are in
  // | real numbers and input that way to make sure nothing is messed up in
  // | an attempt to scale the variance
  // |
  // | The 'mr_f' and 'mr_m' matrices are to scale the mark-recapture abundances
  // | by the recapture probability (i.e. commercial selectivity-at-age)
  // | The youngest age retained in the commercial fishery is age 4, but the
  // | model begins with age 2, thus the '-2'. This is not necessary with the
  // | longline survey, as the youngest there is age 2. 

     dvar_matrix mr_f(1,nyrs_MR,1,nages-2);
     dvar_matrix mr_m(1,nyrs_MR,1,nages-2);
     dvar_matrix mr_f_SRV(1,nyrs_MR,1,nages);
     dvar_matrix mr_m_SRV(1,nyrs_MR,1,nages);
  
     for (i=1;i<=nyrs_MR;i++)
       {
         for(j=3; j<=nages; j++)
           {
             mr_f(i,j-2) = natagef(yrs_MR(i),j)*fish_sel_f(j);
             mr_m(i,j-2) = natagem(yrs_MR(i),j)*fish_sel_m(j);
           }
       }

     for (i=1;i<=nyrs_MR;i++)
       {
             mr_f_SRV(i) = elem_prod(natagef(yrs_MR(i)),srv1_sel_f);
             mr_m_SRV(i) = elem_prod(natagem(yrs_MR(i)),srv1_sel_m);
       }


     for (i=1;i<=nyrs_MR;i++)
       {
         pred_MR_abund(i) =   s_fshyf * 1000 * sum(mr_f(i))+
                              s_fshym * 1000 * sum(mr_m(i));

         pred_MR_abund_SRV(i) =   s_srv1f * 1000 * sum(mr_f_SRV(i))+
                                  s_srv1m * 1000 * sum(mr_m_SRV(i));


        //1000* sum ((s_srv1f *elem_prod(natagef(yrs_MR(i)),srv1_sel_f))+
        //(s_srv1m * elem_prod(natagem(yrs_MR(i)),srv1_sel_m))) =   obs_MR_abund_SRV(i);
                              
       } 
                              
      
  // |-------------------------------------------------------------------------|
  // | ADF&G LONGLINE SURVEY CPUE 1988 - 1996
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance

     for (i=1;i<=nyrs_cpue1;i++)
       {
         pred_cpue1(i) = surv_q1 * sum(elem_prod(natagef(yrs_cpue1(i)),
                                       wt_avg_srv1_f)+  
                                       elem_prod(natagem(yrs_cpue1(i)),
                                       wt_avg_srv1_m));  
       }


  // |-------------------------------------------------------------------------|
  // | ADF&G LONGLINE SURVEY CPUE 1997 - PRESENT
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance
					             
     for (i=1;i<=nyrs_cpue2;i++)
       {
         pred_cpue2(i) = surv_q2 * sum(elem_prod(natagef(yrs_cpue2(i)),
                                       wt_avg_srv1_f)+  
                                       elem_prod(natagem(yrs_cpue2(i)),
                                       wt_avg_srv1_m));  
       }


  // |-------------------------------------------------------------------------|
  // | COMMERCIAL FISHERY CPUE
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance


     for (i=1;i<=nyrs_fshy_cpue;i++)
       {

         pred_fshy_cpue(i) =  comm_q * sum(elem_prod(natagef(yrs_fshy_cpue(i)),wt_fsh_sel_f)
                                           + 
                                           elem_prod(natagem(yrs_fshy_cpue(i)),wt_fsh_sel_m));
       }


  // |-----------------------------------------------------------------------|
  // | AGE COMPS, N, EFFN, SDNR
  // |-----------------------------------------------------------------------|


     //FISHERY
     for (i=1;i<=nyrs_fish_age;i++)
       {
         eac_fishm(i)  = catagem(yrs_fish_age(i)) / sum(catagem(yrs_fish_age(i))) * ageage;   
         eac_fishf(i)  = catagef(yrs_fish_age(i)) / sum(catagef(yrs_fish_age(i))) * ageage;  

         effn_fish_agem(i) = (1-eac_fishm(i))*eac_fishm(i)/norm2(oac_fishm(i)-eac_fishm(i)); 
         effn_fish_agef(i) = (1-eac_fishf(i))*eac_fishf(i)/norm2(oac_fishf(i)-eac_fishf(i)); 
       }


     // SURVEY
     for (i=1;i<=nyrs_srv1_age;i++)
       {
         eac_srv1f(i)  = elem_prod(srv1_sel_f, natagef(yrs_srv1_age(i))*s_srv1f)
                         /sum(elem_prod(natagef(yrs_srv1_age(i))*s_srv1f,srv1_sel_f)) * ageage; 

         eac_srv1m(i)  = elem_prod(srv1_sel_m, natagem(yrs_srv1_age(i))*s_srv1f)
                         /sum(elem_prod(natagem(yrs_srv1_age(i))*s_srv1m,srv1_sel_m)) * ageage;
       }

FUNCTION Population
   recruitment.initialize();
   tot_biomass.initialize();
   Age4_biom.initialize();
   expl_biom.initialize();
   spawn_biom.initialize();

  // |-------------------------------------------------------------------------|
  // | MAKE SURE UNITS ARE CORRECT
  // |-------------------------------------------------------------------------|
  // |*Note conversion from metric tonnes (for total biomass) to straight pounds
  // |  for Age 4 and Exploitable biomass calcs - this is for managers
  // |
  // |  - kg (wt-at-age)  x thousands (natage)  = metric tonnes
  // |  - 1 metric tonne = 2,204 pounds
  // |  - Age 4 & Exploitable = pounds (instead of thousands of pounds)
  // |  - Explotable biomass uses fishery selectivity
  // |  - These estimates are snapshots at the time of the longline survey, which 
  // |    are assumed equivalent to the conditions upon the start of the
  // |    year's commercial fishery


     for (i=styr;i<=endyr;i++)
       {
         recruitment(i) = natage(i,1);
    
         tot_biomass(i) = (natagef(i)*wt_avg_srv1_f)+(natagem(i)*wt_avg_srv1_m); //metric tonnes
    
         Age4_biom(i)   = 2204 * (wt_avg_srv1_f(3,nages) * (natagef(i)(3,nages) *
                                 mfexp(-srv1_fract*Mf))+
                                 wt_avg_srv1_m(3,nages) * (natagem(i)(3,nages) *
                                 mfexp(-srv1_fract*Mm)));

         expl_biom(i)   = 2204 * (wt_avg_fshy_f * elem_prod(natagef(i) *
                                  mfexp(-srv1_fract*Mf),fish_sel_f)+
                                  wt_avg_fshy_m * elem_prod(natagem(i) *
                                  mfexp(-srv1_fract*Mm),fish_sel_m));

    	       
         spawn_biom(i)  = 2204 * (elem_prod(wt_mature_s,natagef(i))*mfexp(-spawn_fract*Mf)); 
       }


FUNCTION Penalties

  penalties.initialize();


    // |------------------------------------------------------------------------|
    // | Recruitment penalty
    // |------------------------------------------------------------------------|
    
       penalties(1)  = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) +
                      size_count(log_rec_dev)*log(sigr);


    // |------------------------------------------------------------------------|
    // | Year 1 abundance-at-age penalty
    // |------------------------------------------------------------------------|
    
       penalties(2)  = norm2(init_pop+1*sig1*sig1/2.)/(2.*square(sig1)) +
                       size_count(init_pop)*log(sig1);
                       

    // |------------------------------------------------------------------------|
    // | Fisheries deviations penalty
    // |------------------------------------------------------------------------|
       int fpen = 0.2;
       if (last_phase())
         {
           fpen = 1;
         }
      
       penalties(3)  = 0.5*log(2*M_PI) + log(fpen) +
                       0.5*(norm2(log_F_devsf)) / (2*square(fpen));

       penalties(4)  = 0.5*log(2*M_PI) + log(fpen) +
                       0.5*(norm2(log_F_devsm)) / (2*square(fpen));


FUNCTION Catch_Likelihood

  catch_like.initialize();

  
  catch_like(1)  +=  0.5*log(2*M_PI)  +log(sigma_catch) +  0.5*(norm2(log(obs_catchf+oo) -
                                     log(pred_catchf+oo))) / (2.*square(sigma_catch));

  catch_like(2)  +=  0.5*log(2*M_PI)  +log(sigma_catch) +  0.5*(norm2(log(obs_catchm+oo) -
                                     log(pred_catchm+oo))) / (2.*square(sigma_catch));



FUNCTION Surv_Likelihood

  surv_like.initialize();

  // |----------------------------------------------------------------------|
  // | Mark-recapture estimates of abundance - normal
  // |----------------------------------------------------------------------|

    for (i=1; i<=nyrs_MR; i++)
      {
   
         surv_like(1) += 0.5*log(2*M_PI) + log(sqrt(obs_MR_var(i))) +
                    0.5*norm2(obs_MR_abund-pred_MR_abund)
                   / (2*obs_MR_var(i));

         surv_like(2) += 0.5*log(2*M_PI) + log(sqrt(obs_MR_var_SRV(i))) +
                    0.5*norm2(obs_MR_abund_SRV-pred_MR_abund_SRV)
                   / (2*obs_MR_var_SRV(i));
                   
      }

     //for (i=1; i<=nyrs_MR; i++)
     //  {
     //    l_var(i)     = log(1. + (square(obs_MR_var(i))/
     //                    square(obs_MR_abund(i))));
     //                     
     //    l_var_SRV(i) = log(1. + (square(obs_MR_var_SRV(i))/
     //                    square(obs_MR_abund_SRV(i))));
     //
     //
     //    surv_like(1) += 0.5*log(2*M_PI) + log(sqrt(l_var(i))) +
     //                    0.5*(norm2(log(obs_MR_abund) - log(pred_MR_abund))) /
     //                    (2*(l_var(i)));
     //
     //    surv_like(2) += 0.5*log(2*M_PI) + log(sqrt(l_var_SRV(i))) +
     //                     0.5*(norm2(log(obs_MR_abund_SRV) - log(pred_MR_abund_SRV))) /
     //                    (2*(l_var_SRV(i)));
     //  }



  // |----------------------------------------------------------------------|
  // | ADF&G survey CPUE - normal 1980 - 1996
  // |----------------------------------------------------------------------|

  for (i=1; i<=nyrs_cpue1; i++)
    {

       surv_like(3) += 0.5*log(2*M_PI) + log(sqrt(obs_cpue1_var(i))+oo) +
                       0.5*(square(obs_cpue1_biom(i)-pred_cpue1(i)))
                       / (2*obs_cpue1_var(i)+oo);
    }


  // |----------------------------------------------------------------------|
  // | ADF&G survey  CPUE - normal 1996- present
  // |----------------------------------------------------------------------|
  
  for (i=1; i<=nyrs_cpue2; i++)
    {

       surv_like(4) += 0.5*log(2*M_PI) + log(sqrt(obs_cpue2_var(i))+oo) +
                       0.5*(square(obs_cpue2_biom(i)-pred_cpue2(i)))
                       / (2*obs_cpue2_var(i)+oo);
    }


  // |----------------------------------------------------------------------|
  // | COMMERCIAL CPUE - normal 
  // |----------------------------------------------------------------------|

  for (i=1; i<=nyrs_fshy_cpue; i++)
    {

       surv_like(5) += 0.5*log(2*M_PI) + log(sqrt(obs_fshy_cpue_var(i))+oo) +
                       0.5*(square(obs_fshy_cpue(i)-pred_fshy_cpue(i)))
                       / (2*obs_fshy_cpue_var(i)+oo);
    }
    


FUNCTION Multinomial_Likelihood

  age_like.initialize();

  // |----------------------------------------------------------------------|
  // | MULTINOMIAL LIKELIHOODS FOR AGE COMPOSITIONS
  // |----------------------------------------------------------------------|

    // |--------------------------------------------------------------------|
    // | Commercial fishery age compositions
    // |--------------------------------------------------------------------|

       for (i=1; i <= nyrs_fish_age; i++)
         {
           age_like(1) -= nsamples_fish_agem(i)*((oac_fishm(i) + oo) *
                                                  log(eac_fishm(i) + oo)) ;
                                                  
           age_like(2) -= nsamples_fish_agef(i)*((oac_fishf(i) + oo) *
                                                  log(eac_fishf(i) + oo)) ; 
         }

    // |--------------------------------------------------------------------|
    // | Longline survey age compositions
    // |--------------------------------------------------------------------|

       for (i=1; i <= nyrs_srv1_age; i++)
         {
           age_like(3) -= nsamples_srv1_agem(i)*((oac_srv1m(i) + oo) *
                                                  log(eac_srv1m(i) + oo)) ;
                                                  
           age_like(4) -= nsamples_srv1_agef(i)*((oac_srv1f(i) + oo) *
                                                  log(eac_srv1f(i) + oo)) ; 
         }


  age_like   -= offset; 


FUNCTION Objective_Function

  Like.initialize();

  // |----------------------------------------------------------------------|
  // | CALL LIKELIHOOD FUNCTIONS
  // |----------------------------------------------------------------------|
  
     Catch_Likelihood();
     Multinomial_Likelihood();
     Surv_Likelihood();
     Penalties();


  // |----------------------------------------------------------------------|
  // | SUM DATA LIKELIHOODS
  // |----------------------------------------------------------------------|

     Like           = sum(catch_like);	
     Like           += sum(surv_like);	
     Like           += sum(age_like);	

     obj_fun    = Like;         // Put here to capture the data likelihood

     obj_fun      += sum(penalties);

     if (last_phase())
        obj_fun += sprpen; 


FUNCTION spr

  SBx=0.;
  Bx.initialize();
  Nspr.initialize();
  N0.initialize();
  sprpen.initialize();

  // |----------------------------------------------------------------------|
  // | SPR RATES - ADDED TO LIKELIHOOD AS 'sprpen'
  // |----------------------------------------------------------------------|
  // | The code below has been modded a bit to match the R script (YPR.R)
  // |  which, in turn, was written to reflect the older Excel
  // |  yield-per-recruit file. Thus, while Nspr(foo,1) can be set equal to
  // |  1, 500, or a bajillion, it is set to 500.
          

  // |----------------------------------------------------------------------|
  // | RECRUITMENT (AGE 2)
  // |----------------------------------------------------------------------|

     for (i=1;i<=8;i++)
       {
         Nspr(i,1) = 500.;
         N0(i,1) = mean(recruitment(1982,endyr-recage));
       }


  // |----------------------------------------------------------------------|
  // | AGES 3 - 42
  // |----------------------------------------------------------------------|

     for (j=2;j<=nages;j++)
       {
         Nspr(1,j) = Nspr(1,j-1) * mfexp(-1.*Mf);
         N0(1,j)   = N0(1,j-1)   * mfexp(-1.*Mf);

       for (i=2; i<=8; i++)
         {
           Nspr(i,j) = Nspr(i,j-1) * mfexp(-1.*(Mf+mF(i-1)
                                                *fish_sel_f(j-1)));
                                                
           N0(i,j)   = N0(i,j-1)   * mfexp(-1.*(Mf+mF(i-1)
                                                *fish_sel_f(j-1)));
         }
       } // close 'j' loop


  // |----------------------------------------------------------------------|
  // | PLUS CLASS
  // |----------------------------------------------------------------------|

      Nspr(1,nages) += Nspr(1,nages) * mfexp(-1.*Mf);
      N0(1,nages)   += N0(1,nages)   * mfexp(-1.*Mf);


     for(i=2; i<=8; i++)
      {
         Nspr(i,nages) += Nspr(i,nages) * mfexp(-1.*(Mf+mF(i-1)*
                                                fish_sel_f(nages)));
         N0(i,nages)   += N0(i,nages)   * mfexp(-1.*(Mf+mF(i-1)*
                                                fish_sel_f(nages)));
       }


  // |----------------------------------------------------------------------|
  // | PROJECTED BIOMASS FOR F LEVELS
  // |----------------------------------------------------------------------|

     SBx(1) +=  2.20462 * Nspr(1)*wt_mature_s;

     
     for (j=1;j<=nages;j++)
       {
         for(i=2; i<=8; i++)
           {
             SBx(i) += 2.20462 * Nspr(i,j)*wt_mature_s(j);
           }

         for(i=1; i<=8; i++)
           {
             Bx(i) +=  N0(i,j)*wt_mature_s(j);
           }
       } // close 'j' loop


  // |----------------------------------------------------------------------|
  // | CALCULATE F LEVELS
  // |----------------------------------------------------------------------|
  
     sprpen    = 100.*square(SBx(8)/SBx(1)-0.6);
     sprpen   += 100.*square(SBx(7)/SBx(1)-0.55);
     sprpen   += 100.*square(SBx(6)/SBx(1)-0.5);
     sprpen   += 100.*square(SBx(5)/SBx(1)-0.45);
     sprpen   += 100.*square(SBx(4)/SBx(1)-0.4);
     sprpen   += 100.*square(SBx(3)/SBx(1)-0.35);
     sprpen   += 100.*square(SBx(2)/SBx(1)-0.3);
  
     B40 = SBx(4)*mean(recruitment(1982,endyr-recage));
     F50 = mF(5);
     F40 = mF(3);
     F35 = mF(2);

     // |-------------------------------------|
     // | NOTE
     // |-------------------------------------|
     // | As a check, mF vector values (F50,
     // |  F40, etc.) match the output from
     // |  'YPR.R' to the fourth decimal place


FUNCTION Projection

  // |------------------------------------------------------------------------|
  // | ABUNDANCE FOR FIRST PROJECTION YEAR (endyr+1)
  // |------------------------------------------------------------------------|


 // |-------------------------------------------------------------------------|
  // | Population Projection (endyr + 1)
  // |-------------------------------------------------------------------------|

    // |-----------------------------------------------------------------------|
    // | Recruitment (age 2)
    // |-----------------------------------------------------------------------|
    // | Mean over recent 15 years, excluding last two model years (variability)

       natage_projf(1)  = 0.5 * mean(mfexp(log_mean_rec
                                + log_rec_dev(endyr-16, endyr-2)));

       natage_projm(1)  = 0.5 * mean(mfexp(log_mean_rec
                                + log_rec_dev(endyr-16, endyr-2)));
                                

    // |-----------------------------------------------------------------------|
    // | AGES 3 - 42
    // |-----------------------------------------------------------------------|

       //for(j=2; j<=nages; j++)
        // {
        //   natage_projf(j)  = natagef(endyr,j-1)*Sf(endyr,j-1);
        //   natage_projm(j)  = natagem(endyr,j-1)*Sm(endyr,j-1); 
        // }  

       natage_projf(2, nages) = ++elem_prod(natagef(endyr)(1,nages-1),
                                  Sf(endyr)(1,nages-1));

       natage_projm(2, nages) = ++elem_prod(natagem(endyr)(1,nages-1),
                                  Sm(endyr)(1,nages-1));

       
    // |-----------------------------------------------------------------------|
    // | PLUS CLASS
    // |-----------------------------------------------------------------------|

       natage_projf(nages)   += natage_projf(nages) * Sf(endyr,nages);
       natage_projm(nages)   += natage_projm(nages) * Sm(endyr,nages);
       

    // |-----------------------------------------------------------------------|
    // | Summations
    // |-----------------------------------------------------------------------|

       for(j=1; j<=nages; j++)
         {
            natage_proj(j) = natage_projf(j)+natage_projm(j);
         }

        for(j=3; j<=nages; j++)
          {
            mr_f_proj(j-2) = natage_projf(j);
            mr_m_proj(j-2) = natage_projm(j);
            }



FUNCTION Get_ABC

  // |-------------------------------------------------------------------------|
  // | Note that the calcs below reflect the ADF&G regulations for sablefish
  // |   and *NOT* standard NOAA federal curves as per Dorn
  // |
  // | Note that fishery selectivities and fishery weights used to compute ABC

     ABC.initialize();
     ABC_all.initialize();
     Age4_biom_proj.initialize();
     expl_abund_proj.initialize();
     expl_biom_proj.initialize();


    for(j=3; j<=nages; j++)
      {
      
     Fself(j-2) = F50*fish_sel_f(j);
     Fselm(j-2) = F50*fish_sel_m(j);
     wt_f_proj(j-2) = wt_avg_fshy_f(j);
     wt_m_proj(j-2) = wt_avg_fshy_m(j);
       }

     Z1 = Mf + Fself;
     S1 = mfexp(-Z1);
     Z2 = Mm + Fselm;
     S2 = mfexp(-Z2);


  // |-------------------------------------------------------------------------|
  // | ALLOWABLE BIOLOGICAL CATCH
  // |-------------------------------------------------------------------------|

     // |----------------------------------------------------------------------|
     // | Projected year (endyr + 1)
     // |----------------------------------------------------------------------|

        
        ABC = 2204 *sum(  ((elem_prod(wt_f_proj,elem_prod(elem_div(Fself, Z1),
                                       elem_prod(1-S1,(mr_f_proj))))) +
                       elem_prod(wt_m_proj,elem_prod(elem_div(Fselm, Z2),
                                       elem_prod(1-S2,(mr_m_proj))))) );


     // |----------------------------------------------------------------------|
     // | All modeled years for comparison with implemented catch levels
     // |----------------------------------------------------------------------|


  // |-------------------------------------------------------------------------|
  // | Mean recruitment
  // |-------------------------------------------------------------------------|

     avg_rec        = mean(recruitment(styr+2,endyr-recage));	 


  // |-------------------------------------------------------------------------|
  // | Trend in spawning biomass relative to last year
  // |-------------------------------------------------------------------------|

     //spbiom_trend   = 0;//spawn_biom(endyr)/spawn_biom(endyr-1);
	

  // |-------------------------------------------------------------------------|
  // | Next year's projected age 4+ biomass at time of survey
  // |-------------------------------------------------------------------------|

    Age4_biom_proj =  2204 * (wt_avg_srv1_f(3,nages) * (natage_projf(3,nages) *
                              mfexp(-srv1_fract*Mf)) +
                              wt_avg_srv1_m(3,nages) * (natage_projm(3,nages) *
                              mfexp(-srv1_fract*Mm)));

  // |-------------------------------------------------------------------------|
  // | Next year's projected exploited abundance at time of survey
  // |-------------------------------------------------------------------------|

        for (int j=3; j<=nages; j++)
        {
        expl_abund_proj_comp(j-2)  = natage_projf(j)*fish_sel_f(j)  +
                                   natage_projm(j)*fish_sel_m(j);
         }



     expl_abund_proj  = 1000 * sum(expl_biom_proj_comp); 

  // |-------------------------------------------------------------------------|
  // | Next year's projected exploited biomass
  // |-------------------------------------------------------------------------|

        for (int j=3; j<=nages; j++)
        {
     expl_biom_proj_comp(j-2) =  (natage_projf(j) * wt_mature_s(j));//fish_sel_f(j) *wt_avg_fshy_f(j)) +
                                   //(natage_projm(j) * fish_sel_m(j)*wt_avg_fshy_m(j));
     }
     expl_biom_proj = 2204 * sum(expl_biom_proj_comp);

	                  





TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=5000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(800000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


RUNTIME_SECTION
  maximum_function_evaluations 5000 5000 5000 5000 5000;
  convergence_criteria 0.0001;


FINAL_SECTION

   // |----------------------------------------------------------------------|
   // | Print run time statistics to the screen.
   // |----------------------------------------------------------------------|
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
        cout<<""<<endl;
	cout<<"--Objective function value: "<<obj_fun<<endl;
        cout<<""<<endl;
	cout<<"--Maximum gradient component: "<<objective_function_value::gmax<<endl;
        cout<<""<<endl;
	cout<<"*******************************************"<<endl;

  // |------------------------------------------------------------------------|
  // | LONG-ASS PARAMETRIC BOOTSTRAP CODE
  // |------------------------------------------------------------------------|


  // |------------------------------------------------------------------------|
  // | CALL IN/OUT FILE AND NAMESPACE
  // |------------------------------------------------------------------------|

     using namespace std;

     ifstream infile;
     ofstream outfile;


  // |------------------------------------------------------------------------|
  // | STRUCTURES
  // |------------------------------------------------------------------------|
  // |
  // |

     const int maxlen=46;
     const int maxline=5000;
     const int maxpar=initial_params::nvarcalc();
     char header_std[maxlen];
     char header_cor[maxline];
     char pname[maxlen];
     int icount;
     int ii;
     int jj;

     int nboot;
     int nwork;
     int ndump;
     int marker;
     int nparameter;
     int index;
     double estimate;
     double stdev;
     const int stuff=570;
     double tmp_cor[stuff];
     double tmp_shit[stuff];

     cout<<" "<<endl;
     cout<<"What's the story, mother?"<<endl;


  // |------------------------------------------------------------------------|
  // | RANDOM NUMBER GENERATOR
  // |------------------------------------------------------------------------|
  // | - uses current system time as random seed

     random_number_generator rng(time(NULL));


  // |------------------------------------------------------------------------|
  // | CONTAINERS
  // |------------------------------------------------------------------------|

     dvector Boot_index_init(1,stuff);
     dvector Boot_stdev_init(1,stuff);

     dvector Boot_index(1,maxpar);
     dvector Boot_stdev(1,maxpar);

     dmatrix Boot_corr(1,stuff,1,stuff);
     dmatrix Boot_correlation(1,maxpar,1,maxpar);
     dmatrix Boot_covariance(1,maxpar,1,maxpar);
     dmatrix Boot_choleski(1,maxpar,1,maxpar);

     dvector Boot_ran(1,maxpar);
     dvector Bootshit(1,maxpar);

     dvector draw(1,maxpar);
     dvector bootN(1,maxpar);
     dmatrix boot(1,nboot,1,maxpar);

     Boot_index_init.initialize();
     Boot_index.initialize();
     Boot_stdev_init.initialize();
     Boot_stdev.initialize();
     Boot_corr.initialize();
     Boot_correlation.initialize();
     Boot_covariance.initialize();
     Boot_choleski.initialize();
     //Boot_ran.initialize();
     draw.initialize();
     bootN.initialize();

     cout<<" "<<endl;
     cout<<"Containers loaded"<<endl;


  // |------------------------------------------------------------------------|
  // | OPEN THE 'MODEL.COR' FILE
  // |------------------------------------------------------------------------|

     infile.open("model.cor");


  // |------------------------------------------------------------------------|
  // | EXISTENCE CHECK
  // |------------------------------------------------------------------------|

     if (!infile.is_open())
       {
         cout<<" "<<endl;
         cout << "I'm sorry Dave... I can't do that..."<<endl;;
         exit(1);
       }


  // |------------------------------------------------------------------------|
  // | REMOVE FIRST TWO LINES OF ADMB CORRELATION FILE
  // |------------------------------------------------------------------------|

     infile.getline(header_cor,maxline,'\n'); //removes first line
     infile.getline(header_cor,maxline,'\n'); //removes second line


  // |------------------------------------------------------------------------|
  // | READ IN PARAMETER VALUE, ST. DEV., AND CORR VALUES
  // |------------------------------------------------------------------------|


     while (!infile.eof())
      {
        infile >> index;
        infile >> pname;
        infile >> estimate;
        infile >> stdev;

        // cout<<"COR file index "<<index<<", parameter name "<<pname<<endl;

        for (icount=1; icount<=index; icount++)
        {
            infile >> tmp_cor[icount];

            ii=index;  //-(Boot_index-1);
            jj=icount; //-(Boot_index-1);

            Boot_stdev_init(ii) = stdev;
            Boot_index_init(ii) = estimate;

            Boot_corr(ii,jj) = tmp_cor[icount];

        }
      }

     infile.close();

     cout<<" "<<endl;
     cout<<"I say take off and nuke it from orbit - only way to be sure"<<endl;


  // |------------------------------------------------------------------------|
  // | TRIM AND MIRROR CORRELATION MATRIX
  // |------------------------------------------------------------------------|

    for(ii=1; ii<=maxpar; ii++)
      {
         Boot_stdev(ii) = Boot_stdev_init(ii);
         Boot_index(ii) = Boot_index_init(ii);

       for(jj=1; jj<=maxpar; jj++)
         {
           Boot_correlation(ii,jj) = Boot_corr(ii,jj);
         }
      }

     for (ii=1; ii<=(maxpar-1); ii++)
       {
        for (jj=(ii+1); jj<=maxpar; jj++)
          {
            Boot_correlation(ii,jj) = Boot_correlation(jj,ii);
          }
       }

     for (ii=1; ii<=maxpar; ii++)
       {
        for (jj=1; jj<=maxpar; jj++)
          {
            Boot_covariance(ii,jj) = Boot_correlation(ii,jj)*Boot_stdev(ii)*Boot_stdev(jj);
          }
       }

     cout<<" "<<endl;
     cout<<"Trim, correlation and covariance implemented."<<endl;


  // |------------------------------------------------------------------------|
  // | CHOLESKI DECOMPOSITION FOR SQUARE ROOT OF CORRELATION
  // |------------------------------------------------------------------------|

     Boot_choleski = choleski_decomp(Boot_covariance);

     cout<<" "<<endl;
     cout<<"Choleski decomp successful; launch sequence initiated"<<endl;


  // |------------------------------------------------------------------------|
  // | MESSAGE OUTPUT
  // |------------------------------------------------------------------------|

     cout<<""<<endl;
     cout<<"M-5 tie in... M-5... working... calling Bootstrap"<<endl;
     cout<<""<<endl;


  // |------------------------------------------------------------------------|
  // | REMOVE OLD BOOTSTRAP RUNS
  // |------------------------------------------------------------------------|

     remove("boot.txt");
     remove("mr.txt");
     remove("N.txt");
     remove("SB.txt");
     remove("cpue1.txt");
     remove("cpue2.txt");
     remove("cpue_fish.txt");
     remove("agef_fish.txt");
     remove("agem_fish.txt");
     remove("agef_srv.txt");
     remove("agem_srv.txt");
     remove("Ff.txt");
     remove("Fm.txt");
     remove("Nf.txt");
     remove("Nm.txt");
     remove("Natagef.txt");
     remove("Natagem.txt");
     remove("mr_f.txt");
     remove("mr_m.txt");


  // |------------------------------------------------------------------------|
  // | MAIN COUNTER
  // |------------------------------------------------------------------------|
  
     void (" 'nboot' itself is located in the control file");

      nwork = 0;             // number of reasonable bootstrap vectors
      ndump = 0;             // number of unreasonable bootstrap vectors

     void (" You get points for better progression quotes ");
     
  // |------------------------------------------------------------------------|



  // |------------------------------------------------------------------------|
  // | BOOTSTRAP LOOP
  // |------------------------------------------------------------------------|
  // |
  // | - Note that the random draws are implemented INSIDE the
  // |   bootstrap loop, and that each call produces a vector
  // |   MAXPAR long which is replicated over each iteration
  // |   of 'nboot'


      while(nwork < nboot)   //Primary loop -
        {                           //Doesn't close
         if(nwork<=nboot)              //until the very end
        {

      draw.initialize();
      draw.fill_randn(rng);         //Random draws ~N(0,1)


      //for (ii=1; ii<=maxpar; ii++)  //Scales N(0,1)
       //{
         //Bootshit(ii) = draw(ii);//  * Boot_stdev(ii);
       //}



  // |------------------------------------------------------------------------|
  // | MULIVARIATE NORMAL
  // |------------------------------------------------------------------------|
  // |
  // | This implements the multivariate normal distribution
  // | by including the correlation between parameters
  // | contained in the choleski matrix. If there were no correlation,
  // | the choleski would contain only zeroes in the off-diagonal,
  // | and condense to a univariate normal
  // |

      //Boot_ran = Boot_choleski * Bootshit;


  // |------------------------------------------------------------------------|
  // | ADD THE RANDOM STANDARD DEVIATION TO THE MLE PARAMETER
  // |------------------------------------------------------------------------|


     bootN.initialize();
     bootN = Boot_index +  Boot_choleski * draw;

     cout<<"My ducks are all in a row"<<endl;
     
  // |------------------------------------------------------------------------|
  // | PARAMETER VALUES
  // |------------------------------------------------------------------------|
  // |
  // | The vector 'bootN' contains all estimated MLE parameters
  // | from the original model, now with random variability added.
  // | This code assigns the elements of the 'bootN' vector their
  // | proper names as per the PARAMETER_SECTION above.
  // |
  // | That way we can simply re-call the primary model functions
  // | already defined above to calculate the new set of bootstrap
  // | derived quantities for each bootstrap iteration

     log_mean_rec.initialize();
     log_mean_y1.initialize();
     log_rec_dev.initialize();
     init_pop.initialize();
     log_surv_q1.initialize();
     log_surv_q2.initialize();
     log_comm_q.initialize();
     log_avg_F.initialize();            
     log_F_devsf.initialize();
     log_F_devsm.initialize();
     mF.initialize();


     log_mean_rec = bootN(1);                  // Mean recruitment;
     log_mean_y1  = bootN(2);                  // Mean Year 1 naa;

     for(ii = 3; ii<=42; ii++)
       {
         init_pop(ii-2) = bootN(ii);          // Year 1 abundance deviations
       }

     for(ii = 43; ii<=78; ii++)
       {
         log_rec_dev(ii+styr-43) = bootN(ii);   // Recruitment deviations
       }
 
     log_surv_q1   =  bootN(79);            // Catchability (survey CPUE 1)
     log_surv_q2   =  bootN(80);            // Catchability (survey CPUE 2)
     log_comm_q    =  bootN(81);            // Catchability (commercial CPUE)

     log_avg_F     =  bootN(82);            // Mean F


      
    for(ii = 83; ii<=118; ii++)
      {
        log_F_devsf(ii - 83+styr) = bootN(ii);      // F deviations, females
      }


    for(ii = 119; ii<=154; ii++)
      {
        log_F_devsm(ii-119 + styr) = bootN(ii);  // F deviations, males
      }


    for(ii = 155; ii<=161; ii++)
      {
        mF(ii-154) = bootN(ii);                // F levels for projected N
      }


     
     // |---------------------------------------------------------------------|
     // | LIMITS ON RANDOM DRAWS AS PER ORIGINAL PARAMETER BOUNDS
     // |---------------------------------------------------------------------|
     // |
     // | TEST DISTRIBUTION OF DRAWS AGAINST MLE DISTRIBUTION FROM MLE MEAN
     // | AND ST. DEV. FROM THE .STD FILE - BOUNDS SHOULD SERVE **ONLY** TO
     // | PROVIDE GENERAL LIMITS TO THE ADMB MINIMIZATION FUNCTION, **NOT**
     // | TO COMPENSATE FOR BAD CODE, ERRORS IN STATISTICAL ASSUMPTIONS, OR
     // | UNINFORMATIVE DATA.
     // |
     // | EXPERIMENT WITH PENALIZED DEVIATIONS FROM PRIOR VALUES AND DIFFERENT
     // | ERROR DISTRIBUTIONS IN THE OBJECTIVE FUNCION IF HARD BOUNDS RESULT
     // | IN A TRUNCATED EXPLORATION OF THE PARAMETER SPACE, AND THEN ALTER THE
     // | BOUNDS ACCORDINGLY
     // |
     // |
    
     if (log_mean_rec     < 0.0   ||
         log_mean_rec     > 10.0  ||
         log_mean_y1      < 0.0   ||
         log_mean_y1      > 8.0  ||
         min(log_rec_dev) < -5.0 ||
         max(log_rec_dev) > 5.0  ||
         min(init_pop)    < -5.0 ||
         max(init_pop)    > 5.0  ||
         log_surv_q1           < -20.0 ||
         log_surv_q1         > 0.0   ||
         log_surv_q2           < -20.0 ||
         log_surv_q2           > 0.0   ||
         log_comm_q          < -20.0 ||
         log_comm_q           > 0.0   ||
         log_avg_F        < -5   ||
         log_avg_F        > -2    )

     {
         //cout<<"Dumping vector"<<endl;
          //go on to the next vector
         ndump++;
         continue;
     }
     else
     {
         nwork++;
         cout<<"Got vector "<<nwork<<endl;
     }


  // |-----------------------------------------------------------|
  // | CALL MODEL FUNCTIONS FROM MLE STRUCTURE ABOVE
  // |-----------------------------------------------------------|
  // |
  // | As we are not estimating anything in the bootstrap, we simply
  // | call the population dynamics functions from the MLE model
  // | but without the objective function. The randomized MLE
  // | parameters are used to calculate the derived quantities
  // | (abundance, spawning biomass, etc.) which are written to
  // | external files for assessment of variance.
  
	
     Mortality();
     Abundance();
     Catch();
     Predicted();
     Population();
     spr();	
     Get_ABC();

  // |-----------------------------------------------------------|
  // | SAID EXTERNAL WRITING
  // |-----------------------------------------------------------|


    ofstream ofs0("boot.txt",ios::app);
    {
        ofs0<<bootN<<endl;
    }

   ofstream ofs1("mr.txt",ios::app);
    {
        ofs1<<pred_MR_abund<<endl;
    }

   ofstream ofs2("N.txt",ios::app);
    {
        ofs2<<tot_n<<endl;
    }

   ofstream ofs3("SB.txt",ios::app);
    {
        ofs3<<spawn_biom<<endl;
    }

   ofstream ofs4("cpue1.txt",ios::app);
    {
        ofs4<<pred_cpue1<<endl;
    }

   ofstream ofs5("cpue2.txt",ios::app);
    {
        ofs5<<pred_cpue2<<endl;
    }

   ofstream ofs6("cpue_fish.txt",ios::app);
    {
        ofs6<<pred_fshy_cpue<<endl;
    }

   ofstream ofs7("agef_fish.txt",ios::app);
    {
        ofs7<<eac_fishf<<endl;
    }

   ofstream ofs8("agem_fish.txt",ios::app);
    {
        ofs8<<eac_fishm<<endl;
    }

   ofstream ofs9("agef_srv.txt",ios::app);
    {
        ofs9<<eac_srv1f<<endl;
    }

   ofstream ofs10("agem_srv.txt",ios::app);
    {
        ofs10<<eac_srv1m<<endl;
    }

   ofstream ofs11("Ff.txt",ios::app);
    {
        ofs11<<Fmortf<<endl;
    }

   ofstream ofs12("Fm.txt",ios::app);
    {
        ofs12<<Fmortm<<endl;
    }

   ofstream ofs13("Nf.txt",ios::app);
    {
        ofs13<<tot_nf<<endl;
    }

   ofstream ofs14("Nm.txt",ios::app);
    {
        ofs14<<tot_nm<<endl;
    }

   ofstream ofs15("nategef.txt",ios::app);
    {
        ofs15<<natagef<<endl;
    }

   ofstream ofs16("nategem.txt",ios::app);
    {
        ofs16<<natagem<<endl;
    }

 

     
   
     //} //if-else loop for 'good_count = 1'
     } // close conditional 'if' g++ loop
     } // close primary g++ loop


 


    cout<<""<<endl;
    cout<<""<<endl;
    cout << "** "<<nboot<<" BOOTSTRAP ITERATIONS ARE COMPLETE **" << endl;
    cout << "** "<<ndump<<" bootstrap vectors dumped **" << endl;
    cout << "** "<<ndump+nboot<<" total iterations completed**" << endl;
    cout << "** "<<(1.0*(ndump*100/(ndump+nboot+oo)))<<"% of the total parameter space ignored"<< endl;
    cout<<""<<endl;

    cout<< "O frabjous day!"<<endl;
    cout<< "The sheep frolic!"<<endl;
    cout<<""<<endl;
    cout<<"        ...moo..."<<endl;
    cout<<"            | "<<endl;
    cout<<"            | "<<endl;
    cout<<"            | "<<endl;
    cout<<"             _.%%%%%%%%%%%%%             " <<endl;
    cout<<"            //-_%%%%%%%%%%%%%            " <<endl;
    cout<<"           (_ %\\%%%%%%%%%%%%%%~             "<<endl;
    cout<<"               %%%%%%%%%%%%%%             "<<endl;
    cout<<"                 %%%%%*%%%%              "<<endl;
    cout<<"            ,,,,,,||,,,,||,,,,,          "<<endl;
    cout<<""<<endl;


GLOBALS_SECTION
	#include <admodel.h>
	#include <string.h>
	#include <time.h>

	#undef REPORT
	#define REPORT(object) report << #object "\n" << setw(8) \
	<< setprecision(4) << setfixed() << object << endl;

	#undef COUT
	#define COUT(object) cout << #object "\n" << setw(6) \
	<< setprecision(3) << setfixed() << object << endl;

	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;

	adstring stripExtension(adstring fileName)
	{
		/*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
		*/
		const int length = fileName.size();
		for (int i=length; i>=0; --i)
		{
			if (fileName(i)=='.')
			{
				return fileName(1,i-1);
			}
		}
		return fileName;
	}
	

REPORT_SECTION
 report<<"Report for the "<<BaseFileName<<" sablefish age-structured model"<<endl;
  report<<""<<endl;
    //likelihoods
  //REPORT(yy);
  //REPORT(aa);
  REPORT(catch_like);
  REPORT(age_like);
  REPORT(surv_like);
  REPORT(penalties);
  report<<""<<endl;
  report<<""<<endl;
  REPORT(ABC)
  report<<""<<endl;
  REPORT(obs_MR_abund)
  REPORT(pred_MR_abund)
  report<<""<<endl;
  REPORT(obs_MR_abund_SRV)
  REPORT(pred_MR_abund_SRV)
  report<<""<<endl;
  REPORT(natagef);
  REPORT(natagem);
  REPORT(spawn_biom);


 //report<<"Report for the "<<BaseFileName<<" sablefish age-structured model"<<endl;
//report<<""<<endl;
//
//report<<"Number parameters estimated"<<endl;
//report<<initial_params::nvarcalc()<<endl;
//report<<" "<<endl;

//report<<""<<endl;
//report << "obj_fun "<< obj_fun <<endl;
//report << "//" <<endl;
//report<<"catch_like"<<endl;
//report<< catch_like<<endl;
//report<<"age_like"<<endl;
//report<< age_like<<endl;
//report<<"surv_like"<<endl;
//report<< surv_like<<endl;
//report<<"penalties"<<endl;
//report<< penalties<<endl;
//report<<" "<<endl;

//report <<"Nf_proj"<<endl;
//report<<natage_projf<<endl;
////report <<"Nm_proj"<<endl;
//report<<natage_projm<<endl;
//report<<""<<endl;
//report <<"ABC: "<<ABC<<endl;
//report << "Exploited abundance projected" << " " <<expl_abund_proj<<endl;
//report << "Exploited biomass projected" << " " <<expl_biom_proj<<endl;
//report << "//" <<endl;
//////////////////////////////////
//report << "Year////////" << yrs_MR//<<endl;
////report << "Predicted MR abund// " << pred_MR_abund<<endl;
//report << "Observed MR abund////" << obs_MR_abund//<<endl<< endl;
//report <<""<<endl;
//////report << "Predicted MR abund SRV//" << pred_MR_abund_SRV<<endl;
//report << "Observed MR abund SRV////" << obs_MR_abund_SRV//<<endl<< endl;
//report <<""<<endl;
//
//report<<"Mean Rec"<<endl;
//report<< avg_rec<<endl;
//report<<" "<<endl;

//report<<"mF"<<endl;
//report<<mF<<endl;
 
//report<<"obs_cpue1"<<endl;
//report<< obs_cpue1_biom<<endl;
//report<<"pred_cpue1"<<endl;
//report<< pred_cpue1<<endl;

//report<<"obs_cpue2"<<endl;
//report<< obs_cpue2_biom<<endl;
//report<<"pred_cpue2"<<endl;
//report<< pred_cpue2<<endl;
//report<<""<<endl;

//report<<"obs_fshy_cpue"<<endl;
//report<< obs_fshy_cpue<<endl;
////report<<"pred_fshy_cpue"<<endl;
//report<< pred_fshy_cpue<<endl;
//report<<""<<endl;


//report<<"natagef"<<endl; 
//report<< natagef<<endl;
//report<<""<<endl;
////report<<"spawn_biom"<<endl; 
//report<< spawn_biom<<endl;
////report<<"Ff"<<endl; 
//report<< Ff<<endl;
//////report<<"Fm"<<endl; 
//report<< Fm<<endl;
//report<<""<<endl;

//report<<"proj_F_catch"<<endl;
//report<<(elem_prod(wt_f_proj,elem_prod(elem_div(Fself, Z1),
////////////////////////////////////// elem_prod(1-S1,(mr_f_proj)))))<<endl;
//report<<"proj_M_catch"<<endl;
//report<<(elem_prod(wt_m_proj,elem_prod(elem_div(Fselm, Z2),
//////////////////////////////////////yy elem_prod(1-S2,(mr_m_proj)))))<<endl;









//FUNCTION Selectivity

  // |-----------------------------------------------------------------------|
  // | FISHERY SELECTIVITIES - keep for future versions
  // |-----------------------------------------------------------------------|
  // | It is assumed that selectivities for the commercial longline fishery
  // | and the ADF&G longline survey are different
  // |


    // |-----------------------------------------------------|
    // | Survey
    // |-----------------------------------------------------|
    
    //for (j=1;j<=nages;j++)
      //{
      //  srv1_sel_m(j) = 1/(1+mfexp(-mfexp(srv1_sel_slope_m)*(j-mfexp(srv1_sel_a50_m))));
      //  srv1_sel_f(j) = 1/(1+mfexp(-mfexp(srv1_sel_slope_f)*(j-mfexp(srv1_sel_a50_f))));
      //}

    //Scale to 1
  
    //srv1_sel_m   = srv1_sel_m / max(srv1_sel_m);         
    //srv1_sel_f   = srv1_sel_f / max(srv1_sel_f);


    // |-----------------------------------------------------|
    // | Fishery
    // |-----------------------------------------------------|

    //for (j=1;j<=nages;j++)
      //{ 
        //fish_sel_m(j) = 1/(1+mfexp(-mfexp(fish_sel_slope_m)*(j-mfexp(fish_sel_a50_m))));
        //fish_sel_f(j) = 1/(1+mfexp(-mfexp(fish_sel_slope_f)*(j-mfexp(fish_sel_a50_f))));
      //}

    //Scale to 1
  
    //fish_sel_m   = fish_sel_m / max(fish_sel_m);       
    //fish_sel_f   = fish_sel_f / max(fish_sel_f);
 
