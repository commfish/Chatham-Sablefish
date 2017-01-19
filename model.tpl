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


 // |---------------------------------------------------------------------------|
 // | DESIGN MATRIX FOR PARAMETER CONTROLS                                      |
 // |---------------------------------------------------------------------------|
 // | - theta_DM -> theta is a vector of estimated parameters.

    init_int    n_theta;
    init_matrix theta_DM(1,n_theta,1,7);
    vector      theta_ival(1,n_theta);
    vector      theta_lb(1,n_theta);
    vector      theta_ub(1,n_theta);
    ivector     theta_phz(1,n_theta);
    ivector     theta_iprior(1,n_theta);
    vector      theta_p1(1,n_theta);
    vector      theta_p2(1,n_theta);
    
    !! theta_ival   = column(theta_DM,1);
    !! theta_lb     = column(theta_DM,2);
    !! theta_ub     = column(theta_DM,3);
    !! theta_phz    = ivector(column(theta_DM,4));
    !! theta_iprior = ivector(column(theta_DM,5));
    !! theta_p1     = column(theta_DM,6);
    !! theta_p2     = column(theta_DM,7);



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
   number M;
   number Mf;
   number srv_slp;
   number fsh_slp;
   
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
          M           = dMiscCont(4);
          Mf          = dMiscCont(5);
          srv_slp     = dMiscCont(6);
          fsh_slp     = dMiscCont(7);
          ph_rec      = dMiscCont(8);
          ph_init     = dMiscCont(9);
          ph_F        = dMiscCont(10);
          ph_spr      = dMiscCont(11);


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
  // | 7) spawn_fract -> spawning month (set to February [2])
  // | 8) srv_fract   -> survey month (set to August [8])
  // | 9) fsh_fract   -> month when fishery takes place (set to October [10])
  // | 10) mr_fract   -> month when fishery mark-recapture occurs
  // |                   NOTE - this is necessary for the M-R in the unified model
  // |                          because the M-R abundance is the MEAN abundance
  // |                          for the fishery (i.e. halfway through). So the
  // |                          mr_fract is 6 weeks later than fsh_fract
               
     init_int      nages
     init_int      styr
     init_int      endyr
     init_int      recage

     init_number   spawn_fract
     init_number   srv_fract
     init_number   fshy_fract
     init_number   mr_fract

     vector        yy(styr,endyr)
     vector        aa(1,nages)


  // |--------------------------------------------------------------------------|
  // | Gear selectivity for fishery and survey (using NOAA values)
  // |--------------------------------------------------------------------------|
  // | Selectivity curves from Dana Hanselman @ NOAA
  
     init_vector	vfish_sel_f(1,nages)
     init_vector	vfish_sel_m(1,nages)
     init_vector	vsrv_sel_f(1,nages)
     init_vector	vsrv_sel_m(1,nages)

  
  // |--------------------------------------------------------------------------|
  // | PHYSIOLOGY                                                               
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) p_mature       -> proportion of females mature-at-age
  // | 2) wt_avg_fshy    -> mean weight-at-age in commercial fishery
  // | 3) wt_avg_srv     -> mean weight-at-age in longline survey
  // | 4) wt_avg_fshy_f  -> mean weight-at-age in commercial fishery, females
  // | 5) wt_avg_fshy_m  -> mean weight-at-age in commercial fishery, males
  // | 6) wt_avg_srv_f   -> mean weight-at-age in longline survey, females
  // | 7) wt_avg_srv_m   -> mean weight-at-age in longline survey, males

     init_vector   p_mature(1,nages)
     init_vector   wt_avg_fshy(1,nages)
     init_vector   wt_avg_srv(1,nages)
     init_vector   wt_avg_fshy_f(1,nages)
     init_vector   wt_avg_fshy_m(1,nages)
     init_vector   wt_avg_srv_f(1,nages)
     init_vector   wt_avg_srv_m(1,nages)


  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL LONGLINE COMMERCIAL CATCH                                   
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) obs_catch     -> total annual catch

     init_vector   obs_catch(styr,endyr)


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
  // |
  // | 5) obs_MR_abund  -> mark-recapture estimates of abundance
  // | 6) obs_MR_var    -> variance of MR estimates

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
  // | 3) nsamples_fish_age   -> sample size
  // | 4) oac_fish            -> age comp data
  // |
  // | Relative sample size = square root of number of fish aged to scale data impact

     init_int      nyrs_fish_age 
     init_vector   yrs_fish_age(1,nyrs_fish_age) 
     init_vector   nsamples_fish_age(1,nyrs_fish_age)   
     init_matrix   oac_fish(1,nyrs_fish_age,1,nages) 



  // |--------------------------------------------------------------------------|
  // | LONGLINE SURVEY AGE COMPOSITION                                          
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) nyrs_srv_age       -> number of years of survey age comps
  // | 2) yrs_srv_age        -> specific years age comps
  // | 3) nsamples_srv_age   -> sample size FEMALES
  // | 4) oac_srv            -> female age comp data
  // |
  // | Relative sample size = square root of number of fish aged to scale data impact

     init_int      nyrs_srv_age 	
     init_ivector  yrs_srv_age(1,nyrs_srv_age)   
     init_vector   nsamples_srv_age(1,nyrs_srv_age) 
     init_matrix   oac_srv(1,nyrs_srv_age,1,nages)  
 
  
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

     vector offset(1,2);       
  
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
     srv_fract   = (srv_fract - 1)   / 12;	// Fraction of year at which survey occurs
     fshy_fract  = (fshy_fract - 1)  / 12;	// Fraction of year at which fishery occurs
     mr_fract    = (mr_fract - 1)    / 12;      // Fraction of year for mark-recapture fit
                                                //   for unified m-r estimates (mean
                                                //   abundance half-way through fishery)



    offset.initialize();


    for (i=1; i<=nyrs_fish_age; i++)
      {
        oac_fish(i)/=sum(oac_fish(i));
        offset(1) -= nsamples_fish_age(i)*((oac_fish(i) + 0.00001)*log(oac_fish(i) + 0.00001));
      }

     for (i=1; i<=nyrs_srv_age; i++)
       {
         oac_srv(i)/=sum(oac_srv(i));
         offset(2) -= nsamples_srv_age(i)*((oac_srv(i) + 0.00001)*log(oac_srv(i) + 0.00001));

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
  // | - log_mean_rec  <- mean recruitment
  // | - log_rec_dev   <- recruitment deviation vector   

     init_bounded_vector   init_pop(1,nages-1,-5,5,ph_init);
     init_bounded_vector   log_rec_dev(styr,endyr,-5,5,ph_rec);

     number mort;


  // |---------------------------------------------------------------------------------|
  // | FISHING MORTALITY PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | - log_F_devs       <- annual fishing mortality deviations
  // |
  // | - Fmort            <- F * selectivity
  // |
  // | - fish_sel_a50     <- age of 50% selectivity fishery
  // |
  // | - vector for fishery selectivity 
  // | - vector for selectivity selectivity

     init_bounded_vector  log_F_devs(styr,endyr,-10.,10.,ph_F);
     vector               Fmort(styr,endyr);
     
     init_bounded_number fish_sel_a50(0,3,ph_init);

     vector	fish_sel(1,nages)
     vector	srv_sel(1,nages)
    

  // |---------------------------------------------------------------------------------|
  // | POPULATION VECTORS
  // |---------------------------------------------------------------------------------|
  // | - tot_n          -> total abundance
  // | - tot_biomass    -> total biomass
  // | - spawn_biomass  -> spawning biomass
  // | - recruitment    -> predicted recruitment
  // |
  // | - pred_catch     -> total annual catch
	
     sdreport_vector  tot_n(styr,endyr)
     sdreport_vector  tot_biomass(styr,endyr)
     sdreport_vector  spawn_biom(styr,endyr);
     sdreport_vector  recruitment(styr,endyr);
     
     vector           pred_catch(styr,endyr);

     vector           wt_mature_s(1,nages);
     vector           wt_mature_f(1,nages);
     vector           wt_fsh_sel(1,nages);
   

  // |---------------------------------------------------------------------------------|
  // | POPULATION MATRICES
  // |---------------------------------------------------------------------------------|
  // | - natage         -> abundance-at-age
  // | - catege         -> catch-at-age
  // |
  // | - Z              -> total mortality-at-age
  // | - F              -> fishing mortality-at-age
  // | - S              -> survival-at-age
  // | 
  // | - eac_fish       -> proportion-at-age fishery 
  // | - eac_srv        -> proportion-at-age survey 

     matrix  natage(styr,endyr,1,nages);
     matrix  catage(styr,endyr,1,nages);
     
     matrix  Z(styr,endyr,1,nages); 
     matrix  F(styr,endyr,1,nages); 
     matrix  S(styr,endyr,1,nages);

     matrix  eac_fish(1,nyrs_fish_age,1,nages); 
     matrix  eac_srv(1,nyrs_srv_age,1,nages); 


  // |---------------------------------------------------------------------------------|
  // | STATISTICAL VECTORS
  // |---------------------------------------------------------------------------------|
  // | - effn_fish_age       -> effective numbers-at-age fishery 
  // | - effn_fish_age       -> effective numbers-at-age fishery 
  // | - sdnr_fish_age       -> SDNR numbers-at-age fishery 
  // | - sdnr_fish_age       -> SDNR numbers-at-age fishery

     vector effn_fish_age(1,nyrs_fish_age)
     vector sdnr_fish_age(1,nyrs_fish_age)

     vector effn_srv_age(1,nyrs_srv_age)
     vector sdnr_srv_age(1,nyrs_srv_age)

     
  // |---------------------------------------------------------------------------------|
  // | MARK-RECAPTURE AND CPUE
  // |---------------------------------------------------------------------------------|
  // | - pred_cpue1          -> Survey CPUE (both sexes) 1988 - 1996
  // | - pred_cpue2          -> Survey CPUE (both sexes) 1997 - present
  // |
  // | - pred_fshy_pue       -> Fishery CPUE (both sexes)
  // |
  // | - pred_MR_abund       -> MR abundance
  // | - pred_MR_abund_SRV   -> MR abundance survey
  // |
  // | - s_srv               -> fraction surviving to time of survey
  // | - s_fshy              -> surviving to time of fishery
  // | - s_mr                -> fraction surviving to time of mark-recapture
  // | - s_spn               -> fraction surviving to spawning
  // |
  // | - q_cpue1             -> q for survey 1 catchability 
  // | - q_cpue2             -> q for survey 2 catchability 
  // | - q_fshy              -> q for fishery catchability both sexes
  // |
  // | NOTE - the presence of 's_mr' differs from the sex-specified model because in
  // |        in this case, we estimate the mark-recapture as a simple mean abundance
  // |        which occurs in the middle of the fishery, as opposed to the abundance
  // |        at the very beginning of the fishery. So we decrement survival until
  // |        the middle of the commercial fishery.

     vector  pred_cpue1(1,nyrs_cpue1); 
     vector  pred_cpue2(1,nyrs_cpue2);

     vector  pred_fshy_cpue(1,nyrs_fshy_cpue);

     sdreport_vector  pred_MR_abund(1,nyrs_MR);
     sdreport_vector  pred_MR_abund_SRV(1,nyrs_MR);

     matrix  mr_f_cpue(1,nyrs_fshy_cpue,1,nages-2);
     matrix  mr_m_cpue(1,nyrs_fshy_cpue,1,nages-2);

     number  s_srv; 
     number  s_fshy; 
     number  s_mr;  
     number  s_spn;  

     sdreport_number  q_cpue1;
     sdreport_number  q_cpue2;
     sdreport_number  q_fshy;

  
  // |--------------------------------------------------------------------------|
  // | SPR / ABC / PROJECTIONS                               
  // |--------------------------------------------------------------------------|
  // | mF                 -> vector of potential F levels
  // | SBx                -> spawning biomass from mF            
  // | Bx                 -> biomass from mF
  // |
  // | Nspr               -> Number of projected spawners
  // | N0                 -> Number projected with F = 0 

     init_bounded_vector mF(1,9,0.01,1,ph_spr)
     vector SBx(1,10)
     vector Bx(1,10)
					       
     matrix Nspr(1,10,1,nages)			
     matrix N0(1,10,1,nages)
	
     matrix N_ABC(endyr+1,endyr+15,1,nages);
     matrix N_OFL(endyr+1,endyr+15,1,nages);
     
     vector FABC(1,nages);
     vector FOFL(1,nages);
     
     vector ZABC(1,nages);
     vector ZOFL(1,nages);
     
     matrix catage_ABC(endyr+1,endyr+15,1,nages);
     matrix catage_OFL(endyr+1,endyr+15,1,nages);
     
     vector catch_ABC(endyr+1,endyr+15);
     vector catch_OFL(endyr+1,endyr+15);

     sdreport_vector spawn_biom_ABC(endyr+1,endyr+15);
     sdreport_vector spawn_biom_OFL(endyr+1,endyr+15);

     sdreport_number B40;
     sdreport_number ABC;
     sdreport_number OFL;
     number F55;
     number F50;
     number F45;
     number F40;

     number stdev_rec;


     

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
     vector  penalties(1,5);
     vector  surv_like(1,5);
     vector  age_like(1,2);
     number  rec_like;
     number  catch_like;
     number  F_mort_regularity;

     number  Like;
     objective_function_value obj_fun;	

						
 LOCAL_CALCS

   wt_mature_f = elem_prod(p_mature/2,wt_avg_fshy);
   wt_mature_s = elem_prod(p_mature/2,wt_avg_srv);

   //wt_fsh_sel = elem_prod(fish_sel,wt_avg_fshy);
   //wt_fsh_sel_m = elem_prod(fish_sel,wt_avg_fshy);


  // |-------------------------------------------------------------------------|
  // | Survival fractions to decrement abundance relative to annual timing
  // |-------------------------------------------------------------------------|
  // |
  // | This will make a difference if distinct natural mortalities are used
  // | between males and females. If both are input at 0.1, then one decrement
  // | for the fishery and another for the survey is all that is required

     s_spn  = mfexp(-1*spawn_fract*M);
     s_srv  = mfexp(-1*srv_fract*M);
     s_fshy  = mfexp(-1*fshy_fract*M);
     s_mr    = mfexp(-1*(mr_fract)*M);



 END_CALCS


PROCEDURE_SECTION
     initializeModelParameters();
     Selectivity();
     Mortality();
            if(DEBUG_FLAG == 1) cout<<"**Mortality**"<<endl;
     Abundance();	
            if(DEBUG_FLAG == 1) cout<<"**Abundance**"<<endl;
     Catch();	
            if(DEBUG_FLAG == 1) cout<<"**Catch**"<<endl;
     Predicted();
            if(DEBUG_FLAG == 1) cout<<"**Predicted**"<<endl;
     Population();
            if(DEBUG_FLAG == 1) cout<<"**Population**"<<endl;
		
     if (last_phase())
       {
         spr();
     	    if(DEBUG_FLAG == 1) cout<<"**spr**"<<endl;
                
         Projection();
            if(DEBUG_FLAG == 1) cout<<"**Projection**"<<endl;

       }

     Objective_Function();


     if(mceval_phase())	
       {
         writePosteriorSamples();

         // |----------------------------------------------
         // | This can either be done here with the
         // | 'evalout' file or below with the calls inside
         // | the "writePosteriorSamples' function.
         // | I prefer using the separate function because
         // | the files are distinct and there is no need
         // | to parse anything in R afterwards
         
         //evalout<<tot_n<<" "<<
         //     pred_MR_abund<<""<<
         //     pred_MR_abund_SRV<<""<<
         //     spawn_biom<<" "<<
         //     recruitment<<" "<<
         //     endl;
       }

FUNCTION void writePosteriorSamples()
	/**
	- This function is only envoked when the -mceval
		command line option is implemented.
	*/
	ofstream ofs0("spbm.ps",ios::app);
	ofs0<<spawn_biom<<endl;

	ofstream ofs1("MR.ps",ios::app);
	ofs1<<pred_MR_abund<<endl;

	ofstream ofs2("MR_srv.ps",ios::app);
	ofs2<<pred_MR_abund_SRV<<endl;

	ofstream ofs3("recruit.ps",ios::app);
	ofs3<<recruitment<<endl;

	ofstream ofs4("tot_n.ps",ios::app);
	ofs4<<tot_n<<endl;

        

FUNCTION void initializeModelParameters()
	//fpen = 0;

	log_mean_rec          = theta(1);
	log_mean_y1           = theta(2);
	log_avg_F             = theta(3);
        log_comm_q            = theta(4);
        log_surv_q1           = theta(5);
        log_surv_q2           = theta(6);

FUNCTION Selectivity

  // |-----------------------------------------------------------------------|
  // | FISHERY SELECTIVITIES - keep for future versions
  // |-----------------------------------------------------------------------|
  // | It is assumed that selectivities for the commercial longline fishery
  // | and the ADF&G longline survey are different. The gear is identical
  // | but behavior differs.
  // |

    // |-----------------------------------------------------|
    // | Survey
    // |-----------------------------------------------------|
    
    for (j=1;j<=nages;j++)
      {
        //srv_sel(j) = 1/(1+mfexp(-srv_slp*(j-mfexp(srv_sel_a50))));
        srv_sel(j) = (vsrv_sel_f(j) + vsrv_sel_m(j))/2;
      }


   //  |-----------------------------------------------------|
   //  | Fishery
   //  |-----------------------------------------------------|

    for (j=1;j<=nages;j++)
      { 
        fish_sel(j) =  1/(1+mfexp(-fsh_slp*(j-mfexp(fish_sel_a50))));
        //fish_sel(j) =  (vfish_sel_f(j) + vfish_sel_m(j))/2;
      }



FUNCTION Mortality
  Fmort.initialize();
  F.initialize();
  Z.initialize();
  S.initialize();

  // |----------------------------------------------------------------------|
  // | ANNUAL DEVIATIONS FOR COMMERCIAL FISHERY
  // |----------------------------------------------------------------------|
  
     Fmort        = mfexp(log_avg_F +  log_F_devs);


  // |----------------------------------------------|
  // | Fishing Mortality matrices
  // |----------------------------------------------|
  // | Note here different mortalities by gender
  // | due to selectivity - input natural mortality
  // | is the same for both as per Johnson and Quinn
  // | (1988) and Hanselman et al. (2014)
  
     for (i=styr; i<=endyr; i++)
       {
             F(i) = (Fmort(i) * fish_sel);
       }		


  // |----------------------------------------------|
  // | Total Mortality & Survival
  // |----------------------------------------------|
     
     Z  = F + M;	
     S  = mfexp(-1.0*Z);
  

FUNCTION Abundance
  natage.initialize();
  tot_n.initialize();

  // |-------------------------------------------------------------------------|
  // | YEAR 1
  // |-------------------------------------------------------------------------|
  // | Equal proportions at recruitment are assumed
  // |
  // | Recruitment variability set to 1.2 as per Sigler et al. (2002) and
  // | Hanselman et al. (2014)

     natage(styr,1)= mfexp(log_rec_dev(styr) + log_mean_rec  
                            + (sigr*sigr)/2)/2;


     for(j=1;j<=nages-1;j++)
       {
         natage(styr,j+1)= mfexp(log_mean_y1 + init_pop(j)  
                                  + (sig1*sig1)/2)/2;
       }
        

  // |-------------------------------------------------------------------------|
  // | Subsequent years
  // |-------------------------------------------------------------------------|
  // | The ++ notation increments index range by 1 so that computations
  // |  for (1,nages-1) correspond to index range (2,nages)

     for (i=styr+1; i<=endyr; i++)
       {
         natage(i,1) = mfexp(log_rec_dev(i) + log_mean_rec + (sigr*sigr)/2)/2;   

       for(j=2; j<=nages; j++)
         {
           natage(i,j)  = natage(i-1,j-1)*S(i-1,j-1);
         }                                                                         
										  
       natage(i,nages)    += natage(i,nages)*S(i,nages);
   
       }


       
FUNCTION Catch
  catage.initialize();
  pred_catch.initialize();

  // |----------------------------------------------------------------------|
  // | COMMERCIAL CATCH-AT-AGE (note that age comps for OF are farther below)
  // |----------------------------------------------------------------------|
  // | Total annual catch from summed catch-at-age
       

     for (i=styr; i<=endyr; i++)
       {
         catage(i)     = elem_div(elem_prod(elem_prod(natage(i),F(i)),(1.-S(i))),Z(i));
         
         pred_catch(i)  = catage(i)*wt_avg_fshy;
       }

 
FUNCTION Predicted
  pred_cpue1.initialize();
  pred_cpue2.initialize();
  pred_fshy_cpue.initialize();
  pred_MR_abund.initialize();
  pred_MR_abund_SRV.initialize();
  eac_fish.initialize();
  eac_srv.initialize();

  // |-------------------------------------------------------------------------|
  // | Catchabilities
  // |-------------------------------------------------------------------------|
  // | Note that catchabilities are assumed *constant* between sexes
  // | Differences in retention are accounted as selectivity

     surv_q1 = mfexp(log_surv_q1);
     surv_q2 = mfexp(log_surv_q2); 
     comm_q  = mfexp(log_comm_q);

     wt_fsh_sel = elem_prod(fish_sel,wt_avg_fshy);


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

     dvar_matrix mr(1,nyrs_MR,1,nages);
     dvar_matrix mr_m(1,nyrs_MR,1,nages);
     dvar_matrix mr_SRV(1,nyrs_MR,1,nages);
     dvar_matrix mr_m_SRV(1,nyrs_MR,1,nages);

  
     for (i=1;i<=nyrs_MR;i++)
       {
         for(j=1; j<=nages; j++)
           {
             mr(i,j) = natage(yrs_MR(i),j) * fish_sel(j);
           }
       }


     for (i=1;i<=nyrs_MR;i++)
       {
         for(j=1; j<=nages; j++)
           {
             mr_SRV(i,j) = natage(yrs_MR(i),j) * srv_sel(j);
           }
       }


     for (i=1;i<=nyrs_MR;i++)
       {
         pred_MR_abund(i) = s_mr * 1000 * sum(mr(i));

         pred_MR_abund_SRV(i) = s_srv * 1000 * sum(mr_SRV(i));
                              
       } 
            
      
  // |-------------------------------------------------------------------------|
  // | ADF&G LONGLINE SURVEY CPUE 1988 - 1996
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance

     for (i=1;i<=nyrs_cpue1;i++)
       {
         pred_cpue1(i) = surv_q1 * sum(elem_prod(natage(yrs_cpue1(i)),
                                       srv_sel));  
       }


  // |-------------------------------------------------------------------------|
  // | ADF&G LONGLINE SURVEY CPUE 1997 - PRESENT
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance
					             
     for (i=1;i<=nyrs_cpue2;i++)
       {
         pred_cpue2(i) = surv_q2 * sum(elem_prod(natage(yrs_cpue2(i)),
                                       srv_sel));  
       }


  // |-------------------------------------------------------------------------|
  // | COMMERCIAL FISHERY CPUE
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance


     for (i=1;i<=nyrs_fshy_cpue;i++)
       {

         pred_fshy_cpue(i) =  comm_q * sum(elem_prod(natage(yrs_fshy_cpue(i)),wt_fsh_sel));
       }


  // |-----------------------------------------------------------------------|
  // | AGE COMPS, N, EFFN, SDNR
  // |-----------------------------------------------------------------------|


     //FISHERY
     for (i=1;i<=nyrs_fish_age;i++)
       {
         eac_fish(i)  = catage(yrs_fish_age(i)) / sum(catage(yrs_fish_age(i))) * ageage;     

         effn_fish_age(i) = (1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i));

                  sdnr_fish_age(i) = sdnr(eac_fish(i),oac_fish(i),
                            double(nsamples_fish_age(i)));
 
       }


     // SURVEY
     for (i=1;i<=nyrs_srv_age;i++)
       {
         eac_srv(i)  = elem_prod(srv_sel, natage(yrs_srv_age(i))*s_srv)
                         /sum(elem_prod(natage(yrs_srv_age(i))*s_srv,srv_sel)) * ageage;

         effn_srv_age(i) = (1-eac_srv(i))*eac_srv(i)/norm2(oac_srv(i)-eac_srv(i));

                  sdnr_srv_age(i) = sdnr(eac_srv(i),oac_srv(i),
                            double(nsamples_srv_age(i)));

       }


FUNCTION Population
   recruitment.initialize();
   tot_biomass.initialize();
   spawn_biom.initialize();
   tot_n.initialize();

  // |-------------------------------------------------------------------------|
  // | MAKE SURE UNITS ARE CORRECT
  // |-------------------------------------------------------------------------|
  // |*Note conversion from metric tonnes (for total biomass) to straight pounds
  // |  for Age 4 and Exploitable biomass calcs - this is for managers
  // |
  // |  - kg (wt-at-age)  x thousands (natage)  = metric tonnes
  // |  - 1 metric tonne = 2,204 pounds


     for (i=styr;i<=endyr;i++)
       {
         tot_n(i)       = sum(natage(i));
         
         recruitment(i) = natage(i,1);
    
         tot_biomass(i) = (natage(i)*wt_avg_srv); 
    	       
         spawn_biom(i)  =  natage(i) * mfexp(-spawn_fract*M) * wt_mature_s;
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
       int fpen = 1;
       if (last_phase())
         {
           fpen = 2;
         }


       penalties(3)  = dnorm(log_avg_F,log(0.09),0.5);
       penalties(4)  = dnorm(log_mean_rec,5,fpen);
       penalties(5)  = dnorm(log_mean_y1,6,fpen);
       

FUNCTION Catch_Likelihood

  catch_like.initialize();

  for(i=styr; i<=endyr; i++)
    {
      catch_like  +=  0.5*log(2*M_PI)  +log(sigma_catch) +  0.5*(square(log(obs_catch(i)+oo) -
                                     log(pred_catch(i)+oo))) / (2.*square(sigma_catch));
    }

FUNCTION Surv_Likelihood

  surv_like.initialize();

  // |----------------------------------------------------------------------|
  // | Mark-recapture estimates of abundance - normal
  // |----------------------------------------------------------------------|

    for (i=1; i<=nyrs_MR; i++)
      {
   
         surv_like(1) += 0.5*log(2*M_PI) + log(sqrt(obs_MR_var(i))) +
                    0.5*(square(obs_MR_abund(i)-pred_MR_abund(i)))
                   / (2*obs_MR_var(i));

         surv_like(2) += 0.5*log(2*M_PI) + log(sqrt(obs_MR_var_SRV(i))) +
                    0.5*(square(obs_MR_abund_SRV(i)-pred_MR_abund_SRV(i)))
                   / (2*obs_MR_var_SRV(i));
                   
      }



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
           age_like(1) -= nsamples_fish_age(i)*((oac_fish(i) + oo) *
                                                  log(eac_fish(i) + oo)) ;

         }

    // |--------------------------------------------------------------------|
    // | Longline survey age compositions
    // |--------------------------------------------------------------------|

       for (i=1; i <= nyrs_srv_age; i++)
         {
           age_like(2) -= nsamples_srv_age(i)*((oac_srv(i) + oo) *
                                                  log(eac_srv(i) + oo)) ;

         }


    age_like(1)   -= offset(1);
    age_like(2)   -= offset(2); 


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

     Like           = catch_like;	
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


  // |----------------------------------------------------------------------|
  // | RECRUITMENT (AGE 8)
  // |----------------------------------------------------------------------|

     for (i=1;i<=10;i++)
       {
         Nspr(i,1) = 1.;
         N0(i,1) = mean(recruitment(1982,endyr-recage));
       }


  // |----------------------------------------------------------------------|
  // | AGES 3 - 42+
  // |----------------------------------------------------------------------|

     for (j=2;j<=nages;j++)  // j loop
       {
         Nspr(1,j) = Nspr(1,j-1)*mfexp(-1.*M);
         N0(1,j)   = N0(1,j-1)*mfexp(-1.*M);

       for (i=2; i<=10; i++)
         {
           Nspr(i,j) = Nspr(i,j-1)*mfexp(-1.*(M+mF(i-1)*fish_sel(j-1)));
           N0(i,j)   = N0(i,j-1)*mfexp(-1.*(M+mF(i-1)*fish_sel(j-1)));
         }
       } // close 'j' loop


  // |----------------------------------------------------------------------|
  // | PLUS CLASS
  // |----------------------------------------------------------------------|

      Nspr(1,nages) += Nspr(1,nages)*mfexp(-1.*M);
      N0(1,nages)   += N0(1,nages)*mfexp(-1.*M);


     for(i=2; i<=10; i++)
      {
         Nspr(i,nages) += Nspr(i,nages)*mfexp(-1.*(M+mF(i-1)*fish_sel(nages)));
         N0(i,nages)   += N0(i,nages)*mfexp(-1.*(M+mF(i-1)*fish_sel(nages)));
       }

  // |----------------------------------------------------------------------|
  // | PROJECTED BIOMASS FOR F LEVELS
  // |----------------------------------------------------------------------|

     SBx(1) +=  Nspr(1)*wt_mature_s*mfexp(-spawn_fract*M);
     

     for (j=1;j<=nages;j++)
       {
         for(i=2; i<=10; i++)
           {
             SBx(i) += Nspr(i,j)*wt_mature_s(j)*
                       mfexp(-spawn_fract*(M+mF(i-1)*fish_sel(j)));
           }

         for(i=1; i<=10; i++)
           {
             Bx(i) +=  N0(i,j)*wt_mature_s(j);
           }
       } // close 'j' loop


  // |----------------------------------------------------------------------|
  // | CALCULATE F LEVELS
  // |----------------------------------------------------------------------|
  
     sprpen    = 100.*square(SBx(10)/SBx(1)-0.7);
     sprpen   += 100.*square(SBx(9)/SBx(1)-0.65);
     sprpen   += 100.*square(SBx(8)/SBx(1)-0.6);
     sprpen   += 100.*square(SBx(7)/SBx(1)-0.55);
     sprpen   += 100.*square(SBx(6)/SBx(1)-0.5);
     sprpen   += 100.*square(SBx(5)/SBx(1)-0.45);
     sprpen   += 100.*square(SBx(4)/SBx(1)-0.4);
     sprpen   += 100.*square(SBx(3)/SBx(1)-0.35);
     sprpen   += 100.*square(SBx(2)/SBx(1)-0.3);
  
     B40 = SBx(4)*mean(recruitment(1982,endyr-recage));
     F55 = mF(6);
     F50 = mF(5);
     F45 = mF(4);
     F40 = mF(3);

     if(DEBUG_FLAG == 1) cout<<"SPR"<<endl;


FUNCTION Projection

  // |------------------------------------------------------------------------|
  // | ABUNDANCE FOR FIRST PROJECTION YEAR (endyr+1)
  // |------------------------------------------------------------------------|


  // |------------------------------------------------------------------------|
  // | AGE ONE
  // |------------------------------------------------------------------------|

     int k;

     if(mceval_phase()) 
       {
         // random_number_generator r(1000);
         // standard deviation of mean recruitment
         stdev_rec   = sqrt(norm2(value(log_rec_dev(1982,endyr-recage))-
                       mean(value(log_rec_dev(1982,endyr-recage))))/
                       (size_count(value(log_rec_dev))-1));
                          
         k=k+1;
            
         // randomized projected recruitment
         N_ABC(endyr+1,1) = mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                       stdev_rec*stdev_rec/2+stdev_rec*randn(k)));

         N_OFL(endyr+1,1) = mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                       stdev_rec*stdev_rec/2+stdev_rec*randn(k)));

       }
          
     else 
       {
         N_ABC(endyr+1,1) = value(mean(recruitment(1982,endyr-recage)));
         N_OFL(endyr+1,1) = value(mean(recruitment(1982,endyr-recage)));
       }
     
  // |------------------------------------------------------------------------|
  // | AGES 8 - 75
  // |------------------------------------------------------------------------|

     for(j=2; j<=nages;j++) 
       {
          N_ABC(endyr+1,j) = natage(endyr,j-1)*S(endyr,j-1);
          N_OFL(endyr+1,j) = natage(endyr,j-1)*S(endyr,j-1);
       }

  // |------------------------------------------------------------------------|
  // | PLUS CLASS
  // |------------------------------------------------------------------------|

     N_ABC(endyr+1,nages) += natage(endyr,nages)*S(endyr,nages);
     N_OFL(endyr+1,nages) += natage(endyr,nages)*S(endyr,nages);

    
  // |------------------------------------------------------------------------|
  // | PROPAGATE FORWARD WITH ABC / OFL MALITIES
  // |------------------------------------------------------------------------|


     // |---------------------------------------------------------------------|
     // | TIER 4 NOAA REGULATIONS
     // |---------------------------------------------------------------------|
     // | FABC = F40
     // | FOFL = F35
     
        for (j=1;j<=nages;j++)
          {  
            FOFL(j) = fish_sel(j) * F50;
            FABC(j) = fish_sel(j) * F55;

            ZABC(j) = FABC(j) + M;
            ZOFL(j) = FOFL(j) + M;
          }


  // |------------------------------------------------------------------------|
  // | AGES TWO THROUGH +
  // |------------------------------------------------------------------------|

 for (i=endyr+2;i<=endyr+15;i++)
       {

         if(mceval_phase()) 
           {
             // random_number_generator r(1000);
             // standard deviation of mean recruitment
             stdev_rec   = sqrt(norm2(value(log_rec_dev(1982,endyr-recage))-
                           mean(value(log_rec_dev(1982,endyr-recage))))/
                           (size_count(value(log_rec_dev))-1));
                          
             k=k+1;
            
             // randomized projected recruitment
             N_ABC(i,1) = mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                   stdev_rec*stdev_rec/2+stdev_rec*randn(k)));

             N_OFL(i,1) = mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                   stdev_rec*stdev_rec/2+stdev_rec*randn(k)));
           }
          
         else 
           {
             N_ABC(i,1) = value(mean(recruitment(1982,endyr-recage)));     
             N_OFL(i,1) = value(mean(recruitment(1982,endyr-recage)));


       }
         
           for(j=2; j<=nages;j++) 
             {
               N_ABC(i,j) = N_ABC(i-1,j-1)* mfexp(-FABC(j-1)-M);
               N_OFL(i,j) = N_OFL(i-1,j-1)* mfexp(-FOFL(j-1)-M);;
             }

         N_ABC(i,nages) += N_ABC(i,nages) * mfexp(-FABC(nages)-M);
         N_OFL(i,nages) += N_OFL(i,nages) * mfexp(-FOFL(nages)-M);

       }


  // |------------------------------------------------------------------------|
  // | CATCH AND BIOMASS
  // |------------------------------------------------------------------------|

     for (i=endyr+1; i<=endyr+15; i++)
       {
         for (j=1; j<=nages; j++)
           {
             catage_ABC(i) = elem_div(elem_prod(elem_prod(N_ABC(i),FABC),
                                     (1.-mfexp(-ZABC))),ZABC);
             catage_OFL(i) = elem_div(elem_prod(elem_prod(N_OFL(i),FOFL),
                                     (1.-mfexp(-ZOFL))),ZOFL);
     
             catch_ABC(i) = catage_ABC(i)*wt_avg_fshy;
             catch_OFL(i) = catage_OFL(i)*wt_avg_fshy;

             spawn_biom_ABC(i) = (N_ABC(i) * mfexp(-spawn_fract * M))
                                           * wt_mature_f;
             spawn_biom_OFL(i) = (N_OFL(i) * mfexp(-spawn_fract * M))
                                           * wt_mature_f;

           }
       }
       
     ABC=catch_ABC(endyr+1);
     OFL=catch_OFL(endyr+1);

     ABC = 2204 * ABC;
     OFL = 2204 * OFL;

    

     if(DEBUG_FLAG == 1) cout<<"Projections"<<endl;
          


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

     //int nboot;
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
     log_F_devs.initialize();
     mF.initialize();


     log_mean_rec =   bootN(1);             // Mean recruitment;
     log_mean_y1  =   bootN(2);             // Mean Year 1 naa;
     log_avg_F    =   bootN(3);             // Mean F
     log_comm_q    =  bootN(4);             // Catchability (commercial CPUE)
     log_surv_q1   =  bootN(5);             // Catchability (survey CPUE 1)
     log_surv_q2   =  bootN(6);             // Catchability (survey CPUE 2)


     for(ii = 7; ii<=46; ii++)
       {
         init_pop(ii-6) = bootN(ii);          // Year 1 abundance deviations
       }

     for(ii = 47; ii<=82; ii++)
       {
         log_rec_dev(ii+styr-47) = bootN(ii);   // Recruitment deviations
       }
 

    for(ii = 83; ii<=118; ii++)
      {
        log_F_devs(ii - 83+styr) = bootN(ii);      // F deviations, females
      }


    for(ii = 119; ii<=154; ii++)
      {
        log_F_devs(ii-119 + styr) = bootN(ii);  // F deviations, males
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
     Population();
     Predicted();

     //spr();	
     //Get_ABC();

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


FUNCTION double sdnr(const dvar_vector& pred,const dvector& obs,double m)
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred);
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;


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

        template<typename T>
        dvar_vector plogis(const dvector x, T location, T scale)
          {
            return(1.0 / (1.0 + exp(-(x-location)/scale)));
          }

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
  REPORT(ABC)
  REPORT(obs_MR_abund)
  REPORT(pred_MR_abund)
  REPORT(obs_MR_abund_SRV)
  REPORT(pred_MR_abund_SRV)
  REPORT(natage);
  REPORT(spawn_biom);
  REPORT(fish_sel);
  REPORT(srv_sel);
  REPORT(pred_catch);
  REPORT(mF);
  REPORT(N0);
  REPORT(Nspr);
  REPORT(SBx);
  REPORT(tot_biomass);
  REPORT(catch_ABC);
  REPORT(obs_cpue1_biom);
  REPORT(pred_cpue1);
  REPORT(obs_cpue2_biom);
  REPORT(pred_cpue2);
  REPORT(obs_fshy_cpue);
  REPORT(pred_fshy_cpue);
  REPORT(wt_mature_s);
  REPORT(mfexp(-spawn_fract*M));
    //SDNR
  REPORT(effn_fish_age);
  REPORT(sdnr_fish_age);
  REPORT(effn_srv_age);
  REPORT(sdnr_srv_age);

 

