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
   number Mm;
   number Mf;
   number srv_slpm;
   number srv_slpf;
   number fsh_slpm;
   number fsh_slpf;
    
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
          srv_slpm    = dMiscCont(6);
          srv_slpf    = dMiscCont(7);
          fsh_slpm    = dMiscCont(8);
          fsh_slpf    = dMiscCont(9);
          ph_rec      = dMiscCont(10);
          ph_init     = dMiscCont(11);
          ph_F        = dMiscCont(12);
          ph_spr      = dMiscCont(13);


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
  // |    
               
     init_int      nages
     init_int      styr
     init_int      endyr
     init_int      recage

     init_number   spawn_fract
     init_number   srv_fract
     init_number   fshy_fract

     vector        yy(styr,endyr)
     vector        aa(1,nages)


  // |--------------------------------------------------------------------------|
  // | Gear selectivity for fishery and survey (using NOAA values)
  // |--------------------------------------------------------------------------|
  // | Selectivity curves from Dana Hanselman @ NOAA
  
     init_vector	vfish_sel_f(1,nages)
     init_vector	vfish_sel_m(1,nages)
     init_vector	srv_sel_f(1,nages)
     init_vector	srv_sel_m(1,nages)

  
  // |--------------------------------------------------------------------------|
  // | PHYSIOLOGY                                                               
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) p_mature       -> proportion of females mature-at-age
  // | 2) wt_avg_fshy_f  -> mean weight-at-age in commercial fishery, females
  // | 3) wt_avg_fshy_m  -> mean weight-at-age in commercial fishery, males
  // | 4) wt_avg_srv_f  -> mean weight-at-age in longline survey, females
  // | 5) wt_avg_srv_m  -> mean weight-at-age in longline survey, males

     init_vector   p_mature(1,nages)
     init_vector   wt_avg_fshy_f(1,nages)
     init_vector   wt_avg_fshy_m(1,nages)
     init_vector   wt_avg_srv_f(1,nages)
     init_vector   wt_avg_srv_m(1,nages)


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
  // | 3) nsamples_fish_agef  -> sample size FEMALES
  // | 4) oac_fishf           -> female age comp data
  // | 5) nsamples_fish_agem  -> sample size MALES
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
  // | 1) nyrs_srv_age       -> number of years of survey age comps
  // | 2) yrs_srv_age        -> specific years age comps
  // | 3) nsamples_srv_agef  -> sample size FEMALES
  // | 4) oac_srvf           -> female age comp data
  // | 5) nsamples_srv_agem  -> sample size MALES
  // | 6) oac_srvm           -> male age comp data
  // |
  // | Relative sample size = square root of number of fish aged to scale data impact

     init_int      nyrs_srv_age 	
     init_ivector  yrs_srv_age(1,nyrs_srv_age)   
     init_vector   nsamples_srv_agef(1,nyrs_srv_age) 
     init_matrix   oac_srvf(1,nyrs_srv_age,1,nages)  
     init_vector   nsamples_srv_agem(1,nyrs_srv_age) 
     init_matrix   oac_srvm(1,nyrs_srv_age,1,nages)  
  
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
     srv_fract   = (srv_fract - 1)   / 12;      // Fraction of year at which survey occurs
     fshy_fract  = (fshy_fract - 1) / 12;	// Fraction of year at which fishery occurs


    offset.initialize();


    for (i=1; i<=nyrs_fish_age; i++)
      {
        oac_fishm(i)/=sum(oac_fishm(i));
        oac_fishf(i)/=sum(oac_fishf(i));
        offset(1) -= nsamples_fish_agem(i)*((oac_fishm(i) + 0.00001)*log(oac_fishm(i) + 0.00001));
        offset(2) -= nsamples_fish_agef(i)*((oac_fishf(i) + 0.00001)*log(oac_fishf(i) + 0.00001));
      }

     for (i=1; i<=nyrs_srv_age; i++)
       {
         oac_srvm(i)/=sum(oac_srvm(i));
         oac_srvf(i)/=sum(oac_srvf(i));
         offset(3) -= nsamples_srv_agem(i)*((oac_srvm(i) + 0.00001)*log(oac_srvm(i) + 0.00001));
         offset(4) -= nsamples_srv_agef(i)*((oac_srvf(i) + 0.00001)*log(oac_srvf(i) + 0.00001));
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


  // |---------------------------------------------------------------------------------|
  // | FISHING MORTALITY PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | - log_F_devsf        <- annual fishing mortality deviations females
  // | - log_F_devsm        <- annual fishing mortality deviations males
  // |
  // | - Fmortf             <- Ff * selectivity
  // | - Fmortm             <- Fm * selectivity
  // |
  // | - fish_sel_a50_f     <- age of 50% selectivity, fishery females
  // | - fish_sel_a50_m     <- age of 50% selectivity, fishery males
  // |
  // | - vector for fishery selectivity females
  // | - vector for fishery selectivity males

     init_bounded_vector  log_F_devsf(styr,endyr,-10.,10.,ph_F);
     init_bounded_vector  log_F_devsm(styr,endyr,-10.,10.,ph_F);
    
     vector               Fmortf(styr,endyr);
     vector               Fmortm(styr,endyr);

     init_bounded_number fish_sel_a50_f(1,3,ph_init);
     init_bounded_number fish_sel_a50_m(1,3,ph_init);

     vector	fish_sel_f(1,nages)
     vector	fish_sel_m(1,nages)
     

  // |---------------------------------------------------------------------------------|
  // | POPULATION VECTORS
  // |---------------------------------------------------------------------------------|
  // | - tot_n          -> total abundance
  // | - tot_nf         -> total abundance females
  // | - tot_nm         -> total abundance males
  // | - tot_biomass    -> total biomass
  // | - spawn_biomass  -> spawning biomass
  // | - recruitment    -> recruitment
  // |
  // | - pred_catch     -> total annual catch
  // | - pred_catchf    -> total annual catch females
  // | - pred_catchm    -> total annual catch males
	
     sdreport_vector  tot_n(styr,endyr)
     sdreport_vector  tot_nf(styr,endyr)
     sdreport_vector  tot_nm(styr,endyr)
     sdreport_vector  tot_biomass(styr,endyr)
     sdreport_vector  spawn_biom(styr,endyr);
     sdreport_vector  recruitment(styr,endyr);
     
     vector           pred_catch(styr,endyr)		
     vector           pred_catchf(styr,endyr)
     vector           pred_catchm(styr,endyr)

     vector           wt_mature_s(1,nages);
     vector           wt_mature_f(1,nages);
     vector           wt_fsh_sel_m(1,nages);
     vector           wt_fsh_sel_f(1,nages);
  

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
  // | - eac_srvf      -> proportion-at-age survey females
  // | - eac_srvm      -> proportion-at-age survey males

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
     matrix  eac_srvf(1,nyrs_srv_age,1,nages); 
     matrix  eac_srvm(1,nyrs_srv_age,1,nages);


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

     vector effn_srv_agef(1,nyrs_srv_age)
     vector effn_srv_agem(1,nyrs_srv_age)
     vector sdnr_srv_agef(1,nyrs_srv_age)
     vector sdnr_srv_agem(1,nyrs_srv_age)
     

  // |---------------------------------------------------------------------------------|
  // | MARK-RECAPTURE AND CPUE
  // |---------------------------------------------------------------------------------|
  // | - pred_cpue1          -> Survey CPUE (both sexes) 1988 - 1996
  // | - pred_cpue2          -> Survey CPUE (both sexes) 1997 - present
  // |
  // | - pred_fshy_pue       -> Fishery CPUE (both sexes)
  // |
  // | - pred_MR_abund       -> MR abundance fishery
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

     vector  pred_cpue1(1,nyrs_cpue1); 
     vector  pred_cpue2(1,nyrs_cpue2);

     vector  pred_fshy_cpue(1,nyrs_fshy_cpue);

     sdreport_vector  pred_MR_abund(1,nyrs_MR);
     sdreport_vector  pred_MR_abund_SRV(1,nyrs_MR);

     number  s_srvf; 
     number  s_fshyf;
     number  s_srvm;  
     number  s_fshym; 

     sdreport_number  q_cpue1;
     sdreport_number  q_cpue2;
     sdreport_number  q_fshy;


  // |--------------------------------------------------------------------------|
  // | SPR / ABC / PROJECTIONS                               
  // |--------------------------------------------------------------------------|
  // | mF                -> vector of potential F levels
  // | SBx               -> spawning biomass from mF            
  // | Bx                -> biomass from mF
  // |
  // | Nspr              -> Number of projected spawners
  // | N0                -> Number projected with F = 0
  // |
  // | N_ABC_f           -> Number of projected females at F_ABC
  // | N_ABC_m           -> Number of projected males at F_ABC
  // | N_OFL_f           -> Number of projected females at F_OFL
  // | N_OFL_m           -> Number of projected males at F_OFL
  // |
  // | FABC_f            -> F * sel of projected females at F_ABC
  // | FABC_m            -> F * sel of projected males at F_ABC
  // | FOFL_f            -> F * sel of projected females at F_OFL
  // | FOFL_m            -> F * sel of projected males at F_OFL
  // | 
  // | ZABC_f            -> Total mortality of projected females at F_ABC
  // | ZABC_m            -> Total mortality of projected males at F_ABC
  // | ZOFL_f            -> Total mortality of projected females at F_OFL
  // | ZOFL_m            -> Total mortality of projected males at F_OFL
  // |
  // | catage_ABC_f      -> Projected catch at age for females at F_ABC
  // | catage_ABC_m      -> Projected catch at age for males at F_ABC
  // | catage_OFL_f      -> Projected catch at age for females at F_OFL
  // | catage_OFL_m      -> Projected catch at age for males at F_OFL
  // |
  // | catch_ABC_f      -> Projected annual catch for females at F_ABC
  // | catch_ABC_m      -> Projected annual catch for males at F_ABC
  // | catch_OFL_f      -> Projected annual catch for females at F_OFL
  // | catch_OFL_m      -> Projected annual catch for males at F_OFL
  // |

     init_bounded_vector mF(1,9,0.01,1,ph_spr)
     vector SBx(1,10)
     vector Bx(1,10)
					       
     matrix Nspr(1,10,1,nages)			
     matrix N0(1,10,1,nages)

     matrix N_ABC_f(endyr+1,endyr+15,1,nages);
     matrix N_ABC_m(endyr+1,endyr+15,1,nages);
     matrix N_OFL_f(endyr+1,endyr+15,1,nages);
     matrix N_OFL_m(endyr+1,endyr+15,1,nages);

     vector FABC_f(1,nages);
     vector FABC_m(1,nages);
     vector FOFL_f(1,nages);
     vector FOFL_m(1,nages)
     
     vector ZABC_f(1,nages);
     vector ZABC_m(1,nages);
     vector ZOFL_f(1,nages);
     vector ZOFL_m(1,nages);
     
     matrix catage_ABC_f(endyr+1,endyr+15,1,nages);
     matrix catage_ABC_m(endyr+1,endyr+15,1,nages);
     matrix catage_OFL_f(endyr+1,endyr+15,1,nages);
     matrix catage_OFL_m(endyr+1,endyr+15,1,nages);

     vector catch_ABC_f(endyr+1,endyr+15);
     vector catch_ABC_m(endyr+1,endyr+15);
     vector catch_OFL_f(endyr+1,endyr+15);
     vector catch_OFL_m(endyr+1,endyr+15);

     matrix spawn_biom_ABC(endyr+1,endyr+15,1,nages);
     matrix spawn_biom_OFL(endyr+1,endyr+15,1,nages);

     sdreport_number B40;
     sdreport_number ABC;
     sdreport_number OFL;
     number ABC_pounds;
     number OFL_pounds;
     number F50;
     number F45;
     number F40;
     number F55;

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
     vector  age_like(1,4);
     number  rec_like;
     vector  catch_like(1,2);
     number  F_mort_regularity;

     number  Like;
     objective_function_value obj_fun;	

						
 LOCAL_CALCS

   wt_mature_f = elem_prod(p_mature,wt_avg_fshy_f);
   wt_mature_s = elem_prod(p_mature,wt_avg_srv_f);
   wt_fsh_sel_f = elem_prod(fish_sel_f,wt_avg_fshy_f);
   wt_fsh_sel_m = elem_prod(fish_sel_m,wt_avg_fshy_m);


  // |-------------------------------------------------------------------------|
  // | Survival fractions to decrement abundance relative to annual timing
  // |-------------------------------------------------------------------------|
  // |
  // | This will make a difference if distinct natural mortalities are used
  // | between males and females. If both are input at 0.1, then one decrement
  // | for the fishery and another for the survey is all that is required

     s_srvf  = mfexp(-1*srv_fract*Mf);
     s_fshyf  = mfexp(-1*fshy_fract*Mf);
     s_srvm  = mfexp(-1*srv_fract*Mm);
     s_fshym  = mfexp(-1*fshy_fract*Mm);

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
         //     tot_nf<<" "<<
         //     tot_nm<<" "<<
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

	ofstream ofs5("tot_nf.ps",ios::app);
	ofs5<<tot_nf<<endl;

	ofstream ofs6("tot_nm.ps",ios::app);
	ofs6<<tot_nm<<endl;
        

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
    
    //for (j=1;j<=nages;j++)
      //{
        //srv_sel_m(j) = 1/(1+mfexp(-mfexp(srv_slpm)*(j-mfexp(srv_sel_a50_m))));
        //srv_sel_f(j) = 1/(1+mfexp(-mfexp(srv_slpf)*(j-mfexp(srv_sel_a50_f))));        
      //}


   //  |-----------------------------------------------------|
   //  | Fishery
   //  |-----------------------------------------------------|

    for (j=1;j<=nages;j++)
      { 
        fish_sel_f(j) =  1/(1+mfexp(-mfexp(fsh_slpf)*(j-mfexp(fish_sel_a50_f))));
        fish_sel_m(j) =  1/(1+mfexp(-mfexp(fsh_slpm)*(j-mfexp(fish_sel_a50_m))));
      }



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

     dvar_matrix mr_f(1,nyrs_MR,1,nages);
     dvar_matrix mr_m(1,nyrs_MR,1,nages);
     dvar_matrix mr_f_SRV(1,nyrs_MR,1,nages);
     dvar_matrix mr_m_SRV(1,nyrs_MR,1,nages);

  
     for (i=1;i<=nyrs_MR;i++)
       {
         for(j=1; j<=nages; j++)
           {
             mr_f(i,j) = natagef(yrs_MR(i),j)*fish_sel_f(j);
             mr_m(i,j) = natagem(yrs_MR(i),j)*fish_sel_m(j);
           }
       }

     for (i=1;i<=nyrs_MR;i++)
       {
             mr_f_SRV(i) = elem_prod(natagef(yrs_MR(i)),srv_sel_f);
             mr_m_SRV(i) = elem_prod(natagem(yrs_MR(i)),srv_sel_m);
       }


     for (i=1;i<=nyrs_MR;i++)
       {
         pred_MR_abund(i) =   s_fshyf * 1000 * sum(mr_f(i))+
                              s_fshym * 1000 * sum(mr_m(i));

         pred_MR_abund_SRV(i) =   s_srvf * 1000 * sum(mr_f_SRV(i))+
                                  s_srvm * 1000 * sum(mr_m_SRV(i));
                              
       } 

                              
      
  // |-------------------------------------------------------------------------|
  // | ADF&G LONGLINE SURVEY CPUE 1988 - 1996
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance

     for (i=1;i<=nyrs_cpue1;i++)
       {                        
         pred_cpue1(i) = surv_q1 * (sum(elem_prod(natagef(yrs_cpue1(i)),
                                       srv_sel_f)) +  
                                    sum(elem_prod(natagem(yrs_cpue1(i)),
                                       srv_sel_m)));  
       }


  // |-------------------------------------------------------------------------|
  // | ADF&G LONGLINE SURVEY CPUE 1997 - PRESENT
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance
					             
     for (i=1;i<=nyrs_cpue2;i++)
       {
         pred_cpue2(i) = surv_q2 * (sum(elem_prod(natagef(yrs_cpue2(i)),
                                       srv_sel_f))+  
                                   sum(elem_prod(natagem(yrs_cpue2(i)),
                                       srv_sel_m)));  
       }


  // |-------------------------------------------------------------------------|
  // | COMMERCIAL FISHERY CPUE
  // |-------------------------------------------------------------------------|
  // | Note that the decrement is not implemented on the CPUE calcs - it is
  // | absorbed into the catchability coefficient and would just add code
  // | with no effect on relative abundance


     for (i=1;i<=nyrs_fshy_cpue;i++)
       {

         pred_fshy_cpue(i) =  comm_q * (sum(elem_prod(natagef(yrs_fshy_cpue(i)),wt_fsh_sel_f))
                                           +
                                        sum(elem_prod(natagem(yrs_fshy_cpue(i)),wt_fsh_sel_m)));
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

         sdnr_fish_agef(i) = sdnr(eac_fishf(i),oac_fishf(i),
                            double(nsamples_fish_agef(i)));
                            
         sdnr_fish_agem(i) = sdnr(eac_fishm(i),oac_fishm(i),
                            double(nsamples_fish_agem(i)));
       }


     // SURVEY
     for (i=1;i<=nyrs_srv_age;i++)
       {
         eac_srvf(i)  = elem_prod(srv_sel_f, natagef(yrs_srv_age(i))*s_srvf)
                         /sum(elem_prod(natagef(yrs_srv_age(i))*s_srvf,srv_sel_f)) * ageage; 

         eac_srvm(i)  = elem_prod(srv_sel_m, natagem(yrs_srv_age(i))*s_srvm)
                         /sum(elem_prod(natagem(yrs_srv_age(i))*s_srvm,srv_sel_m)) * ageage;

         effn_srv_agem(i) = (1-eac_srvm(i))*eac_srvm(i)/norm2(oac_srvm(i)-eac_srvm(i)); 
         effn_srv_agef(i) = (1-eac_srvf(i))*eac_srvf(i)/norm2(oac_srvf(i)-eac_srvf(i));

         sdnr_srv_agef(i) = sdnr(eac_srvf(i),oac_srvf(i),
                            double(nsamples_srv_agef(i)));
                            
         sdnr_srv_agem(i) = sdnr(eac_srvm(i),oac_srvm(i),
                            double(nsamples_srv_agem(i)));
       }

FUNCTION Population
   recruitment.initialize();
   tot_biomass.initialize();
   spawn_biom.initialize();
   tot_nm.initialize();
   tot_nf.initialize();
   tot_n.initialize();

  // |-------------------------------------------------------------------------|
  // | MAKE SURE UNITS ARE CORRECT
  // |-------------------------------------------------------------------------|
  // |*Note conversion from metric tonnes (for total biomass) to straight pounds
  // |  for Age 4 and Exploitable biomass calcs - this is for managers
  // |
  // |  - kg (wt-at-age)  x thousands (natage)  = metric tonnes
  // |  - 1 metric tonne = 2,204 pounds


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
         
     for (i=styr;i<=endyr;i++)
       {
         recruitment(i) = natage(i,1);
    
         tot_biomass(i) = (natagef(i)*wt_avg_srv_f)+(natagem(i)*wt_avg_srv_m); 
    	       
         spawn_biom(i)  = natagef(i) *mfexp(-spawn_fract*Mf) * wt_mature_s; 
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

  catch_like(1)  +=  0.5*log(2*M_PI)  +log(sigma_catch) +  0.5*(square(log(obs_catchf(i)+oo) -
                                     log(pred_catchf(i)+oo))) / (2.*square(sigma_catch));

  catch_like(2)  +=  0.5*log(2*M_PI)  +log(sigma_catch) +  0.5*(square(log(obs_catchm(i)+oo) -
                                     log(pred_catchm(i)+oo))) / (2.*square(sigma_catch));

  }



FUNCTION Surv_Likelihood

  surv_like.initialize();

  // |----------------------------------------------------------------------|
  // | Mark-recapture estimates of abundance - normal
  // |----------------------------------------------------------------------|

    for (i=1; i<=nyrs_MR; i++)
      {
   //
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
  // | ADF&G survey  CPUE - normal 1996 - present
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

       for (i=1; i <= nyrs_srv_age; i++)
         {
           age_like(3) -= nsamples_srv_agem(i)*((oac_srvm(i) + oo) *
                                                  log(eac_srvm(i) + oo)) ;
                                                  
           age_like(4) -= nsamples_srv_agef(i)*((oac_srvf(i) + oo) *
                                                  log(eac_srvf(i) + oo)) ; 
         }


    age_like(1)   -= offset(1);
    age_like(2)   -= offset(2);
    age_like(3)   -= offset(3);
    age_like(4)   -= offset(4); 


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

     Like            = sum(catch_like);	
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
  // | Note that this is only calculated for female spawning biomass
          

  // |----------------------------------------------------------------------|
  // | RECRUITMENT (AGE 2)
  // |----------------------------------------------------------------------|

     for (i=1;i<=10;i++)
       {
         Nspr(i,1) = 1.;
         N0(i,1) = mean(recruitment(1982,endyr-recage))/2;
       }


  // |----------------------------------------------------------------------|
  // | AGES 3 - 42
  // |----------------------------------------------------------------------|

     for (j=2;j<=nages;j++)
       {
         Nspr(1,j) = Nspr(1,j-1) * mfexp(-1.*Mf);
         N0(1,j)   = N0(1,j-1)   * mfexp(-1.*Mf);

       for (i=2; i<=10; i++)
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


     for(i=2; i<=10; i++)
      {
         Nspr(i,nages) += Nspr(i,nages) * mfexp(-1.*(Mf+mF(i-1)*
                                                fish_sel_f(nages)));
         N0(i,nages)   += N0(i,nages)   * mfexp(-1.*(Mf+mF(i-1)*
                                                fish_sel_f(nages)));
       }


  // |----------------------------------------------------------------------|
  // | PROJECTED BIOMASS FOR F LEVELS
  // |----------------------------------------------------------------------|

     //SBx(1) +=  2.20462 * Nspr(1)*wt_mature_s;

     SBx(1) +=  Nspr(1)*wt_mature_s*mfexp(-spawn_fract*Mf);

     
     for (j=1;j<=nages;j++)
       {
         for(i=2; i<=10; i++)
           {
             SBx(i) += Nspr(i,j)*wt_mature_s(j)*
                       mfexp(-spawn_fract*(Mf+mF(i-1)*fish_sel_f(j)));;
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


FUNCTION Projection

  // |--------------------------------------------------------------------------|
  // | ABUNDANCE FOR FIRST PROJECTION YEAR (endyr+1)
  // |--------------------------------------------------------------------------|
     void ("Note that all the jumping up and down that follows can be         ");
     void ("largely circumnagivated with a singe call:                        ");
     void ("ABC =  2204*sum(                                                  ");
     void ("       ((elem_prod(wt_avg_fshy_f,elem_prod(elem_div(Fself, Z1),   ");
     void ("                   elem_prod(1-S1,(natage_projf))))) +            ");
     void ("        elem_prod(wt_avg_fshy_m,elem_prod(elem_div(Fselm, Z2),    ");
     void ("                   elem_prod(1-S2,(natage_projm))))) );           ");
     void (" This way just makes it easier to go step by step and mod things  ");
  //----------------------------------------------------------------------------|

 
  // |-------------------------------------------------------------------------|
  // | Population Projection (endyr + 1)
  // |-------------------------------------------------------------------------|

    // |-----------------------------------------------------------------------|
    // | Recruitment (age 2)
    // |-----------------------------------------------------------------------|
    // | Mean over recent 15 years, excluding last two model years (variability)

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
         N_ABC_f(endyr+1,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                     stdev_rec*stdev_rec/2+stdev_rec*randn(k))));

         N_ABC_m(endyr+1,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                     stdev_rec*stdev_rec/2+stdev_rec*randn(k)))); 

         N_OFL_f(endyr+1,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                     stdev_rec*stdev_rec/2+stdev_rec*randn(k))));

         N_OFL_m(endyr+1,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                     stdev_rec*stdev_rec/2+stdev_rec*randn(k))));
       }
          
     else 
       {
         N_ABC_f(endyr+1,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));
         N_ABC_m(endyr+1,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));
                  
         N_OFL_f(endyr+1,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));
         N_OFL_m(endyr+1,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));

       }

                                

    // |-----------------------------------------------------------------------|
    // | AGES 3 - 42
    // |-----------------------------------------------------------------------|

     for(j=2; j<=nages;j++) 
       {
          N_ABC_f(endyr+1,j) = natagef(endyr,j-1)*Sf(endyr,j-1);
          N_ABC_m(endyr+1,j) = natagem(endyr,j-1)*Sm(endyr,j-1);
                    
          N_OFL_f(endyr+1,j) = natagef(endyr,j-1)*Sf(endyr,j-1);
          N_OFL_m(endyr+1,j) = natagem(endyr,j-1)*Sm(endyr,j-1);

       }

       
    // |-----------------------------------------------------------------------|
    // | PLUS CLASS
    // |-----------------------------------------------------------------------|


          N_ABC_f(endyr+1,nages) = natagef(endyr,nages)*Sf(endyr,nages);
          N_ABC_m(endyr+1,nages) = natagem(endyr,nages)*Sm(endyr,nages);
                    
          N_OFL_f(endyr+1,nages) = natagef(endyr,nages)*Sf(endyr,nages);
          N_OFL_m(endyr+1,nages) = natagem(endyr,nages)*Sm(endyr,nages);



  // |-------------------------------------------------------------------------|
  // | Note that the calcs below reflect the ADF&G regulations for sablefish
  // |   and *NOT* standard NOAA federal curves as per Dorn
  // |
  // | Note that fishery selectivities and fishery weights used to compute ABC

     ABC.initialize();

     for (j=1;j<=nages;j++)
       {  
         FABC_f(j) = fish_sel_f(j) * F50;
         FABC_m(j) = fish_sel_m(j) * F50;
         
         FOFL_f(j) = fish_sel_f(j) * F45;
         FOFL_m(j) = fish_sel_m(j) * F45;

         ZABC_f(j) = FABC_f(j) + Mf;
         ZABC_m(j) = FABC_m(j) + Mm;
         ZOFL_f(j) = FOFL_f(j) + Mf;
         ZOFL_m(j) = FOFL_m(j) + Mm;

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
             N_ABC_f(i,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                   stdev_rec*stdev_rec/2+stdev_rec*randn(k))));

             N_ABC_m(i,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                   stdev_rec*stdev_rec/2+stdev_rec*randn(k)))); 

             N_OFL_f(i,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                   stdev_rec*stdev_rec/2+stdev_rec*randn(k))));

             N_OFL_m(i,1) = 0.5 * (mfexp(value(log(mean(recruitment(1982,endyr-recage)))+
                                   stdev_rec*stdev_rec/2+stdev_rec*randn(k))));
           }
          
         else 
           {
             N_ABC_f(i,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));
             N_ABC_m(i,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));
                  
             N_OFL_f(i,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));
             N_OFL_m(i,1) = 0.5 * (value(mean(recruitment(1982,endyr-recage))));

       }
         
           for(j=2; j<=nages;j++) 
             {
               N_ABC_f(i,j) = N_ABC_f(i-1,j-1)* mfexp(-FABC_f(j-1)-Mf);
               N_ABC_m(i,j) = N_ABC_m(i-1,j-1)* mfexp(-FABC_m(j-1)-Mm);
               
               N_OFL_f(i,j) = N_OFL_f(i-1,j-1)* mfexp(-FOFL_f(j-1)-Mf);
               N_OFL_m(i,j) = N_OFL_m(i-1,j-1)* mfexp(-FOFL_m(j-1)-Mm);
             }


         N_ABC_f(i,nages) += N_ABC_f(i,nages) * mfexp(-FABC_f(nages)-Mf);
         N_ABC_m(i,nages) += N_ABC_m(i,nages) * mfexp(-FABC_m(nages)-Mm);

         N_OFL_f(i,nages) += N_OFL_f(i,nages) * mfexp(-FOFL_f(nages)-Mf);
         N_OFL_m(i,nages) += N_OFL_m(i,nages) * mfexp(-FOFL_m(nages)-Mm);
         
       }


  // |------------------------------------------------------------------------|
  // | CATCH AND BIOMASS
  // |------------------------------------------------------------------------|

     for (i=endyr+1; i<=endyr+15; i++)
       {
         for (j=1; j<=nages; j++)
           {
             catage_ABC_f(i) = elem_div(elem_prod(elem_prod(N_ABC_f(i),FABC_f),
                                       (1.-mfexp(-ZABC_f))),ZABC_f);
             catage_ABC_m(i) = elem_div(elem_prod(elem_prod(N_ABC_m(i),FABC_m),
                                       (1.-mfexp(-ZABC_m))),ZABC_m);
             
             catage_OFL_f(i) = elem_div(elem_prod(elem_prod(N_OFL_f(i),FOFL_f),
                                                  (1.-mfexp(-ZOFL_f))),ZOFL_f);
             catage_OFL_m(i) = elem_div(elem_prod(elem_prod(N_OFL_m(i),FOFL_m),
                                                  (1.-mfexp(-ZOFL_m))),ZOFL_m);

             catch_ABC_f(i) = catage_ABC_f(i) * wt_avg_fshy_f;
             catch_ABC_m(i) = catage_ABC_m(i) * wt_avg_fshy_m;
             
             catch_OFL_f(i) = catage_OFL_f(i) * wt_avg_fshy_f;
             catch_OFL_m(i) = catage_OFL_m(i) * wt_avg_fshy_m;
 
             
             spawn_biom_ABC(i) = (N_ABC_f(i) * mfexp(-spawn_fract * Mf))
                                             * wt_mature_f;
             spawn_biom_OFL(i) = (N_OFL_f(i) * mfexp(-spawn_fract * Mf))
                                             * wt_mature_f;

           }
       }


     ABC = catch_ABC_f(endyr+1)+catch_ABC_m(endyr+1);
     OFL = catch_OFL_f(endyr+1)+catch_OFL_m(endyr+1);

     ABC_pounds = 2204 * (ABC);
     OFL_pounds = 2204 * (OFL);






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
     log_F_devsf.initialize();
     log_F_devsm.initialize();
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
        ofs9<<eac_srvf<<endl;
    }

   ofstream ofs10("agem_srv.txt",ios::app);
    {
        ofs10<<eac_srvm<<endl;
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
  REPORT(ABC_pounds)
  REPORT(obs_MR_abund)
  REPORT(pred_MR_abund)
  REPORT(obs_MR_abund_SRV)
  REPORT(pred_MR_abund_SRV)
  REPORT(pred_catchf);
  REPORT(pred_catchm);
  REPORT(pred_catch);
  REPORT(natagef);
  REPORT(natagem);
  REPORT(mF);
  REPORT(tot_biomass);
  REPORT(catage_ABC_f);
  REPORT(catage_ABC_m);
  REPORT(catch_ABC_f);
  REPORT(catch_ABC_m);
  REPORT(obs_cpue1_biom);
  REPORT(pred_cpue1);
  REPORT(obs_cpue2_biom);
  REPORT(pred_cpue2);
  REPORT(obs_fshy_cpue);
  REPORT(pred_fshy_cpue);
  REPORT(tot_nm);
  REPORT(tot_nf);
  REPORT(tot_n);
  REPORT(wt_mature_s);
  REPORT(mfexp(-spawn_fract*Mf));
      //SDNR
  REPORT(effn_fish_agef);
  REPORT(sdnr_fish_agef);
  REPORT(effn_fish_agem);
  REPORT(sdnr_fish_agem);
  
  REPORT(effn_srv_agef);
  REPORT(sdnr_srv_agef);
  REPORT(effn_srv_agem);
  REPORT(sdnr_srv_agem);


