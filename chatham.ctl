## -------------------------------------------------------------------------- ##
##
## THIS IS THE CONTROL FILE FOR CHATHAM STRAIT SABLEFISH
## 
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
#  DEBUG FLAG (0 = FALSE, 1 = TRUE)
#
	0
#
## -------------------------------------------------------------------------- ##

## -------------------------------------------------------------------------- ##
#  NBOOT for number of bootstraps to run (0 = NOT RUN)
#
	0
#
## -------------------------------------------------------------------------- ##


## —————————————————————————————————————————————————————————————————————————— ##
##                  DESIGN MATRIX FOR PARAMETER CONTROLS                      ##
##  Prior descriptions   Parameter values                                     ##                                                                            
##  -0 uniform           (0,0)                                                ##
##  -1 normal            (p1=mu,p2=sig)                                       ##
##  -2 lognormal         (p1=log(mu),p2=sig)                                  ##
##  -3 beta              (p1=alpha,p2=beta)                                   ##
##  -4 gamma             (p1=alpha,p2=beta)                                   ##
## —————————————————————————————————————————————————————————————————————————— ##
##  init   lower    upper    est   prior
## value   bound    bound    phz    type    p1      p2    # PARAMETER         ##
## —————————————————————————————————————————————————————————————————————————— ##
## - init_int n_theta
   6
#
## - initial population
    5.0     0.00    10.00     1      1     0.0      1   # log_mean_rec
    4.0     0.00    8.00      1      1     0.0      1   # log_mean_y1
#
## - fishing mortality
   -3.5   -5.0     -2.0       2      1     0.0      1   # log_avg_F   
#
## - catchability
   -2.0    -20.0    0.00      2      0     0.0      0   # log_comm_q
   -2.0    -20.0    0.00      2      0     0.0      0   # log_surv_q1
   -2.0    -20.0    0.00      2      0     0.0      0   # log_surv_q2 
## —————————————————————————————————————————————————————————————————————————— ##


## —————————————————————————————————————————————————————————————————————————— ##
##                        OTHER MISCELLANEOUS CONTROLS                        ##
## —————————————————————————————————————————————————————————————————————————— ##
## number of controls to read in.
   9
## Fixed value   # Name         - Description
   0.05          # sigma_catch  - total annual catch variance
   1.2           # sigr         - recruitment variance (Sigler et al. 2002)
   1.2           # sig1         - as sigr
   0.1           # Mm           - male natural mortality (0.1 - Johnson and Quinn, 1988)      
   0.1           # Mf           - female natural mortality (0.1 - Johnson and Quinn, 1988)   
   4             # ph_rec       - recruitment deviations vector phase
   3             # ph_init      - Year 1 deviations vector  phase
   4             # ph_F         - commercial fisheries deviations vector phase
   5             # ph_spr       - spawner-recruit estimates phase

#EOF
42
