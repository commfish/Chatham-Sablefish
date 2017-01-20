# This is the age-structured sablegish model for Chatham Strait in Southeast Alaska 
# DIRECTORY CONTENTS
<hr>  
1.  **master**   
      The ADMB sablefish stock assessment model for Chatham Strait:  
      * sex-specific (distinct abundance, mortality, selectivity for each sex) 
      * uses NOAA longline survey selectivity inputs, male and female
      * calculates sex-specific fishery longline selectivity  
        (accounts for highgrading, which is illegal in the federal fishery but not in the State fishery)  
      * two recapture events for the mark-recapture estimates - fishery and survey  
        These use sex-specific capture rates and selectivities and estimate abundance at the *beginning* of the survey and fishery, respectively  
      * parametric bootstrap
      * mcmc output
2.  **Single-sex**   
      The ADMB sablefish stock assessment model for Chatham Strait:
      * Single-sex structure; no explicit sex parameters or derived quantities 
      * uses AVERAGE NOAA longline survey selectivity input, averaged between male and female
      * calculates a single fishery longline selectivity for both sexes
        (accounts for highgrading, which is illegal in the federal fishery but not in the State fishery)  
      * two recapture events for the mark-recapture estimates - fishery and survey  
        These use general capture rates and selectivities and estimate mean abundance at the *middle* of the survey and fishery, respectively  
        (Note that the survey is so short that the middle and beginning are the same)  
      * parametric bootstrap
      * mcmc output

#SUMMARY
The ADMB stock assessment model in this directory has been under development and not formally implemented to establish harvest rates or any other management quantity. It is, however, mature and ready for testing and implementation for the 2017 fishery season, barring any detected anomalies or problems with functionality.  

The current method of stock assessment in Chatham Strait consists of the following:  

1. Mark-recapture estimates of abundance;  
2. Partition abundance estimates into age-specific cohorts by application of commercial fishery age composition data
3. Increment abundances forward by a single year, with natural mortality = 0.1 and F from the previous year F_50
4. Set youngest age to previous amount (age 4 - assumption of constant recruitment. Not true, of course, but selectivity is so insigificant at that age it is essentially inert)
5. Calculate F_50 from current federal selectivity-at-age vectors and yield-per-recruit calcs
6. Apply F_50 * selectivity to projected abundance-at-age
7. Sum to obtain Allowable Biological Catch for the upcoming fishery  

NOTE: Where no mark-recapture has been conducted in a given year, the estimate of abundance is assumed unchanged from the previous M-R year.  

The age-structured model expands this essential core into a standar stock assessment structure and is intended to:  

a. Continue abundance trends in years for which no mark-recapture work was conducted instead of assuming constant abundance;  
b. Incorporate additional data into the estimate of abundance and relax the assumption that the mark-recapture estimate of abundance is estimated without error;  

# DATA SOURCES:  
1. Commercial fishery:  
    * Total annual catch
    * Age-composition 
    * CPUE  
    * Recapture data from tagged and released fish

2. ADF&G longline survey survey:  
    * Age-composition
    * CPUE
    * Recapture data from tagged and released fish  

# OBJECTIVE FUNCTION  
1. Likelihoods:  
    a. Total annual commercial catch   
    b. Mark-recapture estimates of abundance for both fishery and survey 
    c. Age composition for fishery and survey 
    e. CPUE for fishery and survey 

2. Penalties:  
    a. Annual recruitment deviations  
    b. Year 1 abundance deviations  
    e. Mean recruitment  
    f. Mean year 1 abundance  


#SPECIFIC MODEL NOTES:  

...under construction...
