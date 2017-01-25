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
    
3. Marked fish tagged and released in the pot survey for recapture in the survey and fishery  

# OBJECTIVE FUNCTION  
1. Likelihoods:  
    a. Total annual commercial catch   
    b. Mark-recapture estimates of abundance for both fishery and survey   
    c. Age composition for fishery and survey   
    d. CPUE for fishery and survey   

2. Penalties:  
    a. Annual recruitment deviations  
    b. Year 1 abundance deviations  
    c. Mean recruitment  
    d. Mean year 1 abundance  


#SPECIFIC MODEL NOTES:  

1. Much of the code, especially in the objective function and ABC projections, is written out long-hand; it could be streamlined by using 'dnorm' calls, removing explicit loops, etc. I left it this way to faciliate my own error checking as well as implementation of the retrospective analysis code. Although that is all now resolved, I don't have the time to go through and clean it all up. By all means, feel free to clean this up once you are comfortable with the code workings;  

2. The sample size for age composition has been set to the square root of the sample size scaled to 100 across all sample years. Running the model, however, implements a MANUAL sample size correction: the variances of the residuals in the age-composition data are compared to the assumed input sample sizes by assessing the standard deviation of normalized residuals (Breen at el. 2003). Under the assumption that the normalized residuals are normally distributed, the acceptable limit for SDNR values, following Francis (2011) is given as  

$max(sdnr)<[\chi_{0.95}^2/(m-1)]^{0.5}$

for which m = the number of years in the age-composition data set. If the SNDR for any given year exceeds this limit for a given data set, the input age composition sample size vector is divided by the SDNR vector, and the model re-run with the revised sample sizes. The process is iteratively repeated until the target maximum SDNR value is reached. **This means that the model is initially run using the square root of the sample size as output by the R script, the resulting SDNR examined (in the model.rep file), and if any single element of the vector exceeds the limit, the sample sizes are divided by the SDNR values, AND THE RESULTING REVISED SAMPLE SIZE ARE MANUALLY INSERTED IN THE MODEL.DAT FILE, REPLACING THE SAMPLE SIZE VECTOR (COMMENT THAT LINE OUT).** Using the square root of the sample size scaled to 100 meets this SDNR criteria for several data sets on the first model run, while others needs the manual scaling.  

    Sex-specific model: SDNR limit was met by fishery age composition for both males and females on model run. A single manual iteration was required to bring longline survey age compositions below the limit for the 2015 model.  
    
    Sex-unified model: SDNR limit was met by fishery age composition. Two manual iterations were required to bring the longline survey SDNR below the limit for the 2015 model.   
    
3. The mark-recapture estimates of abundance will differ between the two models. The procedure for determining the number of marks available for recapture and the subsequent estimates of abundance based on the number of recaptures is:

    a. The pot survey catches, tags, and releases a given number of fish in May. Each tagged fish has its length recorded;  
    b. Both models use the NOAA-determined longline survey gear selectivity (as does the traditional stock assessment method);  
    c. The number of tags-at-length are scaled by the gear selectivity at that length. This produces the number of tags available to be recovered by the GEAR in both the survey and fishery, as the gear is identical;  
    d. For the sex-unified model, the scaling selectivity vector is the average between the male and female gear selectivity from NOAA (transformed from selectivity-at-age to selectivity-at-length by the proper R script). For the sex-specific model, the number of pot tags is multiplied by the proportion of males and females tagged and then scaled by the gender-specific selectivity vector;  
    e. Step D will produce two different estimates of the number of tags available to be recovered in a given year;  
    f. Total abundance is calculated as a maximum likelihood fit to a binomial probability:  
        
        1. **Sex-specific model - commercial fishery**  
            a. For each year, the number of recovered tags is reported by date and port surveyed;
            b. Robson Regier (1964) used to define the number of marks necessary for desired precision relative to the total number available (A = 0.25, 1 - alpha = 0.95);  
            c. The recovered tags are partitioned into distinct periods in which each period contains the minimum number of marks. Note - the final period will likely contain fewer than the target; - each period will likely contain more than the minimum number, as dates and vessels are kept intact and tags not divided between them;  
            d. Abundance is estimated as the number of fish available at the BEGINNING of the fishery, with each recapture event (i.e. minimum number of marks partitioned as per c above) tracking the population reduction due to fishing mortality and natural mortality, and each recapture event producing its own estimate of abundance as:
            
             DATA ________________________________________________________________________  
	     M_0   - total number of marked fish AVAILABLE to fishery  
  
	     D     - total number of marks recovered prior to the longline fishery
	             from OUTSIDE Chatham Strait, *NOT* including the longline survey
	             Spreadsheet from Mike Vaughn in Sitka  
	    
	     m_s   - total number of marks recovered in the longline survey  
	   
	     p     - number of longline fishery periods (looped as 'i')
	             NOTE: you must explicitly call this variable not only in the
	                 mr_foo = f(MARK_RECAPTURE...) call, but you must also
	                 explicitly DEFINE it in the function itself. Otherwise,
	                 the function will not loop correctly  
	   
	     mort  - natural mortality. Set to 0.1 and incremented daily (0.1/365)  
	     
	     p_male_p    - % male in pot survey  
	     p_female_p  - % female in pot survey  
	     p_male_s    - % male in longline survey  
	     p_female_s  - % female in longline survey  
	     p_male_f    - % male in longline fishery  
	     p_female_f  - % female in longline fishery  
	    
	     t_i   - time in period i since last recovery event  
	     C_i   - number of recaptured fish examined for marks in period i  
	     m_i   - number of marked individuals in period i  
	     d_i   - number of marks recovered OUTSIDE of Chatham Strait in period i  
	     n_i   - number of fish examined for marks in period i  
	    
	     CALCULATIONS_________________________________________________________________  
	     M_i   - number of marked individuals available for recovery in period i  
	             M_i+1 = [(M_i - m_i)*exp(mort * t_i)] - d_i  
	             M_1   = M_0 * exp(mort * t_1) - D - m_s  
	    
	     r_i   - ratio of marked to examined fish in period i   
	             (m_i / n_i)  
	    
	     pf_i  - probability of recapture for females in period i   
	             (p_female_p * M_i) / Nf_i  
	    
	     pm_i  - probability of recapture for males in period i   
	             (p_male_p * M_i) / Nm_i  
	    
	     ESTIMATED QUANTITIES_________________________________________________________  
	     N     - total estimated abundance (both sexes combined) at the **BEGINNING**
	             of the commercial longline fishery
	             THIS IS YOUR MARK-RECAPTURE ESTIMATE OF ABUNDANCE
	    
	     Nf_i  - estimated abundance of unmarked females in period i  
	             Nf_1 = N * 0.5  
	      
	     Nm_i  - estimated abundance of unmarked males in period i  
	             Nm_1 = N * 0.5  
	    
	     LIKELIHOODS__________________________________________________________________  
	     bin_f - binomial distribution fit to data  
	             BINOMIAL ( m_i * p_female_f | n_i * p_female_f, pf_i)  
	    
	     bin_m - binomial distribution fit to data  
	             BINOMIAL ( m_i * p_male_f | n_i * p_male_f, pm_i)    
	             
	2. **Sex-specific model - longline survey**  
	    a. As above, but with only a single recovery period (the survey is completed within a few days) and the appropriate modifications from fishery to survey data  
	    
	3. **Sex-unified model - commercial fishery**  
	    a. The above calculations are collapsed into a single set of data: number available for recovery, total number recovered, total number recovered outside the fishery or prior to the fishery, etc.  
	    b. The total abundance N calculated here is the AVERAGE abundance over the course of the commercial fishery, meaning abundance at the mid-point of the fishery. Examination of the model code will show that the sex-unified model is fit to abundance half-way through the fishery, whereas the mark-recapture abundance in the sex-specific model is fit to the abundance at the beginning of the fishery;  
	    
	4. **Sex-unified model - longline survey **
	    a. As the sex-specific model, although without sex-specific properties and using the total marks available to the sex-unified structure.  
