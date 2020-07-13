#########################################################################################
### Marcel Kropp, 21.05.2018
### This script estimates VARs following the original contribution of
### Bloom (2009) in his paper entitled 'The Impact of Uncertainty Shocks' and
### corresponds to VAR8 mentioned in the main text of the thesis, meaning
### that all variables are detrended as indicated by Bloom (2009).

### Alternative VAR8-estimations without detrending will be estimated in another
### R-Script (following Jurado et al., 2015, based on Bloom, 2009).

### The data we use is an Excel-file which Nicholas Bloom himself makes available
### via his homepage as part of the available material dedicated to his
### 2009 paper;
### This Excel-file can be downloaded at 
### https://people.stanford.edu/nbloom/sites/default/files/r.zip and 
### comes along bundled together in a zip-file with an R-Script prepared by
### Jesus Fernandez-Villaverde (JFV) that replicates Bloom's (2009) VARs using
### slightly updated data and a longer time-range (the differences between
### Bloom's original used data and JFV's updated data are negligible and
### are mostly limited to the data for Industrial Production having slightly
### changed due to a change in the base-year between Bloom's contrubtion and
### JFV's replication);
### To stick to the original contribution from Bloom (2009), we will use
### JFV's updated data but restrict the range of the data to replicate
### Bloom (2009)'s results using 6 different uncertainty measures:
###         (PART 1) the Bloom-Shocks (marking months with maximum volatility)
###         (PART 2) the VXO/volatility series
###         (PART 3) the MSoC (available as of 1978)
###         (PART 4) the EPU (availbale as of 1985)
###         (PART 5) the GTU (available as of 2004)
###         (PART 6) the Macro Uncertainty Index, for h=1 and h=12 
###         (available as of 1960)
### After restricting the range of the data, the various time-series will
### run from 07/1962-06/2008 (like in Bloom's original contribution from 2009).

### The execution of the below relies on data and objects from 
### '12052018_thesis_master.R'; therefore this R-Script should be run 
### before running this one, to have all necessary objects available.


#########################################
### preparatory steps common to all VARs
#########################################
#####
## Read xlsx data and some basic manipulations
#####

var8_bloom_JFV <- read_excel("VARDATA_UPDATED_JFV.xlsx")

# Drop last observations in 2013 where no VOLATBL-data
# (i.e., VXO/volatility) is available using 'drop_na', 
# drop all rows that go beyond Bloom's original estimations
# (from 2009) by uniting YEAR and MONTH into the helper variable
# 'YEARMONTH',
# change the type to 'yearmon' to be able to filter accordingly,
# filter as we wish (i.e., YEARMONTH <= "Jun 2008")
# and at the end remove the helper-variable
# 'YEARMONTH' again:
var8_bloom_JFV <- var8_bloom_JFV %>% 
                        drop_na %>%
                        unite(YEARMONTH, YEAR, MONTH, remove = FALSE,
                              sep="-") %>%
                        mutate(YEARMONTH = as.yearmon(YEARMONTH)) %>%
                        filter(YEARMONTH <= "Jun 2008") %>%
                        dplyr::select(-YEARMONTH)


#####
## log- and/or HP-transformations
#####
# variables with log-and HP-transformations
cols_log_hp <- c("STOCK", "WAGE", "CPI", "EMPM", "IPM")
for (i in seq_along(cols_log_hp)) {
  var8_bloom_JFV[cols_log_hp[i]] <- log(var8_bloom_JFV[cols_log_hp[i]])
  var_aux <- hpfilter(var8_bloom_JFV[cols_log_hp[i]],freq=129600)
  var8_bloom_JFV[cols_log_hp[i]] <- var_aux$cycle
}

# variables only with HP-transformations
cols_hp <- c("HOURSM", "FFR")
for (i in seq_along(cols_hp)) {
  var_aux <- hpfilter(var8_bloom_JFV[cols_hp[i]],freq=129600)
  var8_bloom_JFV[cols_hp[i]] <- var_aux$cycle
}


## Below will follow six parts that will calculate all VARs and ultimately
## (after some transformations) store all VAR-results (i.e., oirfs including
## Upper and Lower bound of the confidence bands) in a dedicated data-frame.
## This one data-frame will then we feed into a ggplot using 'facets' to
## plot all impulse-response-fuctions we want to plot at once, i.e.,
## automatically split the plot into a matrix of panels (this will
## be Step 7 at the very end).

########################################################################
### PART 1: VAR8_Bloom_2009_HP with Bloom-Shocks (marking months
###         with maximum volatility)
########################################################################
# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_1 <- var8_bloom_JFV

####################
## transformations
####################
# Bloom_Shock: we want to define an indicator function on the 17 events defined
# in the original contribution by Bloom (2009)

# the below executes several steps at once:
# first, we change the name of the column,
# then we set the indicator variable to 1 for months where the
# volatility had reached its peak (according to Table A.1, Bloom, 2009, p. 35),
# in the last step we reorder variables (and only select the ones we
# want for the VARs) and can finally kick off the 
# computation of the VARs (and IRFs)
var8_bloom_JFV_1 <- var8_bloom_JFV_1 %>%
            dplyr::rename(bloom_shock = VOLATBL) %>%
            unite(YEARMONTH, YEAR, MONTH, remove = FALSE,
                  sep="-") %>%
            mutate(YEARMONTH = as.yearmon(YEARMONTH)) %>%
            mutate(bloom_shock = case_when(
              YEARMONTH %in% 
                as.yearmon(c("Oct 1962", #Cuban missile crisis
                              "Nov 1963", #Assassination of JFK
                              "Aug 1966", #Vietnam buildup
                              "Dec 1973", #OPEC I, Arab-Israeli War
                              "Oct 1974", #Franklin National Financial Crisis
                              "Nov 1978", #OPEC II
                              "Mar 1980", #Afghanistan, Iran hostages
                              "Oct 1982", #Monetary cycle turning point
                              "Nov 1987", #Black Monday
                              "Oct 1990", #Gulf War I
                              "Nov 1997", #Asian Crisis
                              "Sep 1998", #Russian, LTCM default
                              "Sep 2001", #9/11 terrorist attack
                              "Sep 2002", #Worldcom and Enron
                              "Feb 2003", #Gulf War II
                              "Mar 2008")) ~ 1, #Credit Crunch
                  TRUE                    ~ 0)) %>%
            # in a last step, we remove all date variables
            dplyr::select(-c(YEARMONTH, YEAR, MONTH)) %>%
            dplyr::select(STOCK, bloom_shock, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)
  
                        
# Check data
# View(var8_bloom_JFV_1)

####################
## VARs (VAR8)
####################
# estimate model (12 lags + constant) + compute impulse-responses:
# because we do not necessarily need the output from the 'var' - function
# as a separate object, we pipe 'var' and 'irf' together:


# var8_bloom_JFV_1.var <- var8_bloom_JFV_1 %>% VAR(p=12, type="const")
# resid <- residuals(var8_bloom_JFV_1.var)


# because we want to plot both the 68% and 95% CIs, we call the VAR-
# function twice and later merge the values for the 68% and 95% CIs together
var8_bloom_JFV_1.var.irfs.68 <- var8_bloom_JFV_1 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "bloom_shock"), 
                        impulse = "bloom_shock", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.68)
# for 


var8_bloom_JFV_1.var.irfs.95 <- var8_bloom_JFV_1 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "bloom_shock"), 
                        impulse = "bloom_shock", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.95)      

# Note: all results are now stored in var8_bloom_JFV_1.var.irfs.68
# and var8_bloom_JFV_1.var.irfs.95, respectively;
# running the below reveals the nested list-structure of the
# 'varirf' - object:
str(var8_bloom_JFV_1.var.irfs.68)

# first preliminary plots
# matplot(var8_bloom_JFV_1.var.irfs.68$irf$bloom_shock[, 3], type='l') 
# lines(var8_bloom_JFV_1.var.irfs.68$Lower$bloom_shock[, 3], col="blue")     
# lines(var8_bloom_JFV_1.var.irfs.68$Upper$bloom_shock[, 3], col="blue")
# lines(var8_bloom_JFV_1.var.irfs.95$Lower$bloom_shock[, 3], col="blue")     
# lines(var8_bloom_JFV_1.var.irfs.95$Upper$bloom_shock[, 3], col="blue") 

# we now need a smart way to extract the relevant components from the
# 'varirf' - object (which we have called 'var8_bloom_JFV_1.var.irfs.68'
# and 'var8_bloom_JFV_1.var.irfs.95'
# in our case), including the values for the oirfs, Lower- and Upper
# bounds of the confidence intervals:
# to do that, we first store the results in a data-frame/tibble
var8_bloom_JFV_1.var.irfs.df <- data.frame(
                                      var8_bloom_JFV_1.var.irfs.68$irf,
                                      var8_bloom_JFV_1.var.irfs.68$Lower,
                                      var8_bloom_JFV_1.var.irfs.68$Upper,
                                      var8_bloom_JFV_1.var.irfs.95$Lower,
                                      var8_bloom_JFV_1.var.irfs.95$Upper)

# and then restructure the data-frame in the way we need it:
# rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for bloom_shock.bloom_shock
var8_bloom_JFV_1.var.irfs.df <- var8_bloom_JFV_1.var.irfs.df %>%
                dplyr::rename(blo.blo.oirf = bloom_shock.bloom_shock, 
                              EMPM.oirf = bloom_shock.EMPM,
                              IPM.oirf = bloom_shock.IPM,
                              EMPM.Lo.68 = bloom_shock.EMPM.1,
                              IPM.Lo.68 = bloom_shock.IPM.1,
                              EMPM.Up.68 = bloom_shock.EMPM.2,
                              IPM.Up.68 = bloom_shock.IPM.2,
                              EMPM.Lo.95 = bloom_shock.EMPM.3,
                              IPM.Lo.95 = bloom_shock.IPM.3,
                              EMPM.Up.95 = bloom_shock.EMPM.4,
                              IPM.Up.95 = bloom_shock.IPM.4) %>%
                dplyr::select(-c(bloom_shock.bloom_shock.1, 
                                 bloom_shock.bloom_shock.2,
                                 bloom_shock.bloom_shock.3, 
                                 bloom_shock.bloom_shock.4))

# to normalize the orthogonalized irfs, we rescale them by dividing
# all values for response = EMPM/IPM by the first entry of
# bloom_shock on itself;
# finally, to rescale, we divide all columns in the data-frame (apart from the
# first one) by the first value of blo.blo.oirf (for n=1),
# and at the same time multiply with 100:

var8_bloom_JFV_1.var.irfs.df.rescaled <- 
            var8_bloom_JFV_1.var.irfs.df[, 2:ncol(var8_bloom_JFV_1.var.irfs.df)] /
            var8_bloom_JFV_1.var.irfs.df[1, 1]*100

# now we perform three 'gather' - commands inside a 'bind_cols'
# - command; the penulatimate command adds a step-counter and references to the
# entire previously created tibble to have the 'response' - variable
# available,
# the last command adds the name of the shock and references to the
# entire previously created tibble to have the 'shock' - variable
# available;

# because we will perform this operation multiple times, below is a function
# which we can hand over to a data-frame which will then apply the necessary 
# manipulations to arrive at the desired result (note that the below
# funtion only works on data-frame with the exact corresponding naming
# of the variables, etc):

create_irf_dataframe <- function(dataframe, shock_name){
                        (bind_cols(
                              dataframe %>% 
                                  gather('IPM.oirf', 
                                         'EMPM.oirf', 
                                         key="response", 
                                         value="oirf") %>%
                                  dplyr::select(response, oirf) %>%
                                  mutate(response = case_when(
                                        response == 'IPM.oirf' ~ 
                                          "% Impact on Production", 
                                        TRUE ~ "% Impact on Employment")),
                              dataframe %>%
                                  gather('IPM.Lo.68', 
                                         'EMPM.Lo.68', 
                                         key="response", 
                                         value="Lo.68") %>%
                                  dplyr::select(Lo.68), 
                              dataframe %>%
                                  gather('IPM.Up.68', 
                                         'EMPM.Up.68', 
                                         key="response", 
                                          value="Up.68") %>%
                                  dplyr::select(Up.68),
                              dataframe %>%
                                gather('IPM.Lo.95', 
                                       'EMPM.Lo.95', 
                                       key="response", 
                                       value="Lo.95") %>%
                                dplyr::select(Lo.95),
                              dataframe %>%
                                gather('IPM.Up.95', 
                                       'EMPM.Up.95', 
                                       key="response", 
                                       value="Up.95") %>%
                                dplyr::select(Up.95)
                                  ) %>%
                              group_by(response) %>% 
                                      mutate(step = row_number())) %>%
                              mutate(shock_name = shock_name)
}


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
var8_irfs_results_blo <- create_irf_dataframe(
                                    var8_bloom_JFV_1.var.irfs.df,
                                  'Bloom-Shock')
            

########################################################################
### PART 2: VAR8_Bloom_2009_HP with actual VXO/volatility series
########################################################################
# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_2 <- var8_bloom_JFV

####################
## transformations
####################
# we only select the variables we
# want for the VARs and can finally kick off the 
# computation of the VARs (and IRFs)
var8_bloom_JFV_2 <- var8_bloom_JFV_2 %>%
                        dplyr::select(STOCK, VOLATBL, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)

# testing an approach to make impulse-responses comparable:
# we standardize the uncertainty-measures:
# var8_bloom_JFV_2 <- var8_bloom_JFV_2 %>%
#                         dplyr::mutate(VOLATBL = 
#                                         (VOLATBL - mean(VOLATBL, 
#                                                         na.rm=TRUE))/
#                                         sd(VOLATBL, na.rm=TRUE))


####################
## VARs (VAR8)
####################
# estimate model (12 lags + constant) + compute impulse-responses:
# because we do not necessarily need the output from the 'var' - function
# as a separate object, we pipe 'var' and 'irf' together:

# the below is just a step in-between to calculate the standard
# deviation of the residuals from the uncertainty measure
var8_bloom_JFV_2.var <- var8_bloom_JFV_2 %>% VAR(p=12, type="const")
resid <- residuals(var8_bloom_JFV_2.var)


# because we want to plot both the 68% and 95% CIs, we call the VAR-
# function twice and later merge the values for the 68% and 95% CIs together
var8_bloom_JFV_2.var.irfs.68 <- var8_bloom_JFV_2 %>% VAR(p=12, type="const") %>%
                          irf(response = c("IPM", "EMPM", "VOLATBL"), 
                          impulse = "VOLATBL", n.ahead = 59, 
                          boot = TRUE, ortho = TRUE, runs=100, ci=0.68)
var8_bloom_JFV_2.var.irfs.95 <- var8_bloom_JFV_2 %>% VAR(p=12, type="const") %>%
                          irf(response = c("IPM", "EMPM", "VOLATBL"), 
                          impulse = "VOLATBL", n.ahead = 59, 
                          boot = TRUE, ortho = TRUE, runs=100, ci=0.95)  


# we first store the results in a data-frame/tibble
var8_bloom_JFV_2.var.irfs.df <- data.frame(
                          var8_bloom_JFV_2.var.irfs.68$irf,
                          var8_bloom_JFV_2.var.irfs.68$Lower,
                          var8_bloom_JFV_2.var.irfs.68$Upper,
                          var8_bloom_JFV_2.var.irfs.95$Lower,
                          var8_bloom_JFV_2.var.irfs.95$Upper)

# and then restructure the data-frame in the way we need it:
# rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for VOLATBL.VOLATBL
var8_bloom_JFV_2.var.irfs.df <- var8_bloom_JFV_2.var.irfs.df %>%
                dplyr::rename(VOLATBL.VOLATBL.oirf = VOLATBL.VOLATBL, 
                              EMPM.oirf = VOLATBL.EMPM,
                              IPM.oirf = VOLATBL.IPM,
                              EMPM.Lo.68 = VOLATBL.EMPM.1,
                              IPM.Lo.68 = VOLATBL.IPM.1,
                              EMPM.Up.68 = VOLATBL.EMPM.2,
                              IPM.Up.68 = VOLATBL.IPM.2,
                              EMPM.Lo.95 = VOLATBL.EMPM.3,
                              IPM.Lo.95 = VOLATBL.IPM.3,
                              EMPM.Up.95 = VOLATBL.EMPM.4,
                              IPM.Up.95 = VOLATBL.IPM.4) %>%                              
                dplyr::select(-c(VOLATBL.VOLATBL.1, 
                                 VOLATBL.VOLATBL.2,
                                 VOLATBL.VOLATBL.3,
                                 VOLATBL.VOLATBL.4))


# to normalize the orthogonalized irfs, we rescale them by dividing
# all values for response = EMPM/IPM by the first entry of
# bloom_shock on itself;
# finally, to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of blo.blo.oirf (for n=1),
# and at the same time multiply with 100:
var8_bloom_JFV_2.var.irfs.df.rescaled <- 
  var8_bloom_JFV_2.var.irfs.df[, 2:ncol(var8_bloom_JFV_2.var.irfs.df)] /
  var8_bloom_JFV_2.var.irfs.df[1, 1]*100

# lastly, we create the desired data-frame using the function
# 'create_irf_dataframe':
var8_irfs_results_vol <- create_irf_dataframe(
                                  var8_bloom_JFV_2.var.irfs.df,
                                  'VXO/volatility')



########################################################################
### PART 3: VAR8_Bloom_2009_HP with Michigan Survey Index
###         following Leduc and Lui (2016)
########################################################################

# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1
var8_bloom_JFV_3 <- var8_bloom_JFV

####################
## transformations
####################
# Merge the data with the data.frame containing the data for the
# Michigan Survey;
# the respective data comes from the data-frame: michigan_survey,
# and the variable of interest i 'vehrn_fb_all'
# first, we have to rename 'YEAR' and 'MONTH' to 'year' and 'month':
var8_bloom_JFV_3 <- var8_bloom_JFV_3 %>%
                    dplyr::rename(year = YEAR,
                                 month = MONTH)

# we perform an inner join to have the data for the michigan_survey in our
# data.frame; we perform an inner join to be only left with a dataset
# that ranges from the start of the availability of the MSoC-index
# until 2008/06; Note that the series for the MSoC starts in
# 02/1978 which is the year that the MSoC switched from monthly to
# quarterly data!
var8_bloom_JFV_3 <- inner_join(x = michigan_survey[, c("vehrn_fb_all",
                                                  "year", "month")], 
                         y = var8_bloom_JFV_3, 
                         by = c("year", "month"))

# because one data point in 01/1978 is missing, we have to remove that
# one row from our new data-frame (so after the below command, the data
# will start 02/1978):
var8_bloom_JFV_3 <- var8_bloom_JFV_3 %>%
                        drop_na


# we reorder variables (and only select the ones we
# want for the VARs)
var8_bloom_JFV_3 <- var8_bloom_JFV_3 %>%
                        dplyr::select(STOCK, vehrn_fb_all, 
                                      FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM) %>%
                        dplyr::rename(Michigan = 
                                        vehrn_fb_all)

# testing an approach to make impulse-responses comparable:
# we standardize the uncertainty-measures:
# var8_bloom_JFV_3 <- var8_bloom_JFV_3 %>%
#                         dplyr::mutate(Michigan = 
#                                         (Michigan - mean(Michigan, 
#                                                         na.rm=TRUE))/
#                                         sd(Michigan, na.rm=TRUE))

####################
## VARs (VAR8)
####################
# estimate model (12 lags + constant) + compute impulse-responses:
# because we do not necessarily need the output from the 'var' - function
# as a separate object, we pipe 'var' and 'irf' together:
# because we want to plot both the 68% and 95% CIs, we call the VAR-
# function twice and later merge the values for the 68% and 95% CIs together
var8_bloom_JFV_3.var.irfs.68 <- var8_bloom_JFV_3 %>% VAR(p=12, type="const") %>%
                         irf(response = c("IPM", "EMPM", "Michigan"), 
                         impulse = "Michigan", n.ahead = 59, 
                         boot = TRUE, ortho = TRUE, runs=100, ci=0.68)
var8_bloom_JFV_3.var.irfs.95 <- var8_bloom_JFV_3 %>% VAR(p=12, type="const") %>%
                         irf(response = c("IPM", "EMPM", "Michigan"), 
                         impulse = "Michigan", n.ahead = 59, 
                         boot = TRUE, ortho = TRUE, runs=100, ci=0.95) 


# we first store the results in data-frames/tibbles
var8_bloom_JFV_3.var.irfs.df <- data.frame(
                         var8_bloom_JFV_3.var.irfs.68$irf,
                         var8_bloom_JFV_3.var.irfs.68$Lower,
                         var8_bloom_JFV_3.var.irfs.68$Upper,
                         var8_bloom_JFV_3.var.irfs.95$Lower,
                         var8_bloom_JFV_3.var.irfs.95$Upper)


# and then restructure the data-frames in the way we need them:
# rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for Michigan.Michigan
var8_bloom_JFV_3.var.irfs.df <- var8_bloom_JFV_3.var.irfs.df %>%
                dplyr::rename(Michigan.Michigan.oirf = Michigan.Michigan, 
                              EMPM.oirf = Michigan.EMPM,
                              IPM.oirf = Michigan.IPM,
                              EMPM.Lo.68 = Michigan.EMPM.1,
                              IPM.Lo.68 = Michigan.IPM.1,
                              EMPM.Up.68 = Michigan.EMPM.2,
                              IPM.Up.68 = Michigan.IPM.2,
                              EMPM.Lo.95 = Michigan.EMPM.3,
                              IPM.Lo.95 = Michigan.IPM.3,
                              EMPM.Up.95 = Michigan.EMPM.4,
                              IPM.Up.95 = Michigan.IPM.4) %>%
                dplyr::select(-c(Michigan.Michigan.1, 
                                 Michigan.Michigan.2,
                                 Michigan.Michigan.3,
                                 Michigan.Michigan.4))

# to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of h1.h1.oirf (for n=1),
# and at the same time multiply with 100;
# as part of the normalization of the oirfs, we have to choose 
# meaningful units; for vol we choose a 15 unit shock
# which amounts to approx four standard deviations of the
# identified error;

var8_bloom_JFV_3.var.irfs.df.rescaled <-
  var8_bloom_JFV_3.var.irfs.df[, 2:ncol(var8_bloom_JFV_3.var.irfs.df)]/
  (var8_bloom_JFV_3.var.irfs.df[1, 1]*100)


# lastly, we create the desired data-frames using the function
# 'create_irf_dataframe':
var8_irfs_results_michigan <- create_irf_dataframe(
                                  var8_bloom_JFV_3.var.irfs.df,
                                  'Michigan')



########################################################################
### PART 4: VAR8_Bloom_2009_HP with Economic Policy Uncertainty
###         following Baker et al (2015)
########################################################################

# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_4 <- var8_bloom_JFV

####################
## transformations
####################
# Merge the data with the data.frame containing the data for the
# EPU;
# the respective data comes from the data-frame: epu_index;
# first, we have to rename 'YEAR' and 'MONTH' to 'year' and 'month':
var8_bloom_JFV_4 <- var8_bloom_JFV_4 %>%
                    dplyr::rename(year = YEAR,
                                 month = MONTH)
# next, we perform an inner join to have the data for the EPU in our
# data.frame; we perform an inner join to be only left with a dataset
# that ranges from the start of the availability of the EPU-index
# until 2008/06;
var8_bloom_JFV_4 <- inner_join(x = epu_index[, c("News_Based_Policy_Uncert_Index",
                                                  "year", "month")], 
                         y = var8_bloom_JFV_4, 
                         by = c("year", "month"))


# we reorder variables (and only select the ones we
# want for the VARs)
var8_bloom_JFV_4 <- var8_bloom_JFV_4 %>%
                        dplyr::select(STOCK, News_Based_Policy_Uncert_Index, 
                                      FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM) %>%
                        dplyr::rename(EPU = News_Based_Policy_Uncert_Index)



####################
## VARs (VAR8)
####################
# estimate model (12 lags + constant) + compute impulse-responses:
# because we do not necessarily need the output from the 'var' - function
# as a separate object, we pipe 'var' and 'irf' together:
# because we want to plot both the 68% and 95% CIs, we call the VAR-
# function twice and later merge the values for the 68% and 95% CIs together
var8_bloom_JFV_4.var.irfs.68 <- var8_bloom_JFV_4 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "EPU"), 
                        impulse = "EPU", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.68)
var8_bloom_JFV_4.var.irfs.95 <- var8_bloom_JFV_4 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "EPU"), 
                        impulse = "EPU", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.95) 


# we first store the results in data-frames/tibbles
var8_bloom_JFV_4.var.irfs.df <- data.frame(
                        var8_bloom_JFV_4.var.irfs.68$irf,
                        var8_bloom_JFV_4.var.irfs.68$Lower,
                        var8_bloom_JFV_4.var.irfs.68$Upper,
                        var8_bloom_JFV_4.var.irfs.95$Lower,
                        var8_bloom_JFV_4.var.irfs.95$Upper)

# and then restructure the data-frames in the way we need them:
# rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for VOLATBL.VOLATBL
var8_bloom_JFV_4.var.irfs.df <- var8_bloom_JFV_4.var.irfs.df %>%
                dplyr::rename(EPU.EPU.oirf = EPU.EPU, 
                              EMPM.oirf = EPU.EMPM,
                              IPM.oirf = EPU.IPM,
                              EMPM.Lo.68 = EPU.EMPM.1,
                              IPM.Lo.68 = EPU.IPM.1,
                              EMPM.Up.68 = EPU.EMPM.2,
                              IPM.Up.68 = EPU.IPM.2,
                              EMPM.Lo.95 = EPU.EMPM.3,
                              IPM.Lo.95 = EPU.IPM.3,
                              EMPM.Up.95 = EPU.EMPM.4,
                              IPM.Up.95 = EPU.IPM.4) %>%
                dplyr::select(-c(EPU.EPU.1, 
                                 EPU.EPU.2,
                                 EPU.EPU.3, 
                                 EPU.EPU.4))

# to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of h1.h1.oirf (for n=1),
# and at the same time multiply with 100;
# as part of the normalization of the oirfs, we have to choose 
# meaningful units; for vol we choose a 15 unit shock
# which amounts to approx four standard deviations of the
# identified error;

var8_bloom_JFV_4.var.irfs.df.rescaled <- 
  var8_bloom_JFV_4.var.irfs.df[, 2:ncol(var8_bloom_JFV_4.var.irfs.df)]/
  (var8_bloom_JFV_4.var.irfs.df[1, 1]*100)


# lastly, we create the desired data-frames using the function
# 'create_irf_dataframe':
var8_irfs_results_epu <- create_irf_dataframe(
                                  var8_bloom_JFV_4.var.irfs.df,
                                  'EPU')

########################################################################
### PART 5: VAR8_Bloom_2009_HP with Macro Uncertainty Index
###         following Jurado et al. (2015);
###         We specifically use the data for h=1 and h=12 to see how
###         they perform compared to each other
########################################################################

# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_6 <- var8_bloom_JFV

####################
## transformations
####################
# Merge the data with the data.frame containing the data for the
# Macro Uncertainty Index (for both h=1 and h=12);
# the respective data comes from the data-frame: macroUncertainty_index;
# first, we have to rename 'YEAR' and 'MONTH' to 'year' and 'month':
var8_bloom_JFV_6 <- var8_bloom_JFV_6 %>%
                    dplyr::rename(year = YEAR,
                                 month = MONTH)
# next, we perform a right-join to have the data for h1 and h12 in our
# data.frame, ranging from 1962-2008 (the same duration as the Bloom, 2009,
# estimations!)
var8_bloom_JFV_6 <- right_join(x = macroUncertainty_index[, c("h1", "h12", "year",
                                                              "month")], 
                         y = var8_bloom_JFV_6, 
                         by = c("year", "month"))


# we reorder variables (and only select the ones we
# want for the VARs)
var8_bloom_JFV_h1 <- var8_bloom_JFV_6 %>%
                        dplyr::select(STOCK, h1, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)
var8_bloom_JFV_h12 <- var8_bloom_JFV_6 %>%
                        dplyr::select(STOCK, h12, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)

####################
## VARs (VAR8) (note that we perform every step for h=1 and h=12)
####################
# estimate models (12 lags + constant) + compute impulse-responses:
# because we do not necessarily need the output from the 'var' - function
# as a separate object, we pipe 'var' and 'irf' together:

var8_bloom_JFV_h1.var <- var8_bloom_JFV_h1 %>% VAR(p=12, type="const")
resid <- residuals(var8_bloom_JFV_h1.var)
var8_bloom_JFV_h12.var <- var8_bloom_JFV_h12 %>% VAR(p=12, type="const")
resid <- residuals(var8_bloom_JFV_h12.var)

var8_bloom_JFV_h1.var.irfs.68 <- var8_bloom_JFV_h1 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "h1"), 
                        impulse = "h1", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.68)
var8_bloom_JFV_h1.var.irfs.95 <- var8_bloom_JFV_h1 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "h1"), 
                        impulse = "h1", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.95)
var8_bloom_JFV_h12.var.irfs.68 <- var8_bloom_JFV_h12 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "h12"), 
                        impulse = "h12", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.68)
var8_bloom_JFV_h12.var.irfs.95 <- var8_bloom_JFV_h12 %>% VAR(p=12, type="const") %>%
                        irf(response = c("IPM", "EMPM", "h12"), 
                        impulse = "h12", n.ahead = 59, 
                        boot = TRUE, ortho = TRUE, runs=100, ci=0.95)

# we first store the results in data-frames/tibbles
var8_bloom_JFV_h1.var.irfs.df <- data.frame(
                        var8_bloom_JFV_h1.var.irfs.68$irf,
                        var8_bloom_JFV_h1.var.irfs.68$Lower,
                        var8_bloom_JFV_h1.var.irfs.68$Upper,
                        var8_bloom_JFV_h1.var.irfs.95$Lower,
                        var8_bloom_JFV_h1.var.irfs.95$Upper)

var8_bloom_JFV_h12.var.irfs.df <- data.frame(
                        var8_bloom_JFV_h12.var.irfs.68$irf,
                        var8_bloom_JFV_h12.var.irfs.68$Lower,
                        var8_bloom_JFV_h12.var.irfs.68$Upper,
                        var8_bloom_JFV_h12.var.irfs.95$Lower,
                        var8_bloom_JFV_h12.var.irfs.95$Upper)

# and then restructure the data-frames in the way we need them:
# rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for VOLATBL.VOLATBL
var8_bloom_JFV_h1.var.irfs.df <- var8_bloom_JFV_h1.var.irfs.df %>%
                dplyr::rename(h1.h1.oirf = h1.h1, 
                              EMPM.oirf = h1.EMPM,
                              IPM.oirf = h1.IPM,
                              EMPM.Lo.68 = h1.EMPM.1,
                              IPM.Lo.68 = h1.IPM.1,
                              EMPM.Up.68 = h1.EMPM.2,
                              IPM.Up.68 = h1.IPM.2,
                              EMPM.Lo.95 = h1.EMPM.3,
                              IPM.Lo.95 = h1.IPM.3,
                              EMPM.Up.95 = h1.EMPM.4,
                              IPM.Up.95 = h1.IPM.4) %>%
                dplyr::select(-c(h1.h1.1, 
                                 h1.h1.2,
                                 h1.h1.3,
                                 h1.h1.4))

var8_bloom_JFV_h12.var.irfs.df <- var8_bloom_JFV_h12.var.irfs.df %>%
                dplyr::rename(h12.h12.oirf = h12.h12, 
                              EMPM.oirf = h12.EMPM,
                              IPM.oirf = h12.IPM,
                              EMPM.Lo.68 = h12.EMPM.1,
                              IPM.Lo.68 = h12.IPM.1,
                              EMPM.Up.68 = h12.EMPM.2,
                              IPM.Up.68 = h12.IPM.2,
                              EMPM.Lo.95 = h12.EMPM.3,
                              IPM.Lo.95 = h12.IPM.3,
                              EMPM.Up.95 = h12.EMPM.4,
                              IPM.Up.95 = h12.IPM.4) %>%
                dplyr::select(-c(h12.h12.1, 
                                 h12.h12.2,
                                 h12.h12.3, 
                                 h12.h12.4))


# to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of h1.h1.oirf (for n=1),
# and at the same time multiply with 100;
# as part of the normalization of the oirfs, we have to choose 
# meaningful units; for vol we choose a 15 unit shock
# which amounts to approx four standard deviations of the
# identified error;

var8_bloom_JFV_h1.var.irfs.df.rescaled <- 
  var8_bloom_JFV_h1.var.irfs.df[, 2:ncol(var8_bloom_JFV_h1.var.irfs.df)]/
  (var8_bloom_JFV_h1.var.irfs.df[1, 1]*100)
var8_bloom_JFV_h12.var.irfs.df.rescaled <- 
  var8_bloom_JFV_h12.var.irfs.df[, 2:ncol(var8_bloom_JFV_h12.var.irfs.df)]/
  (var8_bloom_JFV_h12.var.irfs.df[1, 1]*100)

# lastly, we create the desired data-frames using the function
# 'create_irf_dataframe':
var8_irfs_results_h1 <- create_irf_dataframe(
                                  var8_bloom_JFV_h1.var.irfs.df,
                                  'macro: h1')
var8_irfs_results_h12 <- create_irf_dataframe(
                                  var8_bloom_JFV_h12.var.irfs.df,
                                  'macro: h12')



####################
## plots of orifs
####################
# before plotting, we have to stack the vairous data-frames that we have
# created, together:
var8_irfs_results_all <- bind_rows(var8_irfs_results_blo, var8_irfs_results_vol,
                                  var8_irfs_results_michigan, var8_irfs_results_epu, 
                                  var8_irfs_results_h1, 
                                  var8_irfs_results_h12
                                  )
# we need to create ordered factors, otherwise 'facet_wrap' combines the plots
# in alphabetical order:
var8_irfs_results_all <- var8_irfs_results_all %>%
                          ungroup %>%
                          mutate(shock_name = factor(shock_name, 
                                         ordered = TRUE,
                                         levels=unique(shock_name)),
                                 response = factor(response,
                                          ordered = TRUE,
                                         levels=unique(response)))

# now we can finally plot the irfs:
var8_plot_all_HP_until2008 <- ggplot(data = var8_irfs_results_all) +
  # geom_point() + 
  # geom_line(aes(x = step, y = Up), color="#e80628", 
  #           size=0.5,linetype = 3) +
  # geom_line(aes(x = step, y = Lo), color="#e80628", 
  #           size=0.5, linetype = 3) +
  # geom_ribbon(aes(x = step, ymax=Up, ymin=Lo), 
  #             fill="#706c6c", alpha=.3) +
  # below ribbon is for 68% confidence band
  geom_ribbon(aes(x = step, ymax=Up.68, ymin=Lo.68, fill=shock_name),
              alpha=0.6) +
  # and this one for the 95% confidence band
  geom_ribbon(aes(x = step, ymax=Up.95, ymin=Lo.95, fill=shock_name),
              alpha=.3) +
  # geom_line(aes(x = step, y = oirf), color="black", 
  #           size=0.8) +
  geom_line(aes(x = step, y = oirf), color="black", 
            size=0.5) +
  # geom_point(aes(x = step, y = oirf), color="black", 
  #            size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 60), 
                     breaks = seq(0, 60, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = NULL,  
                     limits = c(-0.008, 0.004),
                     breaks = seq(-0.008, 0.004, by=0.002),
                     minor_breaks = NULL) +
  #ggtitle("% impact on industrial production") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=0.5) +
  theme(axis.text=element_text(size=10),
        plot.title = element_text(size=10, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.position="none",
        #legend.text=element_text(size=14),
        #axis.text.x=element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="black", fill="white", 
                                       size=1.5, linetype="solid")) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 5) + facet_grid(shock_name ~ response,
                                      scales="free_y")

var8_plot_all_HP_until2008

ggsave(file="var8_plot_all_HP_until2008.pdf", width = 210, height = 297, units = "mm")

