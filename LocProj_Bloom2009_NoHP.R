#########################################################################################
### Marcel Kropp, 26.05.2018
### This script estimates Impulse Responses using Jordá's (2005) Local Projection
### Method; 
### range of data: 07/1962 - 06/2009 (following Bloom 2009)
### all variables detrended apart from uncertainty measures (following Bloom 2009)

### Alternative Local Projections without detrending will be estimated in another
### R-Script called ('LocProj_Bloom2009_NoHP)


### While each PART of the below estimations (i.e., whenever we change the
### uncertainty measure in our data-set), we could have refered back to the 
### data.frames we have created for our VAR-estimations,
### we decided to set up a dedicated data.frame for the local projections
### which contains all variables and from which we only select the 
### variables we need one by one;
### the reason for this simplification is that, as opposed to the case 
### of VARs, the variables do not have to enter the regressions in a 
### certain order in order to obey to the recursive scheme as postulated
### via the cholesky decomposition!

### Note: In a revision of this script we should account for the larger
### data-consumption of Local Projections (along leads and lags)!
### This is because, as mentioned in Luca Brugnolini (2018), increasing
### the hoirzons of impulse responses reduces the sample available for the
### estimation itself. 
### But while VAR consumes data only along the lag dimension (p), local
### projections consume data along both the lag(p) and the lead(h) dimensions.

### If our goal is to calculate maximum h=60 steps-ahread forecasts, that 
### amounts to five years.
### Considering that 


## THIS DATA-PREPARATION WILL ONLY BECOME NECESSARY ONCE WE REPLACE
## THE DATA FROM VILLAVERDE WITH OUR OWN DATA!
##--------------------------------------------------------------------------------
##       loading additional_vars.csv;
##       merge with variables from sp500_merge_vxo
##       (most of the data comes from FRED 
##       [https://fred.stlouisfed.org/series/])
##       the employment in manufacturing comes
##       from BLS (see https://data.bls.gov/timeseries/CES3000000001)
##--------------------------------------------------------------------------------

# first we read in the file 'additional_vars.csv'
# additional_vars <- read.csv(file="additional_vars.csv", header=TRUE, sep=",",
#                            stringsAsFactors = FALSE)
# additional_vars <- read.csv(file="additional_vars_ext.csv", header=TRUE, sep=",",
#                            stringsAsFactors = FALSE)
additional_vars <- read.csv(file="additional_vars_ext_long.csv", header=TRUE, sep=",",
                            stringsAsFactors = FALSE)

# next we rename the column names and convert the column 'DATE' to 
# 'year', 'month' and then a numeric 'my' (i.e., the year-month in numeric 
# format)


additional_vars <- rename(additional_vars, 
                          date=DATE, ip=INDPRO, sp500=Close, ffr=FEDFUNDS)

additional_vars$date <- as.Date(additional_vars$date)
# note that the above command converts the date-type to a character again
# (we will deal with this later!)

# next, we create the three variable 'month', 'year' and 'day' using
# the 'separate()' - function from the 'tidyr' - package
additional_vars <- separate(additional_vars, "date", c("year", "month", "day"), sep = "-", 
                            remove=FALSE, convert=TRUE)

# next we create the variable 'my' which is a numerical representation
# of yearmon:
additional_vars <- as.data.frame(additional_vars %>%
                                   mutate(my = year + month/12))
# convert ip to a numeric data-type
additional_vars$ip <- as.numeric(additional_vars$ip)

# and drop a few superfluous columns:
additional_vars <- as.data.frame(additional_vars %>%
                                   dplyr::select(-c( 
                                     Date, Open, 
                                     High, Low, 
                                     Adj.Close, Volume)))


# the variable EMPM has to be retrieved from a separate dataset
# which we downloaded from the homepage of the BLS;
# (see https://data.bls.gov/timeseries/CES3000000001)
# we have to manipulate the data to bring it into our desired format:
# EMPM <- read.xlsx("BLS_EMPM.xlsx", header = TRUE, sheetIndex = 1)
EMPM <- read.xlsx("BLS_EMPM_long.xlsx", header = TRUE, sheetIndex = 1)

# to transform our above data-frame into one column called 'empm', we
# use the 'gather()' - function
# (see https://uc-r.github.io/tidyr for details about how it works):
# In particular, we want to transform our dataframe from wide (where each
# month represents a variable) to gather each month within one
# column and also gather the values associated with each month in 
# a second variable called 'value':
long_EMPM <- EMPM %>% gather(month, value, Jan:Dec)

# next, we transform 'month' to the type numeric
long_EMPM$month <- match(long_EMPM$month, month.abb)

# we rename the 'year' - column:
long_EMPM <- rename(long_EMPM, year = Year, empm = value)

# and create a variable 'my' which we will use to join this variable
# with our data-frames below:
long_EMPM <- as.data.frame(long_EMPM %>%
                             mutate(my = year + month/12))

# we re-order by this newly created 'my' - variable
long_EMPM <- long_EMPM %>% arrange(my)

# next, we join this variable into our data-frame called 'additional_vars':
additional_vars <- inner_join(x = long_EMPM, y = additional_vars, 
                              by = c("my", "year", "month"))
# note that we use a few months of data because the variable 'industrial production'
# in the data.frame 'additional_vars' is only available with a few months of delay
# and hence currently is only availalbe back until Feb/2018.

# our data-frame sp500_merge_vxo contains many superfluous variables
# while the only goal was to construct a consistent volatility - measure
# (measuring uncertainty) and the derivation of the 'Bloom-shocks' 
# thereafter;

# Therefore we create a subset of sp500_merge_vxo which will hold
# our variables of interest with which we will merge 'additional_vars':
# our subset consists of the following variables:
# 'year', 'month', 'my', 'bloom_shock', 'mvol'
sp500_vxo_subset <- as.data.frame(sp500_merge_vxo %>%
                                    dplyr::select(year, month, my, bloom_shock, mvol))


# having sp500_vxo_subset at our disposal, we can now join
# 'additional_vars' and 'sp500_vxo_subset' based on 'my'

regressions <- inner_join(x = sp500_vxo_subset, y = additional_vars, 
                          by = c("my", "year", "month"))


# by specifying an 'inner_join', we make sure to only take months
# where we have data in both datasets that we merged;
# the resulting data-frame 'regressions' can now be used for our
# local projections below!

# we generate a variable 'h' which rungs from 0 to the number of rows-1
# and which resembles the horizon that we are currently looking at
# (when running the Jordá local projections):
regressions <- as.data.frame(regressions %>%
                               mutate(h = seq(0, nrow(regressions)-1)))


# in particular, the variable 'h' will be used in our below loop to know
# exactly where (i.e., in which row) we want to place the coefficients from
# the regressions that we run!
##--------------------------------------------------------------------------------


#########################################
### preparatory steps common to all Local Projections
#########################################
#####
## Read xlsx data and some basic manipulations
#####

# Note: the below initiated data.frame will hold all variables
# that we will use for the stimations below;
# therefore, after reading in the data from the .xlsx-file
# 'VARDATA_UPDATED_JFV.xlsx' as prepared by Jesus Villaverde
# (note: this step will be replaced in a future version of the script
# when we use our own created data.frame!),
# we will immediately merge the data-set with the various data-frames
# holding the uncertainty measues!

locproj8_bloom_JFV <- read_excel("VARDATA_UPDATED_JFV.xlsx")

# we rename the variables 'YEAR' and 'MONTH' to 'year' and
# 'month' to be able to use them below when merging with the
# data.frame containing the uncertainty-measures!
locproj8_bloom_JFV <- locproj8_bloom_JFV %>%
  dplyr::rename(year = YEAR,
                month = MONTH)

# note: the data.frame 'identify_shocks' already contains
# all uncertainty measures we want to merge in here;
# therefore a single right_join is sufficient to have all
# uncertainty measures available in our data.frame
# 'locproj8_bloom_JFV':
locproj8_bloom_JFV <- as.tibble(
  right_join(x = identify_shocks[, 
                                 c("VXO", "EPU", "EPU_Hist", "Michigan",
                                   "Macro1", "Macro12","year",
                                   "month")],
             y = locproj8_bloom_JFV,
             by = c("year", "month"))) %>%
  # we add further steps to the newly created
  # tibble:
  # Bloom_Shock: we want to define an indicator function 
  # on the 17 events defined
  # in the original contribution by Bloom (2009)
  # first, we change the name of the column VOLATBL,
  # then we set the indicator variable to 1 for months 
  # where the volatility had reached its peak (according to 
  # Table A.1, Bloom, 2009, p. 35)
  dplyr::rename(bloom_shock = VOLATBL) %>%
  unite(YEARMONTH, year, month, remove = FALSE,
        sep="-") %>%
  # change the type to 'yearmon' to be able to
  # replace the Bloom-dates accordingly!
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
  # in a last step, we remove all date variables;
  # note that we do not yet remove the variable
  # 'YEARMONTH' but will only do this below
  # in the next stage!
  dplyr::select(-c(year, month))




# Note: in the current setting, the resulting data.frame will range 
# until 05/2013 which is the maximum length of data contained
# in the .xlsx-file "VARDATA_UPDATED_JFV.xlsx" from 
# Jesus Villaverde while the uncertainty-measures in
# 'identify_shocks' all range much longer;
# this will change once we replace "VARDATA_UPDATED_JFV.xlsx" 
# with out own data.frame!

# Note that in the resulting data.frame, we have two series for the
# VXO-data: VXO and VOLATBL;
# VXO is the one we'll reserve for the actual raw VXO series, VOLATBL
# will be overwritten with the bloom-shocks below!

# Drop last observation in 2013 using 'drop_na', 
# for the moment we DO NOT drop all rows that go beyond Bloom's 
# original estimations (from 2009) by uniting YEAR and MONTH into 
# the helper variable 'YEARMONTH' because we still have to figure
# out the estimation window in order to produce results that
# match our VARs!
#  to be able to filter accordingly,
# filter as we wish (i.e., YEARMONTH <= "Jun 2008")
# and at the end remove the helper-variable
# 'YEARMONTH' again:
locproj8_bloom_JFV <- locproj8_bloom_JFV %>% 
  drop_na %>%
  # note that we have muted the filtering
  # for the moment
  # filter(YEARMONTH <= "Jun 2008") %>%
  dplyr::select(-YEARMONTH)

#####
## log- and/or HP-transformations
#####
# variables only with log-transformations
cols_log <- c("STOCK", "WAGE", "CPI", "EMPM", "IPM")
for (i in seq_along(cols_log)) {
  locproj8_bloom_JFV[cols_log[i]] <- log(locproj8_bloom_JFV[cols_log[i]])
}

# no variables with HP-transformations!



# now we perform three 'gather' - commands inside a 'bind_cols'
# - command; the penulatimate command adds a step-counter and references to the
# entire previously created tibble to have the 'response' - variable
# available,
# the last command adds the name of the shock and references to the
# entire previously created tibble to have the 'shock' - variable
# available;



## Below will follow six parts that will calculate all
## local projections and ultimately store the point estimates including
## Upper and Lower bounds of the confidence bands in a dedicated 
## data.frame.
## This one data-frame will then be fed into ggplot using 'facets' to
## plot all impulse-response-fuctions which we have generated via local 
## projections at once, i.e.,
## automatically split the plot into a matrix of panels (this will
## be Step 7 at the very end).



########################################################################
### PART 1: LocProj_Bloom2009_HP with Bloom-Shocks (marking months
###         with maximum volatility)
########################################################################
# as compared to our VAR-estimations, we do not create separate
# data.frames for each estimation; therefore we can jump right ahead
# to the estimations

####################
## transformations
####################

# we only select the variables that are actually necessary for this
# step:
locproj8_bloom_JFV.1 <- locproj8_bloom_JFV %>%
  dplyr::select(c(bloom_shock, IPM, EMPM, HOURSM, CPI, WAGE, FFR,
                  STOCK))

# we want to check whether the procedure below indeed runs all regressions
# correctly; therefore, we create the 12 lags for all 8 variables manually!
locproj8_bloom_JFV.1 <- locproj8_bloom_JFV.1 %>%
  do(setNames(data.frame(., data.table::shift(.$IPM, type="lag", 1:12), 
                         data.table::shift(.$EMPM, type="lag", 1:12),
                         data.table::shift(.$HOURSM, type="lag", 1:12),
                         data.table::shift(.$CPI, type="lag", 1:12),
                         data.table::shift(.$WAGE, type="lag", 1:12),
                         data.table::shift(.$FFR, type="lag", 1:12),
                         data.table::shift(.$STOCK, type="lag", 1:12),
                         data.table::shift(.$bloom_shock, type="lag", 1:12)), 
              c(names(locproj8_bloom_JFV.1), paste0(rep(c("IPM_lag", 
                                                          "EMPM_lag",
                                                          "HOURSM_lag",
                                                          "CPI_lag",
                                                          "WAGE_lag",
                                                          "FFR_lag",
                                                          "STOCK_lag",
                                                          "bloom_shock"), 
                                                        each = 12), 1:12) )))


####################
## Local Projections (with 8 variables)
####################
# our goal now is to loop i from 0 to 60 (i.e., the horizons that we 
# want to consider);
# for the regressions we use 'dynlm' to have easy access to LAG-operators
# within the regressions!



# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.1 %>%
    dplyr::mutate(IPM = data.table::shift(IPM, type="lead", i))
  
  reg <- lm(IPM ~ ., df.reg[, c(1:2, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}


regressions_out_ipm <- regressions_out

# rename the column names
colnames(regressions_out_ipm)[1] <- "IPM.locproj"
colnames(regressions_out_ipm)[2] <- "IPM.Lo.95"
colnames(regressions_out_ipm)[3] <- "IPM.Up.95"
colnames(regressions_out_ipm)[4] <- "IPM.Lo.68"
colnames(regressions_out_ipm)[5] <- "IPM.Up.68"  


# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)

for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.1 %>%
    dplyr::mutate(EMPM = data.table::shift(EMPM, type="lead", i))
  
  reg <- lm(EMPM ~ ., df.reg[, c(1, 3, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}      

regressions_out_empm <- regressions_out
# rename the column names
colnames(regressions_out_empm)[1] <- "EMPM.locproj"
colnames(regressions_out_empm)[2] <- "EMPM.Lo.95"
colnames(regressions_out_empm)[3] <- "EMPM.Up.95"
colnames(regressions_out_empm)[4] <- "EMPM.Lo.68"
colnames(regressions_out_empm)[5] <- "EMPM.Up.68" 

regressions_out_bloom <- data.frame(regressions_out_empm,
                                    regressions_out_ipm)


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
locproj_results_bloom <- create_irf_dataframe_locproj(
  regressions_out_bloom,
  'Bloom-Shock')


########################################################################
### PART 2: LocProj_Bloom2009_HP with actual VXO/volatility series
########################################################################

# we only select the variables that are actually necessary for this
# step:
locproj8_bloom_JFV.2 <- locproj8_bloom_JFV %>%
  dplyr::select(c(VXO, IPM, EMPM, HOURSM, CPI, WAGE, FFR,
                  STOCK))

# we want to check whether the procedure below indeed runs all regressions
# correctly; therefore, we create the 12 lags for all 8 variables manually!
locproj8_bloom_JFV.2 <- locproj8_bloom_JFV.2 %>%
  do(setNames(data.frame(., data.table::shift(.$IPM, type="lag", 1:12), 
                         data.table::shift(.$EMPM, type="lag", 1:12),
                         data.table::shift(.$HOURSM, type="lag", 1:12),
                         data.table::shift(.$CPI, type="lag", 1:12),
                         data.table::shift(.$WAGE, type="lag", 1:12),
                         data.table::shift(.$FFR, type="lag", 1:12),
                         data.table::shift(.$STOCK, type="lag", 1:12),
                         data.table::shift(.$VXO, type="lag", 1:12)), 
              c(names(locproj8_bloom_JFV.2), paste0(rep(c("IPM_lag", 
                                                          "EMPM_lag",
                                                          "HOURSM_lag",
                                                          "CPI_lag",
                                                          "WAGE_lag",
                                                          "FFR_lag",
                                                          "STOCK_lag",
                                                          "VXO"), 
                                                        each = 12), 1:12) )))


####################
## Local Projections (with 8 variables)
####################
# our goal now is to loop i from 0 to 60 (i.e., the horizons that we 
# want to consider);
# for the regressions we use 'dynlm' to have easy access to LAG-operators
# within the regressions!



# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.2 %>%
    dplyr::mutate(IPM = data.table::shift(IPM, type="lead", i))
  
  reg <- lm(IPM ~ ., df.reg[, c(1:2, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}


regressions_out_ipm <- regressions_out

# rename the column names
colnames(regressions_out_ipm)[1] <- "IPM.locproj"
colnames(regressions_out_ipm)[2] <- "IPM.Lo.95"
colnames(regressions_out_ipm)[3] <- "IPM.Up.95"
colnames(regressions_out_ipm)[4] <- "IPM.Lo.68"
colnames(regressions_out_ipm)[5] <- "IPM.Up.68"  


# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)

for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.2 %>%
    dplyr::mutate(EMPM = data.table::shift(EMPM, type="lead", i))
  
  reg <- lm(EMPM ~ ., df.reg[, c(1, 3, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}      

regressions_out_empm <- regressions_out
# rename the column names
colnames(regressions_out_empm)[1] <- "EMPM.locproj"
colnames(regressions_out_empm)[2] <- "EMPM.Lo.95"
colnames(regressions_out_empm)[3] <- "EMPM.Up.95"
colnames(regressions_out_empm)[4] <- "EMPM.Lo.68"
colnames(regressions_out_empm)[5] <- "EMPM.Up.68" 

regressions_out_VXO <- data.frame(regressions_out_empm,
                                  regressions_out_ipm)


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
locproj_results_VXO <- create_irf_dataframe_locproj(
  regressions_out_VXO,
  'VXO')


########################################################################
### PART 3: LocProj_Bloom2009_HP with Michigan Survey Index
###         following Leduc and Lui (2016)
########################################################################
# we only select the variables that are actually necessary for this
# step:
locproj8_bloom_JFV.3 <- locproj8_bloom_JFV %>%
  dplyr::select(c(Michigan, IPM, EMPM, HOURSM, CPI, WAGE, FFR,
                  STOCK))

# we want to check whether the procedure below indeed runs all regressions
# correctly; therefore, we create the 12 lags for all 8 variables manually!
locproj8_bloom_JFV.3 <- locproj8_bloom_JFV.3 %>%
  do(setNames(data.frame(., data.table::shift(.$IPM, type="lag", 1:12), 
                         data.table::shift(.$EMPM, type="lag", 1:12),
                         data.table::shift(.$HOURSM, type="lag", 1:12),
                         data.table::shift(.$CPI, type="lag", 1:12),
                         data.table::shift(.$WAGE, type="lag", 1:12),
                         data.table::shift(.$FFR, type="lag", 1:12),
                         data.table::shift(.$STOCK, type="lag", 1:12),
                         data.table::shift(.$Michigan, type="lag", 1:12)), 
              c(names(locproj8_bloom_JFV.3), paste0(rep(c("IPM_lag", 
                                                          "EMPM_lag",
                                                          "HOURSM_lag",
                                                          "CPI_lag",
                                                          "WAGE_lag",
                                                          "FFR_lag",
                                                          "STOCK_lag",
                                                          "Michigan"), 
                                                        each = 12), 1:12) )))


####################
## Local Projections (with 8 variables)
####################
# our goal now is to loop i from 0 to 60 (i.e., the horizons that we 
# want to consider);
# for the regressions we use 'dynlm' to have easy access to LAG-operators
# within the regressions!



# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.3 %>%
    dplyr::mutate(IPM = data.table::shift(IPM, type="lead", i))
  
  reg <- lm(IPM ~ ., df.reg[, c(1:2, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}


regressions_out_ipm <- regressions_out

# rename the column names
colnames(regressions_out_ipm)[1] <- "IPM.locproj"
colnames(regressions_out_ipm)[2] <- "IPM.Lo.95"
colnames(regressions_out_ipm)[3] <- "IPM.Up.95"
colnames(regressions_out_ipm)[4] <- "IPM.Lo.68"
colnames(regressions_out_ipm)[5] <- "IPM.Up.68"  


# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)

for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.3 %>%
    dplyr::mutate(EMPM = data.table::shift(EMPM, type="lead", i))
  
  reg <- lm(EMPM ~ ., df.reg[, c(1, 3, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}      

regressions_out_empm <- regressions_out
# rename the column names
colnames(regressions_out_empm)[1] <- "EMPM.locproj"
colnames(regressions_out_empm)[2] <- "EMPM.Lo.95"
colnames(regressions_out_empm)[3] <- "EMPM.Up.95"
colnames(regressions_out_empm)[4] <- "EMPM.Lo.68"
colnames(regressions_out_empm)[5] <- "EMPM.Up.68" 

regressions_out_Michigan <- data.frame(regressions_out_empm,
                                       regressions_out_ipm)


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
locproj_results_Michigan <- create_irf_dataframe_locproj(
  regressions_out_Michigan,
  'Michigan')



########################################################################
### PART 4: LocProj_Bloom2009_HP with EPU following
###         Baker et al. (2015)
########################################################################
# we only select the variables that are actually necessary for this
# step:
locproj8_bloom_JFV.4 <- locproj8_bloom_JFV %>%
  dplyr::select(c(EPU_Hist, IPM, EMPM, HOURSM, CPI, WAGE, FFR,
                  STOCK))

# we want to check whether the procedure below indeed runs all regressions
# correctly; therefore, we create the 12 lags for all 8 variables manually!
locproj8_bloom_JFV.4 <- locproj8_bloom_JFV.4 %>%
  do(setNames(data.frame(., data.table::shift(.$IPM, type="lag", 1:12), 
                         data.table::shift(.$EMPM, type="lag", 1:12),
                         data.table::shift(.$HOURSM, type="lag", 1:12),
                         data.table::shift(.$CPI, type="lag", 1:12),
                         data.table::shift(.$WAGE, type="lag", 1:12),
                         data.table::shift(.$FFR, type="lag", 1:12),
                         data.table::shift(.$STOCK, type="lag", 1:12),
                         data.table::shift(.$EPU_Hist, type="lag", 1:12)), 
              c(names(locproj8_bloom_JFV.4), paste0(rep(c("IPM_lag", 
                                                          "EMPM_lag",
                                                          "HOURSM_lag",
                                                          "CPI_lag",
                                                          "WAGE_lag",
                                                          "FFR_lag",
                                                          "STOCK_lag",
                                                          "EPU_Hist"), 
                                                        each = 12), 1:12) )))


####################
## Local Projections (with 8 variables)
####################
# our goal now is to loop i from 0 to 60 (i.e., the horizons that we 
# want to consider);
# for the regressions we use 'dynlm' to have easy access to LAG-operators
# within the regressions!



# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.4 %>%
    dplyr::mutate(IPM = data.table::shift(IPM, type="lead", i))
  
  reg <- lm(IPM ~ ., df.reg[, c(1:2, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}


regressions_out_ipm <- regressions_out

# rename the column names
colnames(regressions_out_ipm)[1] <- "IPM.locproj"
colnames(regressions_out_ipm)[2] <- "IPM.Lo.95"
colnames(regressions_out_ipm)[3] <- "IPM.Up.95"
colnames(regressions_out_ipm)[4] <- "IPM.Lo.68"
colnames(regressions_out_ipm)[5] <- "IPM.Up.68"  


# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)

for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.4 %>%
    dplyr::mutate(EMPM = data.table::shift(EMPM, type="lead", i))
  
  reg <- lm(EMPM ~ ., df.reg[, c(1, 3, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}      

regressions_out_empm <- regressions_out
# rename the column names
colnames(regressions_out_empm)[1] <- "EMPM.locproj"
colnames(regressions_out_empm)[2] <- "EMPM.Lo.95"
colnames(regressions_out_empm)[3] <- "EMPM.Up.95"
colnames(regressions_out_empm)[4] <- "EMPM.Lo.68"
colnames(regressions_out_empm)[5] <- "EMPM.Up.68" 

regressions_out_EPU <- data.frame(regressions_out_empm,
                                  regressions_out_ipm)


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
locproj_results_EPU <- create_irf_dataframe_locproj(
  regressions_out_EPU,
  'EPU')

########################################################################
### PART 5: LocProj_Bloom2009_HP with Macro Uncertainty Index
###         following Jurado et al. (2015);
###         We specifically use the data for h=1 and h=12 to see how
###         they perform compared to each other
########################################################################

# we only select the variables that are actually necessary for this
# step:
locproj8_bloom_JFV.5 <- locproj8_bloom_JFV %>%
  dplyr::select(c(Macro1, Macro12, IPM, EMPM, HOURSM, CPI, WAGE, FFR,
                  STOCK))

# we want to check whether the procedure below indeed runs all regressions
# correctly; therefore, we create the 12 lags for all 8 variables manually!
locproj8_bloom_JFV.5 <- locproj8_bloom_JFV.5 %>%
  do(setNames(data.frame(., data.table::shift(.$IPM, type="lag", 1:12), 
                         data.table::shift(.$EMPM, type="lag", 1:12),
                         data.table::shift(.$HOURSM, type="lag", 1:12),
                         data.table::shift(.$CPI, type="lag", 1:12),
                         data.table::shift(.$WAGE, type="lag", 1:12),
                         data.table::shift(.$FFR, type="lag", 1:12),
                         data.table::shift(.$STOCK, type="lag", 1:12),
                         data.table::shift(.$Macro1, type="lag", 1:12),
                         data.table::shift(.$Macro12, type="lag", 1:12)), 
              c(names(locproj8_bloom_JFV.5), paste0(rep(c("IPM_lag", 
                                                          "EMPM_lag",
                                                          "HOURSM_lag",
                                                          "CPI_lag",
                                                          "WAGE_lag",
                                                          "FFR_lag",
                                                          "STOCK_lag",
                                                          "Macro1.",
                                                          "Macro12."), 
                                                        each = 12), 1:12) )))



####################
## Local Projections (with 8 variables)
####################
# our goal now is to loop i from 0 to 60 (i.e., the horizons that we 
# want to consider);
# for the regressions we use 'dynlm' to have easy access to LAG-operators
# within the regressions!

# for h=1:

# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.5 %>%
    dplyr::mutate(IPM = data.table::shift(IPM, type="lead", i))
  
  reg <- lm(IPM ~ ., df.reg[, c(1, 3, 10:105)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}


regressions_out_ipm <- regressions_out

# rename the column names
colnames(regressions_out_ipm)[1] <- "IPM.locproj"
colnames(regressions_out_ipm)[2] <- "IPM.Lo.95"
colnames(regressions_out_ipm)[3] <- "IPM.Up.95"
colnames(regressions_out_ipm)[4] <- "IPM.Lo.68"
colnames(regressions_out_ipm)[5] <- "IPM.Up.68"  


# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.5 %>%
    dplyr::mutate(EMPM = data.table::shift(EMPM, type="lead", i))
  
  reg <- lm(EMPM ~ ., df.reg[, c(1, 4, 9:104)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}      

regressions_out_empm <- regressions_out
# rename the column names
colnames(regressions_out_empm)[1] <- "EMPM.locproj"
colnames(regressions_out_empm)[2] <- "EMPM.Lo.95"
colnames(regressions_out_empm)[3] <- "EMPM.Up.95"
colnames(regressions_out_empm)[4] <- "EMPM.Lo.68"
colnames(regressions_out_empm)[5] <- "EMPM.Up.68" 

regressions_out_macro1 <- data.frame(regressions_out_empm,
                                     regressions_out_ipm)


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
locproj_results_macro1 <- create_irf_dataframe_locproj(
  regressions_out_macro1,
  'Macro1')

# for h=12:

# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.5 %>%
    dplyr::mutate(IPM = data.table::shift(IPM, type="lead", i))
  
  reg <- lm(IPM ~ ., df.reg[, c(2, 3, 10:93, 106:117)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}


regressions_out_ipm <- regressions_out

# rename the column names
colnames(regressions_out_ipm)[1] <- "IPM.locproj"
colnames(regressions_out_ipm)[2] <- "IPM.Lo.95"
colnames(regressions_out_ipm)[3] <- "IPM.Up.95"
colnames(regressions_out_ipm)[4] <- "IPM.Lo.68"
colnames(regressions_out_ipm)[5] <- "IPM.Up.68"  


# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
  a=numeric(),
  b=numeric(),
  c=numeric(),
  d=numeric(),
  e=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 60
  #i <- 1
  print(i)
  
  # for each i, we take the variable "IPM" and shift it accordingly:
  df.reg <- locproj8_bloom_JFV.5 %>%
    dplyr::mutate(EMPM = data.table::shift(EMPM, type="lead", i))
  
  reg <- lm(EMPM ~ ., df.reg[, c(2, 4, 10:93, 106:117)])
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), 
  # regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-
  # neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] + 0.47*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 0.47*SE_robust[2])
  )
  regressions_out <- rbind(regressions_out, estimates)
  
  # and again add a column 'h'
  regressions_out  <- as.data.frame(regressions_out  %>%
                                      mutate(h = seq(0, 
                                                     nrow(regressions_out)-1)))
  
  rm(reg, SE_robust, estimates)
}      

regressions_out_empm <- regressions_out
# rename the column names
colnames(regressions_out_empm)[1] <- "EMPM.locproj"
colnames(regressions_out_empm)[2] <- "EMPM.Lo.95"
colnames(regressions_out_empm)[3] <- "EMPM.Up.95"
colnames(regressions_out_empm)[4] <- "EMPM.Lo.68"
colnames(regressions_out_empm)[5] <- "EMPM.Up.68" 

regressions_out_macro12 <- data.frame(regressions_out_empm,
                                      regressions_out_ipm)


# creation of data-frame for Bloom_shock
# Achtung: wir haben 'rescaled' ersetzt!
locproj_results_macro12 <- create_irf_dataframe_locproj(
  regressions_out_macro12,
  'Macro12')


####################
## plots of loc_projirfs
####################
# before plotting, we have to stack the vairous data-frames that we have
# created, together:
locproj_irfs_results_all <- bind_rows(locproj_results_bloom, 
                                      locproj_results_VXO,
                                      locproj_results_Michigan,
                                      locproj_results_EPU,
                                      locproj_results_macro1,
                                      locproj_results_macro12)


# we need to create ordered factors, otherwise 'facet_wrap' combines the plots
# in alphabetical order:
locproj_irfs_results_all <- locproj_irfs_results_all %>%
  ungroup %>%
  mutate(shock_name = factor(shock_name, 
                             ordered = TRUE,
                             levels=unique(shock_name)),
         response = factor(response,
                           ordered = TRUE,
                           levels=unique(response)))

# now we can finally plot the irfs:

# next we can plot everything we had plotted above
locproj_plot_all_NoHP_until2008 <- ggplot(data = locproj_irfs_results_all) + 
  # geom_point() + 
  # below ribbon is for 68% confidence band
  geom_ribbon(aes(x=step, ymax=100*(exp(Lo.68)-1), 
                  ymin=100*(exp(Up.68)-1),
                  fill=shock_name), 
              alpha=.6) +
  # and this one for the 95% confidence band
  geom_ribbon(aes(x=step, ymax=100*(exp(Lo.95)-1), 
                  ymin=100*(exp(Up.95)-1),
                  fill=shock_name), 
              alpha=.3) +  
  geom_line(aes(x = step, y = 100*(exp(locproj)-1)), 
            color="black", size=0.5) + 
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, limits = c(0, 60), 
                     breaks = seq(0, 60, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = NULL,
                     # limits = c(-4, 4),
                     # breaks = seq(-4, 4, by=2),
                     minor_breaks = NULL) +
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=0.5) +
  # Change line size
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
  coord_fixed(ratio = 5) + 
  facet_grid(shock_name ~ response,
             scales="free_y")

locproj_plot_all_NoHP_until2008

# and we save the plot to be used in our latex-document
ggsave(file="locproj_plot_all_NoHP_until2008.pdf", width = 210, height = 297, units = "mm")


