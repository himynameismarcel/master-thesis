#########################################################################################
### Marcel Kropp, 31.05.2018
### This script runs the production of summary statistics for the various uncertainty
### measures we have introduced in the text.
### Everything that is produced here will end up in the part about 
### 'Stylized Facts' in the main text.

### The goal of the script is the following:
### We want to produce a mix of the TABLE 1 in Jurado et al. (2016)
### and the Table 2 in Bontempi et al. (2016)

### Jurado et al. (2016) report the following summary statistics:
###       * AR(1) - coef, Half Life
###       * Skewness
###       * Kurtosis
###       * IP-corr(0)
###       * IP-corr(12),
###       * IP-corr(-12)
###       * max IP-corr(k) [k>0] at lag k =
###       * max IP-corr(k) [k<0] at lag k =
### Bontempi et al. (2016) add a few more.

### We merge the summary statistics that these two papers report and make
### our own design out of it. 
### In particular, the below routine will produce four data.frames called
###       * Distribution
###       * Cyclicality
###       * Persistence
###       * Stationarity
### and report the respective figures across all uncertainty measures that we
### consider (VXO, Macro1, Macro12, EPU [speaking of EPU we mean
### the Historical EPU!]) and Financial Uncertainty.

### The individual reported data.frames can then be stacked on top of each
### other.

### We want to reproduce this table to give the reader an overview of the grand
### similarities and differences of the various uncertainty measures.



######################################################################################
# 07.05.2020 (Marcel):
# to calculate the summary statistics, we need the uncertainty-series for which
# we want to calculate the summary statistics;

## Although we have constructed the VXO/stock market volatility-series more
## or less ourselves in the script '20200705_sp500_vxo_bloom_shocks_v0.1.R',
## here we go back to tthe actual series that Bloom has provided in his replication
## files including the column VOLATB;

## All other measures can be extracted from the data.frame 'comparison_measures',
## because we took the data directly from the respective sources (and did not construct
## anything ourselves!)

## The data.frame 'comparison_measures_plot' is produced in the R-script
## '120200705_comparison_measures_v0.2.R' but in the version there comes along
## without the separate variables for MONTH and YEAR!
## Therefore we replicate the generation of 'comparison_measures' here, call the
## resulting data-frame 'uncert_measures.stats' (to avoid confusion) and add in the
## additional variables that we want/need.


## As mentioned above, the dataset that we read in in a preliminary version of this 
## script, is the dataset including the VOLATBL-series provided by Bloom himself, which
## we have stored in 'VARDATA_UPDATED_JFV.xlsx';
## We will replace this step with our own data as soon as we have clarified whether
## our own generated series and the series provided by Bloom is indeed the same!
## Update Marcel (05.07.2020): Ultimately, I have decided to simply stick to the original 
## replication data of Bloom.


###############################
## loading libraries
###############################
if(!exists("foo", mode="function")) source("20200705_load_libraries_v0.1.R")


###############################
## reading in NBER recession dates
###############################
if(!exists("foo", mode="function")) source("20200705_nber_recession_dates_v0.1.R")


###############################
## reading in all other scripts
## that produce the respective series
###############################
# if(!exists("foo", mode="function")) source("20200705_sp500_vxo_bloom_shocks_v0.1.R")
# we do not run the script withthe preparation of the Bloom-series, because
# we take it from the replication data that is readily available
# online!
if(!exists("foo", mode="function")) source("20200705_macro_fin_uncert_v0.1.R")
if(!exists("foo", mode="function")) source("20200705_EPU_baker_et_al_v0.1.R")
if(!exists("foo", mode="function")) source("20200705_additional_vars_v0.1.R")

###############################
## reading in the script that holds all functions
###############################
if(!exists("foo", mode="function")) source("20200628_functions_v0.2.R")
# Note: at the moment this script cannot be run because we still have 
# to fix the functions for the IRFs!

#####################
## preliminaries 1
#####################
# reading in or using the 'VARDATA_BLoom_original.csv' is actually
# not necessary any more because we do not use it
# (we have adequately prepared the file 'sp500_merge_vxo' in another
# script!)
VARDATA_Bloom_original <- read.csv(file="VARDATA_Bloom_original.csv", 
                                   header=TRUE, sep=",")

# we create a variable called 'my'
VARDATA_Bloom_original <- VARDATA_Bloom_original %>%
  mutate(my = YEAR + MONTH/12)

# Now we can start with the construction of the data.frame 'identify_shocks'!
# Note that we here explicitly also include the h12 Macro uncertainty index
# (although the results should be very similar for both h1 and h12!)
uncert_measures.stats <- data.frame("macroUncert_h1" = macroUncertainty_index$h1,
                              "macroUncert_h12" = macroUncertainty_index$h12,
                              "my" = macroUncertainty_index$my,
                              "year" = macroUncertainty_index$year,
                              "month" = macroUncertainty_index$month,
                              "Date" = macroUncertainty_index$Date)

# Then we join in the other uncertainty measures one after the other: 
# (note that we have skipped the addition of the GTU-index!)
#   * followed by the original epu_index,
#   * as well as the historical EPU_index,
#   * then the VXO (called VOLATBL),
#   * in a multi-setp join:
uncert_measures.stats <- as.tibble(right_join(x = epu_index[, 
                                        c("News_Based_Policy_Uncert_Index","my")],
                                        y = uncert_measures.stats,
                                        by = "my") %>%
                               right_join(x = epu_index_historical[, 
                                        c("EPU_Historical","my")],
                                        y = .,
                                        by = "my") %>%
                               # Note:
                               # after an update we can read in 'mvol_VOLATBL'
                               # from the data.frame 'sp500_merge_vxo' and do not
                               # have to rely on the shorter VOLATBL-series
                               # from Bloom's original data!
                               right_join(x = sp500_merge_vxo[, 
                                          c("mvol_VOLATBL","my")],
                                          y = .,
                                          by = "my") %>%
                               # further, we rename all the variables,
                               # select the ones we want to keep
                               dplyr::rename(VXO = mvol_VOLATBL,
                                             EPU = News_Based_Policy_Uncert_Index,
                                             EPU_Hist = EPU_Historical,
                                             #GTU = GTU_US,
                                             Macro1 = macroUncert_h1,
                                             Macro12 = macroUncert_h12) %>%
                               # we further reorder the variables and drop 'my' by 
                               # leaving it out!
                               # update: we need the 'my' variable again further down
                               # below: therefore we should not drop it!
                               dplyr::select(Date, year, month, VXO, EPU, EPU_Hist,
                                             Macro1, Macro12))



#####################
## preliminaries 2
#####################
# Marcel (27.05.2020):
# we add the financial uncertainty measures of Ludvigson et al
# (2018) to it:
uncert_measures.stats <- right_join(x = finUncertainty_index[, c("Uf", "year", "month")], 
                                  y = uncert_measures.stats, 
                                  by = c("year", "month"))
# we put the column 'Uf' last:
uncert_measures.stats <- uncert_measures.stats %>% 
                                  dplyr::select(-Uf,Uf)


# According to the four categories that we want to consider (each producing
# a dedicated data.frame, i.e., for Distribution, Cyclicality, Persistence,
# Stationarity), we split the below function into four parts:


                      #---------------------------------------------------------
                      # (1) DISTRIBUTION:
                      # We want to calculate a few summary statistics for 
                      # each uncertainty measure;
                      # In particular for VXO, EPU_Hist, Macro, FinUnc
                      
                      # To calculate our summary statistics, we want to data.frame 
                      # in a tidy format (i.e., long instead of wide):
                      # For this we select only the variable we want to keep 
                      # (meaning we drop EPU),
                      # and gather the data to get it into a tidy format:
                      
                      uncert_measures.stats.tidy <- uncert_measures.stats %>%
                                              dplyr::select(-EPU) %>%
                                              # we declare 
                                              gather( 
                                                uncert_measure, value, 
                                                -c(Date, year, month), 
                                                na.rm=TRUE) %>%
                                              dplyr::mutate(uncert_measure = 
                                                            as.factor(
                                                            uncert_measure)) %>%
                                              # at the same time we gorup by
                                              # uncert_measure
                                              group_by(uncert_measure) %>%
                                              # and add a column that holds
                                              # the changes of the respective
                                              # uncertainty measures:
                                              dplyr::mutate(change = 
                                                              value - lag(value))
                      
                      
                      # calculate mean and standard deviation across measures
                      # note that we transpose the entire result to have
                      # in a desired format already:
                      summary_stats <- as.data.frame(t(uncert_measures.stats.tidy %>%
                        group_by(uncert_measure) %>%
                        summarise(mean = mean(value), sd = sd(value),
                                  kurtosis = kurtosis(value)-3,
                                  skewness = skewness(value),
                                  kurtosis_change = kurtosis(change, na.rm=TRUE)-3,
                                  skewness_change = skewness(change, na.rm=TRUE),
                                  shapiro = shapiro.test(value)$p.value) %>%
                        dplyr::mutate(coeff_of_variation = sd/mean) %>% 
                        dplyr::select(uncert_measure, 
                                      coeff_of_variation, skewness, 
                                      kurtosis, skewness_change,
                                      kurtosis_change, shapiro)))
                      # we assign new col-names
                      names(summary_stats) <- summary_stats[1, ]
                      # and remove the first row
                      summary_stats <- summary_stats[-1, ]
                      # we transform all data-points to numeric:
                      summary_stats[] <- lapply(summary_stats, 
                                                      function(x) 
                                                        as.numeric(as.character(x)))
                      # we reduce the number of decimals
                      is.num <- sapply(summary_stats, is.numeric)
                      summary_stats[is.num] <- lapply(summary_stats[is.num], round, 2)
                      
                      # and reorder the columns
                      summary_stats <- summary_stats %>%
                                  dplyr::select(VXO, EPU_Hist, 
                                                Macro1, Macro12, Uf) %>%
                                  # at the same time we rename EPU_Hist to EPU
                                  dplyr::rename(EPU = EPU_Hist)

                      
                      #---------------------------------------------------------
                      # (2) CYCLICALITY:
                      # We want to calculate a few cyclicality measures for 
                      # each uncertainty measure;
                      # In particular for VXO, EPU_Hist, Michigan, Macro1, 
                      # Macro12;
                      
                      # The measures we want to calculate are:
                      #     * Downturn/upturn mu-ratios
                      #     * Downturn/upturn sigma-ratios
                      #     * IP-corr(0)
                      #     * IP-corr(12)
                      #     * IP-corr(-12)
                      #     * max IP-corr(k) at lag
                      #     * max IP-corr(k) at lag
                      
                      # for the downturn/upturn ratios we need to read in
                      # the NBER recession dates:
                      # from the script '20200705_nber_recession_dates_v0.1.R',
                      # which we have executed at the start of this script,
                      # the recessions are available in the data.frame
                      # 'nber_recessions' including a dummy 'USREC' 
                      # equal to 1 in case of a recession month;
                      # in addition the data.frame already holds the 
                      # variables 'my', 'year', 'month', etc. 
                      # View(nber_recessions)

                      
                      
                      # we first merge the data.frame 'nber_recessions' with
                      # our data.frame 'identify_shocks'
                      
                      uncert_measures.stats.rec <- inner_join(
                        x = uncert_measures.stats,
                        y = nber_recessions,
                        by = c("year", "month", "Date"))
                      
                      # we are now in a position to calculate the mean 
                      # and standard deviation again across uncertainty
                      # measures but this time also group by
                      # USREC == 1!
                      
                      # Marcel (05.07.2020):
                      # When re-running this script, we ran into the problem of having
                      # two my-variables (my.x and my.y); hence I have reduced it to
                      # one:
                      uncert_measures.stats.rec <- uncert_measures.stats.rec %>%
                        dplyr::select(-my.x) %>%
                        rename(my = my.y)
                      
                      
                      # we bring our data.frame in a tidy format again;
                      # this time including the variable USREC
                      uncert_measures.stats.rec.tidy <- 
                        uncert_measures.stats.rec %>%
                        dplyr::select(-EPU) %>%
                        gather( 
                          uncert_measure, value, 
                          -c(Date, year, month, USREC, my), 
                          na.rm=TRUE) %>%
                        dplyr::mutate(uncert_measure = 
                                        as.factor(
                                          uncert_measure))
                      
                      # next we calculate mean and standard deviation again
                      # but this time group by both uncert_measure and
                      # USREC
                      summary_stats2 <- as.data.frame(t(
                                        uncert_measures.stats.rec.tidy %>%
                                        group_by(uncert_measure, USREC) %>%
                                        summarise(mean = mean(value), 
                                                  sd = sd(value)) %>%
                        # we process the result from above further
                        # to calculate the respective ratios
                        # (downturn/upturn ratios):
                                        group_by(uncert_measure) %>%
                                        summarise(downturn_upturn_ratio_mean=
                                                    mean[USREC==1]/
                                                    mean[USREC==0],
                                                  downturn_upturn_ratio_sd=
                                                    sd[USREC==1]/
                                                    sd[USREC==0])))

                      # we assign new col-names
                      names(summary_stats2) <- summary_stats2[1, ]
                      # and remove the first row
                      summary_stats2 <- summary_stats2[-1, ]
                      # we transform all data-points to numeric:
                      summary_stats2[] <- lapply(summary_stats2, 
                                                function(x) 
                                                  as.numeric(as.character(x)))
                      # we reduce the number of decimals
                      is.num <- sapply(summary_stats2, is.numeric)
                      summary_stats2[is.num] <- lapply(summary_stats2[is.num], 
                                                       round, 2)
                      
                      # and reorder the columns
                      summary_stats2 <- summary_stats2 %>%
                        dplyr::select(VXO, EPU_Hist, 
                                      Macro1, Macro12, Uf) %>%
                        # at the same time we rename EPU_Hist to EPU
                        dplyr::rename(EPU = EPU_Hist)
                      
                      
                      # next, we want to calculate a few measures related
                      # to the correlation of our uncertainty-measures
                      # with industrial production in the spirit of
                      # Jurado et al. (2016);
                      # To be able to calculate these correlation measures, we 
                      # need to add a time-series of industrial production 
                      # into our analysis;
                      # at this point (we will also change that in the other
                      # analyses), we will refer back to 
                      # industrial production (all) instead of 
                      # industrial production in manufacturing, because
                      # industrial production in manufacturing is not 
                      # available as of 1960!
                      
                      # the script '20200507_additional_vars'
                      # treats the preparation of the
                      # data.frame 'additional_vars' which, among others
                      # also contains the variable 'industrial production'
                      # which we need for the below:
                      # (originall additional_vars was dealt with in the
                      # script 21042018_thesis_master.R)
                      
                      # We want to calculate the following:
                      #     * IP-corr(0)
                      #     * IP-corr(12)
                      #     * IP-corr(-12)
                      #     * max IP-corr(k) at lag
                      #     * max IP-corr(k) at lag 
                      #     * EMP-corr(0)
                      #     * EMP-corr(12)
                      #     * EMP-corr(-12)
                      #     * max EMP-corr(k) at lag
                      #     * max EMP-corr(k) at lag
                      #     * for later: we could potentially
                      #       here also add the results for
                      #       employment in manufacturing
                      
                      # According to the information in the TABLE 1
                      # of Jurado et al. (2016), 
                      # IP-corr(k) is the absolute cross-correlation between
                      # a measure of uncertainty ut and 12 month moving
                      # average of industrial production growth in
                      # period t+k, i.e., IP-corr(k) = |corr(ut, delta ln IP t-k)|;
                      # A positive k means uncertainty is correlated with future IP!
                      
                      # First, out of the series for industrial production
                      # and employment in manufacturing
                      # in the data.frame 'additional_vars', we create a
                      # 12 month moving average of industrial production growth
                      # and employment in manufacturing:
                      ip_empm_12mon_MA.growth <- additional_vars %>%
                                # we make a time-series object out of ip
                                dplyr::mutate(ip = zoo(ip, date),
                                              empm = zoo(empm, date)) %>%
                                # and then calcualte the percentage
                                # growth from month to month
                                dplyr::mutate(ip = perc2(ip),
                                              empm = perc2(empm)) %>%
                                # in a last stage we apply a roll-apply
                                dplyr::mutate(ip12MAperc = rollapply(ip, 12, 
                                                             mean,
                                                             fill = NA,
                                                             na.rm = T),
                                              empm12MAperc = rollapply(empm, 12,
                                                             mean,
                                                             fill = NA,
                                                             na.rm  = T)) %>%
                      # see: https://danieljhocking.wordpress.com/2014/12/03/lags-
                      # and-moving-means-in-dplyr/                   
                                # in a last step we restrict the number of columns
                                dplyr::select(year, month, my, ip12MAperc, empm12MAperc)
                      
                      
                      # having the 12 month moving average of ip growth 
                      # and employment growth available,
                      # we can now calculate IP-corr(0), EMP-corr(0) and all other
                      # measures for all uncertainty measures:
                      # to do that, we first merge ip_empm_12mon_MA.growth
                      # with uncert_measures.stats
                      uncert_measures.stats.ip_emp12monMA <- (right_join(
                        x = uncert_measures.stats,
                        y = ip_empm_12mon_MA.growth,
                        by = c("year", "month"))) %>%
                        # and we remove the column 'Date'
                        dplyr::select(-Date)
                      
                      
                      # Marcel (05.07.2020):
                      # When re-running this script, we ran into the problem of having
                      # two my-variables (my.x and my.y); hence I have reduced it to
                      # one:
                      uncert_measures.stats.ip_emp12monMA <- uncert_measures.stats.ip_emp12monMA %>%
                        dplyr::select(-my.x) %>%
                        rename(my = my.y)
                      
                      
                      # we bring our data.frame in a tidy format again;
                      # this time including the variable ip and empm
                      # with their growth rates
                      uncert_measures.stats.ip_emp12monMA.tidy <- 
                        uncert_measures.stats.ip_emp12monMA %>%
                        dplyr::select(-EPU) %>%
                        gather( 
                          uncert_measure, value, 
                          -c(year, ip12MAperc, month, my, empm12MAperc), 
                          na.rm=TRUE) %>%
                        dplyr::mutate(uncert_measure = 
                                        as.factor(
                                          uncert_measure))  %>%
                        # we also drop rows with NAs
                        drop_na
                      
                      # now we are in a position to calculate the cross-
                      # correlations;
                      # for these we make use of the Lag-operator
                      # from the Misc-package because, interestingly,
                      # the others (from dplyr, etc.) did not produce
                      # the desired result for shifts backwards (i.e., leads)
                      summary_stats3 <- as.data.frame(t(
                        uncert_measures.stats.ip_emp12monMA.tidy %>%
                        group_by(uncert_measure) %>%
                        dplyr::summarise(IPcorr0 = cor(ip12MAperc, value,
                                             use = "pairwise.complete.obs"),
                               IPcorr_12 = cor(Hmisc::Lag(ip12MAperc, 12),
                                               value, 
                                               use = "pairwise.complete.obs"),
                               IPcorr12 = cor(Hmisc::Lag(ip12MAperc, -12),
                                              value,
                                              use = "pairwise.complete.obs"),
                               EMPcorr0 = cor(empm12MAperc, value,
                                             use = "pairwise.complete.obs"),
                               EMPcorr_12 = cor(Hmisc::Lag(empm12MAperc, 12),
                                               value, 
                                               use = "pairwise.complete.obs"),
                               EMPcorr12 = cor(Hmisc::Lag(empm12MAperc, -12),
                                              value,
                                              use = "pairwise.complete.obs"))))
                      
                      # we assign new col-names
                      names(summary_stats3) <- summary_stats3[1, ]
                      # and remove the first row
                      summary_stats3 <- summary_stats3[-1, ]
                      # we transform all data-points to numeric:
                      summary_stats3[] <- lapply(summary_stats3, 
                                                 function(x) 
                                                   as.numeric(as.character(x)))
                      # we reduce the number of decimals
                      is.num <- sapply(summary_stats3, is.numeric)
                      summary_stats3[is.num] <- lapply(summary_stats3[is.num], round, 2)
                      
                      # and reorder the columns
                      summary_stats3 <- summary_stats3 %>%
                        dplyr::select(VXO, EPU_Hist, 
                                      Macro1, Macro12, Uf) %>%
                        # at the same time we rename EPU_Hist to EPU
                        dplyr::rename(EPU = EPU_Hist)


                      #---------------------------------------------------------
                      # (3) PERSISTENCE:
                      # We want to calculate a few persistence measures for 
                      # each uncertainty measure;
                      # In particular for VXO, EPU_Hist, Michigan, Macro1, 
                      # Macro12;
                      
                      # The measures we want to calculate are:
                      #     * AR(1) - coefficients
                      #     * based on this result the Half Life
                      # Because we do not need any additional variables
                      # this time, we can go ahead and make use of our
                      # data.frame called: 'uncert_measures.stats'
                      uncert_measures.stats2 <- uncert_measures.stats %>%
                        dplyr::mutate(VXO = ts(VXO, start = c(1960, 7),
                                               frequency = 12),
                                      EPU_Hist = ts(EPU_Hist, start = c(1960, 7),
                                                    frequency = 12),
                                      Macro1 = ts(Macro1, start = c(1960, 7),
                                                   frequency = 12),
                                      Macro12 = ts(Macro12, start = c(1960, 7),
                                                    frequency = 12),
                                      Uf = ts(Macro12, start = c(1960, 7),
                                                   frequency = 12)) %>%
                        # drop columns 'year' and 'month'
                        dplyr::select(-c(year, month, EPU))
                      
                      # not we can loop through all variables and fit
                      # ar- models:
                      # we extract the names of the columns:
                      nam <- names(uncert_measures.stats2)[3:length(uncert_measures.stats2)]
                      # and then loop through the column names:
                      
                      # we initialize an empty vector:
                      my.ar.coef <-  numeric(0)
                      
                      for(i in nam){
                        my.ar <- ar(uncert_measures.stats2[[i]], order.max = 1,
                                    method = "mle", na.action = na.omit)
                        # we extract the coefficient:
                        my.ar.coef <- c(my.ar.coef, my.ar$ar)
                      }
                      
                      # next we store the AR coefficients in a data.frame
                      # and transpose it on the way
                      my.ar.coef <- as.data.frame(t(my.ar.coef),
                                                  col.names = c(nam))
                      
                      # because col.names does not work, we manually change the names:
                      my.ar.coef <- my.ar.coef %>%
                            dplyr::rename(VXO = V1,
                                          EPU = V2, 
                                          Macro1 = V3,
                                          Macro12 = V4,
                                          Uf = V5)
                      
                      # next, we calcualte the half-life of the series
                      # off of the AR-coefficients:
                      # first we transform the columns to a row
                      my.ar.coef.tidy <- my.ar.coef %>%
                                # we declare 
                                gather( 
                                  uncert_measure, ar_coef) 
                      
                      # and then we calculate the half-life:
                      summary_stats4 <- as.data.frame(t(my.ar.coef.tidy %>%
                            dplyr::mutate(half_life = -log(2)/log(ar_coef))))
                      
                      # we assign new col-names
                      names(summary_stats4) <- summary_stats4[1, ]
                      # and remove the first row
                      summary_stats4 <- summary_stats4[-1, ]
                      # we transform all data-points to numeric:
                      summary_stats4[] <- lapply(summary_stats4, 
                                                 function(x) 
                                                   as.numeric(as.character(x)))
                      # we reduce the number of decimals
                      is.num <- sapply(summary_stats4, is.numeric)
                      summary_stats4[is.num] <- lapply(summary_stats4[is.num], round, 2)
                      
                      # and reorder the columns
                      summary_stats4 <- summary_stats4 %>%
                        dplyr::select(VXO, EPU, 
                                      Macro1, Macro12, Uf)
                      
                      
                      #---------------------------------------------------------
                      # (4) Stationarity:
                      # We want to test for the stationarity of our
                      # uncertainty measures!
                      # In particular for VXO, EPU_Hist, Michigan, Macro1, 
                      # Macro12;
                      
                      # we remove Date, year and month from the data.frame:
                      uncert_measures.stats3 <- uncert_measures.stats %>%
                        dplyr::select(-c(Date, year, month, EPU, my)) %>%
                        dplyr::mutate(VXO = ts(VXO, start = c(1960, 7),
                                               frequency = 12),
                                      EPU_Hist = ts(EPU_Hist, start = c(1960, 7),
                                                    frequency = 12),
                                      Macro1 = ts(Macro1, start = c(1960, 7),
                                                  frequency = 12),
                                      Macro12 = ts(Macro12, start = c(1960, 7),
                                                   frequency = 12),
                                      Uf = ts(Macro12, start = c(1960, 7),
                                                   frequency = 12)) %>%
                        # we drop EPU
                        dplyr::rename(EPU = EPU_Hist)
                      # custom function sof rtests
                      kpss.result.level <- matrix(NA, nrow = 1, ncol = 
                                                    ncol(uncert_measures.stats3))
                      colnames(kpss.result.level) <- colnames(uncert_measures.stats3)
                      for (i in 1 : ncol(uncert_measures.stats3)){
                        kpss.result.level[,i ] <- kpss.test(ts(uncert_measures.stats3[,i]))$p.value
                        print(i)
                      }
                      
                      kpss.result.trend <- matrix(NA, nrow = 1, ncol = 
                                                    ncol(uncert_measures.stats3))
                      colnames(kpss.result.trend) <- colnames(uncert_measures.stats3)
                      for (i in 1 : ncol(uncert_measures.stats3)){
                        kpss.result.trend[,i] <- kpss.test(ts(uncert_measures.stats3[,i]), 
                                                           null = "Trend")$p.value
                        print(i)
                      } 
                      
                      summary_stats5 <- as.data.frame(rbind(kpss.result.level,
                                              kpss.result.trend))
                      # we reduce the number of decimals
                      is.num <- sapply(summary_stats5, is.numeric)
                      summary_stats5[is.num] <- lapply(summary_stats5[is.num], round, 2)
                      
                      # we change the row-names
                      row.names(summary_stats5) <- c("kpss.pvalue.level",
                                                     "kpss.pvalue.trend")
                      
                      
                      #---------------------------------------------------------
                      # (5) Combine everything
                      # At the very end, we can combine all
                      # calculated results:
                      # because rbind does not preserve the
                      # row-names, we apply the below custom function:
                      # found here: https://stackoverflow.com/questions/
                      # 14799434/rbind-two-data-frame-preserving-row-order-and-row-names
                      RBIND <- function(datalist) {
                        require(plyr)
                        temp <- rbind.fill(datalist)
                        rownames(temp) <- unlist(lapply(datalist, row.names))
                        temp
                      }
                      summary_stats.all <- 
                        RBIND(list(summary_stats, summary_stats2, 
                              summary_stats3, summary_stats4,
                              summary_stats5))
                      
                      # and export to latex:
                      latextable(summary_stats.all, 
                                 colnames = c("VXO", "Michigan", "EPU", "Macro1",
                                              "Macro12", "Uf"),
                                 rownames = c("Coeff. of Variation",
                                              "Skewness (Levels)",
                                              "Excess Kurtosis (Levels)",
                                              "Skewness (Change)",
                                              "Excess Kurtosis (Change)",
                                              "Shapiro-Worlk (p-value)",
                                              "Downturn/Upturn ratios",
                                              "Downturn/Upturn ratios",
                                              "IP-corr(0)",
                                              "IP-corr(-12)",
                                              "IP-corr(12)",
                                              "EMP-corr(0)",
                                              "EMP-corr(-12)",
                                              "EMP-corr(12)",
                                              "AR(1)",
                                              "Half Life (in months)",
                                              "Level Stationary (p-value)",
                                              "Trend Stationary (p-value)"))


                      
########################################################################
### PART x: correlation matrix between the uncertainty measures
########################################################################
comparison_measures_corr <- uncert_measures.stats[c(5, 7, 8)]

## correlation matrix between the time-series
corr <- cor(comparison_measures_corr, use = "complete.obs")
corr <- round(corr, 2)

stargazer(corr)

                      
                                            
