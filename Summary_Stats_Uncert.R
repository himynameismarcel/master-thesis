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
### consider (VXO, Mchigan, Macro1, Macro12, EPU [speaking of EPU we mean
### the Historical EPU!]).

### The individual reported data.frames can then be stacked on top of each
### other.

### We want to reproduce this table to give the reader an overview of the grand
### similarities and differences of the various uncertainty measures.



### Note: below is a customized function to calculate percentage
### changes of time series objects that adds an NA on top so that
### the resulting vector can be added back in into the respective
### structure (e.g., a data-frame)

### we have called the function 'perc2' because we already defined
### a perc-function somwhere else!
perc2 <- function(x){
  x <- round(diff(x) / x *100, 2)
  x <- c(NA, x)
  return(x)
}


#####################
## preliminaries
#####################

# Note that we have already stored all uncertainty measures 
# (with their complete length) in the script 'Identification_of_Shocks.R'
identify_shocks

# to prevent any unexpected changes to that data.frame (and leave it as is),
# we assign it to 
uncert_measures.stats <- identify_shocks

# According to the four categories that we want to consider (each producing
# a dedicated data.frame, i.e., for Distribution, Cyclicality, Persistence,
# Stationarity), we split the below function into four parts:


                      #---------------------------------------------------------
                      # (1) DISTRIBUTION:
                      # We want to calculate a few summary statistics for 
                      # each uncertainty measure;
                      # In particular for VXO, EPU_Hist, Michigan, Macro1, 
                      # Macro12
                      
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
                                              dplyr::mutate(change = value - lag(value))
                      
                      
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
                                  dplyr::select(VXO, Michigan, EPU_Hist, 
                                                Macro1, Macro12) %>%
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
                      # from the script '12052018_thesis_master.R',
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
                      summary_stats2[is.num] <- lapply(summary_stats2[is.num], round, 2)
                      
                      # and reorder the columns
                      summary_stats2 <- summary_stats2 %>%
                        dplyr::select(VXO, Michigan, EPU_Hist, 
                                      Macro1, Macro12) %>%
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
                      
                      # the script '1205208_thesis_master.R' also
                      # also treates the preparation of the
                      # data.frame 'additional_vars' which, among others
                      # also contains the variable 'industrial production'
                      # which we need for the below:
                      
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
                      # 12 month moving aerage of industrial production growth
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
                        dplyr::select(VXO, Michigan, EPU_Hist, 
                                      Macro1, Macro12) %>%
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
                                      Michigan = ts(Michigan, start = c(1960, 7),
                                                    frequency = 12),
                                      Macro1 = ts(Macro1, start = c(1960, 7),
                                                   frequency = 12),
                                      Macro12 = ts(Macro12, start = c(1960, 7),
                                                    frequency = 12)) %>%
                        # drop columns 'year' and 'month'
                        dplyr::select(-c(year, month, EPU))
                      
                      # not we can loop through all variables and fit
                      # ar- models:
                      # we extract the names of the columns:
                      nam <- names(uncert_measures.stats2)[2:length(uncert_measures.stats2)]
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
                                          Michigan = V3,
                                          Macro1 = V4,
                                          Macro12 = V5)
                      
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
                        dplyr::select(VXO, Michigan, EPU, 
                                      Macro1, Macro12)
                      
                      
                      #---------------------------------------------------------
                      # (4) Stationarity:
                      # We want to test for the stationarity of our
                      # uncertainty measures!
                      # In particular for VXO, EPU_Hist, Michigan, Macro1, 
                      # Macro12;
                      
                      # we remove Date, year and month from the data.frame:
                      uncert_measures.stats3 <- uncert_measures.stats %>%
                        dplyr::select(-c(Date, year, month, EPU)) %>%
                        dplyr::mutate(VXO = ts(VXO, start = c(1960, 7),
                                               frequency = 12),
                                      EPU_Hist = ts(EPU_Hist, start = c(1960, 7),
                                                    frequency = 12),
                                      Michigan = ts(Michigan, start = c(1960, 7),
                                                    frequency = 12),
                                      Macro1 = ts(Macro1, start = c(1960, 7),
                                                  frequency = 12),
                                      Macro12 = ts(Macro12, start = c(1960, 7),
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
                                              "Macro12"),
                                 rownames = c("Coeff. of Variation",
                                              "Skewness",
                                              "Kurtosis",
                                              "Downturn/Upturn ratios",
                                              "Downturn/Upturn ratios",
                                              "IP-corr(0)",
                                              "IP-corr(-12)",
                                              "IP-corr(12)",
                                              "AR(1)",
                                              "Half Life",
                                              "Level Stationary",
                                              "Trend Stationary"))
  
                      
#--------------------------------------------------------------------------------
# based on the data.frame 'uncert_measures.stats.rec' which, besides
# all uncertainty measures also contains a recession dummy, we want
# to reproduce what Kose & Terrores plotted in their 2012 paper
# 'How Does Uncertainty Affect Economic Performance?'

# The authors plot the evolution of uncertainty in the U.S. starting
# three quarters after uncertainty reached its peak during a
# a previous recession. 
# Because we also consider data from the Michigan Survey, we can
# only look at 5 recessions: 1981, 1983, 1991, 2991 and 2007;

# The procedure works as follows:
# (1) we have to standardize all uncertainty-measures in
# uncert_measures.stats.rec;
# (2) we average across all uncertainty-measures to produce
# one uncertainty-series;
# (3) We search for the date with the highest value during a recession
# (4) we extract the trajectories pre- and post those identified
#     peak-points for all recessions!

# preliminaries:
uncert_measures.stats.rec2 <- uncert_measures.stats.rec

                                            
# (1) we standardize all uncertainty measures:
uncert_measures.stats.rec2 <- uncert_measures.stats.rec2 %>%
            dplyr::mutate(VXO = (VXO - mean(VXO, na.rm=TRUE))/sd(VXO, na.rm=TRUE),
                          EPU_Hist  = (EPU_Hist - mean(EPU_Hist, na.rm=TRUE))/
                            sd(EPU_Hist, na.rm=TRUE),
                          Michigan  = (Michigan - mean(Michigan, na.rm=TRUE))/
                            sd(Michigan, na.rm=TRUE),
                          Macro1  = (Macro1 - mean(Macro1, na.rm=TRUE))/
                            sd(Macro1, na.rm=TRUE)) %>%
            # and we remove the original EPU to only keep the historical
            # we also drop Macro12
            dplyr::select(-c(EPU, Macro12)) %>%
            # at the same time we limit the data.frame to the period where
            # we have data for all series:
            drop_na()

# (2) we average across all uncertainty-measures to produce
# one uncertainty-series:
# note: we might get rid of this step and look at all
# measures at once (still to be decided!)
# only: note: we need ONE uncert_measure to decide for the
# peak of uncertainty during a recessionary period;
# therefore we take it with us for now!
uncert_measures.stats.rec2 <- uncert_measures.stats.rec2 %>%
            dplyr::mutate(uncert_mean = 
                            rowMeans(
                              uncert_measures.stats.rec2[, 4:7]))

# (3) We search for the date with the highest value during a recession
# For this, we use the same technique which we have already
# applied above:

                # let us first filter for all rows where the 
                # USREC = 1 to get an 
                # overview:
                as.data.frame(uncert_measures.stats.rec2 %>% 
                                filter(USREC == 1))
                
                # we initialize the grouping variable x as follows:
                x <- 2
                
                # we add an empty column to our data frame that 
                # will hold the respective ID of a recessions:
                uncert_measures.stats.rec2["REC_ID"] <- NA
                
                # then we start the loop:
                for (i in 1:nrow(uncert_measures.stats.rec2)) {
                  if(uncert_measures.stats.rec2$USREC[i] == 1) {
                    # if we detect a shock (i.e., 'bloom_shock' == 1), 
                    # then we have 
                    # to proceed as follows:
                    # we store the 'my' variable in the column 'start' 
                    # and 'end'
                    uncert_measures.stats.rec2[i, 
                              ncol(uncert_measures.stats.rec2)] <- x
                    
                  }else{
                    # no shock
                    x <- x + 1
                  }
                }

                # we can now have a look at the newly created 
                # variable (note that we have suppressed 'NA's for 
                # the creation of the table!):
                
                uncert_measures.stats.rec2 %>%
                                      filter(!is.na(REC_ID))  %>%
                                      group_by(REC_ID) %>%
                                      summarise(Count = n())
                
                # next we group_by REC_ID and extract the respective 
                # start-dates and store the retrieved data to a new 
                # data-frame
                rec_start <- as.data.frame(uncert_measures.stats.rec2 %>%
                                      filter(!is.na(REC_ID))  %>%
                                      group_by(REC_ID) %>%
                                      filter(row_number()==1) %>%
                                      mutate(year_start=year, 
                                             month_start=month, 
                                             my_start = my) %>%
                                      dplyr::select(REC_ID, 
                                                    year_start, 
                                                    month_start, 
                                                    my_start))
                
                # next, we replicate the above query for the end-dates
                # (note that in some scenarios start- and end-dates are 
                # identical!)
                rec_end <- as.data.frame(uncert_measures.stats.rec2 %>%
                                      filter(!is.na(REC_ID))  %>%
                                      group_by(REC_ID) %>%
                                      filter(row_number()==n()) %>%
                                      mutate(year_end=year, 
                                             month_end=month, 
                                             my_end = my) %>%
                                      dplyr::select(REC_ID, year_end, 
                                                    month_end, my_end))
                
                # next, we can merge the two data-frames from above:
                rec_start_end <- merge(x = rec_start, 
                                      y = rec_end, 
                                      by = "REC_ID", all = TRUE, 
                                      na.rm=T)
                # and inspect the resulting data-frame:
                rec_start_end
                
                # we re-arrange the sequence of columns
                rec_start_end <- rec_start_end %>%
                            dplyr::select(REC_ID, year_start, 
                                      month_start, year_end, 
                                      month_end, my_start, 
                                      my_end)
                
                
                # next, we construct a year-month-variable both out of
                # the pair 'year_start' & 'month_start' and 'year_end' 
                # & 'month_end'
                rec_start_end$yearmon_start <- as.yearmon(
                                      rec_start_end$my_start, 
                                      "%Y-%B")
                rec_start_end$yearmon_end <- as.yearmon(
                                      rec_start_end$my_end, 
                                      "%Y-%B")
                # next, we add a helper-column that we need in the 
                # next stage as well:
                rec_start_end$helper_date <- format(
                                      rec_start_end$yearmon_start, 
                                      "%b")
                
                # this means that we can now drop the 'year' and 
                # 'month' variables
                rec_start_end <- rec_start_end %>%
                                      dplyr::select(-c(year_start, 
                                                       year_end, 
                                      month_start, month_end))
                
                # next we add a column that gives us the duration of 
                # the shock:
                rec_start_end <- rec_start_end %>%
                                      mutate(duration = (
                                        rec_start_end$yearmon_end - 
                                          rec_start_end$yearmon_start) 
                                        * 12 + 1) %>%
                                      # and we drop a few helper-
                                      # variables which
                                      # we won't need anymore
                                      dplyr::select(-c(helper_date, 
                                                       yearmon_start,
                                                    yearmon_end))
                
                # next we group_by REC_ID and extract the respective 
                # MAXIMUM VOLATILITY
                # (note that the produced data-frame only has one row per 
                # shock_ID)
                rec_max <- as.data.frame(uncert_measures.stats.rec2 %>%
                                      filter(!is.na(REC_ID))  %>%
                                      group_by(REC_ID) %>%
                                      filter(uncert_mean == 
                                               max(uncert_mean)) %>%
                                      mutate(uncert_max = uncert_mean) %>%
                                      dplyr::select(REC_ID, uncert_max, my))
                
                
                # in a last stage, we merge shocks_max_vol with 
                # rec_start_end:
                rec_start_end <- merge(x = rec_start_end, 
                                      y = rec_max, 
                                      by = "REC_ID", all = TRUE, na.rm=T)

                
                # actually, there is a slight problem, which we fix by
                # adding one month's numeric value to my_end:
                rec_start_end$my_end <- rec_start_end$my_end + (1/12)
                # this makes sure that if we refer to an end-month that
                # we assume that we are talking about the last day of that
                # respective month!
                
                # in the data.frame rec_start_end, the column 'my'
                # marks the month with the highest uncertainty
                # during a recession; we will mark these now by
                # joining with uncert_measures.stats.rec2
                uncert_measures.stats.rec2 <- left_join(
                  x = uncert_measures.stats.rec2,
                  y = rec_start_end[, c("REC_ID", "uncert_max", "my")],
                  by = "my")
                
                # we can then overwrite uncert_max with a 1
                uncert_measures.stats.rec2 <- uncert_measures.stats.rec2 %>% 
                  mutate(uncert_max = ifelse(is.na(uncert_max), 0, 1)) %>%
                  # at the same time we only select and continue the 
                  # columns that we actually need:
                  dplyr::select(my, uncert_mean, 
                                VXO, EPU_Hist, Michigan, Macro1, REC_ID.x,
                                REC_ID.y, uncert_max)
                
                
# (4) now that we have isolated the uncertainy-peaks in each recession,
# we want to extract the trajectories pre- and post those identified
# peak-points for all recessions!

# for this, we first search for the non-NA index-locations in column
# REC_ID.y:
NonNAindex <- which(!is.na(uncert_measures.stats.rec2$REC_ID.y))   
# and then loop through that vector:

# the column index of r'REC_ID.y' is:
col.index <- ncol(uncert_measures.stats.rec2) + 1


# test <- uncert_measures.stats.rec2
# j is the column index because we want to generate a dedicated
# column for each recession
j <- 0


for(i in NonNAindex){
  
        #we extract the value of REC_ID.y in the row i
        value <- uncert_measures.stats.rec2[i, col.index - 2]
        
        # we fill up the rows above and below that index
        # in the data.frame uncert_measures.stats.rec2
        if(i > 45){
          uncert_measures.stats.rec2[i-(3:38), col.index + j] <- "PRE"   
          # prior to recession
          uncert_measures.stats.rec2[i+(3:38), col.index + j] <- "POST"  
          # after recession
        }else{
          # we just go down until i
          uncert_measures.stats.rec2[i-(3:i), col.index + j] <- "PRE"
          uncert_measures.stats.rec2[i+(3:38), col.index + j] <- "POST"
        }
        j <- j + 1
}


# we rename the columns:
uncert_measures.stats.rec2 <- uncert_measures.stats.rec2 %>%
        dplyr::rename(REC1981 = V10,
                      REC1983 = V11,
                      REC1991 = V12,
                      REC2002 = V13,
                      REC2009 = V14)


# to reach a tidy data.format out of the above, we have to proceed
# in multiple steps:
uncert_measures.stats.rec2.tidy <- 
            bind_rows(
        uncert_measures.stats.rec2 %>%
        # we filter for the recession 1981
        filter(REC1981 == "POST" | REC1981 == "PRE") %>%
        # note that we do not take 'uncert_mean'
        # with us anymore because we actually do not need it
        dplyr::select(VXO, EPU_Hist, Michigan, Macro1,
                      REC1981) %>%
        # we rename REC1981 to 'PRE_POST'
        dplyr::rename(PRE_POST = REC1981) %>%
        # and create a new column 'REC_ID' that will hold
        # the year of the respective recession!
        dplyr::mutate(REC_ID = 1981) %>%
        # in a last step we have to use the 'gather' - command
        # to transform our columns to one variable:
        gather( 
          series, value, 
          -c(PRE_POST, REC_ID), na.rm=TRUE),
        
        # RECESSION 1983
        uncert_measures.stats.rec2 %>%
        # we filter for the recession 1981
        filter(REC1983 == "POST" | REC1983 == "PRE") %>%
        # note that we do not take 'uncert_mean'
        # with us anymore because we actually do not need it
        dplyr::select(VXO, EPU_Hist, Michigan, Macro1,
                      REC1983) %>%
        # we rename REC1981 to 'PRE_POST'
        dplyr::rename(PRE_POST = REC1983) %>%
        # and create a new column 'REC_ID' that will hold
        # the year of the respective recession!
        dplyr::mutate(REC_ID = 1983) %>%
        # in a last step we have to use the 'gather' - command
        # to transform our columns to one variable:
        gather( 
          series, value, 
          -c(PRE_POST, REC_ID), na.rm=TRUE),
      
        # RECESSION 1991
        uncert_measures.stats.rec2 %>%
        # we filter for the recession 1981
        filter(REC1991 == "POST" | REC1991 == "PRE") %>%
        # note that we do not take 'uncert_mean'
        # with us anymore because we actually do not need it
        dplyr::select(VXO, EPU_Hist, Michigan, Macro1,
                      REC1991) %>%
        # we rename REC1981 to 'PRE_POST'
        dplyr::rename(PRE_POST = REC1991) %>%
        # and create a new column 'REC_ID' that will hold
        # the year of the respective recession!
        dplyr::mutate(REC_ID = 1991) %>%
        # in a last step we have to use the 'gather' - command
        # to transform our columns to one variable:
        gather( 
          series, value, 
          -c(PRE_POST, REC_ID), na.rm=TRUE),

        # RECESSION 2002
        uncert_measures.stats.rec2 %>%
        # we filter for the recession 1981
        filter(REC2002 == "POST" | REC2002 == "PRE") %>%
        # note that we do not take 'uncert_mean'
        # with us anymore because we actually do not need it
        dplyr::select(VXO, EPU_Hist, Michigan, Macro1,
                      REC2002) %>%
        # we rename REC1981 to 'PRE_POST'
        dplyr::rename(PRE_POST = REC2002) %>%
        # and create a new column 'REC_ID' that will hold
        # the year of the respective recession!
        dplyr::mutate(REC_ID = 2002) %>%
        # in a last step we have to use the 'gather' - command
        # to transform our columns to one variable:
        gather( 
          series, value, 
          -c(PRE_POST, REC_ID), na.rm=TRUE),      

        # RECESSION 2009
        uncert_measures.stats.rec2 %>%
        # we filter for the recession 1981
        filter(REC2009 == "POST" | REC2009 == "PRE") %>%
        # note that we do not take 'uncert_mean'
        # with us anymore because we actually do not need it
        dplyr::select(VXO, EPU_Hist, Michigan, Macro1,
                      REC2009) %>%
        # we rename REC1981 to 'PRE_POST'
        dplyr::rename(PRE_POST = REC2009) %>%
        # and create a new column 'REC_ID' that will hold
        # the year of the respective recession!
        dplyr::mutate(REC_ID = 2009) %>%
        # in a last step we have to use the 'gather' - command
        # to transform our columns to one variable:
        gather( 
          series, value, 
          -c(PRE_POST, REC_ID), na.rm=TRUE)       
      )
      


# now we want to create an ID by REC_ID - Series - group:
uncert_measures.stats.rec2.tidy <- uncert_measures.stats.rec2.tidy %>%
                  group_by(REC_ID, series, PRE_POST) %>%
                  dplyr::mutate(ID = 
                                  case_when((REC_ID == 1981 &
                                            PRE_POST == "PRE") ~ as.numeric(row_number()+12),
                                  TRUE ~ as.numeric(row_number())))

# in a last step we make factors out of the actual
# factor-variables like 'PRE_POST', 'series' and 'REC_ID'
col_names <- c('PRE_POST', 'series', 'REC_ID')

uncert_measures.stats.rec2.tidy <- uncert_measures.stats.rec2.tidy %>%
          ungroup %>%
          dplyr::mutate(REC_ID = factor(REC_ID),
                        PRE_POST = factor(PRE_POST),
                        series = factor(series))

  
  

# in a last step we have to rebase everything to 100:
uncert_measures.stats.rec2.tidy.rebased <- 
                  as.data.frame(uncert_measures.stats.rec2.tidy %>%
                  group_by(REC_ID, series, PRE_POST) %>%
                  mutate(value = value + (100 - value[1L])))
uncert_measures.stats.rec2.tidy.rebased %>%
  dplyr::filter(PRE_POST == "POST" & series == "Michigan") %>%
  ggplot() +
     geom_line(
               aes(x=ID, y=value, colour=REC_ID))

# 


# we add a sequence to each data.frame
test2 <- ts(test2[2:ncol(test2)],
                             start = 1, frequency = 12)
test3 <- ts(test3[2:ncol(test3)],
            start = 1, frequency = 12)

autoplot(test2)
autoplot(test3)


plot(test2$REC2002, type="l")
lines(test2$REC2009, col="red")
lines(test2$REC1983, col="blue")
lines(test2$REC1991, col="green")
lines(test2$VXO, col="red")
lines(test2$Macro1, col="blue")
lines(test2$Michigan, col="green")
lines(test2$EPU_Hist, col="blue")
plot(test3$REC2002, type="l")
lines(test3$REC2009, col="red")
lines(test3$REC1991, col="green")
lines(test3$REC1983, col="blue")
lines(test3$VXO, col="red")
lines(test3$Macro1, col="blue")
lines(test3$Michigan, col="green")
lines(test3$EPU_Hist, col="blue")

                      
########################################################################
### PART x: correlation matrix between the uncertainty measures
########################################################################
comparison_measures_corr <- comparison_measures[c(2, 3, 4, 6)]

## correlation matrix between the time-series
corr <- cor(comparison_measures_corr, use = "complete.obs")
corr <- round(corr, 2)

stargazer(corr)

                      
                                            
