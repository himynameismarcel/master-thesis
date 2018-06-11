## algorithm from: https://stackoverflow.com/questions/22583391/
## peak-signal-detection-in-realtime-timeseries-data


## Haha nice to hear that! There are many ways to improve this algo, 
## so be creative (different treatment up/ down; median instead of mean; 
## robust std; writing the code as a memory-efficient function; 
## threshold margin so the signal doesn't switch too often, etc.).

ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag])
  stdFilter[lag] <- sd(y[0:lag])
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i])
    stdFilter[i] <- sd(filteredY[(i-lag):i])
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

# Data
y <- comparison_measures$Macro

lag       <- 30
threshold <- 2
influence <- 0.03

# Run algo with lag = 30, threshold = 5, influence = 0
result <- ThresholdingAlgo(y,lag,threshold,influence)
seq.y <- 1:length(y)

# Plot result
par(mfrow = c(2,1),oma = c(2,2,0,0) + 0.1,mar = c(0,0,2,1) + 0.2)
plot(1:length(y),y,type="l",ylab="",xlab="") 
lines(1:length(y),result$avgFilter,type="l",col="cyan",lwd=2)
lines(1:length(y),result$avgFilter+threshold*result$stdFilter,type="l",col="green",lwd=2)
lines(1:length(y),result$avgFilter-threshold*result$stdFilter,type="l",col="green",lwd=2)
## add dots to cases where result["signals"] == 1
points(seq.y[result$signals == 1], y[result$signals == 1], col="red", pch=19, cex = 0.5)
plot(result$signals==1,type="S",col="red",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)
#plot(result$signals==-1,type="S",col="black",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)



result["signals"]






### another approach:
### In signal processing, peak detection is often done via wavelet transform. 
### You basically do a discrete wavelet transform on your time series data. 
### Zero-crossings in the detail coefficients that are returned will correspond 
### to peaks in the time series signal. You get different peak amplitudes detected 
### at different detail coefficient levels, which gives you multi-level resolution.

### related to the above:
### this article: https://en.wikipedia.org/wiki/Change_detection
### and: https://stackoverflow.com/questions/12851208/how-to-detect-
### significant-change-trend-in-a-time-series-data/13660660#13660660




### another approach: 
### https://stats.stackexchange.com/questions/139660/detecting-changes-
### in-time-series-r-example

### You could use time series outlier detection to detect changes in time series. 
### Tsay's or Chen and Liu's procedures are popular time series outlier detection 
### methods . See my earlier question on this site.

### Some more background on the tsoutliers package:
### https://jalobe.com/blog/tsoutliers/

### the manual of the package can be found here:
### https://cran.r-project.org/web/packages/tsoutliers/tsoutliers.pdf

### Some more details on Chen and Lui's outlier detection package:
### https://stats.stackexchange.com/questions/104882/detecting-outliers-in-
### time-series-ls-ao-tc-using-tsoutliers-package-in-r-how
### R's tsoutlier package uses Chen and Liu's method for detection outliers. 
### SAS/SPSS/Autobox can also do this. See below for the R code to detect changes 
### in time series.
## Firstly I would like to say a big thank you to the author of the new tsoutliers package which implements Chen and Liu's time series outlier detection which was published in the Journal of the American Statistical Association in 1993 in Open Source software R.

## The package detects 5 different types of outliers iteratively in time series data:

## Additive Outlier (AO)
## Innovation Outlier (IO)
## Level Shift (LS)
## Temporary change (TC)
## Seasonal Level Shift (SLS)

library("tsoutliers")
dat.ts<- ts(comparison_measures$Macro,frequency=1)
data.ts.outliers <- tso(dat.ts)
data.ts.outliers
plot(data.ts.outliers)




#### Another approach 
### 'anomalize': a tidy anomaly detection algorithm that is time-based (built on top of
### tibbletime);
# http://www.business-science.io/code-tools/2018/04/08/introducing-anomalize.html
# The 'anomalize' - package is an open source package and consists of a scalable
# adaption of the Twitter AnomalyDetection package!

# a quick demo: 
# note: we needed to create a dedicated comparison_measures_test - data-frame
# that only consists of a date-column (with actual time-stamps) and 
# a second column which holds the Macro uncertainty-data!!!!
as.tbl(comparison_measures_test) %>%
  time_decompose(Macro) %>%
  anomalize(remainder, alpha=0.9, method = "gesd") %>%
  time_recompose() %>%
  plot_anomalies(time_recomposed = TRUE)

## the above package merges the pros and cons of three packages:
## * Twitter's AnomalyDetection package
## * Rob Hyndman'S forecast::tsoutliers() function available on through the
## forecast package on CRAN
## Javier Lopez-de-Lacalle's package, 'tsoutliers' on CRAN


## Frage: können wir einfach etwas nicht-parametrisches verwenden???
## oder lieber strucchange etc.?



## Another blog-post that is discussing mainly machine-learning algorithms
## https://www.datascience.com/blog/python-anomaly-detection
##  * K-nearest neighbor
##  * Local outlier factor (LOF)
##  * Clustering-Based anomaly detection (unsupervised learning): K-means
##  * Suppoert Vector Machine-Based Anomaly Detection;
## but then the blog-post implements a simple algorithm (MA)


## Another blog-post that is discussing anomaly detection:
## In data mining, anomaly detection (also outlier detection) is the identification
## of items, evens or observations which do not conform to an expected
## pattern or other items in a dataset. Typically the anomalous items will
## translate to some kind of problem such as bank fraud, a structural defect,
## medical problems or errors in a text.
## blog-post: https://bytefish.de/blog/anomaly_detection_with_r/
## and then discusses the twitter AnomalyDetection library
## description of the package: 
## 'AnomalyDetection is an open-source R package to detect anomalies which is robust, 
## from a statistical standpoint, in the presence of seasonality and an
## underlying trend. 
## installation:
install.packages("devtools")
devtools::install_github("twitter/AnomalyDetection")
library(AnomalyDetection)

# Seasonal Hybrid ESD (S-H-ESD) algorithm, which combines 
# seasonal decomposition with robust statistical methods to identify local 
# and global anomalies.
# Detect Anomalies:

# first, we need to transform the date-column in comparison_measures_test
# to POSIXct

# note: we need to have a single column data frame, list or
# vector as the input argument x in order to use anomalydetectorVec;
#
comparison_measures_test_ts <- as.vector(comparison_measures_ts[, 4])

res = AnomalyDetectionVec(comparison_measures_test_ts, 
                          direction='both', period= 12, plot=TRUE)
res$plot
# and we can see that only the extreme outliers are identified:
res$anoms



### Machine learning algorithms (for anomaly detection): 
###       *) time series clustering
###       *) Twitter'S Anomaly Detection package
###       *) Early Warning Signals Toolbox, which includes the earlywarnings package
###           (R)
###       *) deep learning for anomaly detection!



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
### Marcel Kropp, 29.05.2018
### This script discusses the uncertainty measures that we introduce in the text in
### greater detail.
### Everything that is produced here will end up in the part about 
### 'Identification of Shocks' in the main text.

### The goal of the script is the following:
### We want to apply certain peak- or shock-detection-methods to the various
### uncertainty measures we have introduced in the text.

### In particular we want to examine shocks by applying four approaches to the 
### time-series of all uncertainty measures:
###         (1) extend the creation of the 'Bloom-shock' to all uncertainty measures:
###             Bloom defines a shock as those events where the series values are 
###             more than 1.65 stds above the HP-detrended mean of the series.
###         (2) apply the procedure from the 'anomalize' - package to all
###             uncertainty measures
###         (3) apply the procedure from the package 'strcchange' to all
###             uncertainty measures and record the dates where a level shift
###             takes place!
###         (4) transform the uncertainty time-series into percentage changes
###             and use 'tsoutliers' to detect dates with extreme changes!

### All four approaches should produce two DISTINCT data-frame that then are further 
### used to create 
###         * Tables (for the Appendix) and 
###         * facetted time-series plots of all the
###           uncertainty measures and the identified 'shocks'/peaks.

### Because the shape of the data to be used for (*) a Table and (*) to mark episodes
### of unusually high uncertainty in the time-series plots are completely
### distinct, we have to cater for both goals throughout the below step-by-step
### procedure!

### Note:
### For the creation of the above, we restrict all time-series for all uncertainty
### measures that we consider until 06/2008 to work with the same time-span
### that Bloom has allegedly used (according to the Excel-replication-files
### that runs until 06/2008!).

### Note:
### In a final overhaul of this script we should check whether the original data
### we read in from Bloom (2009) with VOLATBL is identical to what we have
### generated ourselves!
### For the Macro Series (h1 and h12), MSoC and EPU there is no issue --> we just
### cut off the the series because even if there were updates we want
### to look back at the most up-to-date data! 

### We have to write a function that takes a vector with the time-series and
### then performs all the necessary steps to produce:
###   * a version of the necessary data that will be used for the TABLE!
###   * a version of the necessary data that will be used for the PLOT! 
### (alternatively we can also write two separate functions, but we will
### so how it goes further down below!)


######################################################################################

### First try:
### If we take the actual series that Bloom has provided in his replication
### files including the column VOLATBL, we should get exactly the same 
### results!

### Note: for the VOLATBL-series we still need a special treatment (once we have made
### sure that our generated VOLATBL is identical to Bloom's VOLATBL the below
### step is not necessary anymore!)
### All other measures can be extracted from the data.frame 'comparison_measures',
### because we took the data directly from the respective sources (and did not construct
### anything ourselves!)

### The data.frame 'comparison_measures' is produced in the R-script
### '12052018_thesis_master.R' but in the version there comes along
### without the separate variables for MONTH and YEAR!
### Therefore we replicate the generation of 'comparison_measures' here, call the
### resulting data-frame 'identify_shocks' (to avoid confusion) and add in the
### additional variables that we want/need.


### As mentioned above, the dataset that we read in in a preliminary version of this 
### script, is the dataset including the VOLATBL-series provided by Bloom himself, which
### we have stored in 'VARDATA_UPDATED_JFV.xlsx';
### We will replace this step with our own data as soon as we have clarified whether
### our own generated series and the series provided by Bloom is indeed the same!

#------------------------------------------------------------------------------------
#####################
## preliminaries
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
identify_shocks <- data.frame("macroUncert_h1" = macroUncertainty_index$h1,
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
#   * and finally the MSoC in a multi-setp join:
identify_shocks <- as.tibble(right_join(x = epu_index[, 
                                 c("News_Based_Policy_Uncert_Index","my")],
                                 y = identify_shocks,
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
                      right_join(x = michigan_survey[, 
                                                  c("veh_share_unc.MA","my")],
                                 y = .,
                                 by = "my") %>%
                      # further, we rename all the variables,
                      # select the ones we want to keep
                      dplyr::rename(VXO = mvol_VOLATBL,
                            EPU = News_Based_Policy_Uncert_Index,
                            EPU_Hist = EPU_Historical,
                            #GTU = GTU_US,
                            Macro1 = macroUncert_h1,
                            Macro12 = macroUncert_h12,
                            Michigan = veh_share_unc.MA) %>%
                      # we further reorder the variables and drop 'my' by 
                      # leaving it out!
                      # update: we need the 'my' variable again further down
                      # below: therefore we should not drop it!
                      dplyr::select(Date, year, month, VXO, EPU, EPU_Hist, Michigan,
                                    Macro1, Macro12))
                      
# filter data to belong to a certain 
# time-range (in our baseline case this
# is the 'Bloom-window' 
# (i.e., July 1962 - June 2008))
# identify_shocks.sub <- identify_shocks  %>%
#                       filter(Date <="2008-06-01" & Date >= "1962-07-01")
# alternative range: (i.e., July 1962 - June 2018))
identify_shocks.sub <- identify_shocks  %>%
                      filter(Date <="2018-06-01" & Date >= "1962-07-01")  

#####################
## function for creation of Bloom-Shock,
## input: data-frame with 'Date', 'year', 'month' and ONE 
##        distinct uncertainty-measure!
##        Note: the uncertainty measure has to enter the function last!
#####################

# first, we generate the data-frames that each only contains a distinct
# subset of time information and the series itself:
identify_shocks_VXO <- identify_shocks.sub %>%
                        dplyr::select(1:4) %>%
                        # we drop potential NAs
                        drop_na %>%
                        # we have to rename the uncertainty measure to 'series'
                        dplyr::rename(series = VXO)
identify_shocks_Macro1 <- identify_shocks.sub %>%
                        dplyr::select(c(1:3, 8)) %>%
                        # we drop potential NAs
                        drop_na %>%
                        # we have to rename the uncertainty measure to 'series'
                        dplyr::rename(series = Macro1)                        
identify_shocks_Macro12 <- identify_shocks.sub %>%
                        dplyr::select(c(1:3, 9)) %>%
                        # we drop potential NAs
                        drop_na %>%
                        # we have to rename the uncertainty measure to 'series'
                        dplyr::rename(series = Macro12)  
identify_shocks_EPU_Hist <- identify_shocks.sub %>%
                        dplyr::select(c(1:3, 6)) %>%
                        # we drop potential NAs
                        drop_na %>%
                        # we have to rename the uncertainty measure to 'series'
                        dplyr::rename(series = EPU_Hist)  
identify_shocks_Michigan <- identify_shocks.sub %>%
                        dplyr::select(c(1:3, 7)) %>%
                        # we drop potential NAs
                        drop_na %>%
                        # we have to rename the uncertainty measure to 'series'
                        dplyr::rename(series = Michigan)  
identify_shocks_EPU <- identify_shocks.sub %>%
                        dplyr::select(c(1:3, 5)) %>%
                        # we drop potential NAs
                        drop_na %>%
                        # we have to rename the uncertainty measure to 'series'
                        dplyr::rename(series = EPU)  

# for series that have some NAs(Michigan and EPU) we
# and need an additional function to fill up missing NAs 
# when merging all the data back in again (took the function from here:
# https://stackoverflow.com/questions/19074163/cbind-is-there-a-way-to-have-
# missing-values-set-to-na)
# the function that let's the user specify where NAs should be added to top
# or bottom:  1 means add NA to top of shorter vector
#             0 means add NA to bottom of shorter vector
# actually, it might be that we do not need the function anymore at the end
# (we'll see!)

my.cbind <- function(x,y,first) {
  
  if(length(y)<length(x[, 1])) {
    if(first==1) y = c(rep(NA, length(x[, 1])-length(y)),y);x=x
    if(first==0) y = c(y,rep(NA, length(x[, 1])-length(y)));x=x
  } 
  
  y <- as.data.frame(y)
  
  return(data.frame(x,y))
  
}


# the above data-frame will now serve as our example to generate the function
# into which we can then feed all other uncertainty-measures!

#-------------------------------------FUNCTION---------------------------------------
# here the function starts:
extractShocks_BLOOM <- function(df, name=name){
  
            nam <- NULL
            helper <- NULL
            mean_cycle <- NULL
            sd_cycle <- NULL
            helper.zoo <- NULL
            change <- NULL

            # (1) detection of shocks (i.e, 0/1)
            
            # detrend the constructed time-series:
            # because we know that our uncertainty-measure will enter last, we 
            # can do the following:
            # extract name for cycle-series which we will generate below:
            nam <- paste(colnames(df[4]), "cycle", sep="_")
            
            # apply hpfilter
            helper <- hpfilter(df[4],freq=129600)
            helper <- helper$cycle
            # store in data-frame:
            df <- cbind(df, helper)
            # rename newly added variable:
            colnames(df)[length(df)] <- nam
            
            
            # next, as mentioned above:
            # according to Bloom (2009) the shocks are chosen as those events 
            # where the series values are more than 1.65 standard deviations 
            # above the HP-detrended mean of the the series.
            # Further: in this setting, each month is being treated as an 
            # independent observation.
            
            # For us this means that we can calculate the mean and standard 
            # deviation of the very last column in our data-frame (which is 
            # the cycle-component):
            # note that we have to store the results in separate object:
            mean_cycle <- mean(df[, 5])
            sd_cycle <- sd(df[, 5])
            
            # Note: Bloom (2009) refers to the threshold which he selects as 
            # the 5% one-tailed significance level TREATING EACH MONTH AS AN 
            # INDEPENDENT OBSERVATION!
            
            # we create a threshold variable and store it in our dataset to 
            # be able to create an indicator variable that allows us to extract 
            # the 'Bloom' - shocks (this time across all variables) which we
            # want to have a look at!
            
            # to name this additional column also accordingly, we conserve our
            # name-object again:
            nam <- paste(colnames(df[4]), "thresh", sep="_")
            df[, nam] <- mean_cycle + 1.65*sd_cycle
            
            # now we can create an indicator-variable
            # (for this we conserve our 'nam' variable again):
            nam <- paste(colnames(df[4]), "shock", sep="_")
            df[, nam] <- ifelse(df[, 5] >
                                df[, 6], 1, 0)
            
            # having finally identified the shocks, we can proceed with 
            # the generation of a data.frame that holds the respective 
            # data 
            #     * for further use in a TABLE
            #     * and the data for further use in PLOTS
            # (see below!)

            

            # (2) creation of data.frame for PLOTS (if solo);
            # below under (3) we will create a data.frame
            # for the case that we have facetted data;
            # Having all the shocks now in our dataset allows us to 
            # create a dedicated data.farme that holds the respective 
            # start- end-dates of the shock-periods. but before we can 
            # do that, we need to CREATE 'start' and 'end'-dates for
            # the respective periods (which we will then ultimately 
            # use for plotting the episodes with shaded regions in our 
            # plots!)
            
            # To retrieve these respective 'start'- and 'end'-dates for the 
            # shock-periods we apply the below procedure that results in a 
            # dedicatd data.frame holding the variables 'start', 'end' for 
            # all shocks!
            
            # let us first filter for all rows where the shock = 1 to get 
            # an overview:
            # (and see whether it corresponds to the identified dates that 
            # Bloom also used himself!)
            df[df[, 7] == 1, c(1, 4, 7)]
            
            # continuation with procedure:
            # we want to give each episode of a shock a unique ID so that
            # we then can easily calculate the max volatility and first
            # volatility by grouping by ID;
            # to accomplish this, we decided for the following procedure:
            
            # we loop through all rows, and first check if the bloom_shock 
            # equals 1 or not:
            
            # we initialize the grouping variable x as follows:
            x <- 2
            
            # we add an empty column to our data frame:
            df["shock_ID"] <- NA
            
            # then we start the loop:
            for (i in 1:nrow(df)) {
              if(df[i, 7] == 1) {
                # if we detect a shock (i.e., 'bloom_shock' == 1), then 
                # we have to record a shock-ID:
                df[i, 8] <- x
                
              }else{
                # no shock
                x <- x + 1
              }
            }
            
            # we can now have a look at the newly created variable
            # (note that we have suppressed 'NA's for the creation 
            # of the table!):
            
            # the below command is just for ourselves to check!
            df %>%
                              filter(!is.na(shock_ID))  %>%
                              group_by(shock_ID) %>%
                              summarise(Count = n())
            # we create a variable 'my':
            df <- df %>%
                              mutate(my = year + month/12)

            
            # because we also want to add the 'change from previous month',
            # before the below step, we have to add an additional variable
            # that displays those changes;
            # specifically, we want to compute the first difference of the
            # series:
            
            # we create another 'nam' - variable:
            nam <- "change"
            
            # we want to add percentage changes as well!
            # note that we have to add
            # an NA at the beginning!
 
            helper.zoo <- as.zoo(df[, 4])
            perc <- function(x) round(diff(x) / x *100, 2)
            change <- perc(helper.zoo)
            
            # finally, we can add the column with the percentage
            # changes:
            df[, nam] <- c(NA, change)
            
            # we continue with a subset of the variables and only look
            # at the rows where we actually have a shock_ID other than 
            # NA!
            # for this we extract the name of the column we want to,
            # among others, filter for:
            
            df.sub <- as.data.frame(df %>%
                              filter(!is.na(shock_ID))  %>%
                              dplyr::select(shock_ID, year, month, my,
                                  series, change))
            
            # we save this version of df.sub because it will be used
            # in one of the steps below to create another data-frame
            # in the style of Bloom (2009), Table A.1 which we want
            # to export as well
            df.sub.blo <- df.sub
            
            # next, we construct a year-month-variable both out of
            # the pair 'year' & 'month';
            # note that we had to add a small correction by '0.083' manually
            # so that the yearmon is created the way we want it!
            df.sub$yearmon <- 
              as.yearmon(df.sub$my-0.083, "%Y-%B")
            
            # next we group_by shock_ID and extract the respective MAXIMUM
            # value of the uncertainty series in the group:
            # (note that the produced data-frame only has one row per shock_ID
            # and that we go back to the data.frame 'df'
            # for this operation:
            
            df.sub.max.vol <- 
                  as.data.frame(df %>%
                                filter(!is.na(shock_ID))  %>%
                                group_by(shock_ID) %>%
                                filter(series == max(series)) %>%
                                mutate(max_vol = series) %>%
                                dplyr::select(shock_ID, max_vol, my) %>%
                                # we add a column called max
                                mutate(max = 'm'))
            
            # note that we replace 'nam' for the name of the
            # column/variable we want to look at:
            df.sub.first.vol <- 
                  as.data.frame(df %>%
                                  filter(!is.na(shock_ID))  %>%
                                  group_by(shock_ID) %>%
                                  filter(row_number()==1) %>%
                                  mutate(first_vol = series) %>%
                                  dplyr::select(shock_ID, first_vol, my)%>%
                                  # we add a column called first
                                  mutate(first = 'f'))


            # now we can merge the three data-frames 
            # and see where we have max
            # and where we have first-values:
            df.sub <- 
                  (right_join(
                  x = df.sub.max.vol[, c("max", "my")], 
                  y = df.sub, 
                              by = "my") %>%
                  right_join(
                  x = df.sub.first.vol[, c("first", "my")], 
                  y = ., 
                              by = "my")) %>%
                  # in a last step we also drop 'year', 'month'
                  dplyr::select(-c(year, month)) %>%
                  # and we rearrange the columns
                  dplyr::select(yearmon, series, first, max,
                                change, shock_ID) %>%
                  # this is a tweak so that we don't paste
                  # NAs together in the 'unite' command below
                  replace_na(list(first = "", max = "")) %>% 
                  unite(label, first, max, sep=" ")

            # here we still should recalculate the shock_ID (starting from 1
            # again!)
            
            # in a last step, we want count the ID-column from scratch again!
            df.sub <- 
              df.sub %>%
              mutate(shock_ID = group_indices(., 
                                factor(shock_ID, levels = unique(shock_ID))))
            
            # lastly, we should round numeric columns in both
            # data-frames!
            # we set the number of decimal place to 2 for numeric columns:
            is.num <- sapply(df.sub, is.numeric)
            df.sub[is.num] <- 
                      lapply(df.sub[is.num], round, 2)
                       
            # in the very last step, we have replace the 
            # name 'series' with the name the user gives to
            # the function as a string!
            colnames(df.sub)[2] <- name

            
            # (3) creation of data.frame for PLOTS (if facetted);
            # above we had already test-wise filtered for all rows where 
            # the shock = 1 to get an overview:
            # this time we perform this operation again and save the
            # result; note that the data.frame 'df' grew in the meantim
            # (with additional variables) that we do not need, there
            # we combine the filter with a select-statement:
            df.dates <- df %>%
              filter(series_shock == 1) %>%
              dplyr::select(my) %>%
              # we create new variables 'Lower' and 'Upper'
              # and fill them with -Ind and + Inf
              # and at the same time we create a variable
              # called 'uncert_measure' containing the string
              # of the respective
              # uncertainty measure:
              # note: the name should be the same like the names
              # we chose for the data.frame 'comparison_measures'
              # so that a join can be performed afterwards more easily!
              dplyr::mutate(Lower = -Inf, Upper = + Inf, uncert_measure = name)
            
            
            # (4) creation of data.frame for TABLES:
            # to retrieve the respective 'start'- and 'end' 
            # - dates for the shock-periods,
            # we apply a procedure to the dataset 
            # 'df' that results in a 
            # separate dedicated data-frame with the following 
            # variables:
            # 'start', 'end', 'max_volatility', 'first_volatility' 
            # for all shock-measures!
            
            # we group_by shock_ID and extract the respective start-dates and
            # store the retrieved data to a new data-frame
            df.start <- as.data.frame(df %>%
                                   filter(!is.na(shock_ID))  %>%
                                   group_by(shock_ID) %>%
                                   filter(row_number()==1) %>%
                                   mutate(year_start=year, month_start=month, 
                                          my_start = my) %>%
                                   dplyr::select(shock_ID, year_start, 
                                                 month_start, my_start))
            
            # next, we replicate the above query for the end-dates
            # (note that in some scenarios start- and end-dates are identical!)
            df.end <- as.data.frame(df %>%
                                  filter(!is.na(shock_ID))  %>%
                                  group_by(shock_ID) %>%
                                  filter(row_number()==n()) %>%
                                  mutate(year_end=year, month_end=month, 
                                         my_end = my) %>%
                                  dplyr::select(shock_ID, year_end, 
                                                month_end, my_end))
            
            
            # # next, we can merge the two data-frames from above:
            df.start.end <- merge(x = df.start, 
                                      y = df.end, by = "shock_ID",
                                      all = TRUE, na.rm=T)
            
            # we re-arrange the sequence of columns
            df.start.end <- df.start.end %>%
                            dplyr::select(shock_ID, my_start, my_end) %>%
                            dplyr::rename(start = my_start,
                                          end = my_end)
            
            # the above data-frame now allows us to plot the episodes of 
            # high volatility into our previous plots (by adding a shaded
            # rectangle!)
            
            # actually, there is a slight problem, which we fix by adding
            # one month's numeric value to my_end:
            # ask Martin about this!
            # df.start.end$end <- 
            #             df.start.end$end + (1/12)
            # this makes sure that if we refer to an end-month that
            # we assume that we are talking about the last day of that
            # respective month!
            
            # in a last step, we want count the ID-column from scratch again!
            df.start.end <- 
              df.start.end %>%
              mutate(shock_ID = group_indices(., 
                              factor(shock_ID, levels = unique(shock_ID))))
            
            # lastly, we should round numeric columns in both
            # data-frames!
            # we set the number of decimal place to 2 for numeric columns:
            is.num <- sapply(df.start.end, is.numeric)
            df.start.end[is.num] <- 
              lapply(df.start.end[is.num], round, 2)
            
            # at the very end we add a column that holds the name of the
            # series we are looking at:
            df.start.end <- df.start.end %>%
                              dplyr::mutate(uncert_measure = name)
            
            
            
            
            
            # (5) creation of data.frame in the style of Bloom (2009)
            # Table A.1 (from Appendix)

            # we can continue working with the df.sub.blo 
            # which we had created above!
            
            df.sub.blo %>%
              filter(!is.na(shock_ID))  %>%
              group_by(shock_ID) %>%
              summarise(Count = n())
            
            # we group_by shock_ID and extract the respective start-dates 
            # and store the retrieved data to a new data-frame
            df.sub.blo.start <- as.data.frame(df.sub.blo %>%
                          filter(!is.na(shock_ID))  %>%
                          group_by(shock_ID) %>%
                          filter(row_number()==1) %>%
                          mutate(year_start=year, 
                                 month_start=month, 
                                 my_start = my) %>%
                          dplyr::select(shock_ID, year_start, 
                                        month_start, my_start))

            # next, we replicate the above query for the end-dates
            # (note that in some scenarios start- and 
            # end-dates are identical!)
            df.sub.blo.end <- as.data.frame(df.sub.blo %>%
                          filter(!is.na(shock_ID))  %>%
                          group_by(shock_ID) %>%
                          filter(row_number()==n()) %>%
                          mutate(year_end=year, 
                                 month_end=month, 
                                 my_end = my) %>%
                          dplyr::select(shock_ID, 
                                        year_end, month_end, 
                                        my_end))
            
            # next, we can merge the two data-frames from above:
            df.blo.start.end <- merge(x = df.sub.blo.start, 
                                      y = df.sub.blo.end, 
                                      by = "shock_ID", all = TRUE, na.rm=T)
            # and inspect the resulting data-frame:
            df.blo.start.end
            
            # we re-arrange the sequence of columns
            df.blo.start.end <- df.blo.start.end %>%
                                dplyr::select(shock_ID, year_start, 
                                        month_start, year_end, 
                                        month_end, my_start, my_end)
            
            
            # next, we construct a year-month-variable both out of
            # the pair 'year_start' & 'month_start' and 'year_end' & 
            # 'month_end'
            df.blo.start.end$yearmon_start <- as.yearmon(
                                    df.blo.start.end$my_start-0.083, "%Y-%B")
            df.blo.start.end$yearmon_end <- as.yearmon(
                                    df.blo.start.end$my_end-0.083, "%Y-%B")
            # next, we add a helper-column that we need in the 
            # next stage as well:
            df.blo.start.end$helper_date <- format(
                                    df.blo.start.end$yearmon_start-0.083, "%b")
            
            # this means that we can now drop the 'year' and 'month' 
            # variables
            df.blo.start.end <- df.blo.start.end %>%
              dplyr::select(-c(year_start, year_end, month_start, 
                               month_end))
            
            # next we add a column that gives us the duration of the 
            # shock:
            df.blo.start.end <- df.blo.start.end %>%
              mutate(duration = (df.blo.start.end$yearmon_end - 
                                   df.blo.start.end$yearmon_start) * 12 + 1) %>%
              dplyr::select(-helper_date)
            
            # next we group_by shock_ID and extract the respective 
            # MAXIMUM
            # VOLATILITY
            # (note that the produced data-frame only has 
            # one row per shock_ID)
            df.blo.max <- as.data.frame(df.sub.blo %>%
                                filter(!is.na(shock_ID))  %>%
                                group_by(shock_ID) %>%
                                filter(series == max(series)) %>%
                                mutate(max_vol = series) %>%
                                dplyr::select(shock_ID, max_vol, my))
            
            # we apply the same date-transformation as above:
            df.blo.max$yearmon_max <- as.yearmon(df.blo.max$my-0.083, "%Y-%B")
            # and we drop 'my':
            df.blo.max <- df.blo.max %>%
              dplyr::select(-my)
            
            # in a last stage, we merge df.blo.max with df.blo.start.end:
            df.blo.start.end.max <- merge(x = df.blo.start.end, 
                                          y = df.blo.max, 
                                      by = "shock_ID", all = TRUE, na.rm=T)
            
            # the above data-frame now allows us to plot the episodes of 
            # high volatility into our previous plots (by adding a shaded
            # rectangle!)
            
            # one addition:
            # knowing that we also have the data.frame df.sub.first.vol
            # that contains the first first volatilites, we add this
            # info to the above data-frame:
            df.blo.start.end.max <- merge(x = df.blo.start.end.max, 
                                          y = df.sub.first.vol, 
                                          by = "shock_ID", 
                                          all = TRUE, na.rm=T)
            
            # the column 'my' now represents the month with the first
            # volatilit:
            # we apply the same date-transformation as above:
            df.blo.start.end.max$yearmon_first <- as.yearmon(
                              df.blo.start.end.max$my-0.083, "%Y-%B")
            # and we drop 'my':
            df.blo.start.end.max <- df.blo.start.end.max %>%
                              dplyr::select(-c(my, first, my_start, my_end)) %>%
                              dplyr::select(shock_ID, yearmon_start,
                                             yearmon_end, duration, 
                                             max_vol, first_vol, 
                                             yearmon_max,
                                             yearmon_first)
            
            # in a very last step, we re-initialize the generation
            # of the shock_IDs:
            # in a last step, we want count the ID-column from scratch again!
            df.blo.start.end.max <- 
              df.blo.start.end.max %>%
              mutate(shock_ID = group_indices(., 
                            factor(shock_ID, levels = unique(shock_ID))))
            
            # in a very last step, we round all numbers to reasonable
            # number of decimals:
            df.blo.start.end.max[is.num] <- 
              lapply(df.blo.start.end.max[is.num], round, 2)
            
            # the below is only necessary if we would plot the above
            # information!
            # actually, there is a slight problem, which we fix by adding
            # one month's numeric value to my_end:
            # shocks_start_end$my_end <- shocks_start_end$my_end + (1/12)
            # this makes sure that if we refer to an end-month that
            # we assume that we are talking about the last day of that
            # respective month!
            
            
            # (6) return objects
            
            # lastly, we return the two data-frames that we want:
            # we want the two data-frames that are returned to 
            # get a 'suffix' according to the 'series'-string chosen
            # as the second argument of the function:
            
            # df1.out <- deparse(substitute(df.start.end))
            # df2.out <- deparse(substitute(df.sub))
            # 
            # df1.out <- paste(df1.out, series, sep="_")
            # df2.out <- paste(df2.out, series, sep="_")
            # 
            # df1.out <- df.start.end
            # df2.out <- df.sub
            
            # return
            return(list(df.sub, df.dates, df.start.end, df.blo.start.end.max))
            
}          
#-----------------------------------ENDOFFUNCTION-----------------------------------#
            
#####################
## running function for creation of Bloom-Shock,
## for all uncertainty-measures we have considered;
## the function 'extractShocks_BLOOM' creates two data-frames
## stored in a list for each of the uncertainty-measures!
#####################
# we apply the function to the various data.frames that we have created above:
VXO <- extractShocks_BLOOM(identify_shocks_VXO, "VXO")
Macro1 <- extractShocks_BLOOM(identify_shocks_Macro1, "Macro1")
Macro <- extractShocks_BLOOM(identify_shocks_Macro1, "Macro") 
Macro12 <- extractShocks_BLOOM(identify_shocks_Macro12, "Macro12") 
EPU_Hist <- extractShocks_BLOOM(identify_shocks_EPU_Hist, "EPU Hist")
Michigan <- extractShocks_BLOOM(identify_shocks_Michigan, "Michigan")
EPU <- extractShocks_BLOOM(identify_shocks_EPU, "EPU")

#####################
## Creation of TABLE
#####################
# finally, we can bind all extracted shocks together
# note that we have to reference the first element of the respective lists!
# Because we know that (after trial and error) that EPU and EPU_Hist
# have a fairly similar outcome, we drop EPU!
BLOOM_Shocks_Table.combined <- (full_join(x = Macro1[[1]], 
                                 y = VXO[[1]], 
                                      by = "yearmon") %>%
                                full_join(x = Macro12[[1]], 
                                      y = ., 
                                      by = "yearmon") %>%
                                full_join(x = EPU_Hist[[1]], 
                                          y = ., 
                                          by = "yearmon") %>%
                                # full_join(x = EPU[[1]], 
                                #           y = ., 
                                #           by = "yearmon") %>%
                                full_join(x = Michigan[[1]], 
                                          y = ., 
                                          by = "yearmon")) %>%
                                replace(is.na(.), "")


# we need to rearrange by yearmon again (after the join)
BLOOM_Shocks_Table.combined <- BLOOM_Shocks_Table.combined %>%
                          dplyr::arrange(yearmon) %>% 
                          # we remove all change - columns 
                          dplyr::select(-c(change, 
                                           change.x.x, 
                                           change.y, 
                                           change.x, 
                                           change.y.y)) %>%
                          # and lastly we have to rename all variables accordingly
                          dplyr::rename(label = label.x.x,
                                        #change = ,
                                        ID= shock_ID.x.x,
                                        #label = label,
                                        #change = ,
                                        ID = shock_ID,
                                        label = label.x,
                                        #change = ,
                                        ID = shock_ID.x,
                                        label = label.y,
                                        #change = ,
                                        ID = shock_ID.y,
                                        #label = label.x,
                                        #change = ,
                                        #shock_ID = shock_ID.x,
                                        label = label.y.y,
                                        #change = ,
                                        ID = shock_ID.y.y)


# we try to export to a latex table
latextable(BLOOM_Shocks_Table.combined, 
           colnames = c("Date", "Michigan", "Label", "ID",
                        "EPU Hist",  "Label", "ID", "Macro12",
                        "Label", "ID", "Macro1", "Label",
                        "ID", "VXO", "Label", "ID"))


# we can also export the smaller table only to be able to 
# compare it to Bloom's Table (2009), Appendix A.1
latextable(VXO[[4]], colnames = c("Shock ID", "Date Begin",
                                  "Date End", "Duration", "Max Volatility",
                                  "First Volatility", "Date Max Vol",
                                  "Date First Vol"))


#####################
## Creation of PLOT (facetted)
#####################
VXO[[3]]
Macro1[[3]]
Macro[[3]]
Macro12[[3]]
EPU_Hist[[3]]
Michigan[[3]]


# we need to stack all second elements of the lists above each other
# to have all shock-periods combined:
shocks_stacked <- bind_rows(VXO[[3]], Macro[[3]], EPU_Hist[[3]], Michigan[[3]])

# we remove the column 'shock_ID' and add the columns
# 'ymin' and 'ymax' holding -Inf and +Inf:
shocks_stacked <- shocks_stacked %>%
                        dplyr::mutate(ymin = -Inf,
                                      ymax = +Inf) %>%
                        dplyr::select(-shock_ID)

# from the script '12052018_thesis_master.R', the object
# comparison_measures contains all time-series we want to
# have a look at;
comparison_measures

# for plotting purpose we restrict the data to the window
# 1962/07 - 2008/06 (the window that Bloom looked at)
# comparison_measures.sub <- comparison_measures %>%
#                       filter(my >= 1962.583 & my <= 2008.500)
comparison_measures.sub <- comparison_measures %>%
                      filter(my >= 1962.583 & my <= 2016.500)

# to bring the dataset into a useful shape for our purposes,
# we 'gather()' the variables with the variable-Names
# together;
# note that before applying the 'gather' - function, we remove
# the series for 'EPU' (because we only want to keep EPU Historical!)
comparison_measures.sub.stacked <- comparison_measures.sub %>%
                        dplyr::select(-EPU) %>%
                                      gather( 
                                      uncert_measure, value, 
                                      -my, na.rm=TRUE)

# note: 
# we have designed the handling of dates in such a way that
# e.g., Dec 2015 results in a numeric value of 2016.00.
# If we now would want to plot a shock that only lasted
# for one month, start and end would be the same and ggplot
# would not pick it up! 
# Further, the duration of shocks would be slightly misread
# because in a scenario like the above, we would think that the
# shock 'only' started with the beginning of 2016 while it
# actually started in ec 2015. 
# But to be consistent with our handling of 'my' and considering
# that if my = 2016.00 the actualy date would refer to Dec 2015,
# we must take into account that we have introduced such a slight
# shift also when looking at the actual time-series;
# therefore, to be consistent, we move all end-dates up
# by 1/12!

shocks_stacked$end <- shocks_stacked$end + (1/12)



# besides the generation of 'shocks_stacked' which holds the PERIODS
# of significant volatility, we also want to plot the actual
# Bloom-Shocks for each measure (i.e., the month with maximum
# volatility):
# therefore, we repeat the above procedure also for the fourth
# data-frame that comes out of our function:
VXO[[4]]
Macro1[[4]]
Macro[[4]]
Macro12[[4]]
EPU_Hist[[4]]
Michigan[[4]]


# we need to stack all second elements of the lists above each other
# to have all shock-periods combined:
max_shocks_stacked <- bind_rows(VXO[[4]] %>%
                        dplyr::mutate(uncert_measure = "VXO"), 
                      Macro[[4]] %>%
                        dplyr::mutate(uncert_measure = "Macro"), 
                      EPU_Hist[[4]] %>%
                        dplyr::mutate(uncert_measure = "EPU Hist"), 
                      Michigan[[4]] %>%
                        dplyr::mutate(uncert_measure = "Michigan")) %>%
                # at the same time we remove all columns which we
                # we don't need further down below
                # (better: only keep the ones we need)
                      dplyr::select(yearmon_max, uncert_measure, max_vol) %>%
                # rename 'yearmon_max' to 'yearmon_max_start'
                      dplyr::rename(yearmon_max_start = yearmon_max) %>%
                # and add a column that is yearmon_max + 1/12 
                # for plotting purposes:
                      dplyr::mutate(yearmon_max_end = yearmon_max_start + 1/24)
                

# we add the columns
# 'ymin' and 'ymax' holding -Inf and +Inf:
max_shocks_stacked <- max_shocks_stacked %>%
                    dplyr::mutate(ymin = -Inf,
                                  ymax = +Inf)


# now we should be able to plot the result:
# note: the below plot marks both entire periods with 'significant'
# volatility;
# 
BLOOM_Shocks_plot_combined <- ggplot(comparison_measures.sub.stacked, aes(x=my,
                                        y=value)) +
  geom_line() + 
  facet_grid(uncert_measure ~ ., scales="free_y") +
  geom_rect(data=shocks_stacked, aes(x = NULL, y=NULL, xmin=start, xmax=end,
                                     ymin=ymin, ymax=ymax), 
                                     alpha=0.2, fill="red") +
  geom_rect(data=max_shocks_stacked, aes(x = NULL, y=NULL, 
                                         xmin=yearmon_max_start, 
                                         xmax=yearmon_max_end,
                                         ymin=ymin, ymax=ymax),
                                        fill="red") +
  # following the latest paper of Jurado et al, 
  # we add dots that identify the shock-events 
  # (in addition to the blue vertical lines!)
  geom_point(data=max_shocks_stacked, aes(x = yearmon_max_start, y=max_vol), 
             alpha=0.9, colour="red") + 
  scale_y_continuous(name = "Uncertainty Measures") +
  scale_x_continuous(name = "Year", 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  theme(legend.position="bottom") + 
  labs(color=NULL) +
  theme(legend.position="none",
        axis.text=element_text(size=10),
        #plot.title = element_text(size=10, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
        #legend.text=element_text(size=14),
        #axis.text.x=element_blank(),
        #plot.margin = unit(c(0,0,0,0), "mm")
        #panel.grid.major = element_blank()) 
        #panel.grid.minor = element_blank())
        #strip.background = element_rect(colour="black", fill="white", 
         #                               size=1.5, linetype="solid"))


BLOOM_Shocks_plot_combined

# a dedicated function to add sublabs to ggplot plots:
# see https://stackoverflow.com/questions/44616530/axis-labels-on-two-
# lines-with-nested-x-variables-year-below-months

# ggsave(file="BLOOM_Shocks_plot_combined.pdf")
ggsave(file="BLOOM_Shocks_plot_combined_all2016.pdf")



# Testing for structural breaks
# lines(confint(break_point, breaks = 1))
# 
# test <- ts(identify_shocks.sub$VXO, start=c(1962, 7), end=c(2008, 6), frequency=12)
# bp.test <- breakpoints(test ~ 1, h = 0.005)
# plot(test)
# par(mfrow = c(1, 1))
# plot(test)
# lines(fitted(break_point, breaks = 1), col = 4)
# lines(fitted(bp.test, breaks = 1), col = 4)
# lines(fitted(bp.test, breaks = 6), col = 4)
# lines(confint(break_point, breaks = 1))
# lines(confint(bp.test, breaks = 1))
# bp.test <- breakpoints(test ~ 1, h = 0.002)
# bp.test <- breakpoints(test ~ 1, h = 0.004)
# plot(test)
# summary(bp.test)
# lines(fitted(bp.test, breaks = 30), col = 4)
# ci.test <- confint(bp.test)
# ci.test
# lines(ci.test)
# 
# test2 <- ts(identify_shocks.sub$Macro1, start=c(1962, 7), end=c(2008, 6), frequency=12)
# 
# 
# # fit an AR(1) model to time series
# test.ar <- ar(test, aic = TRUE, order.max = 1,
#    method = "ols")
# test2.ar <- ar(test2, aic = TRUE, order.max = 1,
#    method = "ols")
# 
# # the coefficients can be extracted with 'ar' from the model:
# h <- - log(2) / log(abs(test.ar$ar))
# h2 <- - log(2) / log(abs(test2.ar$ar))




