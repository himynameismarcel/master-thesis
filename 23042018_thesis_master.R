########################################################################
### Marcel Kropp, 15.04.2018
### This script is separated into 2 big parts:
### PART 1: deals with the construction of the 'Bloom-shock', following 
###        Bloom (2009) thereby making use of 
###         (*) VXO-data from CBOE for the period where data is available 
###             (i.e., 1986 - 2018)
###             and 
###         (*) actual monthly volatilites (for the period pre-1986) 
###             calculated as the monthly standard 
###             deviation of the daily S&P500 index (normalized to the same 
###             mean and variance as the VXO index for when the two
###             series overlap, i.e., from 1986 - 2003).
###         Overall, the combined volatility measure will span from
###         June/1962 until 03/2018.
###         While the ultimate goal is the construction of the 'Bloom-shocks',
###         the script in PART 1 also generates various plots and tables
###         along the way.
###         The respective steps are outlined in detail below.
### PART 2: deals with the estimation of impulse-response-functions
###         following Jordá (2005);
###         Thereby the constructed dataset containing, among others,
###         the 'Bloom-shocks' is extended with additional 
###         regressors and (a) dependent variable(s) representing output
###         (most likely we will refer to industrial production for this!).
###         Following Bloom's (2009) VAR-estimations and the length
###         of our uncertainty-measure (i.e., the volatility series as well
###         as the derived 'Bloom-shocks'), the estimations use data
###         ranging from June/1962 - March/2018.
########################################################################

########################################################################
### PART 1: construction of 'Bloom-shock' following Bloom (2009)
########################################################################

# Note to self: a very nice article to read is:
# https://rpubs.com/bradleyboehmke/data_wrangling

## PART 1 is separated into six parts:
## (1.1)  1.1 handles the sp500 historical data (for which 
##        Bloom (2009) has not exactly declared the source from which 
##        he got the data from).
##        Note that for the longer historical sp500 data, Bloom (2009) 
##        has mentioned Bill Schwert (1989) as source). 
##        We have not changed the data in any respect but would need to 
##        clarify its source because a web-search did not reveal where 
##        the data originally comes from.
##        The data in sp500.csv ranges from 196201 to 200312.
## (1.2)  1.2 handles the raw volatility figures downloaded directly from
##        the homepage of the CBOE (for the purpose of our replication 
##        below, we have downloaded the entire data available from the 
##        homepage of the CBOE).
##        The .csv-file called 'vxo_new.csv' that we read in below, is 
##        constructed out of two files which are available at 
##        http://www.cboe.com/products/vix-index-volatility/vix-options-and
##        -futures/vix-index/vix-historical-data
##        Because Bloom (2009) uses VXO-data (and not the newer VIX), we have 
##        accordingly downloaded the 
##        following two files:
##                * 'Old methodology: VXO data for 1986 - 2003' 
##                  (downloaded as 'vxoarchive.xls')
##                * 'Old methodology: VXO data for 2004 to present 
##                  (Updated Daily) 
##                  ('downloaded as 'vxocurrent.csv')
##      Both above files come along with the variables 'Date', 'Close', 
##      'Open', 'High' and 'Low'.
##      A comparison with the VXO data that Bloom used himself 
##      (file vxo.csv) reveals, that Bloom (2009) used the Close-Price of 
##      the series. 
##      We therefore also only keep the columns 'Date' and 'Close' and 
##      stack the data of the two files from above together to get
##      one VXO volatility series ranging from 1986 - 2018, stored in a 
##      file called 'vxo_new.csv'.
##              Note: we had to apply the following corrections to 
##              'vxo_new.csv':
##                              * clean data for 10/18/2000 (two numbers 
##                                in cell for Close)
##                              * remove empty row between 12/31/1999 
##                                and 01/03/2000
##                              * remove empty row between 6/5/2006 
##                                and 6/6/2006
##              These corrections were, for now, done manually;
##              For a later version of the script it would be great 
##              to have these checks and corrections run 
##              automatically!
## (1.3)  1.3 performs the merge (join) of the two datasets
##        and all related manipulations with the goal to normalize the 
##        returns for the sp500 series for the period where the two series
##        (sp500 volatility and VOX) overlap.
##        At the end of this stage, the dataframe that ultimately contains
##        the monthly volatility figure from 1962 - 2018 is generated.
## (1.4)  1.4 deals with a first visual representation of the volatility 
##        time-series.
##        (Note that we not yet include the actual historical events into
##        our figure [this wil be left to our figure below after we have
##        identified the Bloom-shocks!])
## (1.5)  1.5 finally creates the actual 'Bloom-shocks' (i.e., an indicator
##        variable 0/1). According to Bloom (2009), shocks are selected as
##        those volatilities that are 1.65 SDs above the
##        HP-filtered trend. Accordingly, we will construct the HP-filtered
##        trend classify shocks according to the above stated criterion.
##        In his original paper, Bloom (2009) has identified 17 shocks with
##        this technique; we should expect to detect additional shocks for
##        the post-2009 period until today which we would
##        then include in our analysis!
## (1.6)  In 1.6 we repeat the plot from (1.4) and add vertical bars
##        for the identified 'Bloom-shocks' as well as try to identify 
##        the corresponding political/economic event.

###############################
## loading libraries
###############################
options(stringsAsFactors = F)
# When loading the packages note that according to
# http://stat545.com/block013_plyr-ddply.html
# we should always load 'plyr' before 'dplyr' when loading/using both.
library(plyr)
library(dplyr)
library(tibble)
library(tidyr)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())
library(mFilter)
library(zoo)
library(miscFuncs)
library(dynlm)
library(sandwich)
library(lmtest)
library(xlsx)

# clear workspace
rm(list = ls())

###############################
## (1.1) s&p500 historical data
###############################
## Generate baseline stock volatility data (to merge with
## VXO data)
## based on sp500.csv (data ranges from 1962 - 2003, daily
## data!)

# first we read in the historical sp500 data (for which Bloom
# does not exactly declare where he's
# got it from)
sp500 <- read.csv(file="sp500.csv", header=TRUE, sep=",")
# this file contains the following variables:
#           * caldt.............the calendar date in the format
#                               YYYYMMDD
#           * year..............YYYY
#           * month.............M
#           * day...............D
#           * ym................in the 'readme.txt' - file of 
#                               Bloom (2009) this variable
#                               is declared as 'ym=year+(month-1)/12;
#                               my impression is that it is rather 
#                               'ym=year+month/12.
#           * date..............seems to be linked to 'caldt' in that
#                               is is a
#                               sequency of integers starting at 913
#                               and going up
#                               in steps according to 'caldt'
#           * index............the value of the index (S&P 500)
#           * ri...............the daily return of the index in %

# after reading in, our dataset is stored in a data.frame; this can
# be seen when calling the 'class()' - function
class(sp500)

# str(sp500) reveals that the variables have been assigned either 
# the type 'int' or 'num' (i.e., the date-variable 'caldt' has not 
# been transformed to a date)

# we convert the data frame to a tibble
sp500 <- as_data_frame(sp500)

# next, we calculate the standard deviation for all 'ri' (daily returns of the index in %)
# for each year-month (i.e., BY the variable 'ym'!);
# the most convenient function to do this is the 'ddply' function (see the reference
# here: https://www.rdocumentation.org/packages/plyr/versions/1.8.4/topics/ddply);
# I found this solution suggested here: 
# https://stackoverflow.com/questions/15467219/calculate-group-characteristics-without-ddply-and-merge#comment21890299_15467219
# For a detailed treatment of the 'plyr' - package, see Hadley's article: 
# file:///Users/marcelkropp/Downloads/v40i01.pdf.
# Also, the following blog-post gives an overview: http://stat545.com/block013_plyr-ddply.html.
# In essence: 'ddply()' accepts a data.frame, splits it into pieces based on one or more
# factors, computes on the pieces, then returns the results as a data.frame.
# (the base R functions most relevant to 'ddply()' are 'tapply()', etc.)

sp500 <- ddply(sp500, ~ym, mutate, sdri = sd(ri, na.rm=T))
# Note that we have added 'na.rm=T' to the 'sd()' - function to make sure that
# also for 196207 the standard deviation is calculated which has an NA for the
# the very first value (otherwise all entries for the month 07 would end up
# being NA).

# finally, we order the data frame (actually: tibble) by the variable 'date' 
# (which in our case is a peculiar variable but which fulfills the same purpose
# like the actual 'date' - variable 'caldt')

sp500 <- as_data_frame(sp500[with(sp500, order(date)), ])

###############################
## (1.2) VXO - data from CBOE
###############################
## Here we process the raw volatility figures which we have downloaded directly
## from the CBOE's homepage and which we have stored in 'vxo_new.csv'
## (data ranges from 1986 to 2018, daily data!)

# first we read in the VOX-data
vxo_new <- read.csv(file="vxo_new.csv", header=TRUE, sep=",")
# vxo_new <- read.csv(file="vxo.csv", header=TRUE, sep=",")

# we check out the structure of 'vxo_new'
str(vxo_new)
# as we can see, the variable 'vol' has already been read-in as 'numeric',
# so everything is okay

# next, we want to split the variable 'date' into three variables 'month',
# 'day' and 'year';
# for this, we first need to convert the 'date' - variable (which is currently 
# declared as a factor) into a DATE datatype:
vxo_new$date <- as.Date(vxo_new$date, format = "%m/%d/%Y")
# we change the format of the date variable
vxo_new$date <- format(vxo_new$date, "%m/%d/%Y")
# note that the above command converts the date-type to a character again
# (we will deal with this later!)

# next, we create the three variable 'month', 'year' and 'day' using
# the 'separate()' - function from the 'tidyr' - package
vxo_new <- separate(vxo_new, "date", c("month", "day", "year"), sep = "/", 
                    remove=FALSE, convert=TRUE)
# as we can see, with the option 'convert=TRUE', the values are automatically
# transformed into integers!

# having sliced our date variable into month, day and year,
# we can drop the original variable 'date' (actually, we could
# have dropped 'date' automatically by keeping 'remove=TRUE', but
# we preferred it this way to double-check the output!)

# we convert our date-variable back to a date datatype
vxo_new$date <- as.Date(vxo_new$date, format = "%m/%d/%Y")

# next, we want to apply a consistency check to make sure that there are no trading days
# on Sat or Sun in our data!
# for this, we create two helper variables: 'm2' = 'month' and 'd2'='day', date2=date
vxo_new <- mutate(vxo_new, m2 = month, d2=day, date2=date)
# and then extract the day of the week using the 'wday' component of a 
# 'POSIXlt' - object (because it is numeric, starting on Sunday)
# as opposed to the normal 'weekdays' - function
vxo_new$dow <- as.POSIXlt(vxo_new$date)$wday
# Note: also lubridate sohuld have a 'wday' componsent (check out later!)

# next, Bloom (2009) performs two very strange commands that actually
# do not make much sense (have to ask Martin about this!)

# Above, we have created the variable 'date2' as a copy of 'date';
# In Bloom's stata-code, the created variable 'date' (which was made up of
# month, day and year) is a date in elapsed date format.
vxo_new$date2 <- as.numeric(vxo_new$date2)
# Note: for what follows next, we cannot use the elapsed time values
# that Bloom uses, because 

# we drop the 'day' - variable (which we had actually generated above;
# but because we use a different strategy for the generation of dates
# than Bloom, we could have also easily just skipped it!);
# edit: we decided to keep the 'day' - variable to have a better insight
# into whther or not our 'merge' of the two datasets below worked or not!
# vxo_new <- vxo_new %>% select(-day)

# next, we fill in 9/11 figures using European exchange extrapolation
# the following filter-statement shows us where we would need to input
# the data:
vxo_new %>% filter(date2>11574 & date2 < 11590)

# it is exactly between 2001-09-10 where for six days until 2001-09-17
# there is no data available;
# our goal is to add additional rows to our data-frame, remove weekends
# and then overwrite the vol-value with the value from the European
# exchange!

# first we generate an additional variable that gives us the numbers of the
# row that we are currently in in our data frame:
vxo_new$row <- seq.int(nrow(vxo_new))

# running our above filter-statement again, we know in which rows we want to 
# replicate:
vxo_new %>% filter(date2>11574 & date2 < 11590)

# we replicate row number 3958 that holds the data for 2001-09-10
vxo_new <- rbind(vxo_new, vxo_new[rep(3958, 6), ])
# so in our case, this command creates 9 additional duplicates of 2001-09-10

# we order by date
vxo_new <- vxo_new[with(vxo_new, order(date)), ]

# filtering with the below command would show us that everything
# worked as expected
vxo_new %>% filter(date2>11574 & date2 < 11590)

# we need to replace the dates between 2001-09-10 and 2001-09-17
# with their actual values!
vxo_new$date[vxo_new$row==3958] <- seq(as.Date("2001-09-10"), as.Date("2001-09-16"), by="days")

# we do the same with the the date2-column
vxo_new$date2[vxo_new$row==3958] <- seq(11575, 11581, 1)
# and then re-create the row-column
vxo_new$row <- seq.int(nrow(vxo_new))

# then we replace the value for vol from 2001-09-11 until 2001-09-16 (inclusive)
# with the value from the European exchange
vxo_new$vol[vxo_new$date2>11575 & vxo_new$date2<11582] <- 58.2

# and re-create the dow-variable
vxo_new$dow <- as.POSIXlt(vxo_new$date)$wday
# to be able to drop weekend days
vxo_new <- vxo_new %>%
            filter(dow != 6 & dow != 0)

# and we sort again by date
vxo_new <- vxo_new[with(vxo_new, order(date)), ]

# next, we want to drop months with less than 12 days of data:
# for this purpose we first create a variable 'count' that counts the number
# of days per month
vxo_new <- vxo_new %>% group_by(year, month) %>% mutate(dayc = n())
# we convert it back to a data frame
vxo_new <- as.data.frame(vxo_new)

# in our cae, the below command causes the observations for april 2018
# to be deleted because we only have 8 observations when we downloaded
# the data;
# therefore, depending on when we run the final code, it might be 
# that the month that we are currently in will still be included!
vxo_new <- vxo_new %>%
            filter(dayc>10)

# next, we drop these two helper variables (dayc, date2) because we do not
# need them anymore in what follows:
# (note that date2 was a helper variable which we had used in the expansion
# of our dataset when we wanted to include the data from the European exchanges);
# note that we also drop the variable 'd2' and 'm2' because we could not 
# really figure out the usage of them in Bloom's Stata-code;
# note that because we have a 'clash' of the select-statement of the
# dplyr-package with another package, we have to explicitly tell R to 
# use the version from the 'dplyr' - package!
# alernatively we could get rid of having to type 'dplyr::select' 
# by executing 'select <- dplyr::select'
vxo_new <- vxo_new %>% dplyr::select(-c(dayc, date2, d2, m2))

###############################
## (1.3) merge s&p500 with VXO - data from CBOE
##     and normalize returns over period where series overlap
###############################
# there are many options of how to perform the join between the two
# datasets (for some input see the stackoverflow discussion here:
# https://stackoverflow.com/questions/1299871/how-to-join-merge-data-frames-inner-outer-left-right);

# we decided for the dplyr-way:
# before applying the dplyr-way of a join, we have to make sure that the
# columns by which we want to join are both type numeric
vxo_new$date <- as.numeric(vxo_new$date)

# we add an extension: because we cannot interpret the meaning of 'date' - variable in the sp500
# dataset, we copy the variable, transfer it to a date-datatype and then to a numeric
# date datatype (which counts the elapsed numer of days since 1970!)
# first we copy date to date2 (to free the name 'date' for the actual variable
# that we want to create!)
sp500$date2 <- sp500$date
# then we overwrite date with the values from caldt (which are still integers!)
sp500$date <- sp500$caldt
# we make a date out of 'date'
sp500$date <- as.Date(as.character(sp500$date), "%Y%m%d")
# and then we make a numeric date-variable out of it
sp500$date <- as.numeric(sp500$date)

# the below command performs an outer join of the two datasets
sp500_merge_vxo <- merge(x = vxo_new, y = sp500, by = "date", all = TRUE, na.rm=T)
# an inspection of the new data-frame sp500_merge_vxo reveals:
str(sp500_merge_vxo)
# we have year.y and month.y which are the year and month variables from the sp500 dataset,
# we have year.x and month.x which are the year and month variables from the vxo_new dataset;
# to explain the many NAs (be aware that we performed an outer join),
# we have to keep in mind that the sp500 dataset ranges from 07/1962 until 12/2003 and
#                         that the vxo_new data ranges from 01/1986 until 04/2018;
# this means that in the period 07/1962 until 12/1985 we only have sp500 data,
#                 in the period 01/1986 until 12/2003 we have BOTH sp500 and vxo_new data,
#             and in the period 01/2004 until 04/2018 weonly have vxo_new data!

# to execute the below commands, we need a 'year' - and 'month' - variable
# that spans the entire sp500_merge_vxo dataset!
sp500_merge_vxo <- sp500_merge_vxo %>%
                  mutate(year = case_when(is.na(year.x) ~ year.y,
                          !(is.na(year.x)) ~ year.x))
sp500_merge_vxo <- sp500_merge_vxo %>%
                  mutate(month = case_when(is.na(month.x) ~ month.y,
                           !(is.na(month.x)) ~ month.x))

# after the above merge, the below command drops only 7 observations after
# 1990 that are only in the sp500 dataset (i.e., have no match with vox_new.csv
# after 1990!)
# exactly 7 such cases exist:
sp500_merge_vxo %>%
                  filter(!complete.cases(vol) & year > 1990)

sp500_merge_vxo <- sp500_merge_vxo %>%
                  filter(!(!complete.cases(vol) & year > 1990))

# the next command creates a moving average of the vxo-volatility value.
# this calculation is necessary because currently the vol-values from the
# vxo-dataset are DAILY figures (whereas the sdri-values are 
# the standard deviation of the daily returns over a month!)
sp500_merge_vxo <- as.data.frame((sp500_merge_vxo %>% 
                  group_by(year.x, month.x) %>% mutate(mvol = mean(vol))))

# we check whether indeed the calculation worked:
sp500_merge_vxo %>%
                  filter(complete.cases(vol))

# for the below command note that sdri = standard deviation of the daily returns
# over a month/year (as stated above!)
# the below command annualizes our volatility (standard deviation) of the monthly
# figures (ask Martin why Bloom used the square root of 366???)
sp500_merge_vxo$sdri <- sp500_merge_vxo$sdri * sqrt(366)

# the entire below sequence normalizes the two series (of vol and sdri)
# to have the same mean and variance (excluding the Black-Monday values along
# the way because they are so extreme!)

# to fully understand the sequence of commands below, we have to remind
# ourselves again of the following:
# we have to keep in mind that the sp500 dataset ranges from 07/1962 until 12/2003 and
#                         that the vxo_new data ranges from 01/1986 until 04/2018;
# this means that in the period 07/1962 until 12/1985 we only have sp500 data,
#                 in the period 01/1986 until 12/2003 we have BOTH sp500 and vxo_new data,
#             and in the period 01/2004 until 04/2018 weonly have vxo_new data!

# we create the variable 'miss' that = '1' for cases where 'mvol' is missing (i.e., NA)
#                                 and = '0' for cases where mvol is not missing;
# (remember that 'mvol' is the moving average of the 'vxo' volatility value by
# year and month!)

sp500_merge_vxo <- sp500_merge_vxo %>%
                    mutate(miss = case_when(is.na(mvol) ~ 1,
                                                    !is.na(mvol) ~ 0))
# note that 'mvol' is the moving average of the vxo volatility value
# by year-month

# the next command sets 'miss' to 1 for the months 10/11/12 in 1987
# and the month 1 in 1988 (to remove the period about Black Monday);
# note to self: because we changed the way we calculate
# the values to normalize the two time-series, we have decided to 
# give the series miss the code=2 for Black Monday and the following
# months!

sp500_merge_vxo <- sp500_merge_vxo %>%
                    mutate(miss = replace(miss, 
                                  (year.y==1987 & month.y %in% c(10, 11, 12)) |
                                        (year.y==1988 & month.y==1), 2))

# a subsequent filtering should show that all 0s were set to 2s:              
head(sp500_merge_vxo %>%
                    filter((year.y==1987 & month.y %in% c(10, 11, 12)) |
                     (year.y==1988 & month.y==1)))


# currently the variable 'miss' only consists of 0s and 1s as follows:
#             * '1': for the period 07/1962 - 12/1985
#             * '0': for the period 01/1986 - 03/2018
#             * '2': for the months 10/11/12 in 1987
#                     and the month 1 in 1988 (to remove the period about Black Monday);
# we want to distinguish one more case, namely if 'sdri' is missing;
# 'sdri' originally comes from the sp500 data and is not available as of
# 12/2004 (Bloom used the sp500 data only from 1962 - 2003);
sp500_merge_vxo <- sp500_merge_vxo %>%
                    mutate(miss = replace(miss, 
                        is.na(sdri), NA))

# if we run a filter to check, we'd see that the 0s turned to NAs for the
# variable 'miss'
sp500_merge_vxo %>%
                    filter((is.na(sdri)))

# the below command shows the distribution of the values of 'miss'
# across 0, 1 and NA.
# the value '1' corresponds to cases where we have both vox AND sp500
# data (i.e., the two series overlap!)
sp500_merge_vxo %>%
                    group_by(miss) %>%
                    summarise(Count = n())

# having calculated the moving average (i.e., the mean of the sdri's across
# year/month before, the below command calculates the standard deviation
# of these moving averages across all 'mvol' values for cases where 'mvol'
# is not missing (i.e., miss==0)

# the below code has been commented out because it follows Stata's procedure
# of how to calculate the desired results; in R itself, there is a simpler
# solution!


# # first we generate the variable 'sd_mvol';
# # but because sd_mvol now also contains the value for 'sd_mvol' for rows
# # where 'mvol' is actually missing, we perform another command (see below!)
# sp500_merge_vxo <- sp500_merge_vxo %>%
#                     mutate(sd_mvol = sd(mvol, na.rm=TRUE))
# # this command makes sure that only rows that actually have an 'mvol'
# # value also have the corresponding standard deviation in the column
# # 'sd_mvol'
# sp500_merge_vxo <- sp500_merge_vxo %>%
#                     mutate(sd_mvol = replace(sd_mvol,
#                         is.na(miss) | miss == 1, NA))
# 
# # to check that everything worked as expected, we can compare:
# sp500_merge_vxo %>%
#             group_by(sd_mvol) %>%
#             summarise(Count = n())
# # with
# sp500_merge_vxo %>%
#             group_by(miss) %>%
#             summarise(Count = n())
# 
# # only the rows from 01/1986 - until 12/2003 should now have a value populated
# # in 'sd_mvol':
# # we can check this by filtering as follows:
# sp500_merge_vxo %>%
#             filter(!(is.na(sd_mvol)))
# 
# # next, the below command calculates the standard deviation across all
# # sdri for cases where mvol is not missing (i.e., miss==0), following
# # the same procedure like above!
# sp500_merge_vxo <- sp500_merge_vxo %>%
#             mutate(sd_sdri = sd(sdri, na.rm=TRUE))
# # this command makes sure that only rows that actually have an 'mvol'
# # value also have the corresponding standard deviation in the column
# # 'sd_mvol'
# sp500_merge_vxo <- sp500_merge_vxo %>%
#             mutate(sd_sdri = replace(sd_sdri,
#                            is.na(miss) | miss == 1, NA))
# # to check that everything worked as expected, we can compare:
# sp500_merge_vxo %>%
#             group_by(sd_sdri) %>%
#             summarise(Count = n())
# # with
# sp500_merge_vxo %>%
#             group_by(miss) %>%
#             summarise(Count = n())
# 

# we want two numbers for further calclations:
#                         * sd_mvol = sd(mvol)
#                         * sd_sdri = sd(sdri)

# note that the 'pull' - command is necessary to transform the result
# of the summary to a numeric value!
# we want to calculate 'sd_mvol' and 'sd_sdri' for the period where the two
# series 'mvol' and 'sdri' overlap (so for the years 1987 - 2003)!
sd_mvol <- pull(summarize(sp500_merge_vxo, sd_mvol = sd(mvol[miss==0], na.rm = T)))
sd_sdri <- pull(summarize(sp500_merge_vxo, sd_sdri = sd(sdri[miss==0], na.rm = T)))

# in a later check of the code, we have to investigate whether the
# generation of the variable 'miss' and all the involved manipulations
# actually are necessary!

# next, we multiply 'sdri' by the ratio of the newly created values
# 'sd_mvol' and 'sd_sdri'
sp500_merge_vxo$sdri <- sp500_merge_vxo$sdri*(sd_mvol/sd_sdri)

# here we repeat the above procedure for the mean of both mvol and
# sdri
m_mvol <- pull(summarize(sp500_merge_vxo, m_mvol = mean(mvol, na.rm = T)))
m_sdri <- pull(summarize(sp500_merge_vxo, m_sdri = mean(sdri, na.rm = T)))

# now we replace sdri with sdri + the difference between m_mvol and m_sdri
sp500_merge_vxo$sdri <- sp500_merge_vxo$sdri + m_mvol - m_sdri

# because 'mvol' is missing between 1962 and 1985 (inclusive!), the
# below command fills up any missing 'mvol' values with
# 'sdri' in that period!

sp500_merge_vxo <- sp500_merge_vxo %>% 
                  mutate(mvol = ifelse(miss == 1 & year < 1987, sdri, mvol))

# the below command generates the variable 'my', which is actually the
# same like 'ym' with the difference that it now spans the entire dataset
# (ranging from 1962 to 2018);
# the variable 'ym' originally comes from the sp500 dataset and therefore
# only spans the range of the sp500 data!

sp500_merge_vxo <- as.data.frame(sp500_merge_vxo %>%
                    mutate(my = year + month/12))


# finally, we sort the data by the newly generated variable 'my'
sp500_merge_vxo <- sp500_merge_vxo[with(sp500_merge_vxo, order(my)), ]

# ultimately, for every month we only keep the top-most value!
sp500_merge_vxo <- sp500_merge_vxo %>%
                    group_by(my) %>% slice(1)


###############################
## (1.4) plotting of the data in sp500_merge_vxo
###############################

# in his original Stata-code, Bloom (2009) replaces mvol about 50 with 50 for
# scaling purposes; we decided to skip this step here!
# in our case instead we decided to take the full graph to give the
# reader the opportunity to see all volatilities in context!

# before we plot, note the following:
# the entire construction of the normalization and then
# filling up the missing 'mvol' values with values from 'sdri' means that 
# 'mvol' now contains the entire time series that we want to plot!

# to be able to distinguish those two parts of the series 'mvol', we 
# create the following indicator variable:

sp500_merge_vxo <- sp500_merge_vxo %>% 
                  mutate(mvol_ind = ifelse(year < 1986, "actual", "implied"))
# we can make use of the above created indicator variable in our plot below!

volatility <- ggplot(sp500_merge_vxo, 
                  aes(x = my, y = mvol, colour = mvol_ind)) +
  # geom_point() + 
  geom_line(size=0.8) + 
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Annualized standard deviation (%)", 
                     limits = c(0, 70), 
                     breaks = seq(0, 70, by = 10), 
                     minor_breaks = NULL) + 
  geom_vline(xintercept=1986, linetype="dashed", 
             color = "black", size=0.4) + 
  # theme_minimal() + 
  labs(col=NULL) + 
  theme(legend.position = c(0.935, 0.93), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))
  

volatility

# and we save the plot to be used in our latex-document
ggsave("volatility.pdf")



###############################
## (1.5) detrend the constructed time-series and generate
## Bloom-shocks
###############################
# first we make a time series out of sp500_merge_vxo$mvol;
# we do this using the funtionality from the 'zoo' - package;
# we discovered that this part is actually not necessary!
# Therefore, it is commented out for now!
mvol_ts <- zoo(sp500_merge_vxo$mvol, as.Date(sp500_merge_vxo$date))

# next we run the hp-filter
mvol_hp <- hpfilter(sp500_merge_vxo$mvol, type="lambda", freq=129600)

# next we can add the 'trend' and 'cycle' component of 'mvol_hp' to our
# existing data-frame:
sp500_merge_vxo$mvol_trend <- mvol_hp$trend
sp500_merge_vxo$mvol_cycle <- mvol_hp$cycle

# next we can plot everything we had plotted above
volatility_trend <- ggplot() +
  # geom_point() + 
  geom_line(data = sp500_merge_vxo, aes(x = my, y = mvol_trend), color="red", size=0.8) +
  geom_line(data = sp500_merge_vxo, aes(x = my, y = mvol), color="black", size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Annualized standard deviation (%)", 
                     limits = c(0, 70), 
                     breaks = seq(0, 70, by = 10), 
                     minor_breaks = NULL) + 
  geom_vline(xintercept=1986, linetype="dashed", 
             color = "black", size=0.4) + 
  # theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

volatility_trend 

# and we save the plot to be used in our latex-document
ggsave("volatility_trend.pdf")


volatility_cycle <- ggplot(sp500_merge_vxo, 
                           aes(x = my, y = mvol_cycle)) +
  # geom_point() + 
  geom_line(size=0.8) + 
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "annualized standard deviation (%), deviations from trend", 
                     limits = c(-20, 50), 
                     breaks = seq(-20, 50, by = 10), 
                     minor_breaks = NULL) + 
  # geom_hline(aes(yintercept = 0)) +
  # theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

volatility_cycle

ggsave("volatility_cycle.pdf")


# According to Bloom (2009): The shocks are chosen as those events with stock-market
# volatility more than 1.65 standard deviations above the HP-detrended mean of the
# mean of the stock-market volatility series.
# Further: Each month has been treated as an independent observation.
# For us this means that we can calculate the sample standard deviation 
# as follows:

sd_mvol_cycle <- sd(mvol_hp$cycle)
mean_mvol_cycle <- mean(mvol_hp$cycle)

# we generate a plot to show the isolated shocks according to 
# Bloom's methodology (those events where the red-dashed line
# crosses the volatility-series)!
volatility_cycle_shocks <- ggplot(sp500_merge_vxo, 
                           aes(x = my, y = mvol_cycle)) +
  # geom_point() + 
  geom_line(size=0.8) + 
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "annualized standard deviation (%), deviations from trend", 
                     limits = c(-20, 50), 
                     breaks = seq(-20, 50, by = 10), 
                     minor_breaks = NULL) + 
  # geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=mean_mvol_cycle+sd_mvol_cycle*1.65, linetype="dashed", 
             color = "red", size=0.4) +
  # theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

volatility_cycle_shocks

ggsave("volatility_cycle_shocks.pdf")

# the red dotted line in our plot above marks the following threshold:
# yintercept=mean_mvol_cycle+sd_mvol_cycle*1.65, i.e., it is exactly what
# Bloom (2009) refers to as the threshold which was selected as the
# 5% one-tailed significance level TREATING EACH MONTH AS AN INDEPENDENT
# OBSERVATION!

# we store the threshold value into our dataset to create an indicator
# variable that allows us to extract the 'Bloom' - shocks which we
# want to have a look at
sp500_merge_vxo$thresh <- mean_mvol_cycle+sd_mvol_cycle*1.65

# the above variable now allows us to create an indicator of whether
# or not a certain point is above the threshold value
sp500_merge_vxo <- sp500_merge_vxo %>%
                      mutate(bloom_shock = case_when(mvol_cycle > thresh ~ 1,
                                            mvol_cycle <= thresh ~ 0))


###############################
## (1.6) replicate plot from (4) and add
## periods of 'Bloom-shocks'
###############################

# having the variable bloom_shock in our dataset now allows us to replicate
# our above graph from (4) and add episodes into the graph that correspond to the
# Bloom-shock being 1

# but before we can do that, we need to create 'start' and 'end' - dates for
# the respective periods (which we will then ultimately use for plotting
# the episodes with shaded regions!)

# to retrieve the respective 'start'- and 'end' - dates for the shock-periods,
# we apply a procedure to the dataset 'sp500_merge_vxo' that results in a 
# separate dedicated data-frame with the following variables:
# 'start', 'end', 'max_volatility', 'first_volatility';
# once we have the above four variables, we can easily calculate 'duration_months'
# by subtracting 'start' from 'end'.

# the construction of such a table will allow us
# (1) to automatically generate the respective data even if something in our
#     dataset changes at some point (i.e., it is not hard-coded)
# (2) we can automatically export the information to a latex-table
# to be included in our text by means of e.g., the R-package
# 'miscFuncs' that is capable of exporting matrix-like data
# to latex-tables!

# let us first filter for all rows where the bloom_shock = 1 to get an overview:
as.data.frame(sp500_merge_vxo %>% filter(bloom_shock == 1))

# initially, we had thought of creating a separate data-frame
# and store the loop's data into it;
# but then we decided for a different approach:
# we want to give each episode of a shock a unique ID so that
# we then can easily calculate the max volatility and first
# volatility by grouping by ID (within the full sp500_merge_vxo table!)



# we decided for the following procedure:
# we loop through all rows, and first check if the bloom_shock equals 1
# or not:

# we initialize the grouping variable x as follows:
x <- 2

# we add an empty column to our data frame sp500_merge_vxo:
sp500_merge_vxo["shock_ID"] <- NA

# then we start the loop:
for (i in 1:nrow(sp500_merge_vxo)) {
  if(sp500_merge_vxo$bloom_shock[i] == 1) {
    # if we detect a shock (i.e., 'bloom_shock' == 1), then we have 
    # to proceed as follows:
    # we store the 'my' variable in the column 'start' and 'end'
    sp500_merge_vxo[i, ncol(sp500_merge_vxo)] <- x
    
  }else{
    # no shock
    x <- x + 1
  }
}

# the above method is maybe not yet the most elegant/efficient solution
# but it suffices for now (we can fine-tune it at a later stage!)
# we can now have a look at the newly created variable
# (note that we have suppressed 'NA's for the creation of the table!):

sp500_merge_vxo %>%
          filter(!is.na(shock_ID))  %>%
          group_by(shock_ID) %>%
          summarise(Count = n())

# next we group_by shock_ID and extract the respective start-dates and
# store the retrieved data to a new data-frame
shocks_start <- as.data.frame(sp500_merge_vxo %>%
                                    filter(!is.na(shock_ID))  %>%
                                    group_by(shock_ID) %>%
                                    filter(row_number()==1) %>%
                                    mutate(year_start=year, month_start=month, my_start = my) %>%
                                    dplyr::select(shock_ID, year_start, month_start, my_start))

# next, we replicate the above query for the end-dates
# (note that in some scenarios start- and end-dates are identical!)
shocks_end <- as.data.frame(sp500_merge_vxo %>%
                                    filter(!is.na(shock_ID))  %>%
                                    group_by(shock_ID) %>%
                                    filter(row_number()==n()) %>%
                                    mutate(year_end=year, month_end=month, my_end = my) %>%
                                    dplyr::select(shock_ID, year_end, month_end, my_end))
                                    
# next, we can merge the two data-frames from above:
shocks_start_end <- merge(x = shocks_start, y = shocks_end, by = "shock_ID", all = TRUE, na.rm=T)
# and inspect the resulting data-frame:
shocks_start_end

# we re-arrange the sequence of columns
shocks_start_end <- shocks_start_end %>%
                  dplyr::select(shock_ID, year_start, month_start, year_end, month_end, my_start, my_end)


# next, we construct a year-month-variable both out of
# the pair 'year_start' & 'month_start' and 'year_end' & 'month_end'
shocks_start_end$yearmon_start <- as.yearmon(shocks_start_end$my_start, "%Y-%B")
shocks_start_end$yearmon_end <- as.yearmon(shocks_start_end$my_end, "%Y-%B")
# next, we add a helper-column that we need in the next stage as well:
shocks_start_end$helper_date <- format(shocks_start_end$yearmon_start, "%b")

# this means that we can now drop the 'year' and 'month' variables
shocks_start_end <- shocks_start_end %>%
                  dplyr::select(-c(year_start, year_end, month_start, month_end))

# next we add a column that gives us the duration of the shock:
shocks_start_end <- shocks_start_end %>%
                  mutate(duration = (shocks_start_end$yearmon_end - 
                                       shocks_start_end$yearmon_start) * 12 + 1) %>%
                  dplyr::select(-helper_date)

# next we group_by shock_ID and extract the respective MAXIMUM
# VOLATILITY
# (note that the produced data-frame only has one row per shock_ID)
shocks_max_vol <- as.data.frame(sp500_merge_vxo %>%
            filter(!is.na(shock_ID))  %>%
            group_by(shock_ID) %>%
            filter(mvol == max(mvol)) %>%
            mutate(max_vol = mvol) %>%
            dplyr::select(shock_ID, max_vol, my))

# we apply the same date-transformation as above:
shocks_max_vol$yearmon_max <- as.yearmon(shocks_max_vol$my, "%Y-%B")
# and we drop 'my':
shocks_max_vol <- shocks_max_vol %>%
                dplyr::select(-my)

# in a last stage, we merge shocks_max_vol with shocks_start_end:
shocks_start_end <- merge(x = shocks_start_end, y = shocks_max_vol, 
                          by = "shock_ID", all = TRUE, na.rm=T)

# the above data-frame now allows us to plot the episodes of 
# high volatility into our previous plots (by adding a shaded
# rectangle!)

# actually, there is a slight problem, which we fix by adding
# one month's numeric value to my_end:
shocks_start_end$my_end <- shocks_start_end$my_end + (1/12)
# this makes sure that if we refer to an end-month that
# we assume that we are talking about the last day of that
# respective month!

volatility_cycle_shocks2 <- ggplot(sp500_merge_vxo, 
                                  aes(x = my, y = mvol_cycle)) +
  # geom_point() + 
  geom_line(size=0.8) + 
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "annualized standard deviation (%), deviations from trend", 
                     limits = c(-20, 50), 
                     breaks = seq(-20, 50, by = 10), 
                     minor_breaks = NULL) + 
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=mean_mvol_cycle+sd_mvol_cycle*1.65, linetype="dashed", 
             color = "red", size=0.4) +
  #theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14), panel.grid.major.x = element_blank()) + 
  # note that panel.grid.major.x = element_blank() suppresses vertical grid lines!
  geom_rect(data=shocks_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='red', alpha=0.5)

volatility_cycle_shocks2

ggsave("volatility_cycle_shocks2.pdf")


# having created the data-frame shocks_start_end, we can now also
# export a latex-table with the exact data of the Bloom-shock!
# first, we slightly adapt the shocks - data frame:
# first, we add '%' to 'max_vol' in our data-frame:
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x, format = format, digits = digits, ...), "%")
}

shocks_table <- shocks_start_end %>%
              mutate(max_vol = percent(max_vol), duration = paste(as.character(floor(duration)), "months", sep=" ")) %>%
              dplyr::select(yearmon_start, yearmon_end, duration, yearmon_max, max_vol)
                
# adapt the below command to supply a vector of column-headers!
latextable(shocks_table, dp=3)



########################################################################
### PART 2: construction of impulse-responses following Jordá (2005)
########################################################################

## PART 2 is separated into ....... parts:
## (2.1)  2.1 handles ........
## (1.2)  1.2 handles the .......


###############################
## (2.1) loading additional_vars.csv;
##       merge with variables from sp500_merge_vxo
##       (most of the data comes from FRED 
##       [https://fred.stlouisfed.org/series/])
##       the employment in manufacturing comes
##       from BLS (see https://data.bls.gov/pdq/SurveyOutputServlet)
###############################

# first we read in the file 'additional_vars.csv'
# additional_vars <- read.csv(file="additional_vars.csv", header=TRUE, sep=",",
#                            stringsAsFactors = FALSE)
additional_vars <- read.csv(file="additional_vars_ext.csv", header=TRUE, sep=",",
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

# and drop a few superfuous columns:
additional_vars <- as.data.frame(additional_vars %>%
                    dplyr::select(-c(DATE.1, Date, Open, High, Low, Adj.Close, Volume)))


# the variable EMPM has to be retrieved from a separate dataset
# which we downloaded from the homepage of the BLS;
# we have to manipulate the data to bring it into our desired format:
EMPM <- read.xlsx("BLS_EMPM.xlsx", header = TRUE, sheetIndex = 1)
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
# while the only goal was to consturct a consistent volatility - measure
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
# and we transform the variables which we want to enter the regressions
# as logs to logs:
regressions <- as.data.frame(regressions %>%
              mutate(lip = log(ip),
                     lsp500 = log(sp500),
                     lempm = log(empm)))

# next, following Bloom (2009) we detrend the variables by means
# of the hp-filter:
lip_hp <- hpfilter(regressions$lip, freq=129600)
lsp500_hp <- hpfilter(regressions$lsp500, freq=129600)
ffr_hp <- hpfilter(regressions$ffr, freq=129600)
lempm_hp <- hpfilter(regressions$lempm, freq=129600)

# next we can add the 'cycle' component of the detrended variables to our
# existing data-frame:
# regressions$lip_trend <- lip_hp$trend
regressions$lip_cycle <- lip_hp$cycle
regressions$lsp500_cycle <- lsp500_hp$cycle
regressions$ffr_cycle <- ffr_hp$cycle
regressions$lempm_cycle <- lempm_hp$cycle

# in particular, the variable 'h' will be used in our below loop to know
# exactly where (i.e., in which row) we want to place the coefficients from
# the regressions that we run!

# our goal now is to loop i from 0 to 48 (i.e., the horizons that we 
# want to consider);
# for the regressinos we use 'dynlm' to have easy access to LAG-operators
# within the regressions!
# before we run the regressions, we declare all variables to be 
# of the class 'time - series'
regressions$bloom_shock <- as.ts(regressions$bloom_shock)
regressions$lip_cycle <- as.ts(regressions$lip_cycle)
regressions$lsp500_cycle <- as.ts(regressions$lsp500_cycle)
regressions$ffr_cycle <- as.ts(regressions$ffr_cycle)
regressions$lip <- as.ts(regressions$lip)
regressions$lempm <- as.ts(regressions$lempm)
regressions$lempm_cycle <- as.ts(regressions$lempm_cycle)
# we create emtpy placeholders for the coefficients and
# upper and lower bounds of the regressions in a 
# new data-frame which we call 'regressions_out':
regressions_out <- data.frame(
              b_lip=numeric(),
              lo90_b_lip=numeric(),
              up90_b_lip=numeric(),
              stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  # i <- 2
  print(i)
  
  reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + lsp500_cycle + ffr_cycle +
                 + L(lip, seq.int(1:3)*(-1)), 
                            regressions)
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]))
  regressions_out <- rbind(regressions_out, estimates)
  rm(reg, SE_robust, estimates)
}

# rename the column names
colnames(regressions_out)[1] <- "b_lip"
colnames(regressions_out)[2] <- "lo90_b_lip"
colnames(regressions_out)[3] <- "up90_b_lip"

# and again add a column 'h'
regressions_out <- as.data.frame(regressions_out %>%
                               mutate(h = seq(0, nrow(regressions_out)-1)))

# next we can plot everything we had plotted above
irf1 <- ggplot() +
  # geom_point() + 
  geom_line(data = regressions_out, aes(x = h, 
                      y = 100*(exp(up90_b_lip)-1)), 
                      color="#677066", 
                      size=0.5,linetype = 4) +
  geom_line(data = regressions_out, aes(x = h, 
                      y = 100*(exp(lo90_b_lip)-1)), 
                      color="#677066", size=0.5, 
                        linetype = 4) +
  geom_ribbon(data = regressions_out, aes(x=h, ymax=100*(exp(up90_b_lip)-1), ymin=100*(exp(lo90_b_lip)-1)), fill="pink", alpha=.3) +
  geom_line(data = regressions_out, aes(x = h, y = 100*(exp(b_lip)-1)), color="black", size=0.8) +
  geom_point(data = regressions_out, aes(x = h, y = 100*(exp(b_lip)-1)), color="black", size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = "months", limits = c(0, 60), 
                     breaks = seq(0, 60, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "% impact on industrial production", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  # theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

irf1

# and we save the plot to be used in our latex-document
ggsave("irf1.pdf")


## we generate another irf (this time for employment in manufacturing)

regressions_out <- data.frame(
  b_lip=numeric(),
  lo90_b_lip=numeric(),
  up90_b_lip=numeric(),
  stringsAsFactors=FALSE
)


for (i in 0:60){
  
  # i stands for the horizons that we consider
  #i <- 2
  print(i)
  
  reg <- dynlm(L(lempm_cycle, -i) ~ bloom_shock + lsp500_cycle + ffr_cycle +
                 + L(lempm, seq.int(1:3)*(-1)), 
               regressions)
  #reg <- dynlm(L(lip_cycle, -i) ~ bloom_shock + L(lip, seq.int(1:3)*(-1)), regressions)
  # transformation to newey-west standard errors
  SE_robust <- sqrt(diag(vcovHAC(reg)))
  # good explanation between the vairous 'flavours' of options that the
  # 'sandwich' - package offers:
  # https://stats.stackexchange.com/questions/15608/vcovhc-vcovhac-neweywest-which-function-to-use
  # note: there are slight discrepancies between the above and the
  # below estimation of the standard errors which I still have to
  # figure out!
  # reg_helper <- coeftest(reg, vcov = vcovHAC(reg))
  
  # having estimated the model with HAC-consistent standard errors, we can
  # store the coefficients and upper and lower bounds of the corresponding
  # confidence bands for each horizon:
  estimates <- c(summary(reg)$coefficients[2, 1], 
                 as.numeric(summary(reg)$coefficients[2, 1]+ 1.68*SE_robust[2]),
                 as.numeric(summary(reg)$coefficients[2, 1] - 1.68*SE_robust[2]))
  regressions_out <- rbind(regressions_out, estimates)
  rm(reg, SE_robust, estimates)
}

# rename the column names
colnames(regressions_out)[1] <- "b_lip"
colnames(regressions_out)[2] <- "lo90_b_lip"
colnames(regressions_out)[3] <- "up90_b_lip"

# and again add a column 'h'
regressions_out <- as.data.frame(regressions_out %>%
                                   mutate(h = seq(0, nrow(regressions_out)-1)))

# next we can plot everything we had plotted above
irf2 <- ggplot() +
  # geom_point() + 
  geom_line(data = regressions_out, aes(x = h, y = 100*(exp(up90_b_lip)-1)), color="#677066", size=0.5,linetype = 4) +
  geom_line(data = regressions_out, aes(x = h, y = 100*(exp(lo90_b_lip)-1)), color="#677066", size=0.5, linetype = 4) +
  geom_ribbon(data = regressions_out, aes(x=h, ymax=100*(exp(up90_b_lip)-1), ymin=100*(exp(lo90_b_lip)-1)), fill="pink", alpha=.3) +
  geom_line(data = regressions_out, aes(x = h, y = 100*(exp(b_lip)-1)), color="black", size=0.8) +
  geom_point(data = regressions_out, aes(x = h, y = 100*(exp(b_lip)-1)), color="black", size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = "months", limits = c(0, 60), 
                     breaks = seq(0, 60, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "% impact on employment in manufacturing", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  # theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

irf2

# and we save the plot to be used in our latex-document
ggsave("irf2.pdf")





