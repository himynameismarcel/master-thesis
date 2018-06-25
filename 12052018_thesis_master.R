#########################################################################################
### Marcel Kropp, 15.04.2018
### This script is separated into several big parts:
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
### PART 2: deals with the EPU-index constructed by Baker, Bloom and Davis (2016)
###         Baker, Bloom and Davis (2016) have dedicated an entire web-page
###         to their constructed uncertainty measure where they not
###         only provide several indices for the US, but also several
###         other countries, categories, etc.
###         See the home-page at http://www.policyuncertainty.com/us_monthly.html
###         Because we deal with the US only, we have correspondingly only
###         downloaded data for the US!
### Part 3: deals with the index as constructed by Leduc and Liu (2016)
###         based on the Michigan Survey.
###         The data off of which Leduc and Lui construct their index
###         is available at https://data.sca.isr.umich.edu/tables.php;
###         In particular, the data derived from the questions 
###         'Buying Conditions for Vehicles' and the follow-up
###         questions 'Reasons for Opinions for Buying Conditions for Vehicles'
###         are used.
### Part 4: deals with the GT (Bontempi et al, 2015) and GTU-index
###         (Castelnuovo and Tran, 2017);
###         Both contributions use Google-Trends data to contruct
###         an uncertainty index (with subtle differences in the words
###         they search for; details see below);
###         both indices cover the period 2004 - 2018 (to date),
###         the uncertainty index from Bontempi et al. (2015) 
###         is available from their dedicated website at
###         https://sites.google.com/site/uncertaintyindexbasedongt/;
###         the uncertainty index from Castelnuovo and Tran (2017)
###         is available on Efrem Castelnuovo's homepage at
###         https://sites.google.com/site/efremcastelnuovo/home/publications;
### PART 5: deals with the Macro Uncertainty Index constructed
###         by Jurado et al. (2015)
###         The scripts accompanying their paper is available on their
###         paper's dedicated web presence of the American Economic Review
###         (see https://www.aeaweb.org/articles?id=10.1257/aer.20131193)
###         The constructed macro and firm-uncertainty index is available
###         at Sydney Ludvigson's homepage at
###         https://www.sydneyludvigson.com/data-and-appendixes/;
###         There both the original data used for the paper with the data
###         ranging from 1960 until 2015 as well as updated versions
###         ranging from 1960 until 2017.
### PART 6: deals with the plotting of all our analyzed time-series into
###         one plot (to make them comparable!)
###         the series which we want to plot together at once are
###           *) VIX (used by Bloom)
###           *) Michigan Survey
###           *) EPU (Baker et al, 2016)
###           *) GT and GTU (Bontempi et al. 2015/Castelnuovo and Tran, 2017)
###           *) Macro Uncertainty Index (Jurado et al, 2015)
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
## (1.2)  1.2 handles the raw VXO volatility figures downloaded directly from
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
# install.packages("plyr")
# install.packages("dplyr")
# install.packages("tibble")
# install.packages("tidyr")
# install.packages("lubridate")
# install.packages("ggplot2")
# install.packages("mFilter")
# install.packages("zoo")
# install.packages("miscFuncs")
# install.packages("dynlm")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("xlsx")
# install.packages("readxl")
# install.packages("forecast")
# install.packages("stargazer")
# install.packages("vars")
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("anomalize")
# install.packages("rowr")
# install.packages("xts")
# install.packages("dse")
# install.packages("e1071")
# install.packages("Hmisc")
# install.packages("broom")
# install.packages("purrr")
# install.packages("timetk")
# install.packages("sweep")
# install.packages("tseries")
# install.packages("highfrequency")
# install.packages("data.table")
# install.packages("matlib")

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
library(readxl)
library(forecast)
library(stargazer)
library(vars)
library(grid)
library(gridExtra)
library(ggpubr)
library(anomalize)
library(rowr)
library(xts)
library(dse)
library(e1071)
library(Hmisc)
library(broom)
library(purrr)
library(timetk)
library(sweep)
library(tseries)
library(highfrequency)
library(gtable)
library(data.table)
library(matlib)

# clear workspace
rm(list = ls())

###############################
## (1.1) s&p500 historical data
###############################
# Generate baseline stock volatility data (to merge with
# VXO data)
# based on sp500.csv (data ranges from 1962 - 2003, daily
# data!)

# first we read in the historical sp500 data (for which Bloom
# does not exactly declare the source)
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

# next, we calculate the standard deviation for all 'ri' (daily returns
# of the index in %) for each year-month (i.e., BY the variable 'ym'!);
# the most convenient function to do this is the 'ddply' function (see 
# the reference
# here: https://www.rdocumentation.org/packages/plyr/versions/1.8.4/
# topics/ddply);
# I found this solution suggested here: 
# https://stackoverflow.com/questions/15467219/calculate-group-
# characteristics-without-ddply-and-merge#comment21890299_15467219
# For a detailed treatment of the 'plyr' - package, see Hadley's article: 
# file:///Users/marcelkropp/Downloads/v40i01.pdf.
# Also, the following blog-post gives an overview: http://stat545.com/
# block013_plyr-ddply.html.
# In essence: 'ddply()' accepts a data.frame, splits it into pieces 
# based on one or more
# factors, computes on the pieces, then returns the results as a 
# data.frame
# (the base R functions most relevant to 'ddply()' are 'tapply()', 
# etc.)

sp500 <- ddply(sp500, ~ym, mutate, 
               sdri = sd(ri, na.rm=T))
# Note that we have added 'na.rm=T' to the 'sd()' - function to make 
# sure that also for 196207 the standard deviation is calculated which 
# has an NA for the the very first value (otherwise all entries for the 
# month 07 would end up being NA).

# finally, we order the data frame (actually: tibble) by the variable 
# 'date' (which in our case is a peculiar variable but which fulfills 
# the same purpose like the actual 'date' - variable 'caldt')
sp500 <- as_data_frame(sp500[with(sp500, order(date)), ])

###############################
## (1.2) VXO - data from CBOE
###############################
# Here we process the raw volatility figures which we have 
# downloaded directly
# from the CBOE's homepage and which we have stored in 'vxo_new.csv'
# (data ranges from 1986 to 2018, daily data!)

# first we read in the VOX-data
vxo_new <- read.csv(file="vxo_new.csv", header=TRUE, sep=",")
# vxo_new <- read.csv(file="vxo.csv", header=TRUE, sep=",")

# we check out the structure of 'vxo_new'
str(vxo_new)
# as we can see, the variable 'vol' has already been read-in as 
# 'numeric', so everything is okay

# next, we want to split the variable 'date' into three variables
# 'month', 'day' and 'year';
# for this, we first need to convert the 'date' - variable 
# (which is currently declared as a factor) into a DATE datatype:
vxo_new$date <- as.Date(vxo_new$date, format = "%m/%d/%Y")
# we change the format of the date variable
vxo_new$date <- format(vxo_new$date, "%m/%d/%Y")
# note that the above command converts the date-type to a character 
# again (we will deal with this later!)

# next, we create the three variable 'month', 'year' and 'day' 
# using the 'separate()' - function from the 'tidyr' - package
vxo_new <- separate(vxo_new, "date", c("month", "day", "year"), 
                    sep = "/", 
                    remove=FALSE, convert=TRUE)
# as we can see, with the option 'convert=TRUE', the values are 
# automatically transformed into integers!

# having sliced our date variable into month, day and year,
# we can drop the original variable 'date' (actually, we could
# have dropped 'date' automatically by keeping 'remove=TRUE', but
# we preferred it this way to double-check the output!)

# we convert our date-variable back to a date datatype
vxo_new$date <- as.Date(vxo_new$date, format = "%m/%d/%Y")

# next, we want to apply a consistency check to make sure that 
# there are no trading days
# on Sat or Sun in our data!
# for this, we create two helper variables: 'm2' = 'month' and 
# 'd2'='day', date2=date
vxo_new <- mutate(vxo_new, m2 = month, d2=day, date2=date)
# and then extract the day of the week using the 'wday' component 
# of a 'POSIXlt' - object (because it is numeric, starting on 
# Sunday) as opposed to the normal 'weekdays' - function
vxo_new$dow <- as.POSIXlt(vxo_new$date)$wday
# Note: also the package lubridate should have a 'wday' 
# componsent (check out later!)

# next, Bloom (2009) performs two very strange commands that 
# actually do not make much sense (have to ask Martin about this!)

# Above, we have created the variable 'date2' as a copy of 'date';
# In Bloom's stata-code, the created variable 'date' (which was 
# made up of month, day and year) is a date in elapsed date 
# format.
vxo_new$date2 <- as.numeric(vxo_new$date2)
# Note: for what follows next, we cannot use the elapsed time 
# values that Bloom uses, because 

# we drop the 'day' - variable (which we had actually generated 
# above; but because we use a different strategy for the generation 
# of dates than Bloom, we could have also easily just skipped it!);
# edit: we decided to keep the 'day' - variable to have a better 
# insight into whther or not our 'merge' of the two datasets below 
# worked or not!
# vxo_new <- vxo_new %>% select(-day)

# next, we fill in 9/11 figures using European exchange extrapolation
# the following filter-statement shows us where we would need to input
# the data:
vxo_new %>% filter(date2>11574 & date2 < 11590)

# it is exactly between 2001-09-10 where for six days until 
# 2001-09-17 there is no data available;
# our goal is to add additional rows to our data-frame, remove 
# weekends and then overwrite the vol-value with the value from 
# the European exchange!

# first we generate an additional variable that gives us the 
# numbers of the row that we are currently in in our data frame:
vxo_new$row <- seq.int(nrow(vxo_new))

# running our above filter-statement again, we know in which rows 
# we want to replicate:
vxo_new %>% filter(date2>11574 & date2 < 11590)

# we replicate row number 3958 that holds the data for 2001-09-10
vxo_new <- rbind(vxo_new, vxo_new[rep(3958, 6), ])
# so in our case, this command creates 9 additional duplicates of 
# 2001-09-10

# we order by date
vxo_new <- vxo_new[with(vxo_new, order(date)), ]

# filtering with the below command would show us that everything
# worked as expected
vxo_new %>% filter(date2>11574 & date2 < 11590)

# we need to replace the dates between 2001-09-10 and 2001-09-17
# with their actual values!
vxo_new$date[vxo_new$row==3958] <- seq(as.Date("2001-09-10"), 
                                       as.Date("2001-09-16"), 
                                       by="days")

# we do the same with the the date2-column
vxo_new$date2[vxo_new$row==3958] <- seq(11575, 11581, 1)
# and then re-create the row-column
vxo_new$row <- seq.int(nrow(vxo_new))

# then we replace the value for vol from 2001-09-11 until 
# 2001-09-16 (inclusive) with the value from the European 
# exchange
vxo_new$vol[vxo_new$date2>11575 & vxo_new$date2<11582] <- 58.2

# and re-create the dow-variable
vxo_new$dow <- as.POSIXlt(vxo_new$date)$wday
# to be able to drop weekend days
vxo_new <- vxo_new %>%
            filter(dow != 6 & dow != 0)

# and we sort again by date
vxo_new <- vxo_new[with(vxo_new, order(date)), ]

# next, we want to drop months with less than 12 days of data:
# for this purpose we first create a variable 'count' that counts 
# the number
# of days per month
vxo_new <- vxo_new %>% 
                        group_by(year, month) %>% 
                        mutate(dayc = n())
# we convert it back to a data frame
vxo_new <- as.data.frame(vxo_new)

# in our cae, the below command causes the observations for 
# April 2018 to be deleted because we only have 8 observations 
# when we downloaded the data;
# therefore, depending on when we run the final code, it might be 
# that the month that we are currently in will still be included!
vxo_new <- vxo_new %>%
            filter(dayc>10)

# next, we drop these two helper variables (dayc, date2) because 
# we do not need them anymore in what follows:
# (note that date2 was a helper variable which we had used in 
# the expansion of our dataset when we wanted to include the data 
# from the European exchanges);
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
# https://stackoverflow.com/questions/1299871/how-to-join-merge-
# data-frames-inner-outer-left-right);

# we decided for the dplyr-way:
# before applying the dplyr-way of a join, we have to make sure that the
# columns by which we want to join are both type numeric
vxo_new$date <- as.numeric(vxo_new$date)

# we add an extension: 
# because we cannot interpret the meaning of 'date' - variable in the 
# sp500 dataset, we copy the variable, transfer it to a date-datatype 
# and then to a numeric date datatype (which counts the elapsed 
# numer of days since 1970!)
# first we copy date to date2 (to free the name 'date' for the actual 
# variable that we want to create!)
sp500$date2 <- sp500$date
# then we overwrite date with the values from caldt (which are still 
# integers!)
sp500$date <- sp500$caldt
# we make a 'date' - datatype out of 'date'
sp500$date <- as.Date(as.character(sp500$date), "%Y%m%d")
# and then we make a numeric date-variable out of it
sp500$date <- as.numeric(sp500$date)

# the below command finally performs an outer join of the two datasets
sp500_merge_vxo <- merge(x = vxo_new, y = sp500, 
                         by = "date", all = TRUE, 
                         na.rm=T)
# an inspection of the new data-frame sp500_merge_vxo reveals:
str(sp500_merge_vxo)
# we have year.y and month.y which are the year and month variables 
# from the sp500 dataset,
# we have year.x and month.x which are the year and month variables 
# from the vxo_new dataset;
# to explain the many NAs (be aware that we performed an outer join),
# we have to keep in mind that the sp500 dataset ranges from 07/1962 
#                         until 12/2003 and
#                         that the vxo_new data ranges from 01/1986 
#                         until 04/2018;
# this means that in the period 07/1962 until 12/1985 we only 
#                         have sp500 data,
#                 in the period 01/1986 until 12/2003 we have BOTH 
#                         sp500 and vxo_new data,
#             and in the period 01/2004 until 04/2018 weonly have 
#                         vxo_new data!

# to execute the below commands, we need a 'year' - and 'month' - 
# variable that spans the entire sp500_merge_vxo dataset!
sp500_merge_vxo <- sp500_merge_vxo %>%
                  mutate(year = case_when(is.na(year.x) ~ year.y,
                          !(is.na(year.x)) ~ year.x))
sp500_merge_vxo <- sp500_merge_vxo %>%
                  mutate(month = case_when(is.na(month.x) ~ month.y,
                           !(is.na(month.x)) ~ month.x))

# after the above merge, the below command drops only 7 observations 
# after 1990 that are only in the sp500 dataset (i.e., have no match 
# with vox_new.csv after 1990!)
# exactly 7 such cases exist:
sp500_merge_vxo %>%
                  filter(!complete.cases(vol) & year > 1990)

sp500_merge_vxo <- sp500_merge_vxo %>%
                  filter(!(!complete.cases(vol) & year > 1990))

# the next command creates a moving average of the vxo-volatility value.
# this calculation is necessary because currently the vol-values from 
# the vxo-dataset are DAILY figures (whereas the sdri-values are 
# the standard deviation of the daily returns over a month!)
sp500_merge_vxo <- as.data.frame((sp500_merge_vxo %>% 
                    group_by(year.x, month.x) %>% 
                    mutate(mvol = mean(vol))))

# we check whether indeed the calculation worked:
sp500_merge_vxo %>%
                  filter(complete.cases(vol))

# for the below command note that sdri = standard deviation of 
# the daily returns over a month/year (as stated above!)
# the below command annualizes our volatility (standard deviation)
# of the monthly figures (ask Martin why Bloom used the square
# root of 366???)
sp500_merge_vxo$sdri <- sp500_merge_vxo$sdri * sqrt(366)

# the entire below sequence normalizes the two series (of vol and sdri)
# to have the same mean and variance (excluding the Black-Monday values 
# along the way because they are so extreme!)

# to fully understand the sequence of commands below, we have to remind
# ourselves again of the following:
# we have to keep in mind that the sp500 dataset ranges from 07/1962 
#                         until 12/2003 and
#                         that the vxo_new data ranges from 01/1986 
#                         until 04/2018;
# this means that in the period 07/1962 until 12/1985 we only have 
#                         sp500 data,
#                 in the period 01/1986 until 12/2003 we have BOTH 
#                         sp500 and vxo_new data,
#             and in the period 01/2004 until 04/2018 weonly have 
#                         vxo_new data!

# we create the variable 'miss' that = '1' for cases where 'mvol' 
#                         is missing (i.e., NA)
#                         and = '0' for cases where mvol is 
#                         not missing;
# (remember that 'mvol' is the moving average of the 'vxo' 
# volatility value by year and month!)

sp500_merge_vxo <- sp500_merge_vxo %>%
                    mutate(miss = case_when(
                          is.na(mvol) ~ 1,
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
                    filter((year.y==1987 & month.y 
                            %in% c(10, 11, 12)) |
                     (year.y==1988 & month.y==1)))


# currently the variable 'miss' only consists of 0s and 1s as follows:
#             * '1': for the period 07/1962 - 12/1985
#             * '0': for the period 01/1986 - 03/2018
#             * '2': for the months 10/11/12 in 1987
#                     and the month 1 in 1988 (to remove the period 
#                     about Black Monday);
# we want to distinguish one more case, namely if 'sdri' is missing;
# 'sdri' originally comes from the sp500 data and is not available as 
# of 12/2004 (Bloom used the sp500 data only from 1962 - 2003);
sp500_merge_vxo <- sp500_merge_vxo %>%
                    mutate(miss = replace(miss, 
                        is.na(sdri), NA))

# if we run a filter to check, we'd see that the 0s turned to NAs for the
# variable 'miss'
sp500_merge_vxo %>%
                    filter((is.na(sdri)))

# the below command shows the distribution of the values of 'miss'
# across 0, 1 and NA.
# the value '1' corresponds to cases where we have both VXO AND sp500
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
# we want to calculate 'sd_mvol' and 'sd_sdri' for the period where 
# the two series 'mvol' and 'sdri' overlap (so for the years 1987 - 
# 2003)!
sd_mvol <- pull(dplyr::summarize(sp500_merge_vxo, sd_mvol = 
                            sd(mvol[miss==0], na.rm = T)))
sd_sdri <- pull(dplyr::summarize(sp500_merge_vxo, sd_sdri = 
                            sd(sdri[miss==0], na.rm = T)))

# in a later check of the code, we have to investigate whether the
# generation of the variable 'miss' and all the involved 
# manipulations actually are necessary!

# next, we multiply 'sdri' by the ratio of the newly created values
# 'sd_mvol' and 'sd_sdri'
sp500_merge_vxo$sdri <- sp500_merge_vxo$sdri*(sd_mvol/sd_sdri)

# here we repeat the above procedure for the mean of both mvol and
# sdri
m_mvol <- pull(dplyr::summarize(sp500_merge_vxo, m_mvol = 
                           mean(mvol, na.rm = T)))
m_sdri <- pull(dplyr::summarize(sp500_merge_vxo, m_sdri = 
                           mean(sdri, na.rm = T)))

# now we replace sdri with sdri + the difference between 
# m_mvol and m_sdri
sp500_merge_vxo$sdri <- sp500_merge_vxo$sdri + m_mvol - m_sdri

# because 'mvol' is missing between 1962 and 1985 (inclusive!), the
# below command fills up any missing 'mvol' values with
# 'sdri' in that period!

sp500_merge_vxo <- sp500_merge_vxo %>% 
                  mutate(mvol = ifelse(miss == 1 & 
                                         year < 1987, sdri, mvol))

# the below command generates the variable 'my', which is actually 
# the same like 'ym' with the difference that it now spans the entire 
# dataset (ranging from 1962 to 2018);
# the variable 'ym' originally comes from the sp500 dataset and 
# therefore only spans the range of the sp500 data!

sp500_merge_vxo <- as.data.frame(sp500_merge_vxo %>%
                    mutate(my = year + month/12))


# finally, we sort the data by the newly generated variable 'my'
sp500_merge_vxo <- sp500_merge_vxo[with(sp500_merge_vxo, order(my)), ]

# ultimately, for every month we only keep the top-most value!
sp500_merge_vxo <- sp500_merge_vxo %>%
                    group_by(my) %>% slice(1)


#--------------------------------------------------------------------
#           Reading-in and transformation of
#           NBER RECESSION DATES:
#--------------------------------------------------------------------

# Note: at multiple occasions below we want to add in the NBER 
# recession dates which we retrieved from 
# https://fred.stlouisfed.org/series/USREC?utm
# _source=series_page&utm_medium=related_content&utm_term=other_
# formats&utm_campaign=other_format;
# therefore we read in the data and perform some manipulations to
# be able to seamlessly integrate the recession dates (year-month
# combinations market with 1) into our dataset for the uncertainty
# measures:
nber_recessions <- read.csv(file="NBER_Recessions_USA.csv", 
                            header=TRUE, sep=",")
# next, we split the variable 'DATE' into 'year', 'month' and 'day'
# using # the 'separate()' - function from the 'tidyr' - package:
nber_recessions <- separate(nber_recessions, "DATE", 
                            c("year", "month", "day"), sep = "-", 
                            remove=FALSE, convert=TRUE)
# we drop the 'day' - variable
nber_recessions <- nber_recessions %>%
                  dplyr::select(-day)

# and create the variable 'my'
nber_recessions <- as.data.frame(nber_recessions %>%
                            mutate(my = year + month/12))


# we quickly rename 'DATE' to 'Date'
nber_recessions <- nber_recessions %>%
                  dplyr::rename(Date = DATE) %>%
                  dplyr::mutate(Date = as.Date(Date))

# the procedure that follows is identical to the procedure
# that we use in another script to extract the exact start- 
# and end-dates of the bloom-shocks:

# let us first filter for all rows where the recession indicator = 1 
# to get an overview:
as.data.frame(nber_recessions %>% filter(USREC == 1))


# then we loop through all rows, and first check if the USREC equals 1
# or not:

# we assign nber_recessions to a new data.frame because
# we want to leave 'nber_recessions' as is to be possibly used 
# by other scripts!

nber_recessions2 <- nber_recessions
# we initialize the grouping variable x as follows:
x <- 2

# we add an empty column to our data frame nber_recessions:
nber_recessions2["recession_ID"] <- NA

# then we start the loop:
for (i in 1:nrow(nber_recessions2)) {
  if(nber_recessions2$USREC[i] == 1) {
    # if we detect a shock (i.e., 'USREC' == 1), then we have 
    # to proceed as follows:
    # we store the 'my' variable in the column 'start' and 'end'
    nber_recessions2[i, ncol(nber_recessions2)] <- x
    
  }else{
    # no shock
    x <- x + 1
  }
}


# we filter to see whether our ID-assignment indeed worked:
nber_recessions2 %>%
                    filter(!is.na(recession_ID))  %>%
                    group_by(recession_ID) %>%
                    summarise(Count = n())

# next we group_by recession_ID and extract the respective 
# start-dates and store the retrieved data to a new 
# data-frame
recessions_start <- as.data.frame(nber_recessions2 %>%
                    filter(!is.na(recession_ID))  %>%
                    group_by(recession_ID) %>%
                    filter(row_number()==1) %>%
                    mutate(year_start=year, month_start=month, 
                           my_start = my) %>%
                    dplyr::select(recession_ID, year_start, 
                                  month_start, my_start))

# next, we replicate the above query for the end-dates
# (note that in some scenarios start- and end-dates 
# might be identical!)
recessions_end <- as.data.frame(nber_recessions2 %>%
                    filter(!is.na(recession_ID))  %>%
                    group_by(recession_ID) %>%
                    filter(row_number()==n()) %>%
                    mutate(year_end=year, month_end=month, 
                           my_end = my) %>%
                    dplyr::select(recession_ID, year_end, 
                                  month_end, my_end))

# next, we can merge the two data-frames from above:
recessions_start_end <- merge(x = recessions_start, 
                              y = recessions_end, 
                              by = "recession_ID", 
                              all = TRUE, na.rm=T)
# and inspect the resulting data-frame:
recessions_start_end


# we re-arrange the sequence of columns
recessions_start_end <- recessions_start_end %>%
                  dplyr::select(recession_ID, year_start, 
                              month_start, year_end, 
                              month_end, my_start, 
                              my_end)


# next, we construct a year-month-variable both out of
# the pair 'year_start' & 'month_start' and 'year_end' & 'month_end'
recessions_start_end$yearmon_start <- as.yearmon(
                          recessions_start_end$my_start, "%Y-%B")
recessions_start_end$yearmon_end <- as.yearmon(
                          recessions_start_end$my_end, "%Y-%B")
# next, we add a helper-column that we need in the next stage as well:
recessions_start_end$helper_date <- format(
                          recessions_start_end$yearmon_start, "%b")

# this means that we can now drop the 'year' and 'month' variables
recessions_start_end <- recessions_start_end %>%
                          dplyr::select(-c(year_start, year_end, 
                                           month_start, month_end))

# next we add a column that gives us the duration of the shock:
recessions_start_end <- recessions_start_end %>%
                          mutate(duration = (
                                recessions_start_end$yearmon_end - 
                                recessions_start_end$yearmon_start)
                                * 12 + 1) %>%
                          dplyr::select(-helper_date)


# Lastly: actually, there is a slight problem, which we fix by adding
# one month's numeric value to my_end:
recessions_start_end$my_end <- recessions_start_end$my_end + (1/12)
# this makes sure that if we refer to an end-month that
# we assume that we are talking about the last day of that
# respective month!
#--------------------------------------------------------------------

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
  theme(legend.position = "bottom", axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))  + 
  # we add in recession dates
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=Inf), 
            fill='#606060', alpha=0.5)
  

volatility

# and we save the plot to be used in our latex-document
ggsave("volatility.pdf")



###############################
## (1.5) detrend the constructed time-series and generate
## Bloom-shocks
###############################
# first we have to make a time series out of sp500_merge_vxo$mvol;
# we do this using the funtionality from the 'zoo' - package;
# we discovered that this part is actually not necessary!
# Therefore, it is commented out for now!
mvol_ts <- zoo(sp500_merge_vxo$mvol, 
               as.Date(sp500_merge_vxo$date))

# next we run the hp-filter
mvol_hp <- hpfilter(sp500_merge_vxo$mvol, type="lambda", 
                    freq=129600)

# next we can add the 'trend' and 'cycle' component of 
# 'mvol_hp' to our existing data-frame:
sp500_merge_vxo$mvol_trend <- mvol_hp$trend
sp500_merge_vxo$mvol_cycle <- mvol_hp$cycle

# next we can plot everything we had plotted above
volatility_trend <- ggplot() +
  # geom_point() + 
  geom_line(data = sp500_merge_vxo, aes(x = my, y = mvol_trend), 
            color="red", size=0.8) +
  geom_line(data = sp500_merge_vxo, aes(x = my, y = mvol), 
            color="black", size=0.8) +
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
  scale_y_continuous(name = "annualized standard deviation (%), 
                     deviations from trend", 
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

#------------------------Initial Creation of Bloom-Shocks--------------------#
# Note: In the below, the Bloom-shocks are only generated for the
# VXO series while in the script 'Identification_of_Shocks.R' they are
# generated for all uncertainty-meaures;

# According to Bloom (2009): The shocks are chosen as those events with 
# stock-market volatility more than 1.65 standard deviations above the 
# HP-detrended mean of the mean of the stock-market volatility series.

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
  scale_y_continuous(name = "annualized standard deviation (%), 
                     deviations from trend", 
                     limits = c(-20, 50), 
                     breaks = seq(-20, 50, by = 10), 
                     minor_breaks = NULL) + 
  # geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=mean_mvol_cycle+sd_mvol_cycle*1.65, 
             
             linetype="dashed", 
             color = "red", size=0.4) +
  # theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), 
        axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

volatility_cycle_shocks

ggsave("volatility_cycle_shocks.pdf")

# the red dotted line in our plot above marks the following 
# threshold:
# yintercept=mean_mvol_cycle+sd_mvol_cycle*1.65, i.e., it is 
# exactly what Bloom (2009) refers to as the threshold which 
# was selected as the
# 5% one-tailed significance level TREATING EACH MONTH AS AN 
# INDEPENDENT OBSERVATION!

# we store the threshold value into our dataset to create an 
# indicator variable that allows us to extract the 'Bloom' - 
# shocks which we want to have a look at
sp500_merge_vxo$thresh <- mean_mvol_cycle+sd_mvol_cycle*1.65

# the above variable now allows us to create an indicator of whether
# or not a certain point is above the threshold value
sp500_merge_vxo <- sp500_merge_vxo %>%
                      mutate(bloom_shock = 
                               case_when(
                                 mvol_cycle > thresh ~ 1,
                                 mvol_cycle <= thresh ~ 0))


###############################
## (1.6) replicate plot from (4) and add
## periods of 'Bloom-shocks'
###############################

# having the variable bloom_shock in our dataset now allows us to 
# replicate our above graph from (4) and add episodes into the graph that 
# correspond to the Bloom-shock being 1

# but before we can do that, we need to create 'start' and 'end' - dates 
# for the respective periods (which we will then ultimately use for 
# plotting the episodes with shaded regions!)

# to retrieve the respective 'start'- and 'end' - dates for the 
# shock-periods, we apply a procedure to the dataset 'sp500_merge_vxo' 
# that results in a separate dedicated data-frame with the following 
# variables: 'start', 'end', 'max_volatility', 'first_volatility';
# once we have the above four variables, we can easily calculate 
# 'duration_months'
# by subtracting 'start' from 'end'.

# the construction of such a table will allow us
# (1) to automatically generate the respective data even if something 
#     in our dataset changes at some point (i.e., it is not hard-coded)
# (2) we can automatically export the information to a latex-table
# to be included in our text by means of e.g., the R-package
# 'miscFuncs' that is capable of exporting matrix-like data
# to latex-tables!

# let us first filter for all rows where the bloom_shock = 1 to get an 
# overview:
as.data.frame(sp500_merge_vxo %>% filter(bloom_shock == 1))

# initially, we had thought of creating a separate data-frame
# and store the loop's data into it;
# but then we decided for a different approach:
# we want to give each episode of a shock a unique ID so that
# we then can easily calculate the max volatility and first
# volatility by grouping by ID (within the full sp500_merge_vxo table!)



# we decided for the following procedure:
# we loop through all rows, and first check if the bloom_shock 
# equals 1 or not:

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

# the above method is maybe not yet the most elegant/efficient 
# solution but it suffices for now (we can fine-tune it at a later 
# stage!)
# we can now have a look at the newly created variable
# (note that we have suppressed 'NA's for the creation of the 
# table!):

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
                        mutate(year_start=year, 
                               month_start=month, 
                               my_start = my) %>%
                        dplyr::select(shock_ID, year_start, 
                                      month_start, my_start))

# next, we replicate the above query for the end-dates
# (note that in some scenarios start- and end-dates are identical!)
shocks_end <- as.data.frame(sp500_merge_vxo %>%
                        filter(!is.na(shock_ID))  %>%
                        group_by(shock_ID) %>%
                        filter(row_number()==n()) %>%
                        mutate(year_end=year, 
                               month_end=month, 
                               my_end = my) %>%
                        dplyr::select(shock_ID, year_end, 
                                      month_end, my_end))
                                    
# next, we can merge the two data-frames from above:
shocks_start_end <- merge(x = shocks_start, 
                          y = shocks_end, 
                          by = "shock_ID", all = TRUE, 
                          na.rm=T)
# and inspect the resulting data-frame:
shocks_start_end

# we re-arrange the sequence of columns
shocks_start_end <- shocks_start_end %>%
                  dplyr::select(shock_ID, year_start, 
                                month_start, year_end, 
                                month_end, my_start, my_end)


# next, we construct a year-month-variable both out of
# the pair 'year_start' & 'month_start' and 'year_end' & 'month_end'
shocks_start_end$yearmon_start <- as.yearmon(
                  shocks_start_end$my_start, "%Y-%B")
shocks_start_end$yearmon_end <- as.yearmon(
                  shocks_start_end$my_end, "%Y-%B")
# next, we add a helper-column that we need in the next stage as well:
shocks_start_end$helper_date <- format(
                  shocks_start_end$yearmon_start, "%b")

# this means that we can now drop the 'year' and 'month' variables
shocks_start_end <- shocks_start_end %>%
                  dplyr::select(-c(year_start, year_end, 
                                   month_start, month_end))

# next we add a column that gives us the duration of the shock:
shocks_start_end <- shocks_start_end %>%
                  mutate(duration = (
                    shocks_start_end$yearmon_end - 
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
shocks_start_end <- merge(x = shocks_start_end, 
                          y = shocks_max_vol, 
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
  scale_x_continuous(name = "Year", limits = c(1962, 2019), 
                     breaks = seq(1962, 2019, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "annualized standard deviation (%), 
                     deviations from trend", 
                     limits = c(-20, 50), 
                     breaks = seq(-20, 50, by = 10), 
                     minor_breaks = NULL) + 
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=mean_mvol_cycle+sd_mvol_cycle*1.65, 
             linetype="dashed", 
             color = "red", size=0.4) +
  #theme_minimal() + 
  labs(color=NULL) +
  theme(legend.position = c(0.93, 0.90), 
        axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14), 
        panel.grid.major.x = element_blank()) + 
  # note that panel.grid.major.x = element_blank() suppresses vertical grid lines!
  geom_rect(data=shocks_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, 
                ymax=+Inf), 
            fill='red', alpha=0.5)
  # we also want to add in indications of NBER recessions
  # geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            # aes(xmin=my_start, xmax=my_end, ymin=30, ymax=+Inf), 
            # fill='#606060', alpha=0.5) +   
  # geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            # aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=-10), 
            # fill='#606060', alpha=0.5)

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
              mutate(max_vol = percent(max_vol), 
                     duration = paste(as.character(floor(duration)), 
                                      "months", sep=" ")) %>%
              dplyr::select(yearmon_start, 
                            yearmon_end, 
                            duration, 
                            yearmon_max, 
                            max_vol)
                
# adapt the below command to supply a vector of column-headers!
latextable(shocks_table, dp=3)


########################################################################
### PART 2: loading and plotting of EPU from Baker, Bloom and Davis (2016)
### See also their dedicated homepage with latest updates to the 
### data at: http://www.policyuncertainty.com/us_historical.html
########################################################################

## PART 2 is separated into ....... parts:
## (2.1)  2.1 handles the original series from 1985 - 2018; note that
## on their homepage, Baker et al. (2016) report both old values from 
## their index calculations and new ones;
## (2.2)  2.2 handles the historical data as of 1900 (- 2018); note that
## on their homepage, Baker et al. (2016) report both old values from 
## their index calculations and new ones;


###############################
## (2.1) loading EPU_USA_Bakeretal2016.xlsx;
##       the .xlsx-file was directly downloaded from
##       http://www.policyuncertainty.com/us_historical.html
##       and on the first sheet the columns 'Year', 'Month' and
##       'News_Based_Policy_Uncert_Index' are read in
###############################

## we import the data (sheet 'Main Index' only)
epu_index <- read_excel("EPU_USA_Bakeretal2016.xlsx", 
                        sheet = "Main Index")

## and generate the variable 'my'
## (note that we first have to convert 'Year' to numeric)
epu_index$Year <- as.numeric(epu_index$Year)
## we rename Year and Month to year and month
epu_index <- epu_index %>%
                    dplyr::rename(year = Year,
                                 month = Month)

epu_index <- as.data.frame(epu_index %>%
                      mutate(my = year + month/12))

## we remove the last row that only consists of NAs
epu_index <- epu_index %>%
                        drop_na

## and finally plot the series
epu_index_plot <- ggplot(epu_index, 
                     aes(x = my, 
                         y = News_Based_Policy_Uncert_Index)) +
  # geom_point() + 
  geom_line(size=0.8) + 
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = "Year", limits = c(1985, 2018), 
                     breaks = seq(1985, 2018, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Policy Uncertainty Index", 
                     limits = c(0, 330), 
                     breaks = seq(0, 330, by = 50), 
                     minor_breaks = NULL) +
  theme(legend.position = c(0.935, 0.93), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14)) + 
  annotate("text", x=1987.833, y=200, label="Black Monday") + 
  annotate("text", x=1991.100, y=230, label="Gulf \n War I") + 
  annotate("text", x=1992.917, y=185, label="Clinton \nElection") + 
  annotate("text", x=1998, y=185, label="Russian \n Crisis/LTCM") +
  annotate("text", x=2000.9, y=210, label="Bush Election", angle=90) + 
  annotate("text", x=2001.9, y=290, label="9/11") +
  annotate("text", x=2003.3, y=240, label="Gulf \n War II") + 
  annotate("text", x=2007.5, y=200, label="Stimulus \n Debate") +
  annotate("text", x=2008.7, y=265, label="Lehman \n and TARP") + 
  annotate("text", x=2010.55, y=215, label="Euro \n Crisis") +
  annotate("text", x=2011.6, y=320, label="Debt \n Ceiling \n Dispute") + 
  annotate("text", x=2013, y=250, label="Fiscal Cliff", angle=90) + 
  annotate("text", x=2013.8, y=280, label="Govt. Shutdown", angle=90) + 
  annotate("text", x=2015.667, y=230, label="EU Migration Crisis", angle=90) + 
  annotate("text", x=2016.42, y=280, label="Coup Turkey", angle=90) + 
  annotate("text", x=2016.99, y=300, label="Trump Election", angle=90) + 
  # we also plot the recession dates
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5)

epu_index_plot

# and we save the plot to be used in our latex-document
ggsave("epu_index_plot.pdf")


###############################
## (2.2) loading EPU_USAHistorical_Bakeretal2016.xlsx;
##       the .xlsx-file was directly downloaded from
##       http://www.policyuncertainty.com/us_historical.html
##       and on the first sheet the columns 'Year', 'Month' and
##       'News-Based Historical Economic Policy Uncertainty' 
##       are read in
###############################
## we import the data (sheet 'Historical EPU' only)
epu_index_historical <- read_excel("EPU_USAHistorical_Bakeretal2016.xlsx", 
                                    sheet = "Historical EPU")

## next, we rename the column with the index, Year and Month
## and generate the variable 'my':
epu_index_historical <- epu_index_historical %>%
  dplyr::rename(year = Year,
               month = Month,
      EPU_Historical = 'News-Based Historical Economic Policy Uncertainty')

epu_index_historical <- as.data.frame(epu_index_historical %>%
                        mutate(my = year + month/12))

## we remove the last row that only consists of NAs
epu_index_historical <- epu_index_historical %>%
                        drop_na


# we also want to output the seven largest spikes in a table
epu_largest <- epu_index_historical %>% 
                      filter(EPU_Historical >= 280) %>%
                      mutate(ID = row_number()) %>%
                      dplyr::select(ID, year, month, 
                                    EPU_Historical, my)

# and create a latex-table out of it:
latextable(epu_largest[, 1:4],
           caption="7 Largest Spikes of U.S. EPU Historical",
            colnames = c("ID", "Year", "Month", 
                      "EPU Historical Index"),
            dp = 4)

# because we want to mark the 7 largest spikes in our plot below, 

## and finally plot the series
epu_index_historical_plot <- ggplot(epu_index_historical, 
                     aes(x = my, y = EPU_Historical)) +
  # geom_point() + 
  geom_line(size=0.4) + 
  # we add points to the 7 largest spikes
  geom_point(data=epu_largest, aes(x=my, y=EPU_Historical), 
             colour="red", size=2) + 
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = "Year", limits = c(1900, 2018), 
                     breaks = seq(1900, 2018, by = 10),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Policy Uncertainty Index", 
                     limits = c(0, 400), 
                     breaks = seq(0, 400, by = 50), 
                     minor_breaks = NULL) +
  theme(legend.position = c(0.935, 0.93), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14)) +
  # lastly we add the IDs of the 7 events we have isolated
  geom_text(data = epu_largest, aes(x = my+2, y = EPU_Historical, label = ID), 
            colour="red", size=4) +
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5)

epu_index_historical_plot

# and we save the plot to be used in our latex-document
ggsave("epu_index_historical_plot.pdf")


#----------------------------------------------------------------------------
# further analysis with the categorical EPU data availalbe in the file
# EPU_USA_Categorical_Data.xlsx
categorical_epu_indices <- read_excel("EPU_USA_Categorical_Data.xlsx", 
                                    sheet = "Indices")

str(categorical_epu_indices)
# first we need to change the format of our 'Date' - variable:
categorical_epu_indices <-  categorical_epu_indices %>%
                        mutate(Date = as.Date(
                          as.numeric(Date), origin = "1899-12-30"))

categorical_epu_indices[, 2:13] <- sapply(
                        categorical_epu_indices[, 2:13], as.character)
categorical_epu_indices[, 2:13] <- sapply(
                        categorical_epu_indices[, 2:13], as.numeric)

# create groups based on date-ranges
# for this we use the cut-function:
# specify levels:
# at the same time, at the end of our command we remove rows which are
# beyond 2014/12
levels_epu_groups <- as.Date(c("1985/01/01", "1990/07/01", "1992/01/01", 
                               "2001/09/01", "2003/01/01", "2007/07/01", 
                               "2008/09/01", "2010/01/01", "2013/11/01"))
# specify labels:
labels_epu_groups <- c("Mid-80s to \n Gulf War I", "Gulf War I", 
                       "1990s boom \n to 9/11", "9/11 \n attacks", 
                       "2000s \n boom", "Early \n credit \n crunch", 
                       "Lehman \n collapse & \n  recession", 
                       "Fiscal \n policy \n battles")
categorical_epu_indices <-  categorical_epu_indices %>% 
                  mutate(groups = cut(Date, breaks=levels_epu_groups, 
                                      labels=labels_epu_groups)) %>%
                  filter(Date <= "2014/12/01")

# we slightly rename the column names
categorical_epu_indices <- categorical_epu_indices %>%
  dplyr::rename('Economic_Policy_Uncertainty' = '1. Economic Policy Uncertainty',
                'Fiscal_Policy' = 'Fiscal Policy (Taxes OR Spending)',
                'Taxes' = '3. Taxes',
                'Government_Spending_&_Other' = '4. Government spending',
                'Monetary_Policy' = '2. Monetary policy',
                'Health_Care' = '5. Health care',
                'National_Security' = '6. National security',
                'Regulation' = '8. Regulation',
                'Financial_Regulation' = 'Financial Regulation',
                'Sovereign_Debt_&_Currency_Crises' = '10. Sovereign debt, currency crises',
                'Entitlement_Programs' = '7. Entitlement programs',
                'Trade_Policy' = '9. Trade policy')
#----------------------------------------------------------------------------


########################################################################
### PART 3: loading and construction of Consumer Uncertainty based off
###         of data from the Michigan Survey following Leduc and Liu (2016).
###         The data off of which Leduc and Lui construct their index
###         is available at https://data.sca.isr.umich.edu/subset/subset.php;
###         In particular, the data derived from the questions 
###         'Buying Conditions for Vehicles' (VEH) and the follow-up
###         questions 'Reasons for Opinions for Buying Conditions for Vehicles'
###         (VEHRN) are used.
########################################################################

## PART 3 is separated into ....... parts:
## (3.1)  3.1 handles ........
## (3.2)  3.2 handles the .......


###############################
## (3.1) loading of downloaded data;
##       the michigan_survey.csv-files were directly downloaded from
##       https://data.sca.isr.umich.edu/subset/subset.php
##       by ticking 'VEH', 'VEHRN', 'Produce all response
##       codes, including index scores', 'Monthly', 
##       'Starting year = 1978', 'Ending year = 2018',
##       'All Households (aggregate; same as our tables & charts)'
###############################

michigan_survey <- read.csv(file="michigan_survey.csv", 
                            header=TRUE, sep=";",
                            stringsAsFactors = FALSE)

# corresponding to our selection of data we want to download on the
# website https://data.sca.isr.umich.edu/subset/subset.php, the time-series
# starts in January 1978 (note that we actually have to start in February
# 1978 because the follow-up question has no data in January 1978).

# To know what the variables stand for, see the Codebook at
# https://data.sca.isr.umich.edu/fetchdoc.php?docid=45121.

## make 'yyyymm' - variable a 'Date' - type 
## (instead of int)
michigan_survey$yyyymm <- as.Date(paste0(as.character(
                    michigan_survey$yyyymm), '01'), format='%Y%m%d')
## and rename the variable
michigan_survey <- michigan_survey %>%
                rename(Date = yyyymm)

# next, we create the three variable 'month', 'year' and 'day' using
# the 'separate()' - function from the 'tidyr' - package
michigan_survey <- separate(michigan_survey, "Date", 
                      c("year", "month", "day"), sep = "-", 
                      remove=FALSE, convert=TRUE)

# next we create the variable 'my' which is a numerical representation
# of yearmon:
michigan_survey <- as.data.frame(michigan_survey %>%
                             mutate(my = year + month/12))
# and the variable veh_share_unc which will is calculated
# as vehrn_fb_all/100
michigan_survey <- michigan_survey %>%
                mutate(veh_share_unc = vehrn_fb_all)

# before plotting, we calculate the three-month-MA out of the series
# veh_share_unc like it is done in Leduc and Liu (2016):
michigan_survey <- michigan_survey %>%
                dplyr::mutate(veh_share_unc = zoo(veh_share_unc)) %>%
                dplyr::mutate(veh_share_unc.MA = rollapply(
                  veh_share_unc, width=3, mean, fill = NA, 
                  align = "right", na.rm = TRUE
                )) %>%
                # in a last step we transfer column 'veh_share_unc.MA'
                # back to a numerical data-type:
                dplyr::mutate(veh_share_unc.MA = as.numeric(veh_share_unc.MA))


## and finally plot the series
michigan_plot <- ggplot() +
  geom_line(data = michigan_survey, 
            aes(x = my, 
                y = veh_share_unc.MA), 
            color="black", 
            size=1,linetype = 1) +
  scale_x_continuous(name = "Year", limits = c(1978, 2019), 
                     breaks = seq(1980, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Consumer Uncertainty", 
                     limits = c(0, 15), 
                     breaks = seq(0, 15, by = 5), 
                     minor_breaks = NULL) +
  theme(legend.position = c(0.935, 0.93), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14)) +
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5)

michigan_plot

# and we save the plot to be used in our latex-document
ggsave("michigan_plot.pdf")

########################################################################
### PART 4: loading and plotting of GTU and GT-index
### from Bontempi et al. (2015) and
### Castelnuovo and Tran (2017); 
### See also Efrem Castelnuovo's homepage
### for latest updates to their GTU-index at
### https://sites.google.com/site/efremcastelnuovo/home/publications
########################################################################
## we import the data (sheet 'Foglio1' only)
GTU_index <- read_excel("GTU_indices_EcLetts.xls", 
                                     sheet = 1)


## make 'month' - variable a 'Date' - type 
## (instead of POSIXct)
GTU_index$Month <- as.Date(GTU_index$Month)

## then we rename 'Month' to 'Date'
GTU_index <- rename(GTU_index,
                    Date = Month)

# next, we create the three variable 'month', 'year' and 'day' using
# the 'separate()' - function from the 'tidyr' - package
GTU_index <- separate(GTU_index, "Date", 
                      c("year", "month", "day"), sep = "-", 
                      remove=FALSE, convert=TRUE)

# next we create the variable 'my' which is a numerical representation
# of yearmon:
GTU_index <- as.data.frame(GTU_index %>%
                            mutate(my = year + month/12))


## and finally plot the series
GTU_index_plot <- ggplot() +
  geom_line(data = GTU_index, 
            aes(x = my, 
                y = GTU_US), 
            color="black", 
            size=1,linetype = 1) +
  scale_x_continuous(name = "Year", limits = c(2004, 2016), 
                     breaks = seq(2004, 2016, by = 1),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "GTU index, US", 
                     limits = c(0, 250), 
                     breaks = seq(0, 250, by = 50), 
                     minor_breaks = NULL) +
  theme(legend.position = c(0.935, 0.93), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14))

GTU_index_plot

# and we save the plot to be used in our latex-document
ggsave("GTU_index_plot.pdf")


########################################################################
### PART 5: loading and plotting of macro uncertainty index
### from Jurado et al. (2015); See also Sydney Ludvigson's homepage
### for latest updates to the data at
### https://www.sydneyludvigson.com/data-and-appendixes/
########################################################################

## PART 5 is separated into ....... parts:
## (5.1)  5.1 handles ........
## (5.2)  5.2 handles the .......


###############################
## (5.1) loading MacroUncertaintyToCirculate.xlsx;
##       the .xlsx-file was directly downloaded from
##       https://www.sydneyludvigson.com/data-and-appendixes/
##       and comes along together with another data-source
##       called 'FinancialUncertaintyToCirculate' (not used here!);
##       on the very first sheet called 'Macro Uncertainty' the columns
##       Date, h=1, h=3 and h=12 are read in.
###############################

## we import the data (sheet 'Macro Uncertainty' only)
macroUncertainty_index <- read_excel("MacroUncertaintyToCirculate.xlsx", 
                        sheet = "Macro Uncertainty")

## and rename three columns
macroUncertainty_index <- rename(macroUncertainty_index, 
                          h1 = "h=1",
                          h3 = "h=3",
                          h12 = "h=12")

## make 'Date' - variable a 'Date' - type 
## (instead of POSIXct)
macroUncertainty_index$Date <- as.Date(macroUncertainty_index$Date)

# next, we create the three variable 'month', 'year' and 'day' using
# the 'separate()' - function from the 'tidyr' - package
macroUncertainty_index <- separate(macroUncertainty_index, "Date", 
                            c("year", "month", "day"), sep = "-", 
                            remove=FALSE, convert=TRUE)

# next we create the variable 'my' which is a numerical representation
# of yearmon:
macroUncertainty_index <- as.data.frame(macroUncertainty_index %>%
                                   mutate(my = year + month/12))

# for all three respective series (h=1, h=3 and h=12), we calculate
# the mean and standard deviation to be able to add horizontal
# lines marking 1.65 standard deviations above the mean to the plot
sd_h1 <- sd(macroUncertainty_index$h1)
sd_h3 <- sd(macroUncertainty_index$h3)
sd_h12 <- sd(macroUncertainty_index$h12)
mean_h1 <- mean(macroUncertainty_index$h1)
mean_h3 <- mean(macroUncertainty_index$h3)
mean_h12 <- mean(macroUncertainty_index$h12)

# we also want to add shaded areas to our plot (similar to the Bloom-
# shocks we have added above), that denote episodes where the macro
# uncertainty series exeeds 1.65 standard deviations above its mean
# (we will see that we have far fewer uncertainty episodes than other
# popular proxies for uncertainty) including NBER
# recession dates:

# for this, we store the threshold value into our dataset to 
# create an indicator variable that allows us to extract
# the 'shocks' according to the macro uncertainty series
# (note that we perform this only for the h=1 series!!)
macroUncertainty_index$thresh <- mean_h1+sd_h1*1.65

# the above variable now allows us to create an indicator of whether
# or not a certain point is above the threshold value
macroUncertainty_index <- macroUncertainty_index %>%
            mutate(macro_shock = case_when(h1 > thresh ~ 1,
                                 h1 <= thresh ~ 0))

# having the variable macro_shock in our dataset now allows us to 
# add episodes into the graph that correspond to the
# macro-shock being 1;

# but before we can do that, we need to create 'start' and 'end' - dates for
# the respective periods (which we will then ultimately use for plotting
# the episodes with shaded regions!)

# to retrieve the respective 'start'- and 'end' - dates for the shock-periods,
# we apply a procedure to the dataset 'macroUncertainty_index' that results in a 
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
as.data.frame(macroUncertainty_index %>% filter(macro_shock == 1))


# initially, we had thought of creating a separate data-frame
# and store the loop's data into it;
# but then we decided for a different approach:
# we want to give each episode of a shock a unique ID so that
# we then can easily calculate the max volatility and first
# volatility by grouping by ID (within the full macroUncertainty_index!)

# we decided for the following procedure:
# we loop through all rows, and first check if the macro_shock equals 1
# or not:

# we initialize the grouping variable x as follows:
x <- 2

# we add an empty column to our data frame macroUncertainty_index:
macroUncertainty_index["shock_ID"] <- NA

# then we start the loop:
for (i in 1:nrow(macroUncertainty_index)) {
  if(macroUncertainty_index$macro_shock[i] == 1) {
    # if we detect a shock (i.e., 'bloom_shock' == 1), then we have 
    # to proceed as follows:
    # we store the 'my' variable in the column 'start' and 'end'
    macroUncertainty_index[i, ncol(macroUncertainty_index)] <- x
    
  }else{
    # no shock
    x <- x + 1
  }
}


# the above method is maybe not yet the most elegant/efficient solution
# but it suffices for now (we can fine-tune it at a later stage!)
# we can now have a look at the newly created variable
# (note that we have suppressed 'NA's for the creation of the table!):

macroUncertainty_index %>%
  filter(!is.na(shock_ID))  %>%
  group_by(shock_ID) %>%
  summarise(Count = n())

# next we group_by shock_ID and extract the respective start-dates and
# store the retrieved data to a new data-frame
macro_shocks_start <- as.data.frame(macroUncertainty_index %>%
                    filter(!is.na(shock_ID))  %>%
                    group_by(shock_ID) %>%
                    filter(row_number()==1) %>%
                    mutate(year_start=year, 
                           month_start=month, 
                           my_start = my) %>%
                    dplyr::select(shock_ID, year_start, 
                                  month_start, my_start))


# next, we replicate the above query for the end-dates
# (note that in some scenarios start- and end-dates are identical!)
macro_shocks_end <- as.data.frame(macroUncertainty_index %>%
                    filter(!is.na(shock_ID))  %>%
                    group_by(shock_ID) %>%
                    filter(row_number()==n()) %>%
                    mutate(year_end=year, 
                           month_end=month, 
                           my_end = my) %>%
                    dplyr::select(shock_ID, year_end, 
                                  month_end, my_end))

# next, we can merge the two data-frames from above:
macro_shocks_start_end <- merge(x = macro_shocks_start, 
                                y = macro_shocks_end, 
                                by = "shock_ID", 
                                all = TRUE, na.rm=T)
# and inspect the resulting data-frame:
macro_shocks_start_end

# we re-arrange the sequence of columns
macro_shocks_start_end <- macro_shocks_start_end %>%
                  dplyr::select(shock_ID, 
                                year_start, 
                                month_start, 
                                year_end, 
                                month_end, 
                                my_start, 
                                my_end)


# next, we construct a year-month-variable both out of
# the pair 'year_start' & 'month_start' and 'year_end' & 'month_end'
macro_shocks_start_end$yearmon_start <- as.yearmon(
                      macro_shocks_start_end$my_start, "%Y-%B")
macro_shocks_start_end$yearmon_end <- as.yearmon(
                      macro_shocks_start_end$my_end, "%Y-%B")
# next, we add a helper-column that we need in the next stage as well:
macro_shocks_start_end$helper_date <- format(
                      macro_shocks_start_end$yearmon_start, "%b")


# this means that we can now drop the 'year' and 'month' variables
macro_shocks_start_end <- macro_shocks_start_end %>%
  dplyr::select(-c(year_start, year_end, month_start, month_end))

# next we add a column that gives us the duration of the shock:
macro_shocks_start_end <- macro_shocks_start_end %>%
  mutate(duration = (macro_shocks_start_end$yearmon_end - 
                       macro_shocks_start_end$yearmon_start) 
                              * 12 + 1) %>%
  dplyr::select(-helper_date)

# next we group_by shock_ID and extract the respective MAXIMUM
# VOLATILITY
# (note that the produced data-frame only has one row per shock_ID)
macro_shocks_max_vol <- as.data.frame(macroUncertainty_index %>%
                        filter(!is.na(shock_ID))  %>%
                        group_by(shock_ID) %>%
                        filter(h1 == max(h1)) %>%
                        mutate(max_h1 = h1) %>%
                        dplyr::select(shock_ID, max_h1, my))

# we apply the same date-transformation as above:
macro_shocks_max_vol$yearmon_max <- as.yearmon(
                        macro_shocks_max_vol$my, "%Y-%B")
# and we drop 'my':
macro_shocks_max_vol <- macro_shocks_max_vol %>%
                    dplyr::select(-my)

# in a last stage, we merge macro_shocks_max_vol with 
# macro_shocks_start_end:
macro_shocks_start_end <- merge(x = macro_shocks_start_end, 
                                y = macro_shocks_max_vol, 
                          by = "shock_ID", all = TRUE, na.rm=T)


# actually, there is a slight problem, which we fix by adding
# one month's numeric value to my_end:
macro_shocks_start_end$my_end <- macro_shocks_start_end$my_end + (1/12)
# this makes sure that if we refer to an end-month that
# we assume that we are talking about the last day of that
# respective month!


# the above data-frame now allows us to plot the episodes of 
# high volatility into our previous plots (by adding a shaded
# rectangle!)
# as well as the NBER recession dates

## and finally plot the series
macroUncertainty_index_plot <- ggplot() +
  geom_line(data = macroUncertainty_index, 
                    aes(x = my, 
                        y = h1), color="red", 
            size=1,linetype = 1) +
  geom_line(data = macroUncertainty_index, 
                    aes(x = my, 
                        y = h3),
            color="darkgreen",
            size=1, 
            linetype = 1) +
  geom_line(data = macroUncertainty_index, 
            aes(x = my, 
                y = h12),
            color="blue",
            size=1, 
            linetype = 1) +
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Macro Uncertainty Indices", 
                     limits = c(0.5, 1.3), 
                     breaks = seq(0.5, 1.5, by = 0.2), 
                     minor_breaks = NULL) +
  scale_color_discrete(name = NULL, labels = c("h=1", "h=3", "h=12")) +
  # geom_hline(yintercept=mean_h1+sd_h1*1.65, linetype="dashed", 
  #            color = "blue", size=0.4) +
  # geom_hline(yintercept=mean_h3+sd_h3*1.65, linetype="dashed", 
  #            color = "green", size=0.4) +
  # geom_hline(yintercept=mean_h12+sd_h12*1.65, linetype="dashed", 
  #            color = "red", size=0.4) +
  theme(legend.position = "bottom", axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14)) + 
  # add in macro - shocks (similar to Bloom's approach)
  # geom_rect(data=macro_shocks_start_end, inherit.aes = FALSE,
  #           aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
  #           fill='red', alpha=0.5) + 
  # we also want to add in indications of NBER recessions
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5)
  # geom_rect(data=recessions_start_end, inherit.aes = FALSE,
  #           aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=0.6), 
  #           fill='#606060', alpha=0.5)
  # change ratio of y and x - axis
  # coord_fixed(ratio = 40)
  
macroUncertainty_index_plot



# and we save the plot to be used in our latex-document
ggsave("macroUncertainty_index_plot.pdf")


# having created the data-frame macro_shocks_start_end, we can now also
# export a latex-table with the exact data of the macro-shock!
# first, we slightly adapt the shocks - data frame:
# first, we add '%' to 'max_vol' in our data-frame:
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x, format = format, digits = digits, ...), "%")
}

macro_shocks_table <- macro_shocks_start_end %>%
  mutate(max_h1 = percent(max_h1), 
         duration = paste(as.character(floor(duration)), 
                          "months", sep=" ")) %>%
  dplyr::select(yearmon_start, 
                yearmon_end, 
                duration, 
                yearmon_max, 
                max_h1)

# adapt the below command to supply a vector of column-headers!
latextable(macro_shocks_table, colnames = c("Start", 
                                            "End", "Duration", 
                                            "Maximum Year/Month", 
                                            "Maximum Value"), dp=3)

########################################################################
### PART 6: construction of one time-series plot comparing
### all discussed uncertainty-measures
########################################################################

# To generate our plot with all the series, we first park all the
# time series which we want to plot into one data-frame:

# The necessary data-frames are:
#      * 'macroUncertainty_index' for the measure following Jurado et al. (2015)
#      * 'GTU_index' for the measure following Castelnuovo and Tran (2017)
#      * 'epu_index' for the EPU-index following Baker et al. (2016)
#      * 'sp500_merge_vxo' for the VXO (the cycle-part after detrending
#          following Bloom (2009))
#      * 'michigan_survey' for the consumer uncertainty series
#          following Leduc and Liu (2016)

# We combine the relevant variables from the above data-frames into 
# a new data-frame (note that we only take h=1 from the 
# three constructed uncertainty-measures of Jurado et al.):

# Because all the series that we want to have a look at have different
# sample length, we perform a cascading left-join:

# we start with the following data-frame where we store
# the data from the macro-uncertainty-index:
comparison_measures <- data.frame("macroUncert" = macroUncertainty_index$h1,
                                  "my" = macroUncertainty_index$my)

# then we join in the data one after the other, starting with 
# the data from the GTU_index,
# followed by the original epu_index,
# as well as the historical EPU_index,
# then the VXO (the cycle-part after detrending),
# and finally the MSoC in a multi-setp join:

# before we join in all uncertainty-measures, note the following:
# a quick comparison between Bloom's original data (including his VOLATBL
# measures) and our reconstructed measure revealed a small but consistent
# difference of approx 0.6 between 1962 and 1986.
# So, in order to use Bloom's data for the period 1962 - 1986, we
# replace our 'mvol' measure with Bloom's data in the period
# 1962 - 1986;
# i.e., until 'my' == 1986.083!

VARDATA_Bloom_orig <- read.csv(file="VARDATA_Bloom_original.csv", 
                            header=TRUE, sep=",")

# we create a variable called 'my'
VARDATA_Bloom_orig <- VARDATA_Bloom_orig %>%
                            mutate(my = YEAR + MONTH/12)

# we add the variable 'VOLATBL' from Bloom's original
# data to our data-frame
sp500_merge_vxo <- left_join(x = sp500_merge_vxo,
            y= VARDATA_Bloom_orig[, c("VOLATBL", "my")],
            by = "my")
            
# and then create the variable 'mvol_VOLATBL' which is a merge
# out of our 'mvol' - variable and Bloom's 'VOLATBL',
# whereby 'VOLATBL' ranges from 1962 - 1986,
# and mvol as of 1986 - 2018;
sp500_merge_vxo <- sp500_merge_vxo %>%
                      dplyr::mutate(mvol_VOLATBL = 
                                      ifelse(my <= 1986, VOLATBL,
                                             mvol))


comparison_measures <- right_join(x = GTU_index[, c("GTU_US", "my")], 
                                 y = comparison_measures, 
                                by = "my") %>%
                            right_join(x = epu_index[, 
                                  c("News_Based_Policy_Uncert_Index","my")],
                                  y = .,
                                  by = "my") %>%
                            right_join(x = epu_index_historical[, 
                                  c("EPU_Historical","my")],
                                  y = .,
                                  by = "my") %>%
                            # here we switched 'mvol_cycle' to
                            # 'mvol' - we had wrongly
                            # used 'mvol_cycle'
                            right_join(x = sp500_merge_vxo[, 
                                                c("mvol_VOLATBL","my")],
                                  y = .,
                                  by = "my") %>%
                            right_join(x = michigan_survey[, 
                                                      c("veh_share_unc.MA","my")],
                                  y = .,
                                  by = "my")

# Note that the resulting data.frame ranges until the end of 2017 due to
# the macro uncertainty index having exactly that range!
# (i.e., all other data.frames that we read in afterwards, only run
# until the end of 2017 as well!)

# Rename all the variables:
comparison_measures <- comparison_measures %>%
                      dplyr::rename(VXO = mvol_VOLATBL,
                             EPU = News_Based_Policy_Uncert_Index,
                             'EPU Hist' = EPU_Historical,
                             GTU = GTU_US,
                             Macro = macroUncert,
                             Michigan = veh_share_unc.MA)


# we reorder the varaibles
# and ultimately decided to drop the GTU
comparison_measures <- comparison_measures[c(2, 1, 3, 4, 5, 7)]

# we make a time-series object out of our data-frame
# (and exclude the simple 'EPU')
comparison_measures_ts <- ts(comparison_measures[c(2, 3, 4, 6)],
                             start = c(1960, 7), frequency = 12)


# we plot all the time-series by making use of
# of the forecast-package (that reverts back to ggplot2)
# and add in the NBER-dates;

comparison_plot <- autoplot(comparison_measures_ts, facets = TRUE) + 
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Uncertainty Measures")  +
  labs(color=NULL) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14), 
        panel.grid.major.x = element_blank()) +
  # note that panel.grid.major.x = element_blank() 
  # suppresses vertical grid lines!
  geom_rect(data=recessions_start_end, 
            inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, 
                ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5)

comparison_plot
  
ggsave("comparison_plot.pdf")


# Next, we have to normalize all the series to be able to put
# them all on the same scale:
# we use the 'scale' - function on the entire data-frame:
scaled_comparison_measures_ts <- scale(comparison_measures_ts)


# Note:
# autoplot() will pick a sensible default for the object it's passed. 
# If we want to customize its appearance it's probably better to use 
# the standard ggplot() function.
# To be able to do that the zoo object should be passed trough fortify():
# see https://stackoverflow.com/questions/43068203/unable-changing-
# linetype-ggplot2

# because autoplot does not quite do what we want and we cannot
# apply the 'gather()' - function to a time-series object,
# we apply the scaling also to 'comparison_measures' (without '_ts')
# so that we can use a non-time-series data.frame in our plot below:
# (but this time we perform the scaling manually!)
scaled_comparison_measures <- comparison_measures %>%
              dplyr::rename(EPU_Hist = 'EPU Hist') %>%
              dplyr::mutate(VXO = (VXO - mean(VXO, 
                                             na.rm=TRUE))
                                    /sd(VXO, na.rm=TRUE),
                          EPU_Hist  = (EPU_Hist - mean(EPU_Hist, 
                                             na.rm=TRUE))/
                            sd(EPU_Hist, na.rm=TRUE),
                          Michigan  = (Michigan - mean(Michigan, 
                                             na.rm=TRUE))/
                            sd(Michigan, na.rm=TRUE),
                          Macro  = (Macro - mean(Macro, 
                                             na.rm=TRUE))/
                            sd(Macro, na.rm=TRUE))    %>%
               # and we remove EPU:
               dplyr::select(-EPU)
  
# we bring the data into a tidy data-format:
scaled_comparison_measures <- scaled_comparison_measures %>%
              # note that the 'gather()' - function produces a warning
              # message because 'Michigan' enter as a 'zoo' - object!
              gather(uncert_measure, value, -my, na.rm = TRUE)

ggplot(scaled_comparison_measures, 
       aes(x = my, y = value, colour=uncert_measure,
           linetype = uncert_measure, shape=uncert_measure)) +
    geom_line(size=0.6) +
  geom_point(size=1) + 
  scale_color_manual(values=c("#CC0000",
                              "#2ac113",
                              "#202020", "#0452d1", "#FF69B4")) + 
  scale_y_continuous(name = "Uncertainty Measures") +
  scale_x_continuous(name = "Year", limits = c(1960, 2018), 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  theme(legend.position="bottom") + 
  labs(color=NULL) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        #legend.text=element_text(size=14), 
        panel.grid.major.x = element_blank()) + 
  # note that panel.grid.major.x = element_blank() suppresses vertical grid lines!
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5) + 
  scale_linetype_manual(values = c(1,1,1,1)) +
  scale_shape_manual(values=c(0,1,2,3))

ggsave("comparison_plot_combined.pdf")



