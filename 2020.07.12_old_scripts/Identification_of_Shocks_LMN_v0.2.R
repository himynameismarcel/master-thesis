# This script generates Figure 1 of Ludvigson et al 2018:

# The approach is a combination of our very first script and the 
# script Identification_of_Shocks.R

# we read in the data from the data that Sai Ma provided us with
## -------------------------------------------------------------------
### Reading in Data
## -------------------------------------------------------------------
SVAR.data <- read_excel("Replication_data_Ludvigson.xlsx", 
                        sheet = "Data")
# a quick check shows us that the columns are already correctly
# named, but the type for the 'Date'-column is set to num;
# the Date-column consists of a numeric holding 'yearmonth',
# hence we have to proceed as follows to change that:
SVAR.data$yearmon <- as.yearmon(as.character(SVAR.data$Date), "%Y%m")

# then we extract year and month and store them in two separate
# and distinct columns:
SVAR.data <- SVAR.data %>%
            mutate(year = year(yearmon),
                   month = month(yearmon))



## -------------------------------------------------------------------
# we standardize the two columns Uf and Um
# for this, we simply use the scale-function
## -------------------------------------------------------------------
SVAR.data <- SVAR.data %>%
            mutate(Uf = (Uf - mean(Uf))/sd(Uf),
                   Um = (Um - mean(Um))/sd(Um))
# quick test, whether the standardization worked as intended
mean(SVAR.data[["Uf"]])
sd(SVAR.data[["Uf"]])
mean(SVAR.data[["Um"]])
sd(SVAR.data[["Um"]])

# we create a variable 'my':
SVAR.data <- SVAR.data %>%
  mutate(my = year + month/12)

# for the below function, we create two data-frames which 
# hold Uf or Um and to which we can then apply the function to:
identify_shocks_Um <- SVAR.data %>%
              dplyr::select(1, 7, 8, 3) %>%
              # we drop potential NAs
              drop_na %>%
              # we have to rename the uncertainty measure to 'series'
              dplyr::rename(series = Um)
identify_shocks_Uf <- SVAR.data %>%
              dplyr::select(1, 7, 8, 5) %>%
              # we drop potential NAs
              drop_na %>%
              # we have to rename the uncertainty measure to 'series'
              dplyr::rename(series = Uf)

## -------------------------------------------------------------------
# then we use the approach out of 'Identification_of_Shocks.R' to 
# identify 'shocks' (i.e. episodes where the series goes beyond
# 1.65 std. devs)
## -------------------------------------------------------------------
# because we have actually already implemented a function that
# identifies shocks according to bloom (called 'extractShocks_BLOOM' from
# the script 'Identification_of_Shocks.R'), we can leverage the function 
# also for this application.
# The only difference is, that Bloom constructs the shock as follows:
#   (1) detrending the uncertainty measure using a HP-filter with lambda
#   = 129 600;
#   (2) shocks are then chosen as those events with stock-market volatility
#   more than 1.65 standard deviations above the HP-detrended mean
#   (selected as the 5% one-tailed significance level) of the measures'
#   time series; thereby, each month is being treated as an independent
#   observation;

# Instead, in the application at hand, we standardize the respective
# columns (i.e., mean = 0 and sd = 1) and then simply consider
# cases for which the value of the time series is > 1.65;
# hence, we take the function 'extractShocks_BLOOM' and adjust
# it accordingly:
# here the function starts:

# 20200705, MarceL:
# Note that we have moved the function to our script called 
# '20200628_functions.R'.
# Accordingly we have to source this file 20200628_functions.R'
# on top of this script!

extractShocks_LMN <- function(df, name=name){
  
            nam <- NULL
            mean <- NULL        # in our case: 0
            sd <- NULL          # in our case: 1

            
            # (1) detection of shocks (i.e, 0/1)
            # note that we slightly adjust the function here because
            # we do not want to use the HP-filter:
            # we standardize the column 'series'
            df <- df %>%
                        mutate(series = (series - mean(series))/sd(series))
            
            # we create a threshold variable and store it in our dataset to 
            # be able to create an indicator variable that allows us to extract 
            # the 'Bloom' - shocks (this time across all variables) which we
            # want to have a look at!
            
            # to name this additional column also accordingly, we conserve our
            # name-object again:
            nam <- paste(colnames(df[4]), "thresh", sep="_")
            df[, nam] <- 1.65
            
            # now we can create an indicator-variable
            # (for this we conserve our 'nam' variable again):
            nam <- paste(colnames(df[4]), "shock", sep="_")
            df[, nam] <- ifelse(df[, 4] >
                                  df[, 5], 1, 0)
            
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
            df[df[, 6] == 1, c(1, 4, 6)]
            
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
              if(df[i, 6] == 1) {
                # if we detect a shock (i.e., 'bloom_shock' == 1), then 
                # we have to record a shock-ID:
                df[i, 7] <- x
                
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
            
            
            # we continue with a subset of the variables and only look
            # at the rows where we actually have a shock_ID other than 
            # NA!
            # for this we extract the name of the column we want to,
            # among others, filter for:
            
            df.sub <- as.data.frame(df %>%
                            filter(!is.na(shock_ID))  %>%
                            dplyr::select(shock_ID, year, month, my,
                                          series))
            
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
                            shock_ID) %>%
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
                                        factor(shock_ID, 
                                               levels = unique(shock_ID))))
            
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
            return(list(df.sub, df.dates, df.start.end, 
                        df.blo.start.end.max, df.sub.blo))
  
}          
#-----------------------------------ENDOFFUNCTION-----------------------------------#


Um <- extractShocks_LMN(identify_shocks_Um, "Um") 
Uf <- extractShocks_LMN(identify_shocks_Uf, "Uf") 


## -------------------------------------------------------------------
# Creation of Plot
## -------------------------------------------------------------------
Um[[3]]
Uf[[3]]

# we need to stack all second elements of the lists above each other
# to have all shock-periods combined:
shocks_stacked <- bind_rows(Um[[3]], Uf[[3]])

# we remove the column 'shock_ID' and add the columns
# 'ymin' and 'ymax' holding -Inf and +Inf:
shocks_stacked <- shocks_stacked %>%
                  dplyr::mutate(ymin = -Inf,
                                ymax = +Inf) %>%
                  dplyr::select(-shock_ID)


# SVAR.data contains all time-series we want to
# have a look at;
SVAR.data <- SVAR.data %>%
              dplyr::select(my, Um, Uf)

# to bring the dataset into a useful shape for our purposes,
# we 'gather()' the variables with the variable-Names
# together;
# note that before applying the 'gather' - function, we remove
# the series for 'EPU' (because we only want to keep EPU Historical!)
SVAR.data <- SVAR.data %>%
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
Uf[[4]]
Um[[4]]


# we need to stack all second elements of the lists above each other
# to have all shock-periods combined:
max_shocks_stacked <- bind_rows(Um[[4]] %>%
                                  dplyr::mutate(uncert_measure = "Um"), 
                                Uf[[4]] %>%
                                  dplyr::mutate(uncert_measure = "Uf")) %>%
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


# we just have to do the same again with the very last
# list-items:
Uf[[5]]
Um[[5]]


# we need to stack all second elements of the lists above each other
# to have all shock-periods combined:
all_shocks_stacked <- bind_rows(Um[[5]] %>%
                                  dplyr::mutate(uncert_measure = "Um"), 
                                Uf[[5]] %>%
                                  dplyr::mutate(uncert_measure = "Uf"))

# new facets labels:
supp.labs <- c(`Um`= "Aggregate Macro Uncertainty Um", 
               `Uf`= "Aggregate Macro Uncertainty Uf")


# now we should be able to plot the result:
# note: the below plot marks both entire periods with 'significant'
# volatility;
# 
LMN_Shocks_plot_combined <- ggplot(SVAR.data, aes(x=my,
                                                  y=value,
                                                  color = uncert_measure)) +
  geom_rect(data=recessions_start_end[16:23, ], inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=Inf), 
            fill='#606060', alpha=0.3) + 
  geom_line(size = 1) + 
  facet_grid(uncert_measure ~ .,  labeller = as_labeller(supp.labs), scales="free_y") +
  # geom_rect(data=shocks_stacked, aes(x = NULL, y=NULL, xmin=start, xmax=end,
  #                                    ymin=ymin, ymax=ymax),
  #           alpha=0.2, fill="red") +
  # geom_rect(data=max_shocks_stacked, aes(x = NULL, y=NULL, 
  #                                        xmin=yearmon_max_start, 
  #                                        xmax=yearmon_max_end,
  #                                        ymin=ymin, ymax=ymax),
  #           fill="red") +

  # following the latest paper of Jurado et al, 
  # we add dots that identify the shock-events 
  # (in addition to the blue vertical lines!)
  geom_point(data=all_shocks_stacked, aes(x = my, y=series), 
             alpha=0.9, colour="black", size = 2) + 
  scale_y_continuous(name = "Uncertainty Measures") +
  scale_x_continuous(name = "Year", 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  scale_color_manual(values = c("Um" = "blue", 
                                "Uf" = "red")) + 
  theme(legend.position="bottom") + 
  labs(color=NULL) +
  geom_hline(yintercept=1.65, linetype="dashed", color = "black") + 
  theme(legend.position="none",
        axis.text=element_text(size=10),
        #plot.title = element_text(size=10, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  geom_text(aes(x=2015, y=2, label = "1.65 std"))


LMN_Shocks_plot_combined

# ggsave("LMN_Shocks_plot_combined.pdf")

## -------------------------------------------------------------------
# ultimately, we also want to add the NBER recession dates to the
# plot (which we have also done in 'Identification_of_Shocks.R')
## -------------------------------------------------------------------
## see above! we have taken the NBER-recession dates from another 
## script (where we have prepared them);


## calculation of correlation:
cor(SVAR.data$Uf, SVAR.data$Um)
















