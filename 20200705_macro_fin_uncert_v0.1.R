#########################################################################################
### Marcel Kropp, 15.04.2018
### This script is separated into several big parts:
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



###############################
## loading libraries
###############################
if(!exists("foo", mode="function")) source("20200705_load_libraries_v0.1.R")


###############################
## reading in NBER recession dates
###############################
if(!exists("foo", mode="function")) source("20200705_nber_recession_dates_v0.1.R")


########################################################################
### PART 5: loading and plotting of macro uncertainty index
### from Jurado et al. (2015); See also Sydney Ludvigson's homepage
### for latest updates to the data at
### https://www.sydneyludvigson.com/data-and-appendixes/
########################################################################

## PART 5 is separated into ....... parts:
## (5.1)  5.1 handles the macro uncertainty measures
## (5.2)  5.2 handles the .financial uncertainty measure (both from LMN)


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
                                          mutate(my = year + month/12 - 1/12))

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
  geom_hline(yintercept=mean_h1+sd_h1*1.65, linetype="dashed", 
             color = "blue", size=0.4) +
  geom_hline(yintercept=mean_h3+sd_h3*1.65, linetype="dashed",
             color = "green", size=0.4) +
  geom_hline(yintercept=mean_h12+sd_h12*1.65, linetype="dashed",
             color = "red", size=0.4) +
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
# ggsave("macroUncertainty_index_plot.pdf")


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



###############################
## (5.2) this part reproduces what we did for the macro
## uncertainty measure for the financial uncertainty measure
## of Ludvigson et al (2018)
###############################
## we import the data (sheet 'Macro Uncertainty' only)
finUncertainty_index <- read_excel("Replication_data_Ludvigson.xlsx", 
                                   sheet = "Data")
# a quick check shows us that the columns are already correctly
# named, but the type for the 'Date'-column is set to num;
# the Date-column consists of a numeric holding 'yearmonth',
# hence we have to proceed as follows to change that:
finUncertainty_index$yearmon <- as.yearmon(
  as.character(finUncertainty_index$Date), "%Y%m")

# then we extract year and month and store them in two separate
# and distinct columns:
finUncertainty_index <- finUncertainty_index %>%
  mutate(year = year(yearmon),
         month = month(yearmon)) %>%
  dplyr::select(Date, yearmon, year, month, Uf)


# next we create the variable 'my' which is a numerical representation
# of yearmon:
finUncertainty_index <- as.data.frame(finUncertainty_index %>%
                                        mutate(my = year + month/12 - 1/12))

# we calculate
# the mean and standard deviation to be able to add horizontal
# lines marking 1.65 standard deviations above the mean to the plot
sd_Uf <- sd(finUncertainty_index$Uf)
mean_Uf <- mean(finUncertainty_index$Uf)

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
finUncertainty_index$thresh <- mean_Uf+sd_Uf*1.65

# the above variable now allows us to create an indicator of whether
# or not a certain point is above the threshold value
finUncertainty_index <- finUncertainty_index %>%
  mutate(fin_shock = case_when(Uf > thresh ~ 1,
                               Uf <= thresh ~ 0))

# having the variable fin_shock in our dataset now allows us to 
# add episodes into the graph that correspond to the
# macro-shock being 1;

# but before we can do that, we need to create 'start' and 'end' - dates for
# the respective periods (which we will then ultimately use for plotting
# the episodes with shaded regions!)

# to retrieve the respective 'start'- and 'end' - dates for the shock-periods,
# we apply a procedure to the dataset 'finUncertainty_index' that results in a 
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
as.data.frame(finUncertainty_index %>% filter(fin_shock == 1))


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
finUncertainty_index["shock_ID"] <- NA

# then we start the loop:
for (i in 1:nrow(finUncertainty_index)) {
  if(finUncertainty_index$fin_shock[i] == 1) {
    # if we detect a shock (i.e., 'bloom_shock' == 1), then we have 
    # to proceed as follows:
    # we store the 'my' variable in the column 'start' and 'end'
    finUncertainty_index[i, ncol(finUncertainty_index)] <- x
    
  }else{
    # no shock
    x <- x + 1
  }
}


# the above method is maybe not yet the most elegant/efficient solution
# but it suffices for now (we can fine-tune it at a later stage!)
# we can now have a look at the newly created variable
# (note that we have suppressed 'NA's for the creation of the table!):

finUncertainty_index %>%
  filter(!is.na(shock_ID))  %>%
  group_by(shock_ID) %>%
  summarise(Count = n())

# next we group_by shock_ID and extract the respective start-dates and
# store the retrieved data to a new data-frame
fin_shocks_start <- as.data.frame(finUncertainty_index %>%
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
fin_shocks_end <- as.data.frame(finUncertainty_index %>%
                                  filter(!is.na(shock_ID))  %>%
                                  group_by(shock_ID) %>%
                                  filter(row_number()==n()) %>%
                                  mutate(year_end=year, 
                                         month_end=month, 
                                         my_end = my) %>%
                                  dplyr::select(shock_ID, year_end, 
                                                month_end, my_end))

# next, we can merge the two data-frames from above:
fin_shocks_start_end <- merge(x = fin_shocks_start, 
                              y = fin_shocks_end, 
                              by = "shock_ID", 
                              all = TRUE, na.rm=T)
# and inspect the resulting data-frame:
fin_shocks_start_end

# we re-arrange the sequence of columns
fin_shocks_start_end <- fin_shocks_start_end %>%
  dplyr::select(shock_ID, 
                year_start, 
                month_start, 
                year_end, 
                month_end, 
                my_start, 
                my_end)


# next, we construct a year-month-variable both out of
# the pair 'year_start' & 'month_start' and 'year_end' & 'month_end'
fin_shocks_start_end$yearmon_start <- as.yearmon(
  fin_shocks_start_end$my_start, "%Y-%B")
fin_shocks_start_end$yearmon_end <- as.yearmon(
  fin_shocks_start_end$my_end, "%Y-%B")
# next, we add a helper-column that we need in the next stage as well:
fin_shocks_start_end$helper_date <- format(
  fin_shocks_start_end$yearmon_start, "%b")


# this means that we can now drop the 'year' and 'month' variables
fin_shocks_start_end <- fin_shocks_start_end %>%
  dplyr::select(-c(year_start, year_end, month_start, month_end))

# next we add a column that gives us the duration of the shock:
fin_shocks_start_end <- fin_shocks_start_end %>%
  mutate(duration = (fin_shocks_start_end$yearmon_end - 
                       fin_shocks_start_end$yearmon_start) 
         * 12 + 1) %>%
  dplyr::select(-helper_date)

# next we group_by shock_ID and extract the respective MAXIMUM
# VOLATILITY
# (note that the produced data-frame only has one row per shock_ID)
fin_shocks_max_vol <- as.data.frame(finUncertainty_index %>%
                                      filter(!is.na(shock_ID))  %>%
                                      group_by(shock_ID) %>%
                                      filter(Uf == max(Uf)) %>%
                                      mutate(max_Uf = Uf) %>%
                                      dplyr::select(shock_ID, max_Uf, my))

# we apply the same date-transformation as above:
fin_shocks_max_vol$yearmon_max <- as.yearmon(
  fin_shocks_max_vol$my, "%Y-%B")
# and we drop 'my':
fin_shocks_max_vol <- fin_shocks_max_vol %>%
  dplyr::select(-my)

# in a last stage, we merge fin_shocks_max_vol with 
# fin_shocks_start_end:
fin_shocks_start_end <- merge(x = fin_shocks_start_end, 
                              y = fin_shocks_max_vol, 
                              by = "shock_ID", all = TRUE, na.rm=T)


# actually, there is a slight problem, which we fix by adding
# one month's numeric value to my_end:
fin_shocks_start_end$my_end <- fin_shocks_start_end$my_end + (1/12)
# this makes sure that if we refer to an end-month that
# we assume that we are talking about the last day of that
# respective month!


# the above data-frame now allows us to plot the episodes of 
# high volatility into our previous plots (by adding a shaded
# rectangle!)
# as well as the NBER recession dates

## and finally plot the series
finUncertainty_index_plot <- ggplot() +
  geom_line(data = finUncertainty_index, 
            aes(x = my, 
                y = Uf), color="red", 
            size=1,linetype = 1) +
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Fin Uncertainty Indices", 
                     limits = c(0.5, 1.6), 
                     breaks = seq(0.5, 1.5, by = 0.2), 
                     minor_breaks = NULL) +
  scale_color_discrete(name = NULL, labels = c("Uf")) +
  geom_hline(yintercept=mean_Uf+sd_Uf*1.65, linetype="dashed",
             color = "blue", size=0.4) +
  theme(legend.position = "bottom", axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14)) + 
  # add in fin - shocks (similar to Bloom's approach)
  # geom_rect(data=fin_shocks_start_end, inherit.aes = FALSE,
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

finUncertainty_index_plot



# and we save the plot to be used in our latex-document
# ggsave("macroUncertainty_index_plot.pdf")


# having created the data-frame macro_shocks_start_end, we can now also
# export a latex-table with the exact data of the macro-shock!
# first, we slightly adapt the shocks - data frame:
# first, we add '%' to 'max_vol' in our data-frame:
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x, format = format, digits = digits, ...), "%")
}

fin_shocks_table <- fin_shocks_start_end %>%
  mutate(max_Uf = percent(max_Uf), 
         duration = paste(as.character(floor(duration)), 
                          "months", sep=" ")) %>%
  dplyr::select(yearmon_start, 
                yearmon_end, 
                duration, 
                yearmon_max, 
                max_Uf)

# adapt the below command to supply a vector of column-headers!
latextable(fin_shocks_table, colnames = c("Start", 
                                          "End", "Duration", 
                                          "Maximum Year/Month", 
                                          "Maximum Value"), dp=3)