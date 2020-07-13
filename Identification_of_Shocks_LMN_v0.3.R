# This script generates Figure 1 of Ludvigson et al 2018:

# The approach is a combination of our very first script
# (called '12052018_thesis_master v0.2.R') and the 
# script 'Identification_of_Shocks.R' which was originally
# implemented for the identification of shocks for all the
# uncertainty series in out data!


###############################
## loading libraries
###############################
if(!exists("foo", mode="function")) source("20200705_load_libraries_v0.1.R")

###############################
## reading in NBER recession dates
###############################
if(!exists("foo", mode="function")) source("20200705_nber_recession_dates_v0.1.R")

###############################
## reading in the script 20200628_functions_v0.2.R
###############################
if(!exists("foo", mode="function")) source("20200628_functions_v0.2.R")



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

#------------------------------------------------------------------------
# 20200705, Marcel:
# Note that we have moved the function to our script called 
# '20200628_functions.R'.
# Accordingly we have to source this file 20200628_functions.R'
# on top of this script!
#------------------------------------------------------------------------

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
# actually started in Dec 2015. 
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
               `Uf`= "Aggregate Financial Uncertainty Uf")


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
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "", 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  scale_color_manual(values = c("Um" = "blue", 
                                "Uf" = "red")) + 
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 12, colour = "black")) + 
  labs(color=NULL) +
  geom_hline(yintercept=1.65, linetype="dashed", color = "black") + 
  theme(legend.position="none",
        axis.text=element_text(size=15),
        #plot.title = element_text(size=10, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  geom_text(aes(x=2015, y=2, label = "1.65 std"))


LMN_Shocks_plot_combined

ggsave("LMN_Shocks_plot_combined.pdf")

## -------------------------------------------------------------------
# ultimately, we also want to add the NBER recession dates to the
# plot (which we have also done in 'Identification_of_Shocks.R')
## -------------------------------------------------------------------
## see above! we have taken the NBER-recession dates from another 
## script (where we have prepared them);


## calculation of correlation:
cor(SVAR.data$Uf, SVAR.data$Um)
















