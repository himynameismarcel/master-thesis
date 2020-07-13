#########################################################################################
### Marcel Kropp, 15.04.2018
### This script is separated into several big parts:
### PART 6: deals with the plotting of all our analyzed time-series into
###         one plot (to make them comparable!)
###         the series which we want to plot together at once are
###           *) VIX (used by Bloom)
###           *) Michigan Survey
###           *) EPU (Baker et al, 2016)
###           *) GT and GTU (Bontempi et al. 2015/Castelnuovo and Tran, 2017)
###           *) Macro Uncertainty Index (Jurado et al, 2015)
########################################################################



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
if(!exists("foo", mode="function")) source("20200705_sp500_vxo_bloom_shocks_v0.1.R")
if(!exists("foo", mode="function")) source("20200705_macro_fin_uncert_v0.1.R")
if(!exists("foo", mode="function")) source("20200705_EPU_baker_et_al_v0.1.R")

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

# Marcel (17.05.2020):
# Note: Because we have switched the variable 'my' by 1/12 for both
# the macro uncertainty series and financial uncertainty series 
# from Ludvigson et al, we have to switch them back for the below;
# when we go through this code again, we ahve to change everything
# accordingly so that all series are consistent between each other!

macroUncertainty_index <- macroUncertainty_index %>%
                                           mutate(my = my + 1/12)
finUncertainty_index <- finUncertainty_index %>%
                                         mutate(my = my + 1/12)

comparison_measures <- data.frame("macroUncert" = macroUncertainty_index$h1,
                                  "my" = macroUncertainty_index$my)

# then we join in the data one after the other, starting with 
# the data from the GTU_index,
# followed by the original epu_index,
# as well as the historical EPU_index,
# then the VXO (the cycle-part after detrending),
# and finally the MSoC in a multi-setp join:

comparison_measures <- right_join(x = epu_index_historical[, 
                                      c("EPU_Historical","my")],
                                 y = comparison_measures,
                                 by = "my") %>%
                      # here we switched 'mvol_cycle' to
                      # 'mvol' - we had wrongly
                      # used 'mvol_cycle'
                      right_join(x = sp500_merge_vxo[, 
                                 c("mvol_VOLATBL","my")],
                                 # c("mvol","my")],
                                 y = .,
                                 by = "my") %>%
                      # Marcel (17.05.2020):
                      # below we are adding the financial uncertainty
                      # measure from Ludvigson et al:
                      right_join(x = finUncertainty_index[, 
                                      c("Uf","my")],
                                       y = .,
                                       by = "my")

# Note that the resulting data.frame ranges until the end of 2017 due to
# the macro uncertainty index having exactly that range!
# (i.e., all other data.frames that we read in afterwards, only run
# until the end of 2017 as well!)

# Rename all the variables:
comparison_measures <- comparison_measures %>%
  dplyr::rename(VXO = mvol_VOLATBL,
                'EPU Hist' = EPU_Historical,
                Macro = macroUncert,
                FinUnc = Uf)


# we reorder the varaibles
# and ultimately decided to drop the GTU and Michigan and the 
# News_Based_Policy_Uncert_Index
#comparison_measures <- as_tibble(comparison_measures[c(2, 1, 3, 4, 5, 7)])
#and we added the financial uncertainty measure from jurado et al
comparison_measures <- as_tibble(comparison_measures[c(2, 3, 5, 1, 4)])



# we make a time-series object out of our data-frame
# (and exclude the simple 'EPU')
# Marcel (17.05.2020): Note that we do not include the
# financial uncertainty measure of Ludvigson et al (2018)
# at this stage (we only include it for the production of
# the summary statistics)
comparison_measures_ts <- ts(comparison_measures[c(2, 3, 5)],
                             start = c(1960, 7), frequency = 12)


# we plot all the time-series by making use of
# of the forecast-package (that reverts back to ggplot2)
# and add in the NBER-dates;

comparison_plot <- autoplot(comparison_measures_ts, facets = TRUE) + 
  scale_x_continuous(name = "Year", limits = c(1960, 2020), 
                     breaks = seq(1960, 2020, by = 5),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "")  +
  labs(color=NULL) +
  theme(axis.text=element_text(size=14),
        #axis.title=element_text(size=15,face="bold"),
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

# ggsave("comparison_plot.pdf")


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
                # Michigan  = (Michigan - mean(Michigan, 
                #                    na.rm=TRUE))/
                #   sd(Michigan, na.rm=TRUE),
                Macro  = (Macro - mean(Macro, 
                                       na.rm=TRUE))/
                  sd(Macro, na.rm=TRUE))    
# %>%
# and we remove EPU:
# dplyr::select(-EPU)

# we bring the data into a tidy data-format:
scaled_comparison_measures <- scaled_comparison_measures %>%
  # note that the 'gather()' - function produces a warning
  # message because 'Michigan' enter as a 'zoo' - object!
  gather(uncert_measure, value, -my, na.rm = TRUE)

ggplot(scaled_comparison_measures, 
       aes(x = my, y = value, colour=uncert_measure,
           linetype = uncert_measure, shape=uncert_measure)) +
  geom_line(size=0.8) +
  # geom_point(size=1) + 
  scale_color_manual(values=c("#01D71E",
                              "#FF4DCA",
                              "#0000FF")) + 
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

# ggsave("comparison_plot_combined.pdf")
