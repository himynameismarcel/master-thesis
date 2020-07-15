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
                EPU_Hist = EPU_Historical,
                Macro = macroUncert,
                FinUnc = Uf)


# we reorder the variables
# and ultimately decided to drop the GTU and Michigan and the 
# News_Based_Policy_Uncert_Index
#comparison_measures <- as_tibble(comparison_measures[c(2, 1, 3, 4, 5, 7)])
#and we added the financial uncertainty measure from jurado et al
comparison_measures_plot <- as_tibble(comparison_measures[c(2, 3, 5, 4)])



# we make a time-series object out of our data-frame
# (and exclude the simple 'EPU')
# Marcel (17.05.2020): Note that we do not include the
# financial uncertainty measure of Ludvigson et al (2018)
# at this stage (we only include it for the production of
# the summary statistics)

# Marcel (05.07.2020):
# because using the previous method of 'autoplot' does not
# allow enough flexibility in plotting, we decided to go 
# back to the old-school approach of plotting.
# hence the below code is commented out and instead we 
# follow the same procedure that we use in the script
# 'Identification_of_Shocks_LMN_v0.3.R':

          # comparison_measures_ts <- ts(comparison_measures[c(2, 3, 5)],
          #                              start = c(1960, 7), frequency = 12)
          # 
          # 
          # # we plot all the time-series by making use of
          # # of the forecast-package (that reverts back to ggplot2)
          # # and add in the NBER-dates;
          # 
          # comparison_plot <- autoplot(comparison_measures_ts, facets = TRUE) + 
          #   scale_x_continuous(name = "Year", limits = c(1960, 2020), 
          #                      breaks = seq(1960, 2020, by = 5),
          #                      minor_breaks = NULL) + 
          #   scale_y_continuous(name = "")  +
          #   labs(color=NULL) +
          #   theme(axis.text=element_text(size=14),
          #         #axis.title=element_text(size=15,face="bold"),
          #         legend.text=element_text(size=14), 
          #         panel.grid.major.x = element_blank()) +
          #   # note that panel.grid.major.x = element_blank() 
          #   # suppresses vertical grid lines!
          #   geom_rect(data=recessions_start_end, 
          #             inherit.aes = FALSE,
          #             aes(xmin=my_start, xmax=my_end, 
          #                 ymin=-Inf, ymax=+Inf), 
          #             fill='#606060', alpha=0.5)
          # 
          # comparison_plot

          # ggsave("comparison_plot.pdf")


# with the tibble 'comparison_measures_plot' containing my, VXO,
# Macro and EPU Hist, we perform the following gather_statement:
comparison_measures_plot <- comparison_measures_plot %>%
                      gather(
                        uncert_measure, value,
                        -my, na.rm = TRUE
                      )

# new facets labels:
supp.labs <- c(`Macro`= "Macro Uncertainty Um", 
               `VXO`= "VXO/stock market vol.",
               `EPU_Hist`= "Historical EPU")

# now we should be able to plot the result:
comparison_measures.plot <- ggplot(comparison_measures_plot,
                                   aes(x=my,
                                       y=value,
                                       color=uncert_measure)) +
  geom_rect(data=recessions_start_end[16:23, ], inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=Inf), 
            fill='#606060', alpha=0.3) + 
  geom_line(size = 0.8) + 
  facet_grid(factor(uncert_measure, levels = c("VXO","Macro", "EPU_Hist")) ~ .,  
             labeller = as_labeller(supp.labs), scales="free_y") +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "", 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  scale_color_manual(values = c("Macro" = "#0586ff",
                                "VXO" = "#0fa31e",
                                "EPU_Hist" = "#7a05ff")) +
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 12, colour = "black")) + 
  labs(color=NULL) +
  theme(legend.position="none",
        axis.text=element_text(size=15),
        #plot.title = element_text(size=10, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


comparison_measures.plot

ggsave("comparison_plot.pdf")



# because autoplot does not quite do what we want and we cannot
# apply the 'gather()' - function to a time-series object,
# we apply the scaling also to 'comparison_measures' (without '_ts')
# so that we can use a non-time-series data.frame in our plot below:
# (but this time we perform the scaling manually!)
scaled_comparison_measures <- as_tibble(comparison_measures[c(2, 3, 5, 4)]) %>%
                dplyr::mutate(VXO = (VXO - mean(VXO, 
                                  na.rm=TRUE))/
                                sd(VXO, na.rm=TRUE),
                EPU_Hist  = (EPU_Hist - mean(EPU_Hist, 
                                  na.rm=TRUE))/
                                sd(EPU_Hist, na.rm=TRUE),
                # Michigan  = (Michigan - mean(Michigan, 
                #                    na.rm=TRUE))/
                #   sd(Michigan, na.rm=TRUE),
                Macro  = (Macro - mean(Macro, 
                                  na.rm=TRUE))/
                                sd(Macro, na.rm=TRUE))    


# we bring the data into a tidy data-format:
scaled_comparison_measures <- scaled_comparison_measures %>%
          gather(uncert_measure, 
         value, -my, na.rm = TRUE)

# finally we can generate the plot with the standardized series:
comparison_plot_combined <- ggplot(scaled_comparison_measures, 
       aes(x = my, y = value, colour=uncert_measure,
                            shape=uncert_measure)) +
  geom_rect(data=recessions_start_end, inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=+Inf), 
            fill='#606060', alpha=0.5) + 
  geom_line(size=1.1) +
  # geom_point(size=1) + 
  scale_colour_manual(values=c("#7a05ff",
                              "#0586ff",
                              "#0fa31e"),
                     labels=c("Historical EPU",
                              "Macro Uncertainty Um",
                              "VXO/stock market vol.")) + 
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "", limits = c(1960, 2018), 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  #labs(color=NULL) +
  theme(legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 15)
        #legend.text=element_text(size=14)
        ) + 
  # note that panel.grid.major.x = element_blank() suppresses vertical grid lines!
  scale_linetype_manual(values = c(1,1,1,1)) +
  scale_shape_manual(values=c(0,1,2,3)) + 
  guides(colour = guide_legend(override.aes = list(size = 5),
                               reverse=TRUE))


comparison_plot_combined

# ggsave("comparison_plot_combined.pdf")
