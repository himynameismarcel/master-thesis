#########################################################################################
### Marcel Kropp, 21.05.2018
### This script estimates VARs following the original contribution of
### Bloom (2009) in his paper entitled 'The Impact of Uncertainty Shocks' and
### corresponds to VAR8 mentioned in the main text of the thesis;
### Alternative VAR8-estimations without detrending will be estimated in another
### R-Script (following Jurado et al., 2015, based on Bloom, 2009) and will be 
### parked in the Appendix (and hence also in a another dedicated R-Script).

### The data we use is an Excel-file which Nicholas Bloom himself makes available
### via his homepage as part of the available material dedicated to his
### 2009 paper;
### This Excel-file can be downloaded at 
### https://people.stanford.edu/nbloom/sites/default/files/r.zip and 
### comes along bundled together in a zip-file with an R-Script prepared by
### Jesus Fernandez-Villaverde (JFV) that replicates Bloom's (2009) VARs using
### slightly updated data and a longer time-range (the differences between
### Bloom's original used data and JFV's updated data are negligible and
### are mostly limited to the data for Industrial Production having slightly
### changed due to a change in the base-year between Bloom's contrubtion and
### JFV's replication);
### To stick to the original contribution from Bloom (2009), we will use
### JFV's updated data but restrict the range of the data to replicate
### Bloom (2009)'s results using 5 different uncertainty measures:
###         (PART 1) the Bloom-Shocks (marking months with maximum volatility)
###         (PART 2) the VXO/volatility series
###         (PART 3) the MSoC (available as of 1978)
###         (PART 4) the EPU (availbale as of 1985)
###         (PART 5) the GTU (available as of 2004)
###         (PART 6) the Macro Uncertainty Index, for h=1 and h=12 
###         (available as of 1960)
### After restricting the range of the data, the various time-series will
### run from 07/1962-06/2008 (like in Bloom's original contribution from 2009).

### The execution of the below relies on data and objects from 
### '12052018_thesis_master.R'; therefore this R-Script should be run 
### before running this one, to have all necessary material available

### Note: Because each of the below parts uses a different uncertainty
###       measure for the VARs, each part is self-contained in the sense
###       that it starts again from scratch with the readining-in of the
###       raw data 'VAR_UPDATED_JFV.xlsx'!


#########################################
### preparatory steps common to all VARs
#########################################
#####
## Read xlsx data and some basic manipulations
#####

var8_bloom_JFV <- read_excel("VARDATA_UPDATED_JFV.xlsx")

# Drop last observations in 2013 where no volatbl-data
# (i.e., VXO/volatility) is available using 'drop_na', 
# drop all rows that go beyond Bloom's original estimations
# (from 2009) by uniting YEAR and MONTH into the helper variable
# 'YEARMONTH',
# change the type to 'yearmon' to be able to filter accordingly,
# filter as we wish (i.e., YEARMONTH <= "Jun 2008")
# and at the end remove the helper-variable
# 'YEARMONTH' again:
var8_bloom_JFV <- var8_bloom_JFV %>% 
                        drop_na %>%
                        unite(YEARMONTH, YEAR, MONTH, remove = FALSE,
                              sep="-") %>%
                        mutate(YEARMONTH = as.yearmon(YEARMONTH)) %>%
                        filter(YEARMONTH <= "Jun 2008") %>%
                        dplyr::select(-YEARMONTH)

#####
## log- and/or HP-transformations
#####
# because we want to apply the same transformation to a couple of
# variable, we decide for the dplyr-approach:

# log-and HP-transformations
cols_log_hp <- c("STOCK", "WAGE", "CPI", "EMPM", "IPM")
for (i in seq_along(cols_log_hp)) {
  var8_bloom_JFV[cols_log_hp[i]] <- log(var8_bloom_JFV[cols_log_hp[i]])
  var_aux <- hpfilter(var8_bloom_JFV[cols_log_hp[i]],freq=129600)
  var8_bloom_JFV[cols_log_hp[i]] <- var_aux$cycle
}

# HP-transformations
cols_hp <- c("HOURSM", "FFR")
for (i in seq_along(cols_log)) {
  var_aux <- hpfilter(var8_bloom_JFV[cols_hp[i]],freq=129600)
  var8_bloom_JFV[cols_hp[i]] <- var_aux$cycle
}


########################################################################
### PART 1: VAR8_Bloom_2009 with Bloom-Shocks (marking months
###         with maximum volatility)
########################################################################
# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_1 <- var8_bloom_JFV

####################
## transformations
####################
# Bloom_Shock: we want to define an indicator function on the 17 events defined
# in the original contribution by Bloom (2009)

# note: all below steps could be executed in ONE COMMAND (correct this!)
# first, we change the name of the column
var8_bloom_JFV_1 <- var8_bloom_JFV_1 %>%
                        dplyr::rename(bloom_shock = VOLATBL)

# and then set the indicator variable to 1 for months were the
# volatility had reached its peak (according to Table A.1, Bloom, 2009, p. 35)
var8_bloom_JFV_1 <- var8_bloom_JFV_1 %>%
            unite(YEARMONTH, YEAR, MONTH, remove = FALSE,
                  sep="-") %>%
            mutate(YEARMONTH = as.yearmon(YEARMONTH)) %>%
            mutate(bloom_shock = case_when(
              YEARMONTH %in% 
                as.yearmon(c("Oct 1962", #Cuban missile crisis
                              "Nov 1963", #Assassination of JFK
                              "Aug 1966", #Vietnam buildup
                              "Dec 1973", #OPEC I, Arab-Israeli War
                              "Oct 1974", #Franklin National Financial Crisis
                              "Nov 1978", #OPEC II
                              "Mar 1980", #Afghanistan, Iran hostages
                              "Oct 1982", #Monetary cycle turning point
                              "Nov 1987", #Black Monday
                              "Oct 1990", #Gulf War I
                              "Nov 1997", #Asian Crisis
                              "Sep 1998", #Russian, LTCM default
                              "Sep 2001", #9/11 terrorist attack
                              "Sep 2002", #Worldcom and Enron
                              "Feb 2003", #Gulf War II
                              "Mar 2008")) ~ 1, #Credit Crunch
                  TRUE                    ~ 0)) %>%
            # in a last step, we remove all date variables
            dplyr::select(-c(YEARMONTH, YEAR, MONTH))
  

# Check data
View(var8_bloom_JFV_1)


# we reorder variables (and only select the ones we
# want for the VARs) and can finally kick off the 
# computation of the VARs (and IRFs)
var8_bloom_JFV_1 <- var8_bloom_JFV_1 %>%
                        dplyr::select(STOCK, bloom_shock, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)


####################
## VARs (VAR8)
####################
# estimate model (12 lags + constant)
var8_bloom_JFV_1.var <- VAR(var8_bloom_JFV_1, p=12, type="const")


# next, compute the irfs: we extract the responses for "IPM", "EMPM",
# and "bloom_shock" at once
var8_bloom_JFV_1.var.irfs <- irf(var8_bloom_JFV_1.var, 
                                response = c("IPM", "EMPM", "bloom_shock"), 
                                impulse = "bloom_shock", n.ahead = 36, 
                                boot = TRUE, ortho = TRUE, runs=100)

# Note: all results are now stored in var8_bloom_JFV_1.var.irfs
# running the below reveals the nested list-structure of the
# 'varrirf' - object:
str(var8_bloom_JFV_1.var.irfs)

# first preliminary plots
# matplot(var8_bloom_JFV_1.var.irfs$irf$bloom_shock[, 3], type='l') 
# lines(var8_bloom_JFV_1.var.irfs$Lower$bloom_shock[, 3], col="blue")     
# lines(var8_bloom_JFV_1.var.irfs$Upper$bloom_shock[, 3], col="blue") 

# to normalize the orthogonalized irfs, we rescale them by dividing
# all values for response = EMPM/IPM by the first entry of
# bloom_shock on itself;
# to do that, we first store the results in a data-frame/tibble
var8_bloom_JFV_1.var.irfs.df <- data.frame(var8_bloom_JFV_1.var.irfs$irf,
                                           var8_bloom_JFV_1.var.irfs$Lower,
                                           var8_bloom_JFV_1.var.irfs$Upper)

# and rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for bloom_shock.bloom_shock
var8_bloom_JFV_1.var.irfs.df <- var8_bloom_JFV_1.var.irfs.df %>%
                dplyr::rename(blo.blo.oirf = bloom_shock.bloom_shock, 
                              blo.EMPM.oirf = bloom_shock.EMPM,
                              blo.IPM.oirf = bloom_shock.IPM,
                              blo.EMPM.Lo = bloom_shock.EMPM.1,
                              blo.IPM.Lo = bloom_shock.IPM.1,
                              blo.EMPM.Up = bloom_shock.EMPM.2,
                              blo.IPM.Up = bloom_shock.IPM.2) %>%
                dplyr::select(-c(bloom_shock.bloom_shock.1, 
                                 bloom_shock.bloom_shock.2))

# finally, to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of blo.blo.oirf (for n=1),
# and at the same time multiply with 100:

var8_bloom_JFV_1.var.irfs.df.rescaled <- 
            var8_bloom_JFV_1.var.irfs.df[, 2:ncol(var8_bloom_JFV_1.var.irfs.df)] /
            var8_bloom_JFV_1.var.irfs.df[1, 1]*100

# and lastly, we add a step-counter to the df:
var8_bloom_JFV_1.var.irfs.df.rescaled$step <- 
            seq.int(nrow(var8_bloom_JFV_1.var.irfs.df.rescaled))

# and split our data-set in two separate ones (which is actually not
# necessary anymore because we generate the plots one by one):
var8_bloom_JFV_1.ipm.plot <- var8_bloom_JFV_1.var.irfs.df.rescaled %>%
                        dplyr::select(contains("IPM"))
var8_bloom_JFV_1.empm.plot <- var8_bloom_JFV_1.var.irfs.df.rescaled %>%
                        dplyr::select(contains("EMPM"))


####################
## plots of orifs
####################
# now we can plot both the oirfs (rescaled) for IPM (next) and EMPM (below!):
irf_blo_shock_ipm <- ggplot(data = var8_bloom_JFV_1.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = blo.IPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = blo.IPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=blo.IPM.Lo, ymin=blo.IPM.Up), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = blo.IPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = blo.IPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Bloom-Shock", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("% impact on industrial production") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# oirf for EMPM:
irf_blo_shock_empm <- ggplot(data = var8_bloom_JFV_1.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = blo.EMPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = blo.EMPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=blo.EMPM.Up, ymin=blo.EMPM.Lo), 
            fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = blo.EMPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = blo.EMPM.oirf), color="black", 
            size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("% impact on employment in manufacturing") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# print to device
grid.arrange(irf_blo_shock_ipm, irf_blo_shock_empm, ncol=2)


########################################################################
### PART 2: VAR8_Bloom_2009 with actual VXO/volatility series
########################################################################
# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_2 <- var8_bloom_JFV

####################
## transformations
####################
# we only select the variables we
# want for the VARs and can finally kick off the 
# computation of the VARs (and IRFs)
var8_bloom_JFV_2 <- var8_bloom_JFV_2 %>%
                        dplyr::select(STOCK, VOLATBL, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)


####################
## VARs (VAR8)
####################
# estimate model
var8_bloom_JFV_2.var <- VAR(var8_bloom_JFV_2, p=12, type="const")

# next, compute the irfs: we extract the responses for "IPM", "EMPM",
# and "bloom_shock" at once
var8_bloom_JFV_2.var.irfs <- irf(var8_bloom_JFV_2.var, 
                                 response = c("IPM", "EMPM", "VOLATBL"), 
                                 impulse = "VOLATBL", n.ahead = 36, 
                                 boot = TRUE, ortho = TRUE, runs=100)

# Note: all results are now stored in var8_bloom_JFV_2.var.irfs
# running the below reveals the nested list-structure of the
# 'varrirf' - object:
str(var8_bloom_JFV_2.var.irfs)

# first preliminary plots
# matplot(var8_bloom_JFV_2.var.irfs$irf$VOLATBL[, 3], type='l') 
# lines(var8_bloom_JFV_2.var.irfs$Lower$VOLATBL[, 3], col="blue")     
# lines(var8_bloom_JFV_2.var.irfs$Upper$VOLATBL[, 3], col="blue") 

# to normalize the orthogonalized irfs, we rescale them by dividing
# all values for response = EMPM/IPM by the first entry of
# VOLATBL on itself;
# to do that, we first store the results in a data-frame/tibble
var8_bloom_JFV_2.var.irfs.df <- data.frame(var8_bloom_JFV_2.var.irfs$irf,
                                           var8_bloom_JFV_2.var.irfs$Lower,
                                           var8_bloom_JFV_2.var.irfs$Upper)

# and rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for bloom_shock.bloom_shock
var8_bloom_JFV_2.var.irfs.df <- var8_bloom_JFV_2.var.irfs.df %>%
  dplyr::rename(VOL.VOL.oirf = VOLATBL.VOLATBL, 
                VOL.EMPM.oirf = VOLATBL.EMPM,
                VOL.IPM.oirf = VOLATBL.IPM,
                VOL.EMPM.Lo = VOLATBL.EMPM.1,
                VOL.IPM.Lo = VOLATBL.IPM.1,
                VOL.EMPM.Up = VOLATBL.EMPM.2,
                VOL.IPM.Up = VOLATBL.IPM.2) %>%
  dplyr::select(-c(VOLATBL.VOLATBL.1, 
                   VOLATBL.VOLATBL.2))


# to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of VOL.VOL.oirf (for n=1),
# and at the same time multiply with 100;
# as part of the normalization of the oirfs, we have to choose 
# meaningful units; for vol we choose a 15 unit shock
# which amounts to approx four standard deviations of the
# identified error;

var8_bloom_JFV_2.var.irfs.df.rescaled <- 
  var8_bloom_JFV_2.var.irfs.df[, 2:ncol(var8_bloom_JFV_2.var.irfs.df)]*15 /
  var8_bloom_JFV_2.var.irfs.df[1, 1]*100

# and lastly, we add a step-counter to the df:
var8_bloom_JFV_2.var.irfs.df.rescaled$step <- 
  seq.int(nrow(var8_bloom_JFV_2.var.irfs.df.rescaled))

# and split our data-set in two separate ones:
var8_bloom_JFV_2.ipm.plot <- var8_bloom_JFV_2.var.irfs.df.rescaled %>%
  dplyr::select(contains("IPM"))
var8_bloom_JFV_2.empm.plot <- var8_bloom_JFV_2.var.irfs.df.rescaled %>%
  dplyr::select(contains("EMPM"))


####################
## plots of orifs
####################
# now we can plot both the oirfs (rescaled) for IPM (next) and EMPM (below!):
irf_VOL_shock_ipm <- ggplot(data = var8_bloom_JFV_2.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = VOL.IPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = VOL.IPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=VOL.IPM.Lo, ymin=VOL.IPM.Up), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = VOL.IPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = VOL.IPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "VXO/volatility shock", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# oirf for EMPM:
irf_VOL_shock_empm <- ggplot(data = var8_bloom_JFV_2.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = VOL.EMPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = VOL.EMPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=VOL.EMPM.Up, ymin=VOL.EMPM.Lo), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = VOL.EMPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = VOL.EMPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# print to device
grid.arrange(irf_blo_shock_ipm, irf_blo_shock_empm, 
             irf_VOL_shock_ipm, irf_VOL_shock_empm, ncol=2)


########################################################################
### PART 6: VAR8_Bloom_2009 with Macro Uncertainty Index
###         following Jurado et al. (2015);
###         We specifically use the data for h=1 and h=12 to see how
###         they perform compared to each other
########################################################################

# because we do not want to 'destroy' the data-frame 'var8_bloom_JFV'
# which we have set up above, we create a new data-frame on which we
# perform all operations related to PART 1

var8_bloom_JFV_6 <- var8_bloom_JFV

####################
## transformations
####################
# Merge the data with the data.frame containing the data for the
# Macro Uncertainty Index (for both h=1 and h=12);
# the respective data comes from the data-frame: macroUncertainty_index;
# first, we have to rename 'YEAR' and 'MONTH' to 'year' and 'month':
var8_bloom_JFV_6 <- var8_bloom_JFV_6 %>%
                    dplyr::rename(year = YEAR,
                                 month = MONTH)
# next, we perform a right-join to have the data for h1 and h12 in our
# data.frame, ranging from 1962-2008 (the same duration as the Bloom, 2009,
# estimations!)
var8_bloom_JFV_6 <- right_join(x = macroUncertainty_index[, c("h1", "h12", "year",
                                                              "month")], 
                         y = var8_bloom_JFV_6, 
                         by = c("year", "month"))


# we reorder variables (and only select the ones we
# want for the VARs)
var8_bloom_JFV_6_h1 <- var8_bloom_JFV_6 %>%
                        dplyr::select(STOCK, h1, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)
var8_bloom_JFV_6_h12 <- var8_bloom_JFV_6 %>%
                        dplyr::select(STOCK, h12, FFR, WAGE,
                                      CPI, HOURSM, EMPM, IPM)

####################
## VARs (VAR8) (note that we perform every step for h=1 and h=12)
####################
# estimate models
var8_bloom_JFV_6_h1.var <- VAR(var8_bloom_JFV_6_h1, p=12, type="const")
var8_bloom_JFV_6_h12.var <- VAR(var8_bloom_JFV_6_h12, p=12, type="const")

# next, compute the irfs: we extract the responses for "IPM", "EMPM",
# and "bloom_shock" at once
var8_bloom_JFV_6_h1.var.irfs <- irf(var8_bloom_JFV_6_h1.var, 
                                 response = c("IPM", "EMPM", "h1"), 
                                 impulse = "h1", n.ahead = 36, 
                                 boot = TRUE, ortho = TRUE, runs=100)
var8_bloom_JFV_6_h12.var.irfs <- irf(var8_bloom_JFV_6_h12.var, 
                                    response = c("IPM", "EMPM", "h12"), 
                                    impulse = "h12", n.ahead = 36, 
                                    boot = TRUE, ortho = TRUE, runs=100)

# Note: all results are now stored in var8_bloom_JFV_6_h1.var.irfs and
# var8_bloom_JFV_6_h12.var.irfs, respectively.
# running the below reveals the nested list-structure of the
# 'varrirf' - object:
str(var8_bloom_JFV_6_h1.var.irfs)
str(var8_bloom_JFV_6_h12.var.irfs)

# first preliminary plots
# matplot(var8_bloom_JFV_6_h12.var.irfs$irf$h12[, 3], type='l') 
# lines(var8_bloom_JFV_6_h12.var.irfs$Lower$h12[, 3], col="blue")     
# lines(var8_bloom_JFV_6_h12.var.irfs$Upper$h12[, 3], col="blue") 

# to normalize the orthogonalized irfs, we rescale them by dividing
# all values for response = EMPM/IPM by the first entry of
# bloom_shock on itself;
# to do that, we first store the results in a data-frame/tibble
var8_bloom_JFV_6_h1.var.irfs.df <- data.frame(var8_bloom_JFV_6_h1.var.irfs$irf,
                                              var8_bloom_JFV_6_h1.var.irfs$Lower,
                                              var8_bloom_JFV_6_h1.var.irfs$Upper)
var8_bloom_JFV_6_h12.var.irfs.df <- data.frame(var8_bloom_JFV_6_h12.var.irfs$irf,
                                              var8_bloom_JFV_6_h12.var.irfs$Lower,
                                              var8_bloom_JFV_6_h12.var.irfs$Upper)


# and rename the columns with the values for the Upper- and Lower
# CI to sensible names and at the same time delete
# the Upper and Lower CI bands for bloom_shock.bloom_shock
var8_bloom_JFV_6_h1.var.irfs.df <- var8_bloom_JFV_6_h1.var.irfs.df %>%
  dplyr::rename(h1.h1.oirf = h1.h1, 
                h1.EMPM.oirf = h1.EMPM,
                h1.IPM.oirf = h1.IPM,
                h1.EMPM.Lo = h1.EMPM.1,
                h1.IPM.Lo = h1.IPM.1,
                h1.EMPM.Up = h1.EMPM.2,
                h1.IPM.Up = h1.IPM.2) %>%
  dplyr::select(-c(h1.h1.1, 
                   h1.h1.2))
var8_bloom_JFV_6_h12.var.irfs.df <- var8_bloom_JFV_6_h12.var.irfs.df %>%
  dplyr::rename(h12.h12.oirf = h12.h12, 
                h12.EMPM.oirf = h12.EMPM,
                h12.IPM.oirf = h12.IPM,
                h12.EMPM.Lo = h12.EMPM.1,
                h12.IPM.Lo = h12.IPM.1,
                h12.EMPM.Up = h12.EMPM.2,
                h12.IPM.Up = h12.IPM.2) %>%
  dplyr::select(-c(h12.h12.1, 
                   h12.h12.2))

# to rescale, we divide all columns in the data-rame (apart from the
# first one) by the first value of h1.h1.oirf (for n=1),
# and at the same time multiply with 100;
# as part of the normalization of the oirfs, we have to choose 
# meaningful units; for vol we choose a 15 unit shock
# which amounts to approx four standard deviations of the
# identified error;

var8_bloom_JFV_6_h1.var.irfs.df.rescaled <- 
  var8_bloom_JFV_6_h1.var.irfs.df[, 2:ncol(var8_bloom_JFV_6_h1.var.irfs.df)]*220/
  (var8_bloom_JFV_6_h1.var.irfs.df[1, 1]*100)
var8_bloom_JFV_6_h12.var.irfs.df.rescaled <- 
  var8_bloom_JFV_6_h12.var.irfs.df[, 2:ncol(var8_bloom_JFV_6_h12.var.irfs.df)]*220/
  (var8_bloom_JFV_6_h12.var.irfs.df[1, 1]*100)

# and lastly, we add a step-counter to the df:
var8_bloom_JFV_6_h1.var.irfs.df.rescaled$step <- 
  seq.int(nrow(var8_bloom_JFV_6_h1.var.irfs.df.rescaled))
var8_bloom_JFV_6_h12.var.irfs.df.rescaled$step <- 
  seq.int(nrow(var8_bloom_JFV_6_h12.var.irfs.df.rescaled))

# and split our data-set in two separate ones (which is actually not necessary
# anymore; because we decided to generate the plots step-by-step):
var8_bloom_JFV_6_h1.ipm.plot <- var8_bloom_JFV_6_h1.var.irfs.df.rescaled %>%
  dplyr::select(contains("IPM"))
var8_bloom_JFV_6_h1.empm.plot <- var8_bloom_JFV_6_h1.var.irfs.df.rescaled %>%
  dplyr::select(contains("EMPM"))


####################
## plots of orifs
####################
# now we can plot both the oirfs (rescaled) for IPM (next) and EMPM (below!):
# for h1:
irf_h1_shock_ipm <- ggplot(data = var8_bloom_JFV_6_h1.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = h1.IPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = h1.IPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=h1.IPM.Lo, ymin=h1.IPM.Up), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = h1.IPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = h1.IPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Macro uncertainty: h1", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# for h12:
irf_h12_shock_ipm <- ggplot(data = var8_bloom_JFV_6_h12.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = h12.IPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = h12.IPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=h12.IPM.Lo, ymin=h12.IPM.Up), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = h12.IPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = h12.IPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "Macro uncertainty: h12", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

irf_h1_shock_ipm
irf_h12_shock_ipm

# oirf for EMPM:
# h1:
irf_h1_shock_empm <- ggplot(data = var8_bloom_JFV_6_h1.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = h1.EMPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = h1.EMPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=h1.EMPM.Up, ymin=h1.EMPM.Lo), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = h1.EMPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = h1.EMPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# h12:
irf_h12_shock_empm <- ggplot(data = var8_bloom_JFV_6_h12.var.irfs.df.rescaled) +
  # geom_point() + 
  geom_line(aes(x = step, y = h12.EMPM.Up), color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(aes(x = step, y = h12.EMPM.Lo), color="#e80628", 
            size=0.8, linetype = 3) +
  geom_ribbon(aes(x = step, ymax=h12.EMPM.Up, ymin=h12.EMPM.Lo), 
              fill="#cecaca", alpha=.3) +
  geom_line(aes(x = step, y = h12.EMPM.oirf), color="black", 
            size=0.8) +
  geom_point(aes(x = step, y = h12.EMPM.oirf), color="black", 
             size=0.8) +
  # annotation_custom("decreasing %<->% increasing") +
  scale_x_continuous(name = NULL, 
                     limits = c(0, 36), 
                     breaks = seq(0, 36, by = 12),
                     minor_breaks = NULL) + 
  scale_y_continuous(name = "", 
                     limits = c(-4, 4), 
                     breaks = seq(-4, 4, by=2), 
                     minor_breaks = NULL) + 
  ggtitle("") + 
  # theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
             size=1) +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  # change ratio of y and x - axis
  coord_fixed(ratio = 4)

# print to device
grid.arrange(irf_h1_shock_ipm, irf_h1_shock_empm, 
             irf_h12_shock_ipm, irf_h12_shock_empm, ncol=2)


########################################################################################

# print to pdf
VAR8 <- arrangeGrob(irf_blo_shock_ipm, irf_blo_shock_empm, ncol=2)

ggsave("VAR8.pdf", VAR8)

