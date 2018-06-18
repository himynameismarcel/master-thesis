#########################################################################################
### Marcel Kropp, 26.05.2018
### This script estimates Impulse Responses using Jordá's (2005) Local Projection
### Method; 
### range of data: 07/1962 - 06/2009 (following Bloom 2009)
### all variables detrended apart from uncertainty measures (following Bloom 2009)

### Alternative Local Projections without detrending will be estimated in another
### R-Script.


#########################################
### preparatory steps common to all VARs
#########################################
#####
## Read xlsx data and some basic manipulations
#####

#####
## log- and/or HP-transformations
#####

########################################################################
### PART 1: LocProj_Bloom2009_HP with Bloom-Shocks (marking months
###         with maximum volatility)
########################################################################

####################
## transformations
####################

####################
## VARs (VAR8)
####################




########################################################################
### PART 2: LocProj_Bloom2009_HP with actual VXO/volatility series
########################################################################



########################################################################
### PART 3: LocProj_Bloom2009_HP with Michigan Survey Index
###         following Leduc and Lui (2016)
########################################################################




####################
## plots of orifs
####################







###############################
## (8.3) loading additional_vars.csv;
##       merge with variables from sp500_merge_vxo
##       (most of the data comes from FRED 
##       [https://fred.stlouisfed.org/series/])
##       the employment in manufacturing comes
##       from BLS (see https://data.bls.gov/timeseries/CES3000000001)
###############################

# first we read in the file 'additional_vars.csv'
# additional_vars <- read.csv(file="additional_vars.csv", header=TRUE, sep=",",
#                            stringsAsFactors = FALSE)
# additional_vars <- read.csv(file="additional_vars_ext.csv", header=TRUE, sep=",",
#                            stringsAsFactors = FALSE)
additional_vars <- read.csv(file="additional_vars_ext_long.csv", header=TRUE, sep=",",
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

# and drop a few superfluous columns:
additional_vars <- as.data.frame(additional_vars %>%
                                   dplyr::select(-c( 
                                                    Date, Open, 
                                                    High, Low, 
                                                    Adj.Close, Volume)))


# the variable EMPM has to be retrieved from a separate dataset
# which we downloaded from the homepage of the BLS;
# (see https://data.bls.gov/timeseries/CES3000000001)
# we have to manipulate the data to bring it into our desired format:
# EMPM <- read.xlsx("BLS_EMPM.xlsx", header = TRUE, sheetIndex = 1)
EMPM <- read.xlsx("BLS_EMPM_long.xlsx", header = TRUE, sheetIndex = 1)

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
# while the only goal was to construct a consistent volatility - measure
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
            color="#e80628", 
            size=0.8,linetype = 3) +
  geom_line(data = regressions_out, aes(x = h, 
                                        y = 100*(exp(lo90_b_lip)-1)), 
            color="#e80628", size=0.8, 
            linetype = 3) +
  geom_ribbon(data = regressions_out, aes(x=h, ymax=100*(exp(up90_b_lip)-1), 
                                          ymin=100*(exp(lo90_b_lip)-1)), fill="#cecaca", alpha=.3) +
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
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", size=1) + 
  # Change line size
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

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
  geom_line(data = regressions_out, aes(x = h, y = 100*(exp(up90_b_lip)-1)), color="#e80628", size=0.8,linetype = 3) +
  geom_line(data = regressions_out, aes(x = h, y = 100*(exp(lo90_b_lip)-1)), color="#e80628", size=0.8, linetype = 3) +
  geom_ribbon(data = regressions_out, aes(x=h, ymax=100*(exp(up90_b_lip)-1), ymin=100*(exp(lo90_b_lip)-1)), fill="#cecaca", alpha=.3) +
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
  #theme_minimal() + 
  labs(color=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", size=1) + 
  theme(legend.position = c(0.93, 0.90), axis.text=element_text(size=14),
        axis.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

irf2

# and we save the plot to be used in our latex-document
ggsave("irf2.pdf")





