#########################################################################################
### Marcel Kropp, 15.04.2018
### This script is separated into several big parts:
### PART 2: deals with the EPU-index constructed by Baker, Bloom and Davis (2016)
###         Baker, Bloom and Davis (2016) have dedicated an entire web-page
###         to their constructed uncertainty measure where they not
###         only provide several indices for the US, but also several
###         other countries, categories, etc.
###         See the home-page at http://www.policyuncertainty.com/us_monthly.html
###         Because we deal with the US only, we have correspondingly only
###         downloaded data for the US!


###############################
## loading libraries
###############################
if(!exists("foo", mode="function")) source("20200705_load_libraries_v0.1.R")

###############################
## reading in NBER recession dates
###############################
if(!exists("foo", mode="function")) source("20200705_nber_recession_dates_v0.1.R")

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
# ggsave("epu_index_plot.pdf")


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
# ggsave("epu_index_historical_plot.pdf")


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