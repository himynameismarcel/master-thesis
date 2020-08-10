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
if(!exists("foo", mode="function")) source("Preparation_Price_Of_Gold_v0.1.R")
if(!exists("foo", mode="function")) source("20200628_functions_v0.3.R")



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


return.data <- SVAR.data %>% 
  dplyr::select(S, yearmon) %>%
  dplyr::mutate(S = S/100) %>%
  rename(S1 = S)

# we merge in the data for the price of gold
other.data <- as_tibble(
  merge(x = return.data, 
        y = gold.data, 
        by = "yearmon", all = FALSE, na.rm=T)
)

# and then merge in the data for the log(ip)
other.data <- as_tibble(
  merge(x = other.data, 
        y = SVAR.data[c(4, 6)], 
        by = "yearmon", all = FALSE, na.rm=T)
)

# then we extract year and month and store them in two separate
# and distinct columns:
other.data <- other.data %>%
            mutate(year = year(yearmon),
                   month = month(yearmon))

# we create a variable 'my':
other.data <- other.data %>%
  mutate(my = year + month/12)




# SVAR.data contains all time-series we want to
# have a look at;
other.data <- other.data %>%
              dplyr::select(my, ip, S1, S2)

# to bring the dataset into a useful shape for our purposes,
# we 'gather()' the variables with the variable-Names
# together;
other.data <- other.data%>%
                gather( 
                  series, value, 
                  -my, na.rm=TRUE)


# new facets labels:
supp.labs <- c(`ip`= "log(ind. production)", 
               `S1`= "Aggr. stock market ret.",
               `S2`= "return of gold")


# now we should be able to plot the result:
# note: the below plot marks both entire periods with 'significant'
# volatility;
# 
other_data_plot_combined <- ggplot(other.data, aes(x=my,
                                                  y=value,
                                                  color = series)) +
  geom_rect(data=recessions_start_end[16:23, ], inherit.aes = FALSE,
            aes(xmin=my_start, xmax=my_end, ymin=-Inf, ymax=Inf), 
            fill='#606060', alpha=0.3) + 
  geom_line(size = 1) + 
  facet_grid(series ~ .,  labeller = as_labeller(supp.labs), scales="free_y") +
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
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "", 
                     breaks = seq(1960, 2018, by = 5),
                     minor_breaks = NULL) + 
  scale_color_manual(values = c("ip" = "#576aff", 
                                "S1" = "#cc52ba",
                                "S2" = "#fc9a49")) + 
  theme(legend.position="bottom",
        strip.text.y = element_text(size = 12, colour = "black")) + 
  labs(color=NULL) +
  theme(legend.position="none",
        axis.text=element_text(size=15),
        #plot.title = element_text(size=10, face="bold", hjust = 0.5),
        axis.title=element_text(size=10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


other_data_plot_combined

# ggsave("other_data_plot_combined.pdf")
















