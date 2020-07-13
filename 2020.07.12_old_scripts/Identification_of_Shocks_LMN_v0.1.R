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


## -------------------------------------------------------------------
# then we use the approach out of 'Identification_of_Shocks.R' to 
# identify 'shocks' (i.e. episodes where the series goes beyond
# 1.65 std. devs)
## -------------------------------------------------------------------
# we create a threshold variable and store it in our dataset to 
# be able to create an indicator variable that allows us to extract 
# the periods that are beyond 1.65 std. dev. of the series

SVAR.data <- SVAR.data %>%
            mutate(thresh = 1.65) %>%
            # now we can create an indicator-variable
            mutate(shock_Uf = ifelse(Uf > 1.65, 1, 0),
                   shock_Um = ifelse(Um > 1.65, 1, 0))


## -------------------------------------------------------------------
# creation of a dedicated data.frame that also hold the respective
# start- and end-dates of the shock-periods
## -------------------------------------------------------------------
# To retrieve these respective 'start'- and 'end'-dates for the 
# shock-periods we apply the below procedure that results in a 
# dedicatd data.frame holding the variables 'start', 'end' for 
# all shocks!

# we create a variable 'my':
SVAR.data <- SVAR.data %>%
            mutate(my = year + month/12)
# to bring the dataset into a useful shape for our purposes,
# we 'gather()' the variables with the variable-Names
# together;
# note that before applying the 'gather' - function, we remove
# the series for 'EPU' (because we only want to keep EPU Historical!)
SVAR.data <- SVAR.data %>%
                gather( 
                  uncert_measure, value, 
                  -my, na.rm=TRUE)

# at this stage, we have to split the 


## -------------------------------------------------------------------
# ultimately, we also want to add the NBER recession dates to the
# plot (which we have also done in 'Identification_of_Shocks.R')
## -------------------------------------------------------------------




















