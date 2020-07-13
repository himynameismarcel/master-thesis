## Marcel, 05.07.2020
########################################################################
### PART 2: construction of impulse-responses following Jordá (2005)
########################################################################

## PART 2 is separated into ....... parts:
## (2.1)  2.1 handles ........
## (1.2)  1.2 handles the .......


###############################
## (2.1) loading additional_vars.csv;
##       merge with variables from sp500_merge_vxo
##       (most of the data comes from FRED 
##       [https://fred.stlouisfed.org/series/])
##       the employment in manufacturing comes
##       from BLS (see https://data.bls.gov/pdq/SurveyOutputServlet)
###############################

# first we read in the file 'additional_vars.csv'
# additional_vars <- read.csv(file="additional_vars.csv", header=TRUE, sep=",",
#                            stringsAsFactors = FALSE)
additional_vars <- read.csv(file="additional_vars_ext.csv", header=TRUE, sep=",",
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

# and drop a few superfuous columns:
additional_vars <- as.data.frame(additional_vars %>%
                      dplyr::select(-c(DATE.1, Date, 
                                       Open, High, Low, Adj.Close, Volume)))


# the variable EMPM has to be retrieved from a separate dataset
# which we downloaded from the homepage of the BLS;
# we have to manipulate the data to bring it into our desired format:

# old version when still using the xlsx-package:
# EMPM <- read.xlsx("BLS_EMPM.xlsx", header = TRUE, sheetIndex = 1)

# after switching to the readxl-package:
EMPM <- read_excel("BLS_EMPM.xlsx")
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
