#--------------------------------------------------------------------
#           Reading-in and transformation of
#           NBER RECESSION DATES:
#--------------------------------------------------------------------

# Note: at multiple occasions below we want to add in the NBER 
# recession dates which we retrieved from 
# https://fred.stlouisfed.org/series/USREC?utm
# _source=series_page&utm_medium=related_content&utm_term=other_
# formats&utm_campaign=other_format;
# therefore we read in the data and perform some manipulations to
# be able to seamlessly integrate the recession dates (year-month
# combinations market with 1) into our dataset for the uncertainty
# measures:
nber_recessions <- read.csv(file="NBER_Recessions_USA.csv", 
                            header=TRUE, sep=",")
# next, we split the variable 'DATE' into 'year', 'month' and 'day'
# using # the 'separate()' - function from the 'tidyr' - package:
nber_recessions <- separate(nber_recessions, "DATE", 
                            c("year", "month", "day"), sep = "-", 
                            remove=FALSE, convert=TRUE)
# we drop the 'day' - variable
nber_recessions <- nber_recessions %>%
  dplyr::select(-day)

# and create the variable 'my'
nber_recessions <- as.data.frame(nber_recessions %>%
                                   mutate(my = year + month/12))


# we quickly rename 'DATE' to 'Date'
nber_recessions <- nber_recessions %>%
  dplyr::rename(Date = DATE) %>%
  dplyr::mutate(Date = as.Date(Date))

# the procedure that follows is identical to the procedure
# that we use in another script to extract the exact start- 
# and end-dates of the bloom-shocks:

# let us first filter for all rows where the recession indicator = 1 
# to get an overview:
as.data.frame(nber_recessions %>% filter(USREC == 1))


# then we loop through all rows, and first check if the USREC equals 1
# or not:

# we assign nber_recessions to a new data.frame because
# we want to leave 'nber_recessions' as is to be possibly used 
# by other scripts!

nber_recessions2 <- nber_recessions
# we initialize the grouping variable x as follows:
x <- 2

# we add an empty column to our data frame nber_recessions:
nber_recessions2["recession_ID"] <- NA

# then we start the loop:
for (i in 1:nrow(nber_recessions2)) {
  if(nber_recessions2$USREC[i] == 1) {
    # if we detect a shock (i.e., 'USREC' == 1), then we have 
    # to proceed as follows:
    # we store the 'my' variable in the column 'start' and 'end'
    nber_recessions2[i, ncol(nber_recessions2)] <- x
    
  }else{
    # no shock
    x <- x + 1
  }
}


# we filter to see whether our ID-assignment indeed worked:
nber_recessions2 %>%
  filter(!is.na(recession_ID))  %>%
  group_by(recession_ID) %>%
  summarise(Count = n())

# next we group_by recession_ID and extract the respective 
# start-dates and store the retrieved data to a new 
# data-frame
recessions_start <- as.data.frame(nber_recessions2 %>%
                                    filter(!is.na(recession_ID))  %>%
                                    group_by(recession_ID) %>%
                                    filter(row_number()==1) %>%
                                    mutate(year_start=year, month_start=month, 
                                           my_start = my) %>%
                                    dplyr::select(recession_ID, year_start, 
                                                  month_start, my_start))

# next, we replicate the above query for the end-dates
# (note that in some scenarios start- and end-dates 
# might be identical!)
recessions_end <- as.data.frame(nber_recessions2 %>%
                                  filter(!is.na(recession_ID))  %>%
                                  group_by(recession_ID) %>%
                                  filter(row_number()==n()) %>%
                                  mutate(year_end=year, month_end=month, 
                                         my_end = my) %>%
                                  dplyr::select(recession_ID, year_end, 
                                                month_end, my_end))

# next, we can merge the two data-frames from above:
recessions_start_end <- merge(x = recessions_start, 
                              y = recessions_end, 
                              by = "recession_ID", 
                              all = TRUE, na.rm=T)
# and inspect the resulting data-frame:
recessions_start_end


# we re-arrange the sequence of columns
recessions_start_end <- recessions_start_end %>%
  dplyr::select(recession_ID, year_start, 
                month_start, year_end, 
                month_end, my_start, 
                my_end)


# next, we construct a year-month-variable both out of
# the pair 'year_start' & 'month_start' and 'year_end' & 'month_end'
recessions_start_end$yearmon_start <- as.yearmon(
  recessions_start_end$my_start, "%Y-%B")
recessions_start_end$yearmon_end <- as.yearmon(
  recessions_start_end$my_end, "%Y-%B")
# next, we add a helper-column that we need in the next stage as well:
recessions_start_end$helper_date <- format(
  recessions_start_end$yearmon_start, "%b")

# this means that we can now drop the 'year' and 'month' variables
recessions_start_end <- recessions_start_end %>%
  dplyr::select(-c(year_start, year_end, 
                   month_start, month_end))

# next we add a column that gives us the duration of the shock:
recessions_start_end <- recessions_start_end %>%
  mutate(duration = (
    recessions_start_end$yearmon_end - 
      recessions_start_end$yearmon_start)
    * 12 + 1) %>%
  dplyr::select(-helper_date)


# Lastly: actually, there is a slight problem, which we fix by adding
# one month's numeric value to my_end:
recessions_start_end$my_end <- recessions_start_end$my_end + (1/12)
# this makes sure that if we refer to an end-month that
# we assume that we are talking about the last day of that
# respective month!
#--------------------------------------------------------------------