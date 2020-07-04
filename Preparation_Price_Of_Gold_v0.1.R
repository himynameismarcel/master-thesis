# 12.06.2020
# preparation of the series for the real price of gold:

library(readr)

gold.data <- read_csv("MacroTrends_Data_Download_clean.csv")
gold.data$yearmon <- as.yearmon(gold.data$date, "%Y%m%d")

returns <- as.matrix(diff(log(gold.data$real), lag=1))
gold.data$S2 <- rbind(NA, returns)

# we want the return-series to be in percentage:
gold.data <-  gold.data %>%
          dplyr::mutate(S2 = S2*100) %>%
          # and we drop all other superfluous columns
          dplyr::select(S2, yearmon)

# with this, the data should be in a shape to be 
# merged with the rest of the SVAR.data
# (so we switch to the other data.frame!)