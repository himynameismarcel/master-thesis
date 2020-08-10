#########################################################################################
### Marcel Kropp, 18.07.2020
### This script performs the search for structural shocks following LMN (2019);
### In their 2018 contribution, LMN solely formulated a-priori restrictions, 
### in the update to their paper they perform their algorithm unconstrained and 
### SEARCH the derived structural shocks for the entire sample period for 
### episodes of heightened uncertainty;

###############################
## loading libraries
###############################
if(!exists("foo", mode="function")) source("20200705_load_libraries_v0.1.R")


###############################
## reading in NBER recession dates
###############################
if(!exists("foo", mode="function")) source("20200705_nber_recession_dates_v0.1.R")


###############################
## reading in all other necessary scripts
###############################
if(!exists("foo", mode="function")) source("Preparation_Price_Of_Gold_v0.1.R")
if(!exists("foo", mode="function")) source("20200628_functions_v0.3.R")


### In their paper they write:
### With this in mind, we ask what can be said about the uncertainty shocks themselves
### BEFORE imposing any identifying restrictions. To address this question, we construct
### the unconstrained set, which is based on the reduced-form covariance
### restrictions alone, and then study when big shocks in this set have occurred over
### the sample period, by searching across the unconstrained set for the month in which
### the uncertainty shocks (the structural shocks epsilon_Ft and epsilon_Mt)
### are largest;

### To construct the unconstrained solution set, the authors initialize A_0_inv to be
### the unique lower-triangular Cholesky factor of Sigma_u with non-negative
### diagonal elements, P, and then rotate it by K = 1.5 million random orthogonal
### matrices Q. Each rotation begins by drawing an n x n matrix M.
### Then Q is taken to be the orthonormal matrix in the QR decomposition of M. Since
### A_0_inv = PQ, the procedure imposes the covariance restrictions 
### by construction. Let epsilon = A_0 u_t be the shocks implied by a 
### for given u_t. The moments implied by the covariance structure alone give
### 1.5 million values of A_0_inv, and thus 1.5 unconstrained values of epsilon.





# To perform the above, we simply call the function 'LMN_algorithm' from the 
# script '20200628_functions_v0.3.R':
search.shocks <- LMN_algorithm("NO_constr", 1500000)

# now to extract the epsilon_t_valid, we have to access the third 
# list-item of the function's return-obejct:
epsilon.search <- search.shocks[[3]]


# extracting the macro uncertainty shocks:

        # and then we create data.frames from all data.frames in the list
        # where we column-wise bind together the respective columns that
        # belong together (i.e. macro_h1, financial_h1, etc.)
        epsilon.search.str.macro <- data.frame(do.call(cbind, 
                          lapply(epsilon.search, function(x) x$macro_h1)))
        
        # we extract the dates-column from one of the data-frames in the list:
        dates.search <- epsilon.search$iteration_1$yearmon
        
        # and then give the row-names the dates:
        rownames(epsilon.search.str.macro) <- dates.search
        
        # next we extract the row-name corresponding to every column's maximum 
        # value and create a separate data.frame for this:
        dates.max <- tibble(date.max = as.yearmon(
                            rownames(epsilon.search.str.macro)
                            [apply(epsilon.search.str.macro, 2, which.max)]))
        
        # now we want to create a relative frequency table of the extracted
        # dates in R:
        dates.freq = table(dates.max)   # apply the table function
        dates.relfreq.macro = sort(dates.freq / nrow(dates.max),
                             decreasing = TRUE)
        
# extracting the financial uncertainty shocks:

        # and then we create data.frames from all data.frames in the list
        # where we column-wise bind together the respective columns that
        # belong together (i.e. macro_h1, financial_h1, etc.)
        epsilon.search.str.fin <- data.frame(do.call(cbind, 
                            lapply(epsilon.search, function(x) x$financial_h1)))
        
        # we extract the dates-column from one of the data-frames in the list:
        # we already did this above!
        
        # and then give the row-names the dates:
        rownames(epsilon.search.str.fin) <- dates.search
        
        # next we extract the row-name corresponding to every column's maximum 
        # value and create a separate data.frame for this:
        dates.max <- tibble(date.max = as.yearmon(
          rownames(epsilon.search.str.fin)
          [apply(epsilon.search.str.fin, 2, which.max)]))
        
        # now we want to create a relative frequency table of the extracted
        # dates in R:
        dates.freq = table(dates.max)   # apply the table function
        dates.relfreq.fin = sort(dates.freq / nrow(dates.max),
                             decreasing = TRUE)


# after performing the above for various data-sets, we have the 
# following results at our disposal (for both macro and financial 
# uncertainty):
        dates.relfreq.macro.202004        # stock market data = SP500 data; rest: own built!
        dates.relfreq.macro.201904        # stock market data = SP500 data; rest: own built!
        
        dates.relfreq.macro.CRSP.201504   # stock market data = CRSP; rest: own built!
        dates.relfreq.macro.SP500.201504  # stock market data = SP500; rest: own built!
        # it is completely clear why the above produce the same results, because 
        # for the unconstrained version of the algorithm we do not make any use 
        # of the stock-market variable!
        
        dates.relfreq.macro.LMN.201504    # original LMN-data
        
        # interestingly, using the original LMN-data produced a slightly different
        # distribution than 
        # (i) using the original CRSP-data with the updated uncertainty-series or 
        # (ii) using the S&P500 with the updated uncertainty-series;
        
        # in particular, it is about switching the positions of Dec 1970 and 
        # Nov 1987!
        
        # for the financial uncertainty series we cannot identify any discrepancies
        # for the first four dates with the most maxima!
        
        dates.relfreq.fin.202004  
        dates.relfreq.fin.201904  
        
        dates.relfreq.fin.CRSP.201504  
        dates.relfreq.fin.SP500.201504 
        # it is completely clear why the above produce the same results, because 
        # for the unconstrained version of the algorithm we do not make any use 
        # of the stock-market variable!
        
        dates.relfreq.fin.LMN.201504     
        

# with the above results at our disposal, for the dates with the most
# maxima in the respective shock-series we want to calculate a 
# distribution of the values of the respective shock series for that
# particular day so that we can use our own thresholds for the respective
# shock-restricted estimations (and don't have to rely on the one's of LMN);

# looking at         
dates.relfreq.macro.LMN.201504 # and
dates.relfreq.fin.LMN.201504
# we want to calculate the quantiles for the first two dates

# before being able to filter for the respective dates and then calculate the 
# summary statistics, we have to add the date-column into the data.frame which 
# we had created above since we 'only' have named the respective row-names
# according to the date-columns!
epsilon.search.str.fin$Date <- dates.search
epsilon.search.str.macro$Date <- dates.search

# -------------------------------------------------------------------------------
# financial_uncertainty:
      epsilon.search.str.fin %>%
        # before applying a filter, we have to transform our data.frame 
        # into a tidy data-format:+
        pivot_longer(-Date, values_to = "Value", names_to = "iteration") %>%
        # firs we filter for the two dates of our interest:
        filter(Date == "Sep 2008"| Date == "Oct 1987") %>%
        # after filtering we should group by the two dates:
        group_by(Date) %>%
        # now we need to transform the row into a column
        summarise(q25 = quantile(Value, 0.25),
               q50 = quantile(Value, 0.50),
               q75 = quantile(Value, 0.75),
               q95 = quantile(Value, 0.95))


# -------------------------------------------------------------------------------
# macro_uncertainty:
      epsilon.search.str.macro %>%
        # before applying a filter, we have to transform our data.frame 
        # into a tidy data-format:+
        pivot_longer(-Date, values_to = "Value", names_to = "iteration") %>%
        # firs we filter for the two dates of our interest:
        filter(Date == "Sep 2008"| Date == "Nov 1987"  | Date == "Dec 1970") %>%
        # after filtering we should group by the two dates:
        group_by(Date) %>%
        # now we need to transform the row into a column
        summarise(q25 = quantile(Value, 0.25),
                  q50 = quantile(Value, 0.50),
                  q75 = quantile(Value, 0.75),
                  q95 = quantile(Value, 0.95))

### -------------------------------------------------------------------------------
### Having decided to implement this 'search' - algorithm, we also use this 
### to perform this search-algorithm for the longeest data-sample possible
### to see how the search-algorithm performs over the longer sample