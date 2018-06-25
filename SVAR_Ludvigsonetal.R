#########################################################################################
### Marcel Kropp, 24.06.2018
### This script estimates SVARs following Ludvigson et al (2018)'s algorithm which
### we have described in detail in the main-text;

### In particular, the below runs the SVAR as suggested by Ludvigson et al (2018)
### for their baseline case which is a VAR consisting of
### (U_M, IPM, U_F) 
### where U_M stands for macro uncertainty,
### IPM for a measure of economic activity (here Industrial Production in M.)
### and U_F for a measure of financial uncertainty;

#######################
### Preliminary Steps
#######################
## The preliminary steps involve getting the data into the right format;
## In particular, we need   * U_M,
##                          * IPM, and
##                          * U_F (VXO in our case)
#set.seed(7)

## Reading in the measure for financial uncertainty as constructed
## by Ludvigson et al (2018):
## Note: Also the measure for financial uncertainty by Ludvigson et al. (2018)
## comes along with three forecast horizon: h=1, h=3 and h=12;
## For our replication-exercise, we solely onsider h=1!
financialUncertainty_index <- read_excel("FinancialUncertaintyToCirculate.xlsx", 
                                     sheet = "Financial Uncertainty")

## and rename three columns
financialUncertainty_index <- rename(financialUncertainty_index, 
                                 h1 = "h=1",
                                 h3 = "h=3",
                                 h12 = "h=12")

## make 'Date' - variable a 'Date' - type 
## (instead of POSIXct)
financialUncertainty_index$Date <- as.Date(financialUncertainty_index$Date)

# next, we create the three variable 'month', 'year' and 'day' using
# the 'separate()' - function from the 'tidyr' - package
financialUncertainty_index <- separate(financialUncertainty_index, "Date", 
                                   c("year", "month", "day"), sep = "-", 
                                   remove=FALSE, convert=TRUE)

# next we create the variable 'my' which is a numerical representation
# of yearmon:
financialUncertainty_index <- as.data.frame(financialUncertainty_index %>%
                                          mutate(my = year + month/12))


## The measure for macro uncertainty is stored in the data.frame
## called 'macroUncertainty_index';
## The measure for industrial production in manufacturing is stored
## in the data.frame 'additional_vars';
## To join all three series of U_M, IPM and I_F, we perform the below:
SVAR.data <- data.frame("macroUncert_h1" = macroUncertainty_index$h1,
                              "my" = macroUncertainty_index$my,
                              "year" = macroUncertainty_index$year,
                              "month" = macroUncertainty_index$month,
                              "Date" = macroUncertainty_index$Date)

SVAR.data <- as.tibble(right_join(x = financialUncertainty_index[, 
                                      c("h1","my")],
                                        y = SVAR.data,
                                        by = "my") %>%
                               right_join(x = additional_vars[, 
                                      c("ip","my")],
                                          y = .,
                                          by = "my") %>%
                               # further, we rename a few variables,
                               # select the ones we want to keep
                               dplyr::rename(financial_h1  = h1, 
                                             macro_h1 = macroUncert_h1) %>%
                               # we further reorder the variables
                               dplyr::select(Date, year, month, 
                                             macro_h1, ip, financial_h1)) %>%
                               # and transform ip to the log of ip (lip)
                               dplyr::mutate(lip = log(ip))

# filter data to belong to a certain 
# time-range (in our baseline case this
# is the 'Bloom-window' 
# (i.e., July 1962 - June 2008))
# SVAR.data.sub <- SVAR.data  %>%
#                       filter(Date <="2008-06-01" & Date >= "1962-07-01")
# alternative range: (i.e., July 1962 - June 2018))
# SVAR.data.sub <- SVAR.data  %>%
#   filter(Date <="2018-06-01" & Date >= "1962-07-01")  
SVAR.data.sub <- SVAR.data %>%
              dplyr::select(financial_h1, ip, macro_h1)



# # Preparation and Calculation of returns-data from SP500:
# sp500 <- as.data.frame(additional_vars$sp500)
# colnames(sp500) <- c("sp500_prices")
# 
# sp500_returns <- sp500 %>%
#         dplyr::mutate(sp500_ret = diff(log(sp500_prices)))
                

#######################
### Algorithm
#######################

    ##########
    ## STEP 1:
    ## Estimation of the reduced-form model and initialization of A_0^{-1} as the
    ## unique lower-triangular Cholesky factor P of Sigma_u with non-negative
    ## diagonal elements!
    ##########
    # estimation of model
    my.var <- VAR(SVAR.data.sub, type = "const", p = 6)
    # storing the variance-covariance matrix in Sigma_u
    Sigma_u <- summary(my.var)$covres
    
    # and storing the residuals for each variable in the matrix
    u_t <- residuals(my.var)
    
    # cholesky-decomposition of Sigma_u to retrieve P (
    # the unique lower triangular matrix);
    # because the function 'chol' returns an upper diagonal
    # matrix by default, we have to transpose the result;
    P <- t(chol(Sigma_u))
    
    # we initialize A_0^{-1} as the lower triangular matrix P
    A_0_inv <- P

    # we further initialize the necessary parameters:
    lambda1 <- -0.05
    lambda2 <- 2
    lambda3 <- 0.18
    k1 <- 4
    k2 <- 4
    k3 <- 2
    
    # we initialize a data-frame for the maxG-solutions:
    maxG.df <- data.frame(matrix(ncol = 2, nrow = 1))
    x <- c("maxG", "iteration")
    colnames(maxG.df) <- x
    # and we write the values 0 into the data.frame
    maxG.df[1, 1] <- 0
                
    
    # to check how many draws actually pass the second test so
    # far, we create the object x
    x <- 0
    
    # we further initialize a list for the matrices epsilon_t
    # that pass all constraints
    # and the matrices A_0 that generate those epsilon_t's that
    # pass all the tests
    epsilon_t_valid <- list()
    A_0_valid <- list()
    
    ## the below has to run 1.5 million times:
    for(k in 1:150000){
          #k <- 1
      
          ##########
          ## STEP 2:
          ## rotation of P by Q (i.e., right-multiplication of P with Q;
          ## To get a Q,
          ## the matrix M is drawn from NID(0,1); 
          ## note: NID stands for normally and independently distributed
          ## Q is then taken to be the
          ## orthonormal matrix resulting from the QR-decomposition of M
          ##########
      
          # random draw of matrix M
          # in our case n = 3:
          n <- 3
          M <- matrix(rnorm(n*n,mean=0,sd=1), n, n)
      
          # QR-decomposition of M:
          # note: the QR factorization is based on the result that any full rank
          # matrix can be decomposed into a orthogonal and upper triangular matrix;
          QR <- qr(M)
          Q <- qr.Q(QR)    
      
          ##########
          ## STEP 3:
          ## calculation of A_0^{-1} = P*Q:
          ##########       
          # matrix multiplication to retrieve
          # A_0_inv = P*Q
          A_0_inv <- A_0_inv %*% Q
          A_0_inv
          # for our below computations we also need A_0
          A_0 <- inv(A_0_inv)
    
          ##########
          ## STEP 4:
          ## before actually being able to expose the structural
          ## errors (residuals) \epsilon to the constraints, 
          ## we need to create a matrix with the structural 
          ## residuals implied by the respective A_0_inv which
          ## we have randomly generated above;
          ## for this, we take the matrix u_t which contains the
          ## models' reduced-form residuals, and for each row
          ## (which corresponds to a certain time 't'), multiple
          ## the respective row with A_0_inv
          ##########   
          
          # we take the no. of rows and columns from u_t
          n.rows <- dim(u_t)[1]
          n.cols <- dim(u_t)[2]
          
          # we initialize the matrix epsilon_t to have the
          # same size like u_t (for this we pre-allocate the
          # matrix and then fill it)
          epsilon_t <- matrix(nrow = n.rows, ncol = n.cols)

          # and then we loop through the rows of u_t,
          # and multiply each row with A_0,
          # and write it to epsilon_t:
          for(i in 1:n.rows){
            epsilon_t[i, ] <- A_0 %*% u_t[i, ]
          }
          
          #---------EVENT CONSTRAINTS
          # having epsilon_t at our disposal, we are in a position
          # to impose the EVENT CONSTRAINTS on the 
          # structural residuals epsilon_t:
          # EVENT CONSTRAINTS:
          # As outlined in the text, the event constraints require that
          # identified financial uncertainty shocks have plausible properties
          # during two episodes of heightened financial uncertainty:
          #     * the 1987 stock market crash
          #     * the 2007-2009 financial crisis
          
          
          #---------PRELIMINARIES: EVENT CONSTRAINTS
          # all data in SVAR.data starts in 1960-07-01; the data.frame
          # has 690 rows;
          # the residuals in epsilon_t have 684 rows; we lose 6
          # data points due to the 6 lags in the VAR-model;
          # hence, the residuals in epsilon_t start in 1961-01-01!
          # (6 months later!)
          
          # with this info, we first make a data.frame out of
          # epsilon_t
          epsilon_t <- as.data.frame(epsilon_t)
          # add a date-column
          dates = seq(from = as.Date("1960-12-01"), 
                       to = as.Date("2017-11-01"), by = 'month')
          epsilon_t$date <- dates
          # and we assign sensible col-names
          epsilon_t <- as.tibble(epsilon_t %>%
                dplyr::rename(financial_h1 = V1,
                              lip  =V2,
                              macro_h1 = V3))
          # we further split the 'Date' - variable into three pieces:
          epsilon_t <- separate(epsilon_t, "date", c("year", "month", "day"), 
                              sep = "-", 
                              remove=FALSE, convert=TRUE)
          # and we create a column 'yearmon'
          
          epsilon_t$yearmon <- as.yearmon(paste0(epsilon_t$year, epsilon_t$month), 
                                          "%Y %m")
          # and drop the other columns which we don't need anymore
          epsilon_t <- epsilon_t %>%
            dplyr::select(-c(date, year, month))
          
          
          #---------CHECKS: EVENT CONSTRAINTS        
          # Note: Because the constraints must all be fulfilled (which
          # translates into an 'AND' - condition), they can be tested
          # sequentially!
        
          # test <- as.ts(epsilon_t[, 1])
          # p <- autoplot(test)
          # ggplotly(p)
          
          
          # EVENT CONSTRAINT: FE_1:
          # the financial uncertainty shock that occurred in 
          # October 1987 must be a large one: in particular, it must be
          # beyond four standard deviations:
          if(epsilon_t %>% filter(yearmon == "Oct 1987") %>%
             dplyr::select(financial_h1) >= k1){
            
              print("Passed first test! Continue!")
            
              # EVENT CONSTRAINT: FE_2:
              # postulated that there must be at least one month
              # during the financial crises 2007-2009 for which the financial 
              # uncertainty shock was large and positive!
            
              # note: the function 'any' checks:
              # given a set of logical vectors, is at least one of the 
              # values true?
              
              if(any(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                      yearmon <= "Jun 2009") %>%
                 select(financial_h1) >= k2) > 0){
                 print("Passed second test! Continue!")
                 
                 
                 # EVENT CONSTRAINT: FE_3:
                 # the third and final EVENT CONSTRAINT demands that
                 # any real activity shocks found during the
                 # Great Recession do not take any unusually large
                 # positive values!
                 # In particular, none of the lip-shocks is allowed to
                 # be larger k_3 = 2 to dismiss real activity shocks
                 # that are greater than two standard deviations above
                 # its sample mean during the Great Recession!
                 if(any(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                            yearmon <= "Jun 2009") %>%
                       select(lip) >= k3) > 0){
                  print("Did not pass third test! New Draw!")                 
                 } else {
                   print("Passed third test! Continue!")
                   x <- x + 1
                   
                   
                   # CORRELATION CONSTRAINT: FC_1:
                   # As our 'external' financial variable
                   # we use the SP500 data which we have
                   # prepared above!
                   
                   
                   
                         # Note: out of the pool of A_0's and epsilon_t's that
                         # 'survived' (i.e., passed all tests),
                         # we want to compute the maxG-solution;
                         # To make this procedure computationally efficient,
                         # we have decided for the following procedure:
                         
                         # For the first A_0 that ends up in the solution 
                         # set, we calculate maxG
                         
                         # we need the below to calculate and store 
                         # the max-G-solutions:
                         
                         k1.maxG <-  as.vector(epsilon_t %>% 
                                                 filter(yearmon == "Oct 1987") %>%
                                                  dplyr::select(financial_h1) - k1)
                         
                         # for k2.maxG we subset epsilon_t with the expression
                         # that we used for the if-statement above:
                         helper <- as.matrix(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                                          yearmon <= "Jun 2009") %>%
                                      select(financial_h1) - k2)
                         k2.maxG <- as.vector(helper[which(helper > 0)])
                         
                         # lastly, we check for industrial production
                         # for k2.maxG we subset epsilon_t with the expression
                         # that we used for the if-statement above:
                         k3.maxG <- as.vector(k3 -epsilon_t %>% filter(
                                                      yearmon >= "Dec 2007" &
                                                      yearmon <= "Jun 2009") %>%
                                               select(lip))                   
        
                         
                         # we bind together all values associated with maxG:
                         maxG.all <- as.vector(unlist(c(k1.maxG, k2.maxG, k3.maxG)))
                         # next, we calculate the associated value for 
                         # maxG:
                         # the below is a numeric which holds the calculated
                         # maxG-value
                         maxG.value <- as.numeric(sqrt(t(maxG.all) %*% maxG.all))
                         
                         # then we check whether the currently stored value
                         # is smaller or larger than the new maxG-value:
                         if(maxG.value > maxG.df[1, 1]){
                           # if yes, we replace the entry
                           maxG.df[1, 1] <- maxG.value
                           # and store the iteration number next to it
                           maxG.df[1, 2] <- paste0(c("iteration_"), k)
                         }
                         
      
                         # for draws that passed all tests, we want to store
                         # the respective residuals and the matrix A_0
                         epsilon_t_valid[[length(epsilon_t_valid)+1]] <- epsilon_t
                         # and we rename the entry to know to which
                         # iteration it belongs
                         names(epsilon_t_valid)[length(epsilon_t_valid)] <- 
                           paste0(c("iteration_"), k)
                        
                         A_0_valid[[length(A_0_valid)+1]] <- A_0
                         # and we rename the entry to know to which
                         # iteration it belongs
                         names(A_0_valid)[length(A_0_valid)] <- 
                           paste0(c("iteration_"), k)
      
                         #break
                   
                 }
       
                
              }else{
                print("Did not pass second test! New Draw!")
                k <- k + 1
              }
            
            
          }else{
            print("Did not pass first test! New Draw!")
            k <- k + 1
          }
          
    }      
          
    #-----------Figure 3
    # epsilon_t_valid holds all residuals of the 
    # structural models that passed all tests;
    # A_0_valid holds all impact matrices
    # of the structural models that passed all tests;
    # Figure 3 in Ludvigson et al (2018) displays 
    #     * the date and size of
    #       e_Macro and e_Fin shocks that are at least two standard deviations
    #       above the mean and
    #     * negative e_ipm shocks exceeding two standard deviations for
    #       all solutions in the identified set!
    
    # to get to the desired Figure, we have to loop through all 
    # data.frames stored in 'epsilon_t_valid' and
    #       * standardize the series financial_ha, lip, macro_h1
    #       * make tidy data-frames out of the data-frames
    #       * and finally create indicator-variables with respect
    #         to the thresholds mentioned above!
    # We can achieve all of the above by means of 'lapply';
    epsilon_t_valid.figure3 <- lapply(epsilon_t_valid, function(x) {
            dplyr::mutate(x,
                          financial_h1  = (financial_h1 - mean(financial_h1, 
                            na.rm=TRUE))/
                            sd(financial_h1, na.rm=TRUE), 
                          lip  = (lip - mean(lip, 
                            na.rm=TRUE))/
                            sd(lip, na.rm=TRUE),
                          macro_h1  = (macro_h1 - mean(macro_h1, 
                            na.rm=TRUE))/
                            sd(macro_h1, na.rm=TRUE),
                          # and we transform the variable 'yearmon'
                          my =
                            as.numeric(epsilon_t_maxG$yearmon)+1/12,
                          # To make the calculations easier, we flip the sign of all
                          # entries for lip
                          lip = lip*(-1))

    })
    # next we use 'lapply' together with 'gather' to generate
    # a tidy data-format
    epsilon_t_valid.figure3 <- lapply(epsilon_t_valid.figure3, function(x) {
                        dplyr::select(x,
                                      -c(day, yearmon)) %>%
                          # first we rename the variables
                          dplyr::rename(
                            e_fin = financial_h1,
                            e_ipm = lip,
                            e_macro = macro_h1
                          ) %>%
                          gather(
                            key = series_name,
                            value = series_value,
                            -my)
    })
    # lastly, we create indicator-variables with respect to the
    # thresholds mentioned above!
    epsilon_t_valid.figure3 <- lapply(epsilon_t_valid.figure3, function(x) {
      dplyr::mutate(x,
                indicator = ifelse(series_value > 2, 1, 0))
    })
    
    # with the indicator-variable at our disposal, in each data.frame
    # we filter for indicator == 1
    epsilon_t_valid.figure3 <- lapply(epsilon_t_valid.figure3, function(x) {
      dplyr::filter(x,
                    indicator == 1)
    })
    
    # and finally stack all data.frames together:
    epsilon_t_valid.figure3 <- do.call("rbind", epsilon_t_valid.figure3)
    # and we order by 'my'
    epsilon_t_valid.figure3 <- epsilon_t_valid.figure3 %>%
              dplyr::arrange(my)
    
    # now we are finally in a position to plot figure 3 as well:
    time_series_epsilon_t_largeShocks <- ggplot(data = epsilon_t_valid.figure3) +
      geom_rect(data=recessions_start_end, 
                inherit.aes = FALSE,
                aes(xmin=my_start, xmax=my_end, 
                    ymin=-Inf, ymax=+Inf), 
                fill='#606060', alpha=0.3) + 
      # geom_point(aes(x = my, y = series_value, color=series_name), 
      #           #color="black", 
      #           size=1.1) +
      geom_segment(data=epsilon_t_valid.figure3, aes( 
                                             x=my, 
                                             xend=my,
                                             y=-Inf, yend=series_value),
                   colour = "#0072B2", size = 1.2) +
      geom_hline(yintercept=3,
                 linetype="dashed",
                 size=0.9) +
      scale_x_continuous(name = "Year",  
                         limits = c(1960, 2015),
                         breaks = seq(1960, 2015, by = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(name = NULL, 
                limits = c(2, 7),
                breaks = seq(2, 8, by=2),
                minor_breaks = NULL) +
      labs(color=NULL) +
      theme(axis.text=element_text(size=10),
            plot.title = element_text(size=10, face="bold", hjust = 0.5),
            axis.title=element_text(size=10),
            legend.position="none",
            #legend.text=element_text(size=14),
            #axis.text.x=element_blank(),
            plot.margin = unit(c(1,1,1,1), "mm"),
            #panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(colour="black", fill="white", 
                                            size=1.5, linetype="solid")) + 
      # change ratio of y and x - axis
      coord_fixed(ratio = 5) + facet_grid(series_name ~ ., 
                                          scales="free_y")
    
      time_series_epsilon_t_largeShocks
      ggsave(file="time_series_epsilon_t_largeShocks.pdf")
    
    
    #-----------Figure 2
    # Replication of Figure 2 from Ludvigson et al (2018)
    # We want to get the data into a format so that we can easily
    # plot the maxG-solution with ggplot and facets for the three
    # corresponding series of the maxG-solution;
    # For this purpose, we have to extract the list-entry with the
    # maxG-values from 'epsilon_t_valid'
    maxG.df # we see that the name of the iteration is 'iteration_2106'
    
    # to automate this process, we extract the name of the maxG-solution
    # from maxG.df
    name.iteration.maxG <- maxG.df$iteration
    
    # and then extract from epsilon_t_valid using this value:
    epsilon_t_maxG <- epsilon_t_valid[[name.iteration.maxG]]
    
    # and we need to transform the variable 'yearmon'
    epsilon_t_maxG <- epsilon_t_maxG %>%
                dplyr::mutate(my =
                                as.numeric(epsilon_t_maxG$yearmon)+1/12) %>%
                # and we standardize all residuals in the data.frame
                dplyr::mutate(
                  financial_h1  = (financial_h1 - mean(financial_h1, 
                                  na.rm=TRUE))/
                                  sd(financial_h1, na.rm=TRUE), 
                  lip  = (lip - mean(lip, 
                                  na.rm=TRUE))/
                                  sd(lip, na.rm=TRUE),
                  macro_h1  = (macro_h1 - mean(macro_h1, 
                                  na.rm=TRUE))/
                                  sd(macro_h1, na.rm=TRUE))             

      
    
    # next we have to gather to get to the data-frame that we 
    # actually want for plotting:
    epsilon_t_maxG.tidy <- epsilon_t_maxG %>%
                        dplyr::select(-c(day, yearmon)) %>%
                        # first we rename the variables
                        dplyr::rename(
                              e_fin = financial_h1,
                              e_ipm = lip,
                              e_macro = macro_h1
                        ) %>%
                            gather(
                                  key = series_name,
                                  value = series_value,
                                  -my)
                        
    
    # we need to create ordered factors, otherwise 'facet_wrap' combines the plots
    # in alphabetical order:
    epsilon_t_maxG.tidy <- epsilon_t_maxG.tidy %>%
      ungroup %>%
      mutate(series_name = factor(series_name, 
                                 ordered = TRUE,
                                 levels=unique(series_name)))
    
    
    # now we are in a position to replicate Figure 2 in
    # Ludvigson et al (2018):
    # now we can finally plot the irfs:
    time.series_epsilon_t_maxG <- ggplot(data = epsilon_t_maxG.tidy) +
      geom_rect(data=recessions_start_end, 
                inherit.aes = FALSE,
                aes(xmin=my_start, xmax=my_end, 
                    ymin=-Inf, ymax=+Inf), 
                fill='#606060', alpha=0.3) + 
      geom_line(aes(x = my, y = series_value, color=series_name), 
                color="#0072B2", 
                size=0.8) +
      geom_hline(yintercept=c(-3, 3),
                 linetype="dashed",
                 size=0.4) + 
      # geom_hline(yintercept=sd(series_value), 
      #            linetype="dashed", 
      #            color = "red", size=0.4) + 
      # geom_point(aes(x = step, y = oirf), color="black", 
      #            size=0.8) +
      # annotation_custom("decreasing %<->% increasing") +
      scale_x_continuous(name = "Year",  
                         limits = c(1960, 2015),
                         breaks = seq(1960, 2015, by = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(name = NULL) + 
                         # limits = c(-6, 8),
                         # breaks = seq(-6, 8, by=2),
                         # minor_breaks = NULL) +
      #ggtitle("% impact on industrial production") + 
      # theme_minimal() + 
      labs(color=NULL) +
      # geom_hline(yintercept=0, linetype="dashed", color = "#514e4e", 
      #            size=0.5) +
      theme(axis.text=element_text(size=10),
            plot.title = element_text(size=10, face="bold", hjust = 0.5),
            axis.title=element_text(size=10),
            legend.position="none",
            #legend.text=element_text(size=14),
            #axis.text.x=element_blank(),
            plot.margin = unit(c(1,1,1,1), "mm"),
            #panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(colour="black", fill="white", 
                                            size=1.5, linetype="solid")) + 
      # change ratio of y and x - axis
      coord_fixed(ratio = 5) + facet_grid(series_name ~ ., 
                                          scales="free_y")
    
      time.series_epsilon_t_maxG
      ggsave(file="time_series_epsilon_t_maxG.pdf")
    
      ## The procedure above was initialized with a reduced-form
      ## VAR model stored in my.var;
      ## With this information at hand, a set of A_0's were isolated 
      ## we estimate SVARs for the set of A_0's:
      # my.svar <- SVAR(x = my.var, estmethod = "direct", Amat = A_0_valid[[1]],
      #                 hessian = TRUE, method="BFGS")
      
      

