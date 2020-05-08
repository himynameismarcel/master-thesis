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
### Installing Packages
#######################
library(readxl)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(vars)
library(matlib)
library(forecast)

#######################
### Preliminary Steps
#######################
## The preliminary steps involve getting the data into our R-session
## In particular, we need   * U_M,
##                          * IPM, and
##                          * U_F
#set.seed(7)

## Because we have reached a stage where we want to improve our algorithm 
## and be sure that it is computationally correct, as a first stept, we 
## decided to use the replication-data from Ludvigson et al (2018) together 
## with our algorithm, to rule out, that the result of not getting a single 
## successful run in our simulation is due to the wrong data!

## Both, the macro uncertainty as well as the financial uncertainty measure 
## of Ludvigson et al (2018) come along with three forecast horizon: 
## h=1, h=3 and h=12;
## For our replication-exercise, we always consider h=1 (this is also the 
## series which Mr. Sia has provided us with in the excel-file
## "Replication_data_Ludvigson.xlsx"

SVAR.data <- read_excel("Replication_data_Ludvigson.xlsx", 
                        sheet = "Data")

# a quick check shows us that the columns are already correctly
# named, but the type for the 'Date'-column is set to num;
# the Date-column consists of a numeric holding 'yearmonth',
# hence we have to proceed as follows to change that:
SVAR.data$yearmon <- as.yearmon(as.character(SVAR.data$Date), "%Y%m")

# the data in SVAR.data runs from 07/1960 until 04/2015;
# without further ado, we subset SVAR.data into SVAR.data.sub
# which then only includes the two uncertainty measures (Um, Uf)
# and the time-series for industrial production:
SVAR.data.sub <- SVAR.data %>%
                        dplyr::select(Um, ip, Uf)

# the returns-data is also already available in the Excel-file
# "Replication_data_Ludvigson.xlsx", which is why, in our 
# case, it is also already available in the data.frame SVAR.data;
# for now, we want to store the time-series in a dedicated
# data.frame, but might change this procedure later!
# at the same time, we divide the values by 100 because the original
# series comes in %!
return.data <- SVAR.data %>% 
                     dplyr::select(S, yearmon) %>%
                     dplyr::mutate(S = S/100)


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
    for(k in 1:10000){
          # k <- 3
          print(k)
          print(x)
      
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
          # print(class(u_t))
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
          # all data in SVAR.data starts in 1960-07 and runs until 
          # xx; the data.frame has 690 rows;
          
          # the residuals in epsilon_t have 684 rows; we lose 6
          # data points due to the 6 lags in the VAR-model;
          # hence, the residuals in epsilon_t start in 1961-01-01!
          # (6 months later!)
          
          # Marcel (02.05.2020): The above seems to have been the case
          # when we created this script in June 2018;
          # When trying to re-run this script, we saw that because we 
          # were using the replication-data from Ludvigson et al (2018),
          # which only runs until April 2015,
          # contrary to the above, also when re-running this script,
          # the residuals in epsilon_t have 652 rows;
          # and because we lose 6 data points due to the 6 lags in the VAR-model,
          # the residuals in epsilon_t start only in 1961-01-01!
          # accordingly, we have updated the below code!
          
          # with this info, we first make a data.frame out of
          # epsilon_t
          epsilon_t <- as.data.frame(epsilon_t)
          # add a date-column
          dates = seq(from = as.Date("1961-01-01"), 
                      to = as.Date("2015-04-01"), by = 'month')
          epsilon_t$date <- dates # which has exactly length 652!
          # and we assign sensible col-names
          epsilon_t <- as_tibble(epsilon_t %>%
                dplyr::rename(financial_h1 = V1,
                              lip  = V2,
                              macro_h1 = V3))
          # we further split the 'Date' - variable into three pieces:
          epsilon_t <- separate(epsilon_t, "date", c("year", "month", "day"), 
                              sep = "-", 
                              remove=FALSE, convert=TRUE)
          # and we create the column 'yearmon'
          
          epsilon_t$yearmon <- as.yearmon(paste0(epsilon_t$year, epsilon_t$month), 
                                          "%Y %m")
          # and drop the other columns which we don't need anymore
          epsilon_t <- epsilon_t %>%
                      dplyr::select(-c(date, year, month, day))
          
          
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
          if(abs(epsilon_t %>% dplyr::filter(yearmon == "Oct 1987") %>%
             dplyr::select(financial_h1)) >= k1){
              # Marcel, 04.05.2020: Note: when comparing our code above with the 
              # code from Lugvigson et al (2018), I noticed that there seems to 
              # be a a sign-difference in the error-terms in that the shock-series
              # for the structural errors of F_t seem to be have the exact 
              # opposite sign when compared with the series that the code 
              # of Ludvigson et al (2018) generates; this seems to be limited 
              # to the shock-series of the financial variable only 
              # (since the constraint FE_3 which related to the real variable 
              # seems to be fine again.)
            
            
            
              # print("Passed first test! Continue!")
            
              # EVENT CONSTRAINT: FE_2:
              # postulated that there must be at least one month
              # during the financial crises 2007-2009 for which the financial 
              # uncertainty shock was large and positive!
            
              # note: the function 'any' checks:
              # given a set of logical vectors, is at least one of the 
              # values true?
              
              if(any(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                      yearmon <= "Jun 2009") %>%
                     dplyr::mutate(financial_h1 = financial_h1 * (-1)) %>%
                     # Marcel, 04.05.2020: Note: when comparing our code above with the 
                     # code from Lugvigson et al (2018), I noticed that there seems to 
                     # be a a sign-difference in the error-terms in that the shock-series
                     # for the structural errors of F_t seem to be have the exact 
                     # opposite sign when compared with the series that the code 
                     # of Ludvigson et al (2018) generates; this seems to be limited 
                     # to the shock-series of the financial variable only 
                     # (since the constraint FE_3 which related to the real variable 
                     # seems to be fine again.)
                     dplyr::select(financial_h1) >= k2)){
                
                
                 # print("Passed second test! Continue!")
                 
                 
                 # EVENT CONSTRAINT: FE_3:
                 # the third and final EVENT CONSTRAINT demands that
                 # any real activity shocks found during the
                 # Great Recession do not take any unusually large
                 # positive values!
                 # In particular, none of the lip-shocks is allowed to
                 # be larger k_3 = 2 to dismiss real activity shocks
                 # that are greater than two standard deviations above
                 # its sample mean during the Great Recession!
                 if(sum(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                            yearmon <= "Jun 2009") %>%
                       dplyr::select(lip) <= k3) == 19){
                   
                       # print("Passed third test! Continue!")
                       
    
                       # CORRELATION CONSTRAINTS: PRELIMINARIES
                       # Marcel, 03.05.2020, I hadn't yet implemented
                       # the correlation constraints; so here we go:
                       # Ludvigson et al (2018) declare S_t be a measure of the aggregate
                       # stock market return; u_St is the first order
                       # autoregressive residual for S_t; in their baseline specification
                       # the authors impose restrictions on the correlation of 
                       # u_St with uncertainty shocks;
                       
                       # for our analysis here we get the stock market return from 
                       # the data.frame 'return.data' which we have prepared above
                       # and which contains the variable 'S' 
                       # on which we have to run an AR(1)-estimation
                       # and store the corresponding residuals!
                       # at the same time we have to make sure that the time-series
                       # has the same length as the time-series which we have
                       # used in the analysis of the VAR(6)-system above 
                       # (ore more precisely, the the residual outputs of the 
                       # AR(1)-process and the structural errors of our main
                       # var-model from above have the same length!)
                       
                       # we fit an AR(1)-model to return.data$S:
                       ar_1 <- arima(as.ts(return.data$S), order=c(1, 0, 0))
                       # we store the AR(1)-models residuals (reduced-form residuals):
                       u_St <- residuals(ar_1)
                       
                       # similar to the VAR-data itself, the returns-data started
                       # in 1960-07; with the VAR-data we lost 6 entries due to the
                       # 6 lags we had estimated, with the returns-data we only
                       # estimated an AR(1)-process, hence, the residuals start in 
                       # 1960-08; interestingly, my checks showed that the time-series
                       # u_St has a length of 658 (as many as the original data-series);
                       # interestingly, in the Ludvigson-Code, the residuals from the
                       # AR(1)-estimation only have 657 entries, while the residuals
                       # from our AR(1)-estimation have 658 entires; 
                       # having checked the residuals in the Ludvigson-paper (object U),
                       # the entry corresponding to 'Black Monday' in October 1987
                       # is at position 327, while in our residuals it is at position 
                       # 328; and because Ludvigson et al (2018) only remove the first
                       # five entries for the calculation of the correlation, we,
                       # accordingly have to remove the first 6 entries!
                       # the rationale for the above, comes from the observation,
                       # that in epsilon_t, October 1987 is on position 322, and
                       # we need it to be on position 322 also in u_St!
                       # (hence, we have to remove the first 6 entries!)
                       
                       # next, we calculate the three sample correlations c_M, c_Y and c_F
                       # between the AR(1)-model-residuals u_St and the shocks
                       # u_Mt, u_Yt and u_Ft of the reduced VAR, respectively;
                       # these shocks are the reduced-form errors which are stored 
                       # in u_t;
                       
                       # also, from the detailed explanations above it should be clear that
                       # the respective series (from which we now calculate the correlation)
                       # have the same length.
                       
                       # u_t contains the reduced-form errors of the 
                       # original VAR;
                       u_t <- as.data.frame(u_t)
                       
                       # u_St is currently too long, hence we remove the
                       # first six entries:
                       u_St <- as.ts(tail(as.zoo(u_St), -6))
                       
                       
                       # Marcel (04.04.2020):
                       # For some reason, Ludvigson et al (2018) seem to be using
                       # the strucural errors for the calculation of the correlation;
                       # I have to figure out why, but continue for the 
                       # moment like this!
                       # --> investigate later!
                       
                       # sample correlation between reduced-form residual
                       # u_St and the reduced-form error u_Mt:
                       c_M <- cor(u_St, epsilon_t$macro_h1)
                       
                       # sample correlation between reduced-form residual
                       # u_St aand the reduced-form error u_Yt:
                       c_Y <- cor(u_St, epsilon_t$lip)
                       
                       # sample correlation between reduced-form residual
                       # u_St and the reduced-form error u_Ft:
                       c_F <- cor(u_St, (epsilon_t$financial_h1)*(-1))
                       # Note: the multiplication with (-1) stems from the fact
                       # that we have noticed above, that the residuals for
                       # the financial shock seem to have the opposite sign
                       # as compared to the results of Ludvigson et al (2018):
                       
                       # interestingly, the results for the correlations
                       # between the Ludvigson et al (2018) - code and our
                       # own code also show opposite signs;
                       # although we have compared the inputs (i.e. the series
                       # u_St) and 
                       
                       # after setting u_t to be a data.frame, we have
                       # to convert it to a matrix again for the 
                       # above calculations:
                       u_t <- as.matrix(u_t)
                       
                       # having calculated the above correlations, we are now in 
                       # a possition to proceed with the correlation constraints
                       # one after the other:
                       
                       # CORRELATION CONSTRAINT: FC_1:
                       # FC_1 requires that epsilon_M and epsilon_F (two of the structural
                       # errors) are negatively correlated with S_t;
                       # specifically, each individual correlation has to exceed lambda1
                       # in absolute terms!
                       
                       if((lambda1 - c_M) > 0 & (lambda1 - c_F) > 0){
                         # Marcel (04.04.2020): Interestingly, Ludvigson et al (2018),
                         # have implemented this without an equal sign!
                         
                         # print("Passed fourth test! Continue!")
                         
                         
                         # CORRELATION CONSTRAINT: FC_2:
                         # FC_2 requires that financial uncertainty shocks be more highly
                         # correlated with u_St than macro uncertainty shocks, according to
                         # a magnitued dictated by the lower bound lambda2:
                         if((abs(c_F) - lambda2*abs(c_M)) > 0){
                           # Marcel (04.04.2020): Interestingly, Ludvigson et al (2018),
                           # have implemented this without an equal sign!
                           
                           
                           # print("Passed fifth test! Continue!")
                           
                           
                           # CORRELATION CONSTRAINT: FC_3:
                           # Constraint FC_3 (similar) to FC_1 from above requires
                           # that epslion_M and epsilon_F are negatively correlated 
                           # with S_t;
                           # but FC_1 requires that each individual correlation exceeds 
                           # lambda1,
                           # while FC_3 requires that they collectively exceed lambda3:
                           
                           # for FC_3 we have to calculate c_MF which is (according to the
                           # paper) defined as: c_MF^2 = c_M^2 + c_F^2
                           # hence:
                           c_MF <- sqrt(c_M * c_M + c_F * c_F)
                           if(c_MF - lambda3 > 0){
                             # Marcel (04.04.2020): Interestingly, Ludvigson et al (2018),
                             # have implemented this without an equal sign!
                             
                             # print("Passed sixth test! Continue!")
                             
                             
                             # once all tests are successfully passed, we add 
                             # +1 to our counter-variable x!
                             x <- x + 1
                             
                             
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
                             helper <- as.matrix(epsilon_t %>% filter(
                               yearmon >= "Dec 2007" &
                                 yearmon <= "Jun 2009") %>%
                                 dplyr::select(financial_h1) - k2)
                             k2.maxG <- as.vector(helper[which(helper > 0)])
                             
                             # lastly, we check for industrial production
                             # for k2.maxG we subset epsilon_t with the expression
                             # that we used for the if-statement above:
                             k3.maxG <- as.vector(k3 -epsilon_t %>% filter(
                               yearmon >= "Dec 2007" &
                                 yearmon <= "Jun 2009") %>%
                                 dplyr::select(lip))                   
                             
                             
                             # we bind together all values associated with maxG:
                             maxG.all <- as.vector(unlist(c(k1.maxG, 
                                                            k2.maxG, k3.maxG)))
                             # next, we calculate the associated value for 
                             # maxG:
                             # the below is a numeric which holds the calculated
                             # maxG-value
                             maxG.value <- as.numeric(sqrt(t(maxG.all) %*% 
                                                             maxG.all))
                             
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
                   
                               
                       }else{
                         # print("Did not pass sixth test! New Draw!")
                         k <- k + 1                         
                       }
               
                       
                     }else{
                       # print("Did not pass fifth test! New Draw!")
                       k <- k + 1     
                     }
 
                     
                   }else{
                     # print("Did not pass fourth test! New Draw!")
                     k <- k + 1                    
                   }
                       
                 } else {
                   # print("Did not pass third test! New Draw!")
                   k <- k + 1 
                 }
       
              }else{
                # print("Did not pass second test! New Draw!")
                k <- k + 1
              }
            
          }else{
            # print("Did not pass first test! New Draw!")
            k <- k + 1
          }
          
    }      

    
#######################
### Figure 4 (Replication of Ludvigson et al (2018));
### Calculation of Impulse Response-Functions
#######################

    # To be able to calculate the impulse-response functions by hand
    # and to be able to make use of the impact matrices stored in 
    # 'A_0_valid', we first have to construct the Bi-matrices of the
    # reduced-form VAR out of the coefficients that are stored in my.var 
    # which we have created above!
    
    
    ## STEP 1: extracting reduced-form matrices from
    ## VAR-estimation from above
    #----------------------------------------------
    # because our model consists of 6 lags, we need 6 reduced-form
    # matrices B filled with the coefficients from the model
    # estimation;
    # further, because it is a three-variable VAR (i.e., a VAR(6)-3),
    # we know that each of the matrices has dimension 3x3!
    
    # the most efficient approach would be to loop through 
    # 'my.var$varresult' which, in our case, is a list item with 
    # 3 components (for the three estimated equations), and further
    # holds the coefficients in the nested list-item '$coefficients';
    
    l <- 3 # number of equations
    p <- 6 # number of lags
    
    # for(i in 1:l){
    #   print(sapply(my.var$varresult, `[[`, 1))
    # }
    
    # the below command extracts, in our case, a 3x19-matrix (3x18) + 1 
    # columns vector of constants;
    B_block <- t(sapply(my.var$varresult, `[[`, 1))
    
    # to efficiently create the reduced-form VAR-coefficients (the
    # B'i's out of this block-matrix), we run the following:
    
    
    # to do this we use a function; 'matsplitter()'
    matsplitter<-function(M, r, c) {
      rg <- (row(M)-1)%/%r+1
      cg <- (col(M)-1)%/%c+1
      rci <- (rg-1)*max(cg) + cg
      N <- prod(dim(M))/r/c
      cv <- unlist(lapply(1:N, function(x) M[rci==x]))
      dim(cv)<-c(r,c,N)
      cv
    }       
    
    # another option would be (slightly modified output,
    # padded with NAs)
    # mat_split <- function(M, r, c){
    #   nr <- ceiling(nrow(M)/r)
    #   nc <- ceiling(ncol(M)/c)
    #   newM <- matrix(NA, nr*r, nc*c)
    #   newM[1:nrow(M), 1:ncol(M)] <- M
    #   
    #   div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
    #   matlist <- split(newM, div_k)
    #   N <- length(matlist)
    #   mats <- unlist(matlist)
    #   dim(mats)<-c(r, c, N)
    #   return(mats)
    # }
    
    # finally we store the B_i's to an array
    B_array <- matsplitter(B_block, 3, 3)
    # and convert the output to a list
    B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
    # and the vector of constants is stored in the matrix (a column vector)
    # 'constant'
    constant <- as.matrix(B_block[, ncol(B_block)])
    
    ## STEP 2: recursively calculate corresponding
    ## matrices of the reduced-form MA(infinity) - representation
    ## (note: because we only look 60 steps ahead in our
    ## impulse-response-functions, we also consequently
    ## only need to calculate 60 psi-matrices!)
    #----------------------------------------------
    
    # in the below function,
    # "list" is the list of B-matrices that we have available!
        # in our case (since we have estimated a VAR(6)), we have
        # extracted 6 B-matrices with the coefficients from the
        # VAR-estimation!
    # k is the position of the psi-matrix we are calculating
        # (i.e. k = 10 means that we are calculating the 10th Psi-matrix!)
    # n is the number of variables in the VAR-system (in our case n = 3)
    calculate_psi_matrices <- function(B_list, k, n){
        if(k == 0){
           return(diag(n))       
        }
        if(k == 1){
          return(
            calculate_psi_matrices(B_list, k-1, n)%*%B_list[[1]]
          )        
        }
        if(k == 2){
          return(
            calculate_psi_matrices(B_list, k-1, n)%*%B_list[[1]] + 
            calculate_psi_matrices(B_list, k-2, n)%*%B_list[[2]]
          )        
        }
        if(k == 3){
          return(
            calculate_psi_matrices(B_list, k-1, n)%*%B_list[[1]] + 
            calculate_psi_matrices(B_list, k-2, n)%*%B_list[[2]] + 
            calculate_psi_matrices(B_list, k-3, n)%*%B_list[[3]]
          )        
        }
        if (k == 4){
          return(
            calculate_psi_matrices(B_list, k-1, n)%*%B_list[[1]] + 
            calculate_psi_matrices(B_list, k-2, n)%*%B_list[[2]] + 
            calculate_psi_matrices(B_list, k-3, n)%*%B_list[[3]] + 
            calculate_psi_matrices(B_list, k-4, n)%*%B_list[[4]]
          )        
        }
        if(k == 5){
          return(
            calculate_psi_matrices(B_list, k-1, n)%*%B_list[[1]] + 
            calculate_psi_matrices(B_list, k-2, n)%*%B_list[[2]] + 
            calculate_psi_matrices(B_list, k-3, n)%*%B_list[[3]] + 
            calculate_psi_matrices(B_list, k-4, n)%*%B_list[[4]] + 
            calculate_psi_matrices(B_list, k-5, n)%*%B_list[[5]]
          )        
        }
        if(k >=6){
          return(
            calculate_psi_matrices(B_list, k-1, n)%*%B_list[[1]] + 
            calculate_psi_matrices(B_list, k-2, n)%*%B_list[[2]] + 
            calculate_psi_matrices(B_list, k-3, n)%*%B_list[[3]] + 
            calculate_psi_matrices(B_list, k-4, n)%*%B_list[[4]] + 
            calculate_psi_matrices(B_list, k-5, n)%*%B_list[[5]] + 
            calculate_psi_matrices(B_list, k-6, n)%*%B_list[[6]]
          )        
        }
      }
      
    
    # we let the recursion run and store the calculated
    # Psi in a list:
    # initialize list:

    # because the recursion does not manage more than 20 at once,
    # we split the calculation of all 60 psi-matrices
    # into three chunks
    # (note higher i the more involved the recursion!)
    # note that the below calculates the psi-matrices
    # only until lag 25!
    # for i = 0 we should get back the Identity-Matrix!
    psi_list <- list()
    
    for(i in 0:10){
      #i <- 26
      psi_list[[length(psi_list)+1]] <- calculate_psi_matrices(B_list, i, 3)
      # and we rename the entry to know to which
      # iteration it belongs
      names(psi_list)[length(psi_list)] <- 
        paste0(c("psi_"), i)    
    }  
    
    
    ## STEP 3: having the list 'psi_list' at our
    ## disposal, we can finally calculate the Thetas
    ## (which are the coefficients of the structural 
    ## MA-representation and which will then finally hold
    ## the coefficients for the impulse-response-functions
    ## that we want/need for plotting below!)
    #----------------------------------------------
    ## remember: A_0_valid holds all matrices that passed
    ## our restrictions!
    
    ## this last step will now loop through all matrices stored
    ## in 'A_0_valid' and for each of them, multiply the respective
    ## A_0_valid with the sequence of psi's in 'psi_list';
    ## ultimately, the result will be stored in a list again;
    ## be aware that we have to multiply by the inverse of 
    ## the stored A_0_valid!
    
    # we initialize a list that will hold all Thetas
    # for each A_0_valid - matrix!
    # This means that for each valid A_0 - matrix, 
    # we will get as many Theta-entries as many psi
    # coefficient matrices we have at our disposal!
    Thetas <- list()
    
    for(z in 0:(length(A_0_valid)-1)){
      
        Thetas[[z+1]] <- 
              lapply(psi_list, function(x){
              x %*% inv(A_0_valid[[z+1]])
              })
              names(Thetas)[length(Thetas)] <- paste0(c("A_0_"), z)

    }
    
    #-------------------------------------------------------------------
    # extract coefficients to create impulse-responses:
    # now all that is left is to efficiently extract the respective
    # coefficients;
    # what complicates this proceure is that coefficients that blong
    # together have to be extracted from the exact same position from
    # the respective matrices;
    
    # we want to create 1 df with 763 columns and
    # 25 rows (the length of the list holding the
    # 'Theta' - coefficients) for each of the
    # combinations; we start with Fin_Macro:
    
    # the length of the list Thetas tells us how
    # often we have to execute the below code;

    # the below function extracts all relevant
    # coefficients for ONE POSITION!
    # this means that we need to run the below
    # function for each position in the matrix;
    # in our case 9 times!
    extract.coefs.irf <- function(impulse, response, row, col){
    
        # we initialize a data-frame:
        coefs.df <- data.frame(matrix(ncol = length(Thetas), 
                                      nrow = length(psi_list)))
        
        # # and change its name
        # name_df <- coefs.df
    
        for(m in 1:length(Thetas)){
        # loop to fill columns (through all ~700 elements of the 
        # list 'Theta' that holds for each of the admissable
        # A_0_valid - matrices the respective Thetas in 
        # dedicated coefficient matrices!) 
          
          for(b in 1:11){
          # loop to fill rows (in our case 26); to fill the 
          # coefficients at EACH STEP (hopefully in the future 
          # we'll have more steps!)
            
            coefs.df[b, m] <-  sapply(lapply(Thetas, "[[", b), 
                                      function(x) x[row, col])[m]
            colnames(coefs.df)[m] <- paste0(impulse, "_", response, "_", m)
          
          }
        }
        
        return(as_tibble(coefs.df))
    }
    
    # below we execute the above function 9 times to retrieve
    # the coefficients for each of the series of impulse-
    # response functions
    Fin_Macro.coefs <- 100*extract.coefs.irf("Fin", "Macro", 1, 3)
    Fin_Ipm.coefs <- 100*extract.coefs.irf("Fin", "Ipm", 2, 3)
    Fin_Fin.coefs <- 100*extract.coefs.irf("Fin", "Fin", 3, 3)
    Macro_Macro.coefs <- 100*extract.coefs.irf("Macro", "Macro", 1, 1)
    Macro_Ipm.coefs <- 100*extract.coefs.irf("Macro", "Ipm", 2, 1)
    Macro_Fin.coefs <- 100*extract.coefs.irf("Macro", "Fin", 3, 1)
    Ipm_Macro.coefs <- 100*extract.coefs.irf("Ipm", "Macro", 1, 2)
    Ipm_Ipm.coefs <- 100*extract.coefs.irf("Ipm", "Ipm", 2, 2)
    Ipm_Fin.coefs <- 100*extract.coefs.irf("Ipm", "Fin", 3, 2)
    
    # we add a step-counter to each of the separate data-frames:
    Fin_Macro.coefs$step <- seq.int(nrow(Fin_Macro.coefs))
    Fin_Ipm.coefs$step <- seq.int(nrow(Fin_Ipm.coefs))
    Fin_Fin.coefs$step <- seq.int(nrow(Fin_Fin.coefs))
    Macro_Macro.coefs$step <- seq.int(nrow(Macro_Macro.coefs))
    Macro_Ipm.coefs$step <- seq.int(nrow(Macro_Ipm.coefs))
    Macro_Fin.coefs$step <- seq.int(nrow(Macro_Fin.coefs))
    Ipm_Macro.coefs$step <- seq.int(nrow(Ipm_Macro.coefs))
    Ipm_Ipm.coefs$step <- seq.int(nrow(Ipm_Ipm.coefs))
    Ipm_Fin.coefs$step <- seq.int(nrow(Ipm_Fin.coefs))
    
    # the below command transforms all of the above data-frames
    # into a tidy format for plotting and at the same time
    # binds them together
    plot_SVAR.irfs.all <- bind_rows(
          Macro_Macro.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Macro Uncertainty",
                          response = "Macro Uncertainty"),
          Macro_Ipm.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Macro Uncertainty",
                          response = "Ind. Production"),
          Macro_Fin.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Macro Uncertainty",
                          response = "Fin. Uncertainty"), 
          Ipm_Macro.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Ind. Production",
                          response = "Macro Uncertainty"), 
          Ipm_Ipm.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Ind. Production",
                          response = "Ind. Production"),
          Ipm_Fin.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Ind. Production",
                          response = "Fin. Uncertainty"), 
          Fin_Macro.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
              dplyr::mutate(impulse = "Fin. Uncertainty",
                            response = "Macro Uncertainty"),
          Fin_Fin.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Fin. Uncertainty",
                          response = "Fin. Uncertainty"),
          Fin_Ipm.coefs %>%
            gather(
              key = series_name,
              value = series_value,
              -c(step)) %>%
            # at the same time we add the names for impulse and response
            dplyr::mutate(impulse = "Fin. Uncertainty",
                          response = "Ind. Production")
    )
    
    
    
    # we need to create ordered factors, otherwise 'facet_wrap' combines the plots
    # in alphabetical order:
    plot_SVAR.irfs.all <- plot_SVAR.irfs.all %>%
      ungroup %>%
      mutate(series_name = factor(series_name, 
                                 ordered = TRUE,
                                 levels=unique(series_name)),
             response = factor(response,
                               ordered = TRUE,
                               levels=unique(response)),
             impulse = factor(impulse,
                               ordered = TRUE,
                               levels=unique(impulse)))

    
    
    #-------------maxG-solutions----------------------
    # the below takes care of the maxG-solutions:    
    # to automate this process, we extract the name of the maxG-solution
    # from maxG.df
    name.iteration.maxG <- maxG.df$iteration
    
    # and then extract from epsilon_t_valid using this value:
    A_0_maxG <- A_0_valid[[name.iteration.maxG]]

    # we write the maxG-solution to the collection of Thetas
    Thetas[[length(Thetas)+1]] <- 
      lapply(psi_list, function(x){
        x %*% inv(A_0_maxG)
      })
    names(Thetas)[length(Thetas)] <- paste0(c("A_0_maxG"))
    
    
    
    extract.coefs.irf_maxG <- function(impulse, response, row, col){
      
      # we initialize a data-frame:
      coefs.df <- data.frame(matrix(ncol = 1, 
                                    nrow = length(psi_list)))
      
      # # and change its name
      # name_df <- coefs.df
      
      # we are only interested in the maxG-solution
      m <- length(Thetas)
        # loop to fill columns (through all ~700 elements of the 
        # list 'Theta' that holds for each of the admissable
        # A_0_valid - matrices the respective Thetas in 
        # dedicated coefficient matrices!) 
        
        for(b in 1:11){
          # loop to fill rows (in our case 26); to fill the 
          # coefficients at EACH STEP (hopefully in the future 
          # we'll have more steps!)
          
          coefs.df[b, 1] <-  sapply(lapply(Thetas, "[[", b), 
                                    function(x) x[row, col])[m]
          colnames(coefs.df)[1] <- paste0(impulse, "_", response, "_", "maxG")
          
        }
      
      return(as.tibble(coefs.df))
    }
    
    
    # below we execute the above function 9 times to retrieve
    # the coefficients for each of the maxG - series of impulse-
    # response functions
    Fin_Macro.coefs_maxG <- 100*extract.coefs.irf_maxG("Fin", "Macro", 1, 3)
    Fin_Ipm.coefs_maxG <- 100*extract.coefs.irf_maxG("Fin", "Ipm", 2, 3)
    Fin_Fin.coefs_maxG <- 100*extract.coefs.irf_maxG("Fin", "Fin", 3, 3)
    Macro_Macro.coefs_maxG <- 100*extract.coefs.irf_maxG("Macro", "Macro", 1, 1)
    Macro_Ipm.coefs_maxG <- 100*extract.coefs.irf_maxG("Macro", "Ipm", 2, 1)
    Macro_Fin.coefs_maxG <- 100*extract.coefs.irf_maxG("Macro", "Fin", 3, 1)
    Ipm_Macro.coefs_maxG <- 100*extract.coefs.irf_maxG("Ipm", "Macro", 1, 2)
    Ipm_Ipm.coefs_maxG <- 100*extract.coefs.irf_maxG("Ipm", "Ipm", 2, 2)
    Ipm_Fin.coefs_maxG <- 100*extract.coefs.irf_maxG("Ipm", "Fin", 3, 2)
    
    # we add a step-counter to each of the separate data-frames:
    Fin_Macro.coefs_maxG$step <- seq.int(nrow(Fin_Macro.coefs_maxG))
    Fin_Ipm.coefs_maxG$step <- seq.int(nrow(Fin_Ipm.coefs_maxG))
    Fin_Fin.coefs_maxG$step <- seq.int(nrow(Fin_Fin.coefs_maxG))
    Macro_Macro.coefs_maxG$step <- seq.int(nrow(Macro_Macro.coefs_maxG))
    Macro_Ipm.coefs_maxG$step <- seq.int(nrow(Macro_Ipm.coefs_maxG))
    Macro_Fin.coefs_maxG$step <- seq.int(nrow(Macro_Fin.coefs_maxG))
    Ipm_Macro.coefs_maxG$step <- seq.int(nrow(Ipm_Macro.coefs_maxG))
    Ipm_Ipm.coefs_maxG$step <- seq.int(nrow(Ipm_Ipm.coefs_maxG))
    Ipm_Fin.coefs_maxG$step <- seq.int(nrow(Ipm_Fin.coefs_maxG))
    
    
    # we do not add the maxG-solutions to the big data-frame
    # because want to plot the maxG-solutions separately
    # with another data.frame as source;
    # rather, we make a separate data.frame
    plot_SVAR.irfs.maxG <- bind_rows(
      Macro_Macro.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Macro Uncertainty",
                      response = "Macro Uncertainty"),
      Macro_Ipm.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Macro Uncertainty",
                      response = "Ind. Production"),
      Macro_Fin.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Macro Uncertainty",
                      response = "Fin. Uncertainty"), 
      Ipm_Macro.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Ind. Production",
                      response = "Macro Uncertainty"), 
      Ipm_Ipm.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Ind. Production",
                      response = "Ind. Production"),
      Ipm_Fin.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Ind. Production",
                      response = "Fin. Uncertainty"), 
      Fin_Macro.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Fin. Uncertainty",
                      response = "Macro Uncertainty"),
      Fin_Fin.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Fin. Uncertainty",
                      response = "Fin. Uncertainty"),
      Fin_Ipm.coefs_maxG %>%
        gather(
          key = series_name,
          value = series_value,
          -c(step)) %>%
        # at the same time we add the names for impulse and response
        dplyr::mutate(impulse = "Fin. Uncertainty",
                      response = "Ind. Production")
    )
    
    # we need to create ordered factors, otherwise 'facet_wrap' combines the plots
    # in alphabetical order:
    plot_SVAR.irfs.maxG <- plot_SVAR.irfs.maxG %>%
      ungroup %>%
      mutate(series_name = factor(series_name, 
                                  ordered = TRUE,
                                  levels=unique(series_name)),
             response = factor(response,
                               ordered = TRUE,
                               levels=unique(response)),
             impulse = factor(impulse,
                              ordered = TRUE,
                              levels=unique(impulse)))
    

    # finally we plot the impulse responses:
    impulse.responses_all.SVAR <- 
      ggplot(data=plot_SVAR.irfs.all, aes(x=step, y=series_value)) + 
      geom_line(aes(colour=series_name), alpha=0.6, size=1.5) +
      geom_line(data=plot_SVAR.irfs.maxG, aes(x=step, y=series_value),
                colour="black", size=2, linetype="dashed") +
      labs(color=NULL) + 
      geom_hline(yintercept=0, color = "#514e4e", 
                 size=1) +
      scale_x_continuous(name = NULL) + 
      scale_y_continuous(name = NULL) +
      theme(axis.text=element_text(size=10),
            plot.title = element_text(size=10, face="bold", hjust = 0.5),
            axis.title=element_text(size=10),
            legend.position="none",
            #legend.text=element_text(size=14),
            #axis.text.x=element_blank(),
            plot.margin = unit(c(1,1,1,1), "mm"),
            #panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_rect(colour="black", 
                                            fill="white", 
                                            size=1.5, 
                                            linetype="solid")) + 
      facet_wrap(impulse ~ response,
                 scales="free", strip.position = "top")
      # coord_cartesian(ylim = c(-0.5, 0.5))

    
      impulse.responses_all.SVAR
       
      ggsave(file="impulse_responses_all_SVAR_test.pdf")
       
#######################
### Figure 3 (Replication of Ludvigson et al (2018))
#######################
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
    #       * standardize the series financial_h1, lip, macro_h1
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
                          my = as.numeric(yearmon)+1/12,
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
                            e_macro = macro_h1,
                            e_ipm = lip,
                            e_fin = financial_h1 
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
    
    
#######################
### Figure 2 (Replication of Ludvigson et al (2018))
#######################
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
      
      

