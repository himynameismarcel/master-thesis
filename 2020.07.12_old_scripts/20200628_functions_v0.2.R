## 28.06.2020
## This file contains all the functions that I have created and use
## in the course of the analyses for the Master's Thesis:

## Hence, before running any other script, this script should be executed
## before-hand.


## -------------------------------------
## Functions for Impulse Response Functions
## -------------------------------------
constr.identified_set.incl.maxG <- function(p, scenario)
    #######################
    ### Installing Packages
    #######################ds
    library(readxl)
    library(plyr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(vars)
    library(matlib)
    library(forecast)
    library(Hmisc)
    
    ###-------------------------------
    ### Preliminary Explanations
    ###-------------------------------
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
    
    
    ###-------------------------------
    ### Reading in Data
    ###-------------------------------
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
      dplyr::mutate(S = S/100) %>%
      rename(S1 = S)
    
    # 12.06.2020
    # at this stage, we add in the return-data for the 
    # gold-series which we have prepared in the script 
    # Preparation_Price_Of_Gold.R' and which is stored in
    # the data.frame 'gold.data'
    
    # 28.06.2020 (Marcel)
    # here, we source the script 'Preparation_Price_Of_Gold.R'
    if(!exists("foo", mode="function")) source("Preparation_Price_Of_Gold_v0.1.R")
    
    external.data <- as_tibble(
      merge(x = return.data, 
            y = gold.data, 
            by = "yearmon", all = FALSE, na.rm=T)
    )
    
    ###-------------------------------
    ### Algorithm
    ###-------------------------------
    
    set.seed(1)
    ##----------------------------
    ## STEP 1:
    ## Estimation of the reduced-form model and initialization of A_0^{-1} as the
    ## unique lower-triangular Cholesky factor P of Sigma_u with non-negative
    ## diagonal elements!
    ##----------------------------
    # estimation of model
    my.var <- VAR(SVAR.data.sub, type = "const", p = 6)
    # storing the variance-covariance matrix in Sigma_u
    Sigma_u <- summary(my.var)$covres * 633/652
    # alternative plot of Sigma_u without scientific notation:
    # format(summary(my.var)$covres, scientific=FALSE)
    
    # Marcel (06.05.2020): Note that in my opinion Ludvigson et al (2018)
    # have at this stage mistakenly only divided the sum of squares
    # (i.e. u_t' * u_t) of the residuals by 652 (=658-6) degrees of freedom
    # intsead of 633 (=658-(19+6))! --> Mention to Scharler when I send him
    # my thesis!
    # Marcel (17.05.2020): After studying this issue again, I came to the
    # conclusion, that, as mentioned in https://homepage.univie.ac.at/robert.kunst/var.pdf
    # Ludvigson et al (2018) simply divided by the sample size and didn't
    # apply a degrees-of-freedom correction.
    # Hence, above we have applied a correction factor to translate the
    # OLS-estimator into an MLE-estimate.
    
    # and storing the residuals for each variable in the matrix
    u_t <- residuals(my.var)
    
    # cholesky-decomposition of Sigma_u to retrieve P (
    # the unique lower triangular matrix);
    # because the function 'chol' returns an upper diagonal
    # matrix by default, we have to transpose the result;
    P <- t(chol(Sigma_u))
    
    
    # -----------------------------------------------
    # CORRELATION CONSTRAINTS: PRELIMINARIES
    # -----------------------------------------------
    
    # note that the returns-data starts in Jul 1960
    # and hence has 658 entries;
    
    # the entry corresponding to 'Black Monday' in October 1987
    # is at position 328; we accordingly
    # have to remove the first six entries!
    # the rationale for the above, comes from the observation,
    # that in epsilon_t, October 1987 is on position 322, and
    # we need it to be on position 322 also in the return-series
    # S1 and S2!
    
    external.data <- tail(external.data,-6)
    
    # jetzt sollte es so sein, dass die return-series
    # S1 und S2 genau gleich lang sind und matching
    # positions mit den epsilons aus dem 
    # structural VAR!
    
    
    # we initialize A_0^{-1} as the lower triangular matrix P
    # (note: we have already created P above via
    # P <- t(chol(Sigma_u)))
    A_0_inv <- P
    
    # we further initialize the necessary parameters:
    # k1 <- 4.2 # 4.2 std. dev. above mean in 1987:10 for e_Ft
    # k2 <- 4.6 # 4.6 std. dev. above mean in 2008:09 for e_Ft
    # k3 <- 4.7 # 4.7 std. dev. above mean in 2008:09 for e_Mt
    # k4 <- 4   # 4 std. dev. above mean in 1970:12 for e_Mt
    
    k1 <- 4.0 # 2.8 std. dev. above mean in 1987:10 for e_Ft
    k2 <- 4.1 # 4.6 std. dev. above mean in 2008:09 for e_Ft
    k3 <- 4.2 # 4.7 std. dev. above mean in 2008:09 for e_Mt
    k4 <- 3.9   # 4 std. dev. above mean in 1970:12 for e_Mt    
    
    
    # we initialize a data-frame for the maxG-solutions:
    maxG.df <- data.frame(matrix(ncol = 2, nrow = 1))
    x <- c("maxG", "iteration")
    colnames(maxG.df) <- x
    # and we write the values 0 into the data.frame
    maxG.df[1, 1] <- 0
    
    # date-column for the structural errors
    dates = seq(from = as.Date("1961-01-01"), 
                to = as.Date("2015-04-01"), by = 'month')
    
    # to check how many draws actually pass all the tests,
    # we create the object x
    x <- 0
    
    # 13.06.2020
    # Test-weise eingebaut:
    # a <- 0
    # b <- 0
    # c <- 0
    # d <- 0
    # e <- 0
    # f <- 0
    # g <- 0
    # h <- 0
    
    # 13.06.2020
    # Weil ich auch noch testweise überprüfen möchte, in welcher Konstellation
    # die Bedingungen bestanden bzw. nicht bestanden werden,
    # initialisieren wir an dieser Stelle auch noch eine Liste:
    restrictions.pass_nopass <- list()
    
    
    n <- 3 # stands for the number of variables!
    
    # we further initialize a list for the matrices epsilon_t
    # that pass all constraints
    # and the matrices A_0 that generate those epsilon_t's that
    # pass all the tests
    epsilon_t_valid <- list()
    A_0_valid <- list()
    
    ## the below has to run 1.5 million times:
    for(k in 1:p){
      
      # 13.06.2020
      # Test-weise eingebaut:
      a <- 0
      b <- 0
      c <- 0
      d <- 0
      e <- 0
      f <- 0
      g <- 0
      h <- 0
      
      # k <- 3
      print(k)
      # print(x)
      
      ##----------------------------
      ## STEP 2:
      ## rotation of P by Q (i.e., right-multiplication of P with Q;
      ## To get a Q,
      ## the matrix M is drawn from NID(0,1); 
      ## note: NID stands for normally and independently distributed
      ## Q is then taken to be the
      ## orthonormal matrix resulting from the QR-decomposition of M
      ##----------------------------
      
      # random draw of matrix M
      # in our case n = 3:
      # set.seed(1)
      # n <- 3
      M <- matrix(rnorm(n*n,mean=0,sd=1), n, n)
      
      # QR-decomposition of M:
      # note: the QR factorization is based on the result that a
      # any full rank matrix can be decomposed into a orthogonal 
      # and upper triangular matrix;
      # note that we only need the Q-part (out of the QR-factorization)
      # to continue with the below part!
      QR <- qr(M)
      Q <- qr.Q(QR)    
      
      # Marcel, 06.05.2020: In their own code in the script
      # rand_orthomat.m, Ludvigson et al. (2018) perform an operation
      # which normalizes 
      
      # see: https://stackoverflow.com/questions/38426349/how-to-create-random-
      # orthonormal-matrix-in-python-numpy/38430739#38430739
      # and https://web.eecs.umich.edu/~rajnrao/Acta05rmt.pdf
      # we extract the signs of the matrix R:
      R.sign <- sign(qr.R(QR))
      # and lastly replace the elements to the top-right of the diagonal
      # with zeros:
      R.sign[upper.tri(R.sign)] <- 0
      
      # in a last step we multiply Q with R.sign
      Q <-Q %*% R.sign
      
      
      ##########
      ## STEP 3:
      ## calculation of a possible impact matrix A_0^{-1} = P*Q:
      ##########       
      # matrix multiplication to retrieve
      # A_0_inv = P*Q
      # Note: A_0_inv = P (and P is the unique
      # lower-triangular Cholesky factor)
      A_0_inv <- A_0_inv %*% Q
      # A_0_inv
      
      # Marcel (17.06.2020):
      # On p. 4 of their paper, Ludvigson et al. (2018) write
      # that the restrictions on the sign of A_inv_0 (i.e. 
      # that it is assumed to be positive) follows from combining
      # the unit effect normalization on H with the restriction 
      # sigma_jj >= 0; hence, we too implement this assumption
      # at this stage (we follow the same produced as above
      # when recovering the orthogonal matrices Q):
      A_0_inv.sign <- sign(A_0_inv)
      # and lastly replace the elements to the top-right of the diagonal
      # with zeros:
      A_0_inv.sign[upper.tri(A_0_inv.sign)] <- 0
      A_0_inv.sign[lower.tri(A_0_inv.sign)] <- 0
      
      # in a last step we multiply Q with R.sign
      A_0_inv <-A_0_inv %*% A_0_inv.sign
      
      # note that A_0_inv is a candidate for a 
      # possible impact matrix!
      
      # for our below computations we also need A_0
      A_0 <- inv(A_0_inv)
      
      # Marcel, 06.05.2020:
      # up until here our code and the code of Ludvigson et al (2018)
      # is exactly the same; 
      # in their code, Ludvigson et al. (2018) perform a step called 
      # 'Step 3: give possible impact matrix positive diagonal elements';
      # in my opinion, this step simply involves the multiplication of
      # A_0_inv with the identity-matrix --> hence, I have skipped this
      # step!
      
      ##########
      ## STEP 4:
      ## before actually being able to expose the structural
      ## errors (residuals)\epsilon to the constraints, 
      ## we need to create a matrix with the structural 
      ## residuals implied by the respective A_0_inv which
      ## we have randomly generated above;
      ## for this, we take the matrix u_t which contains the
      ## models' reduced-form residuals, and for each row
      ## (which corresponds to a certain time 't'), multiple
      ## the respective row with A_0_inv
      ##########   
      
      # Marcel (06.06.2020): since we have made the below generation 
      # of the possible structural errors more efficient, we do not 
      # need to manually create epsilon_t and/or set its dimensions
      # anymore! hence the below code is commented out!
      # we take the no. of rows and columns from u_t
      # n.rows <- dim(u_t)[1]
      # n.cols <- dim(u_t)[2]
      
      # we initialize the matrix epsilon_t to have the
      # same size like u_t (for this we pre-allocate the
      # matrix and then fill it)
      # epsilon_t <- matrix(nrow = n.rows, ncol = n.cols)
      
      # and then we loop through the rows of u_t,
      # and multiply each row with A_0,
      # and write it to epsilon_t:
      # print(class(u_t))
      # for(i in 1:n.rows){
      # epsilon_t[i, ] <- A_0 %*% u_t[i, ]
      # }
      
      # The above method, however, is very inefficient!
      # Therefore we want to try out a slightly smarter
      # procedure:
      epsilon_t <- t(apply(t(u_t), 2, function(row) A_0%*%row))
      
      # Marcel, 06.06.2020: Note that we have performed a test
      # to see whether the method for generating epsilon_t in the 
      # code of Ludvigson et al (2018) is the same as here
      # and indeed our test has shown that it is identical!
      # (as a test we have overwritte A_0 with Eposs from the 
      # code of Ludvigson et al (2018) to see if we get the 
      # exact same values!)
      
      
      # -------------------------------------------
      # EVENT CONSTRAINTS
      # -------------------------------------------
      # having epsilon_t at our disposal, we are in a position
      # to impose the EVENT CONSTRAINTS on the 
      # structural residuals epsilon_t:
      # EVENT CONSTRAINTS:
      # As outlined in the text, the event constraints require that
      # identified financial uncertainty shocks have plausible properties
      # during two episodes of heightened financial uncertainty:
      #     * the 1987 stock market crash
      #     * the 2007-2009 financial crisis
      
      
      # -------------------------------------------
      # PRELIMINARIES: EVENT CONSTRAINTS
      # -------------------------------------------
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
      
      # add a date-column which we have created outside of the loop
      epsilon_t$date <- dates # which has exactly length 652!
      # and we assign sensible col-names
      epsilon_t <- as_tibble(epsilon_t %>%
                               # Marcel (06.05.2020): Attention - at this position 
                               # I have accidently flipped the position of financial_h1 and 
                               # macro_h1; this also explains why my output was
                               # upside down!
                               dplyr::rename(macro_h1 = V1,
                                             lip  = V2,
                                             financial_h1 = V3))
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
      
      
      # -------------------------------------------
      # CORRELATION CONSTRAINTS: PRELIMINARIES
      # -------------------------------------------
      
      
      # Marcel (12.06.2020):
      # next, we calculate the four sample correlations c_MS1, c_FS1, 
      # c_MS2, c_FS2;
      
      # also, from the detailed explanations above it should be clear that
      # the respective series (from which we now calculate the correlation)
      # have the same length.
      
      c_M.S1 <- cor(external.data$S1, epsilon_t$macro_h1)
      c_F.S1 <- cor(external.data$S1, epsilon_t$financial_h1)
      c_M.S2 <- cor(external.data$S2, epsilon_t$macro_h1)
      c_F.S2 <- cor(external.data$S2, epsilon_t$financial_h1)
      
      
      # Marcel (06.06.2020): Note to self and also something I should
      # mention in the thesis: in the paper of Ludvigson et al (2018),
      # the authors write >= or <= for all the contraints; in the code 
      # itself they have implemented everything without an = 
      # in the contraint-statements!  
      

      # below, we have to test everything at once, depending on the 
      # particular scenario we are looking at:
      if(scenario == 'ALL'){
        if(
          # EVENT CONSTRAINTS:
          # big-shock events:
          # FE_1: financial_uncertainty
          (epsilon_t %>% dplyr::filter(yearmon == "Oct 1987") %>%
           dplyr::select(financial_h1) >= k1) &
          # ----------------------------------------------------#
          # FE_2: Part 1: financial_uncertainty
          (
            (epsilon_t %>% filter(yearmon == "Sep 2008") %>%
             dplyr::select(financial_h1) >= k2) |
            # FE_2: Part 2: macro_uncertainty
            (epsilon_t %>% filter(yearmon == "Sep 2008") %>%
             dplyr::select(macro_h1) >= k3)
          ) &
          # ----------------------------------------------------#
          # FE_3: macro_uncertainty
          (epsilon_t %>% filter(yearmon == "Dec 1970") %>%
           dplyr::select(macro_h1) >= k4) &
          # ----------------------------------------------------#
          # non-negative constraints:
          # FE_4
          ((sum(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                     yearmon <= "Jun 2009") %>%
                dplyr::select(lip))) <= 0) &
          # ----------------------------------------------------#
          # FE_5: Part 1
          (epsilon_t %>% filter(yearmon == "Oct 1979") %>%
           dplyr::select(financial_h1) >= 0) &             
          # FE_5: Part 2
          (epsilon_t %>% filter(yearmon == "Oct 1979") %>%
           dplyr::select(macro_h1) >= 0) &
          # ----------------------------------------------------#
          # FE_6: Part 1
          (sum(epsilon_t %>% filter(yearmon >= "Jul 2011" & 
                                    yearmon <= "Aug 2011") %>%
               dplyr::select(financial_h1) >= 0) == 2) &             
          # FE_6: Part 2  
          (sum(epsilon_t %>% filter(yearmon >= "Jul 2011" & 
                                    yearmon <= "Aug 2011") %>%
               dplyr::select(macro_h1) >= 0) == 2)  &
          # CORRELATION CONSTRAINTS:
          # ----------------------------------------------------#
          # FC_1
          (c_M.S1 <= 0 & c_F.S1 <= 0) &
          # FC_2
          (c_M.S2 >= 0 & c_F.S2 >= 0)        
      }


      ){
        x <- x + 1
        print(x)
        
        
        #------------------------
        ## Construction of maxG-solution
        #------------------------
        
        # Note: out of the pool of A_0's and epsilon_t's that
        # 'survived' (i.e., passed all tests),
        # we want to compute the maxG-solution;
        # To make this procedure computationally efficient,
        # we have decided for the following procedure:
        
        # For the first A_0 that ends up in the solution 
        # set, we calculate maxG
        
        # we need the below to calculate and store 
        # the max-G-solutions:
        
        # according to LMN(2019):
        # the maxG-solution is defined as the single solution in the 
        # identified set for which the inequalities pertaining
        # to the external variable constraints are collectively
        # maximized as measured by an equally-weighted
        # quadratic norm:
        k.maxG <- (c_M.S1^2)*0.25 + (c_F.S1^2)*0.25 + 
          (c_M.S2^2)*0.25 + (c_F.S2^2)*0.25
        
        # we bind together all values associated with maxG:
        maxG.all <- as.vector(unlist(c(k.maxG)))
        
        # next, we calculate the associated value for 
        # maxG:
        # the below is a numeric which holds the calculated
        # maxG-value
        maxG.value <- as.numeric(maxG.all)
        
        # then we check whether the currently stored value
        # is smaller or larger than the new maxG-value:
        if(maxG.value > maxG.df[1, 1]){
          # if yes, we replace the entry
          maxG.df[1, 1] <- maxG.value
          # and store the iteration number next to it
          maxG.df[1, 2] <- paste0(c("iteration_"), k)
        }
        
        
        #------------------------
        ## identified solution set and calculation of IRFs
        #------------------------
        # for draws of the orthonormal matrix which gives rise to a certain
        # A_0_inv (A_0) that passed all tests, we want to store
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
        
      }
      
    }



matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}   



## -------------------------------------
## Functions for Forward Error Variance Decomposition
## -------------------------------------

mspe.matrix <- function(Thetas, p, h){
  # Thetas[[1]] contains all the Thetas for the first admissable
  # solution for which we have calculated all thetas up until
  # a forecast horizon of h = 60!
  
  # note that the forecast horizon h == 1 corresponds to 
  # Theta_0 and that the expression Thetas[[1]][[1]]
  # extracts Theta_0 for the sequence of the Thetas
  # of the first admissable A_0!
  mspe.mat<- matrix(data = 0, nrow = 3, ncol = 3)
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- Thetas[[p]][[i]] %*% t(Thetas[[p]][[i]])
    mspe.mat <- mspe.mat + helper
    # print(test)
  }
  
  return(mspe.mat)
  
}


contr.forecast.err.mat <- function(Thetas, p, h){
  # Thetas[[1]] contains all the Thetas for the first admissable
  # solution of A_0 for which we have calculated all thetas up until
  # a forecast horizon of h = 200!
  
  # note that the forecast horizon h == 1 corresponds to 
  # Theta_0 and that the expression Thetas[[1]][[1]]
  # extracts Theta_0 for the sequence of the Thetas
  # of the first admissable A_0!
  
  # also, the output of this function is a matrix with
  # 9 numbers that represent the respective contributions of innovations
  # in variable k to the forecast error variance or MSE of the h-step 
  # ahead forecast of variable j:
  contr.forecast.err.matrix <- matrix(data = 0, nrow = 3, ncol = 3)
  
  e1 <- c(1, 0, 0) # Spaltenvektor
  e2 <- c(0, 1, 0) # Spaltenvektor
  e3 <- c(0, 0, 1) # Spaltenvektor
  
  contr.forecast.err.1.1 <- 0
  contr.forecast.err.1.2 <- 0
  contr.forecast.err.1.3 <- 0
  contr.forecast.err.2.1 <- 0
  contr.forecast.err.2.2 <- 0
  contr.forecast.err.2.3 <- 0
  contr.forecast.err.3.1 <- 0
  contr.forecast.err.3.2 <- 0
  contr.forecast.err.3.3 <- 0
  
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e1) %*% Thetas[[p]][[i]] %*% t(t(e1))))^2
    contr.forecast.err.1.1 <- contr.forecast.err.1.1 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e1) %*% Thetas[[p]][[i]] %*% t(t(e2))))^2
    contr.forecast.err.1.2 <- contr.forecast.err.1.2 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e1) %*% Thetas[[p]][[i]] %*% t(t(e3))))^2
    contr.forecast.err.1.3 <- contr.forecast.err.1.3 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e2) %*% Thetas[[p]][[i]] %*% t(t(e1))))^2
    contr.forecast.err.2.1 <- contr.forecast.err.2.1 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e2) %*% Thetas[[p]][[i]] %*% t(t(e2))))^2
    contr.forecast.err.2.2 <- contr.forecast.err.2.2 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e2) %*% Thetas[[p]][[i]] %*% t(t(e3))))^2
    contr.forecast.err.2.3 <- contr.forecast.err.2.3 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e3) %*% Thetas[[p]][[i]] %*% t(t(e1))))^2
    contr.forecast.err.3.1 <- contr.forecast.err.3.2 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e3) %*% Thetas[[p]][[i]] %*% t(t(e2))))^2
    contr.forecast.err.3.2 <- contr.forecast.err.3.2 + helper
    # print(test)
  }
  for(i in 1:h){
    # we loop through the Thetas
    # print(i)
    helper <- ((t(e3) %*% Thetas[[p]][[i]] %*% t(t(e3))))^2
    contr.forecast.err.3.3 <- contr.forecast.err.3.3 + helper
    # print(test)
  }
  
  # in a last step we fill the matrix with the respective
  # contributions:
  
  contr.forecast.err.matrix[1, 1] <- contr.forecast.err.1.1
  contr.forecast.err.matrix[1, 2] <- contr.forecast.err.1.2
  contr.forecast.err.matrix[1, 3] <- contr.forecast.err.1.3
  contr.forecast.err.matrix[2, 1] <- contr.forecast.err.2.1
  contr.forecast.err.matrix[2, 2] <- contr.forecast.err.2.2
  contr.forecast.err.matrix[2, 3] <- contr.forecast.err.2.3
  contr.forecast.err.matrix[3, 1] <- contr.forecast.err.3.1
  contr.forecast.err.matrix[3, 2] <- contr.forecast.err.3.2
  contr.forecast.err.matrix[3, 3] <- contr.forecast.err.3.3
  
  return(contr.forecast.err.matrix)
  
}


FEVD <- function(Thetas, h){
  
  # at the beginning of the function we create a data.frame which
  # will hold the final results of the function-call:
  
  df <- data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, 
                         c(paste0('prop1.1_horizon_',h), 
                           paste0('prop1.2_horizon_',h),
                           paste0('prop1.3_horizon_',h),
                           paste0('prop2.1_horizon_',h),
                           paste0('prop2.2_horizon_',h),
                           paste0('prop2.3_horizon_',h),
                           paste0('prop3.1_horizon_',h),
                           paste0('prop3.2_horizon_',h),
                           paste0('prop3.3_horizon_',h)))))
  
  
  # be aware that if we write h == 1 here, we technically mean
  # a forecast horizon of h == 0;
  # we have to choose this nomenclature, because the matrices
  # in Thetas start at position == 1 and on position == 1
  # is Theta_0 (which corresponds) to h == 0!
  
  # for the 1-step-ahead forecast we get:
  # be aware that in Thetas[[p]][[h]]
  # p stands for the particular solution in the identified set
  # (i.e. each solution contains 200 Thetas (because the forecast
  # horizon is 200 long!)
  # h stands for the forecast horizon
  
  # within the function we want to loop through all Thetas that
  # correspond to every solution in the identified set
  # (currently we have 80 A_0-matrices that pass the test,
  # hence there are also 80 sequences of Theta-matrices
  # for the respective forecast horizons):
  for(i in 1:length(Thetas)){
    
    ## forecast error variance decomposition of y_1t (aka Macro Uncertainty):
    # prop1 is the proportion of variance of the forecast error
    # in y_1t due to shock in macro_uncertainty:
    prop1.1 <- contr.forecast.err.mat(Thetas, i, h)[1, 1] / 
      mspe.matrix(Thetas, 1, h)[1, 1]
    
    # prop2 is the proportion of variance of the forecast error
    # in y_1t due to shock in industrial_production:
    prop1.2 <- contr.forecast.err.mat(Thetas, i, h)[1, 2]  /
      mspe.matrix(Thetas, i, h)[1, 1]
    
    # prop3 is the proportion of variance of the forecast error
    # in y_1t due to shock in financial_uncertainty:
    prop1.3 <- contr.forecast.err.mat(Thetas, i, h)[1, 3] /
      mspe.matrix(Thetas, i, h)[1, 1]
    
    ## forecast error variance decomposition of y_2t (aka industrial production):
    prop2.1 <- contr.forecast.err.mat(Thetas, i, h)[2, 1] / 
      mspe.matrix(Thetas, i, h)[2, 2]
    
    # prop2 is the proportion of variance of the forecast error
    # in y_1t due to shock in industrial_production:
    prop2.2 <- contr.forecast.err.mat(Thetas, i, h)[2, 2] /
      mspe.matrix(Thetas, i, h)[2, 2]
    
    # prop3 is the proportion of variance of the forecast error
    # in y_1t due to shock in financial_uncertainty:
    prop2.3 <- contr.forecast.err.mat(Thetas, i, h)[2, 3] /
      mspe.matrix(Thetas, i, h)[2, 2]
    
    ## forecast error variance decomposition of y_3t (aka financial uncertainty): 
    prop3.1 <- contr.forecast.err.mat(Thetas, i, h)[3, 1] / 
      mspe.matrix(Thetas, i, h)[3, 3]
    
    # prop2 is the proportion of variance of the forecast error
    # in y_1t due to shock in industrial_production:
    prop3.2 <- contr.forecast.err.mat(Thetas, i, h)[3, 2] /
      mspe.matrix(Thetas, i, h)[3, 3]
    
    # prop3 is the proportion of variance of the forecast error
    # in y_1t due to shock in financial_uncertainty:
    prop3.3 <- contr.forecast.err.mat(Thetas, i, h)[3, 3] /
      mspe.matrix(Thetas, i, h)[3, 3]
    
    
    # after each run of the loop, we want to append the results
    # to the data.frame which we have already created ablove:
    # note that we add the results as a new row, meaning that 
    # each row contains the FEVD at the specified forecast
    # horizon for every solution in the identified set!
    
    df[nrow(df) + 1,] = c(prop1.1, prop1.2, prop1.3, 
                          prop2.1, prop2.2, prop2.3,
                          prop3.1, prop3.2, prop3.3)
  }
  
  # at the end we return the filled data.frame
  return(df)
}



extract.table.FEVD <- function(df){
  # we want to know at which forecast horizon the respective FEVD
  # are maximized; therfore we add a column:
  df$h <- seq.int(nrow(df))
  
  # now we want to calculate the maximum (for which we use dplyr again):
  df <- df %>%
    mutate(U_m_Shock.range = abs(U_m_Shock_min - U_m_Shock_max),
           ip_Shock.range = abs(ip_Shock_min - ip_Shock_max),
           U_f_Shock.range = abs(U_f_Shock_min - U_f_Shock_max)
    ) %>%
    # and then we re-order the columns:
    dplyr::select(U_m_Shock_min, U_m_Shock_max, U_m_Shock.range,
                  ip_Shock_min, ip_Shock_max, ip_Shock.range,
                  U_f_Shock_min, U_f_Shock_max, U_f_Shock.range, h) 
  
  # next we take the rows where for U_m shock, ip-Shock and Uf_shock the
  # respective values are maximized:
  helper.max <- cbind(df %>% filter(U_m_Shock.range == 
                                      max(U_m_Shock.range)) %>% 
                        dplyr::select(U_m_Shock_min, U_m_Shock_max, h),
                      df %>% filter(ip_Shock.range == 
                                      max(ip_Shock.range)) %>% 
                        dplyr::select(ip_Shock_min, ip_Shock_max, h),
                      df %>% filter(U_f_Shock.range == max(U_f_Shock.range)) %>%
                        dplyr::select(U_f_Shock_min, U_f_Shock_max, h)  
  )
  # because the names for h are not unique, we use a different way than the dplyr-way:
  names(helper.max) <- c("U_m_Shock_min", "U_m_Shock_max", "h_U_m",
                         "ip_Shock_min", "ip_Shock_max", "h_ip",
                         "U_f_Shock_min", "U_f_Shock_max", "h_U_f")
  # because below we want to stack the data on top of each other, we have to
  # glue the forecast horizons into one particular column:
  helper.max <- transform(helper.max,
                          h=paste0(h_U_m,'_', h_ip,'_', h_U_f))
  # we transform the column 'h' to  character:
  helper.max$h <- as.character(helper.max$h)
  
  helper.max <- helper.max %>%
    dplyr::select(U_m_Shock_min, U_m_Shock_max, 
                  ip_Shock_min, ip_Shock_max,
                  U_f_Shock_min, U_f_Shock_max, 
                  h)
  
  # helper.min.max contains the respective maximum according to the definition of 
  # Ludvigson et al (2018):
  # in a last step, we collect the numbers we want for the respective forecast
  # horizons:
  df.selection <- df[c(1:5, 12, 201), c(1, 2, 4, 5, 7, 8, 10)]
  df.selection[nrow(df.selection)+1, ] <- helper.max    
  
  # in the very last step, we create the ranges as characters:
  df.selection.final <- data.frame("U_m_Shock"=character(),
                                   "ip_Shock"=character(),
                                   "U_f_Shock"=character(),
                                   "h"=character(),
                                   stringsAsFactors = FALSE)
  # the most efficient approach would be to loop through the rows
  # in df.selection:
  for(i in 1:nrow(df.selection)){
    helper <- c(paste0('[', 
                       format(round(df.selection[i, 1], digits=2), 
                              nsmall = 2), ', ', 
                       format(round(df.selection[i, 2], digits=2), 
                              nsmall = 2), ']'),
                paste0('[', format(round(df.selection[i, 3], digits=2), 
                                   nsmall = 2), ', ', 
                       format(round(df.selection[i, 4], digits=2), 
                              nsmall = 2), ']'),
                paste0('[', format(round(df.selection[i, 5], digits=2), 
                                   nsmall = 2), ', ', 
                       format(round(df.selection[i, 6], digits=2), 
                              nsmall = 2), ']'),
                # lastly we add the colum that gives us the respective horizon
                # we are looking at:
                df.selection[i, 7])
    # and then we add helper to the data.frame df.selection.final:
    df.selection.final[nrow(df.selection.final) + 1,] <- helper
  }
  
  return(df.selection.final)
}



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
    
    for(b in 1:61){
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
  
  for(b in 1:61){
    # loop to fill rows (in our case 26); to fill the 
    # coefficients at EACH STEP (hopefully in the future 
    # we'll have more steps!)
    
    coefs.df[b, 1] <-  sapply(lapply(Thetas, "[[", b), 
                              function(x) x[row, col])[m]
    colnames(coefs.df)[1] <- paste0(impulse, "_", response, "_", "maxG")
    
  }
  
  return(as.tibble(coefs.df))
}



## -------------------------------------
## Functions for identification of shocks in a time-series (here: 
## applied to the macro and financial uncertainty series of Ludvigson et al)
## -------------------------------------
extractShocks_LMN <- function(df, name=name){
  
  nam <- NULL
  mean <- NULL        # in our case: 0
  sd <- NULL          # in our case: 1
  
  
  # (1) detection of shocks (i.e, 0/1)
  # note that we slightly adjust the function here because
  # we do not want to use the HP-filter:
  # we standardize the column 'series'
  df <- df %>%
    mutate(series = (series - mean(series))/sd(series))
  
  # we create a threshold variable and store it in our dataset to 
  # be able to create an indicator variable that allows us to extract 
  # the 'Bloom' - shocks (this time across all variables) which we
  # want to have a look at!
  
  # to name this additional column also accordingly, we conserve our
  # name-object again:
  nam <- paste(colnames(df[4]), "thresh", sep="_")
  df[, nam] <- 1.65
  
  # now we can create an indicator-variable
  # (for this we conserve our 'nam' variable again):
  nam <- paste(colnames(df[4]), "shock", sep="_")
  df[, nam] <- ifelse(df[, 4] >
                        df[, 5], 1, 0)
  
  # having finally identified the shocks, we can proceed with 
  # the generation of a data.frame that holds the respective 
  # data 
  #     * for further use in a TABLE
  #     * and the data for further use in PLOTS
  # (see below!)
  
  
  
  # (2) creation of data.frame for PLOTS (if solo);
  # below under (3) we will create a data.frame
  # for the case that we have facetted data;
  # Having all the shocks now in our dataset allows us to 
  # create a dedicated data.farme that holds the respective 
  # start- end-dates of the shock-periods. but before we can 
  # do that, we need to CREATE 'start' and 'end'-dates for
  # the respective periods (which we will then ultimately 
  # use for plotting the episodes with shaded regions in our 
  # plots!)
  
  # To retrieve these respective 'start'- and 'end'-dates for the 
  # shock-periods we apply the below procedure that results in a 
  # dedicatd data.frame holding the variables 'start', 'end' for 
  # all shocks!
  
  # let us first filter for all rows where the shock = 1 to get 
  # an overview:
  # (and see whether it corresponds to the identified dates that 
  # Bloom also used himself!)
  df[df[, 6] == 1, c(1, 4, 6)]
  
  # continuation with procedure:
  # we want to give each episode of a shock a unique ID so that
  # we then can easily calculate the max volatility and first
  # volatility by grouping by ID;
  # to accomplish this, we decided for the following procedure:
  
  # we loop through all rows, and first check if the bloom_shock 
  # equals 1 or not:
  
  # we initialize the grouping variable x as follows:
  x <- 2
  
  # we add an empty column to our data frame:
  df["shock_ID"] <- NA
  
  # then we start the loop:
  for (i in 1:nrow(df)) {
    if(df[i, 6] == 1) {
      # if we detect a shock (i.e., 'bloom_shock' == 1), then 
      # we have to record a shock-ID:
      df[i, 7] <- x
      
    }else{
      # no shock
      x <- x + 1
    }
  }
  
  # we can now have a look at the newly created variable
  # (note that we have suppressed 'NA's for the creation 
  # of the table!):
  
  # the below command is just for ourselves to check!
  df %>%
    filter(!is.na(shock_ID))  %>%
    group_by(shock_ID) %>%
    summarise(Count = n())
  # we create a variable 'my':
  df <- df %>%
    mutate(my = year + month/12)
  
  
  # we continue with a subset of the variables and only look
  # at the rows where we actually have a shock_ID other than 
  # NA!
  # for this we extract the name of the column we want to,
  # among others, filter for:
  
  df.sub <- as.data.frame(df %>%
                            filter(!is.na(shock_ID))  %>%
                            dplyr::select(shock_ID, year, month, my,
                                          series))
  
  # we save this version of df.sub because it will be used
  # in one of the steps below to create another data-frame
  # in the style of Bloom (2009), Table A.1 which we want
  # to export as well
  df.sub.blo <- df.sub
  
  # next, we construct a year-month-variable both out of
  # the pair 'year' & 'month';
  # note that we had to add a small correction by '0.083' manually
  # so that the yearmon is created the way we want it!
  df.sub$yearmon <- 
    as.yearmon(df.sub$my-0.083, "%Y-%B")
  
  # next we group_by shock_ID and extract the respective MAXIMUM
  # value of the uncertainty series in the group:
  # (note that the produced data-frame only has one row per shock_ID
  # and that we go back to the data.frame 'df'
  # for this operation:
  
  df.sub.max.vol <- 
    as.data.frame(df %>%
                    filter(!is.na(shock_ID))  %>%
                    group_by(shock_ID) %>%
                    filter(series == max(series)) %>%
                    mutate(max_vol = series) %>%
                    dplyr::select(shock_ID, max_vol, my) %>%
                    # we add a column called max
                    mutate(max = 'm'))
  
  # note that we replace 'nam' for the name of the
  # column/variable we want to look at:
  df.sub.first.vol <- 
    as.data.frame(df %>%
                    filter(!is.na(shock_ID))  %>%
                    group_by(shock_ID) %>%
                    filter(row_number()==1) %>%
                    mutate(first_vol = series) %>%
                    dplyr::select(shock_ID, first_vol, my)%>%
                    # we add a column called first
                    mutate(first = 'f'))
  
  
  # now we can merge the three data-frames 
  # and see where we have max
  # and where we have first-values:
  df.sub <- 
    (right_join(
      x = df.sub.max.vol[, c("max", "my")], 
      y = df.sub, 
      by = "my") %>%
       right_join(
         x = df.sub.first.vol[, c("first", "my")], 
         y = ., 
         by = "my")) %>%
    # in a last step we also drop 'year', 'month'
    dplyr::select(-c(year, month)) %>%
    # and we rearrange the columns
    dplyr::select(yearmon, series, first, max,
                  shock_ID) %>%
    # this is a tweak so that we don't paste
    # NAs together in the 'unite' command below
    replace_na(list(first = "", max = "")) %>% 
    unite(label, first, max, sep=" ")
  
  # here we still should recalculate the shock_ID (starting from 1
  # again!)
  
  # in a last step, we want count the ID-column from scratch again!
  df.sub <- 
    df.sub %>%
    mutate(shock_ID = group_indices(., 
                                    factor(shock_ID, levels = unique(shock_ID))))
  
  # lastly, we should round numeric columns in both
  # data-frames!
  # we set the number of decimal place to 2 for numeric columns:
  is.num <- sapply(df.sub, is.numeric)
  df.sub[is.num] <- 
    lapply(df.sub[is.num], round, 2)
  
  # in the very last step, we have replace the 
  # name 'series' with the name the user gives to
  # the function as a string!
  colnames(df.sub)[2] <- name
  
  
  # (3) creation of data.frame for PLOTS (if facetted);
  # above we had already test-wise filtered for all rows where 
  # the shock = 1 to get an overview:
  # this time we perform this operation again and save the
  # result; note that the data.frame 'df' grew in the meantim
  # (with additional variables) that we do not need, there
  # we combine the filter with a select-statement:
  df.dates <- df %>%
    filter(series_shock == 1) %>%
    dplyr::select(my) %>%
    # we create new variables 'Lower' and 'Upper'
    # and fill them with -Ind and + Inf
    # and at the same time we create a variable
    # called 'uncert_measure' containing the string
    # of the respective
    # uncertainty measure:
    # note: the name should be the same like the names
    # we chose for the data.frame 'comparison_measures'
    # so that a join can be performed afterwards more easily!
    dplyr::mutate(Lower = -Inf, Upper = + Inf, uncert_measure = name)
  
  
  # (4) creation of data.frame for TABLES:
  # to retrieve the respective 'start'- and 'end' 
  # - dates for the shock-periods,
  # we apply a procedure to the dataset 
  # 'df' that results in a 
  # separate dedicated data-frame with the following 
  # variables:
  # 'start', 'end', 'max_volatility', 'first_volatility' 
  # for all shock-measures!
  
  # we group_by shock_ID and extract the respective start-dates and
  # store the retrieved data to a new data-frame
  df.start <- as.data.frame(df %>%
                              filter(!is.na(shock_ID))  %>%
                              group_by(shock_ID) %>%
                              filter(row_number()==1) %>%
                              mutate(year_start=year, month_start=month, 
                                     my_start = my) %>%
                              dplyr::select(shock_ID, year_start, 
                                            month_start, my_start))
  
  # next, we replicate the above query for the end-dates
  # (note that in some scenarios start- and end-dates are identical!)
  df.end <- as.data.frame(df %>%
                            filter(!is.na(shock_ID))  %>%
                            group_by(shock_ID) %>%
                            filter(row_number()==n()) %>%
                            mutate(year_end=year, month_end=month, 
                                   my_end = my) %>%
                            dplyr::select(shock_ID, year_end, 
                                          month_end, my_end))
  
  
  # # next, we can merge the two data-frames from above:
  df.start.end <- merge(x = df.start, 
                        y = df.end, by = "shock_ID",
                        all = TRUE, na.rm=T)
  
  # we re-arrange the sequence of columns
  df.start.end <- df.start.end %>%
    dplyr::select(shock_ID, my_start, my_end) %>%
    dplyr::rename(start = my_start,
                  end = my_end)
  
  # the above data-frame now allows us to plot the episodes of 
  # high volatility into our previous plots (by adding a shaded
  # rectangle!)
  
  # actually, there is a slight problem, which we fix by adding
  # one month's numeric value to my_end:
  # ask Martin about this!
  # df.start.end$end <- 
  #             df.start.end$end + (1/12)
  # this makes sure that if we refer to an end-month that
  # we assume that we are talking about the last day of that
  # respective month!
  
  # in a last step, we want count the ID-column from scratch again!
  df.start.end <- 
    df.start.end %>%
    mutate(shock_ID = group_indices(., 
                                    factor(shock_ID, 
                                           levels = unique(shock_ID))))
  
  # lastly, we should round numeric columns in both
  # data-frames!
  # we set the number of decimal place to 2 for numeric columns:
  is.num <- sapply(df.start.end, is.numeric)
  df.start.end[is.num] <- 
    lapply(df.start.end[is.num], round, 2)
  
  # at the very end we add a column that holds the name of the
  # series we are looking at:
  df.start.end <- df.start.end %>%
    dplyr::mutate(uncert_measure = name)
  
  
  
  
  
  # (5) creation of data.frame in the style of Bloom (2009)
  # Table A.1 (from Appendix)
  
  # we can continue working with the df.sub.blo 
  # which we had created above!
  
  df.sub.blo %>%
    filter(!is.na(shock_ID))  %>%
    group_by(shock_ID) %>%
    summarise(Count = n())
  
  # we group_by shock_ID and extract the respective start-dates 
  # and store the retrieved data to a new data-frame
  df.sub.blo.start <- as.data.frame(df.sub.blo %>%
                                      filter(!is.na(shock_ID))  %>%
                                      group_by(shock_ID) %>%
                                      filter(row_number()==1) %>%
                                      mutate(year_start=year, 
                                             month_start=month, 
                                             my_start = my) %>%
                                      dplyr::select(shock_ID, year_start, 
                                                    month_start, my_start))
  
  # next, we replicate the above query for the end-dates
  # (note that in some scenarios start- and 
  # end-dates are identical!)
  df.sub.blo.end <- as.data.frame(df.sub.blo %>%
                                    filter(!is.na(shock_ID))  %>%
                                    group_by(shock_ID) %>%
                                    filter(row_number()==n()) %>%
                                    mutate(year_end=year, 
                                           month_end=month, 
                                           my_end = my) %>%
                                    dplyr::select(shock_ID, 
                                                  year_end, month_end, 
                                                  my_end))
  
  # next, we can merge the two data-frames from above:
  df.blo.start.end <- merge(x = df.sub.blo.start, 
                            y = df.sub.blo.end, 
                            by = "shock_ID", all = TRUE, na.rm=T)
  # and inspect the resulting data-frame:
  df.blo.start.end
  
  # we re-arrange the sequence of columns
  df.blo.start.end <- df.blo.start.end %>%
    dplyr::select(shock_ID, year_start, 
                  month_start, year_end, 
                  month_end, my_start, my_end)
  
  
  # next, we construct a year-month-variable both out of
  # the pair 'year_start' & 'month_start' and 'year_end' & 
  # 'month_end'
  df.blo.start.end$yearmon_start <- as.yearmon(
    df.blo.start.end$my_start-0.083, "%Y-%B")
  df.blo.start.end$yearmon_end <- as.yearmon(
    df.blo.start.end$my_end-0.083, "%Y-%B")
  # next, we add a helper-column that we need in the 
  # next stage as well:
  df.blo.start.end$helper_date <- format(
    df.blo.start.end$yearmon_start-0.083, "%b")
  
  # this means that we can now drop the 'year' and 'month' 
  # variables
  df.blo.start.end <- df.blo.start.end %>%
    dplyr::select(-c(year_start, year_end, month_start, 
                     month_end))
  
  # next we add a column that gives us the duration of the 
  # shock:
  df.blo.start.end <- df.blo.start.end %>%
    mutate(duration = (df.blo.start.end$yearmon_end - 
                         df.blo.start.end$yearmon_start) * 12 + 1) %>%
    dplyr::select(-helper_date)
  
  # next we group_by shock_ID and extract the respective 
  # MAXIMUM
  # VOLATILITY
  # (note that the produced data-frame only has 
  # one row per shock_ID)
  df.blo.max <- as.data.frame(df.sub.blo %>%
                                filter(!is.na(shock_ID))  %>%
                                group_by(shock_ID) %>%
                                filter(series == max(series)) %>%
                                mutate(max_vol = series) %>%
                                dplyr::select(shock_ID, max_vol, my))
  
  # we apply the same date-transformation as above:
  df.blo.max$yearmon_max <- as.yearmon(df.blo.max$my-0.083, "%Y-%B")
  # and we drop 'my':
  df.blo.max <- df.blo.max %>%
    dplyr::select(-my)
  
  # in a last stage, we merge df.blo.max with df.blo.start.end:
  df.blo.start.end.max <- merge(x = df.blo.start.end, 
                                y = df.blo.max, 
                                by = "shock_ID", all = TRUE, na.rm=T)
  
  # the above data-frame now allows us to plot the episodes of 
  # high volatility into our previous plots (by adding a shaded
  # rectangle!)
  
  # one addition:
  # knowing that we also have the data.frame df.sub.first.vol
  # that contains the first first volatilites, we add this
  # info to the above data-frame:
  df.blo.start.end.max <- merge(x = df.blo.start.end.max, 
                                y = df.sub.first.vol, 
                                by = "shock_ID", 
                                all = TRUE, na.rm=T)
  
  # the column 'my' now represents the month with the first
  # volatilit:
  # we apply the same date-transformation as above:
  df.blo.start.end.max$yearmon_first <- as.yearmon(
    df.blo.start.end.max$my-0.083, "%Y-%B")
  # and we drop 'my':
  df.blo.start.end.max <- df.blo.start.end.max %>%
    dplyr::select(-c(my, first, my_start, my_end)) %>%
    dplyr::select(shock_ID, yearmon_start,
                  yearmon_end, duration, 
                  max_vol, first_vol, 
                  yearmon_max,
                  yearmon_first)
  
  # in a very last step, we re-initialize the generation
  # of the shock_IDs:
  # in a last step, we want count the ID-column from scratch again!
  df.blo.start.end.max <- 
    df.blo.start.end.max %>%
    mutate(shock_ID = group_indices(., 
                                    factor(shock_ID, levels = unique(shock_ID))))
  
  # in a very last step, we round all numbers to reasonable
  # number of decimals:
  df.blo.start.end.max[is.num] <- 
    lapply(df.blo.start.end.max[is.num], round, 2)
  
  # the below is only necessary if we would plot the above
  # information!
  # actually, there is a slight problem, which we fix by adding
  # one month's numeric value to my_end:
  # shocks_start_end$my_end <- shocks_start_end$my_end + (1/12)
  # this makes sure that if we refer to an end-month that
  # we assume that we are talking about the last day of that
  # respective month!
  
  
  # (6) return objects
  
  # lastly, we return the two data-frames that we want:
  # we want the two data-frames that are returned to 
  # get a 'suffix' according to the 'series'-string chosen
  # as the second argument of the function:
  
  # df1.out <- deparse(substitute(df.start.end))
  # df2.out <- deparse(substitute(df.sub))
  # 
  # df1.out <- paste(df1.out, series, sep="_")
  # df2.out <- paste(df2.out, series, sep="_")
  # 
  # df1.out <- df.start.end
  # df2.out <- df.sub
  
  # return
  return(list(df.sub, df.dates, df.start.end, 
              df.blo.start.end.max, df.sub.blo))
  
}          



## -------------------------------------
## Functions for calculation of summary statistics
## -------------------------------------

### Note: below is a customized function to calculate percentage
### changes of time series objects that adds an NA on top so that
### the resulting vector can be added back in into the respective
### structure (e.g., a data-frame)

### we have called the function 'perc2' because we already defined
### a perc-function somwhere else!
perc2 <- function(x){
  x <- round(diff(x) / x *100, 2)
  x <- c(NA, x)
  return(x)
}

