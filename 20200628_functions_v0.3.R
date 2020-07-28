## 28.06.2020
## This file contains all the functions that I have created and use
## in the course of the analyses for the Master's Thesis:

## Hence, before running any other script, this script should be executed
## before-hand.


## -------------------------------------
## Functions for Impulse Response Functions
## -------------------------------------



matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}  



LMN_algorithm <- function(restr_type, no_of_replications){
  ###-------------------------------
  ### Reading in Data
  ###-------------------------------
  
      # estimation window: 1960 - 2015:04 with CRSP-data
      # SVAR.data <- read_excel("20200718_data_LMN_incl_CRSP_201504_v0.2.xlsx",
      #                         sheet = "Tabelle1")
      
      # estimation window: 1960 - 2015:04 with CRSP-data
      SVAR.data <- read_excel("Replication_data_Ludvigson.xlsx",
                              sheet = "Data")
      
      # estimation window: 1960 - 2015:04 with S&P500-data
      # SVAR.data <- read_excel("20200718_data_LMN_incl_SP500_201504_v0.2.xlsx",
      #                         sheet = "Tabelle1")

      # # estimation window: 1960 - 2020:04
      # SVAR.data <- read_excel("20200718_data_LMN_incl_SP500_202004_v0.2.xlsx", 
      #                         sheet = "Tabelle1")
      
      # estimation window: 1960 - 2019:04
      # SVAR.data <- read_excel("20200718_data_LMN_incl_SP500_201904_v0.2.xlsx",
      #                         sheet = "Tabelle1")

  SVAR.data$yearmon <- as.yearmon(as.character(SVAR.data$Date), "%Y%m")
  
  
  SVAR.data.sub <- SVAR.data %>%
    dplyr::select(Um, ip, Uf)
  
  return.data <- SVAR.data %>% 
    dplyr::select(S, yearmon) %>%
    dplyr::mutate(S = S/100) %>%
    rename(S1 = S)
  
  
  external.data <- as_tibble(
    merge(x = return.data, 
          y = gold.data, 
          by = "yearmon", all = FALSE, na.rm=T)
          )
  
  ###-------------------------------
  ### Algorithm
  ###-------------------------------
  
  # set.seed(999)
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
  
  external.data <- tail(external.data,-6)
  
  # jetzt sollte es so sein, dass die return-series
  # S1 und S2 genau gleich lang sind und matching
  # positions mit den epsilons aus dem 
  # structural VAR!
  
  
  # we initialize A_0^{-1} as the lower triangular matrix P
  # (note: we have already created P above via
  # P <- t(chol(Sigma_u)))
  A_0_inv <- P
  
  # we further initialize the necessary parameters
  
  k1 <- 4.16 
  k2 <- 4.57
  k3 <- 4.73 
  k4 <- 4.05     
  k5 <- 2 
  
      # date-column for the structural errors for estimation-window
      # 1960 - 2015:
      dates = seq(from = as.Date("1961-01-01"),
                  to = as.Date("2015-04-01"), by = 'month')
      
      # date-column for the structural errors for estimation-window
      # 1960 - 2020:04: 
      # dates = seq(from = as.Date("1961-01-01"), 
      #             to = as.Date("2020-04-01"), by = 'month')
      
      # date-column for the structural errors for estimation-window
      # 1960 - 2019:04: 
      # dates = seq(from = as.Date("1961-01-01"), 
      #             to = as.Date("2019-04-01"), by = 'month')
  
  
  # to check how many draws actually pass all the tests,
  # we create the object x
  x <- 0
  
  n <- 3 # stands for the number of variables!
  
  # we further initialize a list for the matrices epsilon_t
  # that pass all constraints
  # and the matrices A_0 that generate those epsilon_t's that
  # pass all the tests
  epsilon_t_valid <- list()
  A_0_valid <- list()
  
  z <- no_of_replications
  
  ## the below has to run 1.5 million times:
  for(k in 1:z){
    
    #k <- 3
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
    
    M <- matrix(rnorm(n*n,mean=0,sd=1), n, n)
    
    QR <- qr(M)
    Q <- qr.Q(QR)    
    
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
    A_0_inv <- A_0_inv %*% Q
    
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
    
    # The above method, however, is very inefficient!
    # Therefore we want to try out a slightly smarter
    # procedure:
    epsilon_t <- t(apply(t(u_t), 2, function(row) A_0%*%row))
    
    
    # -------------------------------------------
    # EVENT CONSTRAINTS
    # -------------------------------------------
    # -------------------------------------------
    # PRELIMINARIES: EVENT CONSTRAINTS
    # -------------------------------------------
    epsilon_t <- as.data.frame(epsilon_t)
    
    # add a date-column which we have created outside of the loop
    epsilon_t$date <- dates # which has exactly length 652!
    # and we assign sensible col-names
    epsilon_t <- as_tibble(epsilon_t %>%
                             dplyr::rename(
                               macro_h1 = V1,
                               lip  = V2,
                               financial_h1 = V3))
    # we further split the 'Date' - variable into three pieces:
    epsilon_t <- separate(epsilon_t, 
                          "date", c("year", "month", "day"), 
                          sep = "-", 
                          remove=FALSE, convert=TRUE)
    # and we create the column 'yearmon'
    
    epsilon_t$yearmon <- as.yearmon(paste0(epsilon_t$year, 
                          epsilon_t$month), 
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
    
    
    # checking all conditions at once:
    if(restr_type == "FC_NN_RA_constr"){
      if(
        # other constraints:
        # real activity constraint
        # FE_4
        (sum(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                  yearmon <= "Jun 2009") %>%
             dplyr::select(lip) < k5) == 19) &
        # ----------------------------------------------------#
        # non-negative constraints
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
             dplyr::select(macro_h1) >= 0) == 2) &
        
        # CORRELATION CONSTRAINTS:
        # ----------------------------------------------------#
        # FC_1
        (c_M.S1 <= 0 & c_F.S1 <= 0) &
        # FC_2
        (c_M.S2 >= 0 & c_F.S2 >= 0)
      ){
        
        
        x <- x + 1
        print(x)
        
        #------------------------
        ## identified solution set and calculation of IRFs
        #------------------------
        
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
    }else if(restr_type == "FC_constr"){
      if(
        # CORRELATION CONSTRAINTS:
        # ----------------------------------------------------#
        # FC_1
        (c_M.S1 <= 0 & c_F.S1 <= 0) &
        # FC_2
        (c_M.S2 >= 0 & c_F.S2 >= 0)
      ){
        x <- x + 1
        print(x)
        
        
        #------------------------
        ## identified solution set and calculation of IRFs
        #------------------------
        
        epsilon_t_valid[[length(epsilon_t_valid)+1]] <- epsilon_t
        names(epsilon_t_valid)[length(epsilon_t_valid)] <- 
          paste0(c("iteration_"), k)
        
        A_0_valid[[length(A_0_valid)+1]] <- A_0
        names(A_0_valid)[length(A_0_valid)] <- 
          paste0(c("iteration_"), k)
        
        }
      }else if(restr_type == "NO_constr"){
        x <- x + 1
        print(x)
        
        #------------------------
        ## identified solution set and calculation of IRFs
        #------------------------
        
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
      }else if(restr_type == "FE_constr"){
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
      # other constraints:
        # real activity constraint
          # FE_4
          (sum(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                      yearmon <= "Jun 2009") %>%
                 dplyr::select(lip) < k5) == 19) &
          # ----------------------------------------------------#
        # non-negative constraints
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
                 dplyr::select(macro_h1) >= 0) == 2)
      ){
        x <- x + 1
        print(x)
        
        #------------------------
        ## identified solution set and calculation of IRFs
        #------------------------
        
        epsilon_t_valid[[length(epsilon_t_valid)+1]] <- epsilon_t
        names(epsilon_t_valid)[length(epsilon_t_valid)] <- 
          paste0(c("iteration_"), k)
        
        A_0_valid[[length(A_0_valid)+1]] <- A_0
        names(A_0_valid)[length(A_0_valid)] <- 
          paste0(c("iteration_"), k)
        }
      }else if(restr_type == "ALL_NO_LEHM_constr"){
        if(
          # EVENT CONSTRAINTS:
          # big-shock events:
          # FE_1: financial_uncertainty
          (epsilon_t %>% dplyr::filter(yearmon == "Oct 1987") %>%
           dplyr::select(financial_h1) >= k1) &
          # ----------------------------------------------------#
          # FE_2: Part 1: financial_uncertainty
          # ----------------------------------------------------#
          # FE_3: macro_uncertainty
          (epsilon_t %>% filter(yearmon == "Dec 1970") %>%
           dplyr::select(macro_h1) >= k4) &
          # ----------------------------------------------------#
          # other constraints:
          # real activity constraint
          # FE_4
          (sum(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                    yearmon <= "Jun 2009") %>%
               dplyr::select(lip) < k5) == 19) &
          # ----------------------------------------------------#
          # non-negative constraints
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
          # ----------------------------------------------------#
          # CORRELATION CONSTRAINTS:
          # ----------------------------------------------------#
          # FC_1
          (c_M.S1 <= 0 & c_F.S1 <= 0) &
          # FC_2
          (c_M.S2 >= 0 & c_F.S2 >= 0)
        ){
          x <- x + 1
          print(x)
          
          #------------------------
          ## identified solution set and calculation of IRFs
          #------------------------
          
          epsilon_t_valid[[length(epsilon_t_valid)+1]] <- epsilon_t
          names(epsilon_t_valid)[length(epsilon_t_valid)] <- 
            paste0(c("iteration_"), k)
          
          A_0_valid[[length(A_0_valid)+1]] <- A_0
          names(A_0_valid)[length(A_0_valid)] <- 
            paste0(c("iteration_"), k)
        }
      }else if(restr_type == "FC_FE_BIG_constr"){
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
          # CORRELATION CONSTRAINTS:
          # ----------------------------------------------------#
          # FC_1
          (c_M.S1 <= 0 & c_F.S1 <= 0) &
          # FC_2
          (c_M.S2 >= 0 & c_F.S2 >= 0)
        ){
          x <- x + 1
          print(x)
          
          #------------------------
          ## identified solution set and calculation of IRFs
          #------------------------
          
          epsilon_t_valid[[length(epsilon_t_valid)+1]] <- epsilon_t
          names(epsilon_t_valid)[length(epsilon_t_valid)] <- 
            paste0(c("iteration_"), k)
          
          A_0_valid[[length(A_0_valid)+1]] <- A_0
          names(A_0_valid)[length(A_0_valid)] <- 
            paste0(c("iteration_"), k)
        }
      }
    
    
    }
    return(list(A_0_valid, my.var, epsilon_t_valid))
}



# after the above test, we need to find an efficient way to 
# loop through the list of the above data.frames and extract
# min and max-values; the best way to do so would be to 
# use a function which we then also can use for the other
# scenarios where we need to extract max- and min-values:
extract.min.max <- function(df1, df2, df3, 
                            df4, df5, df6, 
                            df7, df8, df9){
  irfs.list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9)
  
  irfs.min.max.list <- list()
  # now we can loop through the list:
  for(z in 1:length(irfs.list)){
    
    # for each data.frame we extract min and max
    # first min:
    irfs.min <- apply(irfs.list[[z]], 1, FUN=min)
    # then max:
    irfs.max <- apply(irfs.list[[z]], 1, FUN=max)
    
    # next we combine the two matrices into a new 
    # data.frame:
    
    irfs.min.max.df <- as_tibble(cbind(irfs.min, irfs.max))
    
    # to be sure that we do not mix up the series
    # we decide for the same naming convention as
    # we use for the data.frames with the complete
    # amount of series:
    if(z == 1)
    {names(irfs.min.max.df) <- c("Fin_Macro.min", "Fin_Macro.max")}
    if(z == 2)
    {names(irfs.min.max.df) <- c("Fin_Ipm.min", "Fin_Ipm.max")}
    if(z == 3)
    {names(irfs.min.max.df) <- c("Fin_Fin.min", "Fin_Fin.max")}
    if(z == 4)
    {names(irfs.min.max.df) <- c("Macro_Macro.min", "Macro_Macro.max")}
    if(z == 5)
    {names(irfs.min.max.df) <- c("Macro_Ipm.min", "Macro_Ipm.max")}
    if(z == 6)
    {names(irfs.min.max.df) <- c("Macro_Fin.min", "Macro_Fin.max")}
    if(z == 7)
    {names(irfs.min.max.df) <- c("Ipm_Macro.min", "Ipm_Macro.max")}
    if(z == 8)
    {names(irfs.min.max.df) <- c("Ipm_Ipm.min", "Ipm_Ipm.max")}
    if(z == 9)
    {names(irfs.min.max.df) <- c("Ipm_Fin.min", "Ipm_Fin.max")}
    
    
    # lastly, we append the newly created data.frame
    # to the list irfs.min.max.list.
    irfs.min.max.list[[z]] <- irfs.min.max.df
    # at this stage, we also add a step-counter to each
    # data-frame:
    irfs.min.max.list[[z]]$step <- seq.int(nrow(irfs.min.max.list[[z]]))
    
  }    
  
  # after the loop we return the newly created list
  # that contains all the data.frames:
  return(irfs.min.max.list)
}


generate.irfs.plots <- function(df1, df2, df3, 
                                df4, df5, df6, 
                                df7, df8, df9){
  bind_rows(
    df1 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Macro Uncertainty",
                    response = "Macro Uncertainty"),
    df2 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Macro Uncertainty",
                    response = "Ind. Production"),
    df3 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Macro Uncertainty",
                    response = "Fin. Uncertainty"), 
    df4 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Ind. Production",
                    response = "Macro Uncertainty"), 
    df5 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Ind. Production",
                    response = "Ind. Production"),
    df6 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Ind. Production",
                    response = "Fin. Uncertainty"), 
    df7 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Fin. Uncertainty",
                    response = "Macro Uncertainty"),
    df8 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Fin. Uncertainty",
                    response = "Ind. Production"),
    df9 %>%
      gather(
        key = series_name,
        value = series_value,
        -c(step)) %>%
      # at the same time we add the names for impulse and response
      dplyr::mutate(impulse = "Fin. Uncertainty",
                    response = "Fin. Uncertainty")
  )
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

