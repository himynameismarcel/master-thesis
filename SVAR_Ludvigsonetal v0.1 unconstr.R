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


###--------------------------------------------------------------------------------------------
### Algorithm Only With EVENT CONSTRAINTS
###--------------------------------------------------------------------------------------------

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
library(Hmisc)

###-------------------------------
### Preliminary Explanations
###-------------------------------
###-------------------------------
### Reading in Data
###-------------------------------
SVAR.data <- read_excel("Replication_data_Ludvigson.xlsx", 
                        sheet = "Data")

SVAR.data$yearmon <- as.yearmon(as.character(SVAR.data$Date), "%Y%m")

SVAR.data.sub <- SVAR.data %>%
                        dplyr::select(Um, ip, Uf)

return.data <- SVAR.data %>% 
                     dplyr::select(S, yearmon) %>%
                     dplyr::mutate(S = S/100)


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
    my.var <- VAR(SVAR.data.sub, type = "const", p = 6)
    Sigma_u <- summary(my.var)$covres * 633/652
    u_t <- residuals(my.var)
    P <- t(chol(Sigma_u))
    
    A_0_inv <- P

    lambda1 <- -0.05
    lambda2 <- 2
    lambda3 <- 0.18
    k1 <- 4
    k2 <- 4
    k3 <- 2
    
    dates = seq(from = as.Date("1961-01-01"), 
                to = as.Date("2015-04-01"), by = 'month')
                
    x <- 0
    n <- 3
    
    epsilon_t_valid <- list()
    A_0_valid <- list()
    
    for(k in 1:50000){ # 50000 führt zu 750 bestandenen Fällen
          # print(k)
      
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
                R.sign[upper.tri(R.sign)] <- 0
                Q <-Q %*% R.sign
                
      
          ##########
          ## STEP 3:
          ## calculation of a possible impact matrix A_0^{-1} = P*Q:
          ##########       
          A_0_inv <- A_0_inv %*% Q
          A_0_inv.sign <- sign(A_0_inv)
          A_0_inv.sign[upper.tri(A_0_inv.sign)] <- 0
          A_0_inv.sign[lower.tri(A_0_inv.sign)] <- 0
        
          A_0_inv <-A_0_inv %*% A_0_inv.sign
          
          A_0 <- inv(A_0_inv)
    
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
          epsilon_t <- t(apply(t(u_t), 2, function(row) A_0%*%row))
          
          # -------------------------------------------
          # EVENT CONSTRAINTS
          # -------------------------------------------
          # -------------------------------------------
          # PRELIMINARIES: EVENT CONSTRAINTS
          # -------------------------------------------
          epsilon_t <- as.data.frame(epsilon_t)

          epsilon_t$date <- dates
          epsilon_t <- as_tibble(epsilon_t %>%
                dplyr::rename(macro_h1 = V1,
                              lip  = V2,
                              financial_h1 = V3))
          epsilon_t <- separate(epsilon_t, "date", c("year", "month", "day"), 
                              sep = "-", 
                              remove=FALSE, convert=TRUE)
          
          epsilon_t$yearmon <- as.yearmon(paste0(epsilon_t$year, 
                                                 epsilon_t$month), 
                                          "%Y %m")

          epsilon_t <- epsilon_t %>%
                      dplyr::select(-c(date, year, month, day))


          if(
              # FE_1
              (epsilon_t %>% dplyr::filter(yearmon == "Oct 1987") %>%
                                dplyr::select(financial_h1) > k1) & 
              # FE_2
              (any(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                yearmon <= "Jun 2009") %>%
                                dplyr::select(financial_h1) > k2)) &
              # FE_3
              (sum(epsilon_t %>% filter(yearmon >= "Dec 2007" &
                                yearmon <= "Jun 2009") %>%
                                dplyr::select(lip) < k3) == 19)
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
          
    }
    
#######################
### Figure 4 (Replication of Ludvigson et al (2018));
### Calculation of Impulse Response-Functions for the
### set only containing the EVENT CONSTRAINTS
#######################
    
    l <- 3 # number of equations
    p <- 6 # number of lags

    B_block <- t(sapply(my.var$varresult, `[[`, 1))
    matsplitter<-function(M, r, c) {
      rg <- (row(M)-1)%/%r+1
      cg <- (col(M)-1)%/%c+1
      rci <- (rg-1)*max(cg) + cg
      N <- prod(dim(M))/r/c
      cv <- unlist(lapply(1:N, function(x) M[rci==x]))
      dim(cv)<-c(r,c,N)
      cv
    }       

    B_array <- matsplitter(B_block, 3, 3)
    B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
    constant <- as.matrix(B_block[, ncol(B_block)])
    
    ## STEP 2: 
    ## in our initial approach our idea was to 
    ## (1) first recursively calculate 
        ## the corresponding matrices of the reduced-form MA(infinity)-
        ## representation (the psi-matrices) for each position 1:60 and
    ## (2) then later (see below)
        ## multiply all the psi-matrices of the MA(infinity)-representation
        ## of the reduced-form VAR at each step with the impact matrix inv(A_0)
        ## for each A_0 which is part of the identified solution set to 
        ## calculate the Theta-coefficients of the MA(infinity)-representation
        ## of the structural VAR (SVAR) at each step which then hold the 
        ## values of the impulse responses (IRFs) at the respective 
        ## forecast horizon;
          psi_list <- list()
          
          # Psi0:
          psi_list[[1]] <- diag(3)
          # Psi1:
          psi_list[[2]] <- B_list[[1]]
          # Psi2:
          psi_list[[3]] <- cbind(B_list[[1]], B_list[[2]]) %*%
                            rbind(psi_list[[2]], psi_list[[1]])
          # Psi3:
          psi_list[[4]] <- cbind(B_list[[1]], B_list[[2]], 
                                 B_list[[3]]) %*%
                            rbind(psi_list[[3]], psi_list[[2]], 
                                  psi_list[[1]])  
          # Psi4:
          psi_list[[5]] <- cbind(B_list[[1]], B_list[[2]], 
                                 B_list[[3]], B_list[[4]]) %*%
                            rbind(psi_list[[4]], psi_list[[3]], 
                                  psi_list[[2]], psi_list[[1]]) 
          # Psi5:
          psi_list[[6]] <- cbind(B_list[[1]], B_list[[2]], 
                                 B_list[[3]], B_list[[4]],
                                 B_list[[5]]) %*%
                            rbind(psi_list[[5]], psi_list[[4]], 
                                  psi_list[[3]], psi_list[[2]],
                                  psi_list[[1]])  
          # Psi6:
          psi_list[[7]] <- cbind(B_list[[1]], B_list[[2]], 
                                 B_list[[3]], B_list[[4]],
                                 B_list[[5]], B_list[[6]]) %*%
                            rbind(psi_list[[6]], psi_list[[5]], 
                                  psi_list[[4]], psi_list[[3]],
                                  psi_list[[2]], psi_list[[1]])

          for(i in 8:61){
            psi_list[[i]] <- B_block[, 1:18] %*%
                              rbind(psi_list[[i-1]], psi_list[[i-2]], 
                                    psi_list[[i-3]], psi_list[[i-4]],
                                    psi_list[[i-5]], psi_list[[i-6]])
          }
          
          # at last, we rename the matrices:
          for(i in 0:length(psi_list)-1){
            names(psi_list)[i+1] <- 
              paste0(c("psi_"), i)            
          }
    
    ## STEP 3: having the list 'psi_list' at our
    ## disposal, we can finally calculate the Thetas
    ## (which are the coefficients of the structural 
    ## MA-representation and which will then finally hold
    ## the coefficients for the impulse-response-functions
    ## that we want/need for plotting below!)
    #----------------------------------------------
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
    extract.coefs.irf <- function(impulse, response, row, col){
    
        # we initialize a data-frame:
        coefs.df <- data.frame(matrix(ncol = length(Thetas), 
                                      nrow = length(psi_list)))
        
        # # and change its name
        # name_df <- coefs.df
    
        for(m in 1:length(Thetas)){
          for(b in 1:61){

            coefs.df[b, m] <-  sapply(lapply(Thetas, "[[", b), 
                                      function(x) x[row, col])[m]
            colnames(coefs.df)[m] <- paste0(impulse, "_", response, "_", m)
          
          }
        }
        
        return(as_tibble(coefs.df))
    }

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
    plot_SVAR.irfs.all.FE_only <- plot_SVAR.irfs.all %>%
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

# remove all but one object from work-space    
# rm(list=setdiff(ls(), c("plot_SVAR.irfs.all.FE_only",
#                         "plot_SVAR.irfs.maxG.CONSTR_ALL")))
#     
###--------------------------------------------------------------------------------------------
### Algorithm Only With CORRELATION CONSTRAINTS
###--------------------------------------------------------------------------------------------

###-------------------------------
### Reading in Data
###-------------------------------
SVAR.data <- read_excel("Replication_data_Ludvigson.xlsx", 
                        sheet = "Data")

SVAR.data$yearmon <- as.yearmon(as.character(SVAR.data$Date), "%Y%m")

SVAR.data.sub <- SVAR.data %>%
  dplyr::select(Um, ip, Uf)


return.data <- SVAR.data %>% 
  dplyr::select(S, yearmon) %>%
  dplyr::mutate(S = S/100)


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

my.var <- VAR(SVAR.data.sub, type = "const", p = 6)

Sigma_u <- summary(my.var)$covres * 633/652

u_t <- residuals(my.var)

P <- t(chol(Sigma_u))


# -----------------------------------------------
# CORRELATION CONSTRAINTS: PRELIMINARIES
# -----------------------------------------------
return.data$S_minus_1 <- Hmisc::Lag(return.data$S, -1)

return.data <- head(return.data, -1)

ar.1 <- lm(S_minus_1 ~ S, return.data)

u_St <- residuals(ar.1)


u_St <- as.ts(tail(as.zoo(u_St), -5))

A_0_inv <- P

lambda1 <- -0.05
lambda2 <- 2
lambda3 <- 0.18
k1 <- 4
k2 <- 4
k3 <- 2

dates = seq(from = as.Date("1961-01-01"), 
            to = as.Date("2015-04-01"), by = 'month')

x <- 0

n <- 3

epsilon_t_valid <- list()
A_0_valid <- list()

for(k in 1:10000){ # 10000 führt zu 74 bestandenen Fällen
  # rint(k)

  M <- matrix(rnorm(n*n,mean=0,sd=1), n, n)
  QR <- qr(M)
  Q <- qr.Q(QR)    
  R.sign <- sign(qr.R(QR))
  R.sign[upper.tri(R.sign)] <- 0
  Q <-Q %*% R.sign

  A_0_inv <- A_0_inv %*% Q
  A_0_inv.sign <- sign(A_0_inv)
  A_0_inv.sign[upper.tri(A_0_inv.sign)] <- 0
  A_0_inv.sign[lower.tri(A_0_inv.sign)] <- 0
  

  A_0_inv <-A_0_inv %*% A_0_inv.sign
  
  A_0 <- inv(A_0_inv)

  
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
  epsilon_t <- t(apply(t(u_t), 2, function(row) A_0%*%row))

  # -------------------------------------------
  # EVENT CONSTRAINTS
  # -------------------------------------------

  # -------------------------------------------
  # PRELIMINARIES: EVENT CONSTRAINTS
  # -------------------------------------------
  epsilon_t <- as.data.frame(epsilon_t)
  

  epsilon_t$date <- dates
  epsilon_t <- as_tibble(epsilon_t %>%
                           dplyr::rename(macro_h1 = V1,
                                         lip  = V2,
                                         financial_h1 = V3))
  epsilon_t <- separate(epsilon_t, "date", c("year", "month", "day"), 
                        sep = "-", 
                        remove=FALSE, convert=TRUE)
  
  epsilon_t$yearmon <- as.yearmon(paste0(epsilon_t$year, epsilon_t$month), 
                                  "%Y %m")
  epsilon_t <- epsilon_t %>%
    dplyr::select(-c(date, year, month, day))
  
  
  # -------------------------------------------
  # CORRELATION CONSTRAINTS: PRELIMINARIES
  # -------------------------------------------
  u_t <- as_tibble(u_t)

  c_M <- cor(u_St, epsilon_t$macro_h1)

  c_F <- cor(u_St, (epsilon_t$financial_h1))
  c_MF <- sqrt(c_M * c_M + c_F * c_F)
  u_t <- as.matrix(u_t)
  if(
    # FC_1
    ((lambda1 - c_M) > 0 & (lambda1 - c_F) > 0) &
    # FC_2
    ((abs(c_F) - lambda2*abs(c_M)) > 0) &
    # FC_3
    (c_MF - lambda3 > 0)
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

#######################
### Figure 4 (Replication of Ludvigson et al (2018));
### Calculation of Impulse Response-Functions
#######################

## STEP 1: extracting reduced-form matrices from
## VAR-estimation from above
#----------------------------------------------

l <- 3 # number of equations
p <- 6 # number of lags

B_block <- t(sapply(my.var$varresult, `[[`, 1))

matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}       


B_array <- matsplitter(B_block, 3, 3)
B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
constant <- as.matrix(B_block[, ncol(B_block)])

## STEP 2: 
## in our initial approach our idea was to 
## (1) first recursively calculate 
## the corresponding matrices of the reduced-form MA(infinity)-
## representation (the psi-matrices) for each position 1:60 and
## (2) then later (see below)
## multiply all the psi-matrices of the MA(infinity)-representation
## of the reduced-form VAR at each step with the impact matrix inv(A_0)
## for each A_0 which is part of the identified solution set to 
## calculate the Theta-coefficients of the MA(infinity)-representation
## of the structural VAR (SVAR) at each step which then hold the 
## values of the impulse responses (IRFs) at the respective 
## forecast horizon;

psi_list <- list()

# Psi0:
psi_list[[1]] <- diag(3)
# Psi1:
psi_list[[2]] <- B_list[[1]]
# Psi2:
psi_list[[3]] <- cbind(B_list[[1]], B_list[[2]]) %*%
  rbind(psi_list[[2]], psi_list[[1]])
# Psi3:
psi_list[[4]] <- cbind(B_list[[1]], B_list[[2]], 
                       B_list[[3]]) %*%
  rbind(psi_list[[3]], psi_list[[2]], 
        psi_list[[1]])  
# Psi4:
psi_list[[5]] <- cbind(B_list[[1]], B_list[[2]], 
                       B_list[[3]], B_list[[4]]) %*%
  rbind(psi_list[[4]], psi_list[[3]], 
        psi_list[[2]], psi_list[[1]]) 
# Psi5:
psi_list[[6]] <- cbind(B_list[[1]], B_list[[2]], 
                       B_list[[3]], B_list[[4]],
                       B_list[[5]]) %*%
  rbind(psi_list[[5]], psi_list[[4]], 
        psi_list[[3]], psi_list[[2]],
        psi_list[[1]])  
# Psi6:
psi_list[[7]] <- cbind(B_list[[1]], B_list[[2]], 
                       B_list[[3]], B_list[[4]],
                       B_list[[5]], B_list[[6]]) %*%
  rbind(psi_list[[6]], psi_list[[5]], 
        psi_list[[4]], psi_list[[3]],
        psi_list[[2]], psi_list[[1]])

for(i in 8:61){
  psi_list[[i]] <- B_block[, 1:18] %*%
    rbind(psi_list[[i-1]], psi_list[[i-2]], 
          psi_list[[i-3]], psi_list[[i-4]],
          psi_list[[i-5]], psi_list[[i-6]])
}

# at last, we rename the matrices:
for(i in 0:length(psi_list)-1){
  names(psi_list)[i+1] <- 
    paste0(c("psi_"), i)            
}

## STEP 3: having the list 'psi_list' at our
## disposal, we can finally calculate the Thetas
## (which are the coefficients of the structural 
## MA-representation and which will then finally hold
## the coefficients for the impulse-response-functions
## that we want/need for plotting below!)
#----------------------------------------------
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

extract.coefs.irf <- function(impulse, response, row, col){
  
  # we initialize a data-frame:
  coefs.df <- data.frame(matrix(ncol = length(Thetas), 
                                nrow = length(psi_list)))

  
  for(m in 1:length(Thetas)){
    
    for(b in 1:61){
      coefs.df[b, m] <-  sapply(lapply(Thetas, "[[", b), 
                                function(x) x[row, col])[m]
      colnames(coefs.df)[m] <- paste0(impulse, "_", response, "_", m)
      
    }
  }
  
  return(as_tibble(coefs.df))
}

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
plot_SVAR.irfs.all.FC_only <- plot_SVAR.irfs.all %>%
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


# remove all but two object from work-space    
# rm(list=setdiff(ls(), c("plot_SVAR.irfs.all.FE_only",
#                         "plot_SVAR.irfs.all.FC_only",
#                         "plot_SVAR.irfs.maxG.CONSTR_ALL")))

###--------------------------------------------------------------------------------------------
### Algorithm Only With NO CONSTRAINTS AT ALL
###--------------------------------------------------------------------------------------------

###-------------------------------
### Preliminary Explanations
###-------------------------------
###-------------------------------
### Reading in Data
###-------------------------------
SVAR.data <- read_excel("Replication_data_Ludvigson.xlsx", 
                        sheet = "Data")

SVAR.data$yearmon <- as.yearmon(as.character(SVAR.data$Date), "%Y%m")

SVAR.data.sub <- SVAR.data %>%
  dplyr::select(Um, ip, Uf)

return.data <- SVAR.data %>% 
  dplyr::select(S, yearmon) %>%
  dplyr::mutate(S = S/100)


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
    my.var <- VAR(SVAR.data.sub, type = "const", p = 6)
    Sigma_u <- summary(my.var)$covres * 633/652
    u_t <- residuals(my.var)
    P <- t(chol(Sigma_u))
    
    A_0_inv <- P
    
    lambda1 <- -0.05
    lambda2 <- 2
    lambda3 <- 0.18
    k1 <- 4
    k2 <- 4
    k3 <- 2
    
    dates = seq(from = as.Date("1961-01-01"), 
                to = as.Date("2015-04-01"), by = 'month')
    
    x <- 0
    
    n <- 3
    
    epsilon_t_valid <- list()
    A_0_valid <- list()
    
    for(k in 1:700){
      # print(k)
      
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
      R.sign[upper.tri(R.sign)] <- 0
      Q <-Q %*% R.sign
      
      
      ##########
      ## STEP 3:
      ## calculation of a possible impact matrix A_0^{-1} = P*Q:
      ##########       
      A_0_inv <- A_0_inv %*% Q
      A_0_inv.sign <- sign(A_0_inv)
      A_0_inv.sign[upper.tri(A_0_inv.sign)] <- 0
      A_0_inv.sign[lower.tri(A_0_inv.sign)] <- 0
      
      A_0_inv <-A_0_inv %*% A_0_inv.sign
      
      A_0 <- inv(A_0_inv)
      
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
      epsilon_t <- t(apply(t(u_t), 2, function(row) A_0%*%row))
      
      # -------------------------------------------
      # EVENT CONSTRAINTS
      # -------------------------------------------
      # -------------------------------------------
      # PRELIMINARIES: EVENT CONSTRAINTS
      # -------------------------------------------
      epsilon_t <- as.data.frame(epsilon_t)
      
      epsilon_t$date <- dates
      epsilon_t <- as_tibble(epsilon_t %>%
                               dplyr::rename(macro_h1 = V1,
                                             lip  = V2,
                                             financial_h1 = V3))
      epsilon_t <- separate(epsilon_t, "date", c("year", "month", "day"), 
                            sep = "-", 
                            remove=FALSE, convert=TRUE)
      
      epsilon_t$yearmon <- as.yearmon(paste0(epsilon_t$year, 
                                             epsilon_t$month), 
                                      "%Y %m")
      
      epsilon_t <- epsilon_t %>%
        dplyr::select(-c(date, year, month, day))
        
        
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
    
    #######################
    ### Figure 4 (Replication of Ludvigson et al (2018));
    ### Calculation of Impulse Response-Functions for the
    ### set only containing the EVENT CONSTRAINTS
    #######################
    
    l <- 3 # number of equations
    p <- 6 # number of lags
    
    B_block <- t(sapply(my.var$varresult, `[[`, 1))
    matsplitter<-function(M, r, c) {
      rg <- (row(M)-1)%/%r+1
      cg <- (col(M)-1)%/%c+1
      rci <- (rg-1)*max(cg) + cg
      N <- prod(dim(M))/r/c
      cv <- unlist(lapply(1:N, function(x) M[rci==x]))
      dim(cv)<-c(r,c,N)
      cv
    }       
    
    B_array <- matsplitter(B_block, 3, 3)
    B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
    constant <- as.matrix(B_block[, ncol(B_block)])
    
    ## STEP 2: 
    ## in our initial approach our idea was to 
    ## (1) first recursively calculate 
    ## the corresponding matrices of the reduced-form MA(infinity)-
    ## representation (the psi-matrices) for each position 1:60 and
    ## (2) then later (see below)
    ## multiply all the psi-matrices of the MA(infinity)-representation
    ## of the reduced-form VAR at each step with the impact matrix inv(A_0)
    ## for each A_0 which is part of the identified solution set to 
    ## calculate the Theta-coefficients of the MA(infinity)-representation
    ## of the structural VAR (SVAR) at each step which then hold the 
    ## values of the impulse responses (IRFs) at the respective 
    ## forecast horizon;
    psi_list <- list()
    
    # Psi0:
    psi_list[[1]] <- diag(3)
    # Psi1:
    psi_list[[2]] <- B_list[[1]]
    # Psi2:
    psi_list[[3]] <- cbind(B_list[[1]], B_list[[2]]) %*%
      rbind(psi_list[[2]], psi_list[[1]])
    # Psi3:
    psi_list[[4]] <- cbind(B_list[[1]], B_list[[2]], 
                           B_list[[3]]) %*%
      rbind(psi_list[[3]], psi_list[[2]], 
            psi_list[[1]])  
    # Psi4:
    psi_list[[5]] <- cbind(B_list[[1]], B_list[[2]], 
                           B_list[[3]], B_list[[4]]) %*%
      rbind(psi_list[[4]], psi_list[[3]], 
            psi_list[[2]], psi_list[[1]]) 
    # Psi5:
    psi_list[[6]] <- cbind(B_list[[1]], B_list[[2]], 
                           B_list[[3]], B_list[[4]],
                           B_list[[5]]) %*%
      rbind(psi_list[[5]], psi_list[[4]], 
            psi_list[[3]], psi_list[[2]],
            psi_list[[1]])  
    # Psi6:
    psi_list[[7]] <- cbind(B_list[[1]], B_list[[2]], 
                           B_list[[3]], B_list[[4]],
                           B_list[[5]], B_list[[6]]) %*%
      rbind(psi_list[[6]], psi_list[[5]], 
            psi_list[[4]], psi_list[[3]],
            psi_list[[2]], psi_list[[1]])
    
    for(i in 8:61){
      psi_list[[i]] <- B_block[, 1:18] %*%
        rbind(psi_list[[i-1]], psi_list[[i-2]], 
              psi_list[[i-3]], psi_list[[i-4]],
              psi_list[[i-5]], psi_list[[i-6]])
    }
    
    # at last, we rename the matrices:
    for(i in 0:length(psi_list)-1){
      names(psi_list)[i+1] <- 
        paste0(c("psi_"), i)            
    }
    
    ## STEP 3: having the list 'psi_list' at our
    ## disposal, we can finally calculate the Thetas
    ## (which are the coefficients of the structural 
    ## MA-representation and which will then finally hold
    ## the coefficients for the impulse-response-functions
    ## that we want/need for plotting below!)
    #----------------------------------------------
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
    extract.coefs.irf <- function(impulse, response, row, col){
      
      # we initialize a data-frame:
      coefs.df <- data.frame(matrix(ncol = length(Thetas), 
                                    nrow = length(psi_list)))
      
      # # and change its name
      # name_df <- coefs.df
      
      for(m in 1:length(Thetas)){
        for(b in 1:61){
          
          coefs.df[b, m] <-  sapply(lapply(Thetas, "[[", b), 
                                    function(x) x[row, col])[m]
          colnames(coefs.df)[m] <- paste0(impulse, "_", response, "_", m)
          
        }
      }
      
      return(as_tibble(coefs.df))
    }
    
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
    plot_SVAR.irfs.all.NO_CONSTR <- plot_SVAR.irfs.all %>%
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

# remove all but three objects from workspace
# rm(list=setdiff(ls(), c("plot_SVAR.irfs.all.FE_only",
#                         "plot_SVAR.irfs.all.FC_only",
#                         "plot_SVAR.irfs.all.NO_CONSTR"
#                         )))

## before plotting, we have to add one more column to each of the data
## frames to be able to use that information for the aes in ggplot:
plot_SVAR.irfs.all.FE_only <- plot_SVAR.irfs.all.FE_only %>%
                          mutate(type = "FE_only")
plot_SVAR.irfs.all.FC_only <- plot_SVAR.irfs.all.FC_only %>%
                          mutate(type = "FC_only")                
plot_SVAR.irfs.all.NO_CONSTR <- plot_SVAR.irfs.all.NO_CONSTR %>%
                          mutate(type = "NO_CONSTR") 
plot_SVAR.irfs.all.CONSTR_ALL <- plot_SVAR.irfs.all.CONSTR_ALL %>%
                          mutate(type = "CONSTR_ALL")  
###----------------------------------------------------------------------------------------
### finally we plot the impulse responses:  
###----------------------------------------------------------------------------------------    
    impulse.responses_all.SVAR.COMPLETE <- 
      ggplot(data=plot_SVAR.irfs.all.NO_CONSTR, aes(x=step, y=series_value)) + 
      geom_line(aes(group=series_name, colour=type), alpha=0.01, size=0.9) +
  
      geom_line(data=plot_SVAR.irfs.all.FC_only, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 0.4, size = 0.3) +
      geom_line(data=plot_SVAR.irfs.all.FE_only, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 0.1, size = 0.1) +
      geom_line(data=plot_SVAR.irfs.all.CONSTR_ALL, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 1, size = 0.1) +  
      facet_wrap(impulse ~ response,
             scales="free_y", strip.position = "top") +
      labs(color=NULL) + 
      geom_hline(yintercept=0, color = "#514e4e", 
                 size=1) +
      scale_x_continuous(name = NULL) + 
      scale_y_continuous(name = NULL) +
      theme(axis.text=element_text(size=8),
            plot.title = element_text(size=10, face="bold", hjust = 0.5),
            axis.title=element_text(size=10),
            legend.position="top",
            #legend.text=element_text(size=14),
            #axis.text.x=element_blank(),
            plot.margin = unit(c(1,1,1,1), "mm"),
            #panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            aspect.ratio = 0.9,
            strip.text = element_text(size = 8, margin = margin(0.7, 0.7, 0.7, 0.7, "mm")),
            strip.background = element_rect(colour="black", 
                                            fill="white", 
                                            size=0.5, 
                                            linetype="solid")) + 
      coord_cartesian(ylim = c(-3, 3)) + 
      guides(colour = guide_legend(override.aes = list(size = 5))) + 
      #https://stackoverflow.com/questions/35712062/scale-fill-discrete-does-not-change-label-names
      scale_colour_discrete(
                      labels = c("All Constraints", "Correlation Constraints Only", 
                                 "Event Constraints Only", "No Constraints"))

    
      impulse.responses_all.SVAR.COMPLETE
       
      ggsave(file="impulse_responses_all_SVAR_unconstr_constr.pdf")
      
       

