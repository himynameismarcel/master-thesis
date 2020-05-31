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
library(scales) # to get rid of scientific notation in the plots!

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
    
    for(k in 1:10000){ # 50000 führt zu 750 bestandenen Fällen
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
                  
                  # Marcel, 30.05.2020
                  # in a last step we rename the object A_0_valid so that we can
                  # differentiate it from the other sets created above:
                  A_0_valid.FE_CONSTR.II <- A_0_valid
                  # in a last step, we take the inverse of A_0_valid_no_constr
                  # to have the actual A_0_inv-matrices again:
                  A_0_valid.FE_CONSTR.II <- lapply(A_0_valid.FE_CONSTR.II, inv)
            
          }
          
    }
    

# having created A_0_validFE_CONSTR, we want to extract to sets of
# values: these are A_0_YM (which is the elemtn [2, 1]) and
# A_0_YF (Which is the element [2, 3]);
# to be able to extract these values, we have to run through the list
# and copy out the respective values and store them in a vector:

A_0_valid.FE_CONSTR.YM <- do.call(rbind, 
                                  lapply(A_0_valid.FE_CONSTR.II, function(x) x[1,2]))
colnames(A_0_valid.FE_CONSTR.YM) <- "YM"

A_0_valid.FE_CONSTR.YF <- do.call(rbind, 
                                  lapply(A_0_valid.FE_CONSTR.II, function(x) x[3,2]))
colnames(A_0_valid.FE_CONSTR.YF) <- "YF"
# lastly, we create a new matrix where we add these two column vectors together:
A_0_valid.FE_CONSTR.comb <- as_tibble(cbind(A_0_valid.FE_CONSTR.YM, A_0_valid.FE_CONSTR.YF))
    
    
 

###--------------------------------------------------------------------------------------------
### Plotting the respective values for A_0_YF and A_0_YM
###--------------------------------------------------------------------------------------------
    
    # to be able to use facet, we have to reshape all data.frames:
    A_0_valid.FE_CONSTR.tidy <- A_0_valid.FE_CONSTR.comb   %>%
                            gather(series, data)
    # and ultimately we add another column called NO_CONSTR
    # and muliply the data-series with 100:
    A_0_valid.FE_CONSTR.tidy <- A_0_valid.FE_CONSTR.tidy %>%
                            mutate(TYPE_CONSTR = "EVENT CONSTRAINTS ONLY") %>%
                            mutate(data = data*100)

    
    ggplot(data=A_0_valid.NO_CONSTR.tidy) + 
      geom_histogram(aes(x=data, y=(..count..)/sum(..count..),
                         fill="#288cc9"),
                     color="black", 
                     alpha = 0.8,
                     binwidth=0.1) + 
      geom_histogram(data=A_0_valid.FE_CONSTR.tidy,
                   aes(x=data, y=(..count..)/sum(..count..),
                       fill="#c92222"),
                   color="black", 
                   alpha = 0.8,
                   binwidth=0.04) +
      geom_histogram(data=A_0_valid.ALL_CONSTR.tidy,
                     aes(x=data, y=(..count..)/sum(..count..),
                         fill="black"),
                     color="black",  
                     alpha = 0.8,
                     binwidth=0.009) + 
      scale_x_continuous(name = NULL, 
                         breaks = seq(-0.6, 0.6, by = 0.2),
                         minor_breaks = NULL, labels=comma) + 
      scale_y_continuous(name = NULL) +
      # very useful blog-post about how to include and color 
      # a legend:
      # https://aosmith.rbind.io/2018/07/19/manual-legends-ggplot2/
      scale_fill_identity(guide = "legend",
                          name="",
                          labels=c("No Constraints", 
                                   "Event Constraints (w.o. Lehman Event)",
                                   "All Constraints")) +
      theme(legend.position = "top",aspect.ratio=0.5) +
      facet_wrap(.~series)

      ggsave(file="distribution_impact_matrices_type2.pdf")


