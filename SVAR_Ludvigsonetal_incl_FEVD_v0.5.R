#########################################################################################
### Marcel Kropp, 26.05.2020
### This script calculates the FEVD for the identified set following 
### Ludvigson et al (2018)'s algorithm.

### As outlined at the end of the script 'SVAR_Ludvigsonetal newnew v0.10.R',
### this script does not re-perform the entire algorithm for the generation 
## of the identified set, but rather continues where, in comparison to the 
## calculation of the IRFs in the script 'SVAR_Ludvigsonetal newnew v0.10.R'
## we need to calculate longer horizons for the Psi and Theta matrices!

## hence, below we re-perform only that part of the script 
## 'SVAR_Ludvigsonetal newnew v0.10.R', where we calculate the Psi and 
## Theta matrices and calculate them up until h = 201.
      
          # we initialize the list of Psi coefficient matrices
          psi_list <- list()
          
          # in our case length(B_list) = 6 (the number of lags in our VAR-system)
          # hence, we have to explicitly write down how the first six 
          # Psi's have to be calculated in our system:
          # Note that we have to store the first Psi (which is actually Psi 0),
          # at position 1; later (see below) we will rename the respective Psi's
          # stored in psi_list to match the technical notation in the thesis!
          
          # The approach to calculate the respective Psi's is as follows
          # (this follows directly from the recursion as written down in 
          # the thesis!)
          # Psi0 = identity matrix
          # Psi1 = B1
          # Psi2: B1 * Psi1 + B2 * Psi0
          # Psi3: B1 * Psi2 + B2 * Psi1 + B3 * Psi0
          # Psi4: B1 * Psi3 + B2 * Psi2 + B3 * Psi1 + B4 * Psi0
          # Psi5: B1 * Psi4 + B2 * Psi3 + B3 * Psi3 + B4 * Psi1 + B5 * Psi0
          # Psi6: B1 * Psi5 + B2 * Psi4 + B3 * Psi3 + B4 * Psi2 + B5 * Psi1 + B5 * Psi0
          
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
          
          # after having manually calculated the first 6 Psi-matrices,
          # for the remaining matrices the logic of the recursion remains 
          # the same; but because our VAR-system is a VAR(6)-system with 
          # 6 lags, there are no more B_matrices to take into account, hence
          # the procedure reduces to:
          # (1) B_block[, 1:18] contains all coefficients of the reduced-form VAR
          # which stays the same for every forecast horizon as of 
          # Psi_7 and is multiplied with a block-matrix that for every
          # forecast horizon contains the 
          # previous six Psi's stacked in a block-matrix from top (most recent)
          # to bottom:
          
          # Note that as of here we calculate Psi7 - Psi60 but the loop runs
          # from 8 to 61 because we had to store Psi0 on position 1 in the list
          # (there is no position 0 in lists!)
          
          # Marcel (26.05.2020):
          # Note that for the IRF-Analysis we use 60 steps, hence in the original
          # implementation of the below function we have let the loop run from
          # 8:61; but for the FEVD we want a longer forecast horizon (at least
          # until 100); therefore we let the loop run up until 101 here
          # to be able to extract a longer forecast horizon for the FEVD;
          # for the IRF-analysis the forecast horizon remains at h = 60!
          
          # for(i in 8:61){
          for(i in 8:201){
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
          

    # Marcel, 25.05.2020      
    ###-------------------------------
    ### Calculating the Forecast Error Variance Decomposition
    ###------------------------------- 
          
    # Having read a few ressources about the FEVD, we know how to
    # calculate the FEVD for different horizons;
    # we know that the list 'Theta' holds for each admissable A_0-matrix
    # the corresponding calculated Thetas at various forecast horizons
    # (here from 0 to 60);
    
    # hence, to calculate the MSPE for all admissable solutions in the 
    # list 'Thetas' at the various forecast horizons following Kilian/Lütkepohl
    # or Lütkepohl, respectively, we define the following function:
    
    # the function takes the huge list of Thetas and the forecast 
    # horizon h and the particular solution p for A_0 we are looking at:
          
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'mspe.matrix' was moved to the R-Script 
    # 20200628_functions.R  
    #---------------------------------------------------------------------------
          
    # running e.g. 
    # mspe.matrix(Thetas, 1, 1)
    # we know that this is the 1-step ahead MSE forecast matrix; 
    # as described in the Time Series Book of Lütkepohl, the diagonal
    # elements of this matrix are the MSEs (forecast error variance) of the
    # variables we are interested in!
          
    # if we now want to perform the forecast error variance DECOMPOSITION
    # we need a little more of matrix algebra and we are done:
    # in particular, on the p. 64 in the book of Lütkepohl, the algorithm
    # for calculating the respesctive denominators is described, consisting
    # of elementary vectors;
          
    # so, depending on the particular forecast horizon we are looking at, we
    # could create the following to calculate the respective decomposition:
    
    # this time the output will be a number, not matrix, like above!
    # here, we have to watch out, because if h = 1 (i.e. the one-step ahead
    # forecast), then the sum reduces to one single summand!
    # hence:
    
    # we also should set up a function that calculates the contribution
    # of a shock to the variance of the forecast error (i.e. the numerator
    # in the fractions below!)
          
    # this function will then allow us to calculate the contribution of a 
    # shock to the variance of the forecast error, depending at which
    # horizon we are looking at:
          
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'contr.forecast.err.mat' was moved to the R-Script 
    # 20200628_functions.R  
    #---------------------------------------------------------------------------      
    
    # having the functions 'contr.forecast.err.mat' and 
    # 'mspe.matrix' at our disposal and knowing that we want
    # to calculate the FEVD for h=1, h=12 and h=100 (i.e. infty)
    # for all Theta-sequences in Thetas, we have to embed everything
    # into a function;
    # in particular, the function should calculate the respective
    # FEVD for a particular horizon for ALL Thetas in the identified
    # set (which are stored in Thetas):
          
          
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'FEVD' was moved to the R-Script 
    # 20200628_functions.R  
    #---------------------------------------------------------------------------  
    
    # now we run the above function 'FEVD' for all Thetas (80) in our
    # Thetas-list for the respective forecast-horizons:
    # Marcel (29.05.2020): Because we actually want to transform the below
    # and create a loop that loops from 1:201 (the maximum forecast horizon
    # for which we have calculated the Thetas);
    
    # so instead of separate calls to the function FEVD, we loop from 1:201
    # and store the results (which are data.frames) in a list:
    # FEVD.h1 <- FEVD(Thetas, 1)
    # FEVD.h12 <- FEVD(Thetas, 12)
    # FEVD.h100 <- FEVD(Thetas, 100)
    FEVD.list <- list()
    for(h in 1:length(Thetas[[1]])){
      # the loop runs for each forecast horizon from h=1 to H=201!
      
      # we calculate the FEVD for all forecast horizons across all Thetas;
      # so the list FEVD.list will in the end hold a list of 200 data.frames
      # (because we have calculated the Thetas up until a forecast horizon of
      # 200!)
      # and each data.frame will hold the respective calculated FEVD for
      # all 9 constellations (Because we have 3 variables in the system)
      # for each corresponding A_0-solution! so if we have e.g. 80 A_0-solutions,
      # each data.frame will hold 9 columns (for the respective proportions, i.e.
      # the FEVDs) and the number of rows will be equal to the number of 
      # A-0-solutions (i.e. the number of different Thetas)!
      calculate.FEVD <- FEVD(Thetas, h)
      # 'h' stands for the forecast horizon!
      
      FEVD.list[[h]] <- calculate.FEVD
      names(FEVD.list)[length(FEVD.list)] <- paste0('forecast_horizon_', h)
    }
    
    # if we now as a quick test run
    FEVD.list[[1]]
    # we see that we have calculated the respective proportions 8
    # times because our solution set contains 8 A-0s and hence
    # a difference Theta-sequences for each forecast horizon!


    # with the results for the forecast-horizons at our disposal, we
    # want to proceed like Ludvigson et al (2018) and calculate the ranges
    # of the FEVD for the respective columns:
    # finding the range translates to finding the largest and the smallest
    # values of the respective data.frames;
    # hence, we proceed as follows to calculate the minimum
    # and maximum value for every column in the data.frames:
    
    # instead of performing the below function manually step by 
    # step for each forecast horizon, we should rather loop through 
    # all lists in FEVD.list and calculate the respective data.frames:
    
    # so, the below lops through all data.farmes in FEVD.list, and overwrites
    # the presentation of the information in the data.frames 
    
    # we create a new list to hold the results of the below loop:
    FEVD.min.max.list <- list()
    for(h in 1:length(FEVD.list)){
        # we create a data.frame in the desired shape
        FEVD.min.max <- cbind(
          sapply(FEVD.list[[h]], function(x) min(as.numeric(x))),
          sapply(FEVD.list[[h]], function(x) max(as.numeric(x)))
        )
        # and then we add it to our list:
        FEVD.min.max.list[[h]] <- FEVD.min.max
        # and give it the desired name:
        names(FEVD.min.max.list)[length(FEVD.min.max.list)] <- paste0('forecast_horizon_', h)
    }
    
    # having a look at
    FEVD.min.max.list[[1]]
    # we see that across all possible solutions for the proportions,
    # we have calculated the min and max value (i.e. instead of 
    # 8 different possibilities for the respective solutions because
    # we had 8 different A-0-solutions in our set), now we are only
    # left with ONE particular solution (i.e. min and max) for each 
    # proportion!
    
    
    # the format of the above matrix is perfect for further manipulation
    # to make it look like the table in Ludvigson et al 2018:
    
    # so far we had combined the results across all variables in the 
    # data-frames; now we start splitting, i.e., 
    # for each forecast error variance decomposition of each of the three
    # variables we now create a separate data.frame so that we can easily
    # append the results at the respective forecast horizons:
    
    # but because we want to add all solutions in the list FEVD.min.max 
    # to all the data.frames that we create below, we create three 
    # loops that loop through FEVD.min.max.list and extract the respective
    # data that we need at each stage (i.e. the min and the max)!
    
    # FEVD.macro
    FEVD.macro <- data.frame("U_m_Shock_min"=numeric(),
                             "U_m_Shock_max"=numeric(),
                             "ip_Shock_min"=numeric(),
                             "ip_Shock_max"=numeric(),
                             "U_f_Shock_min"=numeric(),
                             "U_f_Shock_max"=numeric(),
                             stringsAsFactors = FALSE)
    # looping through FEVD.min.max.list (to get all forecast-horizons
    # into the data.frame FEVD.macro);
    # we need to split the min and max for the respective shocks to be 
    # able to easily calculate the maximum afterwards (like in 
    # Ludvigson et al, 2018):
    for(h in 1:length(FEVD.min.max.list)){
      FEVD.macro[nrow(FEVD.macro) + 1,] <- c(FEVD.min.max.list[[h]][1, 1], 
                                             FEVD.min.max.list[[h]][1, 2],
                                             FEVD.min.max.list[[h]][2, 1],
                                             FEVD.min.max.list[[h]][2, 2],
                                             FEVD.min.max.list[[h]][3, 1],
                                             FEVD.min.max.list[[h]][3, 2]
      )
    }
    
    # then we do the same for industrial production and financial_uncertainty:
    # FEVD.ip
    FEVD.ip <- data.frame("U_m_Shock_min"=numeric(),
                             "U_m_Shock_max"=numeric(),
                             "ip_Shock_min"=numeric(),
                             "ip_Shock_max"=numeric(),
                             "U_f_Shock_min"=numeric(),
                             "U_f_Shock_max"=numeric(),
                             stringsAsFactors = FALSE)
    # looping through FEVD.min.max.list (to get all forecast-horizons
    # into the data.frame FEVD.macro);
    # we need to split the min and max for the respective shocks to be 
    # able to easily calculate the maximum afterwards (like in 
    # Ludvigson et al, 2018):
    for(h in 1:length(FEVD.min.max.list)){
      FEVD.ip[nrow(FEVD.ip) + 1,] <- c(FEVD.min.max.list[[h]][4, 1], 
                                       FEVD.min.max.list[[h]][4, 2],
                                       FEVD.min.max.list[[h]][5, 1],
                                       FEVD.min.max.list[[h]][5, 2],
                                       FEVD.min.max.list[[h]][6, 1],
                                       FEVD.min.max.list[[h]][6, 2]
      )
    }
    
    # FEVD.fin
    FEVD.fin <- data.frame("U_m_Shock_min"=numeric(),
                          "U_m_Shock_max"=numeric(),
                          "ip_Shock_min"=numeric(),
                          "ip_Shock_max"=numeric(),
                          "U_f_Shock_min"=numeric(),
                          "U_f_Shock_max"=numeric(),
                          stringsAsFactors = FALSE)
    # looping through FEVD.min.max.list (to get all forecast-horizons
    # into the data.frame FEVD.macro);
    # we need to split the min and max for the respective shocks to be 
    # able to easily calculate the maximum afterwards (like in 
    # Ludvigson et al, 2018):
    for(h in 1:length(FEVD.min.max.list)){
      FEVD.fin[nrow(FEVD.fin) + 1,] <- c(FEVD.min.max.list[[h]][7, 1], 
                                       FEVD.min.max.list[[h]][7, 2],
                                       FEVD.min.max.list[[h]][8, 1],
                                       FEVD.min.max.list[[h]][8, 2],
                                       FEVD.min.max.list[[h]][9, 1],
                                       FEVD.min.max.list[[h]][9, 2]
      )
    }
    
    # the below code has to be transformed into a function so that we 
    # can use it on FEVD.macro, FEVD.ip, and FEVD.fin:
    
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'extract.table.FEVD' was moved to the R-Script 
    # 20200628_functions.R  
    #--------------------------------------------------------------------------- 
    

    # now we can run the above function on all three data.frames
    # FEVD.macro, FEVD.ip and FEVD.fin and hopefully it works:
    FEVD.macro <- extract.table.FEVD(FEVD.macro)
    FEVD.ip <- extract.table.FEVD(FEVD.ip)
    FEVD.fin <- extract.table.FEVD(FEVD.fin)

    # now we combine the results by stacking the data.frames 
    # on top of each other:
    FEVD.all <- rbind(
                            FEVD.macro,
                            FEVD.ip,
                            FEVD.fin
    )
     
             
    # and now we can finally export the table as a latex-table:
    latextable(FEVD.all, 
               colnames = c("Um Shock", "ip Shock", "Uf Shock", "h"),
               rownames = c(
                            "1",
                            "2",
                            "3",
                            "4",
                            "5",
                            "12",
                            "infinity",
                            "h_max",
                            "1",
                            "2",
                            "3",
                            "4",
                            "5",
                            "12",
                            "infinity",
                            "h_max",
                            "1",
                            "2",
                            "3",
                            "4",
                            "5",
                            "12",
                            "infinity",
                            "h_max"))


      

      
