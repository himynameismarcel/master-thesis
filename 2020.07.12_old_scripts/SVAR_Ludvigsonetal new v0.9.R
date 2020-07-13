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
    
    # Marcel (08.05.2020): Because we only need the coefficients and not
    # the constant, we remove the 
    # B_block <- B_block[, 1:18]
    # Note that the above stept is not necessary, because the below
    # automatically only saves the coefficients in B_block and the constant
    # somewhere else!
    
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
    
    
    # finally we store the B_i's to an array
    B_array <- matsplitter(B_block, 3, 3)
    # and convert the output to a list
    B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
    # and the vector of constants is stored in the matrix (a column vector)
    # 'constant'
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
    
    ## This approach turned out to be terribly inefficient for two reasons:
    ## (1) We would first calculate the Psi-matrices (the coefficients of the 
        ## reduced-form MA(infinity) representation) and THEN multiply
        ## the resulting Psi-matrices at each forecast horizon with the inv(A_0)
        ## of the solution-set and (which is even more important!)
    ## (2) In the recursion we would not make use of intermediate results
        ## but calculate the entire expression for the Psi-matrix for the 
        ## respective position from scratch for every forecast horizon;
        ## This approach turned out to be terribly inefficient, because it would
        ## for each position (1:60) start with the recursion from step = 0 
        ## again without making use of previous results that have already 
        ## been calculated along the way;
        ## With this approach, at a certain point the recursions to calculate
        ## the corresponding matrix at the respective forecast horizon became so 
        ## involved and deep that R was not performing well and the loop 
        ## was running for hours!
  
    ## Therefore:
    ## The below approach 
    ## (1) gets rid of the intermediate step below of first calculating the Psi-matrices
        ## (the coefficients of the reduced-form MA(infinity) representation) and 
        ## only THEN multiplying the resulting Psi-matrices at each forecast 
        ## horizon with the inv(A_0); 
        ## instead the loop DIRECTLY 
        ## calculates the THETA-MATRICES that then hold the structural 
        ## MA-infinity weights at each forecast horizon!
        ## and
    ## (2) calculates the Thetas at each forecast horizon by making use of 
        ## the previous results to avoid the recursion having to re-start
        ## the calculation of all previous intermediate-results at each
        ## forecast horizon!

    ## Also, because for the below we have to multiply matrices with 
    ## themsleves multiple times, we define the following operator 
    ## for matrix multiplication:
    # "%^%" <- function(mat,power){
    #   base = mat
    #   out = diag(nrow(mat))
    #   while(power > 1){
    #     if(power %% 2 == 1){
    #       out = out %*% base
    #     }
    #     base = base %*% base
    #     power  = power %/% 2
    #   }
    #   out %*% base
    # }
    # Marcel (08.05.2020): Note: The above function was not necessary
    # anymore after we decided to write down the calculation of the 
    # Psi-matrices differently!
    
    # because we have changed our approach as described above,
    # hence, we could also remove the function 'calculate_psi_matrices'
    # which contained the recursion and simply calculate 
    # all the Psi's of our (reduced-form) VAR-system:
    # (i.e. the MA-infinity weights):
      
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
    plot_SVAR.irfs.all.CONSTR_ALL <- plot_SVAR.irfs.all %>%
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
        
        for(b in 1:61){
          # loop to fill rows (in our case 26); to fill the 
          # coefficients at EACH STEP (hopefully in the future 
          # we'll have more steps!)
          
          coefs.df[b, 1] <-  sapply(lapply(Thetas, "[[", b), 
                                    function(x) x[row, col])[m]
          colnames(coefs.df)[1] <- paste0(impulse, "_", response, "_", "maxG")
          
        }
      
      return(as_tibble(coefs.df))
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
    plot_SVAR.irfs.maxG.CONSTR_ALL <- plot_SVAR.irfs.maxG %>%
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
      ggplot(data=plot_SVAR.irfs.all.CONSTR_ALL, aes(x=step, y=series_value)) + 
      geom_line(aes(colour=series_name), alpha=0.6, size=0.3) +
      geom_line(data=plot_SVAR.irfs.maxG.CONSTR_ALL, aes(x=step, y=series_value),
                colour="black", size=2, linetype="dashed") +
      labs(color=NULL) + 
      geom_hline(yintercept=0, color = "#514e4e", 
                 size=1) +
      scale_x_continuous(name = NULL) + 
      scale_y_continuous(name = NULL) +
      theme(axis.text=element_text(size=8),
            plot.title = element_text(size=10, face="bold", hjust = 0.5),
            axis.title=element_text(size=10),
            legend.position="none",
            # legend.text=element_text(size=5),
            #axis.text.x=element_blank(),
            plot.margin = unit(c(1,1,1,1), "mm"),
            #panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            aspect.ratio = 0.95,
            strip.background = element_rect(colour="black", 
                                            fill="white", 
                                            size=0.5, 
                                            linetype="solid")) + 
      facet_wrap(impulse ~ response,
                 scales="free_y", strip.position = "top") +
      coord_cartesian(ylim = c(-2, 2)) + 
      theme(strip.text = element_text(size = 8, margin = margin(0.7, 0.7, 0.7, 0.7, "mm")))
      

    
      impulse.responses_all.SVAR
       
      # ggsave(file="impulse_responses_all_SVAR.pdf")
      # ggsave(file="impulse_responses_all_SVAR_test.pdf")
      
       
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
                                      -c(yearmon)) %>%
                                      # -c(day, yearmon)) %>%
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
    
    # new facets labels:
    supp.labs2 <- c(`e_fin` = paste("Positive e_fin exceeding", 
                                    "2 standard deviations", sep="\n"),
                    `e_ipm` = paste("Negative e_ipm exceeding",
                                    "2 standard deviations", sep="\n"),
                    `e_macro` = paste("Positive e_m exceeding",
                                    "2 standard deviations", sep="\n"))
    
    # now we are finally in a position to plot figure 3 as well:
    time_series_epsilon_t_largeShocks <- ggplot(data = epsilon_t_valid.figure3) +
      geom_rect(data=recessions_start_end, 
                inherit.aes = FALSE,
                aes(xmin=my_start, xmax=my_end, 
                    ymin=-Inf, ymax=+Inf), 
                fill='#606060', alpha=0.3) + 
      geom_point(aes(x = my, y = series_value),
                size=2, color='#0072B2') +
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
                         breaks = seq(1960, 2015, by = 5),
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
      # coord_fixed(ratio = 5) + 
      coord_cartesian(ylim = c(2, 5)) + 
      facet_grid(series_name ~ ., scales="free_y",
                 labeller = as_labeller(supp.labs2)) + 
      theme(strip.text = element_text(size = 8, margin = margin(.5, 0.3, .5, 0.3, "cm")))
    
      time_series_epsilon_t_largeShocks
      
      #ggsave(file="time_series_epsilon_t_largeShocks.pdf")
    
    
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
    
    # we want to inspect the series of epsilon_t_maxG to get a sense
    # of when the largest values occur:
    epsilon_t_maxG %>% arrange(desc(lip))
    epsilon_t_maxG %>% arrange(lip)
    epsilon_t_maxG %>% arrange(desc(macro_h1))
    epsilon_t_maxG %>% arrange(macro_h1)
    epsilon_t_maxG %>% arrange(desc(financial_h1))
    epsilon_t_maxG %>% arrange(financial_h1)
    
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
                        dplyr::select(-c(yearmon)) %>%
                        # dplyr::select(-c(day, yearmon)) %>%
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
      geom_hline(yintercept=c(-2, 2),
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
                         breaks = seq(1960, 2015, by = 5),
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
      # coord_cartesian(ylim = c(2, 5)) + 
      facet_grid(series_name ~ ., scales="free_y") 
    
      time.series_epsilon_t_maxG
      ggsave(file="time_series_epsilon_t_maxG.pdf")
    
      ## The procedure above was initialized with a reduced-form
      ## VAR model stored in my.var;
      ## With this information at hand, a set of A_0's were isolated 
      ## we estimate SVARs for the set of A_0's:
      # my.svar <- SVAR(x = my.var, estmethod = "direct", Amat = A_0_valid[[1]],
      #                 hessian = TRUE, method="BFGS")
      
      

