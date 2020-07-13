#########################################################################################
### Marcel Kropp, 18.05.2020

## While the script 'SVAR_Ludvigsonetal newnew v0.10.R' runs the LMN-algorithm for the 
## complete set of restrictions, the below calculate three variance to be able
## to produce a figure similar to the one in LMN where the authors compare 
## (i) the case of No Constraints,
## (ii) the case of Correlation Constraints only,
## (iii) and the case of Correlation Constraints + Non-Negative + Real Activity Constraints!
## Note that the last constraints (RA) was taken from LMN (2018) and replaced one 
## of the Non-Negative Constraints in LMN (2019);

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

## Note that after reading in the above scripts, this script can run 
## stand-lone (there are not any other dependencies!)

###--------------------------------------------------------------------------------------------
### Algorithm Only With EVENT CONSTRAINTS (only the 'non-negative constraints' + 'RA constraint')
### in combination with the CORRELATION CONSTRAINTS
###--------------------------------------------------------------------------------------------

FC_NN_RA_constr <- LMN_algorithm("FC_NN_RA_constr", 10000)
    
#######################
### Calculation of Impulse Response-Functions for the
### set only containing the 
### Correlation Constraints + non-negative constr. + RA constr.
#######################
    
    l <- 3 # number of equations
    p <- 6 # number of lags

    B_block <- t(sapply(FC_NN_RA_constr[[2]]$varresult, `[[`, 1))
    
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'matsplitter' was moved to the R-Script 
    # 20200628_functions.R    
    #---------------------------------------------------------------------------  

    B_array <- matsplitter(B_block, 3, 3)
    B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
    constant <- as.matrix(B_block[, ncol(B_block)])
    
    ## STEP 2: 
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
    
    for(z in 0:(length(FC_NN_RA_constr[[1]])-1)){
      
        Thetas[[z+1]] <- 
              lapply(psi_list, function(x){
              x %*% inv(FC_NN_RA_constr[[1]][[z+1]])
              })
              names(Thetas)[length(Thetas)] <- paste0(c("A_0_"), z)

    }
    
    
    #-------------------------------------------------------------------
    # extract coefficients to create impulse-responses:
    # now all that is left is to efficiently extract the respective
    # coefficients;
          
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'extract.coefs.irf' was moved to the R-Script 
    # 20200628_functions.R  
    #--------------------------------------------------------------------------- 

    Fin_Macro.coefs <- 100*extract.coefs.irf("Fin", "Macro", 1, 3)
    Fin_Ipm.coefs <- 100*extract.coefs.irf("Fin", "Ipm", 2, 3)
    Fin_Fin.coefs <- 100*extract.coefs.irf("Fin", "Fin", 3, 3)
    Macro_Macro.coefs <- 100*extract.coefs.irf("Macro", "Macro", 1, 1)
    Macro_Ipm.coefs <- 100*extract.coefs.irf("Macro", "Ipm", 2, 1)
    Macro_Fin.coefs <- 100*extract.coefs.irf("Macro", "Fin", 3, 1)
    Ipm_Macro.coefs <- 100*extract.coefs.irf("Ipm", "Macro", 1, 2)
    Ipm_Ipm.coefs <- 100*extract.coefs.irf("Ipm", "Ipm", 2, 2)
    Ipm_Fin.coefs <- 100*extract.coefs.irf("Ipm", "Fin", 3, 2)
    
    # 14.06.2020
    # Damit wir um die Kollektion der IRFs die Ränder beim
    # Plotten etwas dicker zeichnen können, müssen wir hier
    # an dieser Stelle dedicated data.frames anlegen, in denen
    # wir genau diese max und min-values abspeichern;
    # die zugehörige Funktion findet sich weiter unten in diesem Skript
    # (in letzer Instant wäre es am besten, sämtliche Funktionen
    # in einem dedicated R-Skrpt abzuspeichern, das man einlesen kann
    # )
    
    # Wichtig: der unten stehende Befehl funktioniert nur, wenn
    # die Funktion extract.min.max() vorher eingelesen wurde:
    output.min.max <- extract.min.max(
      Fin_Macro.coefs, Fin_Ipm.coefs, Fin_Fin.coefs, 
      Macro_Macro.coefs, Macro_Ipm.coefs, Macro_Fin.coefs,
      Ipm_Macro.coefs, Ipm_Ipm.coefs, Ipm_Fin.coefs)
    
    
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
    
    
    # 14.06.2020:
    # Da wir bei der Erstellung des data.frames für den Fall 
    # ohne constraints die Funktion 'generate.irfs.plot' gebaut
    # haben, um nicht mehrmals das gleiche aufzurufen zu müssen,
    # können wir auch hier an dieser Stelle auf diese Funktion 
    # zugreifen;
    # Wichtig: diese Funktion sollte am besten schlussendlich
    # zusammmen mit allen anderen von mir gebauten Funktionen
    # in einem separaten R-Skript geparkt und dann eingelesen
    # werden!
    
    plot_SVAR.irfs.all <- generate.irfs.plots(
      Macro_Macro.coefs, Macro_Ipm.coefs, Macro_Fin.coefs, 
      Ipm_Macro.coefs, Ipm_Ipm.coefs, Ipm_Fin.coefs,
      Fin_Macro.coefs, Fin_Ipm.coefs, Fin_Fin.coefs 
    )
    
    # 14.06.2020:
    # the below command does the same as the above command
    # but only for the data.frames with the respective
    # min/max-series:
    
    plot_SVAR.irfs.min.max <- generate.irfs.plots(                        
      output.min.max[[4]],
      output.min.max[[5]],
      output.min.max[[6]],
      output.min.max[[7]],
      output.min.max[[8]],
      output.min.max[[9]],
      output.min.max[[1]], 
      output.min.max[[2]],
      output.min.max[[3]]
    )
    
    
    # we need to create ordered factors, otherwise 'facet_wrap' combines the plots
    # in alphabetical order:
    plot_SVAR.irfs.all.FE_FC <- plot_SVAR.irfs.all %>%
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
    
    # next, we do the exact same for the min.max-series:
    # 14.06.2020
    plot_SVAR.irfs.min.max.FE_FC <- plot_SVAR.irfs.min.max %>%
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

FC_constr <- LMN_algorithm("FC_constr", 3000)
    

#######################
### Calculation of Impulse Response-Functions
### for the set consisting of Correlation Constraints
### only
#######################

## STEP 1: extracting reduced-form matrices from
## VAR-estimation from above
#----------------------------------------------

l <- 3 # number of equations
p <- 6 # number of lags

B_block <- t(sapply(FC_constr[[2]]$varresult, `[[`, 1))

#---------------------------------------------------------------------------
# 28.06.2020 (Marcel)
# Note, the function 'matsplitter' was moved to the R-Script 
# 20200628_functions.R    
#---------------------------------------------------------------------------      


B_array <- matsplitter(B_block, 3, 3)
B_list <- lapply(seq(dim(B_array)[3]), function(x) B_array[ , , x])
constant <- as.matrix(B_block[, ncol(B_block)])

## STEP 2: 
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

for(z in 0:(length(FC_constr[[1]])-1)){
  
  Thetas[[z+1]] <- 
    lapply(psi_list, function(x){
      x %*% inv(FC_constr[[1]][[z+1]])
    })
  names(Thetas)[length(Thetas)] <- paste0(c("A_0_"), z)
  
}


#-------------------------------------------------------------------
# extract coefficients to create impulse-responses:
# now all that is left is to efficiently extract the respective
# coefficients;


#---------------------------------------------------------------------------
# 28.06.2020 (Marcel)
# Note, the function 'extract.coefs.irf' was moved to the R-Script 
# 20200628_functions.R  
#--------------------------------------------------------------------------- 


Fin_Macro.coefs <- 100*extract.coefs.irf("Fin", "Macro", 1, 3)
Fin_Ipm.coefs <- 100*extract.coefs.irf("Fin", "Ipm", 2, 3)
Fin_Fin.coefs <- 100*extract.coefs.irf("Fin", "Fin", 3, 3)
Macro_Macro.coefs <- 100*extract.coefs.irf("Macro", "Macro", 1, 1)
Macro_Ipm.coefs <- 100*extract.coefs.irf("Macro", "Ipm", 2, 1)
Macro_Fin.coefs <- 100*extract.coefs.irf("Macro", "Fin", 3, 1)
Ipm_Macro.coefs <- 100*extract.coefs.irf("Ipm", "Macro", 1, 2)
Ipm_Ipm.coefs <- 100*extract.coefs.irf("Ipm", "Ipm", 2, 2)
Ipm_Fin.coefs <- 100*extract.coefs.irf("Ipm", "Fin", 3, 2)

# 14.06.2020
# Damit wir um die Kollektion der IRFs die Ränder beim 
# Ploten etwas dicker zeichnen können, müssen wir hier
# an dieser Stelle dedicated data.frames anlegen, in denen
# wir genau diese max und min-values abspeichern;
# die zugehörige Funktion findet sich weiter unten in diesem Skript
# (in letzter Instanz wäre es am besten, sämtliche Funktionen
# in einem dedicated R-Skript abzuspeichern, das man einlesen kann
# )

# Wichtig: der unten stehende Befehl funktioniert nur, wenn 
# die Funktion extract.min.max() vorher eingelesen wurde:
output.min.max <- extract.min.max(
            Fin_Macro.coefs, Fin_Ipm.coefs, Fin_Fin.coefs, 
            Macro_Macro.coefs, Macro_Ipm.coefs, Macro_Fin.coefs,
            Ipm_Macro.coefs, Ipm_Ipm.coefs, Ipm_Fin.coefs)

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


# 14.06.2020:
# Da wir bei der Erstellung des data.frames für den Fall 
# ohne constraints die Funktion 'generate.irfs.plot' gebaut
# haben, um nicht mehrmals das gleiche aufzurufen zu müssen,
# können wir auch hier an dieser Stelle auf diese Funktion 
# zugreifen;
# Wichtig: diese Funktion sollte am besten schlussendlich
# zusammmen mit allen anderen von mir gebauten Funktionen
# in einem separaten R-Skript geparkt und dann eingelesen
# werden!

plot_SVAR.irfs.all <- generate.irfs.plots(
            Macro_Macro.coefs, Macro_Ipm.coefs, Macro_Fin.coefs, 
            Ipm_Macro.coefs, Ipm_Ipm.coefs, Ipm_Fin.coefs,
            Fin_Macro.coefs, Fin_Ipm.coefs, Fin_Fin.coefs 
)

# 14.06.2020:
# the below command does the same as the above command
# but only for the data.frames with the respective
# min/max-series:

plot_SVAR.irfs.min.max <- generate.irfs.plots(                        
              output.min.max[[4]],
              output.min.max[[5]],
              output.min.max[[6]],
              output.min.max[[7]],
              output.min.max[[8]],
              output.min.max[[9]],
              output.min.max[[1]], 
              output.min.max[[2]],
              output.min.max[[3]]
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

# next, we do the exact same for the min.max-series:
# 14.06.2020
plot_SVAR.irfs.min.max.FC_only <- plot_SVAR.irfs.min.max %>%
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
### Algorithm Only With NO CONSTRAINTS AT ALL (UNCONSTRAINED SET)
###--------------------------------------------------------------------------------------------

NO_constr <- LMN_algorithm("NO_constr", 500)
    
    #######################
    ### Calculation of Impulse Response-Functions for the
    ### set no constraints at all
    #######################
    
    l <- 3 # number of equations
    p <- 6 # number of lags
    
    B_block <- t(sapply(NO_constr[[2]]$varresult, `[[`, 1))
    
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'matsplitter' was moved to the R-Script 
    # 20200628_functions.R    
    #---------------------------------------------------------------------------   
    
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
    
    for(z in 0:(length(NO_constr[[1]])-1)){
      
      Thetas[[z+1]] <- 
        lapply(psi_list, function(x){
          x %*% inv(NO_constr[[1]][[z+1]])
        })
      names(Thetas)[length(Thetas)] <- paste0(c("A_0_"), z)
      
    }
    
    
    #-------------------------------------------------------------------
    # extract coefficients to create impulse-responses:
    # now all that is left is to efficiently extract the respective
    # coefficients;
    
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'extract.coefs.irf' was moved to the R-Script 
    # 20200628_functions.R  
    #--------------------------------------------------------------------------- 
    
    Fin_Macro.coefs <- 100*extract.coefs.irf("Fin", "Macro", 1, 3)
    Fin_Ipm.coefs <- 100*extract.coefs.irf("Fin", "Ipm", 2, 3)
    Fin_Fin.coefs <- 100*extract.coefs.irf("Fin", "Fin", 3, 3)
    Macro_Macro.coefs <- 100*extract.coefs.irf("Macro", "Macro", 1, 1)
    Macro_Ipm.coefs <- 100*extract.coefs.irf("Macro", "Ipm", 2, 1)
    Macro_Fin.coefs <- 100*extract.coefs.irf("Macro", "Fin", 3, 1)
    Ipm_Macro.coefs <- 100*extract.coefs.irf("Ipm", "Macro", 1, 2)
    Ipm_Ipm.coefs <- 100*extract.coefs.irf("Ipm", "Ipm", 2, 2)
    Ipm_Fin.coefs <- 100*extract.coefs.irf("Ipm", "Fin", 3, 2)
    
    # 14.06.2020
    # Damit wir um die Kollektion der IRFs die Ränder beim 
    # Ploten etwas dicker zeichnen können, müssen wir hier
    # an dieser Stelle dedicated data.frames anlegen, in denen
    # wir genau diese max und min-values abspeichern;
    # als Grundlage dienen die Daten aus obigen data.frames
    # (also Fin_Macro.coefs, Fin_Imp.coefs, usw).
    # Fin_Macro.coefs.max <- Fin_Macro.coefs
    # test <- apply(Fin_Macro.coefs, 1, FUN=min)
    # test2 <- apply(Fin_Macro.coefs, 1, FUN=max) 
    
    # after the above test, we need to find an efficient way to 
    # loop through the list of the above data.frames and extract
    # min and max-values; the best way to do so would be to 
    # use a function which we then also can use for the other
    # scenarios where we need to extract max- and min-values:
    
    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'extract.min.max' was moved to the R-Script 
    # 20200628_functions.R  
    #--------------------------------------------------------------------------- 
    
    # we test the function:
    output.min.max <- extract.min.max(Fin_Macro.coefs, Fin_Ipm.coefs, Fin_Fin.coefs, 
                    Macro_Macro.coefs, Macro_Ipm.coefs, Macro_Fin.coefs,
                    Ipm_Macro.coefs, Ipm_Ipm.coefs, Ipm_Fin.coefs)

    
    # having extracted the min and max across all series, 
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
    # with the full set of ifrs
    # into a tidy format for plotting and at the same time
    # binds them together

    #---------------------------------------------------------------------------
    # 28.06.2020 (Marcel)
    # Note, the function 'generate.irfs.plots' was moved to the R-Script 
    # 20200628_functions.R  
    #--------------------------------------------------------------------------- 
    
    plot_SVAR.irfs.all <- generate.irfs.plots(
                        Macro_Macro.coefs, Macro_Ipm.coefs, Macro_Fin.coefs, 
                        Ipm_Macro.coefs, Ipm_Ipm.coefs, Ipm_Fin.coefs,
                        Fin_Macro.coefs, Fin_Ipm.coefs, Fin_Fin.coefs 
                        )
    
    # 14.06.2020:
    # the below command does the same as the above command
    # but only for the data.frames with the respective
    # min/max-series:
    plot_SVAR.irfs.min.max <- generate.irfs.plots(                        
                            output.min.max[[4]],
                            output.min.max[[5]],
                            output.min.max[[6]],
                            output.min.max[[7]],
                            output.min.max[[8]],
                            output.min.max[[9]],
                            output.min.max[[1]], 
                            output.min.max[[2]],
                            output.min.max[[3]]
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
    
    # next, we do the exact same for the min.max-series:
    # 14.06.2020:
    plot_SVAR.irfs.min.max.NO_CONSTR <- plot_SVAR.irfs.min.max %>%
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
plot_SVAR.irfs.all.FE_FC <- plot_SVAR.irfs.all.FE_FC %>%
                          mutate(type = "FE_FC")
plot_SVAR.irfs.all.FC_only <- plot_SVAR.irfs.all.FC_only %>%
                          mutate(type = "FC_only")                
plot_SVAR.irfs.all.NO_CONSTR <- plot_SVAR.irfs.all.NO_CONSTR %>%
                          mutate(type = "NO_CONSTR") 

# 14.06.2020:
# here we add the columns to the newly created data.frames that
# only hold the respective min/max-series:
plot_SVAR.irfs.min.max.NO_CONSTR <- plot_SVAR.irfs.min.max.NO_CONSTR %>%
                          mutate(type = "NO_CONSTR")
plot_SVAR.irfs.min.max.FC_only <- plot_SVAR.irfs.min.max.FC_only %>%
                          mutate(type = "FC_only")
plot_SVAR.irfs.min.max.FE_FC <- plot_SVAR.irfs.min.max.FE_FC %>%
                          mutate(type = "FE_FC")

###----------------------------------------------------------------------------------------
### finally we plot the impulse responses:  
###----------------------------------------------------------------------------------------    
    impulse.responses_all.SVAR.COMPLETE <- 
      ggplot(data=plot_SVAR.irfs.all.NO_CONSTR, aes(x=step, y=series_value)) + 
      geom_line(aes(group=series_name, colour=type), alpha=0.01, size=0.9) +
      # here we add the borders to the unconstrained set
      geom_line(data=plot_SVAR.irfs.min.max.NO_CONSTR, aes(x=step, y=series_value,
                            group=series_name, colour=type),
                alpha = 1, size = 1, linetype = "dashed") + 
      
      geom_line(data=plot_SVAR.irfs.all.FC_only, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 0.09, size = 0.3) +
      # here we add the borders to the set with CORR Constraints only:
      geom_line(data=plot_SVAR.irfs.min.max.FC_only, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 1, size = 1, linetype = "dashed") + 
  
      geom_line(data=plot_SVAR.irfs.all.FE_FC, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 0.1, size = 0.1) +
      # here we add the borders to the set with CORR Constraints 
      # and the non-negative event constraints:
      geom_line(data=plot_SVAR.irfs.min.max.FE_FC, aes(x=step, y=series_value,
                        group=series_name, colour=type),
                        alpha = 1, size = 1, linetype = "dashed") + 
  
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
            aspect.ratio = 0.95,
            strip.text = element_text(size = 8, 
                                      margin = margin(0.7, 0.7, 0.7, 0.7, "mm")),
            strip.background = element_rect(colour="black", 
                                            fill="white", 
                                            size=0.5, 
                                            linetype="solid")) + 
      coord_cartesian(ylim = c(-3, 3)) + 
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      #https://stackoverflow.com/questions/35712062/scale-fill-discrete-does-not-change-label-names
      scale_colour_discrete(
                      labels = c("FC Only", 
                                 "FC + NN FE + RA FE", 
                                 "No Constraints"))

    
      impulse.responses_all.SVAR.COMPLETE
       
      ## ggsave(file="impulse_responses_all_SVAR_unconstr_constr.pdf")
      
       

