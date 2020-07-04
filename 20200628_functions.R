## 28.06.2020
## This file contains all the functions that I have created and use
## in the course of the analyses for the Master's Thesis:

## Hence, before running any other script, this script should be executed
## before-hand.


## -------------------------------------
## Impulse Response Functions
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



## -------------------------------------
## Forward Error Variance Decomposition
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
