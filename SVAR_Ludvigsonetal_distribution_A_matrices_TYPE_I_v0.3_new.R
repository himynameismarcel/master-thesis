#########################################################################################
### Marcel Kropp, 30.05.2020
### This script estimates SVARs following Ludvigson et al (2019)'s algorithm which
### we have described in detail in the main-text;

### In particular, the below estimates the model under various alternative 
### sets of restrictions.

### Because we have implemented the function 'LMN_algorithm' in the script 
### '20200628_functions_v0.3.R' already, for the below the only thing we needed
### to do was to add two more scenarios to the selection-operator:
### (i) FC Only (which we haven't considered so far)
### (ii) FC without the Lehman Event.

### The case (i) will be relevant for this script, the case of 
### (ii) will be relevant in another script 
### 'SVAR_Ludvigsonetal_distribution_A_matrices_TYPE_II_v0.2_new.R' where
### we replace (i) with (ii) in the comparison of how the restrictions
### affect the solution-set!


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
### Algorithm Only With ALL CONSTRAINTS
###--------------------------------------------------------------------------------------------

## Assuming that the last script we performed is SVAR_Ludvigsonetal newnew v0.10.R
## where we have performed the full set of 1.5 Mio rotations, we should save
## the set of A_0_valid that passed all constraints, before continuing with the below:

A_0_valid.ALL_CONSTR <- A_0_valid
# in a last step, we take the inverse of A_0_valid_no_constr
# to have the actual A_0_inv-matrices again:
A_0_valid.ALL_CONSTR <- lapply(A_0_valid.ALL_CONSTR, inv)


# having created A_0_validALL_CONSTR, we want to extract to sets of
# values: these are A_0_YM (which is the elemtn [2, 1]) and
# A_0_YF (Which is the element [2, 3]);
# to be able to extract these values, we have to run through the list
# and copy out the respective values and store them in a vector:

A_0_valid.ALL_CONSTR.YM <- do.call(rbind, 
                   lapply(A_0_valid.ALL_CONSTR, function(x) x[2,1]))
colnames(A_0_valid.ALL_CONSTR.YM) <- "YM"

A_0_valid.ALL_CONSTR.YF <- do.call(rbind, 
                   lapply(A_0_valid.ALL_CONSTR, function(x) x[2,3]))
colnames(A_0_valid.ALL_CONSTR.YF) <- "YF"
# lastly, we create a new matrix where we add these two column vectors together:
A_0_valid.ALL_CONSTR.comb <- as_tibble(
                    cbind(A_0_valid.ALL_CONSTR.YM, A_0_valid.ALL_CONSTR.YF))


# extracting the miminum values for YM And YF:
A_0_valid.ALL_CONSTR.comb %>%
  mutate(min1 = min(YM),
         min2 = min(YF))



###--------------------------------------------------------------------------------------------
### Algorithm Only With EVENT CONSTRAINTS
###--------------------------------------------------------------------------------------------

FE_constr <- LMN_algorithm("FE_constr", 1000000)           

                  
# Marcel, 30.05.2020
# in a last step we rename the object A_0_valid so that we can
# differentiate it from the other sets created above:
A_0_valid.FE_CONSTR <- FE_constr[[1]]
# in a last step, we take the inverse of A_0_valid_no_constr
# to have the actual A_0_inv-matrices again:
A_0_valid.FE_CONSTR <- lapply(A_0_valid.FE_CONSTR, inv)
            

# having created A_0_validFE_CONSTR, we want to extract to sets of
# values: these are A_0_YM (which is the elemtn [1, 2]) and
# A_0_YF (Which is the element [3, 2]);
# to be able to extract these values, we have to run through the list
# and copy out the respective values and store them in a vector:

A_0_valid.FE_CONSTR.YM <- do.call(rbind, 
                lapply(A_0_valid.FE_CONSTR, function(x) x[2,1]))
colnames(A_0_valid.FE_CONSTR.YM) <- "YM"

A_0_valid.FE_CONSTR.YF <- do.call(rbind, 
                lapply(A_0_valid.FE_CONSTR, function(x) x[2,3]))
colnames(A_0_valid.FE_CONSTR.YF) <- "YF"
# lastly, we create a new matrix where we add these two column vectors together:
A_0_valid.FE_CONSTR.comb <- as_tibble(
              cbind(A_0_valid.FE_CONSTR.YM, A_0_valid.FE_CONSTR.YF))
    

###--------------------------------------------------------------------------------------------
### Algorithm Only With NO CONSTRAINTS AT ALL
###--------------------------------------------------------------------------------------------


### The algorithm with no constraints at all is already implemented in the function 
### 'LMN_algorithm'. Hence, at this stage, we can simply call the function as follows:

NO_constr <- LMN_algorithm("NO_constr", 100000)   


# Marcel, 30.05.2020
# in a last step we rename the object A_0_valid so that we can
# differentiate it from the other sets created above:
A_0_valid.NO_CONSTR <- NO_constr[[1]]
# in a last step, we take the inverse of A_0_valid_no_constr
# to have the actual A_0_inv-matrices again:
A_0_valid.NO_CONSTR <- lapply(A_0_valid.NO_CONSTR, inv)
        

# having created A_0_valid_no_constr, we want to extract to sets of
# values: these are A_0_YM (which is the elemtn [2, 1]) and
# A_0_YF (Which is the element [2, 3]);
# to be able to extract these values, we have to run through the list
# and copy out the respective values and store them in a vector:

A_0_valid.NO_CONSTR.YM <- do.call(rbind, 
                    lapply(A_0_valid.NO_CONSTR, function(x) x[2,1]))
colnames(A_0_valid.NO_CONSTR.YM) <- "YM"

A_0_valid.NO_CONSTR.YF <- do.call(rbind, 
                    lapply(A_0_valid.NO_CONSTR, function(x) x[2,3]))
colnames(A_0_valid.NO_CONSTR.YF) <- "YF"
# lastly, we create a new matrix where we add these two column vectors together:
A_0_valid.NO_CONSTR.comb <- as_tibble(
                    cbind(A_0_valid.NO_CONSTR.YM, A_0_valid.NO_CONSTR.YF))
    
    
 

###--------------------------------------------------------------------------------------------
### Plotting the respective values for A_0_YF and A_0_YM
###--------------------------------------------------------------------------------------------
    
    # to be able to use facet, we have to reshape all data.frames:
    A_0_valid.FE_CONSTR.tidy <- A_0_valid.FE_CONSTR.comb   %>%
                            gather(series, data)
    # and ultimately we add another column called NO_CONSTR
    # and muliply the data-series with 100:
    A_0_valid.FE_CONSTR.tidy <- A_0_valid.FE_CONSTR.tidy %>%
                            mutate(TYPE_CONSTR = "Event Constr. Only") %>%
                            mutate(data = data*100)
    
    A_0_valid.NO_CONSTR.tidy <- A_0_valid.NO_CONSTR.comb   %>%
                            gather(series, data)   
    A_0_valid.NO_CONSTR.tidy <- A_0_valid.NO_CONSTR.tidy %>%
                            mutate(TYPE_CONSTR = "No Constr.") %>%
                            mutate(data = data*100)    
    
    A_0_valid.ALL_CONSTR.tidy <- A_0_valid.ALL_CONSTR.comb   %>%
                            gather(series, data)
    A_0_valid.ALL_CONSTR.tidy <- A_0_valid.ALL_CONSTR.tidy %>%
                            mutate(TYPE_CONSTR = "All Constr.") %>%
                            mutate(data = data*100)       
  
    

    
    ggplot(data=A_0_valid.NO_CONSTR.tidy) + 
      geom_histogram(aes(x=data, y=(..count..)/sum(..count..),
                         fill=TYPE_CONSTR),
                     colour="black",
                     alpha = 1,
                     binwidth=0.145) + 
      geom_histogram(data=A_0_valid.FE_CONSTR.tidy,
                   aes(x=data, y=(..count..)/sum(..count..),
                       fill=TYPE_CONSTR),
                   colour="black",
                   alpha = 0.8,
                   binwidth=0.019) +
      geom_histogram(data=A_0_valid.ALL_CONSTR.tidy,
                     aes(x=data, y=(..count..)/sum(..count..),
                         fill=TYPE_CONSTR),
                     colour="black",
                     alpha = 0.8,
                     binwidth=0.008) + 
      scale_x_continuous(name = NULL
                         # ,
                         # breaks = seq(-0.6, 0.6, by = 0.2),
                         # minor_breaks = NULL, labels=comma
                         ) +
      theme(legend.position="top",
            legend.title = element_blank(),
            legend.text = element_text(size = 15),
            axis.text=element_text(size=15),
            strip.text.x = element_text(size = 15)) +
      scale_y_continuous(name = NULL) +
      # very useful blog-post about how to include and color 
      # a legend:
      # https://aosmith.rbind.io/2018/07/19/manual-legends-ggplot2/
      scale_colour_discrete(labels=c("No Constraints",
                                   "Event Constraints Only",
                                   "All Constraints")) +
      theme(legend.position = "top",aspect.ratio=0.5) +
      facet_wrap(.~series)

      # ggsave(file="distribution_impact_matrices_type1.pdf")



      
###
# Analysis of Recursive Schemes
###
      
# having created A_0_validALL_CONSTR, we want to extract two sets of
# values: these are A_0_YM (which is the elemtn [2, 1]) and
# A_0_YF (Which is the element [2, 3]);
# to be able to extract these values, we have to run through the list
# and copy out the respective values and store them in a vector:

A_0_valid.ALL_CONSTR.FM <- do.call(rbind, 
                                   lapply(A_0_valid.ALL_CONSTR, function(x) x[3,1]))
colnames(A_0_valid.ALL_CONSTR.FM) <- "FM"

A_0_valid.ALL_CONSTR.MY <- do.call(rbind, 
                                   lapply(A_0_valid.ALL_CONSTR, function(x) x[1,2]))
colnames(A_0_valid.ALL_CONSTR.MY) <- "MY"


# lastly, we create a new matrix where we add these two column vectors together:
A_0_valid.ALL_CONSTR.comb.recurs <- as_tibble(cbind(A_0_valid.ALL_CONSTR.FM, 
                                                    A_0_valid.ALL_CONSTR.MY))

# extracting the minimum values
A_0_valid.ALL_CONSTR.comb.recurs %>%
  mutate(min1 = min(FM),
         min2 = min(MY))


# to be able to plot the data, we have to reshape it:
A_0_valid.ALL_CONSTR.tidy.recurs <- A_0_valid.ALL_CONSTR.comb.recurs   %>%
                            gather(series, data)
A_0_valid.ALL_CONSTR.tidy.recurs <- A_0_valid.ALL_CONSTR.tidy.recurs %>%
                            mutate(TYPE_CONSTR = "ALL CONSTRAINTS") %>%
                            mutate(data = data*100)  




# plot
ggplot(data=A_0_valid.ALL_CONSTR.tidy.recurs) + 
  geom_histogram(aes(x=data, y=(..count..)/sum(..count..),
                     fill="#288cc9"),
                 color="black", 
                 alpha = 0.8,
                 binwidth=0.1) + 
  scale_y_continuous(name = NULL) +
  scale_x_continuous(name = NULL) + 
  theme(legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=15),
        strip.text.x = element_text(size = 15)) +
  scale_y_continuous(name = NULL) +
  # very useful blog-post about how to include and color 
  # a legend:
  # https://aosmith.rbind.io/2018/07/19/manual-legends-ggplot2/
  theme(legend.position = "none",aspect.ratio=0.5) +
  facet_wrap(.~series)

#ggsave(file="distribution_impact_matrices_type3.pdf")



