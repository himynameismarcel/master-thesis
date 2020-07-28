#########################################################################################
### Marcel Kropp, 31.05.2020
### This script estimates SVARs following Ludvigson et al (2019)'s algorithm which
### we have described in detail in the main-text;

### In particular, the below estimates the model under various alternative 
### sets of restrictions.

### Because we have implemented the function 'LMN_algorithm' in the script 
### '20200628_functions_v0.3.R' already, for the below the only thing we needed
### to do was to add two more scenarios to the selection-operator:
### (i) FC Only (which we haven't considered so far)
### (ii) FC without the Lehman Event.

### The case (ii) will be relevant for this script, the case of 
### (i) will be relevant in another script 
### 'SVAR_Ludvigsonetal_distribution_A_matrices_TYPE_I_v0.3_new.R' where
### we replace (ii) with (i) in the comparison of how the restrictions
### affect the solution-set!


###--------------------------------------------------------------------------------------------
### Algorithm Only With ALL CONSTRAINTS EXCEPT Lehman Event
###--------------------------------------------------------------------------------------------

ALL_NO_LEHM_constr <- LMN_algorithm("ALL_NO_LEHM_constr", 1000000)           

                  
# Marcel, 30.05.2020
# in a last step we rename the object A_0_valid so that we can
# differentiate it from the other sets created above:
A_0_valid.ALL_NO_LEHM_CONSTR <- ALL_NO_LEHM_constr[[1]]
# in a last step, we take the inverse of A_0_valid_no_constr
# to have the actual A_0_inv-matrices again:
A_0_valid.ALL_NO_LEHM_CONSTR <- lapply(A_0_valid.ALL_NO_LEHM_CONSTR, inv)
            


# having created A_0_validFE_CONSTR, we want to extract to sets of
# values: these are A_0_YM (which is the elemtn [2, 1]) and
# A_0_YF (Which is the element [2, 3]);
# to be able to extract these values, we have to run through the list
# and copy out the respective values and store them in a vector:

A_0_valid.ALL_NO_LEHM_CONSTR.YM <- do.call(rbind, 
                lapply(A_0_valid.ALL_NO_LEHM_CONSTR, function(x) x[2,1]))
colnames(A_0_valid.ALL_NO_LEHM_CONSTR.YM) <- "YM"

A_0_valid.ALL_NO_LEHM_CONSTR.YF <- do.call(rbind, 
                                  lapply(A_0_valid.ALL_NO_LEHM_CONSTR, function(x) x[2,3]))
colnames(A_0_valid.ALL_NO_LEHM_CONSTR.YF) <- "YF"
# lastly, we create a new matrix where we add these two column vectors together:
A_0_valid.ALL_NO_LEHM_CONSTR.comb <- as_tibble(cbind(A_0_valid.ALL_NO_LEHM_CONSTR.YM, 
                                             A_0_valid.ALL_NO_LEHM_CONSTR.YF))
    
    
 

###--------------------------------------------------------------------------------------------
### Plotting the respective values for A_0_YF and A_0_YM
###--------------------------------------------------------------------------------------------
    
    # to be able to use facet, we have to reshape all data.frames:
    A_0_valid.ALL_NO_LEHM_CONSTR.tidy <- A_0_valid.ALL_NO_LEHM_CONSTR.comb   %>%
                            gather(series, data)
    # and ultimately we add another column called NO_CONSTR
    # and muliply the data-series with 100:
    A_0_valid.ALL_NO_LEHM_CONSTR.tidy <- A_0_valid.ALL_NO_LEHM_CONSTR.tidy %>%
                            mutate(TYPE_CONSTR = "All Constr. w.o. LEHMAN") %>%
                            mutate(data = data*100)

    
    ggplot(data=A_0_valid.NO_CONSTR.tidy) + 
      geom_histogram(aes(x=data, y=(..count..)/sum(..count..),
                         fill=TYPE_CONSTR),
                     colour="black",
                     alpha = 1,
                     binwidth=0.145) + 
      geom_histogram(data=A_0_valid.ALL_NO_LEHM_CONSTR.tidy,
                   aes(x=data, y=(..count..)/sum(..count..),
                       fill=TYPE_CONSTR),
                   colour="black",
                   alpha = 0.8,
                   binwidth=0.012) +
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
      scale_colour_discrete(guide = "legend",
                          name="",
                          labels=c("No Constraints", 
                                   "All Constr. (w.o. Lehman Event)",
                                   "All Constraints")) +
      theme(legend.position = "top",aspect.ratio=0.5) +
      facet_wrap(.~series)

      # ggsave(file="distribution_impact_matrices_type2.pdf")


