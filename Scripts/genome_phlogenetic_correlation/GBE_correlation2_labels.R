#Calculations for correlation of traits from Genome evolution in plant pathogenic bacteria metadata.
  #statistics that take into account the strength of phylogenetic signal in the residual variation
#Notes for statistics mentioned in manuscript text:
  #"Amongst the representative genomes in Figure 2, phylogenetically informed statistics revealed evolutionary linkage between larger genome sizes and the presence of T3SS or T6SS with phylogenetic signals (alpha) of 1.25 and 2.05 respectively."
    #These alpha values are produced by Analysis C, logistic regression by phyloglm, and found in the file phyloglm_pairwise_results.csv
#"Amongst the representative genomes chosen for this review (Figure 2), a phylogenetically informed logistic regression revealed that if a pathovar has T6SS, it is 6.5 times more likely to be seedborne (BH adjusted p-value 0.013), with almost no phylogenetic signal (alpha= 89.4)."
  #These results are from Analysis B, lines 182-196, with BH adjusted p-value and alpha value from Analysis C and found in the file phyloglm_pairwise_results.csv
#"Furthermore, logistic regression analysis of the reference genomes (Figure 2) revealed that bacteria that cause disease on a diverse range of plant hosts are 83.5% less likely to be seedborne (BH adjusted p-value 0.095)."
  #These results are found in Analysis B, lines 226-240, and adjusted p-value is from Analysis C and found in file phyloglm_pairwise_results.csv
#"Amongst the reference genomes (Figure 2), T2SS CAZyme repertoire and nutritional lifestyle have a significant relationship (BH adjusted p-value = 0.021) with necrotrophs having a predicted 57.7 more CAZymes than biotrophs."
  #These results are from Analysis A, lines 118-163.

#All calculations are listed below
#Note that this code was written with the assistance of ChatGPT version GPT-5.5

#Load needed libraries
library(ape)
library(phylolm)
library(caper)


#For both binary or multivariate variables: Phylogenetic Logistic Regression, using maximum likelihood
#For continuous variables:Linear Regression, using Phylogenetic Generalized Least Squares

#Load data from excel file into a dataframe
traits <- read_excel("Species_Variables.xlsx")
traits_df <- as.data.frame(traits)

# Binary variables
#Diverse trait: where 1 is affecting multiple plant families and -1 most commonly affecting a single plant family. 
traits_df$diverse      <- as.numeric(traits_df$diverse)
#Seedborne trait: where 1 is known in literature to be seedborne and -1 is not known to be seedborne
traits_df$Seedborne    <- as.numeric(traits_df$Seedborne)
#Insectborne trait: where 1 is known in literature to be insectborne (including nematodes) and -1 is not known to be insectborne 
traits_df$Insectborne  <- as.numeric(traits_df$Insectborne)
#Presence of Type 3 Secretion System: 1 is present and -1 is absent 
traits_df$T3SS         <- as.numeric(traits_df$T3SS)
#Presence of Type 6 Secretion System: 1 is present and -1 is absent
traits_df$T6SS         <- as.numeric(traits_df$T6SS)

# Continuous variables
#genome size of reference in base pairs
traits_df$genome <- as.numeric(traits_df$genome)
#The number of CAZymes identified by run_dbcan3 analysis
traits_df$T2SS   <- as.numeric(traits_df$T2SS)

# Multistate categorical variables
#Cotyledon type of host with -1 being monocots, 1 being dicots and 0 if both are affected
traits_df$Host_cotyledon <- as.numeric(traits_df$Host_cotyledon)
#Lifestyle trait: where -1 is nonvascular, 1 is vascular and 0 is both
traits_df$life_style <- as.numeric(traits_df$life_style)
#Nutritional mode trait: where -1 is necrotrophic, 1 is biotrophic and 0 is hemibiotrophic
traits_df$nutritional_mode <- as.numeric(traits_df$nutritional_mode)

#Load rooted tree with branch lengths
tree <- read.tree("SpeciesTree_rooted_node_labels.txt")
tree$root.edge <- NULL
#add constant to branch lengths to eliminate any branch lengths of 0 
tree$edge.length <- tree$edge.length + 1e-6
#add name labels of species from dataframe
rownames(traits_df) <- traits_df$Species

#Analysis A
#Run linear regression on genome size and all other traits
# 1. Setup variables
y_var <- "genome"
x_vars <- c("diverse", "Host_cotyledon", "life_style", "nutritional_mode", 
            "Seedborne", "Insectborne", "T2SS", "T3SS", "T6SS")

# 2. Initialize results
results <- data.frame(Variable = x_vars, Estimate = NA, P_Value = NA)


comp_data <- comparative.data(
  phy = tree,
  data = traits_df,
  names.col = "Species",
  vcv = TRUE
)

# 3. The Loop
for (i in seq_along(x_vars)) {
  form <- as.formula(paste(y_var, "~", x_vars[i]))
  
  # phylolm is more robust to the "ABNORMAL_TERMINATION" error
  # model = "lambda" estimates the phylogenetic signal via ML
  fit <- try(phylolm(form, data = comp_data$data, phy = comp_data$phy, model = "lambda"), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    s <- summary(fit)$coefficients
    results$Estimate[i] <- s[2, "Estimate"]
    results$P_Value[i]  <- s[2, "p.value"]
    results$Lambda[i] <- fit$optpar
    
  } else {
    message(paste("Failed on variable:", x_vars[i]))
  }
}

# 4. Error Correction
results$FDR_P <- p.adjust(results$P_Value, method = "fdr")
results$BH_P <- p.adjust(results$P_Value, method = "BH")
summary(fit)
print(results)
#results for genome size 
#Variable     Estimate    P_Value     FDR_P
#1          diverse  -70323.0947 0.05167215 0.2527198
#2   Host_cotyledon  114263.7465 0.05615995 0.2527198
#3       life_style  -11868.1380 0.64032826 0.8362569
#4 nutritional_mode   71575.6327 0.22803176 0.6840953
#5        Seedborne   44100.0428 0.40614634 0.7348954
#6      Insectborne   -9491.6344 0.86543033 0.9736091
#7             T2SS     250.8154 0.65042202 0.8362569
#8             T3SS  174817.2745 0.40827521 0.7348954
#9             T6SS   -2803.5147 0.97632691 0.9763269


#run linear regression on number of CAZymes and all other traits
# 1. Setup variables
y_var <- "T2SS"
x_vars <- c("diverse", "Host_cotyledon", "life_style", "nutritional_mode", 
            "Seedborne", "Insectborne", "genome", "T3SS", "T6SS")

# 2. Initialize results
results2 <- data.frame(Variable = x_vars, Estimate = NA, P_Value = NA)

# 3. The Loop
for (i in seq_along(x_vars)) {
  form <- as.formula(paste(y_var, "~", x_vars[i]))
  
  # phylolm is more robust to the "ABNORMAL_TERMINATION" error
  # model = "lambda" estimates the phylogenetic signal via ML
  fit <- try(phylolm(form, data = comp_data$data, phy = comp_data$phy, model = "lambda"), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    s <- summary(fit)$coefficients
    results2$Estimate[i] <- s[2, "Estimate"]
    results2$P_Value[i]  <- s[2, "p.value"]
    results2$Lambda[i] <- fit$optpar
  } else {
    message(paste("Failed on variable:", x_vars[i]))
  }
}

# 4. Error Correction
results2$FDR_P <- p.adjust(results2$P_Value, method = "fdr")
results2$BH_P <- p.adjust(results2$P_Value, method = "BH")

print(results2)
#T2SS
#Variable      Estimate      P_Value        FDR_P         BH_P       
#1          diverse   2.031943034 7.535369e-01 7.535369e-01 7.535369e-01
#2   Host_cotyledon  13.250732679 5.659635e-02 1.273418e-01 1.273418e-01 
#3       life_style   1.677590177 7.502521e-01 7.535369e-01 7.535369e-01 
#4 nutritional_mode -28.874383366 4.814996e-03 2.166748e-02 2.166748e-02 
#5        Seedborne  18.970529845 5.454715e-02 1.273418e-01 1.273418e-01 
#6      Insectborne   9.282547497 3.668380e-01 4.716488e-01 4.716488e-01 
#7           genome   0.000031246 5.386400e-11 4.847760e-10 4.847760e-10 
#8             T3SS -25.417705746 1.871034e-01 2.806551e-01 2.806551e-01 
#9             T6SS -21.234734371 1.146360e-01 2.063449e-01 2.063449e-01 
#Since the estimate for nutritional mode is -28.87, there is a decrease of an average of 28.87 CAZymes per
#genome from biotroph to hemibiotroph and from hemibiotroph to necrotroph, meaning necrtrophs have an average 
#of 57.7 more CAZymes than biotrophs. 

#Analysis B
#For binary and multivariate variables use phyloglm to run Phylogenetic Logistic Regression
#Logistic regression of T3SS and T6SS presence 
T3_T6 <- phyloglm(
  T3SS ~ T6SS,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(T3_T6)
#Parameter estimate(s) from MPLE:
#alpha: 2.861726 
#Coefficients:
# (Intercept)        T6SS 
#-0.9236709   0.2588370 
#p.value = 0.4599

#Logistic regression of Seedborne and T6SS presence
Seed_T6 <- phyloglm(
  Seedborne ~ T6SS,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_T6)
#alpha: 89.44582 
#Coefficients:
#  (Intercept)        T6SS 
#-0.8587212   1.8830015
#p.value = 0.0003773 ***
#log-odds = 1.883, 1e^1.883 = 6.573
#bacteria with a T6SS have an ~6.5x higher chance of being seedborne

#Logistic regression of seedborne and presence of T3SS
Seed_T3 <- phyloglm(
  Seedborne ~ T3SS,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_T3)
#alpha: 4.055587 
#Coefficients:
#  (Intercept)        T3SS 
#-0.1527242  -1.3121003 
#p.value = 0.1023

#logistic regression of seedborne and lifestyle (vascular, nonvascular, both)
Seed_life <- phyloglm(
  Seedborne ~ life_style,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_life)
#alpha: 91.51766 
#Coefficients:
#  Estimate   StdErr z.value p.value
#(Intercept) 0.029371 0.274102  0.1072  0.9147
#life_style  0.077193 0.166079  0.4648  0.6421

#logistic regression for seedborne and host family diversity
Seed_div <- phyloglm(
  Seedborne ~ diverse,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_div)
#alpha: 88.68352 
#Coefficients:
#  Estimate   StdErr z.value   p.value    
#(Intercept)  1.37220  0.37298  3.6790 0.0002341 ***
#  diverse     -1.80150  0.57147 -3.1524 0.0016194 ** 
#log odds -1.8015, e^-1.8015 = 0.165, (0.165-1)x100 = 83.5%
#diverse pathogens are 83.5% less likely to be seedborne

#Logistic regression for seedborne and host cotyledon type
Seed_host <- phyloglm(
  Seedborne ~ Host_cotyledon,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_host)
#alpha: 89.25097 
#Coefficients:
# Estimate  StdErr z.value p.value
#(Intercept)     0.25376 0.29098  0.8721  0.3832
#Host_cotyledon  0.42918 0.28352  1.5137  0.1301

#Logistic regression of seedborne and nutritional mode (biotroph, hemibiotroph, necrotrophs)
Seed_nutr <- phyloglm(
  Seedborne ~ nutritional_mode,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_nutr)
#alpha: 90.36137 
#Coefficients:
# Estimate    StdErr z.value p.value
#(Intercept)      -0.023118  0.277130 -0.0834  0.9335
#nutritional_mode -0.363489  0.335981 -1.0819  0.2793

#logistical regression of T3SS presence and insectborne
Seed_insect <- phyloglm(
  Seedborne ~ Insectborne,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Seed_insect)
#alpha: 91.38999 
#Coefficients:
# Estimate   StdErr z.value p.value  
#(Intercept) -0.19779  0.32031 -0.6175 0.53692  
#Insectborne  0.74646  0.40966  1.8221 0.06844 

#logistical regression of T3SS presence and host family diversity
T3_div <- phyloglm(
  T3SS ~ diverse,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(T3_div)
#alpha: 2.491512 
#Coefficients:
# Estimate   StdErr z.value p.value
#(Intercept) -1.05587  0.89756 -1.1764  0.2394
#diverse     -0.18594  0.31137 -0.5972  0.5504

#logistical regression of T3SS presence and host cotyledon
T3_Host <- phyloglm(
  T3SS ~ Host_cotyledon,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(T3_Host)
#alpha: 1.596597 
#Coefficients:
# Estimate   StdErr z.value p.value
#(Intercept)    -0.69676  1.03466 -0.6734  0.5007
#Host_cotyledon -0.15483  0.17411 -0.8892  0.3739

#Logistical regression of T3SS presence and lifestyle (vascular, nonvascular, and both)
T3_life <- phyloglm(
  T3SS ~ life_style,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(T3_life)
#alpha: 1.731857 
#Coefficients:
#  Estimate   StdErr z.value p.value
#(Intercept) -0.84961  0.96801 -0.8777  0.3801
#life_style   0.13907  0.12276  1.1328  0.2573

#Logistical regression of host diversity and host cotyledon type
Div_Host <- phyloglm(
  diverse ~ Host_cotyledon,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Div_Host)
#alpha: 91.1435 
#Coefficients:
#  Estimate   StdErr z.value   p.value    
#(Intercept)    -1.22044  0.35516 -3.4363 0.0005898 ***
#Host_cotyledon  0.32080  0.34206  0.9379 0.3483194 

#Logistical regression of host diversity and lifestyle (vascular, nonvascular, both)
Div_life <- phyloglm(
  diverse ~ life_style,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Div_life)
#alpha: 60.06954 
#Coefficients:
#Estimate   StdErr z.value  p.value   
#(Intercept) -0.92173  0.32235 -2.8594 0.004245 **
#life_style  -0.25428  0.20470 -1.2422 0.214160 

#Logistical regression of host diversity and nutritional mode (necrotroph, biotroph, hemibiotroph)
Div_nutri<- phyloglm(
  diverse ~ nutritional_mode,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Div_nutri)
#alpha: 89.8575 
#Coefficients:
#  Estimate   StdErr z.value   p.value    
#(Intercept)      -1.20125  0.34829 -3.4490 0.0005626 ***
#  nutritional_mode -1.04039  0.44827 -2.3209 0.0202920 *  

#Logistical regression of host diversity and insectborne
Div_insect<- phyloglm(
  diverse ~ Insectborne,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Div_insect)
#alpha: 63.0742 
#Coefficients:
# Estimate   StdErr z.value  p.value   
#(Intercept) -1.22540  0.39649 -3.0906 0.001997 **
#  Insectborne  0.65104  0.43153  1.5087 0.131379

#Logistical regression of host diversity and T3SS presence
Div_T3<- phyloglm(
  diverse ~ T3SS,
  data = comp_data$data,
  phy = comp_data$phy,
  method = "logistic_MPLE", btol = 12
)
summary(Div_T3)
#alpha: 3.704709 
#Coefficients:
#  Estimate   StdErr z.value p.value  
#(Intercept) -0.98251  0.57601 -1.7057 0.08806 .
#T3SS        -1.98511  1.31543 -1.5091 0.13127

#Analysis C
#Pairwise logistic regression to produce adjusted p-value for multiple tests
#load needed libraries
library(phylolm)
library(dplyr)
library(purrr)
library(tibble)

# List variables to test
vars <- c(
  "Seedborne",
  "nutritional_mode",
  "Insectborne",
  "T3SS",
  "diverse",
  "Host_cotyledon",
  "life_style",
  "T6SS", "T2SS", "genome"
)

# Create all pairwise combinations
combos <- expand.grid(
  response = vars,
  predictor = vars,
  stringsAsFactors = FALSE
)

# Remove self-comparisons
combos <- combos %>%
  filter(response != predictor)

# Function to run phyloglm safely
run_phyloglm <- function(response, predictor, data, phy) {
  
  formula <- as.formula(
    paste(response, "~", predictor)
  )
  
  model <- tryCatch(
    phyloglm(
      formula,
      data = data,
      phy = phy,
      method = "logistic_MPLE",
      btol = 20
    ),
    error = function(e) return(NULL)
  )
  
  # Return NA if model failed
  if (is.null(model)) {
    return(
      tibble(
        response = response,
        predictor = predictor,
        estimate = NA,
        stderr = NA,
        zvalue = NA,
        pvalue = NA,
        alpha = NA
      )
    )
  }
  
  summ <- summary(model)
  
  # Extract coefficient table
  coef_tab <- summ$coefficients
  
  # Second row = predictor
  tibble(
    response = response,
    predictor = predictor,
    estimate = coef_tab[2, "Estimate"],
    stderr = coef_tab[2, "StdErr"],
    zvalue = coef_tab[2, "z.value"],
    pvalue = coef_tab[2, "p.value"],
    alpha = model$alpha
  )
}

# Run all models
results <- pmap_dfr(
  combos,
  ~ run_phyloglm(..1, ..2,
                 data = comp_data$data,
                 phy = comp_data$phy)
)

# Adjust p-values
results <- results %>%
  mutate(
    p_adjust_BH = p.adjust(pvalue, method = "BH"),
    p_adjust_bonf = p.adjust(pvalue, method = "bonferroni")
  )

# View and download results
results
write.csv(results, "phyloglm_pairwise_results.csv", row.names = FALSE)

