############################################
##                                        ##
##    Feature importance analysis         ##
##    of image-enabled flow cytometry     ##
##    data for Figure 1B.                 ##    
##                                        ##
##    Dominik Vonficht & Viktoria Flore   ##
##                08/2024                 ##
##                                        ##
############################################ 

library(rpart)
library(rpart.plot)
library(caret)
library(data.table)
library(xgboost)
library(tidyverse)
library(rattle)
library(tidyr)
library(dplyr)
library(purrr)
library(ggpubr)
library(writexl)

set.seed(127)

## functions
fit_tree <- function(df, r, params, classes){
  ## select parameters
  df <- df %>% filter(status_fine %in% classes, element_name == r)
  df <- df[,colnames(df) %in% c('status_fine', params)]
  df <- df[,!apply(df, 2, function(x) any(is.na(x)))]
  
  ## fix feature names (exclude special characters)
  colnames(df) <- make.names(colnames(df))
  
  ## split into training and test set (70-30)
  idx <- createDataPartition(df$status_fine, p = 0.7)[[1]]
  training <- df[idx,]
  
  ## learn model, 10x CV
  control <- trainControl(method = 'repeatedcv', number = 10, repeats = 3, savePredictions = T)
  model <- train(
    status_fine ~ ., 
    data = training,
    method = 'rpart2',
    trControl = control
  )
  return(list(model = model, index = idx, df = df))
}

## data 
merged_small <- fread("./interactomics_framework/merged_small.csv")

# use events classified as singlets or non-coincident multiplets 
merged_small <- dplyr::filter(merged_small,status_fine!="conincident" & status_fine!="low_quality")
phases_fine <- c("singlet","doublet","triplet","low_quality","quadruplet_more","coincident")
phases_broad <- c("singlet"   , "multiplet",  "undefined",  "coincident")
merged_small <- merged_small %>%
  mutate_at(vars(1:17), as.numeric) %>% as.data.frame()

## models 
models <- list(
  rep1 = list(r = 'Cytostim_1.fcs_6910', subset = colnames(merged_small)[-(18:22)], classes = phases_fine),
  rep2 = list(r = 'Cytostim_5.fcs_16768', subset = colnames(merged_small)[-(18:22)], classes = phases_fine),
  rep3 = list(r = 'Cytostim_2.fcs_34414', subset = colnames(merged_small)[-(18:22)], classes = phases_fine),
  rep4 = list(r = 'Cytostim_3.fcs_46733', subset = colnames(merged_small)[-(18:22)], classes = phases_fine)) %>% 
  map(~ {
    fit_tree(merged_small, .x$r, .x$subset, .x$classes) 
  })

## parameter subset
feature_imp <- models %>% map_df(~ {
  varImp(.x$model)$importance %>%
    as_tibble(rownames = 'feature') %>% 
    arrange(desc(Overall)) %>% 
    mutate(feature = factor(feature, levels = rev(feature)))
}, .id = 'rep')

feature_imp$feature <- as.character(feature_imp$feature)

feature_imp <-  feature_imp %>% dplyr::group_by(feature) %>% 
  dplyr::mutate(median=median(Overall))

## plot feature importance (Figure 1B)
levels <- (feature_imp %>% dplyr::group_by(feature) %>% dplyr::mutate(median=median(Overall)) %>% arrange(desc(median)))$feature %>% unique
p <- ggstripchart(feature_imp %>% mutate(feature = factor(feature, levels = levels)), 
                  x = "feature", y = "Overall", fill = "feature", jitter = 0.2, 
                  size = 6, add.params = list(color = "black"), shape = 21) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  scale_fill_manual(values = pals::kelly()) +
  scale_color_manual(values = pals::kelly()) +
  stat_summary(fun = median, fun.min = median, fun.max = median, 
               geom = "errorbar", color = "black", width = 0.6, size = 0.7)
print(p)

# export source data 
write_xlsx(p$data, path = "/fast/AG_Haas/Viktoria/NatMethods_revision/feature_importance_analysis/Figure1B_sourcedata.xlsx", col_names = T, format_headers = T)


