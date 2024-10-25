if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("mixOmics", force = TRUE)
install.packages("doParallel")
install.packages("plsmod")
install.packages("parsnip")
install.packages("tidymodels")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("rsample")
install.packages("recipes")
install.packages("caret")
install.packages("MASS")
install.packages("magrittr")
library(mixOmics)
library(plsmod)
library(tidyverse)
library(tidymodels)
library(dplyr)
library(ggplot2)
library(tidyr)
library(parsnip)
library(ggplot2)
library(rsample)
library(recipes)
library(caret)
library(magrittr)

#FIRST TRIAL
#read in spectral data
df_spectra1 <- read.csv("~/chloride_prediction/00_rawdata/trial01/230606_svc_reflectance.csv", check.names = F)

#read in chloride data
df_chloride1 <- read.csv("~/chloride_prediction/00_rawdata/trial01/chloridometer_readings.csv")
df_chloride1<-df_chloride1[-which(df_chloride1$svc_id=='no scan'),]

#rename df_spectra and data cleaning
df_spectra1 <- df_spectra1 %>%
  dplyr::rename(svc_id = scan) 
  
#merge spectra with chloride data
df_data1 <- merge(df_chloride1, df_spectra1, by = "svc_id")

#data split into train and test
set.seed(123)
perm_split <- initial_split(df_data1, prop = .7)
perm_train <- training(perm_split)
perm_test <- testing(perm_split)

#5-fold cross validation
set.seed(456)
perm_folds <-  vfold_cv(perm_train, v=5)

#Preparing pls model ---------------------------

#recipe for pls
pls_rec <- recipe(average ~ ., data = perm_train) %>%
  update_role(svc_id, pot.number, genotype, rep, trt,  reading1, reading2, reading3, new_role = 'id') %>%
  update_role(average, new_role = 'outcome') %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors())

#pls tuning specifications
pls_tuning <- parsnip::pls(num_comp = tune()) %>% 
  set_mode("regression") %>% 
  set_engine("mixOmics")

#create tuning workflow
pls_workflow <- workflow() %>% 
  add_recipe(pls_rec) %>% 
  add_model(pls_tuning)

#enable parallel processing to speed up
doParallel::registerDoParallel()

#number of pc to test
comp_grid <- tibble(num_comp = seq(from = 1, to = 20, by = 1))

#pls model tuning
tuned_pls_results <- pls_workflow %>% 
  tune_grid(resamples = perm_folds,
            grid = comp_grid,
            metrics = metric_set(mae, rmse, rsq))

#tuning summary
model_results <- tuned_pls_results %>% 
  collect_metrics()

#plot tuning summary
tuned_pls_results %>% 
  collect_metrics() %>% 
  ggplot(aes(num_comp, mean, col = .metric)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(n.breaks = 20) +
  labs(x = "Number of components",
       y = "Indicator",
       title = "Plot of MAE, RMSE and R-SQ vs number of components for TRAINING dataset, with 10-fold repeated cross validation") +
  facet_grid(.metric ~., scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none")

#select best model
tuned_best <- tuned_pls_results %>% 
  #select_best("rsq")    #absolute best based on either "mae", "rmse", or "rsq"
  select_by_pct_loss(    #selecting least complex model with no more than 5% loss of rmse
    metric = "rmse",
    limit = 5, 1
  )

#update model and workflow with tuned parameters
updated_pls_model <-  parsnip::pls(num_comp = tuned_best$num_comp) %>% 
  set_mode("regression") %>% 
  set_engine("mixOmics")

updated_workflow <- pls_workflow %>% 
  update_model(updated_pls_model)

#create pls model
pls_model <- updated_workflow %>% 
  fit(data = perm_train)

#assessing pls model -------------------------------------------

#extract vip
pls_vip <- pls_model %>% 
  extract_fit_parsnip() %>% 
  tidy()

#plot vip
pls_vip %>% 
  filter(term != "Y", # outcome variable col name
         component == tuned_best$num_comp) %>% 
  group_by(component) %>% 
  ungroup() %>% 
  ggplot(aes(x = as.numeric(term), y = value, fct_reorder(term, value), fill = factor(component))) +
  geom_line(show.legend = F) +
  labs(x = "Wavelength (nm)",
       y = "VIP") +
  theme_bw()


# We can repeat the same process for all five scenarios and identify optimal ncomp values



