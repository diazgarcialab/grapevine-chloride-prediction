install.packages("ggplot2")
install.packages("waves")
install.packages("reshape2")
install.packages("mdatools")
install.packages("scico")
install.packages("ggpubr")
install.packages("doSNOW")
install.packages("pls")
install.packages("doParallel")
install.packages("kernlab")
install.packages("foreach")
install.packages("lme4")
library(doParallel)
library(foreach)
library(ggplot2)
library(waves)
library(reshape2)
library(mdatools)
library(scico)
library(ggpubr)
library(doSNOW)
library(pls)
library(kernlab)
library(lme4)



# reading file from step 01
all_hyper<-read.csv('02_proceseddata/interpolated_signals_ALLFILEs.csv')
rownames(all_hyper)<-all_hyper$X #set row names based on values in the "X" column, making new column
all_hyper<-all_hyper[,-1] #remove the first column (X) from the data frame

# applying all 13 types of spectral pretreatments implemented in the waves package

all_hyper_treated<-pretreat_spectra(df = all_hyper,pretreatment = 1:13)

# reading chloride readings from chloridometer
chr<-rbind(data.frame(read.csv('00_rawdata/trial01/chloridometer_readings.csv'),trial=1), 
           
           data.frame(read.csv('00_rawdata/trial02/chloridometer_readings.csv'),trial=2))

table(substr(chr$svc_id,1,9)) #to see how many samples in each unique svc_id

chr<-chr[-which(chr$svc_id=='no scan'),] #remove no scan


rownames(chr)<-chr$svc_id #set row names of chr to values in svc_id column
str(chr)
chr$rep <- as.factor(chr$rep)
chr$trial <- as.factor(chr$trial)
chr$trt <- as.factor(chr$trt)

# Mixed model
model <- lmer(average ~ trt + (1|genotype)+(1|trial)+(1|genotype:trial), data = chr)
summary(model)

# find data in common from chr and the pre treated data
(f<-intersect(rownames(chr),rownames(all_hyper_treated[[1]])))

# create allcors list to store the correlations between refelence and chloride content
allcors <- list()
for (i in 1:length(all_hyper_treated)) {
  x <- data.frame(chl = chr[f, 'average'], all_hyper_treated[[i]][f, ])
  cors <- numeric()
  for (j in 2:ncol(x)) {
    cors[j - 1] <- cor(x[, 1], x[, j])
  }
  names(cors) <- colnames(all_hyper_treated[[i]])
  allcors[[i]] <- cors
}

alldat <- matrix(NA, ncol(all_hyper), length(allcors))
rownames(alldat) <- colnames(all_hyper)
for (i in 1:length(allcors)) {
  alldat[names(allcors[[i]]), i] <- allcors[[i]]
}

colnames(alldat) <- c('raw data', 'SNV', 'SNV + first derivative', 'SNV + second derivative', 'first derivative', 'second derivative', 'SG', 'SNV + SG', 'gap segment derivative, ws=11', 'SG + 1st derivative, ws=5', 'SG + first derivative, ws=11', 'SG + second derivative, ws=5', 'SG + second derivative, ws=11')
alldat <- data.frame(wavelength = seq(339, 2515, 1), alldat)

# Create a data frame for plotting
alldat2 <- melt(alldat, id.vars = c('wavelength'))
colnames(alldat2)[2:3] <- c('treatment', 'correlation')

# Plot correlation between reflectance and chloride
ggplot(alldat2, aes(x = wavelength, y = abs(correlation))) +
  geom_smooth(se = FALSE) +  
  facet_wrap(~treatment, ncol = 7) +  
  labs(x = "Wavelength (nm)", y = "Pearson's Correlation") +
  theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank()) +
  coord_flip() +
  scale_y_continuous(breaks = c(0.00, 0.25))

# By utilizing the innospectra refelectance data and chloride readings, correlations 
# between reflectance and chlroide content can be calculated using the same procedure
# as described in above codes.

# Set up parallel processing
cores <- 4  # Set the number of cores
cl <- parallel::makeCluster(cores)
registerDoParallel(cl)

# Load necessary libraries on each parallel worker
parallel::clusterEvalQ(cl, {
  library(mdatools)
})

# PLSR models using mdatools
plsr_cor<-numeric(length(all_hyper_treated))
plsr_cor <- foreach (ii = 1:length(all_hyper_treated), .combine = c) %dopar% {
  hyp <- all_hyper_treated[[ii]]
  f<-intersect(rownames(hyp),rownames(chr))
  x<-data.frame(chr=chr[f,'average'],hyp[f,])
  
  # plsr model
  # partitioning training/validation sets
  ix <- sample(1:nrow(x), size = 0.3*nrow(x))
  x1<-x[-ix,]
  x2<-x[ix,]
  
  m<-pls(as.matrix(x1[,-1]),as.matrix(x1$chr),x.test = as.matrix(x2[,-1]),y.test = x2$chr,ncomp = 12)
  
  pred<-m$res$test$y.pred[,1,1]
  obse<-x2[names(pred),1]  
  plsr_cor[ii]<-cor(pred,obse,use = 'complete.obs')
  
}
parallel::stopCluster(cl)
names(plsr_cor)<-colnames(alldat[, 2:14])


# PLSDA using categories
chr1 <- chr[chr$trt != 0, ] #remove trt 0

# Set up parallel processing
cores <- 6  
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(mdatools)
})

# Create empty lists to store results
plsr_misscl <- list()

for (ii in 1:length(all_hyper_treated)){
  hyp<-all_hyper_treated[[ii]]
  (f<-intersect(rownames(hyp),rownames(chr1)))
  x<-data.frame(chr1=chr1[f,'average'],hyp[f,])
  
  x$chr1<-ifelse(x$chr<2538.79 ,'excluder','non-excluder')
  ## plsr model
  ## partitioning training/validation sets
  
  ix <- sample(1:nrow(x), size = 0.3*nrow(x))
  x1<-x[-ix,]
  x2<-x[ix,]
  
  m<-mdatools::plsda(as.matrix(x1[,-1]),as.matrix(x1$chr1),x.test = as.matrix(x2[,-1]),c.test = x2$chr1)
  plsr_misscl[[ii]]<-getConfusionMatrix(m$calres,ncomp = 12)
  
  
}
stopCluster(cl)
names(plsr_misscl)<-names(all_hyper_treated)

# Create empty vectors to store accuracy and specificity
accuracy_vector <- specificity_vector <- numeric(length(plsr_misscl))

# Calculate accuracy and specificity for each model
for (i in seq_along(plsr_misscl)) {
  conf_matrix <- plsr_misscl[[i]]
  
  # Calculate accuracy
  accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  accuracy_vector[i] <- accuracy
  
  # Calculate specificity
  specificity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
  specificity_vector[i] <- specificity
}

# Combine results into a data frame
result_df <- data.frame(
  treatment = names(all_hyper_treated),
  accuracy = accuracy_vector,
  specificity = specificity_vector
)

# Print or use the results as needed
print(result_df)

# PLSR, RF, and SVM models implemented in the waves package
# Set up parallel processing
cores <- 6  
cl <- makeCluster(cores)
registerDoParallel(cl)

# Load necessary libraries
parallel::clusterEvalQ(cl, {
  library(ggplot2)
  library(pls)
  library(randomForest)
  library(e1071)  # For SVM
  library(dplyr)
})

# Define models to train
mdls <- c('pls', 'rf', 'svmLinear')

# Model runner function that returns both residuals and performance
model_runner <- function(hyps, reference_data, treatment, model, tune, num.iterations) {
  hyp <- hyps[[treatment]]
  f <- intersect(rownames(hyp), rownames(reference_data))
  x <- data.frame(unique.id = f, reference = reference_data[f, 'average'], hyp[f, ])
  
  # 70-30 train-test split
  set.seed(123)
  train_idx <- sample(1:nrow(x), size = 0.7 * nrow(x))
  train_data <- x[train_idx, ]
  test_data <- x[-train_idx, ]
  
  # Prepare predictors and target variables
  train_predictors <- as.matrix(train_data[, !(names(train_data) %in% c('unique.id', 'reference'))])
  test_predictors <- as.matrix(test_data[, !(names(test_data) %in% c('unique.id', 'reference'))])
  train_target <- train_data$reference
  test_target <- test_data$reference
  
  # Initialize variables to accumulate metrics
  residuals_list <- list()
  rmse_values <- numeric(num.iterations)
  
  # Perform multiple iterations
  for (iter in 1:num.iterations) {
    if (model == 'pls') {
      mdl <- plsr(train_target ~ train_predictors, ncomp = tune)
      predicted <- as.vector(predict(mdl, newdata = test_predictors, ncomp = tune))
    } else if (model == 'rf') {
      mdl <- randomForest(train_predictors, train_target, ntree = tune)
      predicted <- predict(mdl, test_predictors)
    } else if (model == 'svmLinear') {
      mdl <- svm(train_predictors, train_target, kernel = "linear")
      predicted <- predict(mdl, test_predictors)
    }
    
    # Calculate residuals and RMSE
    residuals <- predicted - test_target
    rmse_values[iter] <- sqrt(mean(residuals^2))
    
    # Store residuals for each iteration
    residuals_list[[iter]] <- data.frame(
      unique.id = test_data$unique.id,
      observed = test_target,
      predicted = predicted,
      residual = residuals,
      iteration = iter
    )
  }
  
  # Average RMSE across iterations
  avg_rmse <- mean(rmse_values)
  
  # Return both residuals and performance as a list
  return(list(
    performance = data.frame(RMSE = avg_rmse, model = model, treatment = treatment),
    residuals = bind_rows(residuals_list)
  ))
}

# Define models, treatments, and parameters
mods <- data.frame(
  treatment = rep(names(all_hyper_treated), each = 3),
  model = mdls,
  tune = c(12, 500, 1)  # Example tuning parameters for PLS, RF, and SVM
)

# Run models in parallel with 100 iterations
superOut <- foreach(i = 1:nrow(mods), .packages = c('pls', 'randomForest', 'e1071', 'dplyr')) %dopar% {
  model_runner(
    hyps = all_hyper_treated,
    reference_data = chr,
    treatment = mods$treatment[i],
    model = mods$model[i],
    tune = mods$tune[i],
    num.iterations = 100  # Set iterations to 100
  )
}

# Stop the parallel cluster
stopCluster(cl)

# Extract performances and residuals from the output
model_performances <- lapply(superOut, function(x) x$performance)
residuals_combined <- lapply(superOut, function(x) x$residuals)

# Combine performances and residuals into data frames
model_pefs <- bind_rows(model_performances)
residuals_df <- bind_rows(residuals_combined)

# Save the combined results
saveRDS(superOut, 'combined_70_30_ncomp12_metricsWAVES.rds')

# Print the outputs
print(model_pefs)
print(residuals_df)

# We can repeat the same process for all five scenarios and for NIR-S-G1 spectral data

# Residuals plot
custom_colors <- c(
  "pls" = "violetred3",
  "rf" = "mediumblue",
  "svmLinear" = "tan4"
)
# Define shapes for each model
custom_shapes <- c(
  "pls" = 16,      # Circle
  "rf" = 15,       # Square
  "svmLinear" = 17 # Triangle
)
# Custom function to split long Pretreatment labels into multiple lines
split_label <- function(label) {
  sapply(label, function(x) {
    paste(strwrap(x, width = 15), collapse = "\n")  # Adjust 'width' as needed
  })
}

# Apply the custom labeller to Pretreatment
pretreatment_labeller <- as_labeller(split_label)

# Updated plot code
residual_plot <- ggplot(residuals_df, aes(x = predicted, y = residuals, color = model, shape = model)) +
  geom_point(alpha = 0.7, size = 1.5) +  # Adjust alpha and size for better visibility
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Adding facet labels for scenarios and pretreatments with custom wrapping
  facet_grid(Pretreatment ~ Scenario, labeller = labeller(
    Scenario = label_both, 
    Pretreatment = pretreatment_labeller
  )) +
  
  # Adjusting axis labels
  labs(
    x = expression('Predicted Chloride Content (' * mu * mol %.% g^-1 * ' of DW)'), 
    y = expression('Standardized Residuals (' * mu * mol %.% g^-1 * ' of DW)')
  ) +
  
  # Setting minimal theme and custom colors/shapes
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  
  # Adjusting theme elements for larger fonts and preventing label cutoff
  theme(
    strip.text.x = element_text(size = 16, face = "bold"),  
    strip.text.y = element_text(size = 13, lineheight = 1.2),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16, face = "bold"),  
    axis.text.x = element_text(size = 12),                  
    axis.text.y = element_text(size = 12),                  
    axis.ticks.x = element_line(color = "black", size = 0.8),  
    axis.ticks.y = element_line(color = "black", size = 0.8),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 14),                # Legend text font size
    legend.title = element_text(size = 16, face = "bold"), # Legend title font size
    legend.position = "bottom"
  ) +
  
  # Adjusting the coordinate limits
  coord_cartesian(xlim = c(-1000, 5000), ylim = c(-5000, 5000))

# Print the plot
print(residual_plot)













