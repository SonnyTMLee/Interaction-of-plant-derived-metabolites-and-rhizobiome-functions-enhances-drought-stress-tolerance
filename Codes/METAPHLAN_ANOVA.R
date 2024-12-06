##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
########Location ANOVA###################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################


library(dplyr)
library(tidyr)
library(ggpubr)

rm(list=ls())

# Set working directory
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_METAPHLAN")

# Read the CSV file (adjust the path as needed)
data <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Population <- as.factor(data$Population)
data$Loc_eco <- as.factor(data$Loc_eco)

# Select metabolite columns (assuming they are all after the first three columns)
metabolite_columns <- names(data)[6:ncol(data)]

# Create a results data frame to store ANOVA results
anova_results <- data.frame(Metabolite = character(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Create a results data frame for significant Metaphlan with Wilcoxon results
final_results <- data.frame(Metabolite = character(),
                            Comparison = character(),
                            p_adjusted = numeric(),
                            stringsAsFactors = FALSE)

# Perform ANOVA and Wilcoxon test for each metabolite
for (metabolite in metabolite_columns) {
  # Remove rows with NA values for the current metabolite
  cleaned_data <- data %>% filter(!is.na(data[[metabolite]]))
  
  # Check if there are enough data points for ANOVA
  if (nrow(cleaned_data) >= 3) {  # At least 3 points needed for ANOVA
    # Perform ANOVA
    model <- aov(cleaned_data[[metabolite]] ~ Location, data = cleaned_data)
    summary_model <- summary(model)
    p_value <- summary_model[[1]][["Pr(>F)"]][1]  # Extract p-value for the Location factor
    
    # Store the ANOVA results
    anova_results <- rbind(anova_results, data.frame(Metabolite = metabolite, p_value = p_value))}}


write.csv(anova_results, file = "Metaphlan_ANOVA_Location.csv", row.names = FALSE)


##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########LOCATION Wilcox##################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
# Select path columns for comparison
path_columns <- names(data)[6:ncol(data)]  # Select Path columns

# Create a results data frame for Wilcoxon results
wilcox_results <- data.frame(Metabolite = character(),
                             Comparison = character(),
                             W = numeric(),
                             p_value = numeric(),
                             Higher = character(),
                             stringsAsFactors = FALSE)

# Perform Wilcoxon tests for each path
for (path in path_columns) {
  # Extract values for each location
  il_values <- data[[path]][data$Location == "IL"]
  h_values <- data[[path]][data$Location == "H"]
  
  # Check if both groups have enough data for the test
  if (length(il_values) > 0 && length(h_values) > 0) {
    # Perform Wilcoxon test
    wilcox_test <- wilcox.test(il_values, h_values)
    
    # Determine which location has the higher median
    higher_location <- ifelse(median(il_values) > median(h_values), "IL", "H")
    
    # Store results
    wilcox_results <- rbind(wilcox_results, data.frame(
      Metabolite = path,
      Comparison = "IL vs H",
      W = wilcox_test$statistic,
      p_value = wilcox_test$p.value,
      Higher = higher_location
    ))
  }
}

# Save Wilcoxon results to CSV
write.csv(wilcox_results, file = "Metaphlan_WILCOX_Location.csv", row.names = FALSE)

# Display results


wilcox_results <- read.csv("Metaphlan_WILCOX_Location.csv", as.is = TRUE)
anova_results_Location <- read.csv("Metaphlan_ANOVA_Location.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = wilcox_results, anova_results_Location, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge2, "FINAL_LOCATION_Merge_Location.csv")


##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########Ecotype ANOVA###################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

rm(list=ls())

data <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Population <- as.factor(data$Population)
data$Loc_eco <- as.factor(data$Loc_eco)

# Select metabolite columns (assuming they are all after the first three columns)
metabolite_columns <- names(data)[6:ncol(data)]

# Create a results data frame to store ANOVA results
anova_results <- data.frame(Metabolite = character(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Create a results data frame for significant Metaphlan with Wilcoxon results
final_results <- data.frame(Metabolite = character(),
                            Comparison = character(),
                            p_adjusted = numeric(),
                            stringsAsFactors = FALSE)

# Perform ANOVA and Wilcoxon test for each metabolite
for (metabolite in metabolite_columns) {
  # Remove rows with NA values for the current metabolite
  cleaned_data <- data %>% filter(!is.na(data[[metabolite]]))
  
  # Check if there are enough data points for ANOVA
  if (nrow(cleaned_data) >= 3) {  # At least 3 points needed for ANOVA
    # Perform ANOVA
    model <- aov(cleaned_data[[metabolite]] ~ Ecotype, data = cleaned_data)
    summary_model <- summary(model)
    p_value <- summary_model[[1]][["Pr(>F)"]][1]  # Extract p-value for the Ecotype factor
    
    # Store the ANOVA results
    anova_results <- rbind(anova_results, data.frame(Metabolite = metabolite, p_value = p_value))}}


write.csv(anova_results, file = "Metaphlan_ANOVA_Ecotype.csv", row.names = FALSE)



##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########Ecotype Wilcox##################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################


# Select metabolite columns (assuming they are all after the first three columns)


# Select path columns for comparison
path_columns <- names(data)[6:ncol(data)]  # Select Path columns

# Create a results data frame for Wilcoxon results
wilcox_results <- data.frame(Metabolite = character(),
                             Comparison = character(),
                             W = numeric(),
                             p_value = numeric(),
                             Higher = character(),
                             stringsAsFactors = FALSE)

# Perform Wilcoxon tests for each path
for (path in path_columns) {
  # Extract values for each Ecotype
  il_values <- data[[path]][data$Ecotype == "W"]
  h_values <- data[[path]][data$Ecotype == "D"]
  
  # Check if both groups have enough data for the test
  if (length(il_values) > 0 && length(h_values) > 0) {
    # Perform Wilcoxon test
    wilcox_test <- wilcox.test(il_values, h_values)
    
    # Determine which Ecotype has the higher median
    higher_Ecotype <- ifelse(median(il_values) > median(h_values), "W", "D")
    
    # Store results
    wilcox_results <- rbind(wilcox_results, data.frame(
      Metabolite = path,
      Comparison = "W vs D",
      W = wilcox_test$statistic,
      p_value = wilcox_test$p.value,
      Higher = higher_Ecotype
    ))
  }
}

# Save Wilcoxon results to CSV
write.csv(wilcox_results, file = "Metaphlan_WILCOX_Ecotype.csv", row.names = FALSE)

# Display results


wilcox_results <- read.csv("Metaphlan_WILCOX_Ecotype.csv", as.is = TRUE)
anova_results_Ecotype <- read.csv("Metaphlan_ANOVA_Ecotype.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = wilcox_results, anova_results_Ecotype, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge2, "FInal_Ecotype.csv")





##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########POPULATION ANOVA#################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

rm(list=ls())

# Read the CSV file (adjust the path as needed)
data <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Population <- as.factor(data$Population)
data$Loc_eco <- as.factor(data$Loc_eco)

# Select metabolite columns (assuming they are all after the first three columns)
metabolite_columns <- names(data)[6:ncol(data)]

# Create a results data frame to store ANOVA results
anova_results <- data.frame(Metabolite = character(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Create a results data frame for significant Metaphlan with Wilcoxon results
final_results <- data.frame(Metabolite = character(),
                            Comparison = character(),
                            p_adjusted = numeric(),
                            stringsAsFactors = FALSE)

# Perform ANOVA and Wilcoxon test for each metabolite
for (metabolite in metabolite_columns) {
  # Remove rows with NA values for the current metabolite
  cleaned_data <- data %>% filter(!is.na(data[[metabolite]]))
  
  # Check if there are enough data points for ANOVA
  if (nrow(cleaned_data) >= 3) {  # At least 3 points needed for ANOVA
    # Perform ANOVA
    model <- aov(cleaned_data[[metabolite]] ~ Population, data = cleaned_data)
    summary_model <- summary(model)
    p_value <- summary_model[[1]][["Pr(>F)"]][1]  # Extract p-value for the Ecotype factor
    
    # Store the ANOVA results
    anova_results <- rbind(anova_results, data.frame(Metabolite = metabolite, p_value = p_value))}}


write.csv(anova_results, file = "Metaphlan_ANOVA_Population.csv", row.names = FALSE)





##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########Population ANOVA#################
####PAIRWISE Wilcox
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

# Sample data
# df <- read.csv("your_data_file.csv")  # Replace with your actual data file

# Step 1: Perform ANOVA for each metabolite
df <- read.csv("summarized_data.csv")  # Replace with your actual file path

metabolite_names <- colnames(df)[grepl("Path", colnames(df))]
anova_results <- list()

for (metabolite in metabolite_names) {
  anova_model <- aov(as.formula(paste(metabolite, "~ Population")), data = df)
  anova_summary <- summary(anova_model)
  anova_results[[metabolite]] <- anova_summary[[1]]$`Pr(>F)`[1]
}


# Step 2: Pairwise Wilcoxon test for significant ANOVA results
wilcox_results <- list()
results_to_write <- data.frame(Metabolite = character(), Group1 = character(), Group2 = character(), P_Value = numeric(), Higher_Group = character(), stringsAsFactors = FALSE)

for (metabolite in metabolite_names) {
  p_value <- anova_results[[metabolite]]
  if (p_value < 0.05) {
    wilcox_test <- pairwise.wilcox.test(df[[metabolite]], df$Population, p.adjust.method = "fdr")
    wilcox_results[[metabolite]] <- wilcox_test
    
    comparisons <- wilcox_test$p.value
    significant_comparisons <- which(comparisons < 0.05, arr.ind = TRUE)
    
    for (comparison in seq_len(nrow(significant_comparisons))) {
      group1 <- rownames(comparisons)[significant_comparisons[comparison, 1]]
      group2 <- colnames(comparisons)[significant_comparisons[comparison, 2]]
      
      # Calculate medians
      median_group1 <- median(df[df$Population == group1, metabolite], na.rm = TRUE)
      median_group2 <- median(df[df$Population == group2, metabolite], na.rm = TRUE)
      
      # Determine which group is higher
      higher_group <- ifelse(median_group1 > median_group2, group1, group2)
      
      results_to_write <- rbind(results_to_write, data.frame(Metabolite = metabolite, Group1 = group1, Group2 = group2, P_Value = comparisons[significant_comparisons[comparison, 1], significant_comparisons[comparison, 2]], Higher_Group = higher_group, stringsAsFactors = FALSE))
    }
  }
}

# Write results to CSV
write.csv(results_to_write, "Metaphlan_WILCOX_Population.csv", row.names = FALSE)

# Print the written results
print(results_to_write)


anova_results_population <- read.csv("Metaphlan_ANOVA_Population.csv", as.is = TRUE)



# Display results


wilcox_results <- read.csv("Metaphlan_WILCOX_Population.csv", as.is = TRUE)
anova_results_Ecotype <- read.csv("Metaphlan_ANOVA_Population.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = wilcox_results, anova_results_Ecotype, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge2, "Final_Population.csv")





##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########Loc_eco ANOVA#################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

rm(list=ls())

data <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Population <- as.factor(data$Population)
data$Loc_eco <- as.factor(data$Loc_eco)

# Select metabolite columns (assuming they are all after the first three columns)
metabolite_columns <- names(data)[6:ncol(data)]

# Create a results data frame to store ANOVA results
anova_results <- data.frame(Metabolite = character(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Create a results data frame for significant Metaphlan with Wilcoxon results
final_results <- data.frame(Metabolite = character(),
                            Comparison = character(),
                            p_adjusted = numeric(),
                            stringsAsFactors = FALSE)

# Perform ANOVA and Wilcoxon test for each metabolite
for (metabolite in metabolite_columns) {
  # Remove rows with NA values for the current metabolite
  cleaned_data <- data %>% filter(!is.na(data[[metabolite]]))
  
  # Check if there are enough data points for ANOVA
  if (nrow(cleaned_data) >= 3) {  # At least 3 points needed for ANOVA
    # Perform ANOVA
    model <- aov(cleaned_data[[metabolite]] ~ Loc_eco, data = cleaned_data)
    summary_model <- summary(model)
    p_value <- summary_model[[1]][["Pr(>F)"]][1]  # Extract p-value for the Ecotype factor
    
    # Store the ANOVA results
    anova_results <- rbind(anova_results, data.frame(Metabolite = metabolite, p_value = p_value))}}


write.csv(anova_results, file = "Metaphlan_ANOVA_Loc_eco.csv", row.names = FALSE)





##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########Loc_eco ANOVA#################
####PAIRWISE Wilcox
##########################################
##########################################
##########################################
##########################################
##########################################


# Sample data
# df <- read.csv("your_data_file.csv")  # Replace with your actual data file

# Step 1: Perform ANOVA for each metabolite
df <- read.csv("summarized_data.csv")  # Replace with your actual file path

metabolite_names <- colnames(df)[grepl("Path", colnames(df))]
anova_results <- list()

for (metabolite in metabolite_names) {
  anova_model <- aov(as.formula(paste(metabolite, "~ Loc_eco")), data = df)
  anova_summary <- summary(anova_model)
  anova_results[[metabolite]] <- anova_summary[[1]]$`Pr(>F)`[1]
}


# Step 2: Pairwise Wilcoxon test for significant ANOVA results
wilcox_results <- list()
results_to_write <- data.frame(Metabolite = character(), Group1 = character(), Group2 = character(), P_Value = numeric(), Higher_Group = character(), stringsAsFactors = FALSE)

for (metabolite in metabolite_names) {
  p_value <- anova_results[[metabolite]]
  if (p_value < 0.05) {
    wilcox_test <- pairwise.wilcox.test(df[[metabolite]], df$Loc_eco, p.adjust.method = "fdr")
    wilcox_results[[metabolite]] <- wilcox_test
    
    comparisons <- wilcox_test$p.value
    significant_comparisons <- which(comparisons < 0.05, arr.ind = TRUE)
    
    for (comparison in seq_len(nrow(significant_comparisons))) {
      group1 <- rownames(comparisons)[significant_comparisons[comparison, 1]]
      group2 <- colnames(comparisons)[significant_comparisons[comparison, 2]]
      
      # Calculate medians
      median_group1 <- median(df[df$Loc_eco == group1, metabolite], na.rm = TRUE)
      median_group2 <- median(df[df$Loc_eco == group2, metabolite], na.rm = TRUE)
      
      # Determine which group is higher
      higher_group <- ifelse(median_group1 > median_group2, group1, group2)
      
      results_to_write <- rbind(results_to_write, data.frame(Metabolite = metabolite, Group1 = group1, Group2 = group2, P_Value = comparisons[significant_comparisons[comparison, 1], significant_comparisons[comparison, 2]], Higher_Group = higher_group, stringsAsFactors = FALSE))
    }
  }
}

# Write results to CSV
write.csv(results_to_write, "Metaphlan_WILCOX_Loc_eco.csv", row.names = FALSE)

# Print the written results
print(results_to_write)


anova_results_Loc_eco <- read.csv("Metaphlan_ANOVA_Loc_eco.csv", as.is = TRUE)



# Display results


wilcox_results <- read.csv("Metaphlan_WILCOX_Loc_eco.csv", as.is = TRUE)
anova_results_Ecotype <- read.csv("Metaphlan_ANOVA_Loc_eco.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = wilcox_results, anova_results_Ecotype, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge2, "Final_Loc_eco.csv")




















rm(list=ls())

# Set working directory
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_HUMANN/HUMANN_ANOVA")

# Read the CSV file (adjust the path as needed)
data <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Population <- as.factor(data$Population)
data$Loc_eco <- as.factor(data$Loc_eco)



data$Population <- as.factor(data$Population)

# Select path columns for comparison
path_columns <- names(data)[6:ncol(data)]  # Select Path columns

# Create a results data frame for Wilcoxon results
wilcox_results <- data.frame(Metabolite = character(),
                             Comparison = character(),
                             W = numeric(),
                             p_value = numeric(),
                             Higher = character(),
                             stringsAsFactors = FALSE)

# Define all combinations of Populations
Populations <- levels(data$Population)

# Perform Wilcoxon tests for each path and each pair of Populations
for (path in path_columns) {
  for (i in 1:(length(Populations) - 1)) {
    for (j in (i + 1):length(Populations)) {
      loc1 <- Populations[i]
      loc2 <- Populations[j]
      
      # Extract values for each Population
      loc1_values <- data[[path]][data$Population == loc1]
      loc2_values <- data[[path]][data$Population == loc2]
      
      # Check if both groups have enough data for the test
      if (length(loc1_values) > 0 && length(loc2_values) > 0) {
        # Perform Wilcoxon test
        wilcox_test <- wilcox.test(loc1_values, loc2_values)
        
        # Determine which Population has the higher median
        higher_Population <- ifelse(median(loc1_values) > median(loc2_values), loc1, loc2)
        
        # Store results
        wilcox_results <- rbind(wilcox_results, data.frame(
          Metabolite = path,
          Comparison = paste(loc1, "vs", loc2),
          W = wilcox_test$statistic,
          p_value = wilcox_test$p.value,
          Higher = higher_Population
        ))
      }
    }
  }
}


wilcox_results$p_adjusted <- p.adjust(wilcox_results$p_value, method = "fdr")
# Save Wilcoxon results to CSV
write.csv(wilcox_results, file = "wilcox_results_Populations.csv", row.names = FALSE)

# Display results
print(wilcox_results)





wilcox_results <- read.csv("HUMANN_WILCOX_Ecotype.csv", as.is = TRUE)
anova_results_Ecotype <- read.csv("HUMANN_ANOVA_Population.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = For_merging, anova_results_Ecotype, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge, "Merge_Population.csv")




##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
#########Loc_Eco ANOVA#################
####PAIRWISE Wilcox
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################


rm(list=ls())

# Set working directory
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_HUMANN/HUMANN_ANOVA")

# Read the CSV file (adjust the path as needed)
data <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Population <- as.factor(data$Population)
data$Loc_eco <- as.factor(data$Loc_eco)

# Select metabolite columns (assuming they are all after the first three columns)
metabolite_columns <- names(data)[6:ncol(data)]


summarise(p_value = wilcox.test(metabolite_columns$Population, Population, exact = FALSE)$p.value)

# Create a results data frame to store ANOVA results
anova_results <- data.frame(Metabolite = character(), 
                            p_value = numeric(),
                            stringsAsFactors = FALSE)

# Create a results data frame for significant metabolites with Wilcoxon results
final_results <- data.frame(Metabolite = character(),
                            Comparison = character(),
                            p_adjusted = numeric(),
                            stringsAsFactors = FALSE)

# Perform ANOVA and Wilcoxon test for each metabolite
for (metabolite in metabolite_columns) {
  # Remove rows with NA values for the current metabolite
  cleaned_data <- data %>% filter(!is.na(data[[metabolite]]))
  
  # Check if there are enough data points for ANOVA
  if (nrow(cleaned_data) >= 3) {  # At least 3 points needed for ANOVA
    # Perform ANOVA
    model <- aov(cleaned_data[[metabolite]] ~ Loc_eco, data = cleaned_data)
    summary_model <- summary(model)
    p_value <- summary_model[[1]][["Pr(>F)"]][1]  # Extract p-value for the Ecotype factor
    
    # Store the ANOVA results
    anova_results <- rbind(anova_results, data.frame(Metabolite = metabolite, p_value = p_value))}}


write.csv(anova_results, file = "HUMANN_ANOVA_Loc_eco.csv", row.names = FALSE)



# Display results


wilcox_results <- read.csv("HUMANN_WILCOX_Ecotype.csv", as.is = TRUE)
anova_results_Ecotype <- read.csv("HUMANN_ANOVA_Ecotype.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = wilcox_results, anova_results_Ecotype, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge2, "Merge_Ecotype.csv")






rm(list=ls())

# Set working directory
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_HUMANN/HUMANN_ANOVA")

# Read the CSV file (adjust the path as needed)
df <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Loc_eco <- as.factor(data$Loc_eco)


# Select path columns for comparison
path_columns <- names(data)[6:ncol(data)]  # Select Path columns

# Sample data

# Step 1: Perform ANOVA for each metabolite
metabolite_names <- colnames(df)[grepl("Path", colnames(df))]  # Identify metabolite columns
anova_results <- list()

for (metabolite in metabolite_names) {
  anova_model <- aov(as.formula(paste(metabolite, "~ Population")), data = df)
  anova_summary <- summary(anova_model)
  anova_results[[metabolite]] <- anova_summary[[1]]$`Pr(>F)`[1]  # Extract p-value
}

# Print ANOVA results
anova_results_df <- data.frame(Metabolite = metabolite_names, P_Value = unlist(anova_results))
print(anova_results_df)

write.csv(anova_results_df, file = "HUMANN_ANOVA_Population.csv", row.names = FALSE)


# Step 2: Pairwise Wilcoxon test for significant ANOVA results
wilcox_results <- list()

for (metabolite in metabolite_names) {
  p_value <- anova_results[[metabolite]]
  if (p_value < 0.05) {  # Check if ANOVA is significant
    wilcox_test <- pairwise.wilcox.test(df[[metabolite]], df$Population, p.adjust.method = "fdr")
    wilcox_results[[metabolite]] <- wilcox_test
  }
}

# Print Wilcoxon results
for (metabolite in names(wilcox_results)) {
  cat("\nMetabolite:", metabolite, "\n")
  print(wilcox_results[[metabolite]])
}

# Optional: Identify significantly higher groups and print results
for (metabolite in names(wilcox_results)) {
  cat("\nMetabolite:", metabolite, "\n")
  comparisons <- wilcox_results[[metabolite]]$p.value
  significant_comparisons <- which(comparisons < 0.05, arr.ind = TRUE)
  
  for (comparison in seq_len(nrow(significant_comparisons))) {
    group1 <- rownames(comparisons)[significant_comparisons[comparison, 1]]
    group2 <- colnames(comparisons)[significant_comparisons[comparison, 2]]
    cat(paste("Significant difference between:", group1, "and", group2, "\n"))
  }
}

# Visualization (Optional)
for (metabolite in metabolite_names) {
  if (metabolite %in% names(wilcox_results)) {
    ggboxplot(df, x = "Population", y = metabolite, 
              title = paste("Boxplot for", metabolite),
              ylab = metabolite) + 
      stat_compare_means(method = "wilcox.test", comparisons = wilcox_results[[metabolite]]$p.value)
  }
}



HUMANN_ANOVA_Population



wilcox_results <- read.csv("HUMANN_WILCOX_Ecotype.csv", as.is = TRUE)
HUMANN_ANOVA_Population <- read.csv("HUMANN_ANOVA_Population.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = For_merging, HUMANN_ANOVA_Population, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge, "HUMANN_ANOVA_Population_Merge.csv")

















library(dplyr)
library(tidyr)
library(ggpubr)


rm(list=ls())

# Set working directory
setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_HUMANN/HUMANN_ANOVA")

# Read the CSV file (adjust the path as needed)
df <- read.csv("summarized_data.csv")  # Replace with your actual file path

# Convert Ecotype to a factor
data$Location <- as.factor(data$Location)
data$Ecotype <- as.factor(data$Ecotype)
data$Loc_eco <- as.factor(data$Loc_eco)



# Step 1: Perform ANOVA for each metabolite
metabolite_names <- colnames(df)[grepl("Path", colnames(df))]  # Identify metabolite columns
anova_results <- list()

for (metabolite in metabolite_names) {
  anova_model <- aov(as.formula(paste(metabolite, "~ Loc_eco")), data = df)
  anova_summary <- summary(anova_model)
  anova_results[[metabolite]] <- anova_summary[[1]]$`Pr(>F)`[1]  # Extract p-value
}

# Print ANOVA results
anova_results_df <- data.frame(Metabolite = metabolite_names, P_Value = unlist(anova_results))
print(anova_results_df)

write.csv(anova_results_df, file = "HUMANN_ANOVA_Loc_eco.csv", row.names = FALSE)



HUMANN_ANOVA_Loc_eco <- read.csv("HUMANN_ANOVA_Loc_eco.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = For_merging, HUMANN_ANOVA_Loc_eco, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge, "HUMANN_ANOVA_Loc_eco_Merge.csv")






# Step 2: Pairwise Wilcoxon test for significant ANOVA results
wilcox_results <- list()

for (metabolite in metabolite_names) {
  p_value <- anova_results[[metabolite]]
  if (p_value < 0.05) {  # Check if ANOVA is significant
    wilcox_test <- pairwise.wilcox.test(df[[metabolite]], df$Loc_eco, p.adjust.method = "fdr")
    wilcox_results[[metabolite]] <- wilcox_test
  }
}

# Print Wilcoxon results
for (metabolite in names(wilcox_results)) {
  cat("\nMetabolite:", metabolite, "\n")
  print(wilcox_results[[metabolite]])
}



r
Copy code
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggpubr)

# Sample data
# df <- read.csv("your_data_file.csv")  # Replace with your actual data file

# Step 1: Perform ANOVA for each metabolite
metabolite_names <- colnames(df)[grepl("Path", colnames(df))]
anova_results <- list()

for (metabolite in metabolite_names) {
  anova_model <- aov(as.formula(paste(metabolite, "~ Loc_eco")), data = df)
  anova_summary <- summary(anova_model)
  anova_results[[metabolite]] <- anova_summary[[1]]$`Pr(>F)`[1]
}


write.csv(anova_results, "anova_results.csv", row.names = FALSE)


# Step 2: Pairwise Wilcoxon test for significant ANOVA results
wilcox_results <- list()
results_to_write <- data.frame(Metabolite = character(), Group1 = character(), Group2 = character(), P_Value = numeric(), Higher_Group = character(), stringsAsFactors = FALSE)

for (metabolite in metabolite_names) {
  p_value <- anova_results[[metabolite]]
  if (p_value < 0.05) {
    wilcox_test <- pairwise.wilcox.test(df[[metabolite]], df$Population, p.adjust.method = "fdr")
    wilcox_results[[metabolite]] <- wilcox_test
    
    comparisons <- wilcox_test$p.value
    significant_comparisons <- which(comparisons < 0.05, arr.ind = TRUE)
    
    for (comparison in seq_len(nrow(significant_comparisons))) {
      group1 <- rownames(comparisons)[significant_comparisons[comparison, 1]]
      group2 <- colnames(comparisons)[significant_comparisons[comparison, 2]]
      
      # Calculate medians
      median_group1 <- median(df[df$Population == group1, metabolite], na.rm = TRUE)
      median_group2 <- median(df[df$Population == group2, metabolite], na.rm = TRUE)
      
      # Determine which group is higher
      higher_group <- ifelse(median_group1 > median_group2, group1, group2)
      
      results_to_write <- rbind(results_to_write, data.frame(Metabolite = metabolite, Group1 = group1, Group2 = group2, P_Value = comparisons[significant_comparisons[comparison, 1], significant_comparisons[comparison, 2]], Higher_Group = higher_group, stringsAsFactors = FALSE))
    }
  }
}

# Write results to CSV
write.csv(results_to_write, "wilcox_results.csv", row.names = FALSE)

# Print the written results
print(results_to_write)




wilcox_results <- read.csv("wilcox_results.csv", as.is = TRUE)


Merge = full_join(x = For_merging, wilcox_results, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge, "wilcox_results_Merge.csv")










# Optional: Identify significantly higher groups and print results
for (metabolite in names(wilcox_results)) {
  cat("\nMetabolite:", metabolite, "\n")
  comparisons <- wilcox_results[[metabolite]]$p.value
  significant_comparisons <- which(comparisons < 0.05, arr.ind = TRUE)
  
  for (comparison in seq_len(nrow(significant_comparisons))) {
    group1 <- rownames(comparisons)[significant_comparisons[comparison, 1]]
    group2 <- colnames(comparisons)[significant_comparisons[comparison, 2]]
    cat(paste("Significant difference between:", group1, "and", group2, "\n"))
  }
}

# Visualization (Optional)
for (metabolite in metabolite_names) {
  if (metabolite %in% names(wilcox_results)) {
    ggboxplot(df, x = "Loc_eco", y = metabolite, 
              title = paste("Boxplot for", metabolite),
              ylab = metabolite) + 
      stat_compare_means(method = "wilcox.test", comparisons = wilcox_results[[metabolite]]$p.value)
  }
}



HUMANN_ANOVA_Loc_eco



wilcox_results <- read.csv("HUMANN_WILCOX_Ecotype.csv", as.is = TRUE)
HUMANN_ANOVA_Loc_eco <- read.csv("HUMANN_ANOVA_Loc_eco.csv", as.is = TRUE)

For_merging <- read.csv("For_merging.csv", as.is = TRUE)

Merge = full_join(x = For_merging, HUMANN_ANOVA_Loc_eco, by = "Metabolite")
Merge2 = full_join(x = Merge, For_merging, by = "Metabolite")
write.csv(Merge, "HUMANN_ANOVA_Loc_eco_Merge.csv")
































install.packages("openxlsx")
rm(list=ls())

library(openxlsx)

setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_Humann_MAGS")


# Load the Excel file
data <- read.xlsx("Combined_24.xlsx")  # Update with your file name


# Initialize an empty data frame to store the modified data
new_data <- data[1, ]  # Start with the first row

# Loop through the data starting from the second row
for (i in 1:nrow(data)) {
  new_data[[length(new_data) + 1]] <- data[i, ]  # Add the current row
  
  # If it's not the first row and the current row does not match the previous one
  if (i > 1 && !all(data[i, ] == data[i - 1, ])) {
    new_data[[length(new_data) + 1]] <- rep(NA, ncol(data))  # Add a blank row
  }
}


# Combine the list into a data frame
new_data <- do.call(rbind, new_data)

# Save the new data to an Excel file
write.xlsx(new_data, "Combined_24_Modified_forDESEQ.xlsx", rowNames = FALSE)  # Replace with your desired file name



setwd("/Users/user/Documents/02_Andropogon/03_SPIPS_Metagenome/00_Humann_MAGS")


# Load the Excel file
data <- read.xlsx("Combined_24.xlsx")


# Initialize an empty list to store the new data

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(openxlsx)
library(tidyr)
library(dplyr)
library(readr)
library(ape)

df <- read.xlsx("Combined_24.xlsx")  # Update with your file name



result <- df %>%
  group_by(PATH) %>%
  do(bind_rows(., tibble(PATH = NA))) %>%
  ungroup()


write.xlsx(result, "Combined_24_Modified_forDESEQ_WITHSPACES.xlsx", rowNames = FALSE)  # Replace with your desired file name

data <- read.xlsx("Combined_24_Modified_forDESEQ_WITHSPACES.xlsx")

summarized_data <- data %>%
  group_by(PATH) %>%
  summarize(across(everything(), sum, na.rm = TRUE), .groups = 'drop')

write.xlsx(summarized_data, "summarized_data.xlsx", rowNames = FALSE)  # Replace with your desired file name



  
  