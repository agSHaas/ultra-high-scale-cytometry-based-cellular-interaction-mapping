# SETUP ------------------------------------------------------------------------
library(jpeg)
library(tidyverse)
library(data.table)

# the following script serves to classify the images from the image-enabled flow 
# cytometry data into singlets, doublets, triplets, quadruplets/multiplets, and 
# coincident doublets/multiplets by manual inspection.

# Authors: Dominik Vonficht, Viktoria Flore 

# RANDOMLY PICK 1000 CELLS PER REPLICATE RECORDING -----------------------------
result_list_final <- readRDS("./interactomics_framework/data/data_CSposCSneg.rds")
set.seed(42)

# Function to randomly sample 1000 rows from a dataframe
sample_1000_rows <- function(df) {
  sampled_df <- df %>%
    sample_n(1000, replace = FALSE)  # Randomly sample 1000 rows
  
  return(sampled_df)
}

# Apply the function to each dataframe in the list
sampled_dataframes_list <- lapply(result_list_final, sample_1000_rows)

combined_df_subset <- bind_rows(sampled_dataframes_list , .id = "element_name")

# Shuffle the rows of combined_df_subset
shuffled_combined_df <- combined_df_subset[sample(nrow(combined_df_subset)),]
sample_indices <- nrow(shuffled_combined_df)

# INTERACTIVE CLASSIFICATION ----------------------------------------------------

# Initialize an empty dataframe to store the results
result_df <- data.frame(image_path = character(0), annotation = character(0))

# Loop through the selected image paths
for (i in 1:sample_indices) {
  df=shuffled_combined_df
  image_path <- df$extended_paths[i]
  image_id <- df$imageid[i]
  sample <- df$element_name[i]
  
  # Load the image using the 'jpeg' library
  img <- readJPEG(image_path)
  
  # Display the image 
  plot(0:2, type = 'n', xlab = '', ylab = '', main = image_path)
  rasterImage(img, 1, 0, 3, 0.4)
  
  # Collect annotation from user 
  annotation <- readline(prompt = "Enter annotation for this image: ")
  
  # Store the result in the result dataframe
  result_df <- rbind(result_df, data.frame(sample_id=sample,imageid=image_id, 
                                           annotation = annotation))
  
  # Close the image display window
  dev.off()
}

# 0 is debris / unidentifiable, 1 is singlet, 2 is doublet, 3 is triplet, 4 quadruplet and more, 5 coincident doublet/multiplet

# results 
random_sampling <- result_df
random_sampling <- random_sampling %>% tidyr::unite("uniq_id",sample_id,imageid,sep="_",remove=F)

# remaining parameters
cell.dat <- shuffled_combined_df
cell.dat <- cell.dat %>% tidyr::unite("uniq_id",element_name,imageid,sep="_",remove=F)

# join
cell.dat <- cell.dat %>% dplyr::left_join(random_sampling %>% select(uniq_id, annotation), by = "uniq_id")

# ANNOTATION -------------------------------------------------------------
cell.dat.anno <- cell.dat %>% 
  dplyr::mutate(status_broad = case_when(
    annotation == 0 ~ "undefined",
    annotation == 1 ~ "singlet",
    annotation %in% c(2, 3, 4) ~ "multiplet",
    annotation == 5 ~ "coincident",
    TRUE ~ ""  # Use an empty string as the default value
  )) %>% 
  dplyr::mutate(status_fine = case_when(
    annotation == 0 ~ "low_quality",
    annotation == 1 ~ "singlet",
    annotation == 2 ~ "doublet",
    annotation == 3 ~ "triplet",
    annotation == 4 ~ "quadruplet_more",
    annotation == 5 ~ "conincident",
    TRUE ~ "" 
  ))

# save
fwrite(cell.dat.anno, file = "./interactomics_framework/data/classified_data_CSposCSneg.csv")