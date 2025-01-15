# SETUP --------------------------------------------------
library(flowCore)
library(CytoML)
library(flowWorkspace)
library(tidyverse)
library(ggrepel)
library(pals)
library(dplyr)
library(purrr) 
library(PICtR)

# the following script serves to important image-enabled flow cytometry data 
# from a flowjo workspace into R and pre-process them as described by and 
# adapted from Schraivogel et al. (2022) 10.1126/science.abj3013, see also
# https://github.com/benediktrauscher/ICS/blob/master/vignettes/flow_data_processing.Rmd

# Authors: Dominik Vonficht, Viktoria Flore, adapted from Daniel Schraivogel

# EXTRACT DATA -----------------------------------------------------------------
dir <- "./interactomics_framework/"
data.dir <- "./interactomics_framework/data/"

# create directory structure if not already present 
for (d in c(data.dir)) {
  ifelse(!dir.exists(d), dir.create(d), FALSE)
}

wsfile <- list.files(data.dir, pattern="CS_all_range.wsp",full=TRUE)
ws <- open_flowjo_xml(wsfile,sample_names_from="sampleNode")

# extract gatingset from group called "samples"
gs <- flowjo_to_gatingset(ws, name = "samples")

# gating tree
plot(gs)

allnodes <- gs_get_pop_paths(gs, path = "auto")
allnodes_small <- allnodes[c(2:4,6,7)]

# determine number of samples
num_samples <- length(gs)

cytoframe_list <- list()  # Initialize an empty list to store cytoframes

for (i in 1:num_samples) {
  replicate <- get_cytoframe_from_cs(gs, i)
  replicate <- as_tibble(exprs(replicate))
  replicate <- replicate %>% mutate(imageid = paste0("id", sprintf("%0.8d", 0:(dim(replicate)[1]-1))))
  
  cytoframe_list[[i]] <- replicate  # Store the transformed_df in the list
}
names(cytoframe_list) <- sampleNames(gs)


# Generate gating results list which contains lists containing image id and one column with FALSE/TRUE for each gate in allnodes_small
gating_results <- list()  # Initialize an empty list

for (i in 1:num_samples) {
  results <- allnodes_small %>%
    lapply(function(x) {
      var <- as_tibble(exprs(gs_pop_get_data(gs, x)[[i]])) %>%
        mutate(new = TRUE)
      var <- left_join(cytoframe_list[[i]], var) %>%
        select(imageid, new) %>%
        mutate(new = ifelse(is.na(new), FALSE, new))
      colnames(var)[2] <- x
      return(var)
    })
  gating_results[[i]] <- results
}

names(gating_results) <-sampleNames(gs)


# compress the lists stored within gating results to columns containing image_id and gating indices containing false/true
gating_indices <- list()
for (i in seq_along(gating_results)) {
  reduced <- purrr::reduce(gating_results[[i]], inner_join)
  gating_indices[[i]] <- reduced
}

names(gating_indices) <-sampleNames(gs)


# Generate a new dataframe. merging the full dataframes stored in cytoframe_list with the newly generated dataframe
# containing the gating information from  the flowjoworkspace
complete_gated_list <- purrr::map2(cytoframe_list, gating_indices, ~ left_join(.x, .y, by = "imageid"))
names(complete_gated_list) <-sampleNames(gs)

# we then need to create tables that can be used as insert for the ICS imageJ plugin to convert .tiff files to .jpg files
# set paths to the image folders containg .tiff raw images
tiffcs1 <- "./interactomics_framework/20230616_1506_Cytostim_1/images/" 
tiffcs2 <- "./interactomics_framework/20230616_1508_Cytostim_2/images/" 
tiffcs3 <- "./interactomics_framework/20230616_1511_Cytostim_3/images/" 
tiffcs4 <-"./interactomics_framework/20230616_1517_Cytostim_4/images/" 
tiffcs5 <- "./interactomics_framework/20230616_1520_Cytostim_5/images/" 
tiffctrl1 <- "./interactomics_framework/20230616_1449_Cytostim_ctrl_1/images/" 
tiffctrl2 <- "./interactomics_framework/20230616_1454_Cytostim_ctrl_2/images/"  
tiffctrl3 <- "./interactomics_framework/20230616_1457_Cytostim_ctrl_3/images/"  
tiffctrl4 <- "./interactomics_framework/20230616_1502_Cytostim_ctrl_4/images/"  

tiff_vector <- c(tiffcs1,tiffcs2,tiffcs3,tiffcs4,tiffcs5,tiffctrl1,tiffctrl2,tiffctrl3,tiffctrl4)

# create an empty list
files_list <- list()

for (i in tiff_vector){
files <- list.files(path = i, full.names = F, pattern = ".tiff$", recursive = TRUE)
files_list[[i]] <- files
}

image_names_list <- purrr::map(files_list, ~ .x %>% as.tibble() %>% 
                          dplyr::mutate(imageid = paste0("id", str_extract(str_extract(value, "[0-9]{8}.tiff"), "[0-9]{8}"))) %>% 
                          dplyr::mutate(path = paste0("../images/", value)) %>% 
                          dplyr::select(-value))

names(image_names_list) <- sampleNames(gs)
complete_gated_list_2 <- purrr::map2(complete_gated_list, image_names_list, ~ left_join(.x, .y, by = "imageid"))


# we only keep values in which waveform <0
complete_gated_list_3 <- purrr::map(complete_gated_list_2, ~ .x %>% dplyr::filter(WaveformPresent < 0 & (complete.cases(path))))

# select specific columns for fiji
fiji_table_list <- purrr::map(complete_gated_list_3, ~ .x %>% dplyr::select(imageid,Cells,Live,`CD45+`,B_T,Myeloid_T,path))

# save elements of fiji_table_list as .csv files for processing in fiji
# Define a function to save a tibble as a CSV file
save_tibble_as_csv <- function(tibble_data, file_name) {
  fwrite(tibble_data, file = file_name)
}

# save each tibble as a CSV file
file_names <- c("CS1.csv","CS2.csv","CS3.csv","CS4.csv","CS5.csv","CS1_ctrl.csv","CS2_ctrl.csv","CS3_ctrl.csv","CS4_ctrl.csv")  # List of file names
setwd(paste0(getwd(), "/cytostim_neg_controls_VF"))
purrr::map2(fiji_table_list, file_names, ~ save_tibble_as_csv(.x, .y))

###########now run fiji plugin for individual csv files
# #paramters in fiji
# channels 1,2,4
# cs1
# minimum: 0.03,0.0,0.0
# maxmimum 0.12,0.7,0.9
#
# crop horizontal: 5


# CS1_ctrl: 27788 images, CS2_ctrl: 20267 images, CS3_ctrl: 27727 images, CS4_ctrl: 28148 images 

# create final dataframes, on which interactive analysis and cell clustering can be run on
jpg_cs1 <- "./interactomics_framework/20230616_1506_Cytostim_1/images-processed/"
jpg_cs2 <- "./interactomics_framework/20230616_1508_Cytostim_2/images-processed/"
jpg_cs3 <- "./interactomics_framework/20230616_1511_Cytostim_3/images-processed/"
jpg_cs4 <- "./interactomics_framework/20230616_1517_Cytostim_4/images-processed/"
jpg_cs5 <- "./interactomics_framework/20230616_1520_Cytostim_5/images-processed/"
jpg_ctrl1 <- "./interactomics_framework/20230616_1449_Cytostim_ctrl_1/images-processed/" 
jpg_ctrl2 <- "./interactomics_framework/20230616_1454_Cytostim_ctrl_2/images-processed/"  
jpg_ctrl3 <- "./interactomics_framework/20230616_1457_Cytostim_ctrl_3/images-processed/"  
jpg_ctrl4 <- "./interactomics_framework/20230616_1502_Cytostim_ctrl_4/images-processed/"  

jpg_vector <- c(jpg_cs1, jpg_cs2, jpg_cs3, jpg_cs4, jpg_cs5, jpg_ctrl1, jpg_ctrl2, jpg_ctrl3, jpg_ctrl4)


# get jpg file names
files_list <- list()
for (i in jpg_vector){
  files <- list.files(path = i, full.names = F, pattern = ".jpg$", recursive = TRUE)
  files_list[[i]] <- files
}

jpg_names_list <- purrr::map(files_list, ~ .x %>% as.tibble() %>% 
                                 dplyr::mutate(imageid = paste0("id", str_extract(str_extract(value, "[0-9]{8}.jpg"), "[0-9]{8}"))) %>% 
                                 dplyr::mutate(path_jpg = paste0("../images-processed/", value)) %>% 
                                 dplyr::select(-value))


# add filenames to dfs
complete_gated_list_4 <- purrr::map2(complete_gated_list_3, jpg_names_list, ~ left_join(.x, .y, by = "imageid"))


# keep only CD45+ cells that have an image
images_CS_cd45 <- map(complete_gated_list_4,~.x %>% dplyr::filter(`CD45+`=="TRUE"))
paths_id_only <- purrr::map(images_CS_cd45, ~.x %>% dplyr::mutate(path_jpg = gsub("\\.\\.", "", path_jpg)))

# Create a list of dataframes by mapping over jpg_vector and paths_id_only1
base_path <- dir

jpg_cs1 <- file.path(base_path, "20230616_1506_Cytostim_1")
jpg_cs2 <- file.path(base_path, "20230616_1508_Cytostim_2")
jpg_cs3 <- file.path(base_path, "20230616_1511_Cytostim_3")
jpg_cs4 <- file.path(base_path, "20230616_1517_Cytostim_4")
jpg_cs5 <- file.path(base_path, "20230616_1520_Cytostim_5")
jpg_ctrl1 <- file.path(base_path, "20230616_1449_Cytostim_ctrl_1")
jpg_ctrl2 <- file.path(base_path, "20230616_1454_Cytostim_ctrl_2")
jpg_ctrl3 <- file.path(base_path, "20230616_1457_Cytostim_ctrl_3")
jpg_ctrl4 <- file.path(base_path, "20230616_1502_Cytostim_ctrl_4")

jpg_vector <- c(jpg_cs1, jpg_cs2, jpg_cs3, jpg_cs4, jpg_cs5, jpg_ctrl1, jpg_ctrl2, jpg_ctrl3, jpg_ctrl4)

modified_paths_list <- map2(jpg_vector, paths_id_only, ~mutate(.y, extended_paths = paste0(.x, path_jpg)))
names(modified_paths_list) <-  names(complete_gated_list_3)

# add FSC_ratio
modified_paths_list <- purrr::map(modified_paths_list, ~.x %>% dplyr::mutate(FSC_ratio=`FSC-A`/`FSC-H`))

# add thresholds for the FSC ratio, FSC-W and FSC-A ----

thresholds_fscratio <- c()
thresholds_fscwidth <- c()
thresholds_fscarea <- c()

for (i in names(modified_paths_list)){
  thresholds_fscratio <- append(thresholds_fscratio, calculateThreshold(hist(modified_paths_list[[i]]$FSC_ratio, breaks=500)))
  thresholds_fscwidth <- append(thresholds_fscwidth, calculateThreshold(hist(modified_paths_list[[i]]$`FSC-W`, breaks=500)))
  thresholds_fscarea <- append(thresholds_fscarea, calculateThreshold(hist(modified_paths_list[[i]]$`FSC-A`, breaks=500)))
}

names(thresholds_fscratio) <- names(modified_paths_list)
names(thresholds_fscwidth) <- names(modified_paths_list)
names(thresholds_fscarea) <- names(modified_paths_list)

# add columns telling whether a cell is above or below different OTSU derived thresholds
# Function to add TRUE/FALSE columns based on the thresholds
add_threshold_columns <- function(df, threshold_fscratio, threshold_fscwidth, threshold_fscarea) {
  df %>%
    mutate(
      above_threshold_fsc_ratio = `FSC_ratio` > threshold_fscratio,
      above_threshold_fsc_width = `FSC-W` > threshold_fscwidth,
      above_threshold_fsc_area = `FSC-A` > threshold_fscarea,
      above_ecc = `Eccentricity (FSC)` > 0.8 & `Eccentricity (SSC)` > 0.8
    )
}

# Create a list of data frames with the added columns for each threshold
result_list <- pmap(
  list(df = modified_paths_list, 
       threshold_fscratio = thresholds_fscratio,
       threshold_fscwidth = thresholds_fscwidth,
       threshold_fscarea = thresholds_fscarea),
  add_threshold_columns
)

# remove 4th CytoStim+ replicate, as quality of images is very low
result_list_final <- result_list[-4]

# save 
saveRDS(result_list_final, file = "./interactomics_framework/data/data_CSposCSneg.rds")
