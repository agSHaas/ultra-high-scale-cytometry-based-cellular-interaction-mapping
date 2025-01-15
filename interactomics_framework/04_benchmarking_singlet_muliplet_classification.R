######################################
##                                  ##
##  BENCHMARKING SINGLET/MULTIPLET  ##
##         CLASSIFICATION           ##
##                                  ##
##  GROUNDTRUTH DATA: CYTOSTIM+     ##
##        S8 IMAGING DATA           ##
##               HD                 ##
##                                  ##
##                                  ##
##        Viktoria, 2024            ##
##                                  ##      
######################################

# SETUP -------------------------------------------------------------------------
# packages 
library(PICtR)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggmagnify)
library(data.table)
library(stringr)
library(autothresholdr) 
library(bluster)
library(data.table)
library(Spectre)
library(flowCore)
library(CytoML) 
library(flowWorkspace)
library(mclust)
library(dbscan)
library(FlowSOM)
library(writexl)
library(rstatix)

# directories
dir <- "./interactomics_framework/"
res.dir <- paste0(dir, "results/")
data.dir <- paste0(dir, "data/")
plot.dir <- paste0(dir, "plots/")

# create directory structure if not already present 
for (d in c(data.dir, res.dir, plot.dir)) {
  ifelse(!dir.exists(d), dir.create(d), FALSE)
}

set.seed(42)

# functions
triangle_thresholding <- function(hist) {
  # code ported from ImageJ's implementation of the triangle algorithm in java 
  # initialise variables 
  min <- 0
  max <- 0 
  min2 <- 0
  min3 <- 0
  dmax <- 0 
  threshold <- 0 
  
  # counts and breaks of the histogram 
  counts <- hist$counts
  breaks <- hist$breaks 
  
  # find first non-zero (minimum) bin of the histogram 
  for (i in seq_along(counts)) {
    if (counts[i] > 0) {
      min <- i
      break
    }
  }
  # move one step back to the first zero point before min
  # note 1-based indexing in R 
  if (min > 1) {
    min <- min - 1
  }
  
  # find last non-zero (minimum) bin of the histogram
  for (i in seq_along(counts)) {
    if (counts[length(counts) - i + 1] > 0) {
      min2 <- length(counts) - i + 1
      break
    }
  }
  # Move one step forward to the first zero point after min2
  # note 1-based indexing in R 
  if (min2 < length(counts)) {
    min2 <- min2 + 1
  }
  
  
  # find the peak 
  for (i in seq_along(counts)) {
    if (counts[i] > dmax) {
      max <- i
      dmax <- counts[i]
    }
  }
  
  
  # two possible thresholds for the two sides of the histogram - find the side
  # where the distance of the peak to the minumun is furthest 
  inverted <- FALSE
  
  if ((max - min) < (min2 - max)) {
    # Reverse the histogram
    inverted <- TRUE
    counts <- rev(counts)
    min <- length(counts) - min2 + 1
    max <- length(counts) - max + 1
  }
  
  # If min and max are the same, return min
  if (min == max) {
    return(min)
  }
  
  # describe the line from the peak of the histogram to the minimum 
  nx <- counts[max] # peak frequency
  ny <- min - max
  d <- sqrt(nx^2 + ny^2)
  nx <- nx / d
  ny <- ny / d
  d <- nx * min + ny * counts[min]
  
  # find the point with the maximum distance from the line connecting peak 
  # and minimum to the histogram 
  threshold <- min # initialize
  distance <- 0
  
  for (i in (min + 1):max) {
    new_distance <- nx * i + ny * counts[i] - d
    if (new_distance > distance) {
      threshold <- i
      distance <- new_distance
    }
  }
  
  # Adjust the split point
  threshold <- threshold - 1
  
  # reverse the histogram / threshold back if needed
  if (inverted) {
    threshold <- length(counts) - threshold + 1
  }
  
  return(breaks[threshold])
}
cm_heatmap <- function(df, pred, xlab, actual = "status_broad", ylab = "actual", sample = "element_name") {
  p <- ggplot(df, aes(x = .data[[pred]], y = .data[[actual]], fill = n)) + 
    geom_tile() + 
    facet_wrap(~ .data[[sample]]) + 
    theme_bw() + 
    coord_equal() + 
    guides(fill = "none") +
    labs(x = xlab, y = ylab) +
    geom_text(aes(label = n), size = 10) +
    scale_fill_distiller(direction = 1) +
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), axis.text.y = element_text(angle = 90))
}
f1_score <- function(df, pred, actual = "status_broad", pos = "multiplet", neg = "singlet", sample = "element_name"){
  res <- data.frame()
  
  for (i in unique(df[[sample]])) {
    # extract TP, TN, FP, FN for each sample
    df_n <- df %>% dplyr::filter(.data[[sample]] == i)
    tp <- df_n %>% dplyr::filter(.data[[actual]] == pos & .data[[pred]] == pos) %>% summarise(sum(n)) %>% pull()
    fp <- df_n %>% dplyr::filter(.data[[actual]] == neg & .data[[pred]] == pos) %>% summarise(sum(n)) %>% pull()
    tn <- df_n %>% dplyr::filter(.data[[actual]] == neg & .data[[pred]] == neg) %>% summarise(sum(n)) %>% pull()
    fn <- df_n %>% dplyr::filter(.data[[actual]] == pos & .data[[pred]] == neg) %>% summarise(sum(n)) %>% pull()
    
    # F1 score
    f1 <- (2*tp) / (2*tp+fp+fn)
    
    # add to results df
    res <- rbind(res, data.frame(pred, i, tp, fp, tn, fn, f1))
  }
  return(res)
}

# DATA -------------------------------------------------------------------------
# load data (S8 data with manual annotations for singlet/multiplets based on images = ground truth)
all <- fread("./interactomics_framework/data/classified_data_CSposCSneg.csv")

# colnames that are easier to handle
all <- all %>% rename(T_APCs = APC_T)
colnames(all) <- colnames(all) %>%
  str_replace_all("-", "_") %>%
  str_replace_all("\\s\\(", "_") %>%
  str_replace_all("\\)", "") %>%
  str_replace_all("\\s", "_")

# remove undefined and coincident events from the start 
# subset on CS+ cause that's what's done in Figure 1
data <- all %>% 
  dplyr::filter(status_broad %in% c("singlet", "multiplet")) %>% 
  dplyr::filter(!str_detect(element_name, "ctrl"))

# OTSU -------------------------------------------------------------------------
otsu_threshold <- calculateThreshold(hist(data$FSC_ratio, breaks = 1000, plot = F))
data <- data %>% 
  mutate(otsu = if_else(FSC_ratio >= otsu_threshold, "multiplet", "singlet"))

# stats
otsu_res <- data %>%
  group_by(element_name, status_broad, otsu) %>% 
  count()

# plot confusion heatmaps
cm_otsu <- cm_heatmap(otsu_res, "otsu", "otsu")
cm_otsu

# plot Otsu
pdf(paste0(plot.dir, Sys.Date(), "_otsu_threshold_histogram.pdf"), height = 5, width = 6)
ggplot(data, aes(x = FSC_ratio, fill = otsu)) + 
  geom_histogram(bins = 1000) + 
  geom_vline(xintercept = otsu_threshold, color = "black", linetype = "dashed") + 
  theme_classic() + 
  scale_fill_manual(values = c("#264896", "#D6D6D6"))
dev.off()

# plot otsu colored by ground truth, nudged 
pdf(paste0(plot.dir, Sys.Date(), "_otsu_threshold_histogram_groundtruth_nudged.pdf"), height = 5, width = 6)
p <- ggplot(data %>% mutate(status_fine = factor(status_fine, levels = c("singlet", "doublet", "triplet", "quadruplet_more"))), 
       aes(x = FSC_ratio, fill = status_fine)) + 
  geom_histogram(bins = 1000, position = "nudge") + 
  geom_vline(xintercept = otsu_threshold, color = "black", linetype = "dashed") + 
  theme_classic() + 
  scale_fill_manual(values = c("#D6D6D6", "#F2911A", "#A490BF", "#4F80AF")) 
print(p)
dev.off()

# export source data for Fig1D
write_xlsx(p$data, path = paste0(plot.dir, "Figure1D_sourcedata.xlsx"), col_names = T, format_headers = T)

# K-MEANS CLUSTERING ----------------------------------------------------------
# just on the FSC ratio 
set.seed(42)
km_ratio <- kmeans(data$FSC_ratio, centers = 2)
data$km_clust_ratio <- as.factor(km_ratio$cluster)

data <- data %>% mutate(km_clust_ratio = if_else(km_clust_ratio == 1, "singlet", "multiplet"))

# plot
ggplot(data, aes(x = FSC_ratio, fill = km_clust_ratio)) + geom_histogram(bins = 1000)
ggplot(data, aes(x = FSC_ratio, fill = km_clust_ratio)) + geom_histogram(bins = 1000) + geom_vline(xintercept = otsu_threshold) # almost corresponds to otsu threshold 

# confusion matrix + heatmap 
km_ratio_res <- data %>% 
  group_by(element_name, status_broad, km_clust_ratio) %>% 
  count() 

cm_km_ratio <- cm_heatmap(km_ratio_res, pred = "km_clust_ratio", xlab = "kmeans FSC_ratio")
cm_km_ratio

# TRIANGLE THRESHOLDING -------------------------------------------------------
triangle_threshold <- triangle_thresholding(hist(data$FSC_ratio, breaks = 1000, plot = F))
data <- data %>% 
  mutate(triangle = if_else(FSC_ratio >= triangle_threshold, "multiplet", "singlet"))

# stats
triangle_res <- data %>%
  group_by(element_name, status_broad, triangle) %>% 
  count()

# plot confusion heatmap
cm_triangle <- cm_heatmap(triangle_res, "triangle", "triangle")
cm_triangle

# MANUAL GATING ---------------------------------------------------------------
# n = 3 people in FlowJo 
ws <- open_flowjo_xml(paste0(data.dir, "CS_all_range_manual_gating_benchmark.wsp"), 
                      sample_names_from="sampleNode") # load workspace

# extract gatingset from group called "samples"
gs <- flowjo_to_gatingset(ws, name = "samples")

# plot gating tree to double-check 
plot(gs)

# create list with all replicates 
cytoframe_list <- list() 
for (i in 1:9) {
  replicate <- get_cytoframe_from_cs(gs, i)
  replicate <- as_tibble(exprs(replicate))
  replicate <- replicate %>% mutate(imageid = paste0("id", sprintf("%0.8d", 0:(dim(replicate)[1]-1))))
  
  cytoframe_list[[i]] <- replicate
}
names(cytoframe_list) <- sampleNames(gs)

# Generate gating results list which contains lists containing image id and one column with FALSE/TRUE for each gate 
gating_results <- list()
for (i in 1:9) {
  results <- gs_get_pop_paths(gs, path = "auto") %>%
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

# compress the lists stored within gating results to columns containing image_id und gating indices containing false true
gating_indices <- list()
for (i in seq_along(gating_results)) {
  reduced <- purrr::reduce(gating_results[[i]], inner_join)
  gating_indices[[i]] <- reduced
}
names(gating_indices) <-sampleNames(gs)

# merging gating information 
complete_gated_list <- purrr::map2(cytoframe_list, gating_indices, ~ left_join(.x, .y, by = "imageid"))
names(complete_gated_list) <- sampleNames(gs)
names(complete_gated_list)

# combine and add element_name column
complete_gated <- bind_rows(complete_gated_list, .id = "element_name")
complete_gated <- complete_gated %>% 
  mutate(uniq_id = str_c(element_name, imageid, sep = "_")) %>% 
  select(c(uniq_id, dbts_1, dbts_2, dbts_3)) %>% 
  mutate(across(2:4, ~ if_else(.x, "multiplet", "singlet"))) %>% 
  rename(manual_1 = dbts_1, manual_2 = dbts_2, manual_3 = dbts_3)

# merge with data 
data <- left_join(data, complete_gated, by = "uniq_id")

# plot confusion heatmaps for each person and calculate F1 score
res_man <- data.frame() 
pdf(paste0(res.dir, Sys.Date(), "_manual_gating_confusion_heatmaps.pdf"), height = 5, width = 5, onefile = T)
for (man in c("manual_1", "manual_2", "manual_3")) { 
  # confusion heatmap
  res <- data %>% group_by(element_name, status_broad, .data[[man]]) %>% count()
  plot(cm_heatmap(res, man, man))
  
  # add F1 score to results df
  res_man <- rbind(res_man, f1_score(res, pred = man))
}
dev.off() 

# IMAGEJ THRESHOLDING METHODS -------------------------------------------------
thresh_methods <- c("Huang", "Intermodes", "IsoData",
                    "Li", "Mean", "RenyiEntropy", 
                    "Shanbhag")

# initialize empty results df
res_imageJ <- data.frame() 
pdf(paste0(res.dir, Sys.Date(), "_imageJ_thresholds_confusion_heatmaps.pdf"), height = 5, width = 5, onefile = T)
for (t in thresh_methods) {
  # find threshold and classify accordingly
  threshold <- auto_thresh(as.integer(data$FSC_ratio), method = t)[1]
  data <- data %>% 
    mutate(!!sym(t) := factor(if_else(FSC_ratio >= threshold, "multiplet", "singlet")))
  
  # confusion heatmap 
  res <- data %>% group_by(element_name, status_broad, .data[[t]]) %>% count()
  plot(cm_heatmap(res, t, t))
  
  # add F1 score to results df
  res_imageJ <- rbind(res_imageJ, f1_score(res, pred = t))
}
dev.off()

# GAUSSIAN MIXTURE MODEL ON FLOW PARAMETERS  ----------------------------------------
library(mclust)

# only on FSC-ratio 
set.seed(42)
gmm_ratio <- Mclust(data$FSC_ratio, G = 2)  # G=2 specifies two clusters
plot(gmm_ratio, what = "classification")

data <- data %>% 
  mutate(gmm_ratio_class = as.factor(gmm_ratio$classification)) %>% 
  mutate(gmm_ratio_class = if_else(gmm_ratio_class == 1, "singlet", "multiplet"))

table(data$gmm_ratio_class, data$status_broad)
ggplot(data, aes(x = FSC_A, y = FSC_H, color = gmm_ratio_class)) + geom_point()

gmm_ratio_res <- data  %>% 
  group_by(element_name, status_broad, gmm_ratio_class) %>% 
  count() 

cm_gmm_ratio <- cm_heatmap(gmm_ratio_res, pred = "gmm_ratio_class", xlab = "gmm_ratio")
cm_gmm_ratio

###############
## COMPARISON
###############

# F1 scores 
f1 <- as.data.frame(rbind(f1_score(otsu_res, "otsu"),
                          f1_score(km_ratio_res, "km_clust_ratio"), 
                          f1_score(triangle_res, "triangle"), 
                          f1_score(gmm_ratio_res, "gmm_ratio_class")))

# merge with F1 scores from ImageJ thresholding and manual gating 
f1 <- rbind(f1, res_imageJ)
f1 <- rbind(f1, res_man %>% group_by(i) %>% mutate(pred = "manual") %>% mutate_at(vars(tp, fp, tn, fn, f1), ~ mean(.)) %>% distinct()) # for manual gating, calculate the mean for each replicate across people 


f1 <- f1 %>% 
  mutate(approach = case_when(
    str_detect(pred, "km|gmm") ~ "unsupervised learning", 
    str_detect(pred, "manual") ~ "manual gating", 
    str_detect(pred, "cluster_type") ~ "louvain",
    .default = "image/histogram thresholding")) %>% 
  mutate(pred_unique = if_else(str_detect(pred, "manual"), str_extract(pred, "(.*)_[0-9]", group = 1), pred)) 

# plot 
levels <- (f1 %>% group_by(pred) %>% mutate(mean = mean(f1)) %>% arrange(desc(mean)))$pred %>% unique()
ggbarplot(f1 %>% mutate(pred = factor(pred, levels = levels)), x = "pred", y = "f1", fill = "approach", add = "mean_sd") + 
  ylim(0, 1) + 
  xlab("Classification Method") + ylab("F1 Score") + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1))

# stats -----
f1 %>% group_by(pred) %>% shapiro_test(f1) # normality can be assumed, but n = 4

f1 %>% anova_test(dv = f1, wid = i, within = pred) # signif 

# pairwise t test 
stats <- f1 %>% 
  pairwise_t_test(f1 ~ pred, ref.group = "otsu", p.adjust.method = "holm", paired = T) %>% 
  add_y_position()

# plots -----
# manuscript plot with diff thresholding methods, not compared to clustering (Extended Data Figure 1A)
levels <- c("otsu", "manual", "km_clust_ratio", "IsoData", "Intermodes", "gmm_ratio_class", "RenyiEntropy", "triangle", "Li", "Shanbhag", "Huang", "Mean") # put otsu in front for comparisons  
stats <- stats %>% mutate(group2 = factor(group2, levels = levels)) %>% arrange(group2) %>% add_y_position(step.increase = 0.3)

# only signif p values 
pdf(paste0(plot.dir, Sys.Date(), "_thresholding_methods_benchmark.pdf"), height = 7, width = 6)
p <- ggbarplot(f1 %>% mutate(pred = factor(pred, levels = levels)), x = "pred", y = "f1", fill = "approach", 
          add = c("mean_sd"),  error.plot = "errorbar", add.params = list(color = "black", shape = 21)) + 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0), limits = c(0, NA)) + 
  xlab("Classification Method") + ylab("F1 Score") + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) + 
  # scale_fill_manual(values = c())
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", color = "black", width = 0.5, size = 0.5) + 
  stat_pvalue_manual(stats, label = "p.adj", hide.ns = T) + 
  geom_jitter(aes(fill = approach), width = 0.2, height = 0, size = 3, shape = 21,  color = "black") + 
  scale_fill_manual(values = c("#BE0032", "#7EA885", "#608DA2"))
print(p)
dev.off()

# export source data
write_xlsx(p$data, paste0(plot.dir, "/ExtendedDataFigure1A_sourcedata.xlsx"))

saveRDS(data, paste0(data.dir, Sys.Date(), "_data_benchmarking.rds"))