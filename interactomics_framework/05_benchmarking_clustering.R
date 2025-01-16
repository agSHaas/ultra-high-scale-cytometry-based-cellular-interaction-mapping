######################################
##                                  ##
##      BENCHMARKING CLUSTERING     ##
##                                  ##
##    GROUNDTRUTH DATA: CYTOSTIM+   ##
##          S8 IMAGING DATA         ##
##                 HD               ##
##                                  ##
##                                  ##
##      Viktoria Flore, 2024        ##
##                                  ##      
######################################

# this script benchmarks several clustering algorithms and feature spaces for 
# clustering within the Interact-omics framework.. Manually classified images 
# from image-enabled flow cytometry serve as ground-truth data. 

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
library(rstatix)
library(clue)
library(flowMeans)
library(Rclusterpp)
library(immunoClust)
library(ComplexHeatmap)
library(tibble)

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
  # initialize variables 
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

cols <- c("T cells" = "#44A998", 
          "other CD45+" = "#342980", 
          "B cells" = "#CC6677",
          "CD33high Myeloid" = "#0E7834", 
          "CD33low Myeloid" = "#509045", 
          "B*T" = "#88CCEE", 
          "My*T" = "#DDCC77", 
          "My*T*B" = "#A74592")

# DATA -------------------------------------------------------------------------
# load data (S8 data with manual annotations for singlet/multiplets based on images = ground truth)
data <- fread(paste0(data.dir, "2024-08-20_data_benchmarking.rds"))

# otsu threshold for later classifying clusters as singlet/multiplet clusters 
otsu_threshold <- calculateThreshold(hist(data$FSC_ratio, breaks = 1000, plot = F))
data <- data %>% 
  mutate(otsu = if_else(FSC_ratio >= otsu_threshold, "multiplet", "singlet"))

####################
##   CLUSTERING   ##
####################

# benchmark only on conventional flow markers including scatter and FSC_ratio
flow_param <- c("SSC_A","SSC_H","FSC_A","FSC_H","Comp_PE_Cy7_A","Comp_BV421_A","Comp_BV605_A", "Comp_APC_A","Comp_FITC_A", "FSC_ratio")

clust_bm <- data %>% 
  dplyr::select(all_of(flow_param),status_fine,status_broad,element_name,uniq_id, otsu) %>% 
  dplyr::mutate_at(vars(flow_param), ~ scale(.)) %>% # scaling and centering
  dplyr::mutate_at(vars(!flow_param), ~ factor(.))

# run UMAP
clust_bm <- run.umap(clust_bm, use.cols = flow_param, umap.seed = 42) # umap with k = 15
clust_bm <- as.data.table(clust_bm)

# what is misclassified after otsu only?
cols <- c("CD33pos" = "palegreen4", 
          "CD19pos" = "firebrick",
          "CD3pos" = "turquoise4",
          "other" = "navy")

pdf(paste0(plot.dir, Sys.Date(), "_missclassified_otsu_only_dotplot_gate.pdf"), width = 6.5, height = 4)
p <- clust_bm %>% 
  dplyr::filter(status_broad == "singlet") %>% 
  left_join(data %>% select(uniq_id, FSC_A, FSC_H) %>% rename(FSC_A_unscaled = FSC_A, FSC_H_unscaled = FSC_H), by = "uniq_id") %>% 
  mutate(pos = factor(case_when(
    Comp_BV421_A > -0.2 ~ "CD33pos", 
    Comp_BV605_A > 0.6 ~ "CD19pos",
    Comp_PE_Cy7_A > 0 ~ "CD3pos", 
    .default = "other"))) %>%
  ggplot(aes(x = FSC_A_unscaled, y = FSC_H_unscaled, color = pos)) + geom_point() + 
  scale_color_manual(values = cols) + 
  geom_abline(slope = 1/otsu_threshold, intercept = 0) + theme_classic()
print(p)
dev.off()

# export data for extended data figure 1B
write_xlsx(p$data, path = paste0(plot.dir, "ExtendedData_Figure1B_sourcedata_dotplot.xlsx"))

for (param in c("Comp_BV421_A", "Comp_PE_Cy7_A", "Comp_BV605_A")) { 
  print(ggscatter(clust_bm, x = param, y = "FSC_A"))
} 

# plot both correctly and misclassified events as bar plots 
pdf(paste0(plot.dir, Sys.Date(), "_classified_otsu_only_barplot_gates_correct_vs_false.pdf"), width = 10, height = 6)
p <- ggbarplot(clust_bm %>% 
                 mutate(classified = factor(if_else(otsu != status_broad, "false", "correct"))) %>% 
                 select(element_name, status_broad, uniq_id, all_of(flow_param), classified) %>% 
                 mutate(element_name = factor(element_name)) %>% 
                 mutate(pos = factor(case_when(
                   Comp_BV421_A > -0.2 ~ "CD33pos", 
                   Comp_BV605_A > 0.6 ~ "CD19pos",
                   Comp_PE_Cy7_A > 0 ~ "CD3pos", 
                   .default = "other"))) %>%
                 group_by(element_name, pos, status_broad, classified, .drop = F) %>%
                 count() %>%
                 group_by(status_broad, element_name, classified, .drop = F) %>% 
                 mutate(total = sum(n), freq = (n / sum(n))) %>% 
                 mutate(pos = factor(pos, levels = c("CD33pos", "CD19pos", "CD3pos", "other"))) %>% 
                 dplyr::filter(status_broad == "singlet"), 
               x = "classified", y = "freq", fill = "pos", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = cols) + 
  theme_classic() +
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data for Extended Data Figure 1B
write_xlsx(p$data, path = paste0(plot.dir, "ExtendedData_Figure1B_sourcedata.xlsx"))


##############################
##                          ##
##  RUN 100 ITERATIONS OF   ##
##  LOUVAIN (igraph)        ##
##  LEIDEN (igraph)         ##
##  HDBSCAN (dbscan)        ##
##  FlowSOM (Spectre)       ##
##                          ##
##############################    

f1_scores <- data.frame() # results df for f1 scores

pdf(paste0(res.dir, Sys.Date(), "_benchmarking_clustering_ratio_cluster_plots.pdf"), height = 6, width = 10, onefile = T)
for (i in 1:100) {
  message(paste0("iteration: ", i))
  
  # save random seed
  seed <- sample(1:1000000, size = 1)
  
  # louvain ----
  set.seed(seed)
  louvain_res <- clust_bm %>% 
    dplyr::select(all_of(flow_param)) %>% 
    bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = 1))) # shared nearest neighborhood graph by default 
  # add to df 
  clust_bm[[paste0("louvain_", seed)]] <- factor(louvain_res)
  
  
  # leiden ----
  set.seed(seed)
  leiden_res <- clust_bm %>% 
    dplyr::select(all_of(flow_param)) %>% 
    bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "leiden", cluster.args = list(resolution_parameter = 1, # shared nearest neighborhood graph by default 
                                                                                       n_iterations = 10, 
                                                                                       objective_function = "modularity"))) # optimize modularity 
  # add to df 
  clust_bm[[paste0("leiden_", seed)]] <- factor(leiden_res)
  
  # hdbscan on UMAP embedding ---- 
  set.seed(seed)
  hdbscan_res <- hdbscan(clust_bm %>% dplyr::select(c(UMAP_X, UMAP_Y)), minPts = 100)
  clust_bm[[paste0("hdbscan_", seed)]] <- factor(hdbscan_res$cluster) # add to df 
  
  
  # flowSOM ----
  clust_bm <- run.flowsom(clust_bm, use.cols = flow_param, meta.k = 20, clust.seed = seed, meta.seed = seed, clust.name = paste0("FlowSOM_", seed), meta.clust.name = paste0("FlowSOM_meta_", seed)) # seeds are set within function
  clust_bm <- clust_bm %>% 
    mutate_at(vars(paste0("FlowSOM_", seed), paste0("FlowSOM_meta_", seed)), ~ factor(.))

  
  for (method in c("louvain", "leiden", "hdbscan", "FlowSOM_meta")) {
    
    # calculate ratio_cluster_plot 
    method_n <- paste0(method, "_", seed)
    
    res <- clust_bm %>% 
      group_by(element_name, .data[[method_n]], otsu, .drop = F) %>% 
      count() %>% 
      group_by(.data[[method_n]], element_name, .drop = F) %>% 
      mutate(total = sum(n), freq = (n / sum(n)))
  
    # select clusters where more than 50 % of events are multiplets
    selected_clusters <- res %>% 
      dplyr::filter(otsu == "multiplet") %>% 
      group_by(.data[[method_n]]) %>% 
      summarize(mean = mean(freq)) %>% 
      dplyr::filter(mean > 0.5) %>%
      pull(.data[[method_n]])
    
    # add cluster type to df 
    clust_bm <- clust_bm %>% 
      mutate(!!sym(paste0("cluster_type_", method_n)) := if_else(.data[[method_n]] %in% selected_clusters, "multiplet", "singlet"))
    
    # plot ratio_cluster_plot and add selected multiplet clusters
    print(ggbarplot(res, x = method_n, y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
        scale_fill_manual(values = c("darkblue", "grey")) + 
        theme_classic() + 
        scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
        labs(title = method_n, subtitle = paste("multiplet clusters: ", paste(selected_clusters, collapse = ","))))
    
    # f1 score 
    confusion_mat <- clust_bm %>% group_by(element_name, status_broad, .data[[paste0("cluster_type_", method_n)]]) %>% count()
    f1 <- f1_score(confusion_mat, pred = paste0("cluster_type_", method_n))
    
    f1_scores <- rbind(f1_scores, f1)
  }  
  # break
}
dev.off()

f1_scores <- f1_scores %>% 
  mutate(method = str_extract(pred, "cluster_type_(.*)_[0-9]", group = 1))

# save clust_bm and f1_scores
saveRDS(clust_bm, paste0(data.dir, "clust_bm.rds"))
saveRDS(f1_scores, paste0(data.dir, "f1_scores.rds"))


##################################
##                              ##
##  add additional clustering   ## 
##  algorithms to benchmark     ##
##                              ##
##################################

# phenograph ----
# export data for PhenoGraph in Python
fwrite(clust_bm[,1:17], paste0(data.dir, "data_benchmark_flowparameters.csv"))

# import phenograph results from python (k=10) 
phenograph <- fread(paste0(data.dir, "data_benchmark_flowparameters_res_phenograph.csv"))
clust_bm <- left_join(clust_bm, phenograph %>% select(-c(all_of(flow_param), status_fine, status_broad, element_name,otsu,UMAP_X, UMAP_Y)), by = "uniq_id")

pdf(paste0(res.dir, Sys.Date(), "_benchmarking_clustering_ratio_cluster_plots_phenograph.pdf"), height = 6, width = 10, onefile = T)
for (col in colnames(clust_bm)) {

  if (str_detect(col, "phenograph")) {
    # calculate ratio_cluster_plot 
    method_n <- col
    
    res <- clust_bm %>% 
      group_by(element_name, .data[[method_n]], otsu, .drop = F) %>% 
      count() %>% 
      group_by(.data[[method_n]], element_name, .drop = F) %>% 
      mutate(total = sum(n), freq = (n / sum(n)))
    
    # select clusters where more than 50 % of events are multiplets
    selected_clusters <- res %>% 
      dplyr::filter(otsu == "multiplet") %>% 
      group_by(.data[[method_n]]) %>% 
      summarize(mean = mean(freq)) %>% 
      dplyr::filter(mean > 0.5) %>%
      pull(.data[[method_n]])
    
    # add cluster type to df 
    clust_bm <- clust_bm %>% 
      mutate(!!sym(paste0("cluster_type_", method_n)) := if_else(.data[[method_n]] %in% selected_clusters, "multiplet", "singlet"))
    
    # plot ratio_cluster_plot and add selected multiplet clusters
    print(ggbarplot(res, x = method_n, y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
            scale_fill_manual(values = c("darkblue", "grey")) + 
            theme_classic() + 
            scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
            labs(title = method_n, subtitle = paste("multiplet clusters: ", paste(selected_clusters, collapse = ","))))
    
    # f1 score 
    confusion_mat <- clust_bm %>% group_by(element_name, status_broad, .data[[paste0("cluster_type_", method_n)]]) %>% count()
    f1 <- f1_score(confusion_mat, pred = paste0("cluster_type_", method_n))
    f1 <- f1 %>% 
      mutate(method = str_extract(pred, "cluster_type_(.*)_[0-9]", group = 1))
    
    f1_scores <- rbind(f1_scores, f1)
  } 
}  
dev.off()

saveRDS(clust_bm, paste0(data.dir, "clust_bm_incl_phenograph.rds"))
saveRDS(f1_scores, paste0(data.dir, "f1_scores_incl_phenograph.rds"))

  
# immunoclust, flowMeans, rclusterpp
pdf(paste0(res.dir, Sys.Date(), "_benchmarking_clustering_ratio_cluster_plots_immunoclust_flowMeans_rclusterpp.pdf"), height = 6, width = 10, onefile = T)
for (i in 1:100) {
  message(paste0("iteration: ", i))
  
  # save random seed
  seed <- sample(1:1000000, size = 1)

  # immunoClust, estimated run time = 10h 
  set.seed(seed)
  immunoclust <- immunoClust::cell.MajorIterationLoop(clust_bm, parameters = flow_param,
                                                      I.buildup=6, I.final=4,
                                                      modelName="mvt2", tol=1e-5, bias=0.3, # mvt2 because of scaled values = potential cutted values at upper/lower edge
                                                      sub.bias=0.3, sub.thres=0.0, sub.tol=1e-4, sub.samples=1500,
                                                      sub.extract=0.8, sub.weights=1, sub.EM="MEt", sub.standardize=TRUE)
  # extract labels
  clust_bm[[paste0("immunoclust_", seed)]] <- factor(immunoclust@label)
  
  
  # flowMeans ----
  set.seed(seed)
  flowmeans_res <- flowMeans(clust_bm, varNames = flow_param, MaxN = NA, NumC = NA, # auto estimation of # clusters 
                             iter.max = 10, nstart = 10, Mahalanobis = TRUE, # standard Mahalanobis distance 
                             Standardize = FALSE, Update = "Mahalanobis", OrthagonalResiduals = TRUE,
                             MaxCovN = NA, MaxKernN = NA, addNoise = TRUE)
  # extract labels 
  clust_bm[[paste0("flowMeans_", seed)]] <- factor(flowmeans_res@Label)
  
  # rclusterpp ----
  set.seed(seed)
  rclusterpp_res <- Rclusterpp.hclust(clust_bm %>% dplyr::select(all_of(flow_param)), method = "ward", distance = "euclidean")
  
  # cut tree at k = 20 (over-clustering) and extract labels 
  clust_bm[[paste0("rclusterpp_", seed)]] <- factor(cutree(rclusterpp_res, k = 20))
  
  # F1 scores ----
  for (method in c("immunoclust", "flowMeans", "rclusterpp")) {
    # calculate ratio_cluster_plot 
    method_n <- paste0(method, "_", seed)
    
    res <- clust_bm %>% 
      group_by(element_name, .data[[method_n]], otsu, .drop = F) %>% 
      count() %>% 
      group_by(.data[[method_n]], element_name, .drop = F) %>% 
      mutate(total = sum(n), freq = (n / sum(n)))
    
    # select clusters where more than 50 % of events are multiplets
    selected_clusters <- res %>% 
      dplyr::filter(otsu == "multiplet") %>% 
      group_by(.data[[method_n]]) %>% 
      summarize(mean = mean(freq)) %>% 
      dplyr::filter(mean > 0.5) %>%
      pull(.data[[method_n]])
    
    # add cluster type to df 
    clust_bm <- clust_bm %>% 
      mutate(!!sym(paste0("cluster_type_", method_n)) := if_else(.data[[method_n]] %in% selected_clusters, "multiplet", "singlet"))
    
    # plot ratio_cluster_plot and add selected multiplet clusters
    print(ggbarplot(res, x = method_n, y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
            scale_fill_manual(values = c("darkblue", "grey")) + 
            theme_classic() + 
            scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
            labs(title = method_n, subtitle = paste("multiplet clusters: ", paste(selected_clusters, collapse = ","))))
    
    # f1 score 
    confusion_mat <- clust_bm %>% group_by(element_name, status_broad, .data[[paste0("cluster_type_", method_n)]]) %>% count()
    f1 <- f1_score(confusion_mat, pred = paste0("cluster_type_", method_n))
    f1 <- f1 %>% 
      mutate(method = str_extract(pred, "cluster_type_(.*)_[0-9]", group = 1))
    
    f1_scores <- rbind(f1_scores, f1)
  }  
  # break
}
dev.off()

# save all clustering results
saveRDS(clust_bm, paste0(data.dir, "clustering_benchmark_incl_all_algorithms.rds"))
saveRDS(f1_scores, paste0(data.dir, "f1_scores_incl_all_algorithms.rds"))


# plots -----
levels <- c("louvain", "phenograph_louvain", "leiden", "phenograph_leiden", "hdbscan", "FlowSOM_meta", "flowMeans", "rclusterpp", "immunoclust") 

cols <- c("louvain" = "dodgerblue3", 
          "phenograph_louvain" = "dodgerblue3", 
          "leiden" = "dodgerblue3", 
          "phenograph_leiden" = "dodgerblue3", 
          "hdbscan" = "dodgerblue3", 
          "FlowSOM_meta" = "dodgerblue3", 
          "flowMeans" = "dodgerblue3", 
          "rclusterpp" = "dodgerblue3", 
          "immunoclust" = "dodgerblue3",
          "Cytostim_1.fcs_6910" = "#B2CCE2", 
          "Cytostim_5.fcs_16768" = "#F5B2AE", 
          "Cytostim_2.fcs_34414"  = "#CCE2C2", 
          "Cytostim_3.fcs_46733" = "#DDCAE3")

# Extended Data Figure 1O
pdf(paste0(plot.dir, Sys.Date(), "_clustering_methods_benchmark.pdf"), height = 7, width = 6)
p <- ggbarplot(f1_scores %>% mutate(method = factor(method, levels = levels)), x = "method", y = "f1", fill = "method", 
                add = c("mean_sd"),  error.plot = "errorbar", add.params = list(color = "black", shape = 21, size = 1)) + 
        scale_y_continuous(breaks = c(0.8, 0.9, 1.0), limits = c(NA, 1.0)) +
        coord_cartesian(ylim=c(0.7, 1.0)) + 
        xlab("Classification Method") + ylab("F1 Score") + 
        theme(axis.text.x = element_text(angle = 65, hjust = 1)) + 
        stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
                     geom = "errorbar", color = "black", width = 0.5, size = 0.5) + 
        geom_jitter(aes(fill = i), width = 0.2, height = 0, size = 1, shape = 21, color = "black") + 
        scale_fill_manual(values = cols) + theme(legend.position = "none")
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "ExtendedDataFigure1O_sourcedata.xlsx"))


#############################
##   DIFF FEATURE SPACES   ##
#############################
# IMPORTANT FEATURES + MARKERS -------

important_features <- c("FSC_A", "FSC_H", "FSC_ratio", "FSC_W", "LightLoss_Blue_A","Size_FSC", "Size_SSC", "SSC_W", "Comp_PE_Cy7_A") # from Fig1 C
features <- c(important_features,"Comp_FITC_A","Comp_APC_A","Comp_BV605_A","Comp_BV421_A") # add remaining fluoresence features 

clust <- data %>% 
  dplyr::select(all_of(features),status_fine, status_broad, element_name,uniq_id, otsu) %>% 
  dplyr::mutate_at(vars(features), ~ scale(.)) %>% # scaling and centering
  dplyr::mutate_at(vars(!features), ~ factor(.))

# run UMAP
clust <- run.umap(clust, use.cols = features, umap.seed = 42) # umap with k = 15
clust <- as.data.table(clust)


f1_scores_import <- data.frame() # empty results df for f1 scores
list <- list() # empty list for consensus classification 

# 100 iterations
pdf(paste0(res.dir, Sys.Date(), "_benchmarking_clustering_ratio_cluster_plots_important_features.pdf"), height = 6, width = 10, onefile = T)
for (i in 1:100) {
  message(paste0("iteration: ", i))
  
  # save random seed
  seed <- sample(1:1000000, size = 1)
  
  # louvain ----
  set.seed(seed)
  louvain_res <- clust %>% 
    dplyr::select(all_of(features)) %>% 
    bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = 1))) # shared nearest neighborhood graph by default 
  # add to df 
  clust[[paste0("louvain_", seed)]] <- factor(louvain_res)
  # add to list for consensus classification
  list[[i]] <- louvain_res

  # calculate ratio_cluster_plot 
  louvain_n <- paste0("louvain_", seed)
  
  res <- clust %>% 
    group_by(element_name, .data[[louvain_n]], otsu, .drop = F) %>% 
    count() %>% 
    group_by(.data[[louvain_n]], element_name, .drop = F) %>% 
    mutate(total = sum(n), freq = (n / sum(n)))
  
  # select clusters where more than 50 % of events are multiplets
  selected_clusters <- res %>% 
    dplyr::filter(otsu == "multiplet") %>% 
    group_by(.data[[louvain_n]]) %>% 
    summarize(mean = mean(freq)) %>% 
    dplyr::filter(mean > 0.5) %>%
    pull(.data[[louvain_n]])
  
  # add cluster type to df 
  clust <- clust %>% 
    mutate(!!sym(paste0("cluster_type_", louvain_n)) := if_else(.data[[louvain_n]] %in% selected_clusters, "multiplet", "singlet"))
  
  # plot ratio_cluster_plot and add selected multiplet clusters
  print(ggbarplot(res, x = louvain_n, y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
          scale_fill_manual(values = c("darkblue", "grey")) + 
          theme_classic() + 
          scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
          labs(title = louvain_n, subtitle = paste("multiplet clusters: ", paste(selected_clusters, collapse = ","))))
  
  # f1 score 
  confusion_mat <- clust %>% group_by(element_name, status_broad, .data[[paste0("cluster_type_", louvain_n)]]) %>% count()
  f1_import <- f1_score(confusion_mat, pred = paste0("cluster_type_", louvain_n))
  
  f1_scores_import <- rbind(f1_scores_import, f1_import)
  # break
}
dev.off()

f1_scores_import <- f1_scores_import %>% 
  mutate(features = "important_features+markers")

# save
saveRDS(f1_scores_import, paste0(data.dir, "f1_scores_import.rds"))


# consensus clustering solution
partition_list <- cl_ensemble(list = lapply(list, as.cl_hard_partition))
set.seed(42)
partition_consensus <- cl_consensus(partition_list)
class <- cl_class_ids(partition_consensus)

clust$consensus <- factor(class)

# plot clustering
ggplot(clust, aes(x = UMAP_X, y = -UMAP_Y, color = consensus)) + geom_point()

# Feature plots 
pdf(paste0(plot.dir, Sys.Date(), "_louvain_important_param_anno_featureplots.pdf"), width = 7, height = 8)
plots = vector('list', length(features))
for(i in seq_along(features)){
  plots[[i]] = ggplot(data = clust, aes(x = UMAP_X, y = -UMAP_Y, color = .data[[features[i]]])) + 
    geom_point() + theme_void() + theme(legend.position = "none") + 
    scale_colour_gradientn(colours = pals::parula(1000)) + 
    labs(title = features[i])
}
print(ggarrange(plotlist = plots))
dev.off()

# CD33-BV421
# CD19-BV605
# HLA.DR-BB515 (FITC)
# CD3-PE-Cy7
# CD45-APC

# anno
clust <- clust %>% 
  dplyr::mutate(celltype = case_when(
    consensus %in% c(2, 14, 11, 15) ~ "T cells",
    consensus  == 9 ~ "B cells",
    consensus  %in% c(6,7) ~ "CD33high Myeloid",
    consensus  == 8 ~ "CD33low Myeloid",
    consensus  %in% c(5,13,10) ~ "other CD45+",
    consensus  == 3 ~ "B*T",
    consensus  == 4 ~ "My*T",
    consensus  == 1 ~ "My*T*B",
    TRUE ~ "")) # Use an empty string as the default value

ggplot(clust, aes(x = UMAP_X, y = -UMAP_Y, color = celltype)) + geom_point()
ggplot(clust, aes(x = UMAP_X, y = -UMAP_Y, color = status_fine)) + geom_point()
ggplot(clust, aes(x = UMAP_X, y = -UMAP_Y, color = otsu)) + geom_point() + 
  scale_color_manual(values = c("#264896", "#D6D6D6"))

# manuscript figures 
# Extended Data Figure 1D
pdf(paste0(plot.dir, Sys.Date(), "_louvain_important_param_anno.pdf"), width = 9, height = 8)
print(ggplot(clust, aes(x = UMAP_X, y = -UMAP_Y, color = celltype)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = cols) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# extended data figure 1E lower panel 
pdf(paste0(plot.dir, Sys.Date(), "_louvain_important_param_ground_truth.pdf"), width = 9, height = 8)
print(ggplot(clust %>% mutate(status_fine = factor(status_fine, levels = c("quadruplet_more", "triplet", "doublet", "singlet"))), 
             aes(x = UMAP_X, y = -UMAP_Y, color = status_fine)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(status_fine)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = c("#4F80AF", "#A490BF", "#F2911A", "#D6D6D6")) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# extended data figure 1E upper panel 
pdf(paste0(plot.dir, Sys.Date(), "_louvain_important_param_otsu_threshold.pdf"), width = 9, height = 8)
print(ggplot(clust, aes(x = UMAP_X, y = -UMAP_Y, color = otsu)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(otsu)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = c("#264896", "#D6D6D6")) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# plot barplot with ground truth anno 
bar <- clust %>% 
  mutate(status_fine = factor(status_fine, levels = c("quadruplet_more", "triplet", "doublet", "singlet"))) %>%
  mutate(celltype = factor(celltype, levels = c("My*T*B", "My*T", "B*T", "CD33high Myeloid", "CD33low Myeloid", "T cells", "B cells", "other CD45+"))) %>%
  group_by(element_name, celltype, status_fine, .drop = F) %>% 
  count() %>%
  group_by(celltype, element_name, .drop = F) %>% 
  mutate(total = sum(n), freq = (n / sum(n)))

# extended data figure 1G
pdf(paste0(plot.dir, Sys.Date(), "_louvain_important_param_barplot_groundtruth.pdf"), width = 10, height = 6) 
p <- ggbarplot(bar, x = "celltype", y = "freq", fill = "status_fine", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = c("#4F80AF", "#A490BF", "#F2911A", "#D6D6D6")) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "ExtendedDataFigure1G_sourcedata.xlsx"))
  
# plot barplot with otsu, extended data figure 1F
otsu_bar <- clust %>% 
  mutate(celltype = factor(celltype, levels = c("My*T*B", "My*T", "B*T", "CD33high Myeloid", "CD33low Myeloid", "T cells", "B cells", "other CD45+"))) %>%
  group_by(element_name, celltype, otsu, .drop = F) %>% 
  count() %>%
  group_by(celltype, element_name, .drop = F) %>% 
  mutate(total = sum(n), freq = (n / sum(n)))

pdf(paste0(plot.dir, Sys.Date(), "_louvain_important_param_barplot_otsu_threshold.pdf"), width = 10, height = 6) 
p <- ggbarplot(otsu_bar, x = "celltype", y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = c("#264896", "#D6D6D6")) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "ExtendedDataFigure1F_sourcedata.xlsx"))

# save
saveRDS(clust, paste0(data.dir, "clust.rds"))

# MARKERS ONLY -----------------------------------------------------------------
markers <- c("Comp_PE_Cy7_A","Comp_FITC_A","Comp_APC_A","Comp_BV605_A","Comp_BV421_A")

clust_markers <- data %>% 
  dplyr::select(all_of(markers),status_fine, status_broad, element_name,uniq_id, otsu) %>% 
  dplyr::mutate_at(vars(markers), ~ scale(.)) %>% # scaling and centering
  dplyr::mutate_at(vars(!markers), ~ factor(.))

# run UMAP
clust_markers <- run.umap(clust_markers, use.cols = markers, umap.seed = 42) # umap with k = 15
clust_markers <- as.data.table(clust_markers)


f1_scores_markers <- data.frame() # empty results df for f1 scores
list_markers <- list() # empty list for consensus classification 

# 100 iterations
pdf(paste0(res.dir, Sys.Date(), "_benchmarking_clustering_ratio_cluster_plots_markers_only.pdf"), height = 6, width = 10, onefile = T)
for (i in 1:100) {
  message(paste0("iteration: ", i))
  
  # save random seed
  seed <- sample(1:1000000, size = 1)
  
  # louvain ----
  set.seed(seed)
  louvain_res <- clust_markers %>% 
    dplyr::select(all_of(markers)) %>% 
    bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = 1))) # shared nearest neighborhood graph by default 
  # add to df 
  clust_markers[[paste0("louvain_", seed)]] <- factor(louvain_res)
  # add to list for consensus classification
  list_markers[[i]] <- louvain_res
  
  # calculate ratio_cluster_plot 
  louvain_n <- paste0("louvain_", seed)
  
  res <- clust_markers %>% 
    group_by(element_name, .data[[louvain_n]], otsu, .drop = F) %>% 
    count() %>% 
    group_by(.data[[louvain_n]], element_name, .drop = F) %>% 
    mutate(total = sum(n), freq = (n / sum(n)))
  
  # select clusters where more than 50 % of events are multiplets
  selected_clusters <- res %>% 
    dplyr::filter(otsu == "multiplet") %>% 
    group_by(.data[[louvain_n]]) %>% 
    summarize(mean = mean(freq)) %>% 
    dplyr::filter(mean > 0.5) %>%
    pull(.data[[louvain_n]])
  
  # add cluster type to df 
  clust_markers <- clust_markers %>% 
    mutate(!!sym(paste0("cluster_type_", louvain_n)) := if_else(.data[[louvain_n]] %in% selected_clusters, "multiplet", "singlet"))
  
  # plot ratio_cluster_plot and add selected multiplet clusters
  print(ggbarplot(res, x = louvain_n, y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
          scale_fill_manual(values = c("darkblue", "grey")) + 
          theme_classic() + 
          scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) + 
          labs(title = louvain_n, subtitle = paste("multiplet clusters: ", paste(selected_clusters, collapse = ","))))
  
  # f1 score 
  confusion_mat <- clust_markers %>% group_by(element_name, status_broad, .data[[paste0("cluster_type_", louvain_n)]]) %>% count()
  f1_markers <- f1_score(confusion_mat, pred = paste0("cluster_type_", louvain_n))
  
  f1_scores_markers <- rbind(f1_scores_markers, f1_markers)
  # break
}
dev.off()

f1_scores_markers <- f1_scores_markers %>% 
  mutate(features = "markers")

saveRDS(f1_scores_markers, paste0(data.dir, "f1_scores_markers.rds"))

# consensus clustering solution
partition_list_markers <- cl_ensemble(list = lapply(list_markers, as.cl_hard_partition))
set.seed(42)
partition_consensus_markers <- cl_consensus(partition_list_markers)
class_markers <- cl_class_ids(partition_consensus_markers)

clust_markers$consensus <- factor(class_markers)

# plot clustering
ggplot(clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = consensus)) + geom_point()

# Feature plots 
pdf(paste0(plot.dir, Sys.Date(), "_louvain_markers_anno_featureplots.pdf"), width = 7, height = 4.5)
plots = vector('list', length(markers))
for(i in seq_along(markers)){
  plots[[i]] = ggplot(data = clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = .data[[markers[i]]])) + 
    geom_point() + theme_void() + theme(legend.position = "none") + 
    scale_colour_gradientn(colours = pals::parula(1000)) + 
    # scale_color_viridis() +
    labs(title = markers[i])
}
print(ggarrange(plotlist = plots))
dev.off()

# subcluster cluster 14
clust_markers$consensus <- as.character(clust_markers$consensus)
set.seed(42)
sub_louvain_res <- clust_markers %>%
  dplyr::filter(consensus == "14") %>% 
  dplyr::select(all_of(markers)) %>% 
  bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = 0.25)))  # shared nearest neighborhood graph by default 

sub_louvain_res <- as.character(sub_louvain_res)

# add to df
sub <- 0
clust_markers$sub.cluster <- NA
for (i in 1:nrow(clust_markers)) { 
  if (clust_markers$consensus[i] == "14") {
    sub <- sub + 1
    clust_markers$sub.cluster[i] <- str_c("14_", sub_louvain_res[sub])
  } else { 
    clust_markers$sub.cluster[i] <- clust_markers$consensus[i]
  }
}
ggplot(clust_markers, aes(x = UMAP_X, y = UMAP_Y, color = sub.cluster)) + geom_point()

# CD33-BV421
# CD19-BV605
# HLA.DR-BB515 (FITC)
# CD3-PE-Cy7
# CD45-APC

# annotation
clust_markers <- clust_markers %>% 
  dplyr::mutate(celltype = case_when(
    sub.cluster %in% c("4", "5", "8", "10", "16", "17") ~ "T cells",
    sub.cluster  == "9" ~ "B cells",
    sub.cluster  %in% c("1","7","15","13") ~ "CD33high Myeloid",
    sub.cluster  == "18" ~ "CD33low Myeloid",
    sub.cluster  %in% c("6","20","19","2") ~ "other CD45+",
    sub.cluster  == "3" ~ "B*T",
    sub.cluster  %in% c("14_1", "14_2") ~ "My*T",
    sub.cluster  == "14_3" ~ "My*T*B",
    TRUE ~ "")) # Use an empty string as the default value

ggplot(clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = celltype)) + geom_point()
ggplot(clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = status_fine)) + geom_point()
ggplot(clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = otsu)) + geom_point() + 
  scale_color_manual(values = c("#264896", "#D6D6D6"))

# manuscript figures 
# extended data figure 1H 
pdf(paste0(plot.dir, Sys.Date(), "_louvain_markers_anno.pdf"), width = 9, height = 8)
print(ggplot(clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = celltype)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = cols) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# extended data figure 1I lower panel
pdf(paste0(plot.dir, Sys.Date(), "_louvain_markers_ground_truth.pdf"), width = 9, height = 8)
print(ggplot(clust_markers %>% mutate(status_fine = factor(status_fine, levels = c("quadruplet_more", "triplet", "doublet", "singlet"))), 
             aes(x = UMAP_X, y = -UMAP_Y, color = status_fine)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(status_fine)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = c("#4F80AF", "#A490BF", "#F2911A", "#D6D6D6")) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# extended data figure 1I upper panel
pdf(paste0(plot.dir, Sys.Date(), "_louvain_markers_otsu_threshold.pdf"), width = 9, height = 8)
print(ggplot(clust_markers, aes(x = UMAP_X, y = -UMAP_Y, color = otsu)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(otsu)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = c("#264896", "#D6D6D6")) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# plot bar plot with ground truth annotation 
bar_markers <- clust_markers %>% 
  mutate(status_fine = factor(status_fine, levels = c("quadruplet_more", "triplet", "doublet", "singlet"))) %>%
  mutate(celltype = factor(celltype, levels = c("My*T*B", "My*T", "B*T", "CD33high Myeloid", "CD33low Myeloid", "T cells", "B cells", "other CD45+"))) %>%
  group_by(element_name, celltype, status_fine, .drop = F) %>% 
  count() %>%
  group_by(celltype, element_name, .drop = F) %>% 
  mutate(total = sum(n), freq = (n / sum(n)))

# extended data figure 1K
pdf(paste0(plot.dir, Sys.Date(), "_louvain_markers_barplot_groundtruth.pdf"), width = 10, height = 6) 
p <- ggbarplot(bar_markers, x = "celltype", y = "freq", fill = "status_fine", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = c("#4F80AF", "#A490BF", "#F2911A", "#D6D6D6")) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "ExtendedDataFigure1K_sourcedata.xlsx"))

# plot barplot with otsu 
otsu_bar_markers <- clust_markers %>% 
  mutate(celltype = factor(celltype, levels = c("My*T*B", "My*T", "B*T", "CD33high Myeloid", "CD33low Myeloid", "T cells", "B cells", "other CD45+"))) %>%
  group_by(element_name, celltype, otsu, .drop = F) %>% 
  count() %>%
  group_by(celltype, element_name, .drop = F) %>% 
  mutate(total = sum(n), freq = (n / sum(n)))

# extended data figure 1J
pdf(paste0(plot.dir, Sys.Date(), "_louvain_markers_barplot_otsu_threshold.pdf"), width = 10, height = 6) 
p <- ggbarplot(otsu_bar_markers, x = "celltype", y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = c("#264896", "#D6D6D6")) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "ExtendedDataFigure1J_sourcedata.xlsx"))

#save
saveRDS(clust_markers, paste0(data.dir, "clust_markers.rds"))


# FLOW PARAMETERS -----------------------------------------------------------------
f1_scores_flow <- f1_scores %>% 
  dplyr::filter(method == "louvain") %>% 
  mutate(features = "flow+scatter+ratio")

# save 
saveRDS(f1_scores_flow, paste0(data.dir, "f1_scores_flow.rds"))

# extract all clustering solutions for louvain based on flow markers 
list_flow <- list()
i <- 0
for (col in colnames(clust_bm)) { 
  if (str_detect(col, "^louvain")) {
    i <- i + 1
    list_flow[[i]] <- clust_bm %>% pull(col)
  }
}

# consensus clustering 
partition_list_flow <- cl_ensemble(list = lapply(list_flow, as.cl_hard_partition))
set.seed(42)
partition_consensus_flow <- cl_consensus(partition_list_flow)
class_flow <- cl_class_ids(partition_consensus_flow)

clust_bm$consensus_louvain <- as.character((class_flow))

# plot consensus clustering
ggplot(clust_bm, aes(x = UMAP_X, y = UMAP_Y, color = consensus_louvain)) + geom_point()

# Feature plots 
pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_anno_featureplots.pdf"), width = 7, height = 5)
plots = vector('list', length(flow_param))
for(i in seq_along(flow_param)){
  plots[[i]] = ggplot(data = clust_bm, aes(x = UMAP_X, y = UMAP_Y, color = .data[[flow_param[i]]])) + 
    geom_point() + theme_void() + theme(legend.position = "none") + 
    scale_colour_gradientn(colours = pals::parula(1000)) + 
    # scale_colour_viridis() + 
    labs(title = flow_param[i])
}
print(ggarrange(plotlist = plots))
dev.off()

# subcluster cluster 12
set.seed(42)
sub_louvain_res <- clust_bm %>%
  dplyr::filter(consensus_louvain == "12") %>% 
  dplyr::select(all_of(flow_param)) %>% 
  bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = 0.25)))  # shared nearest neighborhood graph by default 
  
sub_louvain_res <- as.character(sub_louvain_res)

# add to df
sub <- 0
clust_bm$sub.cluster <- NA
for (i in 1:nrow(clust_bm)) { 
  if (clust_bm$consensus_louvain[i] == "12") {
    sub <- sub + 1
    clust_bm$sub.cluster[i] <- str_c("12_", sub_louvain_res[sub])
  } else { 
    clust_bm$sub.cluster[i] <- clust_bm$consensus_louvain[i]
  }
}

# CD33-BV421
# CD19-BV605
# HLA.DR-BB515 (FITC)
# CD3-PE-Cy7
# CD45-APC

# annotation
clust_bm <- clust_bm %>% 
  dplyr::mutate(celltype = case_when(
    sub.cluster %in% c("9", "11", "14", "15") ~ "T cells",
    sub.cluster  == "2" ~ "B cells",
    sub.cluster  %in% c("7", "1", "5") ~ "CD33high Myeloid",
    sub.cluster  == "8" ~ "CD33low Myeloid",
    sub.cluster  %in% c("3", "13", "6") ~ "other CD45+",
    sub.cluster  == "10" ~ "B*T",
    sub.cluster  == "12_1" ~ "My*T",
    sub.cluster  == "12_2" ~ "My*T*B",
    TRUE ~ "")) # Use an empty string as the default value

# manuscript figures 
# Figure 1F
pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_anno.pdf"), width = 9, height = 8)
print(ggplot(clust_bm, aes(x = UMAP_X, y = UMAP_Y, color = celltype)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = cols) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

louvain_cols <- c(pals::tol.groundcover())
names(louvain_cols) <- NULL
pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_clusters.pdf"), width = 9, height = 8)
print(ggplot(clust_bm, aes(x = UMAP_X, y = UMAP_Y, color = consensus_louvain)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(consensus_louvain)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = louvain_cols) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

pdf(paste0(plot.dir, Sys.Date(), "_ratio_cluster_plot_flow_louvain.pdf"), height = 7, width = 20)
print(clust_bm %>%
        group_by(consensus_louvain, otsu) %>%
        count(otsu) %>%
        mutate(consensus_louvain = factor(consensus_louvain, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))) %>% 
        ggplot(aes(x = consensus_louvain, y = n, fill = otsu)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = rev(c("#D6D6D6", "#264896"))) +
        xlab("Clusters") +
        labs(y = "Count") +
        ggtitle("Cluster Analysis") +
        theme_classic())
dev.off()

# selected dbts clusters
dist <- clust_bm %>% group_by(consensus_louvain, otsu) %>% 
  count(otsu) %>% ungroup() %>% group_by(consensus_louvain) %>% 
  mutate(ratio_ratio = n/sum(n)) %>% dplyr::filter(otsu == 
                                                     "multiplet")
cutoff <- quantile(dist$ratio_ratio, probs = 0.8)
cluster_to_use <- as.vector(pull(dist[dist$ratio_ratio > 
                                        cutoff, ], .data[["consensus_louvain"]]))

# figure 1F right, lower
pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_ground_truth.pdf"), width = 9, height = 8)
print(ggplot(clust_bm %>% mutate(status_fine = factor(status_fine, levels = c("quadruplet_more", "triplet", "doublet", "singlet"))), 
             aes(x = UMAP_X, y = UMAP_Y, color = status_fine)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(status_fine)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = c("#4F80AF", "#A490BF", "#F2911A", "#D6D6D6")) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# figure 1F right, upper
pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_otsu_threshold.pdf"), width = 9, height = 8)
print(ggplot(clust_bm, aes(x = UMAP_X, y = UMAP_Y, color = otsu)) + 
        ggrastr::geom_point_rast(size = 3.5, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(otsu)), size = 2, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = c("#264896", "#D6D6D6")) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# plot barplot with ground truth anno 
bar_flow <- clust_bm %>% 
  mutate(status_fine = factor(status_fine, levels = c("quadruplet_more", "triplet", "doublet", "singlet"))) %>%
  mutate(celltype = factor(celltype, levels = c("My*T*B", "My*T", "B*T", "CD33high Myeloid", "CD33low Myeloid", "T cells", "B cells", "other CD45+"))) %>%
  group_by(element_name, celltype, status_fine, .drop = F) %>% 
  count() %>%
  group_by(celltype, element_name, .drop = F) %>% 
  mutate(total = sum(n), freq = (n / sum(n)))

# figure 1H
pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_barplot_groundtruth.pdf"), width = 10, height = 6)
p <- ggbarplot(bar_flow, x = "celltype", y = "freq", fill = "status_fine", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = c("#4F80AF", "#A490BF", "#F2911A", "#D6D6D6")) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "Figure1H_sourcedata.xlsx"))

# plot barplot with otsu, figure 1G
otsu_bar_flow <- clust_bm %>% 
  mutate(celltype = factor(celltype, levels = c("My*T*B", "My*T", "B*T", "CD33high Myeloid", "CD33low Myeloid", "T cells", "B cells", "other CD45+"))) %>%
  group_by(element_name, celltype, otsu, .drop = F) %>% 
  count() %>%
  group_by(celltype, element_name, .drop = F) %>% 
  mutate(total = sum(n), freq = (n / sum(n)))

pdf(paste0(plot.dir, Sys.Date(), "_louvain_flow_param_barplot_otsu_threshold.pdf"), width = 10, height = 6)
p <- ggbarplot(otsu_bar_flow, x = "celltype", y = "freq", fill = "otsu", add = "mean_sd", error.plot = "errorbar") + 
  scale_fill_manual(values = c("#264896", "#D6D6D6")) + 
  theme_classic() + theme(legend.position = "none") + 
  scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
print(p)
dev.off()

# export source data 
write_xlsx(p$data, paste0(plot.dir, "Figure1G_sourcedata.xlsx"))

# heatmap 
cols_replicate <- c("Cytostim_1.fcs_6910" = "#B2CCE2", 
                    "Cytostim_5.fcs_16768" = "#F5B2AE", 
                    "Cytostim_2.fcs_34414"  = "#CCE2C2", 
                    "Cytostim_3.fcs_46733" = "#DDCAE3")

hm <- clust_bm %>%
  group_by(element_name, celltype) %>% 
  summarise(across(flow_param, ~ mean(.))) %>% # mean expression per celltype
  ungroup() %>%
  dplyr::select(all_of(flow_param), celltype, element_name)

# row annotation
rowAnn <- HeatmapAnnotation(df = dplyr::select(hm, c(element_name, celltype)),
                            which = 'row',
                            col = list(celltype=cols,
                                       element_name=cols_replicate))

pdf(paste0(plot.dir, Sys.Date(), "_heatmap_flow_parameters.pdf"), width = 12, height = 14)
p <- (Heatmap(as.matrix(hm %>% dplyr::select(all_of(flow_param))),
        col=pals::coolwarm(100),
        name = "z-score",
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D2",
        column_title = "", 
        row_title = "",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 9),show_row_names = F,
        cluster_rows = T,
        cluster_columns = T,
        right_annotation = rowAnn,
        heatmap_legend_param = list(legend_gp = gpar(fontsize = 7)),
        show_parent_dend_line = FALSE))
print(p)      
dev.off()

# export source data 
dat <- cbind(as.data.frame(p@matrix),  as.data.frame(rowAnn@anno_list$element_name@fun@var_env$value), as.data.frame(rowAnn@anno_list$celltype@fun@var_env$value))
write_xlsx(dat, paste0(plot.dir, "ExtendedDataFigure1M_sourcedata.xlsx"))

#save
saveRDS(clust_bm, paste0(data.dir, "clust_bm.rds"))

# plot otsu colored by ground truth, nudged, split by cluster anno -> extended data figure 1N
pdf(paste0(plot.dir, Sys.Date(), "_otsu_threshold_histogram_groundtruth_nudged_split_clustertype.pdf"), height = 5, width = 6)
p <- ggplot(left_join(data, clust_bm %>% mutate(consensus_cluster_type = if_else(celltype %in% c("My*T*B", "My*T", "B*T"), "multiplet_cluster", "singlet_cluster")) %>% dplyr::select(uniq_id, consensus_cluster_type), by = "uniq_id") %>%
         mutate(status_fine = factor(status_fine, levels = c("singlet", "doublet", "triplet", "quadruplet_more"))),
       aes(x = FSC_ratio, fill = status_fine)) +
  geom_histogram(bins = 1000, position = "nudge") +
  geom_vline(xintercept = otsu_threshold, color = "black", linetype = "dashed") +
  theme_classic() + facet_wrap(~consensus_cluster_type ) +
  scale_fill_manual(values = c("#D6D6D6", "#F2911A", "#A490BF", "#4F80AF")) +
  scale_y_continuous(expand = c(0,0))
print(p)
dev.off()

# export source data
write_xlsx(p$data, paste0(plot.dir, "ExtendedDataFigure1N_sourcedata.xlsx"))

# compare ----
f1_scores_otsu <- f1_score(otsu_res, "otsu") %>% 
  mutate(features = "otsu")
all_f1_scores <- rbind(f1_scores_flow %>% dplyr::select(-method), f1_scores_markers, f1_scores_import, f1_scores_otsu)

# save 
saveRDS(all_f1_scores, paste0(data.dir, "all_f1_scores_clustering.rds"))

# plot
cols <- c("otsu" = "#768EC3", 
          "markers" = "#768EC3", 
          "flow+scatter+ratio" = "#768EC3", 
          "important_features+markers" = "#768EC3", 
          "Cytostim_1.fcs_6910" = "#B2CCE2", 
          "Cytostim_5.fcs_16768" = "#F5B2AE", 
          "Cytostim_2.fcs_34414"  = "#CCE2C2", 
          "Cytostim_3.fcs_46733" = "#DDCAE3")

levels <- c("otsu", "markers", "flow+scatter+ratio", "important_features+markers")
pdf(paste0(plot.dir, Sys.Date(), "_clustering_methods_benchmark.pdf"), height = 7, width = 6)
p <- ggbarplot(all_f1_scores %>% mutate(features = factor(features, levels = levels)), x = "features", y = "f1", fill = "features", 
          add = c("mean_sd"),  error.plot = "errorbar", add.params = list(color = "black", shape = 21, size = 1)) + 
  scale_y_continuous(breaks = c(0.8, 0.9, 1.0), limits = c(NA, 1.0)) +
  coord_cartesian(ylim=c(0.7, 1.0)) + 
  xlab("Classification Method") + ylab("F1 Score") + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", color = "black", width = 0.5, size = 0.5) + 
  geom_jitter(aes(fill = i), width = 0.2, height = 0, size = 1, shape = 21,  color = "black") + 
  scale_fill_manual(values = cols)
print(p)
dev.off()

# export source data
write_xlsx(p$data, paste0(plot.dir, "Figure1E_sourcedata.xlsx"))

pdf(paste0(plot.dir, Sys.Date(), "_clustering_methods_benchmark_boxplots.pdf"), height = 4.5, width = 6)
print(ggplot(all_f1_scores %>% mutate(features = factor(features, levels = levels)), aes(x = features, y = f1)) + 
        geom_boxplot(aes(fill = features)) + 
        geom_jitter(aes(color = i), alpha = 0.2) + 
        scale_y_continuous(breaks = c(0.8, 0.9, 1.0), limits = c(0.70, NA)) + 
        scale_fill_manual(values = cols) + 
        scale_color_manual(values = cols) + 
        theme_classic()
      )
dev.off()

########################
##                    ##
##   ADJ RAND INDEX   ##
##                    ##
##   Lea Jopp-Saile   ##
##                    ##
########################

# repeat adjusted rand index, but on consensus clustering solutions
image <- clust %>% dplyr::select(all_of(features), status_fine, status_broad, element_name, uniq_id, UMAP_X, UMAP_Y)
conventional <- clust_bm %>% dplyr::select(all_of(flow_param), status_fine, status_broad, element_name, uniq_id, UMAP_X, UMAP_Y)

set.seed(42)
list_image <- list()
list_conv <- list()

# resolutions 0 to 2
for (res in seq(0,2,by=0.1)) {
  message(paste0("resolution: ", res))
  
  # image parameter ----
  for (i in 1:100) { # 100 iterations per resolution
    louvain_res <- image %>% 
      dplyr::select(all_of(features)) %>% 
      bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = res))) # shared nearest neighborhood graph by default 
    
    # add to list for consensus classification
    list_image[[i]] <- louvain_res
  }
  
  # consensus clustering 
  partition_list <- cl_ensemble(list = lapply(list_image, as.cl_hard_partition))
  partition_consensus <- cl_consensus(partition_list)
  class <- cl_class_ids(partition_consensus)
  
  # add to df
  image[[paste0("consensus_", res)]] <- as.character((class))
  
  # conventional flow parameter ----
  for (i in 1:100) { # 100 iterations per resolution
    louvain_res <- conventional %>% 
      dplyr::select(all_of(flow_param)) %>% 
      bluster::clusterRows(NNGraphParam(k=5, cluster.fun = "louvain", cluster.args = list(resolution = res))) 
    
    # add to list for consensus classification
    list_conv[[i]] <- louvain_res
  }
  
  # consensus clustering 
  partition_list <- cl_ensemble(list = lapply(list_conv, as.cl_hard_partition))
  partition_consensus <- cl_consensus(partition_list)
  class <- cl_class_ids(partition_consensus)
  
  # add to df
  conventional[[paste0("consensus_", res)]] <- as.character((class))
}

# join consensus solutions across resolutions
df <- left_join(conventional %>% dplyr::select(contains("consensus")) %>% rownames_to_column(), 
                image %>% dplyr::select(contains("consensus")) %>% rownames_to_column(), 
                by = "rowname")

# calculate rand index for each resolution
rand_index <- unlist(lapply(grep("consensus", colnames(image), value = T), function(i) {
  rand <- mclust::adjustedRandIndex(df[[paste0(i,".x")]], df[[paste0(i,".y")]])
}
))


rand_df <- data.frame(x = seq(0,2, by = 0.1), y = rand_index)

saveRDS(rand_df, paste0(data.dir, "rand_df.rds"))

# plot
pdf(paste0(plot.dir, Sys.Date(), "_adjusted_rand_index_consensus.pdf"), height = 4.3, width = 5)
p <- ggplot(rand_df, aes(x,y)) + 
  geom_line(aes(colour = "adjusted_rand_index"), size = 2) +
  ylim(0,1) +
  scale_color_manual(values = "#A2195B") +
  geom_hline(yintercept = mean(rand_index), col = "black") +
  xlab("Cluster resolution") + ylab("Adjusted Rand Index") +
  theme_classic()
print(p)
dev.off()

# export source data
write_xlsx(p$data, paste0(plot.dir, "Figure1I_sourcedata.xlsx"))


# HEATMAP ----------------------------------------------------------------------
ari_df <- lapply(seq(2,22), function(x){
  colnames(df)[x]
  unlist(lapply(seq(23,43), function(y){
    colnames(df)[y]
    rand <-  mclust::adjustedRandIndex(df[[x]], df[[y]])
  }))
})

ari_df <- do.call(cbind,ari_df)

colnames(ari_df) <- paste0("res_", seq(0,2, by=0.1))
rownames(ari_df) <- paste0("res_", seq(0,2, by=0.1))

library(pheatmap)
pdf(paste0(plot.dir, Sys.Date(), "_adjusted_rand_index_consensus_heatmap.pdf"), height = 6.4, width = 7)
pheatmap(ari_df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#352A87", "#33B7A0", "#F9FB0E"))(100),
         annotation_legend = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Heatmap of Adjusted Rand Index")
dev.off()
#"ARI Interpretation: <0 = Disagreement, 0 = Random, >0 to <0.10 = Poor, 0.10 to <0.40 = Fair, 0.40 to <0.75 = Good, >0.75 = Excellent"
