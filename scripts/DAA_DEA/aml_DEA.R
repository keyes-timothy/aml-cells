#aml_DEA.R
#Description: performs differential expression and differential abundance analysis on aml data
#input: 
#output: 
#side-effects: 
#optimizations: 

#libraries
library(edgeR)
library(samr)
library(multcomp)
library(lme4)
library(gridExtra)
library(viridis)
library(Rtsne)
library(Rtsne.multicore)
library(umap)
library(tidyverse)
library(rlang)

#set up output directories
daa.output <- file.path(OUTPUT.DIRECTORY, "daa_results")
dea.output <- file.path(OUTPUT.DIRECTORY, "dea_results")


#read in cluster summary data
load("/Users/tkeyes/Desktop/scg results/AML_cluster_summaries.RData") #change to Box location

#read in expression.matrix
load(file.path(OUTPUT.DIRECTORY, "data", "AML_matrix_clustered.RData"))


################################           Diagnostic Vs. Relapse Samples           #######################################


paired.samples <- c("PARBIU", "PASGWH", "PARWXU","PARZUU", 
                    "PASVZC", "PASXVC", "PATGIY", "PATHIW", 
                    "PATJSP", "PATLHB")

#####   With literally everything
#form master matrices - m x n matrix with m patient/condition pairs and n features 
#(aka channel centroids within each cluster)
form_dxrx_master <- function(abundances = NULL, centroids = NULL, derived = NULL, stims = NULL, cluster_method = NULL){ 
  master.matrix <- left_join(x = purrr::pluck(abundances, cluster_method ), 
                             y = purrr::pluck(centroids, cluster_method )) %>% 
    left_join(y = purrr::pluck(derived, cluster_method)) %>% 
    left_join(y = purrr::pluck(stims, cluster_method )) %>% 
    dplyr::select(-num.cells, -total.cells) %>% 
    dplyr::filter(condition != "Healthy") %>% 
    dplyr::filter(patient %in% paired.samples) %>% 
    gather(key = "variable", value = "value", -patient, -condition, -!!as.name(cluster_method)) %>% ungroup()
  master.matrix$combined.var <- paste(master.matrix$variable, master.matrix[[cluster_method]], sep = "_")
  master.matrix <- master.matrix %>% dplyr::select(patient, combined.var, value, condition) %>%
    tidyr::spread(key = combined.var, value = value)
  return(master.matrix)
}

#save all feature matrices for all clustering methods that I have 
dx.rx.matrices <- list()

for(i in 1:length(abundance.by.pat.con)){ #for each clustering method, form a master_matrix and store in the dx.rx.matrices list
  cluster_method <- names(abundance.by.pat.con)[[i]]
  new.matrix <- form_dxrx_master(abundances = abundance.by.pat.con, 
                                 centroids = centroids.by.pat.con, 
                                 derived = derived.metrics.by.pat.con, 
                                 stims = stims.by.pat.con, 
                                 cluster_method = cluster_method)
  new.matrix$patient[new.matrix$patient == "PASTRP"] <- "PASRTP"
  dx.rx.matrices[[cluster_method]] <- new.matrix
}


#read in clinical metadata
clinical.matrix <- read_csv(file.path(OUTPUT.DIRECTORY, "..", "AML_metadata.csv"))
clinical.matrix$patient <- clinical.matrix$Patient
clinical.matrix$Patient <- NULL

sam.matrices <- map(dx.rx.matrices, function(data) inner_join(x = data, y = clinical.matrix))
sam.matrices$phenograph.metacluster <- NULL

#apply samr function to Dx-Rx pairs
sam_results <- list()

aml_sam <- function(sam.matrix = NULL){ 
  data <- sam.matrix %>% dplyr::select(contains("_"))
  feature.names <- colnames(data)
  response <- c(1:nrow(sam.matrix), 1:nrow(sam.matrix)) %>% sort() %>% 
    purrr::modify_at(.at = seq(from = 1, to = length(.), by = 2), .f = function(x) -x)
  result <- SAM(x = t(as.matrix(data)), y = response, 
                resp.type = "Two class paired", nperms = 5000, 
                geneid = feature.names, genenames = feature.names)
}

SAMs <- map(.x = sam.matrices, .f = aml_sam)


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


###### With only abundances
form_dxrx_master <- function(abundances = NULL, centroids = NULL, derived = NULL, stims = NULL, cluster_method = NULL){ 
  master.matrix <- purrr::pluck(abundances, cluster_method) %>%
    dplyr::select(-num.cells, -total.cells) %>% 
    dplyr::filter(condition != "Healthy") %>% 
    dplyr::filter(patient %in% paired.samples) %>% 
    gather(key = "variable", value = "value", -patient, -condition, -!!as.name(cluster_method)) %>% ungroup()
  master.matrix$combined.var <- paste(master.matrix$variable, master.matrix[[cluster_method]], sep = "_")
  master.matrix <- master.matrix %>% dplyr::select(patient, combined.var, value, condition) %>%
    tidyr::spread(key = combined.var, value = value, fill = 0)
  return(master.matrix)
}

dx.rx.matrices <- list()

for(i in 1:length(abundance.by.pat.con)){ 
  cluster_method <- names(abundance.by.pat.con)[[i]]
  new.matrix <- form_dxrx_master(abundances = abundance.by.pat.con, 
                                 centroids = tibble(), 
                                 derived = derived.metrics.by.pat.con, 
                                 stims = tibble(), 
                                 cluster_method = cluster_method)
  new.matrix$patient[new.matrix$patient == "PASTRP"] <- "PASRTP"
  dx.rx.matrices[[cluster_method]] <- new.matrix
}

#combine with clinical data
sam.matrices <- map(dx.rx.matrices, function(data) inner_join(x = data, y = clinical.matrix))
#sam.matrices$phenograph.metacluster <- NULL

sam_results <- list()

aml_sam <- function(sam.matrix = NULL){ 
  data <- sam.matrix %>% dplyr::select(contains("_"))
  feature.names <- colnames(data)
  response <- c(1:nrow(sam.matrix), 1:nrow(sam.matrix)) %>% sort() %>% 
     purrr::modify_at(.at = seq(from = 1, to = length(.), by = 2), .f = function(x) -x)
  result <- SAM(x = t(as.matrix(data)), y = response, 
                resp.type = "Two class paired", nperms = 5000, 
                geneid = feature.names, genenames = feature.names)
}

SAMs <- map(.x = sam.matrices, .f = aml_sam)

#save a character vector of cluster numbers that are expanded in the Rx samples
phenograph.expanded.clusters <- SAMs$phenograph.metacluster$siggenes.table$genes.up %>% as_tibble() %>% 
  dplyr::filter(`q-value(%)` < 0.05) %>% dplyr::select(`Gene ID`) %>% 
  purrr::pluck(1) %>% str_split_fixed(pattern = "_", n = 2) %>% as_tibble() %>% 
  purrr::pluck(2)

#make heatmap for overall mah-cluster abundance
heatmap.matrix <- sam.matrices$Mah.cluster %>% dplyr::select(contains("cell.percentage")) %>% as.matrix()
row.names(heatmap.matrix) <- 1:nrow(heatmap.matrix)
row.annotations <- data.frame(patient = sam.matrices$Mah.cluster$patient,
                              condition = sam.matrices$Mah.cluster$condition)


my.colnames <- colnames(heatmap.matrix) %>% str_split_fixed(pattern = "_", n = 2) %>% 
  as_tibble() %>% pluck(2) 
colnames(heatmap.matrix) <- my.colnames

my.order <- CLASSIFIER.POPULATIONS
heatmap.matrix <- heatmap.matrix[,my.order]

annotation_colors <- list(condition = c(Dx = "skyblue", Rx = "red"))
my.heatmap <- pheatmap(mat = heatmap.matrix, 
                       annotation_row = row.annotations, 
                       cluster_cols = FALSE, cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       annotation_colors = annotation_colors,
                       gaps_row = c(seq(from = 2, to = nrow(heatmap.matrix), 
                                        by = 2)), 
                       cell_width = 7, cell_height = 5)

ggsave(filename = "Mah_cluster_overall.pdf", 
       plot = my.heatmap, device = "pdf", path = daa.output, width = 6, 
       height = 5)

#make heatmap for overall phenograph cluster abundance
heatmap.matrix <- sam.matrices$phenograph.metacluster %>% dplyr::select(contains("cell.percentage"))  %>% as.matrix()
row.names(heatmap.matrix) <- 1:nrow(heatmap.matrix)
row.annotations <- data.frame(patient = sam.matrices$phenograph.metacluster$patient,
                              condition = sam.matrices$phenograph.metacluster$condition)
  

my.colnames <- colnames(heatmap.matrix) %>% str_split_fixed(pattern = "_", n = 2) %>% 
  as_tibble() %>% pluck(2) 
colnames(heatmap.matrix) <- my.colnames

my.order <- as.integer(my.colnames) %>% sort.int() %>% as.character()
heatmap.matrix <- heatmap.matrix[,my.order]

annotation_colors <- list(condition = c(Dx = "skyblue", Rx = "red"))
my.heatmap <- pheatmap(mat = heatmap.matrix, 
                       annotation_row = row.annotations, 
                       cluster_cols = FALSE, cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       annotation_colors = annotation_colors,
                       gaps_row = c(seq(from = 2, to = nrow(heatmap.matrix), 
                                        by = 2)), 
                       cell_width = 7, cell_height = 5)

ggsave(filename = "phenoclusters_overall.pdf", 
       plot = my.heatmap, device = "pdf", path = daa.output, width = 10, 
       height = 5)

my.maxima <- map_dbl(1:nrow(heatmap.matrix), function(x) max(heatmap.matrix[x,]))
my.minima <- map_dbl(1:nrow(heatmap.matrix), function(x) min(heatmap.matrix[x,]))
my.num.clusters <- map_int(1:nrow(heatmap.matrix), function(x) sum(heatmap.matrix[x,] > 0))
average.num.clusters <- mean(my.num.clusters)

#make heatmap for expanded phenograph cluster abundance
heatmap.matrix <- heatmap.matrix %>% as_tibble() %>% 
  dplyr::select(which(colnames(heatmap.matrix) %in% phenograph.expanded.clusters)) %>% as.matrix()
row.names(heatmap.matrix) <- 1:nrow(heatmap.matrix)

my.heatmap <- pheatmap(mat = (heatmap.matrix*100), 
                       annotation_row = row.annotations, 
                       cluster_cols = FALSE, cluster_rows = FALSE,
                       show_rownames = FALSE, 
                       annotation_colors = annotation_colors,
                       gaps_row = c(seq(from = 2, to = nrow(heatmap.matrix), 
                                        by = 2)), 
                       cell_width = 3, cell_height = 3, 
                       angle_col = 0, display_numbers = FALSE)

ggsave(filename = "phenoclusters_expanded.pdf", 
       plot = my.heatmap, device = "pdf", path = daa.output, width = 6, 
       height = 5)

my.maxima <- map_dbl(1:nrow(heatmap.matrix), function(x) max(heatmap.matrix[x,]))
my.minima <- map_dbl(1:nrow(heatmap.matrix), function(x) min(heatmap.matrix[x,]))
my.num.clusters <- map_int(1:nrow(heatmap.matrix), function(x) sum(heatmap.matrix[x,] > 0))
average.num.clusters <- mean(my.num.clusters)


#make a faceted boxplot (with points representing each patient) of all metacluster abundances across patients 


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


#####  With only signaling - have to play with this because of all the missing values...not super functional atm
form_dxrx_master <- function(abundances = NULL, centroids = NULL, derived = NULL, stims = NULL, cluster_method = NULL){ 
  master.matrix <- purrr::pluck(stims, cluster_method) %>%
    dplyr::filter(condition != "Healthy") %>% 
    dplyr::filter(patient %in% paired.samples) %>% 
    gather(key = "variable", value = "value", -patient, -condition, -!!as.name(cluster_method)) %>% ungroup()
  master.matrix$combined.var <- paste(master.matrix$variable, master.matrix[[cluster_method]], sep = "_")
  master.matrix <- master.matrix %>% dplyr::select(patient, combined.var, value, condition) %>%
    tidyr::spread(key = combined.var, value = value)
  return(master.matrix)
}

dx.rx.matrices <- list()

for(i in 1:length(abundance.by.pat.con)){ 
  cluster_method <- names(abundance.by.pat.con)[[i]]
  new.matrix <- form_dxrx_master(abundances = abundance.by.pat.con, 
                                 centroids = centroids.by.pat.con, 
                                 derived = derived.metrics.by.pat.con, 
                                 stims = stims.by.pat.con, 
                                 cluster_method = cluster_method)
  new.matrix$patient[new.matrix$patient == "PASTRP"] <- "PASRTP"
  dx.rx.matrices[[cluster_method]] <- new.matrix
}

#combine with clinical data
sam.matrices <- map(dx.rx.matrices, function(data) inner_join(x = data, y = clinical.matrix))
pheno.sam.matrix <- sam.matrices$phenograph.metacluster 
sam.matrices$phenograph.metacluster <- NULL

sam_results <- list()

aml_sam <- function(sam.matrix = NULL){ 
  data <- sam.matrix %>% dplyr::select(contains("_"))
  feature.names <- colnames(data)
  response <- c(1:nrow(sam.matrix), 1:nrow(sam.matrix)) %>% sort() %>% 
    purrr::modify_at(.at = seq(from = 1, to = length(.), by = 2), .f = function(x) -x)
  result <- SAM(x = t(as.matrix(data)), y = response, 
                resp.type = "Two class paired", nperms = 5000, 
                geneid = feature.names, genenames = feature.names)
}

SAMs <- map(.x = sam.matrices, .f = aml_sam)

#visualize difference in signaling in Mah.clusters
dim_red_matrix <- expression.matrix %>% dplyr::filter(patient %in% paired.samples) %>% 
  group_by(Mah.cluster) %>% 
  sample_n(size = 1000, replace = F) %>% ungroup()

tSNE <- Rtsne(X = dim_red_matrix[1:41], dims = 2)
tSNE.data <- tSNE$Y %>% as_tibble()
tSNE.data <- cbind(tSNE.data, dim_red_matrix)
names(tSNE.data)[1:2] <- paste0("tSNE", 1:2)

tSNE.data <- tSNE.data %>% 
  mutate_at(.vars = vars(contains(match = "cluster")), .funs = as.factor)

dir.create(file.path(dea.output, "Mah.cluster tSNE plots"))
for(i in 1:length(colnames(tSNE.data))){ 
  channel.name <- colnames(tSNE.data)[i]
  channel <- sym(colnames(tSNE.data)[i])
  tSNE.plot <- ggplot(data = tSNE.data,
                      mapping = aes(x = tSNE1, y = tSNE2, color = !!channel)) + 
    geom_point() + theme_bw() 
  if(is.numeric(tSNE.data[,channel.name])){
    tSNE.plot <- tSNE.plot + 
      scale_color_viridis(limits = c(NA, 2), trans = "sqrt")
  }
  ggsave(filename = paste0(channel, "_tSNE.pdf"), path = file.path(dea.output, "Mah.cluster tSNE plots"), 
         plot = tSNE.plot, device = "pdf", height = 6, width = 7)
}

#signaling plots
faceted.tSNE.data <- tSNE.data %>% 
  dplyr::select(Mah.cluster, tSNE1, tSNE2, condition, 
                one_of(c(SIGNALING.MARKERS)))

#function to use with purrr::map()
my.ggplot <- function(data = NULL, variable = NULL){ 
  variable <- sym(variable)
  my.plot <- ggplot(data = data, 
                    mapping = aes(x = tSNE1, y = tSNE2, color = !!variable)) + 
    geom_point(size = 1, alpha = 0.7) + theme_minimal() + facet_wrap(~condition)
  if(as_string(variable) %in% SIGNALING.MARKERS){ 
    my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.7), trans = "sqrt")
  }
  my.plot <- my.plot + theme(legend.position="none") + 
    annotate(geom = "label", x = 0, y = 50, label = as_string(variable), 
             label.size = 0, size = 5) +
    xlab("") + ylab("") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  return(my.plot)
}

#produce a list of ggplot objects to read into grid.arrange...
plot.list <- map(.x = colnames(faceted.tSNE.data), 
                 .f = my.ggplot, data = faceted.tSNE.data) 
names(plot.list) <- colnames(faceted.tSNE.data)
plot.list$tSNE1 <- NULL
plot.list$tSNE2 <- NULL

plot.list$`Cluster Type` <- NULL

final.plot <- do.call("grid.arrange", args = c(plot.list, 
                                               ncol = ceiling(sqrt(length(plot.list))), 
                                               nrow = ceiling(sqrt(length(plot.list)))))
ggsave(filename = paste0("Mah_signaling_faceted_plot.pdf"), 
       plot = final.plot, device = "pdf", width = 8, height = 8, path = dea.output)



#pERK plot + pSTAT3 plot
dim_red_matrix <- expression.matrix %>% dplyr::filter(patient %in% paired.samples) %>% 
  group_by(Mah.cluster, stimulation, condition) %>% 
  sample_n(size = 500, replace = F) %>% ungroup()

tSNE <- Rtsne(X = dim_red_matrix[1:41], dims = 2)
tSNE.data <- tSNE$Y %>% as_tibble()
tSNE.data <- cbind(tSNE.data, dim_red_matrix)
names(tSNE.data)[1:2] <- paste0("tSNE", 1:2)

tSNE.data <- tSNE.data %>% 
  mutate_at(.vars = vars(contains(match = "cluster")), .funs = as.factor)


#signaling plots
faceted.tSNE.data <- tSNE.data %>% 
  dplyr::select(Mah.cluster, tSNE1, tSNE2, condition, stimulation, 
                one_of(c(SIGNALING.MARKERS)))

my.ggplot <- function(data = NULL, variable = NULL){ 
  variable <- sym(variable)
  my.plot <- ggplot(data = data, 
                    mapping = aes(x = tSNE1, y = tSNE2, color = !!variable)) + 
    geom_point(size = 1, alpha = 0.7) + theme_minimal() + facet_wrap(stimulation~condition)
  if(as_string(variable) %in% SIGNALING.MARKERS){ 
    my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.2), trans = "sqrt")
  }
  my.plot <- my.plot + theme(legend.position="none") + 
    xlab("") + ylab("") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  return(my.plot)
}

#signal plot
pERK.plot <- my.ggplot(data = faceted.tSNE.data, variable = "pErk")

ggsave(path = dea.output, filename = "pERK_plot.pdf", 
       device = "pdf", width = 8, height = 8)

pSTAT3.plot <- my.ggplot(data = faceted.tSNE.data, variable = "pSTAT3")

ggsave(path = dea.output, filename = "pSTAT3_plot.pdf", 
       device = "pdf", width = 8, height = 8)

#cluster identity plot
faceted.tSNE.data$my.cluster <- ifelse(faceted.tSNE.data$Mah.cluster == 5, 
                                       yes = "MEP", no = "Not MEP")

variable <- sym("my.cluster")
Mah.cluster.plot <- ggplot(data = faceted.tSNE.data, 
                  mapping = aes(x = tSNE1, y = tSNE2, color = !!variable)) + 
  geom_point(size = 1, alpha = 0.7) + theme_minimal() + facet_wrap(~condition)
if(as_string(variable) %in% SIGNALING.MARKERS){ 
  my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.7), trans = "sqrt")
}
Mah.cluster.plot <- Mah.cluster.plot + scale_color_manual(values = c("red", "grey"))
Mah.cluster.plot <- Mah.cluster.plot + theme(legend.position="none") +  
  xlab("") + ylab("") + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
Mah.cluster.plot

ggsave(plot = Mah.cluster.plot, path = dea.output, filename = "MEP_plot.pdf", 
       device = "pdf", width = 8, height = 5)



#SAM with surface markers

form_dxrx_master <- function(abundances = NULL, centroids = NULL, derived = NULL, stims = NULL, cluster_method = NULL){ 
  master.matrix <- purrr::pluck(centroids, cluster_method) %>%
    dplyr::filter(condition != "Healthy") %>% 
    dplyr::filter(patient %in% paired.samples) %>% 
    gather(key = "variable", value = "value", -patient, -condition, -!!as.name(cluster_method)) %>% ungroup()
  master.matrix$combined.var <- paste(master.matrix$variable, master.matrix[[cluster_method]], sep = "_")
  master.matrix <- master.matrix %>% dplyr::select(patient, combined.var, value, condition) %>%
    tidyr::spread(key = combined.var, value = value)
  return(master.matrix)
}

dx.rx.matrices <- list()

for(i in 1:length(abundance.by.pat.con)){ 
  cluster_method <- names(abundance.by.pat.con)[[i]]
  new.matrix <- form_dxrx_master(abundances = abundance.by.pat.con, 
                                 centroids = centroids.by.pat.con, 
                                 derived = derived.metrics.by.pat.con, 
                                 stims = stims.by.pat.con, 
                                 cluster_method = cluster_method)
  new.matrix$patient[new.matrix$patient == "PASTRP"] <- "PASRTP"
  dx.rx.matrices[[cluster_method]] <- new.matrix
}

#combine with clinical data
sam.matrices <- map(dx.rx.matrices, function(data) inner_join(x = data, y = clinical.matrix))
pheno.sam.matrix <- sam.matrices$phenograph.metacluster 
sam.matrices$phenograph.metacluster <- NULL

sam_results <- list()

aml_sam <- function(sam.matrix = NULL){ 
  data <- sam.matrix %>% dplyr::select(contains("_"))
  feature.names <- colnames(data)
  response <- c(1:nrow(sam.matrix), 1:nrow(sam.matrix)) %>% sort() %>% 
    purrr::modify_at(.at = seq(from = 1, to = length(.), by = 2), .f = function(x) -x)
  result <- SAM(x = t(as.matrix(data)), y = response, 
                resp.type = "Two class paired", nperms = 5000, 
                geneid = feature.names, genenames = feature.names)
}



SAMs <- map(.x = sam.matrices, .f = aml_sam)








#SAM with signaling markers only in expanded phenograph clusters
pheno_aml_sam <- function(sam.matrix = NULL){ 
  data <- sam.matrix %>% dplyr::select(contains("_"))
  my.tags <- paste0("_", phenograph.expanded.clusters)
  feature.names <- colnames(data)
  data <- data %>% dplyr::select(matches(paste(my.tags, collapse = "|")))
  response <- c(1:nrow(sam.matrix), 1:nrow(sam.matrix)) %>% sort() %>% 
    purrr::modify_at(.at = seq(from = 1, to = length(response), by = 2), .f = function(x) -x)
  result <- SAM(x = t(as.matrix(data)), y = response, 
                resp.type = "Two class paired", nperms = 5000, 
                geneid = feature.names, genenames = feature.names)
}

pheno.SAM <- pheno_aml_sam(sam.matrix = pheno.sam.matrix)


#tSNE plots with expanded phenograph clusters




#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################



####################################         Using Edge-R         ##########################################

library(edgeR)

#form an m x n matrix, where m is the number of clusters for that clustering method and n is the number of samples being analyzed.
#In this matrix, the [i,j]th entry encodes the count of cells in the ith cluster for the jth sample. 
my.cells <- abundance.by.pat.con$Mah.cluster$num.cells
total.cells <- abundance.by.pat.con$Mah.cluster$total.cells

counts.table <- abundance.by.pat.con$Mah.cluster %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  dcast(formula = Mah.cluster ~ patient + condition, 
        value.var = "num.cells", 
        fill = 0)

counts.matrix <- as.matrix(counts.table[,2:ncol(counts.table)])
row.names(counts.matrix) <- counts.table$Mah.cluster
colnames(counts.matrix) <- rep(x = c("Dx", "Rx"), times = ncol(counts.matrix)/2)
conditions <- rep(x = c("Dx", "Rx"), times = ncol(counts.matrix)/2)

dge <- DGEList(counts = counts.matrix)

design <- model.matrix(~factor(conditions))

y <- estimateDisp(dge, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=2)

sampled.matrix <- expression.matrix %>% ungroup() %>% 
  group_by(Mah.cluster) %>% sample_n(size = 5000)

umap <- umap(d = sampled.matrix[1:41])
tSNE <- Rtsne(X = sampled.matrix[1:41])

tSNE.mult <- Rtsne.multicore(X = sampled.matrix[1:41], num_threads = 10)

#plots 

#settings
par(mfrow=c(1,2))
layout(matrix(c(1,1,1,2,2,2,3), nrow=1,ncol=7))

#first plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(length(CLASSIFIER.POPULATIONS))
plot(umap$layout[,1:2], pch=21, 
     bg=alpha(col=cols[sampled.matrix$Mah.cluster],0.6), 
     col="black", xlab="", ylab="", lwd=0.5, main="Mah clusters")

#second plot
log2fc <- res$table$logFC
log2fc[res$table$PValue>0.05] <- 0
colpal <- colorRampPalette(c('red', "grey", 'blue'))
cols <- colpal(100)[as.numeric(cut(log2fc,breaks = 100))]
plot(umap$layout[,1:2], pch=21, bg=alpha(col=cols[sampled.matrix$Mah.cluster],0.6), col="grey", xlab="", ylab="", lwd=0.5, main="colored by log2FC")

legend_image <- as.raster(matrix(colpal(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'log2FC')
text(x=1.5, y = seq(0,1,l=3), labels = c(round(min(log2fc), digits = 2),0,round(max(log2fc),digits=2)))
rasterImage(legend_image, 0, 0, 1, 1)

###now with tSNE 
#settings
par(mfrow=c(1,2))
layout(matrix(c(1,1,1,2,2,2,3), nrow=1,ncol=7))

#first plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(length(CLASSIFIER.POPULATIONS))
plot(tSNE$Y, pch=21, 
     bg=alpha(col=cols[sampled.matrix$Mah.cluster],0.6), 
     col="black", xlab="", ylab="", lwd=0.5, main="Mah clusters")

#second plot
log2fc <- res$table$logFC
log2fc[res$table$PValue>0.05] <- 0
colpal <- colorRampPalette(c('red', "grey", 'blue'))
cols <- colpal(100)[as.numeric(cut(log2fc,breaks = 100))]
plot(tSNE$Y, pch=21, bg=alpha(col=cols[sampled.matrix$Mah.cluster],0.6), col="grey", xlab="", ylab="", lwd=0.5, main="colored by log2FC")

legend_image <- as.raster(matrix(colpal(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'log2FC')
text(x=1.5, y = seq(0,1,l=3), labels = c(round(min(log2fc), digits = 2),0,round(max(log2fc),digits=2)))
rasterImage(legend_image, 0, 0, 1, 1)



#################         using phenograph clusters       ###################
my.cells <- abundance.by.pat.con$phenograph.metacluster$num.cells
total.cells <- abundance.by.pat.con$phenograph.metacluster$total.cells

counts.table <- abundance.by.pat.con$phenograph.metacluster %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  dcast(formula = phenograph.metacluster ~ patient + condition, 
        value.var = "num.cells", 
        fill = 0)

counts.matrix <- as.matrix(counts.table[,2:ncol(counts.table)])
row.names(counts.matrix) <- counts.table$phenograph.metacluster
colnames(counts.matrix) <- rep(x = c("Dx", "Rx"), times = ncol(counts.matrix)/2)
conditions <- rep(x = c("Dx", "Rx"), times = ncol(counts.matrix)/2)

dge <- DGEList(counts = counts.matrix)

design <- model.matrix(~factor(conditions))

y <- estimateDisp(dge, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=2)

sampled.matrix <- expression.matrix %>% ungroup() %>% 
  group_by(phenograph.metacluster) %>% sample_n(size = 500)

umap <- umap(d = sampled.matrix[1:41])

tSNE <- Rtsne(X = sampled.matrix[1:41], 
              dims = 3, perplexity = 70)



#plots 

par(mfrow=c(1,2))
layout(matrix(c(1,1,1,2,2,2,3), nrow=1,ncol=7))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(length(unique(expression.matrix$phenograph.metacluster)))
plot(umap$layout[,1:2], pch=21, 
     bg=alpha(col=cols[sampled.matrix$phenograph.metacluster],0.6), 
     col="black", xlab="", ylab="", lwd=0.1, main="Phenograph clusters")

#second plot

log2fc <- res$table$logFC
log2fc[res$table$PValue>0.05] <- 0
colpal <- colorRampPalette(c('red', "grey", 'blue'))
cols <- colpal(100)[as.numeric(cut(log2fc,breaks = 100))]
plot(umap$layout[,1:2], pch=21, bg=alpha(col=cols[sampled.matrix$phenograph.metacluster],0.8), 
     col="black", xlab="", ylab="", lwd=0.2, main="colored by log2FC")

legend_image <- as.raster(matrix(colpal(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'log2FC')
text(x=1.5, y = seq(0,1,l=3), labels = c(round(min(log2fc), digits = 2),0,round(max(log2fc),digits=2)))
rasterImage(legend_image, 0, 0, 1, 1)


##Now with tSNE

par(mfrow=c(1,2))
layout(matrix(c(1,1,1,2,2,2,3), nrow=1,ncol=7))

#first plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(length(unique(expression.matrix$phenograph.metacluster)))
plot(tSNE$Y, pch=21, 
     bg=alpha(col=cols[sampled.matrix$phenograph.metacluster],0.6), 
     col="black", xlab="", ylab="", lwd=0.1, main="Phenograph clusters")

#second plot

log2fc <- res$table$logFC
log2fc[res$table$PValue>0.05] <- 0
colpal <- colorRampPalette(c('red', "grey", 'blue'))
cols <- colpal(100)[as.numeric(cut(log2fc,breaks = 100))]
plot(tSNE$Y, pch=21, bg=alpha(col=cols[sampled.matrix$phenograph.metacluster],0.8), 
     col="grey", xlab="", ylab="", lwd=0.2, main="colored by log2FC")

legend_image <- as.raster(matrix(colpal(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'log2FC')
text(x=1.5, y = seq(0,1,l=3), labels = c(round(min(log2fc), digits = 2),0,round(max(log2fc),digits=2)))
rasterImage(legend_image, 0, 0, 1, 1)

my.clusters <- which(res$table$PValue<0.05) 
table(expression.matrix$phenograph.metacluster)

plotly.frame <- cbind(tSNE$Y, 
                      sampled.matrix$phenograph.metacluster) %>% as_tibble()
colnames(plotly.frame) <- c("tSNE1", "tSNE2", "tSNE3", "phenograph.metacluster")

my.clusters <- which(res$table$PValue<0.05)
my.clusters <- row.names(counts.matrix)[my.clusters]
plotly.frame$`cluster type` <- ifelse(plotly.frame$phenograph.metacluster %in% my.clusters, 
                                      yes = "expanded", 
                                      no = "not expanded") 

p <- plot_ly(plotly.frame, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, color = ~`cluster type`, 
             colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'tSNE1'),
                      yaxis = list(title = 'tSNE2'),
                      zaxis = list(title = 'tSNE3')))


#############################################################################################################
#############################################################################################################
#############################################################################################################


#interrogate what the significantly different clusters actually are
my.clusters <- which(res$table$PValue<0.05)
my.clusters <- row.names(counts.matrix)[my.clusters]

#kind of hack-y, but used to interrogate phenograph clusters more until I understand the 
#modeling component 
my.clusters <- phenograph.expanded.clusters

#composition of clusters by patient 
cluster_totals <- expression.matrix %>% ungroup() %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  dplyr::filter(phenograph.metacluster %in% my.clusters) %>% 
  group_by(phenograph.metacluster, condition) %>% 
  summarize(cluster.total.cells = n())

patient_totals <- expression.matrix %>% ungroup() %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  ungroup() %>% 
  group_by(patient, condition) %>% 
  summarize(patient.total.cells = n())

#Rx - Dx 
cluster.differences <- cluster_totals %>% 
  tidyr::spread(key = condition, value = cluster.total.cells, fill = 0) %>% 
  dplyr::transmute(difference = Rx - Dx)

#expanded.clusters <- cluster.differences$phenograph.metacluster[cluster.differences$difference > 0]
#shrunken.clusters <- cluster.differences$phenograph.metacluster[cluster.differences$difference < 0]

cluster.proportions <- expression.matrix %>% ungroup() %>%
  dplyr::filter(patient %in% paired.samples) %>% 
  dplyr::filter(phenograph.metacluster %in% my.clusters) %>% 
  group_by(patient, condition, phenograph.metacluster) %>% 
  summarize(cluster.cells = n()) %>% left_join(cluster_totals) %>% 
  dplyr::mutate(proportion = cluster.cells/cluster.total.cells)
cluster.proportions$phenograph.metacluster <- 
  factor(as.character(cluster.proportions$phenograph.metacluster), 
         levels = as.character(sort(unique(cluster.proportions$phenograph.metacluster))))

missing.clusters <- cluster.proportions %>% 
  dplyr::filter(condition == "Dx")
missing.clusters <- my.clusters[!(my.clusters %in% unique(missing.clusters$phenograph.metacluster))]
for(cluster in missing.clusters){
  new.row <- tibble(patient = "PARBIU", condition = "Dx", 
                    phenograph.metacluster = cluster,
                    cluster.cells = 0, 
                    cluster.total.cells = 1, 
                    proportion = 0)
  cluster.proportions <- bind_rows(cluster.proportions, new.row)
}
cluster.proportions$phenograph.metacluster <- factor(cluster.proportions$phenograph.metacluster, 
                                                     levels = as.character(sort(as.numeric(unique(cluster.proportions$phenograph.metacluster)))))

#proportion of clusters by patient in Dx and Rx
proportion.plot <- ggplot(data = cluster.proportions, 
                          mapping = aes(x = phenograph.metacluster, 
                                        y = proportion, 
                                        fill = patient)) + 
  geom_bar(stat = "identity") + facet_wrap(~condition, scales = "free_x") + 
  theme_bw() + xlab("Phenograph Metacluster")
proportion.plot

ggsave(filename = "patient_proportion_plot.pdf", device = "pdf", 
       plot = proportion.plot, 
       width = 10, height = 6, 
       path = daa.output)

proportion.plot <- ggplot(data = cluster.proportions, 
                          mapping = aes(x = condition, 
                                        y = proportion, 
                                        fill = patient)) + 
  geom_bar(stat = "identity") + facet_wrap(~phenograph.metacluster, scales = "free_x") + 
  theme_bw() + xlab("Phenograph Metacluster") + scale_x_discrete(drop = F)
proportion.plot

ggsave(filename = "patient_proportion_plot_2.pdf", device = "pdf", 
       plot = proportion.plot, 
       width = 6, height = 10, 
       path = daa.output)

#summed percentage of each patient in each cluster by Dx and Rx
cluster.proportions <- cluster.proportions %>% left_join(patient_totals) %>% 
  mutate(`Proportion of Patient Cells` = cluster.cells/patient.total.cells)

missing.clusters <- cluster.proportions %>% 
  dplyr::filter(condition == "Dx")
missing.clusters <- my.clusters[!(my.clusters %in% unique(missing.clusters$phenograph.metacluster))]
for(cluster in missing.clusters){
  new.row <- tibble(patient = "PARBIU", 
                    phenograph.metacluster = cluster, condition = "Dx",
                    cluster.cells = 0, 
                    cluster.total.cells = 0, 
                    proportion = 0, patient.total.cells = 0, 
                    `Proportion of Patient Cells` = 0)
  cluster.proportions <- rbind(cluster.proportions, new.row)
}

proportion.plot <- ggplot(data = cluster.proportions, 
                          mapping = aes(x = phenograph.metacluster, 
                                        y = `Proportion of Patient Cells`, 
                                        fill = patient)) + 
  geom_bar(stat = "identity") + facet_wrap(~condition, scales = "free_x") + 
  theme_bw() + xlab("Phenograph Metacluster")
proportion.plot

ggsave(filename = "summed_patient_proportion_plot_expanded_only.pdf", device = "pdf", 
       plot = proportion.plot, 
       width = 10, height = 6, 
       path = daa.output)

#summed percentage of cells in the expanded clusters vs. non-expanded clusters by Dx and Rx

#recalulate cluster totals from all clusters, not just expanded ones 
cluster_totals <- expression.matrix %>% ungroup() %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  group_by(phenograph.metacluster, condition) %>% 
  summarize(cluster.total.cells = n())

patient_totals <- expression.matrix %>% ungroup() %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  ungroup() %>% 
  group_by(patient, condition) %>% 
  summarize(patient.total.cells = n())

cluster.proportions <- expression.matrix %>% ungroup() %>%
  dplyr::filter(patient %in% paired.samples) %>% 
  group_by(patient, condition, phenograph.metacluster) %>% 
  summarize(cluster.cells = n()) %>% left_join(cluster_totals) %>% left_join(patient_totals) %>% 
  dplyr::mutate(cluster.proportion = cluster.cells/cluster.total.cells, 
                `Proportion of Patient Cells` = cluster.cells/patient.total.cells)

cluster.proportions$phenograph.metacluster <- 
  factor(as.character(cluster.proportions$phenograph.metacluster), 
         levels = as.character(sort(unique(cluster.proportions$phenograph.metacluster))))

#new categorical variable
cluster.proportions$`Cluster Type` <- ifelse(test = cluster.proportions$phenograph.metacluster %in% my.clusters, 
                                             yes = "Expanded", no = "Not Expanded")

proportion.plot <- ggplot(data = cluster.proportions, 
                          mapping = aes(x = patient, 
                                        y = `Proportion of Patient Cells`, 
                                        fill = `Cluster Type`)) + 
  geom_bar(stat = "identity") + facet_wrap(~condition, scales = "free_x") + 
  theme_bw() + xlab("Phenograph Metacluster") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
proportion.plot

ggsave(filename = "summed_patient_expanded_versus_nonexpanded.pdf", device = "pdf", 
       plot = proportion.plot, 
       width = 10, height = 6, 
       path = daa.output)

#composition of patients by cluster

#only expanded
patient_totals <- expression.matrix %>% ungroup %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  dplyr::filter(phenograph.metacluster %in% my.clusters) %>% 
  group_by(patient, condition) %>% 
  summarize(total.cells = n())

patient.proportions <- expression.matrix %>% ungroup() %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  dplyr::filter(phenograph.metacluster %in% my.clusters) %>% 
  group_by(patient, condition, phenograph.metacluster) %>% 
  summarize(patient.cells = n()) %>% left_join(patient_totals) %>% 
  dplyr::mutate(proportion = patient.cells/total.cells)
patient.proportions$phenograph.metacluster <- 
  factor(as.character(patient.proportions$phenograph.metacluster), 
         levels = as.character(sort(unique(patient.proportions$phenograph.metacluster))))
patient.proportions$patient <- as.factor(as.character(patient.proportions$patient))

my.colors <- colorRampPalette(colors = brewer.pal(11, name = "Spectral"))(length(unique(patient.proportions$phenograph.metacluster)))

proportion.plot <- ggplot(data = patient.proportions, 
                          mapping = aes(x = patient, 
                                        y = proportion, 
                                        fill = phenograph.metacluster)) + 
  geom_bar(stat = "identity") + facet_wrap(~condition, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_fill_manual(values = my.colors) + scale_x_discrete(drop = F)
proportion.plot

ggsave(filename = "patient_proportion_plot_expanded_only.pdf", device = "pdf", 
       plot = proportion.plot, 
       width = 12, height = 8, path = daa.output)

#all clusters
patient_totals <- expression.matrix %>% ungroup %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  group_by(patient, condition) %>% 
  summarize(total.cells = n())

patient.proportions <- expression.matrix %>% ungroup() %>% 
  dplyr::filter(patient %in% paired.samples) %>% 
  group_by(patient, condition, phenograph.metacluster) %>% 
  summarize(patient.cells = n()) %>% left_join(patient_totals) %>% 
  dplyr::mutate(proportion = patient.cells/total.cells)
patient.proportions$phenograph.metacluster <- 
  factor(as.character(patient.proportions$phenograph.metacluster), 
         levels = as.character(sort(unique(patient.proportions$phenograph.metacluster))))
patient.proportions$patient <- as.factor(as.character(patient.proportions$patient))

my.colors <- colorRampPalette(colors = brewer.pal(length(unique(patient.proportions$phenograph.metacluster)),
                                                  name = "Spectral"))(length(unique(patient.proportions$phenograph.metacluster)))

proportion.plot <- ggplot(data = patient.proportions, 
                          mapping = aes(x = patient, 
                                        y = proportion, 
                                        fill = phenograph.metacluster)) + 
  geom_bar(stat = "identity") + facet_wrap(~condition, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) + 
  scale_fill_manual(values = my.colors) + scale_x_discrete(drop = F)
proportion.plot

ggsave(filename = "patient_proportion_plot_all.pdf", device = "pdf", 
       plot = proportion.plot, 
       width = 12, height = 8, path = daa.output)

#dimensionality reduction for all clusters with expanded ones indicated 
dim_red_matrix <- expression.matrix %>% dplyr::filter(patient %in% paired.samples) %>% 
  group_by(phenograph.metacluster) %>% 
  dplyr::filter(phenograph.metacluster != 4) %>%
  sample_n(size = 500, replace = TRUE) %>% ungroup() %>% unique()

tSNE <- Rtsne(X = dim_red_matrix[1:41], dims = 2)
tSNE.data <- tSNE$Y %>% as_tibble()
tSNE.data <- cbind(tSNE.data, dim_red_matrix)
names(tSNE.data)[1:2] <- paste0("tSNE", 1:2)

tSNE.data <- tSNE.data %>% 
  mutate_at(.vars = vars(contains(match = "cluster")), .funs = as.factor)
tSNE.data$`Cluster Type` <- as.factor(ifelse(test = (tSNE.data$phenograph.metacluster %in% my.clusters),
                                             yes = "Expanded", no = "Not Expanded"))

dir.create(file.path(daa.output, "tSNE plots expanded"))
for(i in 1:length(colnames(tSNE.data))){ 
  channel.name <- colnames(tSNE.data)[i]
  channel <- sym(colnames(tSNE.data)[i])
  tSNE.plot <- ggplot(data = tSNE.data,
                      mapping = aes(x = tSNE1, y = tSNE2, color = !!channel)) + 
    geom_point() + theme_bw() 
  if(is.numeric(tSNE.data[,channel.name])){
    tSNE.plot <- tSNE.plot + 
      scale_color_viridis(limits = c(NA, 2), trans = "sqrt")
  }
  ggsave(filename = paste0(channel, "_tSNE.pdf"), path = file.path(daa.output, "tSNE plots"), 
         plot = tSNE.plot, device = "pdf", height = 6, width = 7)
}

#dimensionality reduction for only the clusters that have been found to be significant 
dim_red_matrix <- expression.matrix %>% dplyr::filter(patient %in% paired.samples) %>% 
  dplyr::filter(phenograph.metacluster %in% my.clusters) %>% group_by(phenograph.metacluster) %>% 
  sample_n(size = 500, replace = TRUE) %>% ungroup() %>% unique()

tSNE <- Rtsne(X = dim_red_matrix[1:41], dims = 2)
tSNE.data <- tSNE$Y %>% as_tibble()
tSNE.data <- cbind(tSNE.data, dim_red_matrix)
names(tSNE.data)[1:2] <- paste0("tSNE", 1:2)

tSNE.data <- tSNE.data %>% mutate_at(.vars = vars(contains(match = "cluster")), .funs = as.factor)
tSNE.data$`Cluster Type` <- as.factor(ifelse(test = (tSNE.data$phenograph.metacluster %in% my.clusters),
                                           yes = "Expanded", no = "Not Expanded"))

#create a bunch of separate tSNE plots
dir.create(file.path(daa.output, "tSNE plots expanded"))
for(i in 1:length(colnames(tSNE.data))){ 
  channel <- sym(colnames(tSNE.data)[i])
  tSNE.plot <- ggplot(data = tSNE.data,
                      mapping = aes(x = tSNE1, y = tSNE2, color = !!channel)) + 
    geom_point() + theme_bw() 
  if(is.numeric(tSNE.data[,channel.name])){
    tSNE.plot <- tSNE.plot + 
      scale_color_viridis(limits = c(NA, 2), trans = "sqrt")
  }
  ggsave(filename = paste0(channel, '_tSNE_expanded.pdf'), path = file.path(daa.output, "tSNE plots expanded"),
         plot = tSNE.plot, device = "pdf", height = 6, width = 7)
}

#create a giant faceted plot of just the LSC-associated markers

LSC.markers <- c("CD123", "CD47", "CD90", "CD93", "CD117", "CD123", "CD34", "CD38", "TIM-3") 

faceted.tSNE.data <- tSNE.data %>% 
  dplyr::select(phenograph.metacluster, `Cluster Type`, tSNE1, tSNE2, one_of(LSC.markers)) 


#function to use with purrr::map()
my.ggplot <- function(data = NULL, variable = NULL){ 
  variable <- sym(variable)
  my.plot <- ggplot(data = data, 
                    mapping = aes(x = tSNE1, y = tSNE2, color = !!variable)) + 
    geom_point(size = 1, alpha = 0.7) + theme_minimal()
  if(as_string(variable) %in% LSC.markers){ 
    my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.5), trans = "sqrt")
  }
  my.plot <- my.plot + theme(legend.position="none") + 
    annotate(geom = "label", x = 0, y = 50, label = as_string(variable))
  return(my.plot)
}

#produce a list of ggplot objects to read into grid.arrange...
plot.list <- map(.x = colnames(faceted.tSNE.data), .f = my.ggplot, data = faceted.tSNE.data) 
names(plot.list) <- colnames(faceted.tSNE.data)
plot.list$tSNE1 <- NULL
plot.list$tSNE2 <- NULL

plot.list$`Cluster Type` <- NULL

final.plot <- do.call("grid.arrange", args = c(plot.list, ncol = 3, nrow = 3))
ggsave(filename = paste0("stem_cell_marker_faceted_plot.pdf"), 
       plot = final.plot, device = "pdf", width = 8, height = 8, path = dea.output)

#### with all lineage markers 

faceted.tSNE.data <- tSNE.data %>% 
  dplyr::select(phenograph.metacluster, 
                `Cluster Type`, tSNE1, tSNE2, 
                one_of(c(SURFACE.MARKERS, TRX.FACTORS, "MPO")))

#function to use with purrr::map()
my.ggplot <- function(data = NULL, variable = NULL){ 
  variable <- sym(variable)
  my.plot <- ggplot(data = data, 
                    mapping = aes(x = tSNE1, y = tSNE2, color = !!variable)) + 
    geom_point(size = 1, alpha = 0.7) + theme_minimal()
  if(as_string(variable) %in% c(SURFACE.MARKERS, TRX.FACTORS, "MPO")){ 
    my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.7), trans = "sqrt")
  }
  my.plot <- my.plot + theme(legend.position="none") + 
    annotate(geom = "label", x = 0, y = 50, label = as_string(variable), 
             label.size = 0, size = 2) +
    xlab("") + ylab("") + 
    theme(axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())
  return(my.plot)
}

#produce a list of ggplot objects to read into grid.arrange...
plot.list <- map(.x = colnames(faceted.tSNE.data), 
                 .f = my.ggplot, data = faceted.tSNE.data) 
names(plot.list) <- colnames(faceted.tSNE.data)
plot.list$tSNE1 <- NULL
plot.list$tSNE2 <- NULL

plot.list$`Cluster Type` <- NULL

final.plot <- do.call("grid.arrange", args = c(plot.list, 
                                               ncol = ceiling(sqrt(length(plot.list))), 
                                               nrow = ceiling(sqrt(length(plot.list)))))
ggsave(filename = paste0("all_markers_faceted_plot.pdf"), 
       plot = final.plot, device = "pdf", width = 10, height = 10, path = dea.output)


#### with all signaling markers 


faceted.tSNE.data <- tSNE.data %>% 
  dplyr::select(phenograph.metacluster, 
                `Cluster Type`, tSNE1, tSNE2, 
                one_of(c(SIGNALING.MARKERS)))

#function to use with purrr::map()
my.ggplot <- function(data = NULL, variable = NULL){ 
  variable <- sym(variable)
  my.plot <- ggplot(data = data, 
                    mapping = aes(x = tSNE1, y = tSNE2, color = !!variable)) + 
    geom_point(size = 1, alpha = 0.7) + theme_minimal()
  if(as_string(variable) %in% SIGNALING.MARKERS){ 
    my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.7), trans = "sqrt")
  }
  my.plot <- my.plot + theme(legend.position="none") + 
    annotate(geom = "label", x = 0, y = 50, label = as_string(variable), 
             label.size = 0, size = 5) +
    xlab("") + ylab("") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  return(my.plot)
}

#produce a list of ggplot objects to read into grid.arrange...
plot.list <- map(.x = colnames(faceted.tSNE.data), 
                 .f = my.ggplot, data = faceted.tSNE.data) 
names(plot.list) <- colnames(faceted.tSNE.data)
plot.list$tSNE1 <- NULL
plot.list$tSNE2 <- NULL

plot.list$`Cluster Type` <- NULL

final.plot <- do.call("grid.arrange", args = c(plot.list, 
                                               ncol = ceiling(sqrt(length(plot.list))), 
                                               nrow = ceiling(sqrt(length(plot.list)))))
ggsave(filename = paste0("signaling_faceted_plot.pdf"), 
       plot = final.plot, device = "pdf", width = 8, height = 8, path = dea.output)


#fanagle the legends 
#legend <- cowplot::get_legend(plot.list[[1]] + theme(legend.position = "right"))


#create directory for UMAP plots
# dir.create(file.path("~", "Desktop", "UMAP plots"))
# umap <- umap(d = dim_red_matrix[1:41], config = my.umap)
# umap.data <- umap$layout %>% as_tibble()
# umap.data <- cbind(umap$layout, dim_red_matrix)
# names(umap.data)[1:3] <- paste0("umap", 1:3)
# 
# umap.data <- umap.data %>% mutate_at(.vars = vars(contains(match = "cluster")), .funs = as.factor)
# 
# umap.data$`Cluster Type` <- as.factor(ifelse(test = (umap.data$phenograph.metacluster %in% expanded.clusters),
#                                              yes = "Expanded", no = "Shrunken"))
# 
# for(i in 1:length(colnames(umap.data))){ 
#   channel <- sym(colnames(umap.data)[i])
#   umap.plot <- ggplot(data = umap.data,
#                       mapping = aes(x = umap1, y = umap2, color = !!channel)) + 
#     geom_point() + theme_bw() 
#   ggsave(filename = file.path("~", "Desktop","UMAP plots", paste0(channel, "_", "umap_plot.pdf")), 
#          plot = umap.plot, device = "pdf")
# }
# 
# #giant faceted plot
# faceted.umap.data <- umap.data %>% 
#   dplyr::select(phenograph.metacluster, `Cluster Type`, umap1, umap2, one_of(LSC.markers))
# 
# #function to use with purrr::map()
# my.ggplot <- function(data = NULL, variable = NULL){ 
#   variable <- sym(variable)
#   my.plot <- ggplot(data = data, 
#                     mapping = aes(x = umap1, y = umap2, color = !!variable)) + 
#     geom_point(size = 1) + theme_minimal()
#   if(as_string(variable) %in% LSC.markers){ 
#     my.plot <- my.plot + scale_color_viridis(limits = c(NA, 1.5), trans = "sqrt")
#   }
#   my.plot <- my.plot + theme(legend.position="none") + 
#     annotate(geom = "label", x = 0, y = 8, label = as_string(variable))
#   return(my.plot)
# } 
# 
# plot.list <- map(.x = colnames(faceted.umap.data), .f = my.ggplot, data = faceted.umap.data) 
# names(plot.list) <- colnames(faceted.umap.data)
# plot.list$umap1 <- NULL
# plot.list$umap2 <- NULL
# 
# final.plot <- do.call("grid.arrange", args = c(plot.list, ncol = 3, nrow = 3))
# ggsave(filename = file.path("~", "Desktop", "umap_faceted_plot.pdf"), 
#        plot = final.plot, device = "pdf")


#Differential Expression Analysis between clusters using SAM 
#####  With only signaling - have to play with this because of all the missing values...not super functional atm
form_dxrx_master <- function(abundances = NULL, centroids = NULL, derived = NULL, stims = NULL, cluster_method = NULL){ 
  master.matrix <- purrr::pluck(stims, cluster_method) %>%
    dplyr::filter(condition != "Healthy") %>% 
    dplyr::filter(patient %in% paired.samples) %>% 
    gather(key = "variable", value = "value", -patient, -condition, -!!as.name(cluster_method)) %>% ungroup()
  master.matrix$combined.var <- paste(master.matrix$variable, master.matrix[[cluster_method]], sep = "_")
  master.matrix <- master.matrix %>% dplyr::select(patient, combined.var, value, condition) %>%
    tidyr::spread(key = combined.var, value = value)
  return(master.matrix)
}

dx.rx.matrices <- list()

for(i in 1:length(abundance.by.pat.con)){ 
  cluster_method <- names(abundance.by.pat.con)[[i]]
  new.matrix <- form_dxrx_master(abundances = abundance.by.pat.con, 
                                 centroids = centroids.by.pat.con, 
                                 derived = derived.metrics.by.pat.con, 
                                 stims = stims.by.pat.con, 
                                 cluster_method = cluster_method)
  new.matrix$patient[new.matrix$patient == "PASTRP"] <- "PASRTP"
  dx.rx.matrices[[cluster_method]] <- new.matrix
}

#combine with clinical data
sam.matrices <- map(dx.rx.matrices, function(data) inner_join(x = data, y = clinical.matrix))
pheno.sam.matrix <- sam.matrices$phenograph.metacluster 
sam.matrices$phenograph.metacluster <- NULL

sam_results <- list()

aml_sam <- function(sam.matrix = NULL){ 
  sam.matrix <- sam.matrix %>% ungroup()
  data <- as.matrix(sam.matrix[4:ncol(sam.matrix)])
  feature.names <- colnames(data)
  response <- as.numeric(as.factor(pluck(sam.matrix,1)))
  result <- SAM(x = t(as.matrix(data)), y = response, 
                resp.type = "Multiclass", nperms = 5000, 
                geneid = feature.names, genenames = feature.names)
}

SAMs <- map(.x = centroids.by.pat.con, .f = aml_sam)

#heatmaps of surface expression
annotation_colors <- list(Mah.cluster = getPalette(n = length(unique(expression.matrix$Mah.cluster))))

heatmap.matrix <- centroids.by.pat.con$Mah.cluster %>% ungroup() %>% group_by(Mah.cluster) %>% 
  summarize_all(mean) 
row_names <- heatmap.matrix$Mah.cluster
heatmap.matrix <- as.matrix(heatmap.matrix[,4:ncol(heatmap.matrix)])
row.names(heatmap.matrix) <- row_names


my.heatmap <- pheatmap(mat = heatmap.matrix, 
                       cluster_cols = TRUE, cluster_rows = FALSE,
                       show_rownames = TRUE, 
                       annotation_colors = annotation_colors, 
                       cell_width = 7, cell_height = 5, scale = "column")
ggsave(filename = "Mah_cluster_marker_heatmap.pdf", path = dea.output, device = "pdf", 
       width = 10, height = 6, plot = my.heatmap)

#histograms?
annotation_colors <- list(Mah.cluster = getPalette(n = length(unique(expression.matrix$phenograph.metacluster))))

heatmap.matrix <- centroids.by.pat.con$phenograph.metacluster %>% ungroup() %>% group_by(phenograph.metacluster) %>% 
  summarize_all(mean) 
row_names <- heatmap.matrix$phenograph.metacluster
heatmap.matrix <- as.matrix(heatmap.matrix[,4:ncol(heatmap.matrix)])
row.names(heatmap.matrix) <- row_names

my.heatmap <- pheatmap(mat = heatmap.matrix, 
                       cluster_cols = TRUE, cluster_rows = TRUE,
                       show_rownames = TRUE, 
                       annotation_colors = annotation_colors, 
                       cell_width = 7, cell_height = 5, scale = "column")
ggsave(filename = "phenograph_cluster_marker_heatmap.pdf", path = dea.output, device = "pdf", 
       width = 12, height = 10, plot = my.heatmap)


#plot all clusters with DEA between conditions

#Plot all clusters with DEA between 

