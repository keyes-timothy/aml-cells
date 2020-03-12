###AML_clustering.R
###implement clustering algorithms for all AML data 
#input: expression matrix of m cells by n channels
#output: none
#side-effects: save expression matrix with columns indicating which clusters cells belong to. Clustering methods 
#              used include Phenograph, FlowSOM, and k-means
#optimizations: make phenograph fun faster? 

message("Setting up clustering run...")
#libraries
library(Rphenograph)
library(FlowSOM)
library(foreach)
library(doParallel)
library(purrr)
#Xshift
#hierarchical

#####   directories   #####

clustering.output <- file.path(OUTPUT.DIRECTORY, "clustering_results")

###########      phenograph      ############
num_cores <- 30 #can be changed

expression.matrix$Mah.cluster <- expression.matrix$MahID
expression.matrix$MahID <- NULL
surface.cols <- colnames(expression.matrix) %in% SURFACE.MARKERS
signaling.cols <- colnames(expression.matrix) %in% SIGNALING.MARKERS
mah.dist.cols <- dplyr::contains(match = "Mah-", vars = colnames(expression.matrix))
cluster.cols <- dplyr::contains(match = "cluster", vars = colnames(expression.matrix))
metadata.cols <- dplyr::one_of(c("stimulation", "condition", "patient"), .vars = colnames(expression.matrix))
other.cols <- !(colnames(expression.matrix) %in% c(SURFACE.MARKERS, 
                                                  colnames(expression.matrix)[mah.dist.cols],
                                                  colnames(expression.matrix)[cluster.cols], 
                                                  SIGNALING.MARKERS,
                                                  colnames(expression.matrix)[metadata.cols]))

#normalize channels - divide each cell's measurement of a given protein by the 99.5th percentile of healthy
message("Normalizing data")

max.healthies <- expression.matrix %>% ungroup() %>% dplyr::filter(condition == "Healthy") %>% 
  dplyr::select(c(which(surface.cols), which(signaling.cols), which(other.cols))) %>% 
  apply(MARGIN = 2, FUN = quantile, probs = c(0.995))

for(i in 1:length(max.healthies)){ 
  expression.matrix[,names(max.healthies)[[i]]] <- expression.matrix[,names(max.healthies)[[i]]]/max.healthies[[i]]
}

#use the diagnostic plots function here too, to see how rescaling has changed the histograms for each channel. 
# FILL IN 
# FILL IN 
# FILL IN 

#define function to phenograph cluster a single df taken from a list of dataframes 
phenoPatient <- function(my.df = NULL, chosen.cols = surface.cols){ 
  message("Starting phenograph analysis of ", 
          unique(my.df$patient), " ",
          unique(my.df$condition), " ", 
          unique(my.df$stimulation))
  my.phenograph <- Rphenograph(data = my.df[,which(chosen.cols)], k = 15)
  my.clusters <- membership(my.phenograph[[2]]) %>% as.numeric()
  if(length(my.clusters) != nrow(my.df)){ 
    my.string <- as.character(as.numeric(1:nrow(my.df)))
    bad.cells <- which(!(my.string %in% names(membership(my.phenograph[[2]]))))
    bad.cell.clusters <- rep(NA, length(bad.cells)) 
    names(bad.cell.clusters) <- as.character(bad.cells)
    my.clusters <- c(membership(my.phenograph[[2]]), bad.cell.clusters)
    my.clusters <- my.clusters[order(as.numeric(names(my.clusters)))] %>% as.numeric() #check that this works
    message("Had to add NA phenoClusters in ", 
            unique(my.df$patient), " ",
            unique(my.df$condition), " ", 
            unique(my.df$stimulation))
  }
  my.df <- cbind(my.df, my.clusters)
  my.result = my.df
  #filename <- paste(my.df$patient[[1]], my.df$condition[[1]], my.df$stimulation[[1]], sep = "_")
  #filename <- paste0(filename, ".RDS")
  #saveRDS(object = my.result, file = file.path(clustering.output, filename))
  message("Returning phenograph analysis of ", 
          unique(my.df$patient), " ",
          unique(my.df$condition), " ", 
          unique(my.df$stimulation))
  return(list(my.result))
}

#run phenograph on each patient individually (in parallel)
message("Running PhenoGraph on all data")

df.list <- expression.matrix %>% dplyr::ungroup() %>% 
  group_by(patient, condition, stimulation) %>% 
  group_split() #divide up expression.matrix into data collected from each patient, condition, and stim

my.clust <- makeCluster(num_cores, outfile = "")
registerDoParallel(my.clust)

my.phenoResult <- foreach(fcs.df = df.list, .combine = c, 
                          .packages = c("dplyr", "Rphenograph")) %dopar% 
  phenoPatient(my.df = fcs.df, chosen.cols = surface.cols)

stopCluster(my.clust)

message("Consolidating PhenoData from separate cores")
#get a data matrix in which each row is a cluster identified in the single-patient phenograph runs
pheno.medians <- my.phenoResult %>% dplyr::bind_rows() %>% dplyr::select(-Mah.cluster) %>% 
  group_by(patient, condition, stimulation, my.clusters) %>% 
  summarize_all(median)

message("Obtaining Phenograph Metaclusters")
#run phenograph on clusters from each patient to yield metaclusters
my.meta.phenograph <- pheno.medians %>% dplyr::ungroup() %>% 
                      dplyr::select(which(colnames(.) %in% SURFACE.MARKERS)) %>% 
  Rphenograph(k = 15)
my.metaclusters <- membership(my.meta.phenograph[[2]]) %>% as.vector()
pheno.medians <- cbind(pheno.medians, my.metaclusters = my.metaclusters)
metacluster.table <- pheno.medians %>% dplyr::select(patient, condition, stimulation, 
                                                     my.clusters = my.clusters, metacluster = my.metaclusters) %>% 
  group_split()

#adding the phenograph metacluster assignment to expression.matrix
expression.matrix$phenograph.metacluster <- rep(0, times = nrow(expression.matrix)) #remove me? 

my.clust <- makeCluster(num_cores, outfile = "")
registerDoParallel(my.clust)

expression.matrix <- foreach(lookup.table=metacluster.table, fcs.df = my.phenoResult, 
                             .packages = "dplyr") %dopar% {
  dictionary <- lookup.table %>% dplyr::select(my.clusters, metacluster)
  fcs.df <- left_join(x = fcs.df, y = dictionary, by = "my.clusters")
  fcs.df
} %>% dplyr::bind_rows()

stopCluster(my.clust)

#phenograph cluster altogether 
#giant.phenograph <- Rphenograph(data = expression.matrix[,which(surface.cols)], k = 15)
#expression.matrix$phenograph.cluster <- membership(giant.phenograph[[2]]) %>% as.vector()

expression.matrix$phenograph.metacluster <- expression.matrix$metacluster
expression.matrix$metacluster <- NULL
expression.matrix$my.clusters <- NULL
#save data structure
setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("expression.matrix"), 
     file = "AML_matrix_Phenograph.RData")
setwd(CODE.DIRECTORY)

rm(df.list, metacluster.table, my.clust, my.meta.phenograph, my.phenoResult, pheno.medians)

#############     FlowSOM    ###############
message("FlowSOM clustering...") 
expression.matrix$patient <- as.factor(expression.matrix$patient)
expression.matrix$stimulation <- as.factor(expression.matrix$stimulation)
expression.matrix$condition <- as.factor(expression.matrix$condition)
my.FSO <- list(data = data.matrix(expression.matrix), 
               compensate = FALSE, 
               spillover = NULL, 
               transform = FALSE, 
               scale = NULL, 
               prettyColnames = colnames(expression.matrix))
my.SOM <- BuildSOM(fsom = my.FSO, colsToUse = which(surface.cols))

#are there SOM arguments that I want to pass to ... in the function above?
my.MST <- BuildMST(my.SOM, tSNE = FALSE)

# PlotStars(my.MST, backgroundValues = as.factor(my.MST$map$medianValues[,"condition"]))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.MST$map$medianValues[,"condition"]))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.MST$map$medianValues[,"patient"]))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.MST$map$medianValues[,"patient"]))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.MST$map$medianValues[,"phenograph.cluster"]))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.MST$map$medianValues[,"phenograph.cluster"]))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

expression.matrix$FlowSOM.cluster <- my.MST$map$mapping[,1]

# my.flowSOM.modes <- expression.matrix %>% dplyr::group_by(FlowSOM.cluster) %>%
#   summarize_all(Mode) 
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$condition))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$condition))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$patient))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$patient))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$phenograph.cluster))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$phenograph.cluster))
# 
# PlotStars(my.MST, backgroundValues = as.factor(my.flowSOM.modes$Mah.cluster))
# PlotStars(my.MST, view = "tSNE", backgroundValues = as.factor(my.flowSOM.modes$Mah.cluster))

message("Finding FlowSOM metaclusters...")
flowSOM.metaclusters <- MetaClustering(data = my.MST$map$codes, 
                                       method = "metaClustering_consensus", 
                                       max = 50)
flowSOM.metaclusters.cells <- flowSOM.metaclusters[my.MST$map$mapping[,1]]

expression.matrix$FlowSOM.metacluster <- flowSOM.metaclusters.cells

setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("expression.matrix"), file = "AML_matrix_FlowSOM.RData")

rm(my.FSO, my.MST, my.SOM, flowSOM.metaclusters, flowSOM.metaclusters.cells)

#################################       k-means       ###################################
message("K-means clustering...")

my.kmeans <- kmeans(x = expression.matrix[,which(surface.cols)], centers = 200, iter.max = 200,  nstart = 1)
kmeans.clusters <- my.kmeans$cluster
expression.matrix$kmeans.cluster <- kmeans.clusters

message("Saving final expression matrix in ", OUTPUT.DIRECTORY)
setwd(file.path(OUTPUT.DIRECTORY, "data"))
save(list = c("expression.matrix"), file = "AML_matrix_clustered.RData")

rm(my.clust, my.kmeans,
   kmeans.clusters, my.metaclusters, expression.matrix)

message("Clustering analysis completed!")

