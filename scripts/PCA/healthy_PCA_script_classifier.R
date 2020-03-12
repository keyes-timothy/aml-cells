#####################  PCA analysis of AML data   ###########################

#Set-up dataset
test.matrix <- expression.matrix
rm(expression.matrix)

#create directory for plots
dir.create(PCA.output, showWarnings = FALSE)

#separate healthy and cancer data
message("Setting up data...") 
healthy.demographics <- dplyr::filter(test.matrix, condition=="Healthy") %>% 
                        select(condition, stimulation, patient, MahID)
healthy.data <- dplyr::filter(test.matrix, condition=="Healthy") %>% 
                        select(-condition, -stimulation, -patient, -cisplatin) %>% 
                        select(-contains("Mah"))
cancer.demographics <- dplyr::filter(test.matrix, condition!="Healthy") %>% 
                       select(condition, stimulation, patient, MahID)
cancer.data <- dplyr::filter(test.matrix, condition!="Healthy") %>% 
               select(-condition, -stimulation, -patient, -cisplatin) %>% 
               select(-contains("Mah"))
healthy.data.scaled <- scale(healthy.data, center = TRUE, scale = TRUE)
cancer.data.scaled <- scale(cancer.data, center = TRUE, scale = TRUE)

#perform healthy PCA and extract healthy values and healthy coordinate system
message("Plotting PCA for healthy samples...")
healthy.PCA <- prcomp(healthy.data.scaled, scale = FALSE)
healthy.loadings <- healthy.PCA$rotation
healthy.PCs <- healthy.PCA$x
healthy.PC.importance <- summary(healthy.PCA)$importance

######plot PCAs

#set general themes
axis.label.style <- element_text(face = "bold", size = 24)
axis.style <- element_text(size = 14)
PCA.theme <- theme(axis.title.x = axis.label.style, axis.title.y = axis.label.style, 
                   axis.text = axis.style, legend.title = element_text(face = "bold", size = 16), 
                   legend.text = element_text(size = 14)) + theme_bw()

#set up healthy PCA ggplot dataframe
healthy.PCA.plot <- cbind(healthy.PCs, healthy.demographics, healthy.data.scaled)

min.cell.number <- healthy.PCA.plot %>% count(patient) %>% ungroup() %>% select(n) %>% min()
message(paste0("The minimum cell number in all patients is ", min.cell.number))

healthy.PCA.centroids <- healthy.PCA.plot %>% group_by(MahID) %>% 
  dplyr::select(contains("PC")) %>% 
  summarize_all(.funs = mean)


#actual plots
healthy.PCA.plot <- healthy.PCA.plot %>% group_by(patient, MahID) %>% sample_n(size = 1000, replace = FALSE)
setwd(PCA.output)
healthy.PCA.df <- healthy.PCA.plot
save(list = c("healthy.PCA.df"), file = "healthy_PCA_df_classified.RData")
rm(healthy.PCA.df)

healthy.plot.patient <- ggplot(data = healthy.PCA.plot %>% group_by(patient) %>% sample_n(size = 1000), mapping = aes(x = PC1, y = PC2, color = patient)) + 
  geom_point(alpha = 0.5)
healthy.plot.patient <- healthy.plot.patient + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
healthy.plot.patient
ggsave(file = "Patient_healthy_PCA.pdf", plot = healthy.plot.patient, device = "pdf", path = PCA.output, 
       height = 7, width = 8, units = "in")

healthy.plot.CD34 <- ggplot(data = healthy.PCA.plot %>% group_by(patient) %>% sample_n(size = 1000), mapping = aes(x = PC1, y = PC2, color = CD34)) + 
  geom_point(alpha = 0.5)
healthy.plot.CD34 <- healthy.plot.CD34 + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
healthy.plot.CD34
ggsave(file = "CD34_healthy_PCA.pdf", plot = healthy.plot.CD34, device = "pdf", path = PCA.output, 
       height = 7, width = 8, units = "in")


healthy.PCA.plot$MahID <- as.factor(healthy.PCA.plot$MahID)
healthy.PCA.centroids$MahID <- as.factor(healthy.PCA.centroids$MahID)
levels(healthy.PCA.plot$MahID) <- CLASSIFIER.POPULATIONS
levels(healthy.PCA.centroids$MahID) <- CLASSIFIER.POPULATIONS
healthy.plot.populations <- healthy.PCA.plot %>% ungroup() %>% 
  group_by(MahID) %>% sample_n(500, replace = FALSE) %>% ungroup()

healthy.plot.pops <- ggplot(data = healthy.plot.populations, mapping = aes(x = PC1, y = PC2, color = MahID)) + 
  geom_point() + 
  geom_point(data = healthy.PCA.centroids, mapping = aes(x = PC1, y = PC2), color = "black", size = 3) + 
  geom_text_repel(data = healthy.PCA.centroids, 
            mapping = aes(x = PC1, y = PC2, label = MahID), 
            color = "black", fontface = "bold", nudge_y = 0.2)#, position = position_jitter(width = 0.1, height = 0.2))
healthy.plot.pops <- healthy.plot.pops + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme + guides(color=guide_legend(title="Subpopulation", override.aes = list(size = 3, alpha = 1)))
healthy.plot.pops

ggsave(file = "populations_healthy_PCA.pdf", plot = healthy.plot.pops, 
       device = "pdf", path = PCA.output, 
       height = 7, width = 8, units = "in")

# healthy.plot.CD117 <- ggplot(data = healthy.PCA.plot, mapping = aes(x = PC1, y = PC2, color = CD117)) + geom_point()
# healthy.plot.CD117 <- healthy.plot.CD117 + xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + ylab(paste("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)"))
# ggsave(file = "CD117_healthy_PCA.tif", plot = healthy.plot.CD117, 
#        device = "pdf", path = PCA.output)

healthy.plot.all <- ggplot(data = healthy.PCA.plot, 
                           mapping = aes(x = PC1, y = PC2, color = patient)) + 
  #geom_point(alpha = 0.5) + 
  geom_density_2d() + 
  facet_wrap(~patient) + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme + 
  theme(legend.position = "none")
healthy.plot.all
ggsave(file = "healthy_PCA_patient.pdf", plot = healthy.plot.all, device = "pdf", path = PCA.output, 
       height = 8, width = 10, units = "in")

#umap 
healthy.umap.all <- healthy.PCA.plot %>% ungroup() %>% dplyr::select(contains(match = "PC", ignore.case = FALSE)) %>% umap()
umap.coordinates <- healthy.umap.all$layout %>% as_tibble(); colnames(umap.coordinates) <- c("UMAP1", "UMAP2")
healthy.PCA.plot <- cbind(healthy.PCA.plot %>% ungroup(), umap.coordinates) %>% as_tibble()

healthy.plot.all.umap <- ggplot(data = healthy.PCA.plot, 
                           mapping = aes(x = UMAP1, y = UMAP2, color = patient)) + 
  #geom_point(alpha = 0.5) + 
  geom_density_2d() + 
  facet_wrap(~patient) + 
  PCA.theme + 
  theme(legend.position = "none")
healthy.plot.all.umap
ggsave(file = "healthy_umap_patient.pdf", plot = healthy.plot.all.umap, device = "pdf", path = PCA.output, 
       height = 8, width = 10, units = "in")

#classified populations
healthy.PCA.plot$MahID <- as.factor(healthy.PCA.plot$MahID)
levels(healthy.PCA.plot$MahID) <- CLASSIFIER.POPULATIONS
healthy.plot.all2 <- ggplot(data = healthy.PCA.plot, 
                           mapping = aes(x = PC1, y = PC2, color = MahID)) + 
  #geom_point(alpha = 0.5) + 
  geom_density_2d(alpha = 0.7) + 
  facet_wrap(~patient) + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
healthy.plot.all2
ggsave(file = "healthy_PCA_classified_pops.pdf", plot = healthy.plot.all2, device = "pdf", 
       path = PCA.output, height = 8, width = 10, units = "in")

#umap
healthy.plot.all.umap2 <- ggplot(data = healthy.PCA.plot, 
                            mapping = aes(x = UMAP1, y = UMAP2, color = MahID)) + 
  geom_point(alpha = 0.5) + 
  #geom_density_2d(alpha = 0.7) + 
  facet_wrap(~patient) + 
  PCA.theme
healthy.plot.all.umap2
ggsave(file = "healthy_umap_classified_pops.pdf", plot = healthy.plot.all.umap2, device = "pdf", 
       path = PCA.output, height = 8, width = 10, units = "in")

#sampled populations
healthy.PCA.plot2 <- healthy.PCA.plot %>% group_by(patient, MahID) %>% sample_n(250, replace = TRUE) #change this value
healthy.plot.all3 <- ggplot(data = healthy.PCA.plot2, 
                            mapping = aes(x = PC1, y = PC2, color = MahID)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~patient) + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme + guides(color=guide_legend(title="Subpopulation", override.aes = list(size = 3, alpha = 1)))
healthy.plot.all3

ggsave(file = "healthy_PCA_classified_pops_sampled.pdf", plot = healthy.plot.all3, device = "pdf",
       path = PCA.output, height = 8, width = 10, units = "in")


#######Which markers are driving the healthy PCA? 
message("Analyzing PCA results and quality...")
#general quality of PCA
healthy.PC.proportions.explained <- healthy.PC.importance[2,1:9]
healthy.PC.proportions.explained.plot <- ggplot(mapping = aes(x = as.factor(names(healthy.PC.proportions.explained)), y = healthy.PC.proportions.explained)) +  
                                                geom_bar(stat = "identity") + xlab("Components") + ylab("Proportion Variance Explained")
ggsave(file = "Healthy_Screeplot.pdf", plot = healthy.PC.proportions.explained.plot, device = "pdf", 
       path = PCA.output)


#Markers with highest relevance to PC1 and PC2 
healthy.PC1.rotation <- healthy.PCA$rotation[,1] %>% abs() %>% sort(decreasing = TRUE)
healthy.PC2.rotation <- healthy.PCA$rotation[,2] %>% abs() %>% sort(decreasing = TRUE)

#PC1
healthy.PC1.markers.plot <- ggplot(data = NULL, mapping = aes(x = factor(names(healthy.PC1.rotation), levels = names(healthy.PC1.rotation)), y = healthy.PC1.rotation)) + 
      geom_bar(stat = "identity") + xlab("Marker") + ylab("PC1 Loading") + theme(axis.text.x  = element_text(angle = 60, hjust = 1, size = 11))
ggsave(file = "healthy.PC1.pdf", plot = healthy.PC1.markers.plot, device = "pdf", path = PCA.output, width = 18, height = 8, units = "cm")

#PC2
healthy.PC2.markers.plot <- ggplot(data = NULL, mapping = aes(x = factor(names(healthy.PC2.rotation), levels = names(healthy.PC2.rotation)), y = healthy.PC2.rotation)) + 
  geom_bar(stat = "identity") + xlab("Marker") + ylab("PC2 Loading") + theme(axis.text.x  = element_text(angle = 60, hjust = 1, size = 11))
ggsave(file = "healthy.PC2.pdf", plot = healthy.PC2.markers.plot, device = "pdf", path = PCA.output, width = 18, height = 8, units = "cm")

#Markers with highest variance in healthy and cancer
healthy.variances <- healthy.data %>% 
  summarize_all(.funs = var) %>% as.numeric()
names(healthy.variances) <- colnames(healthy.data)
  healthy.variances <- healthy.variances %>% sort(decreasing = TRUE)

healthy.variances.plot <- ggplot(data = NULL, mapping = aes(x = factor(names(healthy.variances), levels = names(healthy.variances)), y = as.numeric(healthy.variances))) + 
  geom_bar(stat = "identity") + xlab("Marker") + ylab("Variance") + theme(axis.text.x  = element_text(angle = 60, hjust = 1, size = 11))
ggsave(file = "healthy.variances.pdf", plot = healthy.variances.plot, device = "pdf", path = PCA.output, width = 18, height = 8, units = "cm")

cancer.variances <- cancer.data %>%
  summarize_all(.funs = var) %>% as.numeric()
names(cancer.variances) <- colnames(cancer.data)
cancer.variances <- cancer.variances %>% sort(decreasing = TRUE)

cancer.variances.plot <- ggplot(data = NULL, mapping = aes(x = factor(names(cancer.variances), levels = names(cancer.variances)), y = as.numeric(cancer.variances))) + 
  geom_bar(stat = "identity") + xlab("Marker") + ylab("Variance") + theme(axis.text.x  = element_text(angle = 60, hjust = 1, size = 11))
ggsave(file = "cancer.variances.pdf", plot = cancer.variances.plot, device = "pdf", path = PCA.output, width = 18, height = 8, units = "cm")

#quantitative representation of which markers are driving each PC in healthy
message("Making healthy biplots...")
library(ggfortify)
healthy.PCA$x <- as.data.frame(healthy.PCA$x) %>% sample_n(7500, replace = TRUE) %>% as.matrix()
healthy.most.PC1 <- as.data.frame(healthy.PCA$rotation[,"PC1"]) %>% abs() %>% as.matrix()
healthy.most.PC1 <- rownames(healthy.most.PC1)[order(healthy.most.PC1, decreasing = TRUE)]
healthy.most.PC1 <- healthy.most.PC1[1:10]
healthy.most.PC2 <- as.data.frame(healthy.PCA$rotation[,"PC2"]) %>% abs() %>% as.matrix()
healthy.most.PC2 <- rownames(healthy.most.PC2)[order(healthy.most.PC2, decreasing = TRUE)]
healthy.most.PC2 <- healthy.most.PC2[1:10]
healthy.scree.PCs <- unique(c(healthy.most.PC1, healthy.most.PC2))
healthy.PCA$rotation <- healthy.PCA$rotation[healthy.scree.PCs,]
healthy.PCA$scale <- healthy.PCA$scale[healthy.scree.PCs]
healthy.PCA$center <- healthy.PCA$center[healthy.scree.PCs]
healthy.PC.drivers <- autoplot(healthy.PCA, loadings = TRUE, colour = "red", loadings.colour = 'blue', alpha = 0.2,
                               loadings.label = TRUE, loadings.label.fontface = "bold", 
                               loadings.label.size = 7, loadings.label.colour = "black") + 
  PCA.theme
ggsave(file = "Healthy_PCA_drivers.pdf", plot = healthy.PC.drivers, device = "pdf", path = PCA.output, 
       width = 10, height = 10, units = "in") 

message("Performing 3D plotting for healthy samples...")
# 3D plotting - healthy
library(plot3D)
setwd(PCA.output)
for (i in seq(from = 0, to = 360, by = 10)){ 
    pdf(file = paste0("healthy3D_CD34_", i, ".pdf"), width = 10, height = 10)
    healthy.3D <- scatter3D(x = healthy.PCA.plot$PC1, y = healthy.PCA.plot$PC2, 
                            z = healthy.PCA.plot$PC3, 
                            colvar = healthy.PCA.plot$CD34,
                            xlab = "PC1", 
                            ylab = "PC2", 
                            zlab = "PC3", 
                            pch = 19, cex = 0.5, phi = 30, theta = i)
    healthy.3D
    dev.off()
}

#classified in 3D
for (i in seq(from = 0, to = 360, by = 10)){ 
  pdf(file = paste0("pops_healthy3D_", i, ".pdf"), width = 15, height = 10)
  healthy.3D <- scatter3D(x = healthy.PCA.plot$PC1, y = healthy.PCA.plot$PC2, 
                          z = healthy.PCA.plot$PC3, 
                          colkey = list(labels = CLASSIFIER.POPULATIONS, 
                                        at = c(1, 2, 3, 4, 5, 6, 7, 8,9), 
                                        dist = -0.1),
                          colvar = as.numeric(healthy.PCA.plot$MahID),
                          xlab = "PC1", 
                          ylab = "PC2", 
                          zlab = "PC3", 
                          pch = 19, cex = 0.5, phi = 30, theta = i)
  healthy.3D
  dev.off()
}

setwd(file.path(OUTPUT.DIRECTORY, "data"))

message("Saving important data files...") 
save(list = c("healthy.demographics", "healthy.data", "healthy.data.scaled"), file = "healthy_data_classified.RData")
save(list = c("cancer.demographics", "cancer.data", "cancer.data.scaled", "healthy.loadings", "healthy.PCA.plot"), file = "cancer_data_classified.RData")

rm(list = c("healthy.demographics", "healthy.data", "healthy.data.scaled",
            "cancer.demographics", "cancer.data", "cancer.data.scaled", 
            "healthy.loadings", "healthy.PCA.plot", "healthy.variances",
            "healthy.most.PC1", "healthy.most.PC2", "healthy.PC.drivers", 
            "healthy.PC1.rotation", "healthy.PC2.rotation", "min.cell.number", 
            "healthy.PC.proportions.explained", "healthy.PC.proportions.explained.plot",
            "healthy.3D", "healthy.plot.all", 
            "healthy.plot.CD34", "healthy.plot.patient", 
            "test.matrix", "healthy.scree.PCs", "cancer.variances", "healthy.PC.importance", 
            "healthy.PC1.markers.plot", "healthy.PC2.markers.plot", "healthy.PCA", 
            "healthy.PCs", "healthy.variances.plot", "cancer.variances.plot", 
            "healthy.PCA.centroids", "healthy.PCA.plot2", "healthy.plot.all2", 
            "healthy.plot.all3", "healthy.plot.pops", "healthy.plot.populations"))

setwd(CODE.DIRECTORY)
message("Healthy PCA analysis completed!")



# plot_ly(healthy.PCA.plot %>% group_by(MahID) %>% sample_n(size = 250), 
#         x = ~PC1, y = ~PC2, z = ~PC3, color = ~MahID, 
#         marker = list(size = 4)) %>%
#   add_markers() %>%
#   layout(scene = list(xaxis = list(title = 'PC1'),
#                       yaxis = list(title = 'PC2'),
#                       zaxis = list(title = 'PC3')))
