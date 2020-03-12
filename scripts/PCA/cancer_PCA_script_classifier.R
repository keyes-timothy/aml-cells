#####################  PCA analysis of AML data   ###########################

#Set-up dataset
setwd(file.path(OUTPUT.DIRECTORY, "data"))
load("cancer_data_classified.RData")

#PCA with cancer cells
message("Plotting PCA for cancer samples...")
cancer.PCA <- prcomp(cancer.data.scaled, scale = FALSE) 
cancer.loadings <- cancer.PCA$rotation 
cancer.PCs <- cancer.PCA$x
cancer.PC.importance <- summary(cancer.PCA)$importance

######Plot PCAs

#set general themes
axis.label.style <- element_text(face = "bold", size = 24)
axis.style <- element_text(size = 14)
PCA.theme <- theme(axis.title.x = axis.label.style, axis.title.y = axis.label.style, 
                   axis.text = axis.style, legend.title = element_text(face = "bold", size = 16), 
                   legend.text = element_text(size = 14)) + theme_bw()
min.cell.number.cancer <- cbind(cancer.PCs, cancer.demographics) %>% group_by(patient) %>% count(patient) %>% ungroup() %>% 
                          select(n) %>% min()
cancer.demographics$CD117 <- cancer.data$CD117

#calculate classified centroids
cancer.PCA.centroids <- cbind(cancer.PCs, cancer.demographics) %>% group_by(MahID) %>% 
  dplyr::select(contains("PC")) %>% 
  summarize_all(.funs = mean)

#Actual plots
message("Plotting cancer.PCA.plot")
cancer.PCA.plot <- cbind(cancer.PCs, cancer.demographics) %>% 
  group_by(patient, condition, MahID) %>% 
  sample_n(size = (300), replace = TRUE)
cancer.PCA.plot$patient <- as.factor(cancer.PCA.plot$patient)

setwd(file.path(OUTPUT.DIRECTORY, "data"))
cancer.PCA.df_classified <- cancer.PCA.plot
save(list = c("cancer.PCA.df_classified"), file = "cancer_PCA_df_classified.RData")
rm(cancer.PCA.df_classified)

cancer.plot <- ggplot(data = cancer.PCA.plot, 
                      mapping = aes(x=PC1, y=PC2, color = patient)) + 
  #geom_point(size = 1, alpha = 0.4)
  geom_density_2d()
cancer.plot <- cancer.plot + 
               xlab(paste0("PC1 (", round(cancer.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
               ylab(paste0("PC2 (", round(cancer.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
               facet_wrap(patient~condition) + PCA.theme + 
  theme(legend.position = "none")
ggsave(file = "cancer_PCA_plot.pdf", plot = cancer.plot, device = "pdf", 
       path = PCA.output, height = 10, width = 12, units = "in")

#classified populations
message("Plotting classified populations")
cancer.PCA.plot <- cbind(cancer.PCs, cancer.demographics) %>% as_tibble()
cancer.PCA.plot$patient <- as.factor(cancer.PCA.plot$patient)
cancer.PCA.plot$MahID <- as.factor(cancer.PCA.plot$MahID)
levels(cancer.PCA.plot$MahID) <- c("Unclassified", CLASSIFIER.POPULATIONS)
cancer.PCA.centroids$MahID <- as.factor(cancer.PCA.centroids$MahID)
levels(cancer.PCA.centroids$MahID) <- c("Unclassified", CLASSIFIER.POPULATIONS)
cancer.PCA.plot <- cancer.PCA.plot %>% dplyr::filter(MahID != "Unclassified") %>% 
  group_by(MahID) %>% sample_n(size = 1000, replace = TRUE)

#faceted populations by patient and stimulation
cancer.plot.pops <- ggplot(data = cancer.PCA.plot, 
                           mapping = aes(x=PC1, y=PC2, color = MahID)) + 
  #geom_point(size = 1, alpha = 0.4) 
  geom_density_2d()

cancer.plot.pops <- cancer.plot.pops + 
  xlab(paste0("PC1 (", round(cancer.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(cancer.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  facet_wrap(patient~condition) + PCA.theme + 
  guides(color=guide_legend(title="Subpopulation", override.aes = list(size = 3, alpha = 1)))
ggsave(file = "cancer_PCA_plot_pops.pdf", plot = cancer.plot.pops, 
       device = "pdf", path = PCA.output, height = 10, width = 14, units = "in")

#not faceted
cancer.plot.pops2 <- ggplot(data = cancer.PCA.plot, 
                           mapping = aes(x=PC1, y=PC2, color = MahID)) + 
  #geom_point(size = 1, alpha = 0.4) + 
  geom_density_2d(alpha = 0.6) + 
  geom_point(data = cancer.PCA.centroids, 
             mapping = aes(x = PC1, y = PC2), color = "black", size = 3) + 
  geom_text_repel(data = cancer.PCA.centroids, 
            mapping = aes(x = PC1, y = PC2, label = MahID), 
            color = "black", fontface = "bold", nudge_y = 0.2)

cancer.plot.pops2 <- cancer.plot.pops2 + 
  xlab(paste0("PC1 (", round(cancer.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(cancer.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme + 
  guides(color=guide_legend(title="Subpopulation", 
                            override.aes = list(size = 3, alpha = 1)))

ggsave(file = "cancer_PCA_plot_pops2.pdf", plot = cancer.plot.pops2, 
       device = "pdf", path = PCA.output, height = 7, width = 8, units = "in")

#Express cancer using healthy coordinate system - healthy trajectory is in black
message("Plotting cancer samples under the healthy coordinate system...")
cancer.under.healthy.PCs <- as.matrix(cancer.data.scaled) %*% as.matrix(healthy.loadings) 
cancer.under.healthy.plot <- cbind(cancer.under.healthy.PCs, cancer.demographics) %>% 
  as_tibble() 
cancer.under.healthy.plot$MahID <- as.factor(cancer.under.healthy.plot$MahID)
levels(cancer.under.healthy.plot$MahID) <- c("Unclassified", CLASSIFIER.POPULATIONS)

cancer.under.healthy.plot$Subpopulation <- cancer.under.healthy.plot$MahID
healthy.PCA.plot$Subpopulation <- as.factor(healthy.PCA.plot$MahID)
levels(healthy.PCA.plot$MahID) <- c("Unclassified",  CLASSIFIER.POPULATIONS)
healthy.PCA.plot$CD34 <- NULL
healthy.PCA.plot$CD38 <- NULL
healthy.PCA.plot <- ungroup(healthy.PCA.plot)
cancer.under.healthy.plot2 <- rbind(cancer.under.healthy.plot, healthy.PCA.plot) %>% dplyr::filter(!(MahID == "Unclassified")) %>%
  group_by(patient, MahID) %>% sample_n(size = 100, replace = TRUE) %>% as_tibble()

cancer.u.healthy.plot <- ggplot(data = cancer.under.healthy.plot2, 
                                mapping = aes(x=PC1, y=PC2, color = Subpopulation)) +
  PCA.theme + 
  geom_point(size = 1, alpha = 0.4) + facet_wrap(patient~condition) + 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))



# print("Adding healthy panels to cancer plots...")
# healthy.PCA.plot <- sample_n(healthy.PCA.plot, size = 1000, replace = TRUE)
# cancer.u.healthy.plot <- cancer.u.healthy.plot +
#                          geom_point(data = healthy.PCA.plot, 
#                                      color = healthy.PCA.plot$Subpopulation, alpha = 0.5, size = 0.5) + 
#   PCA.theme
 ggsave(file = "cancer_under_healthy_plot.pdf", 
        plot = cancer.u.healthy.plot, device = "pdf", 
        path = PCA.output, 
        height = 12, width = 18, units = "in")
 
#not faceted 
 cancer.under.healthy.plot3 <- rbind(cancer.under.healthy.plot, healthy.PCA.plot) %>% dplyr::filter(!(MahID == "Unclassified")) %>%
   group_by(MahID) %>% sample_n(size = 500, replace = TRUE) %>% unique() %>% as_tibble()
 
 cancer.u.healthy.plot3 <- ggplot(data = cancer.under.healthy.plot3, 
                                 mapping = aes(x=PC1, y=PC2, color = Subpopulation)) +
   PCA.theme + 
   geom_point(size = 3, alpha = 0.5) + 
   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
 
 ggsave(file = "cancer_under_healthy_plot3.pdf", 
        plot = cancer.u.healthy.plot3, device = "pdf", 
        path = PCA.output, 
        height = 8, width = 10, units = "in")
 

#quantitative representation of which markers are driving each PC in cancer
message("Making cancer biplots...")
library(ggfortify)
#same in cancer 
cancer.PCA$x <- as.data.frame(cancer.PCA$x) %>% sample_n(15000, replace = TRUE) %>% as.matrix()
cancer.most.PC1 <- as.data.frame(cancer.PCA$rotation[,"PC1"]) %>% abs() %>% as.matrix()
cancer.most.PC1 <- rownames(cancer.most.PC1)[order(cancer.most.PC1, decreasing = TRUE)]
cancer.most.PC1 <- cancer.most.PC1[1:10]
cancer.most.PC2 <- as.data.frame(cancer.PCA$rotation[,"PC2"]) %>% abs() %>% as.matrix()
cancer.most.PC2 <- rownames(cancer.most.PC2)[order(cancer.most.PC2, decreasing = TRUE)]
cancer.most.PC2 <- cancer.most.PC2[1:10]
cancer.scree.PCs <- unique(c(cancer.most.PC1, cancer.most.PC2))
cancer.PCA$rotation <- cancer.PCA$rotation[cancer.scree.PCs,]
cancer.PCA$scale <- cancer.PCA$scale[cancer.scree.PCs]
cancer.PCA$center <- cancer.PCA$center[cancer.scree.PCs]
cancer.PC.drivers <- autoplot(cancer.PCA, loadings = TRUE, colour = "red", loadings.colour = 'blue', 
                              loadings.label = TRUE, loadings.label.fontface = "bold", 
                              loadings.label.size = 7, loadings.label.colour = "black", alpha = 0.2, 
                              loadings.label.repel = TRUE) + PCA.theme
ggsave(file = "Cancer_PCA_drivers.tif", plot = cancer.PC.drivers, device = "pdf", path = PCA.output, 
       width = 10, height = 10, units = "in")


message("Plotting 3D Cancer plots...")
library(plot3D)
setwd(PCA.output)
# 3D plotting - cancer
cancer.under.healthy.plot <- cancer.under.healthy.plot %>% ungroup() %>% 
  group_by(MahID) %>% sample_n(size = 500, replace = TRUE) 
for (i in seq(from = 0, to = 360, by = 10)){ 
  tiff(filename = paste0("cancer3D_", i, ".tiff"), width = 10, height = 10, units = "in", res = 300)
  cancer.u.healthy.3D <- scatter3D(x = cancer.under.healthy.plot$PC1, 
                                    y = cancer.under.healthy.plot$PC2, 
                                    z = cancer.under.healthy.plot$PC3, 
                                   colkey = list(labels = CLASSIFIER.POPULATIONS, 
                                                 at = c(1, 2, 3, 4, 5, 6, 7, 8,9), 
                                                 dist = -0.05),
                                   colvar = as.numeric(cancer.under.healthy.plot$MahID),
                                   xlab = "PC1", 
                                   ylab = "PC2", 
                                   zlab = "PC3", 
                                   pch = 19, cex = 0.5, phi = 30, theta = i)
    cancer.u.healthy.3D
    dev.off()
}

#remove variables that are not needed
rm(list = c("cancer.data", "cancer.data.scaled", "cancer.demographics", 
            "cancer.loadings", "cancer.PC.drivers", "cancer.PC.importance", 
            "cancer.PCA", "cancer.PCA.plot", "cancer.PCs", "cancer.plot", 
            "cancer.u.healthy.3D", "cancer.u.healthy.plot", "cancer.under.healthy.PCs", 
            "cancer.under.healthy.plot", "cancer.most.PC1", "cancer.most.PC2", 
            "cancer.scree.PCs", "healthy.PCA.plot", "healthy.loadings", 
            "cancer.PCA.centroids", "cancer.plot.pops", "cancer.plot.pops2", 
            "cancer.under.healthy.plot2"))

message("Cancer PCA analysis completed!")

