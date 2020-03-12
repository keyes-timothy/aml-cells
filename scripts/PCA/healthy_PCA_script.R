#####################  PCA analysis of AML data   ###########################

#Set-up dataset
test.matrix <- expression.matrix
rm(expression.matrix)

#create directory for plots
dir.create(PCA.output, showWarnings = FALSE)

#separate healthy and cancer data
message("Setting up data...") 
healthy.demographics <- dplyr::filter(test.matrix, condition=="Healthy") %>% select(condition, stimulation, patient)
healthy.data <- dplyr::filter(test.matrix, condition=="Healthy") %>% select(-condition, -stimulation, -patient, -cisplatin)
cancer.demographics <- dplyr::filter(test.matrix, condition!="Healthy") %>% select(condition, stimulation, patient)
cancer.data <- dplyr::filter(test.matrix, condition!="Healthy") %>% select(-condition, -stimulation, -patient, -cisplatin)
healthy.data.scaled <- scale(healthy.data, center = TRUE, scale = TRUE)
cancer.data.scaled <- scale(cancer.data, center = TRUE, scale = TRUE)

#perform healthy PCA, extract healthy values and healthy coordinate system
message("Plotting PCA for healthy samples...")
healthy.PCA <- prcomp(healthy.data.scaled, scale = FALSE)
healthy.loadings <- healthy.PCA$rotation
healthy.PCs <- healthy.PCA$x
healthy.PC.importance <- summary(healthy.PCA)$importance

#plot healthy
axis.label.style <- element_text(face = "bold", size = 24)
axis.style <- element_text(size = 14)
PCA.theme <- theme(axis.title.x = axis.label.style, axis.title.y = axis.label.style, 
                   axis.text = axis.style, legend.title = element_text(face = "bold", size = 16), 
                   legend.text = element_text(size = 14)) + theme_bw()

healthy.PCA.plot <- cbind(healthy.PCs, healthy.demographics)
healthy.PCA.plot$CD34 <- healthy.data$CD34
healthy.PCA.plot$CD38 <- healthy.data$CD38
healthy.PCA.plot$CD117 <- healthy.data$CD117


#saving healthy PCA data for plotting locally 
message("Sampling healthy PCA")
setwd(PCA.output)
healthy.PCA.df <- healthy.PCA.plot %>% group_by(patient, stimulation, condition) %>% 
  sample_n(size = 2000, replace = TRUE)
save(list = c("healthy.PCA.df"), file = "sampled_healthy_PCA_df.RData")
rm(healthy.PCA.df)

min.cell.number <- healthy.PCA.plot %>% count(patient) %>% ungroup() %>% select(n) %>% min()
healthy.PCA.plot <- healthy.PCA.plot %>% group_by(patient) %>% sample_n(size =2000) 

healthy.plot.patient <- ggplot(data = healthy.PCA.plot, mapping = aes(x = PC1, y = PC2, color = patient)) + 
  geom_point(alpha = 0.5)
healthy.plot.patient <- healthy.plot.patient + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
ggsave(file = "Patient_healthy_PCA.pdf", plot = healthy.plot.patient, device = "pdf", path = PCA.output, height = 7, width = 8, units = "in")

healthy.plot.CD34 <- ggplot(data = healthy.PCA.plot, mapping = aes(x = PC1, y = PC2, color = CD34)) + geom_point()
healthy.plot.CD34 <- healthy.plot.CD34 + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
ggsave(file = "CD34_healthy_PCA.pdf", plot = healthy.plot.CD34, device = "pdf", path = PCA.output, 
       height = 7, width = 8, units = "in")

healthy.plot.CD38 <- ggplot(data = healthy.PCA.plot, mapping = aes(x = PC1, y = PC2, color = CD38)) + geom_point()
healthy.plot.CD38 <- healthy.plot.CD38 + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
ggsave(file = "CD38_healthy_PCA.pdf", plot = healthy.plot.CD38, 
       device = "pdf", path = PCA.output, height = 7, width = 8, units = "in")

healthy.plot.CD117 <- ggplot(data = healthy.PCA.plot, mapping = aes(x = PC1, y = PC2, color = CD117)) + geom_point()
healthy.plot.CD117 <- healthy.plot.CD117 + 
  xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
  ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
  PCA.theme
ggsave(file = "CD117_healthy_PCA.pdf", plot = healthy.plot.CD117, 
       device = "pdf", path = PCA.output, height = 7, width = 8, units = "in")

healthy.plot.all <- ggplot(data = healthy.PCA.plot, 
                           mapping = aes(x = PC1, y = PC2, color = patient)) + 
                    geom_point(alpha = 0.5) + facet_wrap(~patient) + theme(legend.position = "none", 
                                                                           axis.title = element_text(face = "bold", size = 17), 
                                                                           axis.text = element_text(size = 12)) + 
                    xlab(paste0("PC1 (", round(healthy.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
                    ylab(paste0("PC2 (", round(healthy.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) 
ggsave(file = "healthy_PCA_all.pdf", plot = healthy.plot.all, device = "pdf", path = PCA.output, height = 6, width = 6, units = "in")


#######Which markers are driving the healthy PCA? 
message("Analyzing PCA results and quality...")
#general quality of PCA
healthy.PC.proportions.explained <- healthy.PC.importance[2,1:9]
healthy.PC.proportions.explained.plot <- ggplot(mapping = aes(x = as.factor(names(healthy.PC.proportions.explained)), y = healthy.PC.proportions.explained)) +  
                                                geom_bar(stat = "identity") + xlab("Components") + ylab("Proportion Variance Explained")
ggsave(file = "Healthy_Screeplot.pdf", plot = healthy.PC.proportions.explained.plot, device = "pdf", path = PCA.output)


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
healthy.PCA$x <- as.data.frame(healthy.PCA$x) %>% sample_n(5000, replace = TRUE) %>% as.matrix()
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
    pdf(file = paste0("healthy3D_vertical_", i, ".pdf"), width = 10, height = 10)
    healthy.3D <- scatter3D(x = healthy.PCA.plot$PC1, 
                            y = healthy.PCA.plot$PC2, 
                            z = healthy.PCA.plot$PC3, 
                            colvar = healthy.PCA.plot$CD34,
                            xlab = "PC1", 
                            ylab = "PC2", 
                            zlab = "PC3", 
                            clab = "CD34 Expression", 
                            pch = 19,
                            cex = 0.5,
                            phi = i,
                            theta = 220)
    healthy.3D
    dev.off()
}

for (i in seq(from = 0, to = 360, by = 10)){ 
  pdf(file = paste0("healthy3D_", i, ".pdf"), width = 10, height = 10)
  healthy.3D <- scatter3D(x = healthy.PCA.plot$PC1, 
                          y = healthy.PCA.plot$PC2, 
                          z = healthy.PCA.plot$PC3, 
                          colvar = healthy.PCA.plot$CD34,
                          xlab = "PC1", 
                          ylab = "PC2", 
                          zlab = "PC3", 
                          clab = "CD34 Expression", 
                          pch = 19,
                          cex = 0.5,
                          phi = 30,
                          theta = i)
  healthy.3D
  dev.off()
}

# for (i in seq(from = 0, to = 360, by = 10)){ 
#   pdf(filename = paste0("CD117_healthy3D_", i, ".pdf"), width = 10, height = 10, units = "in", res = 300)
#   healthy.3D <- scatter3D(x = healthy.PCA.plot$PC1, y = healthy.PCA.plot$PC2, z = healthy.PCA.plot$PC3, colvar = healthy.PCA.plot$CD117,
#                           pch = 19, cex = 0.5, phi = 30, theta = i)
#   healthy.3D
#   dev.off()
# }

setwd(file.path(OUTPUT.DIRECTORY, "data"))

message("Saving important data files...") 
save(list = c("healthy.demographics", "healthy.data", "healthy.data.scaled"), file = "healthy_data.RData")
save(list = c("cancer.demographics", "cancer.data", "cancer.data.scaled", "healthy.loadings", "healthy.PCA.plot"), file = "cancer_data.RData")

setwd(CODE.DIRECTORY)

rm(list = c("healthy.demographics", "healthy.data", "healthy.data.scaled",
            "cancer.demographics", "cancer.data", "cancer.data.scaled", 
            "healthy.loadings", "healthy.PCA.plot", "healthy.variances",
            "healthy.most.PC1", "healthy.most.PC2", "healthy.PC.drivers", 
            "healthy.PC1.rotation", "healthy.PC2.rotation", "min.cell.number", 
            "healthy.PC.proportions.explained", "healthy.PC.proportions.explained.plot",
            "healthy.3D", "healthy.plot.all", "healthy.plot.CD117", 
            "healthy.plot.CD34", "healthy.plot.CD38", "healthy.plot.patient", 
            "test.matrix", "healthy.scree.PCs", "cancer.variances", "healthy.PC.importance", 
            "healthy.PC1.markers.plot", "healthy.PC2.markers.plot", "healthy.PCA", 
            "healthy.PCs", "healthy.variances.plot", "cancer.variances.plot"))

message("Healthy PCA analysis completed!")
