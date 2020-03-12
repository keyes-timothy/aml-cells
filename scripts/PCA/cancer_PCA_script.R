#####################  PCA analysis of AML data   ###########################

#Set-up dataset
setwd(file.path(OUTPUT.DIRECTORY, "data"))
load("cancer_data.RData")
cancer.demographics$patient[cancer.demographics$patient == "PASTRP"] <- "PASRTP"

#PCA with cancer cells
message("Plotting PCA for cancer samples...")
cancer.PCA <- prcomp(cancer.data.scaled, scale = FALSE) 
cancer.loadings <- cancer.PCA$rotation 
cancer.PCs <- cancer.PCA$x
cancer.PC.importance <- summary(cancer.PCA)$importance

#plot cancer
axis.label.style <- element_text(face = "bold", size = 24)
axis.style <- element_text(size = 14)
PCA.theme <- theme(axis.title.x = axis.label.style, axis.title.y = axis.label.style, 
                   axis.text = axis.style, legend.title = element_text(face = "bold", size = 16), 
                   legend.text = element_text(size = 14)) + theme_bw()

min.cell.number.cancer <- cbind(cancer.PCs, cancer.demographics) %>% group_by(patient, condition) %>% count() %>% ungroup() %>% 
                          select(n) %>% min()
cancer.demographics$CD117 <- cancer.data$CD117
cancer.PCA.plot <- cbind(cancer.PCs, cancer.demographics) %>% 
  group_by(patient, condition) %>% 
  sample_n(size = (500), replace = TRUE) #can be varied
cancer.PCA.plot$patient <- as.factor(cancer.PCA.plot$patient)

#Save cancer PCA dataframe for plotting locally
message("Saving cancer PCA dataframe for local plotting")
setwd(PCA.output)
cancer.PCA.df <- cancer.PCA.plot
save(list = c("cancer.PCA.df"), file = "sampled_cancer_PCA_df.RData")
rm(cancer.PCA.df)

setwd(OUTPUT.DIRECTORY)

cancer.plot <- ggplot(data = cancer.PCA.plot, mapping = aes(x=PC1, y=PC2, color = patient)) + geom_point(size = 1, alpha = 0.4)
cancer.plot <- cancer.plot + 
               xlab(paste0("PC1 (", round(cancer.PC.importance["Proportion of Variance", "PC1"]*100, digits = 2), "%)")) + 
               ylab(paste0("PC2 (", round(cancer.PC.importance["Proportion of Variance", "PC2"]*100, digits = 2), "%)")) + 
               facet_wrap(patient~condition) + theme(legend.position = "none") + 
               PCA.theme
ggsave(file = "cancer_PCA_plot.pdf", plot = cancer.plot, device = "pdf", path = PCA.output, 
       height = 9, width = 11, units = "in")

#Express cancer using healthy coordinate system - healthy trajectory is in black
min.cell.number <- healthy.PCA.plot %>% count(patient) %>% ungroup() %>% select(n) %>% min()
message("Plotting cancer samples under the healthy coordinate system...")
cancer.under.healthy.PCs <- as.matrix(cancer.data.scaled) %*% as.matrix(healthy.loadings) 
cancer.under.healthy.plot <- cbind(cancer.under.healthy.PCs, cancer.demographics) %>% 
                            as_tibble() %>% group_by(patient, condition) %>%
                            sample_n(size = 500, replace = TRUE) 
cancer.under.healthy.plot$patient <- as.factor(cancer.under.healthy.plot$patient)
cancer.u.healthy.plot <- ggplot(data = cancer.under.healthy.plot, mapping = aes(x=PC1, y=PC2, color = condition)) +
                         geom_point(size = 1, alpha = 0.4) + facet_wrap(~patient) +
                         guides(colour = guide_legend(override.aes = list(size=5)))
message("Adding healthy panels to cancer plots...")
cancer.u.healthy.plot <- cancer.u.healthy.plot +
                         geom_point(data = sample_n(healthy.PCA.plot, size = 500, replace = TRUE), 
                                     color = "black",alpha = 0.5, size = 0.5) + PCA.theme #+  
                         #theme(legend.position = "none")
ggsave(file = "cancer_under_healthy_plot.pdf", plot = cancer.u.healthy.plot, device = "pdf", path = PCA.output, height = 10, width = 12, units = "in")

#quantitative representation of which markers are driving each PC in cancer
message("Making cancer biplots...")
library(ggfortify)
#same in cancer 
cancer.PCA$x <- as.data.frame(cancer.PCA$x) %>% sample_n(10000, replace = TRUE) %>% as.matrix()
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
ggsave(file = "Cancer_PCA_drivers.pdf", plot = cancer.PC.drivers, device = "pdf", path = PCA.output, 
       width = 10, height = 10, units = "in")


message("Plotting 3D Cancer plots...")
library(plot3D)
setwd(PCA.output)
# 3D plotting - cancer
cancer.under.healthy.plot <- cancer.under.healthy.plot %>% group_by(patient) %>% sample_n(size = 100, replace = TRUE) 
for (i in seq(from = 0, to = 360, by = 10)){ 
  pdf(file = paste0("cancer3D_", i, ".pdf"), width = 10, height = 10)
  cancer.u.healthy.3D <- scatter3D(x = cancer.under.healthy.plot$PC1, 
                                    y = cancer.under.healthy.plot$PC2, 
                                    z = cancer.under.healthy.plot$PC3, 
                                    colvar = as.numeric(as.factor(cancer.under.healthy.plot$patient)), 
                                   xlab = "PC1", 
                                   ylab = "PC2", 
                                   zlab = "PC3", 
                                   colkey = FALSE, 
                                    pch = 19, 
                                    cex = 0.5, 
                                    phi = 30, 
                                    theta = i)
    cancer.u.healthy.3D
    dev.off()
}

for (i in seq(from = 0, to = 360, by = 10)){ 
  pdf(file = paste0("cancer3D_vertical_", i, ".pdf"), width = 10, height = 10)
  cancer.u.healthy.3D <- scatter3D(x = cancer.under.healthy.plot$PC1, 
                                   y = cancer.under.healthy.plot$PC2, 
                                   z = cancer.under.healthy.plot$PC3, 
                                   colvar = as.numeric(as.factor(cancer.under.healthy.plot$patient)), 
                                   xlab = "PC1", 
                                   ylab = "PC2", 
                                   zlab = "PC3", 
                                   colkey = FALSE, 
                                   pch = 19, 
                                   cex = 0.5, 
                                   phi = i, 
                                   theta = 220)
  cancer.u.healthy.3D
  dev.off()
}

rm(list = c("cancer.data", "cancer.data.scaled", "cancer.demographics", 
            "cancer.loadings", "cancer.PC.drivers", "cancer.PC.importance", 
            "cancer.PCA", "cancer.PCA.plot", "cancer.PCs", "cancer.plot", 
            "cancer.u.healthy.3D", "cancer.u.healthy.plot", "cancer.under.healthy.PCs", 
            "cancer.under.healthy.plot", "cancer.most.PC1", "cancer.most.PC2", 
            "cancer.scree.PCs", "healthy.PCA.plot", "healthy.loadings"))

message("Cancer PCA analysis completed!")

