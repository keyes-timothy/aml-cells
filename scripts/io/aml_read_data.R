# aml_read_data.R 



########################


#set up output directory
reader_output <- file.path(OUTPUT.DIRECTORY, "Reader_Plots")
dir.create(reader.output, showWarnings = FALSE)

#indicate some important working directory information
message("Setting Directory Information")
setwd(OUTPUT.DIRECTORY)
AML.directory <- file.path(getwd(), "..", "AML_data_raw_myeloid_compensated_corrected")

#list data files
file.list <- list.files(AML.directory)

#read files
message("Reading in data...")
setwd(file.path(CODE.DIRECTORY, "read_data"))
source("cytof_clean_names.R")
source("cytof_exprsExtract.R")
source("cytofAsinh.R")
source("cytof_exprsMerge.R")
source("clean_protein_name.R")

setwd(AML.directory)
expression.matrix <- cytof_exprsMerge(fcsFiles = file.list, 
                                      transformMethod = "cytofAsinh", 
                                      mergeMethod = "all", verbose = TRUE) %>% 
  cytof.clean.names() %>% as_tibble()

##########clean the data that was read in 
message("Cleaning column and row names...")

#clean protein names
needsCleaning <- grep(pattern = "_", x = colnames(expression.matrix))
new.colnames <- colnames(expression.matrix) 
new.colnames[needsCleaning] <- sapply(new.colnames[needsCleaning], clean.protein.name)
colnames(expression.matrix) <- new.colnames
rm(needsCleaning)
rm(new.colnames)

#remove channels that we don't care about, like the barcoding channels
message("Removing irrelevant channels")
indexIrrelevant <- grep(pattern = "<", x = colnames(expression.matrix))
expression.matrix <- dplyr::select(expression.matrix, -indexIrrelevant)
rm(indexIrrelevant)

#fix the stimulation column, which currently has both the condition and stimulation in it
stim.col.wrong <- grep(pattern = " ", x = expression.matrix$stimulation)
split.stim <- sapply(expression.matrix$stimulation[stim.col.wrong], strsplit, split = " ") %>% as.data.frame %>% t #roundabout way to split each part of the current stimulation variable
expression.matrix$condition <- rep("Healthy", times = nrow(expression.matrix))
expression.matrix$stimulation[stim.col.wrong] <- split.stim[,2]
expression.matrix$condition[stim.col.wrong] <- split.stim[,1]
rm(split.stim)
rm(stim.col.wrong)

message("Saving Files")

#save a sampled version of the giant dataframe for testing locally 
dir.create(file.path(OUTPUT.DIRECTORY, "data"))
setwd(file.path(OUTPUT.DIRECTORY, "data"))
sampled.matrix <- expression.matrix %>% group_by(patient, stimulation, condition) %>%
  sample_n(size = 100)
save(list = c("sampled.matrix"), file = "AML_sampled_matrix.RData")

#save giant data matrix to local directory
save(list = c("expression.matrix"), file = "AML_expression_matrix.RData")
setwd(CODE.DIRECTORY)


########provide some simple figures regarding the quality of the data 
cell.counts <- expression.matrix %>% group_by(patient, stimulation, condition) %>% 
  summarize(`Number of Cells` = n())
message("Plotting quality report...")

cell.counts.plot <- ggplot(data = cell.counts, mapping = aes(x = stimulation, y = `Number of Cells`, 
                                                             fill = condition)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~patient) + 
  xlab("Stimulation") + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9))
ggsave(filename = "cell.counts.pdf", plot = cell.counts.plot, device = "pdf", path = reader.output, 
       width = 10, height = 7, units = "in") 

rm(list = c("cell.counts", "cell.counts.plot", "expression.matrix"))

