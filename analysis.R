library(TCGAbiolinks)
library(SummarizedExperiment)

##Refeerence https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html
##https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html#TCGAanalyze_DEA__TCGAanalyze_LevelTab:_Differential_expression_analysis_(DEA)
#https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html


query.exp.hg38 <- GDCquery(
    project = "TCGA-LIHC", 
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification",workflow.type = "STAR - Counts")

GDCdownload(query.exp.hg38)

expdat <- GDCprepare(query = query.exp.hg38,save = TRUE, save.filename = "exp.rda")

BRCAMatrix <- assay(expdat,"unstranded") 
# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCA.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(expdat)
write.csv(BRCA.RNAseq_CorOutliers, file = "TCGA-LIHC_STARCounts_unstranded.csv", row.names = FALSE)

metadata<-as.data.frame(expdat@colData@listData)
write.csv(metadata, file = "TCGA-LIHC_metadata.csv", row.names = FALSE)




samplesDown <- getResults(query.exp.hg38,cols=c("cases"))



#############-------------------------------------------------------------------__###########

# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gprofiler2")
library("genefilter")



GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

TCGAbiolinks:::getProjectSummary("TCGA-LIHC")

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")

lihc_res = getResults(query_TCGA) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(lihc_res) # columns present in the table

head(lihc_res$sample_type) # first 6 types of tissue.

summary(factor(lihc_res$sample_type)) # summary of distinct tissues types present in this study


query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.type = "Gene Expression Quantification",
  data.category="Transcriptome Profiling",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_TCGA)

tcga_data = GDCprepare(query_TCGA)

dim(tcga_data)

colnames(colData(tcga_data))

table(tcga_data@colData$vital_status)

table(tcga_data@colData$tumor_stage)

sample_meetadata=as.data.frame(tcga_data@colData)

sample_meetadata <- sample_meetadata %>%unnest(colnames(sample_meetadata), names_sep = "_")
# Then, you can save the modified data frame to a CSV file
write.csv(sample_meetadata, file = "TCGA-LIHC.csv", row.names = FALSE)

table(tcga_data@colData$definition)

table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)
dim(assay(tcga_data))     # gene expression matrices.

head(assay(tcga_data)[,1:10]) # expression of first 6 genes and first 10 samples

head(rowData(tcga_data))     # ensembl id and gene id of the first 6 genes.
