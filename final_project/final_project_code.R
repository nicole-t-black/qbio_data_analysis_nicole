# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html

### Load Packages ###
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
if(!require(limma)) BiocManager::install("limma")
if(!require(SummarizedExperiment)) BiocManager::install("SummarizedExperiment")
if(!require(matrixStats)) BiocManager::install("matrixStats")
if(!require(DESeq2)) BiocManager::install("DESeq2")
if(!require(pheatmap)) BiocManager::install("pheatmap")

library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(matrixStats)
library(DESeq2)
library(pheatmap)

################################################################################

### Load in HTSeq Counts from TCGA ###
query_GDC <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
GDCdownload(query_GDC)
sum_exp_GDC <- GDCprepare(query_GDC)

transcriptomics_GDC <- assays(sum_exp_GDC)$"HTSeq - Counts"

### Pre-processing ###
clinical_GDC_no_NA_mask <- !is.na(colData(sum_exp_GDC)$paper_age_at_initial_pathologic_diagnosis)
clinical_GDC <- colData(sum_exp_GDC)[ clinical_GDC_no_NA_mask, ] #filters out patients with NA age

clinical_GDC$age_category = ifelse(clinical_GDC$paper_age_at_initial_pathologic_diagnosis < 40, "Young", 
                                       ifelse(clinical_GDC$paper_age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age category column, categorized by young, mid, old
clinical_GDC$age_category <- factor(clinical_GDC$age_category, levels=c("Young", "Mid", "Old")) #categorical column for age category
clinical_GDC$paper_BRCA_Subtype_PAM50 <- factor(clinical_GDC$paper_BRCA_Subtype_PAM50, levels=c("Her2", "LumA", "LumB", "Basal", "Normal")) #categorical column for PAM50

transcriptomics_GDC <- transcriptomics_GDC[rowMeans(transcriptomics_GDC) >= 10, clinical_GDC_no_NA_mask] #filters out genes where mean counts < 10, filters out patients with NA age

# Saves Transcriptomics and Clinical Data #

write.csv(transcriptomics_GDC, "final_project_data/transcriptomics_GDC.csv")
write.csv(clinical_GDC, "final_project_data/clinical_GDC.csv")

### DESeq2 adjusted by PAM50 subtype ###

dds <- DESeqDataSetFromMatrix(countData = transcriptomics_GDC, colData = clinical_GDC, design = ~age_category +paper_BRCA_Subtype_PAM50) #creates a DESeqDataSet object from the transcriptomic and patient data matrices, using age category and PAM50 subtype as factors
dds_obj <- DESeq(dds) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values

deseq_results <- results(dds_obj, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values, contrast=c specifies the comparison for the fold change, with age_category_median as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

deseq_results$FoldChange <- 2^deseq_results$log2FoldChange #creates column in results with FoldChange

# Volcano plot #
padj_threshold <- 0.05 #sets significance threshold
log2FC_threshold <- 1.0 #sets fold change threshold
jpeg("final_project_figures/DESeq_volcano_plot.jpg")
plot(x= deseq_results$log2FoldChange, y= -log10(deseq_results$padj), main= "Differentially Expressed Genes (Adjusted for PAM50 Subtype)", xlab= "log2(Fold Change)", ylab="-log10(p-adjusted value)") #plots significance (-log10(padj)) vs expression (log2FoldChange)
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

# Polishes and Saves Results #
deseq_results_no_NA_mask <- !is.na(deseq_results$padj)
deseq_results <- deseq_results[deseq_results_no_NA_mask, ]

gene_information <- rowData(sum_exp_GDC)
deseq_results$CommonGeneName <- gene_information[rownames(deseq_results), 2] #adds column with common gene name

deseq_results_all <- deseq_results
deseq_results <- deseq_results[deseq_results$padj < padj_threshold, ] #filters for only significant (padj < 0.05) genes

deseq_results_up_reg <- deseq_results[deseq_results$log2FoldChange > log2FC_threshold, ] #significantly UP regulated genes
deseq_results_down_reg <- deseq_results[deseq_results$log2FoldChange < -log2FC_threshold, ] #significantly DOWN regulated genes

deseq_results_up_reg <- deseq_results_up_reg[order(deseq_results_up_reg$log2FoldChange, decreasing = TRUE), ] #sorts significantly up-regulated results by log2FoldChange
deseq_results_down_reg <- deseq_results_down_reg[order(deseq_results_down_reg$log2FoldChange), ] #sorts significantly down-regulated results by log2FoldChange

write.csv(deseq_results, "final_project_data/deseq_results.csv")
write.csv(deseq_results_up_reg, "final_project_data/deseq_results_up_reg.csv")
write.csv(deseq_results_down_reg, "final_project_data/deseq_results_down_reg.csv")

### Kaplan-Meier - DESeq2 ###
for (i in 1:5) { #creates KM plots for top 5 up and down regulated genes
  ## UP Regulated Genes ##
  
  # Prepares Variables for Specific Gene #
  km_gene_name <- deseq_results_up_reg[i, 8]
  km_ENSG <- rownames(deseq_results_up_reg)[i]
  km_legend <- paste(km_gene_name, "Expression Level")
  km_filename <- paste("final_project_figures/deseq_km_", km_gene_name,"_up_", i, ".jpg", sep="")
  
  # Separates High, Low, No Expression Levels #
  km_zero_mask <- transcriptomics_GDC[km_ENSG, ] == 0
  km_zero_count <- sum(km_zero_mask)
  km_average <- sum(transcriptomics_GDC[km_ENSG, ])/(ncol(transcriptomics_GDC)-km_zero_count)
  clinical_GDC$km_expression = ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] > km_average, "High", ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] == 0, "No Counts", "Low"))
  
  # Creates and Saves KM Plot #
  TCGAanalyze_survival( clinical_GDC, "km_expression", legend= km_legend, filename= km_filename)
  
  ## DOWN Regulated Genes ##
  
  # Prepares Variables for Specific Gene #
  km_gene_name <- deseq_results_down_reg[i, 8]
  km_ENSG <- rownames(deseq_results_down_reg)[i]
  km_legend <- paste(km_gene_name, "Expression Level")
  km_filename <- paste("final_project_figures/deseq_km_", km_gene_name,"_down_", i, ".jpg", sep="")
  
  # Separates High, Low, No Expression Levels #
  km_zero_mask <- transcriptomics_GDC[km_ENSG, ] == 0
  km_zero_count <- sum(km_zero_mask)
  km_average <- sum(transcriptomics_GDC[km_ENSG, ])/(ncol(transcriptomics_GDC)-km_zero_count)
  clinical_GDC$km_expression = ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] > km_average, "High", ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] == 0, "No Counts", "Low"))
  
  # Creates and Saves KM Plot #
  TCGAanalyze_survival( clinical_GDC, "km_expression", legend= km_legend, filename= km_filename)
}

################################################################################

### Gather CPTAC data (from python code) ###
clinical_CPTAC <- read.csv("final_project_data/clinical_CPTAC.csv", header = TRUE)
rownames(clinical_CPTAC) <- clinical_CPTAC$Patient_ID

proteomics_CPTAC <- read.csv("final_project_data/proteomics_CPTAC.csv", header = TRUE)
rownames(proteomics_CPTAC) <- proteomics_CPTAC$Patient_ID
proteomics_CPTAC$Patient_ID <- NULL

### Filter CPTAC data -- pre-processing ###
clinical_CPTAC$Age.in.Years <- clinical_CPTAC[ ,4]/12 #add Age in Years column

clinical_CPTAC$age_category <- ifelse(clinical_CPTAC$Age.in.Years < 40, "Young", 
                                      ifelse(clinical_CPTAC$Age.in.Years >= 60, "Old", "Mid")) #add age categories, old/mid/young

clinical_CPTAC_no_NA_mask <- !is.na(clinical_CPTAC$Age.in.Month)
clinical_CPTAC <- clinical_CPTAC[ clinical_CPTAC_no_NA_mask , ] #filter out patients with NA age

clinical_CPTAC_no_mid_mask <- clinical_CPTAC$age_category != "Mid"
clinical_CPTAC <-clinical_CPTAC[clinical_CPTAC_no_mid_mask, ] #filter out patients with mid age-category

proteomics_CPTAC <- proteomics_CPTAC[clinical_CPTAC_no_NA_mask , ]
proteomics_CPTAC <- proteomics_CPTAC[clinical_CPTAC_no_mid_mask, ] #filter out NA and mid age-categories

proteomics_CPTAC <- t(proteomics_CPTAC) #transpose proteomics dataframe --> rows = proteins, columns = patients

### LIMMA ###

# Analysis #
age_factor <- factor( clinical_CPTAC$age_category, levels=c("Young", "Old")) #create categorical variables
PAM50_factor <- factor( clinical_CPTAC$PAM50, levels=c("Her2", "LumA", "LumB", "Basal", "Normal")) #create vategorical variables
design_matrix <- model.matrix(~0 +age_factor +PAM50_factor, clinical_CPTAC) #design matrix, containing categorical variables (young/old, PAM50 for adjustment)

fit <- lmFit(proteomics_CPTAC, design_matrix) #limma lmFit(data, design_matrix) function

contrast_matrix <- makeContrasts(Y_O="age_factorYoung-age_factorOld", levels=design_matrix) #contrast matrix to express the different conditions

fit2 <- contrasts.fit(fit, contrast_matrix) #save the result of the contrasts.fit(fit, contrasts) 

fit3 <- eBayes(fit2) #eBayes to smooth error

limma_results <- topTable(fit3, adjust="BH", n=Inf) #topTable to get statistics for top differentially expressed

# Volcano Plot #
padj_threshold <- 0.05 #sets significance threshold
log2FC_threshold <- 1.0 #sets fold change threshold
jpeg("final_project_figures/CPTAC_volcano_plot.jpg")
plot(x= limma_results$logFC, y= -log10(limma_results$adj.P.Val), main= "Differentially Expressed Proteins (Adjusted for PAM50 Subtype)", xlab="log2(Fold Change)", ylab="-log10(p-adjusted value)" ) #plots significance (-log10(padj)) vs expression (log2FoldChange)
#plot(x= filter$logFC, y= -log10(filter$finalFDR) )
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

# Polishes and Saves Results #
limma_results <- limma_results[limma_results$adj.P.Val < padj_threshold , ] #filters for only significant (padj < 0.05) proteins

limma_results_up_reg <- limma_results[limma_results$logFC > log2FC_threshold, ] #significantly UP regulated proteins
limma_results_down_reg <- limma_results[limma_results$logFC < -log2FC_threshold, ] #significantly DOWN regulated proteins

limma_results_up_reg <- limma_results_up_reg[order(limma_results_up_reg$logFC, decreasing = TRUE), ] #sorts significantly up-regulated results by log2FoldChange
limma_results_down_reg <- limma_results_down_reg[order(limma_results_down_reg$logFC), ] #sorts significantly down-regulated results by log2FoldChange

write.csv(limma_results, "final_project_data/limma_results.csv")
write.csv(limma_results_up_reg, "final_project_data/limma_results_up_reg.csv")
write.csv(limma_results_down_reg, "final_project_data/limma_results_down_reg.csv")

### Kaplan-Meier - limma ###
limma_results_up_proteins <- rownames(limma_results_up_reg) #gets names of up regulated proteins
limma_results_down_proteins <- rownames(limma_results_down_reg) #gets names of down regulated proteins

for (protein in limma_results_up_proteins) { #removes proteins with NAs
  if (grepl("NA", protein, fixed=TRUE) == TRUE) {
    limma_results_up_proteins <- limma_results_up_proteins[limma_results_up_proteins != protein]
  }
}

for (protein in limma_results_down_proteins) { #removes proteins with NAs
  if (grepl("NA", protein, fixed=TRUE) == TRUE) {
    limma_results_down_proteins <- limma_results_down_proteins[limma_results_down_proteins != protein]
  }
}

limma_results_up_proteins_length <- length(limma_results_up_proteins)
limma_results_down_proteins_length <- length(limma_results_down_proteins)

if (limma_results_up_proteins_length >= 5) {
  for (i in 1:5) { #creates KM plots for top 5 up regulated proteins
    
    # Prepares Variables for Specific Gene #
    km_protein_name <- limma_results_up_proteins[i]
    km_ENSG <- rownames(deseq_results_all)[deseq_results_all$CommonGeneName == km_protein_name]
    km_legend <- paste(km_protein_name, "Expression Level")
    km_filename <- paste("final_project_figures/limma_km_", km_protein_name,"_up_", i, ".jpg", sep="")
    
    # Separates High, Low, No Expression Levels #
    km_zero_mask <- transcriptomics_GDC[km_ENSG, ] == 0
    km_zero_count <- sum(km_zero_mask)
    km_average <- sum(transcriptomics_GDC[km_ENSG, ])/(ncol(transcriptomics_GDC)-km_zero_count)
    clinical_GDC$km_expression = ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] > km_average, "High", ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] == 0, "No Counts", "Low"))
    
    # Creates and Saves KM Plot #
    TCGAanalyze_survival( clinical_GDC, "km_expression", legend= km_legend, filename= km_filename)
  }
} else if (limma_results_up_proteins_length == 0){
    for (i in 1:limma_results_up_proteins_length) {
      print(limma_results_up_proteins_length)
  }
} else {
  for (i in 1:limma_results_up_proteins_length) { #creates KM plots all up regulated proteins
    
    # Prepares Variables for Specific Gene #
    km_protein_name <- limma_results_up_proteins[i]
    km_ENSG <- rownames(deseq_results_all)[deseq_results_all$CommonGeneName == km_protein_name]
    km_legend <- paste(km_protein_name, "Expression Level")
    km_filename <- paste("final_project_figures/limma_km_", km_protein_name,"_up_", i, ".jpg", sep="")
    
    # Separates High, Low, No Expression Levels #
    km_zero_mask <- transcriptomics_GDC[km_ENSG, ] == 0
    km_zero_count <- sum(km_zero_mask)
    km_average <- sum(transcriptomics_GDC[km_ENSG, ])/(ncol(transcriptomics_GDC)-km_zero_count)
    clinical_GDC$km_expression = ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] > km_average, "High", ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] == 0, "No Counts", "Low"))
    
    # Creates and Saves KM Plot #
    TCGAanalyze_survival( clinical_GDC, "km_expression", legend= km_legend, filename= km_filename)
  }
}

if (limma_results_down_proteins_length >= 5) {
  for (i in 1:5) { #creates KM plots for top 5 up regulated proteins
    
    # Prepares Variables for Specific Gene #
    km_protein_name <- limma_results_diwn_proteins[i]
    km_ENSG <- rownames(deseq_results_all)[deseq_results_all$CommonGeneName == km_protein_name]
    km_legend <- paste(km_protein_name, "Expression Level")
    km_filename <- paste("final_project_figures/limma_km_", km_protein_name,"_down_", i, ".jpg", sep="")
    
    # Separates High, Low, No Expression Levels #
    km_zero_mask <- transcriptomics_GDC[km_ENSG, ] == 0
    km_zero_count <- sum(km_zero_mask)
    km_average <- sum(transcriptomics_GDC[km_ENSG, ])/(ncol(transcriptomics_GDC)-km_zero_count)
    clinical_GDC$km_expression = ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] > km_average, "High", ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] == 0, "No Counts", "Low"))
    
    # Creates and Saves KM Plot #
    TCGAanalyze_survival( clinical_GDC, "km_expression", legend= km_legend, filename= km_filename)
  }
} else if (limma_results_down_proteins_length == 0){
  for (i in 1:limma_results_down_proteins_length) {
    print(limma_results_down_proteins_length)
  }
} else {
  for (i in 1:limma_results_down_proteins_length) { #creates KM plots all up regulated proteins
    
    # Prepares Variables for Specific Gene #
    km_protein_name <- limma_results_down_proteins[i]
    km_ENSG <- rownames(deseq_results_all)[deseq_results_all$CommonGeneName == km_protein_name]
    km_legend <- paste(km_protein_name, "Expression Level")
    km_filename <- paste("final_project_figures/limma_km_", km_protein_name,"_down_", i, ".jpg", sep="")
    
    # Separates High, Low, No Expression Levels #
    km_zero_mask <- transcriptomics_GDC[km_ENSG, ] == 0
    km_zero_count <- sum(km_zero_mask)
    km_average <- sum(transcriptomics_GDC[km_ENSG, ])/(ncol(transcriptomics_GDC)-km_zero_count)
    clinical_GDC$km_expression = ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] > km_average, "High", ifelse(transcriptomics_GDC[km_ENSG, clinical_GDC$barcode] == 0, "No Counts", "Low"))
    
    # Creates and Saves KM Plot #
    TCGAanalyze_survival( clinical_GDC, "km_expression", legend= km_legend, filename= km_filename)
  }
}
