if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
if(!require(limma)) BiocManager::install("limma")
if(!require(SummarizedExperiment)) BiocManager::install("SummarizedExperiment")
if(!require(matrixStats)) BiocManager::install("matrixStats")
if(!require(DESeq2)) BiocManager::install("DESeq2")

library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(matrixStats)
library(DESeq2)

# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

### Load in HTSeq Counts from TCGA ###
GDC_query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(GDC_query)
GDC_sum_exp <- GDCprepare(GDC_query)

transcriptomics_GDC <- assays(GDC_sum_exp)$"HTSeq - Counts"

### Pre-processing ###
patients_no_NA_mask <- !is.na(colData(GDC_sum_exp)$paper_age_at_initial_pathologic_diagnosis)
patient_data_GDC <- colData(GDC_sum_exp)[ patients_no_NA_mask, ] #filters out patients with NA age

patient_data_GDC$age_category = ifelse(patient_data_GDC$paper_age_at_initial_pathologic_diagnosis < 40, "Young", 
                             ifelse(patient_data_GDC$paper_age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age category column, categorized by young, mid, old
patient_data_GDC$age_category <- factor(patient_data_GDC$age_category, levels=c("Young", "Mid", "Old")) #categorical column for age category
patient_data_GDC$paper_BRCA_Subtype_PAM50 <- factor(patient_data_GDC$paper_BRCA_Subtype_PAM50, levels=c("Her2", "LumA", "LumB", "Basal", "Normal")) #categorical column for PAM50

transcriptomics_GDC <- transcriptomics_GDC[rowMeans(transcriptomics_GDC) >= 10, patients_no_NA_mask] #filters out genes where mean counts < 10, filters out patients with NA age

### DESeq2 adjusted by PAM50 subtype ###

dds_with_adjustment <- DESeqDataSetFromMatrix(countData = transcriptomics_GDC, colData = patient_data_GDC, design = ~age_category +paper_BRCA_Subtype_PAM50) #creates a DESeqDataSet object from the transcriptomic and patient data matrices, using age category and PAM50 subtype as factors
dds_with_adjustment_obj <- DESeq(dds_with_adjustment) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values
resultsNames(dds_with_adjustment_obj) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young"

results_with_adjustment <- results(dds_with_adjustment_obj, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values, contrast=c specifies the comparison for the fold change, with age_category_median as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

results_with_adjustment$FoldChange <- 2^results_with_adjustment$log2FoldChange #creates column in results with FoldChange

# Volcano plot #
padj_threshold <- 0.05 #sets significance threshold
log2FC_threshold <- 1.0 #sets fold change threshold
jpeg("final_project_data/DESeq_Volcano_Plot_With_Adjustment.jpg")
plot(x= results_with_adjustment$log2FoldChange, y= -log10(results_with_adjustment$padj), main= "Differentially Expressed Genes (Adjusted for PAM50 Subtype)", xlab= "log2(Fold Change)", ylab="-log10(p-adjusted value)") #plots significance (-log10(padj)) vs expression (log2FoldChange)
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

# Saves Results #
results_with_adjustment_no_NA_mask <- !is.na(results_with_adjustment$padj)
results_with_adjustment <- results_with_adjustment[results_with_adjustment_no_NA_mask, ]

results_significant_padj_with_adjustment <- results_with_adjustment[results_with_adjustment$padj < padj_threshold, ] #filters for only significant (padj < 0.05) genes

results_sig_up_regulated_with_adjustment <- results_significant_padj_with_adjustment[results_significant_padj_with_adjustment$log2FoldChange > log2FC_threshold, ] #significantly UP regulated genes
results_sig_down_regulated_with_adjustment <- results_significant_padj_with_adjustment[results_significant_padj_with_adjustment$log2FoldChange < -log2FC_threshold, ] #significantly DOWN regulated genes

gene_information <- rowData(GDC_sum_exp)
results_sig_up_regulated_with_adjustment$CommonGeneName <- gene_information[rownames(results_sig_up_regulated_with_adjustment), 2] #adds column with common gene names
results_sig_down_regulated_with_adjustment$CommonGeneName <- gene_information[rownames(results_sig_down_regulated_with_adjustment), 2] #adds column with common gene names
results_with_adjustment$CommonGeneName <- gene_information[rownames(results_with_adjustment), 2] #adds column with common gene name

write.csv(results_sig_up_regulated_with_adjustment, "final_project_data/results_sig_up_regulated_with_adjustment.csv")
write.csv(results_sig_down_regulated_with_adjustment, "final_project_data/results_sig_down_regulated_with_adjustment.csv")
write.csv(results_with_adjustment, "final_project_data/DESeq_Data_With_Adjustment.csv")

### Kaplan-Meier Plots ###
results_sig_up_regulated_with_adjustment_sorted <- results_sig_up_regulated_with_adjustment[order(results_sig_up_regulated_with_adjustment$log2FoldChange, decreasing = TRUE), ] #sorts significantly up-regulated results by log2FoldChange
results_sig_down_regulated_with_adjustment_sorted <- results_sig_down_regulated_with_adjustment[order(results_sig_down_regulated_with_adjustment$log2FoldChange), ] #sorts significantly down-regulated results by log2FoldChange

for (i in 1:5) { #creates KM plots for top 5 up and down regulated genes
  ## UP Regulated Genes ##
  
  # Prepares Variables for Specific Gene #
  gene_name <- results_sig_up_regulated_with_adjustment_sorted[i, 8]
  ENSG <- rownames(results_sig_up_regulated_with_adjustment_sorted)[i]
  legend <- paste(gene_name, "Expression Level")
  filename_pdf <- paste("final_project_data/survival_expression_of_",gene_name,".pdf",sep="")
  filename_jpg<- paste("final_project_data/survival_expression_of_",gene_name,".jpg",sep="")
  
  # Separates High, Low, No Expression Levels #
  zero_mask <- transcriptomics_GDC[ENSG, ] == 0
  zero_count <- sum(zero_mask)
  average <- sum(transcriptomics_GDC[ENSG, ])/(ncol(transcriptomics_GDC)-zero_count)
  patient_data_GDC$kaplan_meier_expression = ifelse(transcriptomics_GDC[ENSG, patient_data_GDC$barcode] > average, "High", ifelse(transcriptomics_GDC[ENSG, patient_data_GDC$barcode] == 0, "No Counts", "Low"))
  
  # Creates and Saves KM Plot #
  TCGAanalyze_survival( patient_data_GDC, "kaplan_meier_expression", legend= legend, filename= filename_pdf)
  TCGAanalyze_survival( patient_data_GDC, "kaplan_meier_expression", legend= legend, filename= filename_jpg)

  ## DOWN Regulated Genes ##
  
  # Prepares Variables for Specific Gene #
  gene_name <- results_sig_down_regulated_with_adjustment_sorted[i, 8]
  ENSG <- rownames(results_sig_down_regulated_with_adjustment_sorted)[i]
  legend <- paste(gene_name, "Expression Level")
  filename_pdf <- paste("final_project_data/survival_expression_of_",gene_name,".pdf",sep="")
  filename_jpg<- paste("final_project_data/survival_expression_of_",gene_name,".jpg",sep="")
  
  # Separates High, Low, No Expression Levels #
  zero_mask <- transcriptomics_GDC[ENSG, ] == 0
  zero_count <- sum(zero_mask)
  average <- sum(transcriptomics_GDC[ENSG, ])/(ncol(transcriptomics_GDC)-zero_count)
  patient_data_GDC$kaplan_meier_expression = ifelse(transcriptomics_GDC[ENSG, patient_data_GDC$barcode] > average, "High", ifelse(transcriptomics_GDC[ENSG, patient_data_GDC$barcode] == 0, "No Counts", "Low"))
  
  # Creates and Saves KM Plot #
  TCGAanalyze_survival( patient_data_GDC, "kaplan_meier_expression", legend= legend, filename= filename_pdf)
  TCGAanalyze_survival( patient_data_GDC, "kaplan_meier_expression", legend= legend, filename= filename_jpg)
}

# from limma #

limma_results_proteins_kaplan_meier <- list("ESR1", "PTGR2")

for (protein in limma_results_proteins_kaplan_meier) { #creates KM plots for top 5 up and down regulated genes
  ## UP Regulated Genes ##
  
  # Prepares Variables for Specific Gene #
  gene_name <- protein
  ENSG <- rownames(results_with_adjustment)[results_with_adjustment$CommonGeneName == protein]
  legend <- paste(gene_name, "Expression Level")
  filename_pdf <- paste("final_project_data/survival_expression_of_",gene_name,".pdf",sep="")
  filename_jpg<- paste("final_project_data/survival_expression_of_",gene_name,".jpg",sep="")
  
  # Separates High, Low, No Expression Levels #
  zero_mask <- transcriptomics_GDC[ENSG, ] == 0
  zero_count <- sum(zero_mask)
  average <- sum(transcriptomics_GDC[ENSG, ])/(ncol(transcriptomics_GDC)-zero_count)
  patient_data_GDC$kaplan_meier_expression = ifelse(transcriptomics_GDC[ENSG, patient_data_GDC$barcode] > average, "High", ifelse(transcriptomics_GDC[ENSG, patient_data_GDC$barcode] == 0, "No Counts", "Low"))
  
  # Creates and Saves KM Plot #
  TCGAanalyze_survival( patient_data_GDC, "kaplan_meier_expression", legend= legend, filename= filename_pdf)
  TCGAanalyze_survival( patient_data_GDC, "kaplan_meier_expression", legend= legend, filename= filename_jpg)
}
