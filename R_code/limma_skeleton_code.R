# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html

### Load Packages ###
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
if (!require(limma)) BiocManager::install("limma")
if(!require(SummarizedExperiment)) BiocManager::install("SummarizedExperiment")
if(!require(matrixStats)) BiocManager::install("matrixStats")

library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(matrixStats)

### Gather CPTAC data (from python code) ###
clinical_CPTAC <- read.csv("../python_code/data/clinical_data.csv", header = TRUE)
rownames(clinical_CPTAC) <- clinical_CPTAC$Patient_ID

proteomics_CPTAC <- read.csv("../python_code/data/protein_data.csv", header = TRUE)
rownames(proteomics_CPTAC) <- proteomics_CPTAC$Patient_ID
proteomics_CPTAC$Patient_ID <- NULL

### Filter CPTAC data -- pre-processing ###
clinical_CPTAC$Age.in.Years <- clinical_CPTAC[ ,4]/12 #add Age in Years column

clinical_CPTAC$age_category <- ifelse(clinical_CPTAC$Age.in.Years < 40, "Young", 
                                ifelse(clinical_CPTAC$Age.in.Years >= 60, "Old", "Mid")) #add age categories, old/mid/young

clinical_no_NA_mask <- !is.na(clinical_CPTAC$Age.in.Month)
clinical_CPTAC <- clinical_CPTAC[ clinical_no_NA_mask , ] #filter out patients with NA age

clinical_no_mid_mask <- clinical_CPTAC$age_category != "Mid"
clinical_CPTAC <-clinical_CPTAC[clinical_no_mid_mask, ] #filter out patients with mid age-category

proteomics_CPTAC <- proteomics_CPTAC[clinical_no_NA_mask , ]
proteomics_CPTAC <- proteomics_CPTAC[clinical_no_mid_mask, ] #filter out NA and mid age-categories
  
proteomics_CPTAC_t <- t(proteomics_CPTAC) #transpose proteomics dataframe --> rows = proteins, columns = patients
  
### LIMMA ###

# Analysis #
f_age <- factor( clinical_CPTAC$age_category, levels=c("Young", "Old")) #create categorical variables
f_PAM50 <- factor( clinical_CPTAC$PAM50, levels=c("Her2", "LumA", "LumB", "Basal", "Normal")) #create vategorical variables
design <- model.matrix(~0+ f_age +f_PAM50, clinical_CPTAC) #design matrix, containing categorical variables (young/old, PAM50 for adjustment)

fit <- lmFit(proteomics_CPTAC_t, design) #limma lmFit(data, design_matrix) function

cont.matrix <- makeContrasts(Y_O="f_ageYoung-f_ageOld", levels=design) #contrast matrix to express the different conditions

fit2 <- contrasts.fit(fit, cont.matrix) #save the result of the contrasts.fit(fit, contrasts) 
  
fit3 <- eBayes( fit2 ) #eBayes to smooth error

limma_results_with_adjustment <- topTable( fit3, adjust="BH", n=Inf) #topTable to get statistics for top differentially expressed

# Volcano Plot #
padj_threshold <- 0.05 #sets significance threshold
log2FC_threshold <- 1.0 #sets fold change threshold
jpeg("final_project_data/CPTAC_Volcano_Plot_With_Adjustment.jpg")
plot(x= limma_results_with_adjustment$logFC, y= -log10(limma_results_with_adjustment$adj.P.Val), main= "Differentially Expressed Proteins (Adjusted for PAM50 Subtype)", xlab="log2(Fold Change)", ylab="-log10(p-adjusted value)" ) #plots significance (-log10(padj)) vs expression (log2FoldChange)
#plot(x= filter$logFC, y= -log10(filter$finalFDR) )
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

# Save Results #
write.table(limma_results_with_adjustment, "final_project_data/CPTAC_limma.txt", sep = "\t",
            row.names = F, quote = F)

write.csv(limma_results_with_adjustment, "final_project_data/CPTAC_limma_results_with_adjustment.csv")


limma_results_significant_padj_with_adjustment <- limma_results_with_adjustment[limma_results_with_adjustment$adj.P.Val < padj_threshold , ] #filters for only significant (padj < 0.05) proteins

limma_results_sig_up_regulated_with_adjustment <- limma_results_significant_padj_with_adjustment[limma_results_significant_padj_with_adjustment$logFC > log2FC_threshold, ] #significantly UP regulated proteins
limma_results_sig_down_regulated_with_adjustment <- limma_results_significant_padj_with_adjustment[limma_results_significant_padj_with_adjustment$logFC < -log2FC_threshold, ] #Dsignificantly OWN regulated proteins

write.csv(limma_results_sig_up_regulated_with_adjustment, "final_project_data/results_limma_sig_up_regulated_with_adjustment.csv")
write.csv(limma_results_sig_down_regulated_with_adjustment, "final_project_data/results_limma_sig_down_regulated_with_adjustment.csv")

### Kaplan-Meier Plot Prep (code to run is with other KM Plots) ###

limma_results_sig_up_regulated_with_adjustment_sorted <- limma_results_sig_up_regulated_with_adjustment[order(limma_results_sig_up_regulated_with_adjustment$logFC, decreasing = TRUE), ] #sorts significantly up-regulated results by log2FoldChange
limma_results_sig_down_regulated_with_adjustment_sorted <- limma_results_sig_down_regulated_with_adjustment[order(limma_results_sig_down_regulated_with_adjustment$logFC), ] #sorts significantly down-regulated results by log2FoldChange

limma_results_proteins_kaplan_meier <- rownames(limma_results_sig_up_regulated_with_adjustment_sorted) + rownames(limma_results_sig_down_regulated_with_adjustment_sorted)

