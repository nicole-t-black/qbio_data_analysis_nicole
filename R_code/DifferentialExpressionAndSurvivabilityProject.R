if (!require(DESeq2)) BiocManager::install("DESeq2")
install.packages('matrixStats')
library(matrixStats)
library(TCGAbiolinks)
library(DESeq2)

# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

### Load in HTSeq Counts ###

 library(SummarizedExperiment)
 query <- GDCquery(project = "TCGA-BRCA",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts")
 GDCdownload(query)
 sum_exp <- GDCprepare(query)
 
#access the actual counts. counts is genes in rows X patients in columns.
counts <- assays(sum_exp)$"HTSeq - Counts" #creates counts array (?) with number of counts for each gene

#We need to remove patients with unknown ages.
patients_no_NA_mask <- !is.na(colData(sum_exp)$paper_age_at_initial_pathologic_diagnosis) #FALSE if NA, TRUE if has age information

#access the patient_data from coldata
patient_data <- colData(sum_exp)[ patients_no_NA_mask, ] #creates patient_data data frame with the sum_exp data only for TRUE age patients

#### Pre-processing ####

counts <- counts[rowMeans(counts) >= 10, patients_no_NA_mask] #rewrites counts with counts data for genes where mean >= 10 is TRUE and age (patients_no_NA_mask) is TRUE

patient_age <- patient_data[ ,"paper_age_at_initial_pathologic_diagnosis"]
median_age <- median(patient_age)

patient_data$age_category = ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < 40, "Young", 
                             ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age_category column, categorized by young, mid, old

patient_data$age_category_median = ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < median_age, "Young", "Old") #adds age_category column separated by the median age

patient_data$age_category <- factor(patient_data$age_category, levels=c("Young", "Mid", "Old")) #same column (age_category) but now a categorical variable, with levels Young Mid Old

patient_data$age_category_median <- factor(patient_data$age_category_median, levels=c("Young", "Old")) #same column (age_category) but now a categorical variable, with levels Young Old

###DESeq2 with young/old (ignoring mid)###

dds1 <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~age_category) #creates a DESeqDataSet object from the counts and patient_data matrices, using age_category as the condition
dds1_obj <- DESeq(dds1) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values
resultsNames(dds1_obj) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young" (why no Mid vs Old?)

results1 <- results(dds1_obj, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values
  #contrast=c specifies the comparison for the fold change, with age_category as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

head(results1) #view results

results1$FoldChange <- 2^results1$log2FoldChange #creates column in results with FoldChange

write.csv(results1, "../data/DESeq_Data_1.csv")

##Volcano plot with young/old (ignoring mid)##

padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("../data/DESeq_Volcano_Plot_1.jpg")
plot(x= results1$log2FoldChange, y= -log10(results1$padj), main= "Differentially Expressed Genes between Young and Old Patients", xlab= "log2(Fold Change)", ylab="-log10(p-adjusted value)") #plots significance (-log10(padj)) vs expression (log2FoldChange)
#abline() plots straight lines on an R plot.
#v argument is for a vertical line, h argument is for a horizontal line. col argument is color
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

##Tables of Up/Down Regulated Genes with young/old (ignoring mid)##

results_significant_padj_1 <- results1[results1$padj < padj_threshold, ] #filters for only significant (padj < 0.05) genes, creates table

results_sig_up_regulated_1 <- results_significant_padj_1[results_significant_padj_1$log2FoldChange > log2FC_threshold, ] #UP regulated significant genes
results_sig_down_regulated_1 <- results_significant_padj_1[results_significant_padj_1$log2FoldChange < -log2FC_threshold, ] #DOWN regulated significant genes


gene_information <- rowData(sum_exp)

results_sig_up_regulated_1$CommonGeneName <- gene_information[rownames(results_sig_up_regulated_1), 2] #adds column with common gene names to results_sig_up_regulated
results_sig_down_regulated_1$CommonGeneName <- gene_information[rownames(results_sig_down_regulated_1), 2] #adds column with common gene names to results_sig_down_regulated

write.csv(results_sig_up_regulated_1, "../data/results_sig_up_regulated_1.csv") #saves results_sig_up_regulated into .csv file
write.csv(results_sig_down_regulated_1, "../data/results_sig_down_regulated_1.csv")#saves results_sig_down_regulated into .csv file


###DESeq2 with young old defined by median###

dds2 <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~age_category_median) #creates a DESeqDataSet object from the counts and patient_data matrices, using age_category as the condition
dds2_obj <- DESeq(dds2) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values
resultsNames(dds2_obj) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young" (why no Mid vs Old?)

results2 <- results(dds2_obj, contrast=c("age_category_median", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values
#contrast=c specifies the comparison for the fold change, with age_category_median as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

head(results2) #view results

results2$FoldChange <- 2^results2$log2FoldChange #creates column in results with FoldChange

write.csv(results2, "../data/DESeq_Data_2.csv")

##Volcano plot with young old defined by median##

padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("../data/DESeq_Volcano_Plot_2.jpg")
plot(x= results2$log2FoldChange, y= -log10(results2$padj), main= "Differentially Expressed Genes between Young (< 56 y/o) and Old (>= 56 y/o) Patients", xlab= "log2(Fold Change)", ylab="-log10(p-adjusted value)") #plots significance (-log10(padj)) vs expression (log2FoldChange)
#abline() plots straight lines on an R plot.
#v argument is for a vertical line, h argument is for a horizontal line. col argument is color
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

##Tables of Up/Down Regulated Genes with young old defined by median##

results_significant_padj_2 <- results2[results2$padj < padj_threshold, ] #filters for only significant (padj < 0.05) genes, creates table

results_sig_up_regulated_2 <- results_significant_padj_2[results_significant_padj_2$log2FoldChange > log2FC_threshold, ] #UP regulated significant genes
results_sig_down_regulated_2 <- results_significant_padj_2[results_significant_padj_2$log2FoldChange < -log2FC_threshold, ] #DOWN regulated significant genes


gene_information <- rowData(sum_exp)

results_sig_up_regulated_2$CommonGeneName <- gene_information[rownames(results_sig_up_regulated_2), 2] #adds column with common gene names to results_sig_up_regulated
results_sig_down_regulated_2$CommonGeneName <- gene_information[rownames(results_sig_down_regulated_2), 2] #adds column with common gene names to results_sig_down_regulated

write.csv(results_sig_up_regulated_2, "../data/results_sig_up_regulated_2.csv") #saves results_sig_up_regulated into .csv file
write.csv(results_sig_down_regulated_2, "../data/results_sig_down_regulated_2.csv")#saves results_sig_down_regulated into .csv file


###Kaplan-Meier Plots###

#CSN3#

#gene_counts_CSN3 <- counts["ENSG00000171209",]
#gene_counts_CSN3 <- sort(gene_counts_CSN3, decreasing = TRUE, na.last = NA)
#length_gene_CSN3 <- length(gene_counts_CSN3)
#high_bound_gene_CSN3 <- gene_counts_CSN3[length_gene_CSN3/3] #counts of gene at border of high and mid
#mid_bound_gene_CSN3 <- gene_counts_CSN3[length_gene_CSN3-length_gene_CSN3/3] #counts of gene at border of mid and low
#patient_data$gene_level_CSN3 <- ifelse(counts["ENSG00000171209",]>=high_bound_gene_CSN3 ,"High", ifelse(counts["ENSG00000171209",] < mid_bound_gene_CSN3, "Low", "Mid"))

#TCGAanalyze_survival( patient_data, "gene_level_CSN3", legend="CSN3 Expression Level", filename="../data/survival_expression_of_CSN3_Thirds.pdf")

CSN3_zero_mask <- counts["ENSG00000171209", ] == 0
CSN3_zero_count <- sum(CSN3_zero_mask)
CSN3_average <- sum(counts["ENSG00000171209", ])/(ncol(counts)-CSN3_zero_count)

CSN3_median <- rowMedians(counts, rows = "ENSG00000171209")
patient_data$CSN3_expression = ifelse(counts["ENSG00000171209", patient_data$barcode] > CSN3_average, "High", ifelse(counts["ENSG00000171209", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "CSN3_expression", legend="CSN3 Expression Level", filename="../data/survival_expression_of_CSN3.pdf")
TCGAanalyze_survival( patient_data, "CSN3_expression", legend="CSN3 Expression Level", filename="../data/survival_expression_of_CSN3.jpg")

#CSN2#
CSN2_zero_mask <- counts["ENSG00000135222", ] == 0
CSN2_zero_count <- sum(CSN2_zero_mask)
CSN2_average <- sum(counts["ENSG00000135222", ])/(ncol(counts)-CSN2_zero_count)
patient_data$CSN2_expression = ifelse(counts["ENSG00000135222", patient_data$barcode] > CSN2_average, "High", ifelse(counts["ENSG00000135222", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "CSN2_expression", legend="CSN2 Expression Level", filename="../data/survival_expression_of_CSN2.pdf")
TCGAanalyze_survival( patient_data, "CSN2_expression", legend="CSN2 Expression Level", filename="../data/survival_expression_of_CSN2.jpg")

#CSN1S1#
CSN1S1_zero_mask <- counts["ENSG00000126545", ] == 0
CSN1S1_zero_count <- sum(CSN1S1_zero_mask)
CSN1S1_average <- sum(counts["ENSG00000126545", ])/(ncol(counts)-CSN1S1_zero_count)
patient_data$CSN1S1_expression = ifelse(counts["ENSG00000126545", patient_data$barcode] > CSN1S1_average, "High", ifelse(counts["ENSG00000126545", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "CSN1S1_expression", legend="CSN1S1 Expression Level", filename="../data/survival_expression_of_CSN1S1.pdf")
TCGAanalyze_survival( patient_data, "CSN1S1_expression", legend="CSN1S1 Expression Level", filename="../data/survival_expression_of_CSN1S1.jpg")

#CARTPT#
CARTPT_zero_mask <- counts["ENSG00000164326", ] == 0
CARTPT_zero_count <- sum(CARTPT_zero_mask)
CARTPT_average <- sum(counts["ENSG00000164326", ])/(ncol(counts)-CARTPT_zero_count)
patient_data$CARTPT_expression = ifelse(counts["ENSG00000164326", patient_data$barcode] > CARTPT_average, "High", ifelse(counts["ENSG00000164326", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "CARTPT_expression", legend="CARTPT Expression Level", filename="../data/survival_expression_of_CARTPT.pdf")
TCGAanalyze_survival( patient_data, "CARTPT_expression", legend="CARTPT Expression Level", filename="../data/survival_expression_of_CARTPT.jpg")

#LACRT#
LACRT_zero_mask <- counts["ENSG00000135413", ] == 0
LACRT_zero_count <- sum(LACRT_zero_mask)
LACRT_average <- sum(counts["ENSG00000135413", ])/(ncol(counts)-LACRT_zero_count)
patient_data$LACRT_expression = ifelse(counts["ENSG00000135413", patient_data$barcode] > LACRT_average, "High", ifelse(counts["ENSG00000135413", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "LACRT_expression", legend="LACRT Expression Level", filename="../data/survival_expression_of_LACRT.pdf")
TCGAanalyze_survival( patient_data, "LACRT_expression", legend="LACRT Expression Level", filename="../data/survival_expression_of_LACRT.jpg")

#CHGB#
CHGB_zero_mask <- counts["ENSG00000089199", ] == 0
CHGB_zero_count <- sum(CHGB_zero_mask)
CHGB_average <- sum(counts["ENSG00000089199", ])/(ncol(counts)-CHGB_zero_count)
patient_data$CHGB_expression = ifelse(counts["ENSG00000089199", patient_data$barcode] > CHGB_average, "High", ifelse(counts["ENSG00000089199", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "CHGB_expression", legend="CHGB Expression Level", filename="../data/survival_expression_of_CHGB.pdf")
TCGAanalyze_survival( patient_data, "CHGB_expression", legend="CHGB Expression Level", filename="../data/survival_expression_of_CHGB.jpg")

#PRSS2#
PRSS2_zero_mask <- counts["ENSG00000275896", ] == 0
PRSS2_zero_count <- sum(PRSS2_zero_mask)
PRSS2_average <- sum(counts["ENSG00000275896", ])/(ncol(counts)-PRSS2_zero_count)
patient_data$PRSS2_expression = ifelse(counts["ENSG00000275896", patient_data$barcode] > PRSS2_average, "High", ifelse(counts["ENSG00000275896", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "PRSS2_expression", legend="PRSS2 Expression Level", filename="../data/survival_expression_of_PRSS2.pdf")
TCGAanalyze_survival( patient_data, "PRSS2_expression", legend="PRSS2 Expression Level", filename="../data/survival_expression_of_PRSS2.jpg")

#MAGEA3#
MAGEA3_zero_mask <- counts["ENSG00000221867", ] == 0
MAGEA3_zero_count <- sum(MAGEA3_zero_mask)
MAGEA3_average <- sum(counts["ENSG00000221867", ])/(ncol(counts)-MAGEA3_zero_count)
patient_data$MAGEA3_expression = ifelse(counts["ENSG00000221867", patient_data$barcode] > MAGEA3_average, "High", ifelse(counts["ENSG00000221867", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "MAGEA3_expression", legend="MAGEA3 Expression Level", filename="../data/survival_expression_of_MAGEA3.pdf")
TCGAanalyze_survival( patient_data, "MAGEA3_expression", legend="MAGEA3 Expression Level", filename="../data/survival_expression_of_MAGEA3.jpg")

#CYP2A7#
CYP2A7_zero_mask <- counts["ENSG00000198077", ] == 0
CYP2A7_zero_count <- sum(CYP2A7_zero_mask)
CYP2A7_average <- sum(counts["ENSG00000198077", ])/(ncol(counts)-CYP2A7_zero_count)
patient_data$CYP2A7_expression = ifelse(counts["ENSG00000198077", patient_data$barcode] > CYP2A7_average, "High", ifelse(counts["ENSG00000198077", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "CYP2A7_expression", legend="CYP2A7 Expression Level", filename="../data/survival_expression_of_CYP2A7.pdf")
TCGAanalyze_survival( patient_data, "CYP2A7_expression", legend="CYP2A7 Expression Level", filename="../data/survival_expression_of_CYP2A7.jpg")

#FDCSP#
FDCSP_zero_mask <- counts["ENSG00000181617", ] == 0
FDCSP_zero_count <- sum(FDCSP_zero_mask)
FDCSP_average <- sum(counts["ENSG00000181617", ])/(ncol(counts)-FDCSP_zero_count)
patient_data$FDCSP_expression = ifelse(counts["ENSG00000181617", patient_data$barcode] > FDCSP_average, "High", ifelse(counts["ENSG00000181617", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "FDCSP_expression", legend="FDCSP Expression Level", filename="../data/survival_expression_of_FDCSP.pdf")
TCGAanalyze_survival( patient_data, "FDCSP_expression", legend="FDCSP Expression Level", filename="../data/survival_expression_of_FDCSP.jpg")

#BRCA1#
BRCA1_zero_mask <- counts["ENSG00000012048", ] == 0
BRCA1_zero_count <- sum(BRCA1_zero_mask)
BRCA1_average <- sum(counts["ENSG00000012048", ])/(ncol(counts)-BRCA1_zero_count)
patient_data$BRCA1_expression = ifelse(counts["ENSG00000012048", patient_data$barcode] > BRCA1_average, "High", ifelse(counts["ENSG00000012048", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "BRCA1_expression", legend="BRCA1 Expression Level", filename="../data/survival_expression_of_BRCA1.pdf")
TCGAanalyze_survival( patient_data, "BRCA1_expression", legend="BRCA1 Expression Level", filename="../data/survival_expression_of_BRCA1.jpg")

#TP53#
TP53_zero_mask <- counts["ENSG00000141510", ] == 0
TP53_zero_count <- sum(TP53_zero_mask)
TP53_average <- sum(counts["ENSG00000141510", ])/(ncol(counts)-TP53_zero_count)
patient_data$TP53_expression = ifelse(counts["ENSG00000141510", patient_data$barcode] > TP53_average, "High", ifelse(counts["ENSG00000141510", patient_data$barcode] == 0, "No Expression", "Low"))

TCGAanalyze_survival( patient_data, "TP53_expression", legend="TP53 Expression Level", filename="../data/survival_expression_of_TP53.pdf")
TCGAanalyze_survival( patient_data, "TP53_expression", legend="TP53 Expression Level", filename="../data/survival_expression_of_TP53.jpg")

##Violin Plots - Code Does Not Run##
counts_transpose <- t(counts)

counts_transpose$barcodes <- rownames(counts_transpose)

lst <- list(counts_transpose)
do.call(rbind,lst)

counts_transpose$age_category <- patient_data[counts_transpose$barcodes, age_category] #adds age_category column, categorized by young, mid, old
counts$age_category <- patient_data$age_category
patient_data$CSN3_counts = counts["ENSG00000171209", patient_data$barcode]

library(ggplot2)

ggplot(data=patient_data, aes(x=age_category, y=CSN3_counts), fill=age_category) + 
  geom_violin()

##Box Plots - Code Does Not Run##

patient_data$CSN3_counts = counts["ENSG00000171209",]

patient_data$CSN3_counts_log = sapply(counts["ENSG00000171209",], log10)

# Boxplots by age
ggplot(CSN3_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for CNS3 by Age Category") + geom_boxplot()

plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$SERPINB13_counts)

