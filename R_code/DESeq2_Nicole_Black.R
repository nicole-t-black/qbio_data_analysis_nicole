if (!require(DESeq2)) BiocManager::install("DESeq2")

library(TCGAbiolinks)
library(DESeq2)

# REFERENCE: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# If testing on the 6 patients, get the patient barcodes
# clinical_file <- read.csv("../data/tcga_brca_six_example_clinical.csv")
# barcodes <- as.character( clinical_file$barcode )

###### Load in your HTSeq Counts #######
#Option A: Use GDCquery #see below, loaded HTSeq Counts into sum_exp
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

#patient_data_sample <- head(patient_data)
#write.csv(patient_data_sample, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/patient_data_sample.csv")

##### Preprocess your data #####
#How many genes are in counts? 56,602 genes

counts <- counts[rowMeans(counts) >= 10, ] #rewrites counts with counts data for genes where mean >= 10 is TRUE

counts <- counts[ , patients_no_NA_mask] #rewrites counts with only patients where age (patients_no_NA_mask) is TRUE

#counts_sample <- counts[1:20, ]
#write.csv(counts_sample, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/counts_sample.csv")

#We need to add an age_category column to our patient data
patient_data$age_category = ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < 40, "Young", 
                             ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age_category column, categorized by young, mid, old

#Next, we need to make age_category a "factor".
patient_data$age_category <- factor(patient_data$age_category, levels=c("Young", "Mid", "Old")) #same column (age_category) but now a categorical variable, with levels Young Mid Old

####### Now for actual analysis!! #######
dds <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~age_category) #creates a DESeqDataSet object from the counts and patient_data matrices, using age_category as the condition
dds_obj <- DESeq(dds) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values
resultsNames(dds_obj) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young" (why no Mid vs Old?)

results <- results(dds_obj, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values
  #contrast=c specifies the comparison for the fold change, with age_category as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

head(results) #look at the results

#Notice, each gene has a log2FoldChange and a padj value. This is what we are interested in!
#For clarification, please add a FoldChange column by computing 2^log2FoldChange column
results$FoldChange <- 2^results$log2FoldChange #creates column in results with FoldChange

#Save ALL your results to a csv file
write.csv(results, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/DESeq_Data.csv")

####### Interpreting results ########

#We often visualize results via a "volcano plot"
padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/DESeq_Volcano_Plot.jpg")
plot(x= results$log2FoldChange, y= -log10(results$padj) ) #plots significance (-log10(padj)) vs expression (log2FoldChange)
#abline() plots straight lines on an R plot.
#v argument is for a vertical line, h argument is for a horizontal line. col argument is color
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

######## Look at your volcano plot and answer the following questions ########

#What does each dot on the plot represent?
      #One gene, it's expression
#Why might we have two vertical line and only one horizontal line?
      #Vertical lines represent cutoffs for UP or DOWN expression, horizontal is for significance)
#Why are we plotting the -log10 of the adjusted p values rather than the actual adjusted p values?

#padj not -log10 plot
padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/DESeq_Volcano_Plot_Padj_Test.jpg")
plot(x= results$log2FoldChange, y= results$padj ) #plots significance (padj) vs expression (log2FoldChange)
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(padj_threshold), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

#We want to separate between genes that are UP regulated in young (higher expression in young patients),
#                                         and genes that are DOWN regulated in young (lower expression in young patients)
#If the log2FoldChange is POSITIVE, the expression is higher in young
#If the log2FoldChange is NEGATIVE, the expression is lower in young and higher in old
#What does the log2FoldChange equal, if the expression is the same in young and old patients?
      #Log2FoldChange = 0, or near to 0

results_significant_padj <- results[results$padj < padj_threshold, ] #filters for only significant (padj < 0.05) genes, creates table

results_sig_up_regulated <- results_significant_padj[results_significant_padj$log2FoldChange > log2FC_threshold, ] #UP regulated significant genes
results_sig_down_regulated <- results_significant_padj[results_significant_padj$log2FoldChange < -log2FC_threshold, ] #DOWN regulated significant genes

#How could you get the same results using the absolute value of the log2FoldChange?
      #Not sure, I get that you can easily get the UP regulated, but not down...?

#Notice that the gene names are in the ENSG00000#### format. This is the ensembl_gene_id format.
#we probably want the "common" name of the gene.
gene_information <- rowData(sum_exp)

results_sig_up_regulated$CommonGeneName <- gene_information[rownames(results_sig_up_regulated), 2] #adds column with common gene names to results_sig_up_regulated
results_sig_down_regulated$CommonGeneName <- gene_information[rownames(results_sig_down_regulated), 2] #adds column with common gene names to results_sig_down_regulated

write.csv(results_sig_up_regulated, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/results_sig_up_regulated.csv") #saves results_sig_up_regulated into .csv file
write.csv(results_sig_down_regulated, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/results_sig_down_regulated.csv")#saves results_sig_down_regulated into .csv file

##################################################################
#As we have touched on, there are many other variables that may influence the results between young and old
#You may want to ADJUST for this confounding variables like PAM50 subtype, and what is called "Histology"
#Breast cancer subtypes and histology (ductal vs. lobular, feel free to Google!), affect the patient's "omics" data.
#Repeat the above analysis, but add the below modifications.
#Compare your results with adjustment and without adjustment.
#Are there greater or fewer genes significant with the adjustment?
#Are the genes that are signficant the same?

#check that all variables are not "NA"
patients_no_NA_mask <- ( !is.na(colData(sum_exp)$paper_age_at_initial_pathologic_diagnosis)
                        & !is.na(colData(sum_exp)$paper_BRCA_Subtype_PAM50)
                        & !is.na(colData(sum_exp)$paper_BRCA_Pathology)
                        & !colData(sum_exp)$paper_BRCA_Pathology == "NA" )

patient_data <- colData(sum_exp)[ patients_no_NA_mask, ] #recreates patient_data data frame with the sum_exp data only for TRUE age, BRCA pathology, and BRCA subtypes

counts <- assays(sum_exp)$"HTSeq - Counts" #manually remove existing counts in console, then use this to recreate full counts array (?) with number of counts for each gene

counts <- counts[rowMeans(counts) >= 10, ] #rewrites counts with counts data for genes where mean >= 10 is TRUE
counts <- counts[ , patients_no_NA_mask ] #rewrites counts to remove all NA patients (NA for age, BRCA subtype, BRCA pathology)

patient_data$age_category = ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < 40, "Young", 
                                   ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age_category column, categorized by young, mid, old

#all columns must be FACTORS
patient_data$age_category <- factor( patient_data$age_category, levels=c("Young", "Mid", "Old") )
patient_data$paper_BRCA_Subtype_PAM50 <- factor( patient_data$paper_BRCA_Subtype_PAM50, levels=c("Her2","LumA","LumB","Basal","Normal") )
patient_data$paper_BRCA_Pathology <- factor( patient_data$paper_BRCA_Pathology, levels=c("IDC","Other","Mixed","ILC") )

####### Now for actual analysis part 2!! #######
dds_with_adjustment <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~paper_BRCA_Pathology+ paper_BRCA_Subtype_PAM50 +age_category) #creates a DESeqDataSet object from the counts and patient_data matrices, using BRCA_pathology, BRCA_subtype, and age_category as the conditions
dds_obj_with_adjustment <- DESeq(dds_with_adjustment) #runs DESeq on the DESeqDataSet object, returning results tables with log^2 fold, padj, etc. values
resultsNames(dds_obj_with_adjustment) #lists the coefficients, "intercept", "age_category_Mid_vs_Young", "age_category_Old_vs_Young" (why no Mid vs Old?)

results_with_adjustment <- results(dds_obj_with_adjustment, contrast=c("age_category", "Young",'Old')) #extracts analysis table with log^2 fold changes, standard errors, test stats, p-vales, and adjusted p-values
#contrast=c specifies the comparison for the fold change, with age_category as the name of the factor, "Young" as the numerator, and "Old" as the denominator"

head(results_with_adjustment) #look at the results

#Notice, each gene has a log2FoldChange and a padj value. This is what we are interested in!
#For clarification, please add a FoldChange column by computing 2^log2FoldChange column
results_with_adjustment$FoldChange <- 2^results_with_adjustment$log2FoldChange #creates column in results with FoldChange

#Save ALL your results to a csv file
write.csv(results_with_adjustment, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/DESeq_Data_with_Adjustment.csv")

####### Interpreting results ########

#We often visualize results via a "volcano plot"
padj_threshold <- 0.05 #significance threshold
log2FC_threshold <- 1.0 #differential expression threshold
jpeg("/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/DESeq_Volcano_Plot_with_Adjustment.jpg")
plot(x= results_with_adjustment$log2FoldChange, y= -log10(results_with_adjustment$padj) ) #plots significance (-log10(padj)) vs expression (log2FoldChange)
#abline() plots straight lines on an R plot.
#v argument is for a vertical line, h argument is for a horizontal line. col argument is color
abline(v=c(log2FC_threshold, -log2FC_threshold), h= c(-log10(padj_threshold)), col="green") #plots vertical lines at +/- 1 expression and horizontal at significance cutoff of 0.05
dev.off()

#stuck because it says there's an NA??? but I can't find it :(
results_with_adjustment_significant_adjp <- results_with_adjustment[results_with_adjustment$padj > padj_threshold, ] #filters for only significant (padj > 0.05) genes, creates table

results_with_adjustment_sig_up_regulated <- results_with_adjustment_significant_adjp[results_with_adjustment_significant_adjp$log2FoldChange > log2FC_threshold, ] #UP regulated significant genes
results_with_adjustment_sig_down_regulated <- results_with_adjustment_significant_adjp[results_with_adjustment_significant_adjp$log2FoldChange < -log2FC_threshold, ] #DOWN regulated significant genes

gene_information <- rowData(sum_exp)

results_with_adjustment_sig_up_regulated$CommonGeneName <- gene_information[rownames(results_with_adjustment_sig_up_regulated), 2] #adds column with common gene names to results_sig_up_regulated
results_with_adjustment_sig_down_regulated$CommonGeneName <- gene_information[rownames(results_with_adjustment_sig_down_regulated), 2] #adds column with common gene names to results_sig_down_regulated

write.csv(results_with_adjustment_sig_up_regulated, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/results_with_adjustment_sig_up_regulated.csv") #saves results_sig_up_regulated into .csv file
write.csv(results_with_adjustment_sig_down_regulated, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/data/results_with_adjustment_sig_down_regulated.csv")#saves results_sig_down_regulated into .csv file

###Kaplan-Meier Plots##

#CSN3#
gene_counts_CSN3 <- counts["ENSG00000171209",]
gene_counts_CSN3 <- sort(gene_counts_CSN3, decreasing = TRUE, na.last = NA)
length_gene_CSN3 <- length(gene_counts_CSN3)
high_bound_gene_CSN3 <- gene_counts_CSN3[length_gene_CSN3/3] #counts of gene at border of high and mid
mid_bound_gene_CSN3 <- gene_counts_CSN3[length_gene_CSN3-length_gene_CSN3/3] #counts of gene at border of mid and low
patient_data$gene_level_CSN3 <- ifelse(counts["ENSG00000171209",]>=high_bound_gene_CSN3 ,"High", ifelse(counts["ENSG00000171209",] < mid_bound_gene_CSN3, "Low", "Mid"))

TCGAanalyze_survival( patient_data, "gene_level_CSN3", filename="survival_expression_of_CSN3_Thirds.pdf")

CSN3_average <- sum(counts["ENSG00000171209", ])/ncol(counts)
patient_data$CSN3_expression = ifelse(counts["ENSG00000171209", patient_data$barcode] < CSN3_average, "Low", "High") #adds age_category column, categorized by young, mid, old

TCGAanalyze_survival( patient_data, "CSN3_expression", filename="survival_expression_of_CSN3.pdf")

#CSN2#
gene_counts_CSN2 <- counts["ENSG00000135222",]
gene_counts_CSN2 <- sort(gene_counts_CSN2, decreasing = TRUE, na.last = NA)
length_gene_CSN2 <- length(gene_counts_CSN2)
high_bound_gene_CSN2 <- gene_counts_CSN2[length_gene_CSN2/3] #counts of gene at border of high and mid
mid_bound_gene_CSN2 <- gene_counts_CSN2[length_gene_CSN2-length_gene_CSN2/3] #counts of gene at border of mid and low
patient_data$gene_level_CSN2 <- ifelse(counts["ENSG00000135222",]>=high_bound_gene_CSN3 ,"High", ifelse(counts["ENSG00000135222",] < mid_bound_gene_CSN2, "Low", "Mid"))

TCGAanalyze_survival( patient_data, "gene_level_CSN2", filename="survival_expression_of_CSN2_Thirds.pdf")

CSN2_average <- sum(counts["ENSG00000135222", ])/ncol(counts)
patient_data$CSN2_expression = ifelse(counts["ENSG00000135222", patient_data$barcode] < CSN2_average, "Low", "High") #adds age_category column, categorized by young, mid, old

TCGAanalyze_survival( patient_data, "CSN2_expression", filename="survival_expression_of_CSN2.pdf")

#CSN1S1#
gene_counts_CSN1S1 <- counts["ENSG00000126545",]
gene_counts_CSN1S1 <- sort(gene_counts_CSN1S1, decreasing = TRUE, na.last = NA)
length_gene_CSN1S1 <- length(gene_counts_CSN1S1)
high_bound_gene_CSN1S1 <- gene_counts_CSN1S1[length_gene_CSN1S1/3] #counts of gene at border of high and mid
mid_bound_gene_CSN1S1 <- gene_counts_CSN1S1[length_gene_CSN1S1-length_gene_CSN1S1/3] #counts of gene at border of mid and low
patient_data$gene_level_CSN1S1 <- ifelse(counts["ENSG00000126545",]>=high_bound_gene_CSN1S1 ,"High", ifelse(counts["ENSG00000126545",] < mid_bound_gene_CSN1S1, "Low", "Mid"))

TCGAanalyze_survival( patient_data, "gene_level_CSN1S1", filename="survival_expression_of_CSN1S1_Thirds.pdf")

CSN1S1_average <- sum(counts["ENSG00000126545", ])/ncol(counts)
patient_data$CSN1S1_expression = ifelse(counts["ENSG00000126545", patient_data$barcode] < CSN1S1_average, "Low", "High") #adds age_category column, categorized by young, mid, old

TCGAanalyze_survival( patient_data, "CSN1S1_expression", filename="survival_expression_of_CSN1S1.pdf")

#CARTPT#
gene_counts_CARTPT <- counts["ENSG00000164326",]
gene_counts_CARTPT <- sort(gene_counts_CARTPT, decreasing = TRUE, na.last = NA)
length_gene_CARTPT <- length(gene_counts_CARTPT)
high_bound_gene_CARTPT <- gene_counts_CARTPT[length_gene_CARTPT/3] #counts of gene at border of high and mid
mid_bound_gene_CARTPT <- gene_counts_CARTPT[length_gene_CARTPT-length_gene_CARTPT/3] #counts of gene at border of mid and low
patient_data$gene_level_CARTPT <- ifelse(counts["ENSG00000164326",]>=high_bound_gene_CARTPT ,"High", ifelse(counts["ENSG00000164326",] < mid_bound_gene_CARTPT, "Low", "Mid"))

TCGAanalyze_survival( patient_data, "gene_level_CARTPT", filename="survival_expression_of_CARTPT_Thirds.pdf")

CARTPT_average <- sum(counts["ENSG00000164326", ])/ncol(counts)
patient_data$CARTPT_expression = ifelse(counts["ENSG00000164326", patient_data$barcode] < CARTPT_average, "Low", "High") #adds age_category column, categorized by young, mid, old

TCGAanalyze_survival( patient_data, "CARTPT_expression", filename="survival_expression_of_CARTPT.pdf")

#LACRT#
gene_counts_LACRT <- counts["ENSG00000135413",]
gene_counts_LACRT <- sort(gene_counts_LACRT, decreasing = TRUE, na.last = NA)
length_gene_LACRT <- length(gene_counts_LACRT)
high_bound_gene_LACRT <- gene_counts_LACRT[length_gene_LACRT/3] #counts of gene at border of high and mid
mid_bound_gene_LACRT <- gene_counts_LACRT[length_gene_LACRT-length_gene_LACRT/3] #counts of gene at border of mid and low
patient_data$gene_level_LACRT <- ifelse(counts["ENSG00000135413",]>=high_bound_gene_LACRT ,"High", ifelse(counts["ENSG00000135413",] < mid_bound_gene_LACRT, "Low", "Mid"))

TCGAanalyze_survival( patient_data, "gene_level_LACRT", filename="survival_expression_of_LACRT_Thirds.pdf")

LACRT_average <- sum(counts["ENSG00000135413", ])/ncol(counts)
patient_data$LACRT_expression = ifelse(counts["ENSG00000135413", patient_data$barcode] < LACRT_average, "Low", "High") #adds age_category column, categorized by young, mid, old

TCGAanalyze_survival( patient_data, "LACRT_expression", filename="survival_expression_of_LACRT.pdf")

##Violin Plots##
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

##Box Plots##

patient_data$CSN3_counts = counts["ENSG00000171209",]

patient_data$CSN3_counts_log = sapply(counts["ENSG00000171209",], log10)

# Boxplots by age
ggplot(CSN3_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for CNS3 by Age Category") + geom_boxplot()

plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$SERPINB13_counts)

