if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks) #load TCGAbiolinks package

#add barcodes argument to query if you want to run on your local machine for smaller files downloaded
#barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
 #         "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
  #            "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
#barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
 #                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")


#######    Group 1: RNASeq     ############
 library(SummarizedExperiment)
 query <- GDCquery(project = "TCGA-BRCA",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts")
 GDCdownload(query) #only need this line of code once to download the data
 sum_exp <- GDCprepare(query)
# Create a tutorial on SummarizedExperiment

# Boxplots by age
# Add a new column to colData called "age_category"
# If age_at_initial_pathologic_diagnosis is < 40, define patient as "Young" in new column
# If age_at_initial_pathologic_diagnosis is >= 60, define patient as "Old" in new column
# Other patients (between 40 and 60), define as "Mid"
# Choose 3 genes of interest from the paper presentations
# Create 3 different boxplots with age_category on x-axis, counts on the y-axis by repeating the below code for each gene
# remember to give your plot a title and informative axis labels
# png("boxplot_genename.png")
# boxplot(FILL IN HERE)
# *Feel free to add lines here that format your boxplot*
# dev.off()
# Use rsync to copy your create pngfile to local machine for viewing

#######    Group 2: clinical   ###########
 install.packages("survival") #install survival package (only need to do once)
 install.packages("survminer") #install survminer package (only need to do once)
 install.packages("arsenal") #install arsenal package (only need to do once)
 library(survival) #load survival package
 library(survminer) #load survminer package
 library(arsenal) #load arsenal package
# library(dplyr) #unsure if I need this, or what exactly it does
 clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")  #searches GDC for clinical TCGA-BRCA data and loads to clin_query
 GDCdownload( clin_query ) #downloads clin_query data onto local machine, (only need to do once)
 clinic <- GDCprepare_clinic(clin_query, clinical.info="patient") #prepares clin_query data into SE
 names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up" #formatting change in days_to_last_follow_up column

 clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 40, "Young", 
                              ifelse(clinic$age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age_category column, categorized by young, mid, old

 clinic$pathologic_stage_numeric = 
          ifelse(clinic$stage_event_pathologic_stage == "Stage I", "Stage I",
          ifelse(clinic$stage_event_pathologic_stage == "Stage IA", "Stage I",
          ifelse(clinic$stage_event_pathologic_stage == "Stage IB", "Stage I",
                 ifelse(clinic$stage_event_pathologic_stage == "Stage II", "Stage II",
                 ifelse(clinic$stage_event_pathologic_stage == "Stage IIA", "Stage II",
                 ifelse(clinic$stage_event_pathologic_stage == "Stage IIB", "Stage II",
                        ifelse(clinic$stage_event_pathologic_stage == "Stage III","Stage III",
                        ifelse(clinic$stage_event_pathologic_stage == "Stage IIIA","Stage III",
                        ifelse(clinic$stage_event_pathologic_stage == "Stage IIIB","Stage III",
                        ifelse(clinic$stage_event_pathologic_stage == "Stage IIIC","Stage III",
                               ifelse(clinic$stage_event_pathologic_stage == "Stage IV", "Stage IV", "Stage X"))))))))))) #adds pathologic_stage_numeric column, categorized by stage (w/o A, B, C, etc.)


 install.packages("tableone") #install package tableone (only need to do once)
 library(tableone) #load package tableone
 clinic_summary <- CreateTableOne(data = clinic) #creates summary of clinic data and stores in clinic_summary
 
 subtypes <- TCGAquery_subtype(tumor = "BRCA") #creates new "subtypes" data table with only BRCA data
 subtypes$age_category = ifelse(subtypes$age_at_initial_pathologic_diagnosis < 40, "Young", 
                                ifelse(subtypes$age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age_category column, categorized by young, mid, old
 
# How many patients are in subtypes vs clinic? Why? #NB: 1087-subtypes, 1174-clinic, unsure why
 
table_arse <- tableby(age_category ~ (pathologic_stage) + (`mRNA Clusters`) + (BRCA_Pathology),
          data=subtypes,
          numeric.test="kwt", cat.test = "chisq", numeric.stats = c("Nmiss", "meansd"),
          total=FALSE) #creates a table of numeric tests, y value: age; x values: pathologic stage, mRNA clusters, BRCA pathology
df <- as.data.frame(summary(table_arse, text=TRUE, pfootnote=TRUE)) #saves above data as a data frame
write.csv(df, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/R_code/table_arse.csv", row.names=FALSE) #writes data frame into .csv file

#creates Kaplan-Meier plots from clinic data, dependent variables: "stage_event_pathologic_stage", "age_category", etc.
TCGAanalyze_survival( clinic, "stage_event_pathologic_stage", filename="pathologic_stage.pdf")
TCGAanalyze_survival( clinic, "age_category", filename="survival_age_category.pdf")
TCGAanalyze_survival( clinic, "pathologic_stage_numeric", filename="pathologic_stage_simple.pdf")
# rsync to copy and view on local computer #NB: don't know how to use rsync

overall_survival <- as.integer( ifelse( is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death) ) #creates overall_survival value, either days to death (if patient is deceased), or days to last follow up (if patient is alive)
clinic$overall_survival <- overall_survival #adds overall_survival column with data from overall_survival values
clinic$death_event <- ifelse(clinic$vital_status == "Alive", 0,1) #adds death_event column where 0 denotes alive and 1 denotes death event
# colnames(clinic) #use if want to visually see affect of ^^
cox_fit <- coxph(Surv(overall_survival, death_event)~age_at_initial_pathologic_diagnosis, data=clinic) #creates cox plot of survival rate dependent on age at diagnosis
jpeg("cox_plot_age_continuous.jpg") #creates jpeg of cox plot
ggadjustedcurves(cox_fit, data=clinic) #NB: unsure what this does, also removed "cox"
dev.off()

#######    Group 3: MAF   ###########
# BiocManager::install("maftools")
# library(maftools)
# mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="") #choose a pipeline
# after running ^^, navigate to the saved csv file. Open the csv file with below Code
# only need to query once. For repeating code, you can just read in the saved dataframe you created
# maf_dataframe <- read.csv("PATH/FILENAME.csv")
# plotmafSummary(maf = maf_dataframe, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# pdf("maf_summary.pdf")
# dev.off()
# Use rsync to copy onto local machine and view
# Create figures from Part 7: https://bioconductor.riken.jp/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#7_Visualization
