if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks) #load library NB

#add barcodes argument to query if you want to run on your local machine for smaller files downloaded
barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
          "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
              "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
                      "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")


#######    Group 1: RNASeq     ############
# library(SummarizedExperiment)
# query <- GDCquery(project = "TCGA-BRCA",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type = "HTSeq - Counts")
# GDCdownload(query) #only need this line of code once to download the data
# sum_exp <- GDCprepare(query)
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
# install.packages("survival") #install packages, one time only NB
# install.packages("survminer") #install packages, one time only NB
# install.packages("arsenal") #install packages, one time only NB
 library(survival) #load packages NB
 library(survminer) #load packages NB
 library(arsenal) #load packages NB
# clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")  #searches GDC and loads data to "clin query", only need to do once NB
# GDCdownload( clin_query ) #should only need this command once. This downloads the files onto your system. #only need to do once NB
# clinic <- GDCprepare_clinic(clin_query, clinical.info="patient") #prepares info from clin query, only need to do once NB
# names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up" #formatting change, I think, only need to do once, NB
 
# clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 40, "Young", ifelse(clinic$age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age category, only need to do once NB
#  clinic$pathologic_stage_numeric = 
 #         ifelse(clinic$stage_event_pathologic_stage == "Stage I", "Stage I",
  #        ifelse(clinic$stage_event_pathologic_stage == "Stage IA", "Stage I",
   #       ifelse(clinic$stage_event_pathologic_stage == "Stage IB", "Stage I",
    #             ifelse(clinic$stage_event_pathologic_stage == "Stage II", "Stage II",
     #            ifelse(clinic$stage_event_pathologic_stage == "Stage IIA", "Stage II",
      #           ifelse(clinic$stage_event_pathologic_stage == "Stage IIB", "Stage II",
       #                 ifelse(clinic$stage_event_pathologic_stage == "Stage III","Stage III",
        #                ifelse(clinic$stage_event_pathologic_stage == "Stage IIIA","Stage III",
         #               ifelse(clinic$stage_event_pathologic_stage == "Stage IIIB","Stage III",
          #              ifelse(clinic$stage_event_pathologic_stage == "Stage IIIC","Stage III",
           #                    ifelse(clinic$stage_event_pathologic_stage == "Stage IV", "Stage IV", "Stage X"))))))))))) #whole mess...like actually please help NB


# install.packages("tableone") #install package one time NB
 library(tableone) #laods package NB
# clinic_summary <- CreateTableOne(data = clinic) #creates summary, once NB
 
# subtypes <- TCGAquery_subtype(tumor = "BRCA") #creates new table, once NB
# subtypes$age_category = ifelse(subtypes$age_at_initial_pathologic_diagnosis < 40, "Young", ifelse(subtypes$age_at_initial_pathologic_diagnosis >= 60, "Old", "Mid")) #adds age category, only need to do once NB
 
# How many patients are in subtypes vs clinic? Why? 1087-subtypes, 1174-clinic
 
# BELOW CODE WILL NOT RUN. This is an example from another dataset
# Consider what you need to change for it to run on clinic
# Hint: age_category stays the same. "Oncotree.Code" was the name of a column
# Option to use either subtypes or clinic dataframe

 
table_arse <- tableby(age_category ~ (pathologic_stage) + (`mRNA Clusters`) + (BRCA_Pathology),
          data=subtypes,
          numeric.test="kwt", cat.test = "fe", numeric.stats = c("Nmiss", "meansd"),
          simulate.p.value = TRUE, total=FALSE
          )
df <- as.data.frame(summary(table_arse, text=TRUE, pfootnote=TRUE))
write.csv(df, "/Users/nicoleblack/Desktop/d/qbio_data_analysis_nicole_local/qbio_data_analysis_nicole/GDCdata/table_arse.csv", row.names=FALSE)

# Modify below code for different variables of interest
# The "gender" can be changed to any of the column names of the clinic dataframe
# Look at the dataframe via the str(clinic) command

TCGAanalyze_survival( clinic, "stage_event_pathologic_stage", filename="pathologic_stage.pdf")
TCGAanalyze_survival( clinic, "age_category", filename="survival_age_category.pdf")
TCGAanalyze_survival( clinic, "age_category", filename="survival_age_category.pdf")
TCGAanalyze_survival( clinic, "pathologic_stage_numeric", filename="pathologic_stage_simple.pdf")
# use ls to confirm the file was created
# rsync to copy and view on local computer

overall_survival <- as.integer( ifelse( is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death) ) #overall survival value, including data from alive patients (assumes death upon last visitation) NB
clinic$overall_survival <- overall_survival #new column with overall survival info NB
clinic$death_event <- ifelse(clinic$vital_status == "Alive", 0,1) #needed to change some formatting issues NB
# colnames(clinic) #use if want to visually see affect of ^^
cox_fit <- coxph(Surv(overall_survival, death_event)~age_at_initial_pathologic_diagnosis, data=clinic)
jpeg("cox_plot_age_continuous.jpg") #can't find it NB
ggadjustedcurves(cox_fit, data=clinic) #what does this do? also removed "cox" NB
dev.off()
# Use rsync to copy figure onto local system and view #how to use rsync NB

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
