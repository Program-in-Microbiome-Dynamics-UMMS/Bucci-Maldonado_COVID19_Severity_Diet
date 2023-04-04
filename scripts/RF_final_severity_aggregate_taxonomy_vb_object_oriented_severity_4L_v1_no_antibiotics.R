rm(list=ls())
library(tune)
library(workflows)
library(tidymodels)
library(tictoc)
library(phyloseq)
library(DESeq2)
library(dplyr)
library(yingtools2)
library(data.table)
#library(ifultools)
library(stringr)
library(ggplot2)
library(gtools)
library(tidyr)
library(ggprism)
library(ranger)
library(caret)
library(logr)

source('./ml_model_class.R')
source('./plot_importance.R')

get_vst_data <- function(phy, var_name, o_relative = F, o_glom = F, glom_level = NULL){
  # not run:
  # uncomment for debug:
  #phy <- phy_stl
  #o_relative <- F
  #o_glom <- F
  #glom_level <- "Genus"
  #var_name = "severity_remdesivir"
  
  ps_stl <- phy
  #ps_stl <- subset_samples(ps_stl, get(var_name) %in% c("1","0"))
  ps_stl <-  prune_taxa(taxa_sums(ps_stl)>0,ps_stl)
  # Total number of samples:
  nsamples(ps_stl)
  stl_sm_dt <- data.frame(sample_data(ps_stl))
  # Now lets focus on the stool samples:
  ps_stl
  
  ### Compute Diversity
  stl_dt <- get.samp(ps_stl,stats = T) %>% as.data.frame()
  phy_deseq <- ps_stl 
  
  prevdf = apply(X = otu_table(phy_deseq),
                 MARGIN = ifelse(taxa_are_rows(phy_deseq), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(phy_deseq),
                      tax_table(phy_deseq)@.Data)
  prevdf$OTU <-  rownames(prevdf)
  prevdf <- prevdf[order(prevdf$Prevalence,decreasing = T), ]
  prevdf$OTU <- factor(as.character(prevdf$OTU) , levels = as.character(prevdf$OTU))
  prevdf$Prop <- prevdf$TotalAbundance/sum(prevdf$TotalAbundance)
  
  prev_frac<-0.05
  prev_cutoff <- prev_frac*nsamples(phy_deseq) # Cut-off
  ab_cutoff <- 1e-4# Cut-off
  # Prevalence
  prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff, ]
  # Abundance
  prevdf_fil <- prevdf_fil[prevdf_fil$Prop >= ab_cutoff, ]
  most_prevalent_taxa<-unique(rownames(prevdf_fil))
  most_prevalent_taxa
  
  phy_fil <- prune_taxa(most_prevalent_taxa,phy_deseq)
  
  if (o_relative == F){
    dds <- phyloseq_to_deseq2(phy_fil, design=formula(paste("~",var_name))) #replace this with any sample variable(s)
    dds <- estimateSizeFactors(dds,"poscounts")
    dds <- estimateDispersions(dds)
    dds <- DESeq(dds,fitType= "local")
    
    vst_dt <- getVarianceStabilizedData(dds)
    dt_meta <-  data.frame(sample_data(phy_fil))
    dt_meta$SampleID  <- rownames(dt_meta)
    dt_meta <- dt_meta[match(colnames(vst_dt),dt_meta$SampleID),]
    
    #make phyloseq object 
    dt_meta$sample  <- rownames(dt_meta)
    rownames(dt_meta) <- NULL
    phy_vst <- phyloseq(otu_table(vst_dt,taxa_are_rows = T),tax_table(phy_fil), set.samp(dt_meta))
    phy_vst
  } else{
    dt_meta <-  data.frame(sample_data(phy_fil))
    dt_meta$SampleID  <- rownames(dt_meta)
    dt_meta <- dt_meta[match(sample_names(phy_fil),dt_meta$SampleID),]
    #make phyloseq object 
    dt_meta$sample  <- rownames(dt_meta)
    rownames(dt_meta) <- NULL
    phy_vst <- transform_sample_counts(phy_fil, function(x) x/sum(x))
    vst_dt <- otu_table(phy_vst)
  }
  
  if (o_glom == T){
    TTABLE <- as.data.frame(tax_table(phy_vst))
    TTABLE$evalue <- NULL
    TTABLE$pident <- NULL
    TMP <- TTABLE$Kingdom
    TTABLE[,which(names(TTABLE)[1:ncol(TTABLE)]!=glom_level)] <- NULL
    TTABLE <- cbind(TMP, TTABLE)
    lnames<-c(names(TMP), names(TTABLE))
    TTABLE <- tax_table(TTABLE)
    rownames(TTABLE) <- taxa_names(phy_vst)
    colnames(TTABLE) <- lnames
    tax_table(phy_vst) <- TTABLE
    phy_vst <- tax_glom(phy_vst, taxrank = glom_level)
    taxa_names(phy_vst) <- tax_table(phy_vst)[,glom_level]
    vst_dt <- otu_table(phy_vst)
  }
  
  X_mic <- data.frame(t(vst_dt))
  #t_data_all <- cbind(X_mic, Y_class_all = ifelse(dt_meta$Severity == "severe","Y","N"))
  t_data_all <- cbind(X_mic, Y_class_all = ifelse(dt_meta[,var_name] == "1","1","0"))
  
  vst_list <- list(phy_vst,t_data_all)
  names(vst_list) <- c("phy_vst","rf_mat")
  return(vst_list)
}

# Create a directory to save figures and tables
mainDir <- "../../results"
subDir <- "RF_final_severity_4L_aggregate_taxonomy_object_oriented_T"
dir.create(file.path(mainDir, subDir,Sys.Date()), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,Sys.Date(),sep="/")


phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Synergistetes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Euryarchaeota" = "#8c8c8c")


o_RELAB <- T

# Import phyloseq object for TNG and STL samples
# Complete observations
# choose taxonomic_level
##### STL samples #####
phy_stl_0 <-  readRDS("../../data/phy_stl_imp_w_antibiotics-classes_clinical_checked-with-new_severity.rds")

stemp <- as.data.frame(sample_data(phy_stl_0))
stemp$Severity <- as.character(stemp$Severity)
stemp$severity_remdesivir <- as.character(stemp$severity_remdesivir)
rbind(stemp$Severity, stemp$severity_remdesivir)
stemp$Severity[stemp$Severity=="severe"] <- "1"  
stemp$Severity[stemp$Severity!="1"] <- "0"  
stemp$Severity <- factor(stemp$Severity)
sample_data(phy_stl_0) <- sample_data(stemp)
phy_stl <- phy_stl_0

table(stemp$Severity, stemp$vancomycin_1)
table(stemp$Severity, stemp$zosyn_1)
table(stemp$Severity, stemp$severity_remdesivir)

# Use the blast assigned name along with SV ID
taxa_names(phy_stl) <- paste0(taxa_names(phy_stl),"_",make.names(as.character(tax_table(phy_stl)[,"Species"])))
phy_stl
stl_tk <- taxa_names(phy_stl)[which(tax_table(phy_stl)[,3]=="Bacteria")]
phy_stl <- prune_taxa(stl_tk, phy_stl)
phy_stl

# only keeping vanco and zosy positive
samp_mod <- sample_data(phy_stl)
samp_mod$sample_type <- NULL
samp_mod$patient_id <- NULL
samp_mod$SampleID <-  NULL
samp_mod$severity_remdesivir <-  NULL
samp_mod$Fatality <-  NULL
samp_mod$Fatality_new <-  NULL
samp_mod$Height <-  NULL
samp_mod$classification_dexamethasone <-  NULL

#remove antibiotics classes
samp_mod <- samp_mod[,-grep("cefepime",names(samp_mod))]
samp_mod <- samp_mod[,-grep("vancomycin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("cefpodoxime",names(samp_mod))]
samp_mod <- samp_mod[,-grep("azithromycin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("ceftazidime",names(samp_mod))]
samp_mod <- samp_mod[,-grep("levofloxacin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("bactrim",names(samp_mod))]
samp_mod <- samp_mod[,-grep("augmentin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("daptomycin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("meropenem",names(samp_mod))]
samp_mod <- samp_mod[,-grep("ceftriaxone",names(samp_mod))]
samp_mod <- samp_mod[,-grep("unasyn",names(samp_mod))]
samp_mod <- samp_mod[,-grep("amikacin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("zosyn",names(samp_mod))]
samp_mod <- samp_mod[,-grep("cefazolin",names(samp_mod))]
names(samp_mod)
sample_data(phy_stl) <- samp_mod
write.csv(samp_mod,paste0(fig_folder,"/severity_metadata_stool.csv"))

list_vst_stl_ASV_relative  <- get_vst_data(phy = phy_stl,var_name = "Severity", o_relative = o_RELAB, o_glom = F, glom_level = "Genus")
list_vst_stl_Species_relative  <- get_vst_data(phy = phy_stl,var_name = "Severity", o_relative = o_RELAB, o_glom = T, glom_level = "Species")
list_vst_stl_Genus_relative  <- get_vst_data(phy = phy_stl,var_name = "Severity", o_relative = o_RELAB, o_glom = T, glom_level = "Genus")
#list_vst_stl_Family_relative  <- get_vst_data(phy = phy_stl,var_name = "severity_remdesivir", o_relative = o_RELAB, o_glom = T, glom_level = "Family")

phy_vst_stl <- list_vst_stl_ASV_relative[["phy_vst"]]
mic_data_stl_ASV <- list_vst_stl_ASV_relative[["rf_mat"]]
mic_data_stl_Species <- list_vst_stl_Species_relative[["rf_mat"]]
mic_data_stl_Genus <- list_vst_stl_Genus_relative[["rf_mat"]]
#mic_data_stl_Family <- list_vst_stl_Family_relative[["rf_mat"]]

# aggregate sets at different taxonomic levels
#mic_data_stl <- cbind(mic_data_stl_ASV[,"Y_class_all"],
#                      mic_data_stl_ASV[,2:ncol(mic_data_stl_ASV)-1], 
#                      mic_data_stl_Species[,2:ncol(mic_data_stl_Species)-1],
#                      mic_data_stl_Genus[,2:ncol(mic_data_stl_Genus)-1],
#                      mic_data_stl_Family[,2:ncol(mic_data_stl_Family)-1])

mic_data_stl <- cbind(mic_data_stl_ASV[,"Y_class_all"],
                      mic_data_stl_ASV[,2:ncol(mic_data_stl_ASV)-1])
names(mic_data_stl)[1] <- "class"
names(mic_data_stl)[2:length(names(mic_data_stl))]<-paste0("mic_",names(mic_data_stl)[2:length(names(mic_data_stl))])

mic_data_stl_Species <- cbind(mic_data_stl_Species[,"Y_class_all"],
                              mic_data_stl_Species[,2:ncol(mic_data_stl_Species)-1])
names(mic_data_stl_Species)[1] <- "class"
names(mic_data_stl_Species)[2:length(names(mic_data_stl_Species))]<-paste0("mic_",names(mic_data_stl_Species)[2:length(names(mic_data_stl_Species))])

mic_data_stl_Genus <- cbind(mic_data_stl_Genus[,"Y_class_all"],
                            mic_data_stl_Genus[,2:ncol(mic_data_stl_Genus)-1])
names(mic_data_stl_Genus)[1] <- "class"
names(mic_data_stl_Genus)[2:length(names(mic_data_stl_Genus))]<-paste0("mic_",names(mic_data_stl_Genus)[2:length(names(mic_data_stl_Genus))])

# Now combine sample data 
samp_stl <-  data.frame(sample_data(phy_vst_stl))
samp_stl$Severity <- NULL
samp_stl$SampleID <- NULL
rownames(samp_stl)
rownames(mic_data_stl)
names(samp_stl) <- paste0('cc_',names(samp_stl))
names(samp_stl)

str(samp_stl)
t_data_stl <- cbind(samp_stl,mic_data_stl)
# Lets see how many rows gets removed if we remove NAs
t_data_stl_comp <- t_data_stl[complete.cases(t_data_stl),]

# cc and microbiome for STL Species
t_data_stl_Species <- cbind(samp_stl,mic_data_stl_Species)
t_data_stl_Species <- t_data_stl_Species[complete.cases(t_data_stl_Species),]

t_data_stl_Genus <- cbind(samp_stl,mic_data_stl_Genus)
t_data_stl_Genus <- t_data_stl_Genus[complete.cases(t_data_stl_Genus),]

#t_data_stl_Family <- cbind(samp_stl,mic_data_stl_Family)
#t_data_stl_Family <- t_data_stl_Family[complete.cases(t_data_stl_Family),]

###### TNG samples #####
phy_tng_0 <-  readRDS("../../data/phy_tng_w_antibiotics-classes_clinical_checked-with-new_severity.rds")

stemp <- as.data.frame(sample_data(phy_tng_0))
stemp$Severity <- as.character(stemp$Severity)
stemp$Severity[stemp$Severity=="severe"] <- "1"  
stemp$Severity[stemp$Severity!="1"] <- "0"  
stemp$Severity <- factor(stemp$Severity)
sample_data(phy_tng_0) <- sample_data(stemp)
phy_tng <- phy_tng_0

table(stemp$Severity, stemp$vancomycin_1)
table(stemp$Severity, stemp$zosyn_1)

# Use the blast assigned name along with SV ID
taxa_names(phy_tng) <- paste0(taxa_names(phy_tng),"_",make.names(as.character(tax_table(phy_tng)[,"Species"])))
phy_tng
tng_tk <- taxa_names(phy_tng)[which(tax_table(phy_tng)[,3]=="Bacteria")]
phy_tng <- prune_taxa(tng_tk, phy_tng)
phy_tng

# only keeping vanco and zosy positive
samp_mod <- sample_data(phy_tng)
samp_mod$sample_type <- NULL
samp_mod$patient_id <- NULL
samp_mod$SampleID <-  NULL
samp_mod$severity_remdesivir <-  NULL
samp_mod$Fatality <-  NULL
samp_mod$Fatality_new <-  NULL
samp_mod$Height <-  NULL
samp_mod$classification_dexamethasone <-  NULL

#remove antibiotics classes
samp_mod <- samp_mod[,-grep("cefepime",names(samp_mod))]
samp_mod <- samp_mod[,-grep("vancomycin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("cefpodoxime",names(samp_mod))]
samp_mod <- samp_mod[,-grep("azithromycin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("ceftazidime",names(samp_mod))]
samp_mod <- samp_mod[,-grep("levofloxacin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("bactrim",names(samp_mod))]
samp_mod <- samp_mod[,-grep("augmentin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("daptomycin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("meropenem",names(samp_mod))]
samp_mod <- samp_mod[,-grep("ceftriaxone",names(samp_mod))]
samp_mod <- samp_mod[,-grep("unasyn",names(samp_mod))]
samp_mod <- samp_mod[,-grep("amikacin",names(samp_mod))]
samp_mod <- samp_mod[,-grep("zosyn",names(samp_mod))]
samp_mod <- samp_mod[,-grep("cefazolin",names(samp_mod))]
names(samp_mod)
sample_data(phy_tng) <- samp_mod
write.csv(samp_mod,paste0(fig_folder,"/severity_metadata_tng.csv"))

list_vst_tng_ASV_relative  <- get_vst_data(phy_tng, var_name = "Severity", o_relative = o_RELAB, o_glom = F, glom_level = "Genus")
list_vst_tng_Species_relative  <- get_vst_data(phy_tng, var_name = "Severity", o_relative = o_RELAB, o_glom = T, glom_level = "Species")
list_vst_tng_Genus_relative  <- get_vst_data(phy_tng, var_name = "Severity", o_relative = o_RELAB, o_glom = T, glom_level = "Genus")
#list_vst_tng_Family_relative  <- get_vst_data(phy_tng, var_name = "severity_remdesivir", o_relative = o_RELAB, o_glom = T, glom_level = "Family")

phy_vst_tng <- list_vst_tng_ASV_relative[["phy_vst"]]
mic_data_tng_ASV <- list_vst_tng_ASV_relative[["rf_mat"]]
mic_data_tng_Species <- list_vst_tng_Species_relative[["rf_mat"]]
mic_data_tng_Genus <- list_vst_tng_Genus_relative[["rf_mat"]]
#mic_data_tng_Family <- list_vst_tng_Family_relative[["rf_mat"]]

# aggregate sets at different taxonomic levels
#mic_data_stl <- cbind(mic_data_stl_ASV[,"Y_class_all"],
#                      mic_data_stl_ASV[,2:ncol(mic_data_stl_ASV)-1], 
#                      mic_data_stl_Species[,2:ncol(mic_data_stl_Species)-1],
#                      mic_data_stl_Genus[,2:ncol(mic_data_stl_Genus)-1],
#                      mic_data_stl_Family[,2:ncol(mic_data_stl_Family)-1])

mic_data_tng <- cbind(mic_data_tng_ASV[,"Y_class_all"],
                      mic_data_tng_ASV[,2:ncol(mic_data_tng_ASV)-1])
names(mic_data_tng)[1] <- "class"
names(mic_data_tng)[2:length(names(mic_data_tng))]<-paste0("mic_",names(mic_data_tng)[2:length(names(mic_data_tng))])

mic_data_tng_Species <- cbind(mic_data_tng_Species[,"Y_class_all"],
                              mic_data_tng_Species[,2:ncol(mic_data_tng_Species)-1])
names(mic_data_tng_Species)[1] <- "class"
names(mic_data_tng_Species)[2:length(names(mic_data_tng_Species))]<-paste0("mic_",names(mic_data_tng_Species)[2:length(names(mic_data_tng_Species))])

mic_data_tng_Genus <- cbind(mic_data_tng_Genus[,"Y_class_all"],
                            mic_data_tng_Genus[,2:ncol(mic_data_tng_Genus)-1])
names(mic_data_tng_Genus)[1] <- "class"
names(mic_data_tng_Genus)[2:length(names(mic_data_tng_Genus))]<-paste0("mic_",names(mic_data_tng_Genus)[2:length(names(mic_data_tng_Genus))])

# Now combine sample data 
samp_tng <-  data.frame(sample_data(phy_vst_tng))
samp_tng$Severity <- NULL
samp_tng$SampleID <- NULL
rownames(samp_tng)
rownames(mic_data_tng)
names(samp_tng) <- paste0('cc_',names(samp_tng))
names(samp_tng)

str(samp_tng)
t_data_tng <- cbind(samp_tng,mic_data_tng)
# Lets see how many rows gets removed if we remove NAs
t_data_tng_comp <- t_data_tng[complete.cases(t_data_tng),]

# cc and microbiome for tng Species
t_data_tng_Species <- cbind(samp_tng,mic_data_tng_Species)
t_data_tng_Species <- t_data_tng_Species[complete.cases(t_data_tng_Species),]

t_data_tng_Genus <- cbind(samp_tng,mic_data_tng_Genus)
t_data_tng_Genus <- t_data_tng_Genus[complete.cases(t_data_tng_Genus),]

#t_data_tng_Family <- cbind(samp_tng,mic_data_tng_Family)
#t_data_tng_Family <- t_data_tng_Family[complete.cases(t_data_tng_Family),]


### run the models #####
CC_DAT <- t_data_stl_comp[,c(grep("class",names(t_data_stl_comp)),grep('cc_',names(t_data_stl_comp)))]

data_list <- list(CC_DAT, t_data_stl_comp, t_data_stl_Species, t_data_stl_Genus,
                  t_data_tng_comp, t_data_tng_Species, t_data_tng_Genus)
names(data_list) <- c("CC_DAT","CC_STL_ASV","CC_STL_Species","CC_STL_Genus",
                      "CC_TNG_ASV","CC_TNG_Species","CC_TNG_Genus")

#data_list <- list(t_data_stl_comp)
#names(data_list) <- c("CC_STL_ASV")

#data_list <- list(cc_dt_stl)
#names(data_list) <- c("CC_only")
ml_res_holder <-  R6Class("results_container",
                          public = list(
                            f1_list = NULL,
                            cfm_list = NULL,
                            roc_dt_list = NULL,
                            importance_list = NULL,
                            lime_list = NULL,
                            f1_score_list = NULL
                          )

)
ml_res_holder_list <- list()


seed_array <- seq(112,122)
# loop over all the data modalities
start_list <- 1
end_list <- 7
iseed <- 1
is <- 1
for (is in seq(1,length(seed_array))){
  iseed <- seed_array[is]
  set.seed(iseed)
  out <- tryCatch(
    {
      ml_res <- ml_res_holder$new()
      f1_list<-list()
      cfm_list <- list()
      roc_dt_list <- list()
      importance_list <- list()
      lime_list <- list()
      f1_score_list <- list()
      ilist <- 1
      #for (ilist in seq(1,length(data_list))){
      for (ilist in seq(start_list,end_list)){  
        ### load data in #### 
        dat_to_use <- data_list[[ilist]]
        dat_to_use
        dat_to_use$class
        #dat_to_use$class[dat_to_use$class=="Y"] <- "1"
        #dat_to_use$class[dat_to_use$class=="N"] <- "0"
        
        #### first run the boruta using LOO CV #####
        mlo <- ml_data$new()
        mlo$fulldata <- dat_to_use
        mlo$fulldata
        folds<- vfold_cv(mlo$fulldata, v=nrow(mlo$fulldata),repeats=1) 
        vimp_list<-list()
        len<-1
        mlo_boruta_cv_list <- list()
        for(len in 1:nrow(folds)){
          mlo_tmp <- mlo
          mlo_tmp$do_rsample_split_fulldata(v_in = nrow(mlo_tmp$fulldata), repeats_in = 1, fold_id = len)
          mlo_tmp$run_boruta(myseed = iseed, doTrace = T, maxRuns = 100, opt_encoding = F, is_splitted_data = T)
          mlo_tmp$boruta_importance
          mlo_boruta_cv_list[[len]] <- mlo_tmp
          vimp_list[[len]] <- mlo_tmp$boruta_importance
        }
        sel_species <-unique(unlist(vimp_list))
        
        ####### after getting Boruta feature selection run LOO CV for model accuracy estimation  ######
        mlo <- ml_data$new()
        mlo$fulldata <- dat_to_use
        mlo$fulldata
        len <- 1
        f1_score_list_cv <- list()
        lime_list_cv <- list()
        for(len in 1:nrow(folds)){
          split<-folds$splits[[len]]
          mlo_tmp <- mlo
          mlo_tmp$fulldata <- mlo_tmp$fulldata[,c("class",sel_species)]
          mlo_tmp$do_rsample_split_fulldata(v_in = nrow(mlo_tmp$fulldata), repeats_in = 1, fold_id = len)
          mlo_tmp$run_rf_classification(myseed = iseed, ntrees = 1000, ncores = 15, opt_encoding = F, is_splitted_data = T, o_run_pimp = F, o_run_lime = T)
          f1_score_list_cv[[len]] <- mlo_tmp$rfc_f1
          lime_list_cv[[len]] <- mlo_tmp$rfc_lime_explanation
        }  
        
        ##### F1 score and R0C  ######
        f1_dt <-  do.call("rbind",f1_score_list_cv)
        f1_dt$data <- names(data_list)[ilist]
        f1_list[[ilist]] <- f1_dt
        write.csv(f1_dt, paste0(fig_folder,"/f1_",names(data_list)[ilist],".csv"))
        lime_list[[ilist]] <- lime_list_cv
        lime_list[[ilist]]$name <- names(data_list)[ilist]
        
        roc_var_dt <- roc(f1_dt$test,f1_dt$pred,plot=TRUE,smooth = F)
        roc_dt <- data.frame(
          tpp=roc_var_dt$sensitivities*100, ## tpp = true positive percentage
          fpp=(1 - roc_var_dt$specificities)*100, ## fpp = false positive precentage
          #data = f1_dt,
          AUC=as.numeric(gsub("Are.*:","",roc_var_dt$auc)))
        roc_dt$variable <- names(data_list)[ilist]
        roc_dt_list[[ilist]] <- roc_dt
        
        p_roc <- ggplot(roc_dt)+
          geom_path(aes(x= fpp,y = tpp),alpha = 1,size = 1)+
          geom_text(aes(x= fpp,y = tpp,label = AUC))+
          theme_classic()
        p_roc
        
        f1_score <- as.data.frame(F1_Score(y_pred=f1_dt$pred, y_true =f1_dt$test,positive = 1))
        names(f1_score) <-"F1"
        conf_rf <- confusionMatrix(as.factor(f1_dt$pred),as.factor(f1_dt$test), positive = "1")
        conf_rf$data_name <- names(data_list)[ilist]
        f1_score$data_name <- names(data_list)[ilist]
        f1_score_list[[ilist]] <- f1_score
        cfm_list[[ilist]] <- conf_rf
        
        ####### run the full final model   #############
        # Now run random forest model on full data for importance display
        mlo_all_data <- ml_data$new()
        mlo_all_data$fulldata <- dat_to_use
        mlo_all_data$fulldata
        mlo_all_data$fulldata <- mlo_all_data$fulldata[,c("class",sel_species)]
        mlo_all_data$train_data <- mlo_all_data$fulldata
        mlo_all_data$run_rf_classification(myseed = iseed, ntrees = 1000, ncores = 15, opt_encoding = F, is_splitted_data = F, o_run_pimp = F)
        mlo_all_data$rfc_model$importance
        importance_tmp <- as.data.frame(mlo_all_data$rfc_model$variable.importance)
        names(importance_tmp)[1] <- "importance"
        importance_tmp$data_name <- names(data_list)[ilist]
        importance_list[[ilist]] <- importance_tmp[order(importance_tmp$importance,decreasing = T),]
      }
      ml_res$f1_list <- f1_list
      ml_res$cfm_list <- cfm_list
      ml_res$roc_dt_list <- roc_dt_list
      ml_res$importance_list <- importance_list
      ml_res$lime_list <- lime_list
      ml_res$f1_score_list <- f1_score_list
      ml_res
    },
    error = function(cond) {
      print(iseed)
      print
      return(NA)
    }
  )
  ml_res_holder_list[[is]] <- out  
}

if (o_RELAB == F){
  save(phy_stl, phy_tng, phy_vst_stl, phy_vst_tng, data_list, ml_res_holder_list, fig_folder,file = paste(fig_folder,"/","ml_classification_severity_REDEMSVIR_SEVERITY_batch_vst.Rdata", sep=""))
}else{
  save(phy_stl, phy_tng, phy_vst_stl, phy_vst_tng, data_list, ml_res_holder_list, fig_folder,file = paste(fig_folder,"/","ml_classification_severity_REDEMSVIR_SEVERITY_batch_relab.Rdata", sep=""))
}




