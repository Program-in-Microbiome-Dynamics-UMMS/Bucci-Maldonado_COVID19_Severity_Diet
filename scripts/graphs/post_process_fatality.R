rm(list=ls())
library(tune)
library(workflows)
library(tidymodels)
library(tictoc)
library(DESeq2)
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(gtools)
library(tidyr)
library(ggprism)
library(ranger)
library(caret)
library(pROC)
library(plyr)

##### load data ##### 
#fileName_1 <- "../../results/RF_final_fatality_aggregate_taxonomy_object_oriented_relab_F/2022-07-14/ml_classification_severity_REDEMSVIR_SEVERITY_batch_relab.Rdata"
#fileName_1 <- "../../results/RF_final_fatality_aggregate_taxonomy_object_oriented_relab_F/2022-07-14/ml_classification_severity_REDEMSVIR_SEVERITY_batch_vst.Rdata"
#fileName_1 <- "../../results/RF_final_fatality_aggregate_taxonomy_object_oriented_relab_F_with_severity_in/2022-07-15/ml_classification_severity_REDEMSVIR_SEVERITY_batch_vst.Rdata"
fileName_1 <- "../../results/FATALITY-RF_final_fatality_aggregate_taxonomy_object_oriented_relab_F_with_severity_in/2022-07-15/ml_classification_severity_REDEMSVIR_SEVERITY_batch_vst.Rdata"

load(fileName_1)

#analysis_name <- "Severity-4L"
analysis_name <- "Fatality"

phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Synergistetes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Euryarchaeota" = "#8c8c8c")

##### F1 score plot #####
i<-1
res_f1_list <- list()
for (i in seq(1,length(ml_res_holder_list))){
  condition <- length(ml_res_holder_list[[i]])
  if (condition>1){
    tmp <- ml_res_holder_list[[i]]
    tmp_f1<-tmp$f1_list  
    tmp_f1 <- do.call('rbind',tmp_f1)
    tmp_f1$nseed <- i
    tmp_f1$data_mod <- paste0(tmp_f1$data, "-", tmp_f1$nseed)
    res_f1_list[[i]] <- tmp_f1
  }
}
res_f1_list
pred_dt <- do.call("rbind",res_f1_list)
data_var <- unique(pred_dt$data_mod)
roc_var_list  <- list()
var <- data_var[1]
for(var in data_var){
  pred_dt_sel <-  pred_dt[pred_dt$data_mod ==var,]
  roc_var_dt <- roc(pred_dt_sel$test,pred_dt_sel$pred,plot=TRUE,smooth = F)
  roc_dt <- data.frame(
    tpp=roc_var_dt$sensitivities*100, ## tpp = true positive percentage
    fpp=(1 - roc_var_dt$specificities)*100, ## fpp = false positive precentage
    data = var,
    data_0 = unique(pred_dt_sel$data),
    AUC=as.numeric(gsub("Are.*:","",roc_var_dt$auc)))
  roc_var_list[[var]] <- roc_dt
}
roc_comb <-  do.call("rbind",roc_var_list)
roc_comb$AUC <- as.character(signif(100*roc_comb$AUC, digits = 4))
roc_comb_sub <- roc_comb[roc_comb$tpp==100,]
roc_comb_sub$AUC <- as.numeric(roc_comb_sub$AUC)
p_roc <- ggplot(roc_comb_sub)+
  geom_boxplot(aes(x= data_0, y = AUC, fill = data_0, group = data_0),
              alpha = 0.5, size = 1)+
  geom_jitter(aes(x= data_0, y = AUC, fill = data_0, group = data_0),
              alpha = 0.8, size = 3, shape = 21)+
  stat_boxplot(aes(x= data_0, y = AUC, fill = data_0, group = data_0), geom='errorbar',)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  xlab("")+
  ylab("F1 score")+
  theme_prism()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle(analysis_name)
  
print(p_roc)
write.csv(roc_comb_sub,paste0(fig_folder, "/F1_score_",analysis_name,".csv"))
ggsave(paste0(fig_folder, "/F1_score_",analysis_name,".pdf"), p_roc,width=7, height =5)

roc_comb_sub_2 <- roc_comb_sub[grepl('TNG',roc_comb_sub$data)==F,] 

p_roc_2 <- ggplot(roc_comb_sub_2)+
  geom_boxplot(aes(x= data_0, y = AUC, fill = data_0, group = data_0),
               alpha = 0.5, size = 1)+
  geom_jitter(aes(x= data_0, y = AUC, fill = data_0, group = data_0),
              alpha = 0.8, size = 4, shape = 21)+
  stat_boxplot(aes(x= data_0, y = AUC, fill = data_0, group = data_0), geom='errorbar',)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  xlab("")+
  ylab("F1 score")+
  theme_prism()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(y = "prism_offset_minor")+
  scale_y_continuous(limits = c(50, 90))+
  ggtitle(analysis_name)

print(p_roc_2)

ggsave(paste0('./', "/F1_score_STL_only_",analysis_name,".pdf"), p_roc_2,width=5, height =5)


roc_comb_sub_3 <- roc_comb_sub[grepl('STL',roc_comb_sub$data)==F,] 

p_roc_3 <- ggplot(roc_comb_sub_3)+
  geom_boxplot(aes(x= data_0, y = AUC, fill = data_0, group = data_0),
               alpha = 0.5, size = 1)+
  geom_jitter(aes(x= data_0, y = AUC, fill = data_0, group = data_0),
              alpha = 0.8, size = 4, shape = 21)+
  stat_boxplot(aes(x= data_0, y = AUC, fill = data_0, group = data_0), geom='errorbar',)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  xlab("")+
  ylab("F1 score")+
  theme_prism()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  guides(y = "prism_offset_minor")+
  scale_y_continuous(limits = c(50, 90))+
  ggtitle(analysis_name)

print(p_roc_3)

ggsave(paste0('./', "/F1_score_TNG_only_",analysis_name,".pdf"), p_roc_3,width=5, height =5)



mod1<-lm(AUC~data_0 , data = roc_comb_sub)
summary(mod1)
anova_mod <- aov(mod1)
TukeyHSD(anova_mod)


roc_comb_sub_summary <- ddply(roc_comb_sub, ~data_0, summarise,
                              mean_AUC=mean(AUC),
                              n = length(AUC),
                              #se = sd(AUC)/n
                              se = sd(AUC))
p_roc_error <- ggplot(roc_comb_sub_summary)+
  geom_errorbar(aes(x= data_0, y = mean_AUC, ymax = mean_AUC+se, ymin = mean_AUC-se,
                    color = data_0, group = data_0),
               alpha = 1, size = 1, width = 0.2)+
  geom_point(aes(x= data_0, y = mean_AUC,
                    fill = data_0, group = data_0),
                alpha = 1, size = 4, width = 0.2, shape = 21)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  xlab("")+
  ylab("F1 score")+
  theme_prism()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle(analysis_name)
print(p_roc_error)
write.csv(roc_comb_sub_summary,paste0(fig_folder, "/F1_score_errorbar_",analysis_name,".csv"))
ggsave(paste0(fig_folder, "/F1_score_errorbar_",analysis_name,".pdf"), p_roc_error, width=7, height =5)

##### Importance analysis #####
i<-1
res_importance_list <- list()
for (i in seq(1,length(ml_res_holder_list))){
  condition <- length(ml_res_holder_list[[i]])
  if (condition>1){
    tmp <- ml_res_holder_list[[i]]
    tmp_importance<-tmp$importance_list 
    j<-1
    for (j in seq(1,length(tmp_importance))){
      tmp_importance[[j]]$var_name <- rownames(tmp_importance[[j]])
    }
    tmp_importance <- do.call('rbind',tmp_importance)
    tmp_importance$nseed <- i
    res_importance_list[[i]] <- tmp_importance
  }
}
res_importance_list

res_importance <- do.call('rbind',res_importance_list)
res_importance$var_name[grep('Race',res_importance$var_name)] <- "Race/Ethnicity"
res_importance$var_name[grep('vanco',res_importance$var_name)] <- "vancomycin"
res_importance$var_name[grep('zosyn',res_importance$var_name)] <- "zosyn"
res_importance$var_name[grep('meropen',res_importance$var_name)] <- "meropenem"
res_importance$var_name[grep('cefazolin',res_importance$var_name)] <- "cefazolin"
res_importance$var_name[grep('Antibiotics',res_importance$var_name)] <- "antibiotics"
res_importance$var_name[grep('heart_',res_importance$var_name)] <- "heart failure"
res_importance$var_name[grep('cc_',res_importance$var_name)] <- gsub('cc_','',res_importance$var_name[grep('cc_',res_importance$var_name)])
     
res_importance_summary <- ddply(res_importance, ~data_name+var_name, summarise,
                                mean_importance = mean(importance),
                                median_importance = median(importance),
                                n = length(importance),
                                se = sd(importance),
                                mad_importance = mad(importance))


reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}
res_importance_summary <- res_importance_summary[res_importance_summary$mean_importance > 0,]
res_importance_summary$var_name_mod <- res_importance_summary$var_name
res_importance_summary$var_name_mod <- gsub("\\_", " ", res_importance_summary$var_name_mod)
res_importance_summary$var_name_mod <- gsub("X\\.", "", res_importance_summary$var_name_mod)
res_importance_summary$var_name_mod <- gsub("\\.", " ", res_importance_summary$var_name_mod)
res_importance_summary$var_name_mod <- gsub("  ", " ", res_importance_summary$var_name_mod)
res_importance_summary$data_name_mod <- res_importance_summary$data_name
res_importance_summary$data_name_mod <- gsub("\\_", " ", res_importance_summary$data_name_mod)
res_importance_summary <- res_importance_summary[res_importance_summary$var_name_mod!="xxxx",]

g_importance <- ggplot(res_importance_summary, aes(reorder_within(var_name_mod, median_importance, data_name_mod), median_importance)) +
  geom_errorbar(size=1, alpha =1, aes(color = data_name, ymin=median_importance-mad_importance, ymax=median_importance+mad_importance), width = 0.3) +
  geom_point(size = 4, shape = 21, aes(fill=data_name)) +
  scale_x_reordered() +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(data_name_mod ~ ., scales = "free") +
  coord_flip() +
  theme_prism() +
  theme(panel.grid.major.y = element_blank())+
  ylab("Median Importance")+
  xlab('Predictor')
  #guides(y = "prism_offset_minor")+
  ggtitle(analysis_name)
print(g_importance)
write.csv(res_importance_summary,paste0(fig_folder, "/importance",analysis_name,".csv"))
ggsave(paste0('./', "/Importance_",analysis_name,".pdf"), g_importance, width=20, height =22)


res_importance_summary_sub <- res_importance_summary[grepl('TNG',res_importance_summary$data_name_mod)==F,]

g_importance_sub <- ggplot(res_importance_summary_sub, aes(reorder_within(var_name_mod, median_importance, data_name_mod), median_importance)) +
  geom_errorbar(size=1, alpha =1, aes(color = data_name, ymin=median_importance-mad_importance, ymax=median_importance+mad_importance), width = 0.3) +
  geom_point(size = 4, shape = 21, aes(fill=data_name)) +
  scale_x_reordered() +
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(data_name_mod ~ ., scales = "free", nrow = 1) +
  coord_flip() +
  theme_prism() +
  theme(panel.grid.major.y = element_blank())+
  ylab("Median Importance")+
  xlab('Predictor')+
#guides(y = "prism_offset_minor")+
ggtitle(analysis_name)
print(g_importance_sub)
#write.csv(res_importance_summary,paste0(fig_folder, "/importance",analysis_name,".csv"))
ggsave(paste0('./', "/Importance_",analysis_name,".pdf"), g_importance_sub, width=25, height =10)


##### boxplot abundances of important features #####
tchoose <- "CC_STL_Species" 
res_importance_summary_temp <- res_importance_summary[res_importance_summary$data_name==tchoose,]
features_to_plot <- res_importance_summary_temp$var_name
features_to_plot
features_to_plot_micro <- features_to_plot[grep("mic_", features_to_plot)]
xdata <- data_list[tchoose]
xdata<- xdata$CC_STL_Species[,c("class", features_to_plot_micro)]
x_data_m <- melt(xdata)

x_data_m$variable_mod <- x_data_m$variable
x_data_m$variable_mod <- gsub("\\_", " ", x_data_m$variable_mod)
x_data_m$variable_mod <- gsub("X\\.", "", x_data_m$variable_mod)
x_data_m$variable_mod <- gsub("\\.", " ", x_data_m$variable_mod)
x_data_m$variable_mod <- gsub("  ", " ", x_data_m$variable_mod)
x_data_m$variable_mod <- x_data_m$variable_mod
x_data_m$variable_mod <- gsub("\\_", " ", x_data_m$variable_mod)
x_data_m$variable_mod <- gsub("mic ", "", x_data_m$variable_mod)
x_data_m <- x_data_m[x_data_m$variable_mod!="xxxx",]

violing_microbiome <- ggplot(x_data_m, aes(class, value, fill = class))+
  geom_jitter(shape =21, size =4)+
  geom_violin(alpha=0.5, na.rm = T, outlier.colour = NA)+
  facet_wrap(~variable_mod, scales = "free")+
  ylab("Relative Abundance")+
  xlab("")+
  scale_fill_manual(values=c("dodgerblue4","firebrick1"))+
  theme_prism()
ggsave(paste0(fig_folder, "/Abundance_microbiome_","CC_STL_Species",analysis_name,".pdf"), violing_microbiome, width=22, height =12)


# SIRUS analysis
library(sirus)
tchoose <- "CC_STL_Species"
res_importance <- do.call('rbind',res_importance_list)
res_importance_sub <- res_importance[res_importance$data_name==tchoose,]
boruta_sp <- unique(res_importance_sub$var_name)
t_data <- data_list$CC_STL_Species
t_data_sub <-  t_data[,c(boruta_sp,"class")]
train_X <-   t_data_sub[,!colnames(t_data_sub) %in% "class"]
y <-  ifelse(t_data_sub$class == "1",1,0)
cv.grid <- sirus.cv(
  data = train_X,
  y,
  type = "classif",
  nfold = 3,
  ncv = 2,
  num.rule.max = 25,
  q = 10,
  discrete.limit = 10,
  num.trees.step = 1000,
  alpha = 0.05,
  num.trees = NULL,
  num.threads = NULL,
  replace = TRUE,
  sample.fraction = NULL,
  verbose = TRUE,
  seed = 100
)

plot.error <- sirus.plot.cv(cv.grid)$error
plot(plot.error)

## fit SIRUS
sirus.m <- sirus.fit(data = train_X , y, p0 = cv.grid$p0.pred,seed = 1000)
sirus.m$num.trees
sel_rules <-  sirus.print(sirus.m, digits = 3)[-1]

length(sel_rules)

library(stringr)
str_dt  <-  str_split_fixed(sel_rules,pattern = "then", n = 2)
str_dt

# Unique taxa in all rules
sirus_taxa <-  gsub("<.*","",str_dt[,1])
sirus_taxa <-  trimws(gsub("if ","",sirus_taxa))
sirus_taxa <- gsub(' in \\{1\\}','',sirus_taxa)
sirus_taxa <- gsub(' in \\{0\\}','',sirus_taxa)
sirus_taxa

str_dt_2 <-  str_split_fixed(str_dt[,2],pattern = "else", n = 2)
str_dt_2[,1] <- 100* as.numeric(trimws(gsub(" \\(.*","",str_dt_2[,1]) ))
str_dt_2[,2] <- 100* as.numeric(trimws(gsub(" \\(.*","",str_dt_2[,2]) ))

str_dt <-  data.frame(str_dt[,1],str_dt_2)
names(str_dt) <-  c("rule","pr","pnr")
str_dt$rule <- trimws(gsub("if ","",str_dt$rule))

# Conscise rule for each taxa
all_rule_dt <-  str_dt
all_rule_dt_microbiome <- all_rule_dt[grep('mic',all_rule_dt$rule),]
all_rule_dt_microbiome
all_rule_dt_microbiome$Taxa <- gsub("<.*","",all_rule_dt_microbiome$rule) 
all_rule_dt_microbiome$Taxa <- gsub(' ','',all_rule_dt_microbiome$Taxa)

ir <- 1
pdf(paste0(fig_folder, "/sirus_results","CC_STL_Species",".pdf"), width=8, height =8)
for (ir in seq(1,nrow(all_rule_dt_microbiome))){
  isalone <- grepl('&',all_rule_dt_microbiome$Taxa[ir])
  if (isalone == F){
    rule_label <-  paste0("if ",
                          all_rule_dt_microbiome$rule[ir],
                          " then  Ps  = ",
                          as.numeric(all_rule_dt_microbiome$pr[ir]),"%", " else Pns = ",
                          as.numeric(all_rule_dt_microbiome$pnr[ir]),"%"  )
    xdata <- data_list["CC_STL_Genus"]
    xdata<- xdata$CC_STL_Genus[,c("class", all_rule_dt_microbiome$Taxa[ir])]
    x_data_m <- melt(xdata)
    x_data_m$rule <- all_rule_dt_microbiome$rule[ir]
    x_data_m$pr <- all_rule_dt_microbiome$pr[ir]
    x_data_m$pnr <- all_rule_dt_microbiome$pnr[ir]
    x_data_m$thres <-  as.numeric(gsub(".*< ","",x_data_m$rule))
    
    violing_microbiome <- ggplot(x_data_m, aes(class, value, fill = class))+
      geom_jitter(shape =21, size =4)+
      geom_violin(alpha=0.5, na.rm = T, outlier.colour = NA)+
      geom_hline(data = x_data_m, 
                 aes(yintercept= thres),color = "black", linetype = "dashed",size = 1)+
      ylab("Abundance")+
      xlab("")+
      scale_fill_manual(values=c("dodgerblue4","firebrick1"))+
      ggtitle(rule_label)+
      theme_prism()+
      theme(plot.title = element_text(size = 10, face = "bold"))
    print(violing_microbiome)
  }
}
dev.off()

#  heatmap of SIRUS identified bacteria
pred_var <- "Fatality_new"
sig_microbes <- all_rule_dt_microbiome$Taxa
sig_microbes <- gsub('mic_','',sig_microbes)
phy_deseq <- phy_stl
dds <- phyloseq_to_deseq2(phy_deseq, ~ 1) 
dds <- estimateSizeFactors(dds,"poscounts")
dds <- estimateDispersions(dds)
dds <- DESeq(dds,fitType= "local")
vst_dt <- getVarianceStabilizedData(dds)
dt_meta <-  data.frame(sample_data(phy_deseq))
met  <- dt_meta
mat  <- vst_dt 
vst_dt[vst_dt<0] <- 0
met_samp <- met
met_samp$sample  <- rownames(met_samp)
rownames(met_samp) <- NULL
mat <-  mat[rownames(mat) %in% sig_microbes,]
order_res <- tax_table(phy_deseq)[taxa_names(phy_deseq) %in% sig_microbes,]

split_cols<-  met[,pred_var]
split_cols <- factor(split_cols, 
                     levels= c("0","1"))
levels(split_cols) <- c("0","1")

status_col <-  c("blue","red")
names(status_col) <- levels(split_cols)

library(gtools)
library(ComplexHeatmap)
ha_column = HeatmapAnnotation(Status =  met[,pred_var],
                              col=list(Status = status_col))

library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))
teddy_cols <- c('white','#fecc5c','#fd8d3c','#f03b20','#bd0026')
jet.colors <-c("white", "blue", "#007FFF", "cyan",
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
ht_cols <- jet.colors
mat2 <-  as.data.frame(mat) 

library(circlize)
split_rows <-  order_res[,"Order"]
split_rows <- factor(split_rows, levels= unique(split_rows))
dt_col_phylum <- data.frame(phylum_col)
dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
comb_dt <- as.data.frame(order_res)
shades_num <- comb_dt %>%
  select(Phylum,Order) %>%
  group_by(Phylum) %>%
  mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique %>% as.data.frame()

shades_num <- merge(shades_num, dt_col_phylum, by = "Phylum")
shades_num$phylum_col <-  as.character(shades_num$phylum_col)

mat3 <- mat2
rownames(mat3) <- gsub("\\_"," ", rownames(mat3))
rownames(mat3) <- gsub("X\\."," ", rownames(mat3))
rownames(mat3) <- gsub("\\.\\."," ", rownames(mat3))
rownames(mat3) <- gsub("\\."," ", rownames(mat3))

ht1 = Heatmap(as.matrix(mat3), name = "vst", column_title = NA, 
              top_annotation = ha_column,
              col = ht_cols,
              cluster_column_slices = F,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 10,fontface = "bold"), 
              row_title_rot = 0,
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = T, 
              row_names_gp = gpar(fontsize = 9),
              row_names_max_width = max_text_width(rownames(mat2), gp = gpar(fontsize = 12)),
              show_column_names = F, show_row_names = T,
              column_names_side = "bottom", na_col="white",
              show_row_dend = F,
              border = T)

pdf(paste0(fig_folder,"/heatmap_severity_",names(data_list)[i],".pdf"), width = 10, height = 6)
draw(ht1)
dev.off()
    
# print SIRUS rules to text 
library("stargazer")
stargazer(all_rule_dt,            
          summary = FALSE,
          type = "text",
          out = "data_stargazer_txt.txt")
write.csv(all_rule_dt,"./sirus_fatality.csv")

# x_data_m$variable_mod <- x_data_m$variable
# x_data_m$variable_mod <- gsub("\\_", " ", x_data_m$variable_mod)
# x_data_m$variable_mod <- gsub("X\\.", "", x_data_m$variable_mod)
# #res_importance_summary$var_name_mod <- gsub("\\..", " ", res_importance_summary$var_name_mod)
# x_data_m$variable_mod <- gsub("\\.", " ", x_data_m$variable_mod)
# x_data_m$variable_mod <- gsub("  ", " ", x_data_m$variable_mod)
# x_data_m$variable_mod <- x_data_m$variable_mod
# x_data_m$variable_mod <- gsub("\\_", " ", x_data_m$variable_mod)
# x_data_m$variable_mod <- gsub("mic ", "", x_data_m$variable_mod)
# x_data_m <- x_data_m[x_data_m$variable_mod!="xxxx",]
# 
# violing_microbiome <- ggplot(x_data_m, aes(class, value, fill = class))+
#   geom_jitter(shape =21, size =4)+
#   geom_violin(alpha=0.5, na.rm = T, outlier.colour = NA)+
#   facet_wrap(~variable_mod, scales = "free")+
#   ylab("Relative Abundance")+
#   xlab("")+
#   scale_fill_manual(values=c("dodgerblue4","firebrick1"))+
#   theme_prism()
# 
# 
# # Plot lime
# imod <- 2
# lime_df_all <- do.call("rbind",lime_list[[imod]])
# lime_dat <-  lime_df_all
# #lime_dat$Data <-  names(data_list)[i]
# #lime_model_list[[names(data_list)[i]]] <- lime_dat
#   
# # True labels data
# t_data <- data_list[[imod]]
# t_data_r2 <- t_data
# true_labels_dt <-  data.frame(case = rownames(t_data_r2),true_label = t_data_r2$Y_class_all)
#   
# lime_dt <- lime_df_all
# explanation <- lime_dt
# explanation$data <-  NULL
# explanation$prediction <- NULL
# explanation <- unique(data.frame(explanation))
# explanation <- explanation[,c("case","label" ,"feature_weight","feature_desc" )]
#   
# str(explanation)
# explanation <- merge(explanation, true_labels_dt, by = "case")
# names(explanation)
# explanation$case_true <- paste0(explanation$case," ","("," ",explanation$true_label," ",")")
#   
#   #explanation <- explanation[grep("Phasco",explanation$feature_desc),]
# explanation$feature_desc <- gsub("\\_"," ", explanation$feature_desc)
# explanation$feature_desc <- gsub("X\\."," ", explanation$feature_desc)
# explanation$feature_desc <- gsub("\\.\\."," ", explanation$feature_desc)
# p <- ggplot(explanation, aes(case_true, feature_desc)) + geom_tile(aes(fill = feature_weight)) + 
#     # scale_x_discrete("Case", expand = c(0, 0)) + scale_y_discrete("Feature", expand = c(0, 0)) +
#     scale_fill_gradient2("Feature\nweight", low = "firebrick", mid = "#f7f7f7", high = "steelblue") + 
#     theme_prism() +
#     theme(panel.border = element_rect(fill = NA,colour = "grey60", size = 1), panel.grid = element_blank(), 
#           legend.position = "right", axis.text.x = element_text(angle = 45, 
#                                                                 hjust = 1, vjust = 1)) +
#     theme(axis.text.y = element_text(color = "grey20", size = 15, face = "bold"))+
#     facet_wrap(~label,scales = "free_x")+
#     ylab("Features")+
#     xlab("Truth")
# pdf(paste0(fig_folder,"/LIME_severity_",names(data_list)[i],".pdf"), width = 18, height = 10)
# print(p)
# dev.off()
# write.csv(explanation, paste0(fig_folder,"/LIME_severity_",names(data_list)[i],".csv"))
#   
#   explanation
#   u_rule <- unique(explanation$feature_desc)
#   u_rule
#   lime_freq_list <- list()
#   irule <- 1
#   for (irule in seq(1,length(u_rule))){
#     explanation_sub <- explanation[which(explanation$feature_desc==u_rule[irule]),]  
#     explanation_sub 
#     correct_vals <- explanation_sub$label == explanation_sub$true_label
#     correct_vals 
#     explanation_sub_true <- explanation_sub[which(correct_vals==T),]
#     true_contribution <- length(which(explanation_sub_true$feature_weight !=0))/length(unique(explanation_sub$case)) 
#     total_contribution <- length(which(explanation_sub$feature_weight !=0))/length(unique(explanation$case)) 
#     explanation_sub_true_yes <- explanation_sub_true[explanation_sub_true$label=="Y",]
#     explanation_sub_true_no <- explanation_sub_true[explanation_sub_true$label=="N",]
#     if (nrow(explanation_sub_true_yes)>0){
#       sign_association <- sign(explanation_sub_true_yes$feature_weight[1])
#     }else if(nrow(explanation_sub_true_no)>0){
#       sign_association <- -1 * sign(explanation_sub_true_no$feature_weight[1])
#       }
#     lime_freq_list[[irule]] <-c(u_rule[irule],true_contribution,total_contribution,sign_association)
#   }
#   lime_freq_list
#   lime_freq_list <- do.call('rbind',lime_freq_list)
#   write.csv(lime_freq_list, paste0(fig_folder,"/LIME_severity_frequency_list_",names(data_list)[i],".csv"))
#   
#   # Now run random forest model on full data for importance display
#   library(ranger)
#   library(caret)
#   t_data_r2$Y_class_all <- as.factor(t_data_r2$Y_class_all)
#   split <- t_data_r2
#   rf_model <- ranger(Y_class_all~.,data=t_data_r2,importance = "permutation",
#                       probability = T,
#                       num.trees = 10000,
#                       replace=FALSE)
#   #rf_model_full_data<-run_rf_with_Bayesian_Optimization(split, is_cv_run = F) # USING FULL MODEL
#   #rf_model <- rf_model_full_data$rf_model
#   # Add RF list for each model
#   rf_final_list[[names(data_list)[i]]] <- rf_model
#   pred_ranger <- rf_model$predictions
#   res_pred = ifelse(pred_ranger[,"Y"] >= 0.5,"Y","N")
#   t_data_all$Y_class_all <- as.factor(t_data_all$Y_class_all)
#   conf_rf <- confusionMatrix(factor(res_pred),t_data_all$Y_class_all,positive = "Y")
#   
#   # Takes a list of variable importance
#   # Prints out the variable importance 
#   pimp_plot<-plot_importance(rf_model)
#   print(pimp_plot)
#   library(ggthemes)
#   ggsave(paste0(fig_folder, "/importance_Severity_",names(data_list)[i],".pdf"), pimp_plot,width=10, height =8)
#   
#   # if microbioms is included make heatmap and plot of significant bacteria
#   if (grepl('STL',names(data_list)[i]) || grepl('TNG',names(data_list)[i])){
#     mod <- rf_model
#     varimp <- mod$variable.importance
#     dat <- data.frame(variable=names(varimp),importance=varimp)
#     dat <- dat[dat$importance>0,]
#     sig_microbes <- dat$variable[grep('SV',dat$variable)]
#     sig_microbes
#     if (length(sig_microbes)>0){
#       # Heatmap for the genes:
#       # VST phyloseq for Heatmap later
#       #phy_vst <- phy_gene_sel
#       if (length(grep('STL',names(data_list)[i])>0)){
#         phy_deseq <- phy_stl
#       }else{
#         phy_deseq <- phy_tng
#       }
#       phy_deseq <- subset_samples(physeq = phy_deseq, Severity!="none")
#       dds <- phyloseq_to_deseq2(phy_deseq, ~ 1) #replace this with any sample variable(s)
#       dds <- estimateSizeFactors(dds,"poscounts")
#       dds <- estimateDispersions(dds)
#       dds <- DESeq(dds,fitType= "local")
#       vst_dt <- getVarianceStabilizedData(dds)
#       dt_meta <-  data.frame(sample_data(phy_deseq))
#       met  <- dt_meta
#       mat  <- vst_dt 
#       vst_dt[vst_dt<0] <- 0
#       met_samp <- met
#       met_samp$sample  <- rownames(met_samp)
#       rownames(met_samp) <- NULL
#       
#       mat <-  mat[rownames(mat) %in% sig_microbes,]
#       order_res <- tax_table(phy_deseq)[taxa_names(phy_deseq) %in% sig_microbes,]
#       
#       split_cols<-  met$Severity
#       split_cols <- factor(split_cols, 
#                            levels= c("moderate","severe"))
#       levels(split_cols) <- c("moderate","severe")
#       
#       status_col <-  c("blue","red")
#       names(status_col) <- levels(split_cols)
#       
#       library(gtools)
#       library(ComplexHeatmap)
#       ha_column = HeatmapAnnotation(Status =  met$Severity,
#                                     col=list(Status = status_col))
#       
#       library(RColorBrewer)
#       colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
#                                                    "RdYlBu")))
#       teddy_cols <- c('white','#fecc5c','#fd8d3c','#f03b20','#bd0026')
#       jet.colors <-c("white", "blue", "#007FFF", "cyan",
#                      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
#       ht_cols <- jet.colors
#       mat2 <-  as.data.frame(mat) 
#       
#       library(circlize)
#       split_rows <-  order_res[,"Order"]
#       split_rows <- factor(split_rows, levels= unique(split_rows))
#       dt_col_phylum <- data.frame(phylum_col)
#       dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
#       comb_dt <- as.data.frame(order_res)
#       shades_num <- comb_dt %>%
#         select(Phylum,Order) %>%
#         group_by(Phylum) %>%
#         mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique %>% as.data.frame()
#       
#       shades_num <- merge(shades_num, dt_col_phylum, by = "Phylum")
#       shades_num$phylum_col <-  as.character(shades_num$phylum_col)
#       mat3 <- mat2
#       rownames(mat3) <- gsub("\\_"," ", rownames(mat3))
#       rownames(mat3) <- gsub("X\\."," ", rownames(mat3))
#       rownames(mat3) <- gsub("\\.\\."," ", rownames(mat3))
#       rownames(mat3) <- gsub("\\."," ", rownames(mat3))
#       
#       ht1 = Heatmap(as.matrix(mat3), name = "vst", column_title = NA, 
#                     top_annotation = ha_column,
#                     col = ht_cols,
#                     #left_annotation = ha1,
#                     cluster_column_slices = F,
#                     row_names_side = "left",
#                     row_title_gp = gpar(fontsize = 10,fontface = "bold"), 
#                     row_title_rot = 0,
#                     column_split = split_cols,
#                     show_parent_dend_line = F,
#                     width=2, cluster_columns = T, 
#                     row_names_gp = gpar(fontsize = 9),
#                     row_names_max_width = max_text_width(rownames(mat2), gp = gpar(fontsize = 12)),
#                     show_column_names = F, show_row_names = T,
#                     column_names_side = "bottom", na_col="white",
#                     show_row_dend = F,
#                     border = T)
#       #pdf("GSVA_raw.pdf",width = 14,height = 8)
#       pdf(paste0(fig_folder,"/heatmap_severity_",names(data_list)[i],".pdf"), width = 10, height = 6)
#       draw(ht1)
#       dev.off()
#       
#       phy_deseq_plot <- phy_deseq
#       otu_table(phy_deseq_plot) <- otu_table(vst_dt, taxa_are_rows = T)
#       phy_deseq_plot <- prune_taxa(sig_microbes,phy_deseq_plot)
#       ps_sig_m <- psmelt(phy_deseq_plot)
#       
#       ps_sig_m$OTU2 <- ps_sig_m$OTU
#       ps_sig_m$OTU2 <- gsub("\\_"," ", ps_sig_m$OTU2)
#       ps_sig_m$OTU2 <- gsub("X\\."," ", ps_sig_m$OTU2)
#       ps_sig_m$OTU2 <- gsub("\\.\\."," ", ps_sig_m$OTU2)
#       ps_sig_m$OTU2 <- gsub("\\."," ", ps_sig_m$OTU2)
#       
#       gsig <- ggplot()+
#         geom_jitter(data=ps_sig_m,aes(x=get(pred_var),y=Abundance, 
#                                       group=get(pred_var), 
#                                       fill=get(pred_var)), size =4 , alpha=0.5, shape = 21)+
#         geom_violin(data=ps_sig_m,aes(x=get(pred_var),y=Abundance, 
#                                       group=get(pred_var), 
#                                       fill=get(pred_var)),alpha=0.5, outlier.shape = NA)+
#         scale_color_manual(name=pred_var,values=c("blue","red"))+
#         scale_fill_manual(name=pred_var,values=c("blue","red"))+
#         facet_wrap(~OTU2,scales = "free")+
#         ylab("Abundance")+
#         xlab(pred_var)+
#         theme_prism()+
#         theme(axis.text.x = element_blank()) 
#       pdf(paste0(fig_folder,"/boxplot_exp_severity_",names(data_list)[i],".pdf"), width = 25, height = 12)
#       print(gsig)
#       dev.off()
#       write.csv(ps_sig_m[,c("Severity","OTU","Abundance")], paste0(fig_folder,"/boxplot_severity_",names(data_list)[i],".csv"))
#       
#       library(plyr)
#       ps_sig_m_summarised <- ddply(ps_sig_m,~OTU+get(pred_var),
#                                    summarise, mean_Abundance=mean(Abundance))
#       names(ps_sig_m_summarised)[2]<-pred_var
#       
#       dat_csv<-dat
#       dat_csv$direction <- "severe"
#       irow <- 1
#       dat_csv <- dat_csv[grep("SV",dat_csv$variable),]
#       for (irow in seq(1,nrow(dat_csv))){
#         #print(irow)
#         idx <- which(ps_sig_m_summarised$OTU==dat_csv$variable[irow])
#         tmp <- ps_sig_m_summarised[idx,]
#         difference <- tmp$mean_Abundance[tmp$Severity=="severe"] - tmp$mean_Abundance[tmp$Severity=="moderate"]
#         #print(difference)
#         if (difference <0){
#           dat_csv$direction[irow] <- "moderate"
#         }
#       }
#       write.csv(dat_csv, paste0(fig_folder,"/importance_severity_",names(data_list)[i],".csv"))
#       write.csv(order_res, paste0(fig_folder,"/taxonomy_significant_",names(data_list)[i],".csv"))
#     }
#   }
#   
#   if (grepl('CC',names(data_list)[i])){
#     mod <- rf_model
#     varimp <- mod$variable.importance
#     dat <- dat[dat$importance>0,]
#     dat <- data.frame(variable=names(varimp),importance=varimp)
#     sig_CC <- dat$variable[which(grepl('SV_',dat$variable)==F)]
#     sig_CC
#     if (length(sig_CC)>0){
#       cc_dat <- t_data_all[,c(sig_CC,"Y_class_all")]
#       cc_dat$Y_class_all <- as.character(cc_dat$Y_class_all)
#       cc_dat$Y_class_all[cc_dat$Y_class_all=="N"] <- "moderate"
#       cc_dat$Y_class_all[cc_dat$Y_class_all=="Y"] <- "severe"
#       names(cc_dat)[which(names(cc_dat)=="Y_class_all")] <- "Severity"
#       is.fact <- sapply(cc_dat, is.factor)
#       if (length(is.fact)>1){
#         cc_dat_fact <- cc_dat[,c(names(cc_dat)[which(is.fact==T)],"Severity")]
#         cc_dat_fact_m <- melt(cc_dat_fact, id.vars = "Severity")
#         cc_dat_fact_m$variable2 <- cc_dat_fact_m$variable
#         cc_dat_fact_m$variable2 <- gsub("\\_"," ", cc_dat_fact_m$variable2)
#         cc_dat_fact_m$variable2 <- gsub("X\\."," ", cc_dat_fact_m$variable2)
#         cc_dat_fact_m$variable2 <- gsub("\\.\\."," ", cc_dat_fact_m$variable2)
#         cc_dat_fact_m$variable2 <- gsub("\\."," ", cc_dat_fact_m$variable2)
#         cc_dat_fact_m_v2 <- cc_dat_fact_m
#         # array_to_separate <- c('Race','Obesity')
#         # isep <- 1
#         # sep_tk <- c()
#         # for (isep in seq(1,length(array_to_separate))){
#         #   idsep <- grep(array_to_separate[isep],cc_dat_fact_m$variable)
#         #   if (length(idsep)>0){ 
#         #     cc_dat_fact_m_v2 <- cc_dat_fact_m_v2[-idsep,]
#         #     sep_tk <- c(sep_tk,idsep)
#         #     }
#         # }
#         #grep('Obesity',cc_dat_fact_m$variable),]
#         # remove not binary (Race and Obesity and plot separately)
#         g_stack_cc_factor <- ggplot(cc_dat_fact_m_v2, aes(fill=Severity, x=value)) + 
#           geom_bar(position="stack")+
#           facet_wrap(~variable2, scales="free")+
#           scale_fill_manual(name=pred_var,values=c("blue","red"))+
#           theme_prism()
#         g_stack_cc_factor
#         pdf(paste0(fig_folder,"/barplot_severity_CC_factors_",names(data_list)[i],".pdf"), width = 16, height = 8)
#         g_stack_cc_factor
#         dev.off()
#         write.csv(cc_dat_fact, paste0(fig_folder,"/barplot_severity_CC_factors_",names(data_list)[i],".csv"))
#         
#         # if (sep_tk > 0 ){
#         #   cc_dat_fact_m_v3 <-cc_dat_fact_m[sep_tk,]
#         #   cc_dat_fact_m_v3 <- cc_dat_fact_m_v3[cc_dat_fact_m_v3$value == 1,]
#         # }
#       }
#       if (length(which(is.fact==F))>1){
#         cc_dat_num <- cc_dat[,c(names(cc_dat)[which(is.fact==F)])]
#         cc_dat_num_m <- melt(cc_dat_num, id.vars = "Severity")
#         cc_dat_num_m$variable2 <- cc_dat_num_m$variable
#         cc_dat_num_m$variable2 <- gsub("\\_"," ", cc_dat_num_m$variable2)
#         cc_dat_num_m$variable2 <- gsub("X\\."," ", cc_dat_num_m$variable2)
#         cc_dat_num_m$variable2 <- gsub("\\.\\."," ", cc_dat_num_m$variable2)
#         cc_dat_num_m$variable2 <- gsub("\\."," ", cc_dat_num_m$variable2)
#         g_box_cc_num <- ggplot(cc_dat_num_m, aes(x=Severity, y=value, fill = Severity)) + 
#           geom_boxplot(alpha=0.8)+
#           geom_jitter(shape =21, size =3)+
#           facet_wrap(~variable2, scales="free")+
#           scale_fill_manual(name=pred_var,values=c("blue","red"))+
#           theme_prism()
#         g_box_cc_num
#         pdf(paste0(fig_folder,"/box_severity_CC_numeric_",names(data_list)[i],".pdf"), width = 5, height = 5)
#         g_box_cc_num
#         dev.off()
#         write.csv(cc_dat_num_m, paste0(fig_folder,"/box_severity_CC_numeric_",names(data_list)[i],".csv"))
#       }
#       }
#   }
# }  
# 
# # Compare predictions using imputation and clean data
# library(pROC)
# f_dt <- do.call("rbind",f1_list)
# 
# data_var <- unique(f_dt$data)
# f1_var <- c()
# cf_list <- list()
# for(var in data_var){
#   
#   f_dt_sel <-  f_dt[f_dt$data ==var,]
#   f1_score <- F1_Score(y_pred =f_dt_sel$ypred, y_true =f_dt_sel$ytrue,positive = 1)
#   cf <- confusionMatrix(factor(f_dt_sel$ypred),factor(f_dt_sel$ytrue),positive = "1")
#   print(cf)
#   cf_list[[var]] <- cf$byClass
#   f1_var <- c(f1_var,f1_score)
# }
# 
# cf_dt <- data.frame(do.call("rbind",cf_list))
# cf_dt$Data <- data_var
# 
# # Save the confusion matrix metrics
# write.csv(cf_dt, paste0(fig_folder,"/Confusion_matrix.csv"))
# 
# pred_dt <- do.call("rbind",pred_dt_list)
# 
# data_var <- unique(pred_dt$data)
# roc_var_list  <- list()
# var <- data_var[1]
# for(var in data_var){
#   
#   pred_dt_sel <-  pred_dt[pred_dt$data ==var,]
#   roc_var_dt <- roc(pred_dt_sel$test,pred_dt_sel$pred,plot=TRUE,smooth = F)
#   
#   roc_dt <- data.frame(
#     tpp=roc_var_dt$sensitivities*100, ## tpp = true positive percentage
#     fpp=(1 - roc_var_dt$specificities)*100, ## fpp = false positive precentage
#     data = var,
#     AUC=as.numeric(gsub("Are.*:","",roc_var_dt$auc)))
#   
#   roc_var_list[[var]] <- roc_dt
# }
# 
# roc_comb <-  do.call("rbind",roc_var_list)
# 
# roc_comb$AUC <- as.character(signif(100*roc_comb$AUC, digits = 4))
# 
# auc_text_dt <-  unique(roc_comb[,c("data","AUC")])
# auc_text_dt$fpp <- c(75,75,75,75,75)
# auc_text_dt$tpp <- seq(50,30,-5)
# auc_text_dt$AUC <- paste0(auc_text_dt$data," : ",auc_text_dt$AUC)
# 
# roc_comb[1,]
# roc_comb$data <- gsub("AB\\_","",roc_comb$data)
# roc_comb$data <- gsub("CC_only","CC",roc_comb$data)
# 
# roc_comb[1,]
# auc_text_dt$data <- gsub("AB\\_","",auc_text_dt$data)
# auc_text_dt$data <- gsub("CC_only","CC",auc_text_dt$data)
# auc_text_dt$AUC <- gsub("AB\\_","",auc_text_dt$AUC)
# auc_text_dt$AUC <- gsub("CC_only","CC",auc_text_dt$AUC)
# 
# p_roc <- ggplot(roc_comb)+
#   geom_path(aes(x= fpp,y = tpp,color = data,group = data),alpha = 1,size = 1)+
#   geom_text(data = auc_text_dt,aes(x= fpp,y = tpp,color = data,group = data, label = AUC))+
#   scale_color_brewer(palette = "Set1")+
#   
#   
#   
#   #stat_summary(fun.data = "mean_cl_boot", geom = "line",
#   #             colour = "blue")+
#   xlab("False Positive Percentage")+
#   ylab("True Positive Percentage")+
#   geom_abline(intercept =0 , slope = 1, linetype = 'dashed')+
#   coord_fixed()+
#   theme_prism()+
#   labs(color = "Data")
# 
# print(p_roc)
# write.csv(roc_comb,paste0(fig_folder, "/ROC_curve_Severity.csv"))
# ggsave(paste0(fig_folder, "/ROC_curve_Severity.pdf"), p_roc,width=5, height =5)
# 
# # plots for paper
# results_dir <- "../results/RF_final_Sev/2020-12-16/"
# importance_stl <- read.csv(paste0(results_dir,"importance_severity_AB_STL.csv"))
# importance_stl
# 
# library(ggthemes)
# library(ggprism)
# dat <- importance_stl
# dat$variable <- gsub("\\_"," ", dat$variable)
# dat$variable <- gsub("X\\."," ", dat$variable)
# dat$variable <- gsub("\\.\\."," ", dat$variable)
# dat$variable <- gsub("\\."," ", dat$variable)
# dat <-  dat[order(dat$importance,decreasing = T),]
# dat <- dat[which(dat$importance>0),]
# g_imp_stl <- ggplot(dat, aes(x=reorder(variable,importance), y=importance, fill = direction))+
#   geom_bar(stat="identity", position="dodge", color = "black")+ coord_flip()+
#   ylab("Variable Importance")+
#   scale_fill_manual(values = c("blue","red"))+
#   theme_prism()+
#   xlab("Variables")
# g_imp_stl
# 
# importance_tng <- read.csv(paste0(results_dir,"importance_severity_AB_TNG.csv"))
# importance_tng
# 
# dat <- importance_tng
# dat$variable <- gsub("\\_"," ", dat$variable)
# dat$variable <- gsub("X\\."," ", dat$variable)
# dat$variable <- gsub("\\.\\."," ", dat$variable)
# dat$variable <- gsub("\\."," ", dat$variable)
# dat <-  dat[order(dat$importance,decreasing = T),]
# dat <- dat[which(dat$importance>0),]
# g_imp_tng <- ggplot(dat, aes(x=reorder(variable,importance), y=importance, fill = direction))+
#   geom_bar(stat="identity", position="dodge", color = "black")+ coord_flip()+
#   ylab("Variable Importance")+
#   scale_fill_manual(values = c("blue","red"))+
#   theme_prism()+
#   xlab("Variables")
# 
# g_imp_tng
# 
# library(ggpubr)
# library(ggrepel)
# combined_1 <- ggarrange(g_imp_stl, g_imp_tng, nrow = 1, ncol = 2)
# #combined_1 <- ggarrange(p_roc, pimp_plot, nrow = 1, ncol = 2, align = 'h')
# pdf(paste0(fig_folder, "/combined_stl_tng_importance.pdf"), combined_1, width=20, height =10)
# combined_1
# dev.off()
# 
# 
# # lime frequency stl
# stl_lime_freq <- read.csv(paste0(fig_folder,"/LIME_severity_frequency_list_",names(data_list)[2],".csv"))
# stl_lime_freq$V5<-"moderate"
# stl_lime_freq$V5[stl_lime_freq$V4==1]<-"severe"
# 
# g_freq_lime_stl <- ggplot()+
#   geom_point(data=stl_lime_freq[stl_lime_freq$V3>0.2,], aes(x=V3, y=V2, fill=V5), shape= 21, size=4,alpha=1)+
#   geom_point(data=stl_lime_freq[stl_lime_freq$V3<=0.2,], aes(x=V3, y=V2), fill="grey",shape= 21, size=4,alpha=0.5)+
#   geom_text_repel(data=stl_lime_freq[stl_lime_freq$V3>0.2,], aes(x=V3, y=V2, label=V1))+
#   ylab("Frequency correct labeling")+
#   scale_fill_manual(values = c("blue","red"))+
#   theme_prism()+
#   xlab("Frequency rule selected")
# g_freq_lime_stl
# 
# # lime frequency stl
# tng_lime_freq <- read.csv(paste0(fig_folder,"/LIME_severity_frequency_list_",names(data_list)[3],".csv"))
# tng_lime_freq$V5<-"moderate"
# tng_lime_freq$V5[stl_lime_freq$V4==1]<-"severe"
# 
# g_freq_lime_tng <- ggplot()+
#   geom_point(data=tng_lime_freq[tng_lime_freq$V3>0.2,], aes(x=V3, y=V2, fill=V5), shape= 21, size=4,alpha=1)+
#   geom_point(data=tng_lime_freq[tng_lime_freq$V3<=0.2,], aes(x=V3, y=V2), fill="grey",shape= 21, size=4,alpha=0.5)+
#   geom_text_repel(data=tng_lime_freq[tng_lime_freq$V3>0.2,], aes(x=V3, y=V2, label=V1))+
#   ylab("Frequency correct labeling")+
#   scale_fill_manual(values = c("blue","red"))+
#   theme_prism()+
#   xlab("Frequency rule selected")
# g_freq_lime_tng
# 
# combined_2 <- ggarrange(g_freq_lime_stl, g_freq_lime_tng, nrow = 1, ncol = 2)
# #combined_1 <- ggarrange(p_roc, pimp_plot, nrow = 1, ncol = 2, align = 'h')
# pdf(paste0(fig_folder, "/combined_stl_tng_lime.pdf"), combined_2, width=15, height =8)
# combined_2
# dev.off()
# 
