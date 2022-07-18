################################### PACKAGES, HELPER FUNCTIONS, LIST OF TRAITS ###################################
library(ggplot2)
library(ggrepel)
library(readxl)
library(googlesheets4)
library(stringr)
library(latex2exp)
library(Cairo)
library(ggtext)
library(stringi)
library(patchwork)
library(reshape2)
library(mltools)
library(data.table)
library(dplyr)

theme_bhr <- function(){ 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15,  color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          legend.position = "None",
          legend.direction = "horizontal",
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(fill = "white"))
}

null_color = "gray70"
colors_emphasis = c("Diabetes" = "#C95693",
                    "Townsend" = null_color,
                    "Asthma" = null_color,
                    "Osteoarthritis" = null_color,
                    "FIS" = null_color,
                    "Birth weight" = null_color,
                    "ID time" = null_color,
                    "Neuroticism" = "#7B03D1",
                    "Alcohol frequency" = null_color,
                    "BMI" = "#D9B407",
                    "Cancer" = "#7B8FEA",
                    "LDL" = "#B20C00",
                    "Platelet count" = "#F36E04",
                    "ALT" = null_color,
                    "FEV1" = null_color,
                    "Calcium" = null_color,
                    "HbA1c" = null_color,
                    "IGF-1" = null_color,
                    "SBP" = "#67D720",
                    "Height" = "#227E18",
                    "BMD" = null_color,
                    "RBC count" = "#43BBA7",
                    "NA" = null_color)

emphasis_traits = names(colors_emphasis)[colors_emphasis != null_color]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#########################################################################################################
###Figure 1A
#path = "/Users/daniel/Desktop/rare_h2/ms/genebass_summary_statistics_null/"
#sumstats_plof <- readRDS(paste0(path,"bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1.ms.munged.Rds"))
#sumstats_plof = sumstats_plof[sumstats_plof$phenotype_key == "30780NA",c("gene",      "N", "burden_score",    "w_t_beta", "overdispersion", "trait_type", "phenotype_key", "chromosome", "gene_position")]
#write.csv(sumstats_plof, "outputs/sumstats_plof_ldl.csv")
#sumstats_synonymous <- readRDS(paste0(path,"bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0_high1e-05_group1.ms.munged.Rds"))
#sumstats_synonymous = sumstats_synonymous[sumstats_synonymous$phenotype_key == "30780NA",c("gene",      "N", "burden_score",    "w_t_beta", "overdispersion", "trait_type", "phenotype_key", "chromosome", "gene_position")]
#write.csv(sumstats_synonymous, "outputs/sumstats_synonymous_ldl.csv")

# sumstats_plof = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST1"))
# sumstats_synonymous = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST2"))
# 
# sumstats_plof$gamma_sq = sumstats_plof$w_t_beta^2/sumstats_plof$burden_score
# sig_chisq = (sumstats_plof$N - 1) * sumstats_plof$gamma_sq > qchisq(p = 1 - (0.05/nrow(sumstats_plof)), df = 1)
# plof_sumstats <- sumstats_plof[!sig_chisq & sumstats_plof$burden_score < 0.001,]
# bins = bin_data(plof_sumstats$burden_score, bins=25, binType = "quantile",returnDT = TRUE)
# plof_sumstats$group  = bin_data(plof_sumstats$burden_score, bins=25, binType = "quantile",)
# plof_sumstats_gammasq_bin <- sapply(levels(plof_sumstats$group), function(x) mean(plof_sumstats$gamma_sq[plof_sumstats$group == x]))
# plof_sumstats_normw_bin <- sapply(levels(plof_sumstats$group), function(x) mean(plof_sumstats$burden_score[plof_sumstats$group == x]))
# 
# sumstats_synonymous$gamma_sq = sumstats_synonymous$w_t_beta^2/sumstats_synonymous$burden_score
# sig_chisq = (sumstats_synonymous$N - 1) * sumstats_synonymous$gamma_sq > qchisq(p = 1 - (0.05/nrow(syn_sumstats)), df = 1)
# syn_sumstats <- sumstats_synonymous[sumstats_synonymous$burden_score < 0.001,]
# syn_sumstats$gammasq = syn_sumstats$w_t_beta^2/syn_sumstats$burden_score
# syn_sumstats$group  = bin_data(syn_sumstats$burden_score, bins=25, binType = "quantile",)
# syn_sumstats_gammasq_bin <- sapply(levels(syn_sumstats$group), function(x) mean(syn_sumstats$gammasq[syn_sumstats$group == x]))
# syn_sumstats_normw_bin <- sapply(levels(syn_sumstats$group), function(x) mean(syn_sumstats$burden_score[syn_sumstats$group == x]))
# 
# n = 359350
# 
# plot <- ggplot()+
#   geom_point(mapping = aes(x = plof_sumstats_normw_bin, y = plof_sumstats_gammasq_bin/(1/n)), size = 3, color = "#BC321B")+
#   geom_smooth(mapping = aes(x = plof_sumstats_normw_bin, y = plof_sumstats_gammasq_bin/(1/n)),method = "lm", se = FALSE, color = "#BC321B")+
#   geom_point(mapping = aes(x = syn_sumstats_normw_bin, y = syn_sumstats_gammasq_bin/(1/n)),size = 3, color = "grey70")+
#   geom_smooth(mapping = aes(x = syn_sumstats_normw_bin, y = syn_sumstats_gammasq_bin/(1/n)),method = "lm", color = "grey70", se = FALSE)+
#   labs(x = "Mean burden score", y = TeX(r'(Mean $\chi^2$ statistic)'))+
#   theme_bhr()+
#   scale_x_continuous(labels = function(x) format(x,digits = 1,scientific = TRUE,drop0trailing = TRUE))+
#   theme(plot.margin = margin(2,15,2,2))+
#   annotate("text", x = 6.7e-4, y = 3.7e-6/(1/n), size = 6, color = "#BC321B", label = "pLoF")+
#   annotate("text", x = 8.2e-4, y = 3e-6/(1/n), size = 6, color = "grey70", label = "Synonymous")+
#   annotate("text", x = 5.3e-4, y = 3.3e-6/(1/n), size = 6, color = "black", label = "Slope = Squared \n per-allele effect size")


###Figure 1C (simulations)
data = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST3"))
data$tag = paste0(data$model,"_",data$af_group)
data = data[data$tag %in% c("realistic_1", "realistic_2", "strong_selection_1", "small_N_1"),]
data$tag_display = ifelse(data$tag == "realistic_1", "Realistic",
                          ifelse(data$tag == "realistic_2", "High MAF",
                                 ifelse(data$tag == "small_N_1", "Small N",
                                        ifelse(data$tag == "strong_selection_1", "Strong selection", NA))))
data$tag_display = factor(data$tag_display, levels = c("Realistic", "Small N", "High MAF", "Strong selection"))
data$true_h2 = factor(data$true_h2, levels = c("0.005", "0"))
data$true_h2_label <- ifelse(data$true_h2 == 0, "Null (h2 = 0%)","Non-null (h2 = 0.5%)")

ggplot(data, aes(x = factor(true_h2), y = estimated_h2, fill = factor(true_h2_label)))+
  facet_grid(~tag_display, switch = "both")+
  geom_boxplot(size = 1)+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 0.005, linetype = "dashed")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(fill = "white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x = "", y = "Burden heritability", title = "")+
  scale_fill_manual(values = c("#BC321B", "gray70"))+
  scale_y_continuous(labels = function(x) paste0(x*100, "%"))

###Figure 2A (variant frequencies)
data = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST5"))[,c(1,3,4)]

data_common = data[data$Frequency_name == "common",]
data_common$fraction <- data_common$Count/sum(data_common$Count)
data_ur = data[data$Frequency_name == "ultra_rare",]
data_ur$fraction <- data_ur$Count/sum(data_ur$Count)
data_rare = data[data$Frequency_name == "rare",]
data_rare = aggregate(data_rare$Count, by=list(Functional_group=data_rare$Functional_group, Frequency_name = data_rare$Frequency_name), FUN=sum)
data_rare$fraction = data_rare$x/sum(data_rare$x)
colnames(data_rare) <- c("Functional_group", "Frequency_name", "Count", "fraction")
data = rbind(data_common, data_rare, data_ur)

data$display_names <- ifelse(data$Frequency_name == "common", "Common", 
                             ifelse(data$Frequency_name == "rare", "Rare",
                                    ifelse(data$Frequency_name == "ultra_rare", "Ultra-Rare", NA)))

data$Functional_group <- ifelse(data$Functional_group == "pLoF", "pLoF",
                                ifelse(data$Functional_group == "missense_benign", "Missense: Benign",
                                       ifelse(data$Functional_group == "missense_notbenign", "Missense: Pathogenic",
                                       ifelse(data$Functional_group == "synonymous", "Synonymous", NA))))

data$display_names <- factor(data$display_names, levels=c("Ultra-Rare",
                                              "Rare",
                                              "Common"))

data$Functional_group <- factor(data$Functional_group, levels=c("pLoF",
                                                 "Missense: Pathogenic",
                                                 "Missense: Benign",
                                                 "Synonymous"))

my.labels <- c(paste0("Ultra rare\n(",format(sum(data_ur$Count), big.mark = ","), " variants)"),
               paste0("Rare\n(",format(sum(data_rare$Count), big.mark = ","), " variants)"), 
               paste0("Common\n(",format(sum(data_common$Count), big.mark = ","), " variants)"))

Fig2A <- ggplot(data, aes(x = display_names, y = fraction, fill = Functional_group))+
  geom_bar(position = "stack", stat = "identity", color = "black", size = 1)+
  labs(x = "", y = "Proportion of variants", fill = "Variant group")+
  theme(text=element_text(family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20, color = "black"),
        axis.title=element_text(size=20,  color = "black"),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        #legend.position = c(1.2,-0.15),
        legend.position = "None",
        legend.direction = "horizontal",
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(fill = "white"))+
  scale_fill_manual(values=c("#BC321B", "darkorange1", "goldenrod1", "gray70"))+
  geom_hline(yintercept = 0)+
  scale_x_discrete(labels= my.labels)+
  scale_y_continuous(breaks = c(0, 0.5, 1))


###Figure 2B (heritability by ff group)
data = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST6"))
phenotype_file <- read.csv("~/rv_h2/reference_files/ms_phenotype_sheet.csv")
phenotype_file = phenotype_file[phenotype_file$phenotype_core == 1,]
data = data[data$phenotype_key %in% phenotype_file$phenotype_key,]
data$frequency <- substrRight(data$summary_statistic, 6)
data$functional_group <- do.call(rbind.data.frame, stri_split_fixed(str = substring(data$summary_statistic, 42), pattern = "_", n = 2))[,1]
data_ur <- data[data$frequency == "group1",c("phenotype_key", "bhr_h2", "functional_group")]
data_ur$frequency <- "Ultra rare"

data_rare <- data[data$frequency %in% c("group2", "group3"),c("phenotype_key", "bhr_h2", "frequency", "functional_group")]
data_rare = aggregate(data_rare$bhr_h2, by=list(functional=data_rare$functional_group, trait=data_rare$phenotype_key), FUN=sum)
data_rare = data_rare[c("trait", "x", "functional")]
colnames(data_rare) <- c("phenotype_key", "bhr_h2", "functional_group")
data_rare$frequency <- "Rare"
data = rbind(data_ur, data_rare)
data$bhr_h2_percent <- data$bhr_h2*100

data$functional_group <- factor(data$functional_group, levels=c("pLoF",
                                                                "missense-notbenign",
                                                                "missense-benign",
                                                                "synonymous"))
data$frequency <- factor(data$frequency, levels = c("Ultra rare", "Rare"))

Fig2B <- ggplot(data, aes(y = bhr_h2, fill = functional_group))+
  geom_boxplot(size = 1)+
  facet_grid(~frequency)+
  geom_hline(yintercept = 0, size = 1)+
  coord_cartesian(ylim = c(-0.1/100,1/100))+
  labs(x = "Variant group", y = "Burden heritability", fill = "Variant category")+
  theme(text=element_text(family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20,  color = "black"),
        axis.title=element_text(size=20,  color = "black"),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20),
        legend.position = "None",
        legend.direction = "horizontal",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(fill = "white"))+
  scale_fill_manual(values=c("#BC321B", "darkorange1", "goldenrod1", "gray70"))+
  scale_y_continuous(breaks = c(0, 0.005, .01), labels = function(x) paste0(x*100, "%"))


#Figure 2C: Comparison with LDSC heritability
ldsc_h2 = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY", 
                                     sheet = "ST8"))
bhr_h2 = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", 
                                              sheet = "ST7"))
bhr_h2$display_name = ifelse(bhr_h2$display_name == "Standing height", "Height",bhr_h2$display_name)
bhr_h2$emphasis = bhr_h2$display_name
bhr_h2$emphasis[!(bhr_h2$emphasis %in% emphasis_traits)] <- ""
bhr_h2 = left_join(bhr_h2, ldsc_h2)

Fig2C <- ggplot(bhr_h2, mapping = aes(x = ldsc_h2, y = aggregated_h2, 
                              ymin = aggregated_h2 - aggregate_h2_se, ymax = aggregated_h2  + aggregate_h2_se,
                              xmin = ldsc_h2 - ldsc_h2_se, xmax = ldsc_h2 + ldsc_h2_se,
                              label = emphasis,
                              color = emphasis))+
  geom_point(size = 3)+
  geom_smooth(method = "lm", color = "black", se = FALSE)+
  geom_errorbar(width = 0)+
  geom_errorbarh(height = 0)+
  geom_text_repel(max.overlaps = 200, force = 100, size = 5,point.padding = 1,seed = 10332134)+
  labs(x = "Common-variant heritability (SE)", y = "Total burden heritability (SE)")+
  scale_color_manual(values = colors_emphasis)+
  guides(color = "none")+
  theme_bhr()

####Figure 2D
phenotype_file <- data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST4"))
phenotype_file = phenotype_file[phenotype_file$phenotype_core == 1,]

bhr_h2 = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST6"))
bhr_h2$emphasis = bhr_h2$display_name
bhr_h2$emphasis[!(bhr_h2$emphasis %in% emphasis_traits)] <- ""
bhr_h2 = bhr_h2[bhr_h2$phenotype_key %in% phenotype_file$phenotype_key,]

bhr_h2_core = bhr_h2[bhr_h2$summary_statistic %in% c("bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",
                                                                                  "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0_high1e-05_group1"),]

Fig2D = ggplot(data = bhr_h2_core,
       mapping = aes(x = n_bhr*intercept,
                     y = lambda_gc,
                     color = summary_statistic))+
  geom_point(size = 3)+
  geom_abline()+
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_manual(values = c("#BC321B","gray70"))+
  guides(color = "none")+
  theme_bhr()+
  labs(x = "Scaled Intercept", y = "Lambda GC")+
  annotate("text", x = 1.02, y = 1.4, size = 5, color = "#BC321B", label = "Ultra-Rare pLoF")+
  annotate("text", x = 1.02, y = 1.37, size = 5, color = "gray70", label = "Ultra-Rare Synonymous")+
  scale_y_continuous(breaks=seq(1,1.4,0.1))+
  scale_x_continuous(breaks=seq(1,1.4,0.1))

plot_grid(Fig2A, Fig2B, Fig2C, Fig2D, ncol = 2, labels = c('A', 'B', 'C', 'D'), align = "h")

#####Figure 3A
ordered_sigdf = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY", 
                                                     sheet = "ST9"))
ordered_sigdf$labelgene[ordered_sigdf$proportion_bhr_explained < 0.05] <- ""

colors_emphasis_Fig3A = colors_emphasis
colors_emphasis_Fig3A[colors_emphasis_Fig3A == "gray70"] <- "black"


Fig3A <- ggplot(data = ordered_sigdf)+
  geom_col(mapping = aes(x = factor(display_name, levels = unique(display_name[order(frac_sig_winnerscursecorr)])), 
                         y = proportion_bhr_explained,
                         color = display_name),
           position = position_stack(),
           fill = "white")+
  geom_errorbar(mapping = aes(x = factor(display_name, levels = unique(display_name[order(frac_sig_winnerscursecorr)])),
                              ymin = frac_sig_winnerscursecorr - frac_sig_se, 
                              ymax = frac_sig_winnerscursecorr + frac_sig_se ), width = 0,
                color = "gray70",
                alpha = 0.5)+
  geom_text(mapping = aes(x = factor(display_name, levels = unique(display_name[order(frac_sig_winnerscursecorr)])),
                          y = proportion_bhr_explained,
                          label = labelgene),
            position = position_stack(0.5),
            size = 3.5)+
  coord_flip(ylim = c(0, 0.65))+
  labs(y = "Fraction of burden heritability explained by significant genes (SE)", x = "")+
  theme_bhr()+
  theme(axis.text.y = element_text(color = as.vector(colors_emphasis_Fig3A)[match(unique(ordered_sigdf$display_name[order(ordered_sigdf$frac_sig_winnerscursecorr)]),names(colors_emphasis_Fig3A))]))+
  ylim(-0.02,1)+
  scale_color_manual(values = colors_emphasis_Fig3A)+
  guides(fill = "none",  color = "none")
ggsave("3A.pdf", Fig3A,width = 10, height =5, dpi = 1200, device = cairo_pdf)

#Figure 3B
ordered_sigdf_common = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY", 
                                                            sheet = "ST11"))

ordered_sigdf_common$display_name = phenotype_file$display_name[match(ordered_sigdf_common$phenotype_key,phenotype_file$phenotype_key)]

colors_emphasis_Fig3B = colors_emphasis
colors_emphasis_Fig3B[colors_emphasis_Fig3B == "gray70"] <- "black"

Fig3B <- ggplot(data = ordered_sigdf_common)+
  geom_col(mapping = aes(x = factor(display_name, levels = unique(display_name[order(fraction_HESS_sig)])), 
                         y = fraction_HESS,
                         color = display_name),
           position = position_stack(),
           fill = "white")+
  coord_flip(ylim = c(0, 0.65))+
  labs(y = "Fraction of common-variant heritability explained by significant SNPs", x = "")+
  theme_bhr()+
  theme(axis.text.y = element_text(color = as.vector(colors_emphasis_Fig3B)[match(unique(ordered_sigdf_common$display_name[order(ordered_sigdf_common$fraction_HESS_sig)]),names(colors_emphasis_Fig3B))]))+
  ylim(-0.02,1)+
  scale_color_manual(values = colors_emphasis_Fig3B)+
  guides(fill = "none",  color = "none")

ggsave("3B.pdf", Fig3B,width = 10, height =5, dpi = 1200, device = cairo_pdf)

# ###Figure 3x (SUPPLEMENT)
# HESS_output = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY", 
#                                                      sheet = "ST12"))
# 
# gnomad_information <- data.frame(fread("reference_files/gnomad.v2.1.1.lof_metrics.by_gene.txt",
#                                        select = c("gene","gene_id", "chromosome", "start_position",	"end_position", "pLI")))
# gnomad_information = gnomad_information[gnomad_information$chromosome != "X",]
# gnomad_information = gnomad_information[gnomad_information$chromosome != "Y",]
# gnomad_information_counts = data.frame(table(gnomad_information$gene))
# gnomad_information_counts = gnomad_information_counts[gnomad_information_counts$Freq == 1,]
# gene_information <- gnomad_information[gnomad_information$gene %in% gnomad_information_counts$Var1,]
# gene_information$midpoint <- (gene_information$start_position + gene_information$end_position)/2
# 
# 
# get_corresponding_HESSh2_sigdf <- function(assn) {
#   print(assn)
#   gene = ordered_sigdf$gene[assn] 
#   hessoutput = HESS_output[HESS_output$trait == ordered_sigdf$phenotype_key[assn],]
#   start = gene_information$start_position[gene_information$gene_id == gene]
#   end = gene_information$end_position[gene_information$gene_id == gene]
#   chr = as.numeric(gene_information$chromosome[gene_information$gene_id == gene])
#   midpoint = as.numeric(gene_information$midpoint[gene_information$gene_id == gene])
#   block = which(hessoutput$chr == chr & hessoutput$start < midpoint & hessoutput$end > midpoint)
#   
#   print(block)
#   if (length(hessoutput$local_h2g[block]) > 0) {return(hessoutput$local_h2g[block])} else {return(NA)}
# }
# 
# ordered_sigdf$corresponding_HESS_h2 = sapply(1:length(ordered_sigdf$phenotype_key),get_corresponding_HESSh2_sigdf)
# 
# ordered_sigdf$total_HESS_h2 = ordered_sigdf_common$HESS_h2[match(ordered_sigdf$phenotype_key,ordered_sigdf_common$phenotype_key)]
# 
# ordered_sigdf$corresponding_HESS_proportion = ordered_sigdf$corresponding_HESS_h2/ordered_sigdf$total_HESS_h2
# ordered_sigdf = ordered_sigdf[ordered_sigdf$proportion_bhr_explained > 0,]
# Fig3C <- ggplot(data = ordered_sigdf, 
#        mapping = aes(x = corresponding_HESS_proportion, 
#                      y = proportion_bhr_explained,
#                      color = display_name))+
#   geom_point(size = 3)+
#   geom_abline()+
#   geom_abline(intercept = 1, color = "grey50", linetype = "dashed")+
#   geom_abline(intercept = 2, color = "grey50", linetype = "dashed")+
#   geom_abline(intercept = -1, color = "grey50", linetype = "dashed")+
#   geom_abline(intercept = -2, color = "grey50", linetype = "dashed")+
#   theme_bhr()+
#   scale_color_manual(values = colors_emphasis)+
#   scale_x_log10()+
#   scale_y_log10()+
#   labs(x = "Fraction of common-variant heritability",
#        y = "Fraction of burden heritability")+
#   annotate("text", x = 0.2, y = 0.35, size = 6, color = "black", label = "1x")+
#   annotate("text", x = 0.0015, y = 0.35, size = 6, color = "black", label = "100x")+
#   annotate("text", x = 0.2, y = 4e-3, size = 6, color = "black", label = "0.01x")
# 
# ggsave("3C.pdf", Fig3C,width = 5, height =5, dpi = 1200, device = cairo_pdf)

###Figure 3C
bhr_sig_associations = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST9"))
bhr_sig_associations = bhr_sig_associations[!duplicated(bhr_sig_associations[ , c("phenotype_key")]),][c("phenotype_key", "frac_sig_winnerscursecorr", "frac_sig_se")]

amm_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST13"))
amm_enrichments$trait_set_key <- paste0(amm_enrichments$phenotype_key,"_",amm_enrichments$gene_set_display)
amm_enrichments = amm_enrichments[amm_enrichments$trait_set_key %in% c("30620NA_Sig_ALT",
                                                                       "21001NA_Sig_BMI",
                                                                       "30680NA_Sig_calcium",
                                                                       "2453NA_Sig_cancer",
                                                                       "3063NA_Sig_FEV1",
                                                                       "30750NA_Sig_HbA1c",
                                                                       "50NA_Sig_height",
                                                                       "30770NA_Sig_IGF1",
                                                                       "30780NA_Sig_LDL",
                                                                       "30080NA_Sig_platelet_count",
                                                                       "30010NA_Sig_RBC_count"),c("phenotype_key",  "fraction_h2", "fraction_h2_se")]
colnames(amm_enrichments) <- c("phenotype_key",  "amm_fraction_h2", "amm_fraction_h2_se")
merged <- merge(bhr_sig_associations, amm_enrichments, by = "phenotype_key")
phenotypes = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST4"))[c("phenotype_key", "display_name")]
merged = merge(merged, phenotypes, by = "phenotype_key")
colnames(merged) = c("phenotype_key", "bhr_fraction_h2", "bhr_fraction_h2_se", "amm_fraction_h2", "amm_fraction_h2_se", "display_name") 
merged$bhr_fraction_h2_se_display = ifelse(merged$display_name %in% emphasis_traits, merged$bhr_fraction_h2_se, 0)
merged$amm_fraction_h2_se_display = ifelse(merged$display_name %in% emphasis_traits, merged$amm_fraction_h2_se, 0)

merged$display_only = merged$display_name
merged$display_only[!(merged$display_name %in% emphasis_traits)] <- ""

Fig3C <- ggplot(merged, aes(x = amm_fraction_h2, y = bhr_fraction_h2,
                 xmin = amm_fraction_h2 - amm_fraction_h2_se_display, xmax = amm_fraction_h2 + amm_fraction_h2_se_display,
                 ymin = bhr_fraction_h2 - bhr_fraction_h2_se_display, ymax = bhr_fraction_h2 + bhr_fraction_h2_se_display,
                 label = display_only,
                 color = display_name))+
  geom_point(size = 3)+
  geom_errorbar(width = 0)+
  geom_errorbarh(height = 0)+
  labs(x = "Fraction of common-variant heritability (SE)", y = "Fraction of burden heritability (SE)")+
  geom_abline(slope = 1, linetype = "dashed")+
  geom_text_repel(force_pull = .1, box.padding = 1.25, size = 5, family = "Arial", max.overlaps = 1000)+
  theme_bhr()+
  scale_color_manual(values=colors_emphasis)+
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6))

###Figure 3D
amm_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST13"))
amm_enrichments$trait_set_key <- paste0(amm_enrichments$phenotype_key,"_",amm_enrichments$gene_set_display)
amm_enrichments = amm_enrichments[amm_enrichments$trait_set_key %in% c("2453NA_Sig_cancer",
                                                                       "2453NA_Cancer_TSG",
                                                                       "2453NA_Cancer_oncogenes"),
                                  c("phenotype_key",  "gene_set_display",  "fraction_h2", "fraction_h2_se")]
amm_enrichments$gene_set_display <- ifelse(amm_enrichments$gene_set_display == "Cancer_TSG", "cosmic_tsg",
                                           ifelse(amm_enrichments$gene_set_display == "Cancer_oncogenes", "cosmic_oncogene",
                                                  ifelse(amm_enrichments$gene_set_display == "Sig_cancer", "siggene_2453NA", NA)))
amm_enrichments = amm_enrichments[c("gene_set_display", "fraction_h2", "fraction_h2_se")]
amm_enrichments$model <- "Common-variant heritability"
colnames(amm_enrichments) <- c("gene_set", "fraction_h2", "fraction_h2_se", "model")

bhr_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST14"))
bhr_enrichments = bhr_enrichments[bhr_enrichments$display_name == "Cancer" & bhr_enrichments$gene_set_display %in% c("Cancer_oncogenes", "Cancer_TSG", "Sig_cancer"),]
bhr_enrichments = bhr_enrichments[c("gene_set", "fraction_h2", "fraction_h2_se")]
bhr_enrichments$model <- "Burden heritability"

merged <- rbind(amm_enrichments, bhr_enrichments)
merged$gene_set_display_figure <- ifelse(merged$gene_set == "cosmic_oncogene", "Oncogenes\n (n = 101)",
                                                  ifelse(merged$gene_set == "siggene_2453NA", "Exome\nsignificant\n (n = 7)",
                                                         ifelse(merged$gene_set == "cosmic_tsg", "Tumor\nsuppressor\n (n = 172)", NA)))
merged$model <- factor(merged$model, levels=c("Common-variant heritability", "Burden heritability"))

Fig3D <- ggplot(merged, aes(x = reorder(gene_set_display_figure, fraction_h2), y = fraction_h2, ymin = fraction_h2 - fraction_h2_se, ymax = fraction_h2 + fraction_h2_se,
                       fill = model))+
  geom_col(color = "black", position=position_dodge())+
  geom_errorbar(width = 0, size = 1, color = "black", position=position_dodge(.9))+
  coord_flip()+
  geom_hline(yintercept = 0)+
  theme_bhr()+
  labs(x = "", y = "Fraction of cancer heritability (SE)")+
  scale_fill_manual(values = c("grey70", "#7B8FEA"))+
  annotate("text", x = 1.2, y = 0.4, size = 6, color = "#7B8FEA", label = "Burden heritability")+
  annotate("text", x = 1, y = 0.4, size = 6, color = "grey70", label = "Common-variant heritability")

combined <- plot_grid(Fig3A, Fig3C, Fig3B, Fig3D, ncol = 2, labels = c('A', 'C', 'B', 'D'), align = 'hv')

###Figure 4A (common vs rare tissue enrichments)
bhr_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST14"))
bhr_enrichments$set_trait <- paste0(bhr_enrichments$gene_set_display, "_",bhr_enrichments$display_name)

bhr_enrichments = bhr_enrichments[bhr_enrichments$set_trait %in% c("Cortex_Alcohol frequency", "Cortex_Birth weight", "Cortex_BMI","Cortex_FIS","Cortex_ID time","Cortex_Neuroticism","Cortex_Townsend",
                                                                   "GABAergic_neuron_Alcohol frequency", "GABAergic_neuron_Birth weight", "GABAergic_neuron_BMI","GABAergic_neuron_FIS","GABAergic_neuron_ID time","GABAergic_neuron_Neuroticism","GABAergic_neuron_Townsend",
                                                                   "Glutamatergic_neurons_Alcohol frequency", "Glutamatergic_neurons_Birth weight", "Glutamatergic_neurons_BMI","Glutamatergic_neurons_FIS","Glutamatergic_neurons_ID time","Glutamatergic_neurons_Neuroticism","Glutamatergic_neurons_Townsend",
                                                                   "Erythrocyte_HbA1c","Erythrocyte_RBC count",
                                                                   "Hepatocyte_ALT","Hepatocyte_IGF-1","Hepatocyte_LDL",
                                                                   "Liver_ALT","Liver_IGF-1","Liver_LDL",
                                                                   "Megakaryocytes_Platelet count",
                                                                   "Whole_blood_Platelet count","Whole_blood_RBC count"),]
bhr_enrichments = bhr_enrichments[c("gene_set_display", "display_name", "set_trait", "enrichment", "enrichment_se")]
colnames(bhr_enrichments) <- c("gene_set_display", "display_name"          ,            "set_trait", "bhr_enrichment", "bhr_enrichment_se")


amm_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST13"))
amm_enrichments$set_trait <- paste0(amm_enrichments$gene_set_display, "_",amm_enrichments$display_name)
amm_enrichments = amm_enrichments[amm_enrichments$set_trait %in% bhr_enrichments$set_trait,]
amm_enrichments = amm_enrichments[c("set_trait", "Enrichment_mean", "Enrichment_se")]
colnames(amm_enrichments) <- c("set_trait", "amm_enrichment_mean", "amm_enrichment_se")
merged <- merge(bhr_enrichments, amm_enrichments, by = "set_trait")
merged = merged[merged$bhr_enrichment < 10,]
merged = merged[merged$amm_enrichment_mean < 10,]
merged = merged[(merged$bhr_enrichment - 1.96*merged$bhr_enrichment_se) > 0,]
merged = merged[(merged$amm_enrichment_mean - 1.96*merged$amm_enrichment_se) > 0,]
merged$display_label <- ifelse(merged$set_trait == "Megakaryocytes_Platelet count", "Trait: Platelet count\nGene set: Platelets",
                               ifelse(merged$set_trait == "Liver_LDL", "Trait: LDL\nGene set: Liver",
                                      ifelse(merged$set_trait == "Erythrocyte_RBC count", "Trait: RBC Count\nGene set: RBCs",
                                             ifelse(merged$set_trait == "Glutamatergic_neurons_BMI", "Trait: BMI\nGene set: Glutamatergic neurons",""))))

merged$bhr_enrichment_se_display = ifelse(merged$display_label != "", merged$bhr_enrichment_se, 0)
merged$amm_enrichment_se_display = ifelse(merged$display_label != "",  merged$amm_enrichment_se, 0)

null_color = "gray70"
colors_emphasis_f4 = c("Glutamatergic_neurons_BMI" = "#D9B407",
                       "Erythrocyte_RBC count" = "#43BBA7",
                       "Liver_LDL" = "#B20C00",
                       "Megakaryocytes_Platelet count" = "#F36E04",
                       "Cortex_Alcohol frequency" = null_color,
                       "Cortex_BMI"= null_color,
                       "Cortex_FIS"= null_color,
                       "Cortex_ID time"= null_color,
                       "Cortex_Neuroticism"= null_color,
                       "Cortex_Townsend"= null_color,
                       "Erythrocyte_HbA1c"= null_color,
                       "GABAergic_neuron_Alcohol frequency"= null_color,
                       "GABAergic_neuron_Birth weight"= null_color,
                       "GABAergic_neuron_BMI"= null_color,
                       "GABAergic_neuron_ID time"= null_color,
                       "GABAergic_neuron_Neuroticism"= null_color,
                       "GABAergic_neuron_Townsend"= null_color,
                       "Glutamatergic_neurons_Alcohol frequency"= null_color,
                       "Glutamatergic_neurons_Birth weight"= null_color,
                       "Glutamatergic_neurons_FIS"= null_color,
                       "Glutamatergic_neurons_ID time"= null_color,
                       "Glutamatergic_neurons_Neuroticism"= null_color,
                       "Glutamatergic_neurons_Townsend"= null_color,
                       "Liver_ALT"= null_color,
                       "Liver_IGF-1"= null_color,
                       "Whole_blood_Platelet count"= null_color)

f4a <- ggplot(merged, aes(x = amm_enrichment_mean, y = bhr_enrichment, 
                   xmin = amm_enrichment_mean - amm_enrichment_se_display, xmax = amm_enrichment_mean + amm_enrichment_se_display,
                   ymin = bhr_enrichment - bhr_enrichment_se_display, ymax = bhr_enrichment + bhr_enrichment_se_display,
                   color = set_trait,  label = display_label))+
  geom_errorbar(width = 0)+
  geom_errorbarh(height = 0)+
  geom_point(size = 5)+
  ggrepel::geom_text_repel(force_pull = .2, box.padding = 4, size = 5, family = "Arial", max.overlaps = 100)+
  labs(x = "Common-variant enrichment (SE)", y = "Burden heritability enrichment (SE)", color='Gene set')+
  theme(text=element_text(family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"),
        legend.position = "None",
        legend.direction = "horizontal",
        strip.text.x = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "white"))+
  geom_abline(slope = 1, linetype = "dashed")+
  scale_x_continuous(breaks = seq(1, 8, 1))+
  scale_y_continuous(breaks = seq(1, 8, 1))+
  scale_color_manual(values=colors_emphasis_f4)+
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 8))

###Figure 4B (common vs. rare constraint enrichment)
amm_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST13"))
amm_enrichments = amm_enrichments[amm_enrichments$gene_set_display == "oe1",]
amm_enrichments = amm_enrichments[c("phenotype_key", "Enrichment_mean", "Enrichment_se")]
colnames(amm_enrichments) <- c("phenotype_key", "amm_enrichment_mean", "amm_enrichment_se")

bhr = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST6"))
bhr = bhr[bhr$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",]
bhr = bhr[c("summary_statistic", "phenotype_key", "display_name", "bhr_enrichment_oe1", "bhr_enrichment_oe1_se")]
merged <- merge(amm_enrichments, bhr, by = "phenotype_key")
merged$bhr_enrichment_oe1_se_display = ifelse(merged$display_name %in% emphasis_traits, merged$bhr_enrichment_oe1_se, 0)
merged$amm_enrichment_se_display = ifelse(merged$display_name %in% emphasis_traits, merged$amm_enrichment_se, 0)

merged$display_only = merged$display_name
merged$display_only[!(merged$display_name %in% c(emphasis_traits, "FIS", "Townsend"))] <- ""

median(merged$bhr_enrichment_oe1/merged$amm_enrichment_mean)

f4b <- ggplot(merged, aes(x = amm_enrichment_mean, y = bhr_enrichment_oe1, xmin = amm_enrichment_mean - amm_enrichment_se_display, xmax = amm_enrichment_mean + amm_enrichment_se_display,
                   ymin = bhr_enrichment_oe1 - bhr_enrichment_oe1_se_display, ymax = bhr_enrichment_oe1 + bhr_enrichment_oe1_se_display,
                   color = display_name, label = display_only))+
  geom_errorbar(width = 0)+
  geom_errorbarh(height = 0)+
  geom_point(size = 5)+
  ggrepel::geom_text_repel(force_pull = .1, box.padding = 1.7, size = 5, family = "Arial")+
  labs(x = "Common-variant heritability enrichment\n in constrained genes (SE)", y = "Burden heritability enrichment\n in constrained genes (SE)", color='Gene set')+
  theme(text=element_text(family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15, face = "bold"),
        legend.title = element_text(size=15, face = "bold"),
        legend.position = "None",
        legend.direction = "horizontal",
        strip.text.x = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "white"))+
  geom_abline(slope = 1, linetype = "dashed")+
  scale_x_continuous(breaks = seq(1, 8, 1))+
  scale_y_continuous(breaks = seq(1, 8, 1))+
  scale_color_manual(values=colors_emphasis)+
  coord_cartesian(ylim = c(0, 8), xlim = c(0, 8))

###Figure 4C (common vs rare enrichment by constraint bins)
phenotype_file <- data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST4"))
phenotype_file = phenotype_file[phenotype_file$phenotype_core == 1,]
bhr = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST6"))
bhr = bhr[bhr$phenotype_key %in% phenotype_file$phenotype_key,]
bhr$frequency <- substrRight(bhr$summary_statistic, 6)
bhr$functional_group <- do.call(rbind.data.frame, stri_split_fixed(str = substring(bhr$summary_statistic, 42), pattern = "_", n = 2))[,1]
bhr = bhr[bhr$frequency == "group1" & bhr$functional_group == "pLoF",c("bhr_enrichment_oe1","bhr_enrichment_oe2","bhr_enrichment_oe3","bhr_enrichment_oe4", "bhr_enrichment_oe5")]
bhr = data.frame(melt(bhr))
bhr$set <- str_sub(bhr$variable, -3)
bhr = bhr[c("set", "value")]
colnames(bhr) <- c("gene_set", "enrichment")
bhr$model <- "Ultra-rare"

amm_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST13"))
amm_enrichments = amm_enrichments[amm_enrichments$gene_set_display %in% c("oe1", "oe2", "oe3", "oe4", "oe5"), c("gene_set_display", "Enrichment_mean")]
amm_enrichments$model <- "Common"
colnames(amm_enrichments) <- c("gene_set", "enrichment", "model")
merged <- rbind(bhr, amm_enrichments)
merged$gene_set = str_sub(merged$gene_set, -1)

f4c <- ggplot(merged, aes(x = factor(gene_set), y = enrichment, fill = factor(model)))+
  geom_boxplot(size = 1)+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  labs(x = "Constraint quintiles", y = "Heritability enrichment")+
  theme(text=element_text(family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=15,  color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15),
        legend.title = element_text(size=15),
        #legend.position = "None",
        legend.direction = "horizontal",
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(fill = "white"))+
  scale_fill_manual(values = c("white", "firebrick3"))+
  coord_cartesian(ylim = c(-.5, 8))+
  scale_y_continuous(breaks = seq(0, 8, 1))

###Figure 4d
phenotype_file <- data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST4"))
phenotype_file = phenotype_file[phenotype_file$phenotype_core == 1,]
bhr = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST6"))
bhr = bhr[bhr$phenotype_key %in% phenotype_file$phenotype_key,]
bhr = bhr[bhr$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1",]

bhr$display_only = bhr$display_name
bhr$display_only[!(bhr$display_name %in% c(emphasis_traits, "FIS", "Townsend"))] <- ""
bhr$sign = ifelse(bhr$mu_genome > 0 & bhr$display_only != "", "+", 
                  ifelse(bhr$mu_genome < 0 & bhr$display_only != "", "-", ""))
bhr$display_only_sign = ifelse(bhr$sign != "", paste0(bhr$display_name," (",bhr$sign,")"), "")

f4d <- ggplot(bhr, aes(x = bhr_enrichment_oe1, y = abs(mu_genome), color = display_only, label = display_only_sign))+
  geom_point(size = 3)+
  geom_text_repel(force_pull = .1, box.padding = 1.7, size = 5, family = "Arial")+
  theme_bhr()+
  labs(x = "Burden heritability enrichment\n in constrained genes (SE)", y = "Absolute mean effect size\n of ultra-rare pLoF alleles")+
  scale_color_manual(values=colors_emphasis)+
  scale_x_continuous(breaks = seq(1, 8, 1))
  
plot_grid(f4a, f4b, f4c, f4d, ncol = 2, labels = c('A', 'B', 'C', 'D'))

###Figure 5
plof_missense_compare_bhr_df = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST15", na = c("NA")))
plof_missense_compare_bhr_df_sig = plof_missense_compare_bhr_df[abs(plof_missense_compare_bhr_df$plof_h2/plof_missense_compare_bhr_df$plof_h2_se) > 1.96 & 
                                                                  abs(plof_missense_compare_bhr_df$missense_h2/plof_missense_compare_bhr_df$missense_h2_se) > 1.96,]

plof_missense_compare_bhr_df_sig$order = rank(plof_missense_compare_bhr_df_sig$rg)

Fig5A = ggplot(data = plof_missense_compare_bhr_df_sig,
       mapping = aes(x = factor(display_name, levels = display_name[order(rg)]), 
                     y = rg, 
                     ymin = rg - rg_se, 
                     ymax = rg + rg_se,
                     fill = display_name))+
  geom_col(color = "black")+
  geom_errorbar(width = 0, color = "black")+
  geom_hline(mapping = aes(yintercept = mean(rg)), linetype = "dashed", color = "black")+
  coord_flip()+
  theme_bhr()+
  labs(y = "pLoF-missense\nBurden genetic correlation (SE)")+
  scale_y_continuous(breaks = c(0,0.5,1))+
  scale_fill_manual(values = colors_emphasis)+
  theme(axis.title.x = element_text(vjust=10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = as.vector(colors_emphasis_Fig3B)[match(unique(plof_missense_compare_bhr_df_sig$display_name[order(plof_missense_compare_bhr_df_sig$rg)]),names(colors_emphasis_Fig3B))]))
  

bhr_rg = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST16"))
names_df = names(bhr_rg)
bhr_rg = as.data.frame(sapply(1:ncol(bhr_rg), function(x) unlist(bhr_rg[,x])))
bhr_rg[,c(3:6)] = sapply(c(3:6), function(x) as.numeric(bhr_rg[,x]))
names(bhr_rg)  = names_df
subset_traits = c("30780NA","30870NA","20002highcholesterol","30680NA","30600NA","20127NA","2090NA","3062NA","21001NA","20002osteoarthritis")
subset_traits_names = c("LDL","Triglycerides","High Cholesterol","Calcium","Albumin","Neuroticism","Depression","FVC","BMI","Osteoarthritis")
bhr_rg_subset = bhr_rg[bhr_rg$trait1 %in% subset_traits & bhr_rg$trait2 %in% subset_traits,]

cormat = matrix(data = NA, nrow = length(subset_traits), ncol = length(subset_traits))
zmat = matrix(data = NA, nrow = length(subset_traits), ncol = length(subset_traits))
ldscmat = matrix(data = NA, nrow = length(subset_traits), ncol = length(subset_traits))
zmat_ldsc = matrix(data = NA, nrow = length(subset_traits), ncol = length(subset_traits))

pairs = combn(10,2)
for (pair in 1:ncol(pairs)){
  traits = subset_traits[pairs[,pair]]
  index = which(bhr_rg_subset$trait1 %in% traits & bhr_rg_subset$trait2 %in% traits)
  cormat[pairs[1,pair],pairs[2,pair]] = bhr_rg_subset$bhr_rg[index]
  cormat[pairs[2,pair],pairs[1,pair]] = bhr_rg_subset$bhr_rg[index]
  
  ldscmat[pairs[1,pair],pairs[2,pair]] = bhr_rg_subset$ldsc_rg[index]
  ldscmat[pairs[2,pair],pairs[1,pair]] = bhr_rg_subset$ldsc_rg[index]
  
  zmat[pairs[1,pair],pairs[2,pair]] = bhr_rg_subset$bhr_rg[index]/bhr_rg_subset$bhr_rg_se[index]
  zmat[pairs[2,pair],pairs[1,pair]] = bhr_rg_subset$bhr_rg[index]/bhr_rg_subset$bhr_rg_se[index]
  
  zmat_ldsc[pairs[1,pair],pairs[2,pair]] = bhr_rg_subset$ldsc_rg[index]/bhr_rg_subset$ldsc_rg_se[index]
  zmat_ldsc[pairs[2,pair],pairs[1,pair]] = bhr_rg_subset$ldsc_rg[index]/bhr_rg_subset$ldsc_rg_se[index]
  
}
diag(cormat) <- 0


rownames(cormat) <- subset_traits_names
colnames(cormat) <- subset_traits_names

reorder_cormat <- function(cormat,source){
  # Use correlation between variables as distance
  dd <- as.dist((1-source)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
ldscmat <- reorder_cormat(ldscmat,cormat)
zmat  <- reorder_cormat(zmat,cormat)
zmat_ldsc  <- reorder_cormat(zmat_ldsc,cormat)

cormat <- reorder_cormat(cormat,cormat)
cormat[lower.tri(cormat)] <- ldscmat[lower.tri(ldscmat)]
zmat[lower.tri(zmat)] <- zmat_ldsc[lower.tri(zmat_ldsc)]

# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)
melted_cormat$z <- melt(zmat)$value
melted_cormat$sig = sapply(abs(melted_cormat$z), function(x) if(!is.na(x) & x> 1.96) {return("*")} else {return(NA)})

melted_cormat$value2 <- melted_cormat$value
melted_cormat$value2[melted_cormat$value2 > 1] <- 1
melted_cormat$value2[melted_cormat$value2 < -1] <- -1
melted_cormat$value2[melted_cormat$value2 == 0] <- NA

demingmodel = odregress(bhr_rg$ldsc_rg[!is.na(bhr_rg$ldsc_rg)],
                        bhr_rg$bhr_rg[!is.na(bhr_rg$ldsc_rg)])

Fig5C <- ggplot(data = bhr_rg,
                mapping = aes(x = ldsc_rg, y = bhr_rg))+
  geom_point(alpha = 0.25)+
  geom_abline(intercept = demingmodel$coeff[2], slope = demingmodel$coeff[1], linetype = "dashed")+
  geom_abline()+
  theme_bhr()+
  labs(x = "Common-variant genetic correlation", y = "Burden genetic correlation")+
  scale_y_continuous(breaks=seq(-1,1,0.5))+
  scale_x_continuous(breaks=seq(-1,1,0.5))+
  xlim(-1.05,1.05)+
  ylim(-1.05,1.05)+
  theme(axis.title.x = element_text(vjust=10))

# Create a ggheatmap
Fig5B <- ggplot(melted_cormat, aes(Var2, Var1, fill = value2))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name=expression(r[g])) +
   geom_text(aes(label = sig), vjust = 0.77)+
  coord_fixed()+
  theme(text=element_text(family="Arial"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 9, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=10),
        legend.title = element_text(size=15),
        legend.position = "right",
        legend.direction = "vertical",
        strip.background = element_rect(fill = "white"))

#plot_grid(Fig5A, Fig5B, Fig5C, ncol = 3, labels = c('A', 'B', 'C'))

Fig5 = (Fig5A | Fig5B |Fig5C)+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 20))

ggsave("5.pdf", Fig5,width = 15, height =5, dpi = 1200, device = cairo_pdf)

###Figure 6A
scz_bp_output = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=337549921", sheet = "ST17"))
bhr_h2_df = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST18"))
scz_bp_output$class = ifelse(scz_bp_output$class == "Syn", "Synonymous", scz_bp_output$class)
phenotype_data = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST4"))
bhr_h2_df = bhr_h2_df[bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key,]

Fig6A <- ggplot()+
  geom_pointrange(data = scz_bp_output,
                  mapping = aes(x = factor(class, levels = c("pLoF","Missense (MPC >2)", "Synonymous")), 
                                y = bhr_h2, 
                                ymin = bhr_h2 - bhr_h2_se, 
                                ymax = bhr_h2 + bhr_h2_se),
                  color = "white",
                  position = position_dodge2(width = 0.5),
                  size = 0.5)+
  geom_violin(data = bhr_h2_df[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1" &
                                 bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1],],
              mapping = aes(x = factor(rep("pLoF",nrow(bhr_h2_df[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1" &
                                                                   bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1],])),
                                       levels = c("pLoF","Missense (MPC >2)", "Synonymous")),
                            y = bhr_h2),
              color = "gray70",
              bw = 0.0025, 
              size = 2)+
  geom_violin(data = bhr_h2_df[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1" &
                                 bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1],],
              mapping = aes(x = factor(rep("Missense (MPC >2)",nrow(bhr_h2_df[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_missense-notbenign_low0_high1e-05_group1" &
                                                                                bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1],])),
                                       levels = c("pLoF","Missense (MPC >2)", "Synonymous")),
                            y = bhr_h2),
              color = "gray70",
              bw = 0.0025,
              size = 2)+
  geom_violin(data = bhr_h2_df[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0_high1e-05_group1" &
                                 bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1],],
              mapping = aes(x = factor(rep("Synonymous",nrow(bhr_h2_df[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_synonymous_low0_high1e-05_group1" &
                                                                                bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1],])),
                                       levels = c("pLoF","Missense (MPC >2)", "Synonymous")),
                            y = bhr_h2),
              color = "gray70",
              bw = 0.0025,
              size = 2)+
  geom_pointrange(data = scz_bp_output,
                  mapping = aes(x = factor(class, levels = c("pLoF","Missense (MPC >2)", "Synonymous")), 
                                y = bhr_h2, 
                                ymin = bhr_h2 - bhr_h2_se, 
                                ymax = bhr_h2 + bhr_h2_se, 
                                color = factor(dx, levels = c("SCZ","BP"))),
                  position = position_dodge2(width = 0.5), size = 0.75)+
  geom_text(data = scz_bp_output,
            mapping = aes(x = factor(class, levels = c("pLoF","Missense (MPC >2)", "Synonymous")), 
                          y = bhr_h2 + bhr_h2_se+0.0025, 
                          color = factor(dx, levels = c("SCZ","BP")),
                          label = c("*","*","","*","","")),
            position = position_dodge2(width = 0.5), size = 8)+
  geom_hline(yintercept = 0)+
  theme_bhr()+
  labs(x = "", y = "Burden heritability (SE)")+
  scale_color_manual(values = c("#AE8344","#068E90"))+
  guides(color = "none")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))+
  scale_y_continuous(breaks = c(0, 0.01, .02), labels = function(x) paste0(x*100, "%"))+
  annotate("text", x = "Synonymous", y = 0.025, size = 6, color = "#AE8344", label = "Schizophrenia")+
  annotate("text", x = "Synonymous", y = 0.0225, size = 6, color = "#068E90", label = "Bipolar disorder")+
  annotate("text", x = "Synonymous", y = 0.02, size = 6, color = "gray70", label = "UKB (22 traits)")

ggsave("Fig6_top.pdf", Fig6A, width = 10, height = 5, dpi = 1200, device = cairo_pdf)

###Figure 6B
amm_enrichments = data.frame(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1AwvDRKEUJ6EtUbvnvSV34JkPV0Y7ozJHw-vDJOW2PpY/edit#gid=881218879", sheet = "ST13"))
amm_enrichments = amm_enrichments[amm_enrichments$gene_set_display == "oe1",]
amm_enrichments = amm_enrichments[c("phenotype_key", "Enrichment_mean", "Enrichment_se")]
colnames(amm_enrichments) <- c("phenotype_key", "amm_enrichment_mean", "amm_enrichment_se")

enrichmentcompare_df <- data.frame(enrich = c(scz_bp_output$enrichment_constrained[scz_bp_output$dx == "SCZ" & scz_bp_output$class == "pLoF"],
                                              scz_bp_output$enrichment_constrained[scz_bp_output$dx == "BP" & scz_bp_output$class == "pLoF"],
                                              3.01881758,
                                              2.735787707,
                                              mean(bhr_h2_df$bhr_enrichment_oe1[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1" & bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1]]),
                                              mean(amm_enrichments$amm_enrichment_mean)),
                                   se = c(scz_bp_output$enrichment_constrained_se[scz_bp_output$dx == "SCZ" & scz_bp_output$class == "pLoF"],
                                          scz_bp_output$enrichment_constrained_se[scz_bp_output$dx == "BP" & scz_bp_output$class == "pLoF"],
                                          0.282939424,
                                          0.329771062,
                                          sd(bhr_h2_df$bhr_enrichment_oe1[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1" & bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1]])/sqrt(length(bhr_h2_df$bhr_enrichment_oe1[bhr_h2_df$summary_statistic == "bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1" & bhr_h2_df$phenotype_key %in% phenotype_data$phenotype_key[phenotype_data$phenotype_core ==1]])),
                                          sd(amm_enrichments$amm_enrichment_mean)/sqrt(length(amm_enrichments$amm_enrichment_mean))),
                                   dx = c("Schizophrenia","Bipolar Disorder","Schizophrenia","Bipolar Disorder","UKB (22 traits)","UKB (22 traits)"),
                                   type = c("URV PTV (BHR)","URV PTV (BHR)","Common (AMM)","Common (AMM)", "URV PTV (BHR)","Common (AMM)"),
                                   label = c("Ultra-rare\npLoF","Ultra-rare\npLoF","Common","Common","Ultra-rare\npLoF","Common"),
                                   unambiguous = c("SCZURV","SCZCommon","BPURV","BPCommon","UKBURV","UKBCommon"))

Fig6B = ggplot(data = enrichmentcompare_df, mapping = aes(x = factor(dx,levels = c("Schizophrenia","Bipolar Disorder","UKB (22 traits)")), 
                                                          group = label, 
                                                          label = label,
                                                          y = enrich, 
                                                          ymin = enrich - se, 
                                                          ymax = enrich + se, 
                                                          color = unambiguous,
                                                          fill = unambiguous))+
  geom_col(position = position_dodge(), color = "black")+ 
  geom_text(position = position_dodge(width = 0.9), aes(x=factor(dx, levels = c("Schizophrenia","Bipolar Disorder","UKB (22 traits)")), y=1.25), color = "black")+
  scale_fill_manual(values = c("#4CB7CB","#E3C38A","#4CB7CB","#E3C38A","gray70","gray70"))+
  geom_linerange(position = position_dodge(0.9),  color  = "black")+
  theme_bhr()+
  guides(color = "none", fill = "none")+
  labs(x = "", y = "Heritability enrichment\nin constrained genes (SE)")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))

Fig6 = (Fig6A | Fig6B) + plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 20))

ggsave("Fig6.pdf", Fig6, width = 15, height = 5, dpi = 1200, device = cairo_pdf)
