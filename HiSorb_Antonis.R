# Load packages
library(tidyverse)
library(viridis)
library(ggpubr)
library(mixOmics)
library(gridExtra)
library(patchwork)
library(readxl)
library(tidymodels)

# Set working directory
setwd("~/OneDrive - University of California, San Diego Health/Projects/HiSorb_Antonis")

# Get data
polar_gavin <- read_excel("Polar - 28ADvs33C Gavin3.xlsx", col_names = TRUE) %>% as.data.frame() %>% 
  dplyr::filter(!(ExperimentNo %in% c("P17", "P66", "C47", "C52")))

nopolar_gavin <- read_excel("Non-polar - 28ADvs33C Gavin3.xlsx", col_names = TRUE) %>% as.data.frame() %>%
  dplyr::filter(ExperimentNo != "P17")

# Extract metabolites
polar_gavin_met <- polar_gavin %>% dplyr::select(`1-Propene‚ 2-methoxy-`:last_col())
nopolar_gavin_met <- nopolar_gavin %>% dplyr::select(`Diisobutyl cellosolve`:last_col())

# Check metabolite recovery
polar_gavin_met_tot <- rowSums(polar_gavin_met) %>% as.data.frame() %>% rename("TotalArea" = ".")
nopolar_gavin_met_tot <- rowSums(nopolar_gavin_met) %>% as.data.frame() %>% rename("TotalArea" = ".")

ggdensity(polar_gavin_met_tot, "TotalArea")
ggdensity(nopolar_gavin_met_tot, "TotalArea")

# Total area normalisation
polar_gavin_met_ra <- polar_gavin_met/polar_gavin_met_tot$TotalArea
nopolar_gavin_met_ra <- nopolar_gavin_met/nopolar_gavin_met_tot$TotalArea

####################
#### ANALYSIS ######
####################

# PCA - POLAR
PCA_polar_gavin_met <- mixOmics::pca(polar_gavin_met_ra, ncomp = 2, center = TRUE, scale = TRUE)
PCA_polar_gavin_met_scores <- data.frame(PCA_polar_gavin_met$variates$X, dplyr::select(polar_gavin, ExperimentNo:CancerSubtype))

PCA_plot <- PCA_polar_gavin_met_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "CancerSubtype", alpha = 0.6,
            xlab = paste("PC1 (", round(PCA_polar_gavin_met$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_polar_gavin_met$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_polar_gavin_met_scores %>% group_by(CancerSubtype) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (CancerSubtype)), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PCA - NON POLAR
PCA_nopolar_gavin_met <- mixOmics::pca(nopolar_gavin_met_ra, ncomp = 2, center = TRUE, scale = TRUE)
PCA_nopolar_gavin_met_scores <- data.frame(PCA_nopolar_gavin_met$variates$X, dplyr::select(nopolar_gavin, ExperimentNo:CancerSubtype))

PCA_plot <- PCA_nopolar_gavin_met_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "CancerSubtype", alpha = 0.6,
            xlab = paste("PC1 (", round(PCA_nopolar_gavin_met$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_nopolar_gavin_met$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_nopolar_gavin_met_scores %>% group_by(CancerSubtype) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (CancerSubtype)), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Combine polar and non polar data in a single dataset.
common_samples_gavin <- polar_gavin %>% dplyr::filter(TBID %in% nopolar_gavin$TBID) %>% dplyr::select(TBID)

polar_gavin_tounite <- polar_gavin %>% dplyr::filter(TBID %in% common_samples_gavin$TBID) %>% arrange(TBID)
nopolar_gavin_tounite <- nopolar_gavin %>% dplyr::filter(TBID %in% common_samples_gavin$TBID) %>% arrange(TBID)

meta_gavin_unite <- polar_gavin_tounite %>% dplyr::select(ExperimentNo:CancerSubtype)

polar_gavin_met_unite <- polar_gavin_tounite %>% dplyr::select(`1-Propene‚ 2-methoxy-`:last_col())
nopolar_gavin_met_unite <- nopolar_gavin_tounite %>% dplyr::select(`Diisobutyl cellosolve`:last_col())

overlap <- inner_join(as.data.frame(colnames(polar_gavin_met_unite)), as.data.frame(colnames(nopolar_gavin_met_unite)), 
                      by = c("colnames(polar_gavin_met_unite)" = "colnames(nopolar_gavin_met_unite)"))

# 45 metabolites are overlapping --> keep everything, change variable names
colnames(nopolar_gavin_met_unite) <- paste0("nopolar_", colnames(nopolar_gavin_met_unite))
gavin_unite <- cbind(polar_gavin_met_unite, nopolar_gavin_met_unite)

# Total area normalisation
gavin_met_tot <- rowSums(gavin_unite) %>% as.data.frame() %>% rename("TotalArea" = ".")
gavin_met_ra <- gavin_unite/gavin_met_tot$TotalArea

# PCA
PCA_gavin_met <- mixOmics::pca(gavin_met_ra, ncomp = 2, center = TRUE, scale = TRUE)
PCA_gavin_met_scores <- data.frame(PCA_gavin_met$variates$X, dplyr::select(polar_gavin_tounite, ExperimentNo:CancerSubtype))

PCA_plot <- PCA_gavin_met_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "CancerSubtype", alpha = 0.6,
            xlab = paste("PC1 (", round(PCA_gavin_met$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_gavin_met$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_gavin_met_scores %>% group_by(CancerSubtype) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (CancerSubtype)), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_gavin <- mixOmics::plsda(gavin_met_ra, meta_gavin_unite$CancerSubtype, ncomp = 3, scale = TRUE)
PLSDA_gavin_scores <- data.frame(PLSDA_gavin$variates$X, meta_gavin_unite)

PLSDA_gavin_plot <- PLSDA_gavin_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "CancerSubtype", alpha = 0.6,
            xlab = paste("Component 1 (", round(PLSDA_gavin$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_gavin$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), title = "PLS-DA - Adenocarcinoma vs Control") +
  geom_point(data = PLSDA_gavin_scores %>% group_by(CancerSubtype) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = CancerSubtype), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

#ggsave(plot = PLSDA_gavin_plot, filename = "PLSDA_plot.svg", device = "svg", dpi = "retina", width = 4, height = 4)

Loadings_gavin <- plotLoadings(PLSDA_gavin, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

perf_plsda_gavin <- perf(PLSDA_gavin, validation = "loo", progressBar = FALSE, auc = TRUE) 
plot(perf_plsda_gavin, legend = FALSE)

VIPs_gavin <- as.data.frame(mixOmics::vip(PLSDA_gavin))
VIPs_gavin_filter <- dplyr::filter(VIPs_gavin, VIPs_gavin$comp1 > 1.5)
VIPs_gavin_filter$ID <- rownames(VIPs_gavin_filter)
VIPs_gavin_select <- VIPs_gavin_filter %>% dplyr::select(ID, comp1)
VIPs_gavin_Load <- VIPs_gavin_select %>% left_join(Loadings_gavin, by = c("ID" = "rowname"))

#write.csv(VIPs_gavin_Load, "VIPs_PLSDA.csv", fileEncoding="UTF-16LE")

# Select metabolites of interest for final model
gavin_select <- gavin_met_ra %>% dplyr::select("2-Pentanone", "Hexanal", "3-Hexanone","1‚3‚5-Cycloheptatriene‚ 3‚7‚7-trimethyl-")

# Metabolites ID have been double checked with spectra and one has been reclassified as p-Cymene
# Names are "2-Pentanone", "Hexanal", "3-Hexanone", "p-Cymene" but having some problem with new R and getting error when 
# generating plots. I will remove the number from the names
colnames(gavin_select) <- c("Pentanone", "Hexanal", "Hexanone", "Cymene")

# Save table to be exported to MetaboAnalyst too
gavin_select$Class <- meta_gavin_unite$CancerSubtype
#write_csv(gavin_select, "Metabolites_interest.csv")

# Check classification performance of model generated using only the 4 metabolites of interest
set.seed(123)
# Split data in training and test
data_split <- initial_split(gavin_select, strata = "Class")
train_data <- training(data_split)
test_data <- testing(data_split)

train_data_metabolite <- train_data %>% dplyr::select(-Class)
train_data_metadata <- train_data %>% dplyr::select(Class)

test_data_metabolite <- test_data %>% dplyr::select(-Class) %>% as.matrix()
test_data_metadata <- test_data %>% dplyr::select(Class)
test_data_metadata$Class <- as.factor(test_data_metadata$Class)

# PLSDA - TRAIN
PLSDA_gavin <- mixOmics::plsda(train_data_metabolite, train_data_metadata$Class, ncomp = 2, scale = TRUE)
PLSDA_gavin_scores <- data.frame(PLSDA_gavin$variates$X, train_data_metadata)

PLSDA_gavin_plot <- PLSDA_gavin_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Class", alpha = 0.6,
            xlab = paste("Component 1 (", round(PLSDA_gavin$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_gavin$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_gavin_scores %>% group_by(Class) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Class), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_gavin <- plotLoadings(PLSDA_gavin, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

perf_plsda_gavin <- perf(PLSDA_gavin, validation = "loo", progressBar = FALSE, auc = TRUE)
plot(perf_plsda_gavin, legend = FALSE)

# Generate ROC using TEST data
roc_plot <- auroc(PLSDA_gavin, roc.comp = 1, newdata = test_data_metabolite, outcome.test = test_data_metadata$Class,
                  title = "ROC Curve - Control vs Adenocarcinoma")
roc_plot$graph.Comp1 + theme(legend.position="none") # save as pdf

# Check influence of antibiotics and chemotherapy on metabolites of interest
gavin_check <- gavin_select %>% dplyr::mutate(Chemo = meta_gavin_unite$ChemoTherapy) %>%
  dplyr::mutate(Antibiotics = meta_gavin_unite$Antibiotics) %>%
  dplyr::mutate(Diabetes = meta_gavin_unite$Diabetes) 

gavin_check$Chemo <- gsub("0", "No", gavin_check$Chemo)
gavin_check$Chemo <- gsub("1", "Yes", gavin_check$Chemo)
gavin_check$Antibiotics <- gsub("0", "No", gavin_check$Antibiotics)
gavin_check$Antibiotics <- gsub("1", "Yes", gavin_check$Antibiotics)
gavin_check$Diabetes <- gsub("0", "No", gavin_check$Diabetes)
gavin_check$Diabetes <- gsub("1", "Yes", gavin_check$Diabetes)

gavin_chemo <- list()

for (i in c("Pentanone", "Hexanal", "Hexanone", "Cymene")) {
  
  plot_box <- gavin_check %>% 
    dplyr::filter(Class == "AD") %>%
    ggboxplot(x = "Chemo", y = i, title = i, add = "jitter", legend = "none",
              add.params = list(color = "Chemo", alpha = 0.8, size = 1),
              ylab = "Relative Abundance") + stat_compare_means(label = "p.format") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 6),
          plot.title = element_text(size = 7, hjust = 0.5))
  
  gavin_chemo[[i]] <- plot_box  # add each plot into plot list
  
}

boxplots_chemo <- wrap_plots(gavin_chemo, nrow = 1)

#ggsave(boxplots_chemo, filename = "Boxplots_Chemo.svg", device = "svg", dpi = "retina", width = 6, height = 3, limitsize = FALSE)

gavin_ab <- list()

for (i in c("Pentanone", "Hexanal", "Hexanone", "Cymene")) {
  
  plot_box <- gavin_check %>%
    ggboxplot(x = "Antibiotics", y = i, title = i, add = "jitter", legend = "none",
              add.params = list(color = "Antibiotics", alpha = 0.8, size = 1),
              ylab = "Relative Abundance") + stat_compare_means(label = "p.format") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 6),
          plot.title = element_text(size = 7, hjust = 0.5))
  
  gavin_ab[[i]] <- plot_box  # add each plot into plot list
  
}

boxplots_antibiotics <- wrap_plots(gavin_ab, nrow = 1)

#ggsave(boxplots_antibiotics, filename = "Boxplots_Antibiotics.svg", device = "svg", dpi = "retina", width = 6, height = 3, limitsize = FALSE)

gavin_diabetes <- list()

for (i in c("Pentanone", "Hexanal", "Hexanone", "Cymene")) {
  
  plot_box <- gavin_check %>% dplyr::filter(Class == "AD") %>%
    ggboxplot(x = "Diabetes", y = i, title = i, add = "jitter", legend = "none",
              add.params = list(color = "Diabetes", alpha = 0.8, size = 1),
              ylab = "Relative Abundance") + stat_compare_means(label = "p.format") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 6),
          plot.title = element_text(size = 7, hjust = 0.5))
  
  gavin_diabetes[[i]] <- plot_box  # add each plot into plot list
  
}

boxplots_diabetes <- wrap_plots(gavin_diabetes, nrow = 1)

#ggsave(boxplots_diabetes, filename = "Boxplots_Diabetes.svg", device = "svg", dpi = "retina", width = 6, height = 3, limitsize = FALSE)



#########################
### CHECK OSMOLARITY  ###
#########################

# Osmolarity normalisation
polar_gavin_met_os <- polar_gavin_met/polar_gavin$Osmolality
nopolar_gavin_met_os <- nopolar_gavin_met/nopolar_gavin$Osmolality

# PCA - POLAR
PCA_polar_gavin_met <- mixOmics::pca(polar_gavin_met_os, ncomp = 2, center = TRUE, scale = TRUE)
PCA_polar_gavin_met_scores <- data.frame(PCA_polar_gavin_met$variates$X, dplyr::select(polar_gavin, ExperimentNo:CancerSubtype))

PCA_polar_gavin_os_plot <- PCA_polar_gavin_met_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "CancerSubtype", alpha = 0.6,
            xlab = paste("PC1 (", round(PCA_polar_gavin_met$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_polar_gavin_met$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_polar_gavin_met_scores %>% group_by(CancerSubtype) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (CancerSubtype)), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PCA - NON POLAR 
PCA_nopolar_gavin_met <- mixOmics::pca(nopolar_gavin_met_os, ncomp = 2, center = TRUE, scale = TRUE)
PCA_nopolar_gavin_met_scores <- data.frame(PCA_nopolar_gavin_met$variates$X, dplyr::select(nopolar_gavin, ExperimentNo:CancerSubtype))

PCA_nopolar_gavin_os_plot <- PCA_nopolar_gavin_met_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "CancerSubtype", alpha = 0.6,
            xlab = paste("PC1 (", round(PCA_nopolar_gavin_met$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_nopolar_gavin_met$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_nopolar_gavin_met_scores %>% group_by(CancerSubtype) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (CancerSubtype)), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Osmolarity normalization does not appear to work well

# Plot total area and osmolarity to check correlation
polar_gavin_met_os <- cbind(polar_gavin_met_tot, polar_gavin$Osmolality)
nopolar_gavin_met_os <- cbind(nopolar_gavin_met_tot, nopolar_gavin$Osmolality)

a <- polar_gavin_met_os %>% 
  ggscatter("TotalArea", "polar_gavin$Osmolality", add = "reg.line", conf.int = TRUE, color = "dodgerblue4",
            xlab = "Total Peak Area", ylab = "Osmolarity", add.params = list(color = "black"), title = "Gavin - Polar") + 
  stat_cor(method = "pearson", label.x = 1400, label.y = 1100)

b <- nopolar_gavin_met_os %>% 
  ggscatter("TotalArea", "nopolar_gavin$Osmolality", add = "reg.line", conf.int = TRUE, color = "firebrick4",
            xlab = "Total Peak Area", ylab = "Osmolarity", add.params = list(color = "black"), title = "Gavin - Non Polar") + 
  stat_cor(method = "pearson", label.x = 600, label.y = 1100)

# No correlation
#ggsave(ggarrange(a, b), filename = "Osmolarity.svg", device = "svg", dpi = "retina", width = 8, height = 5)
