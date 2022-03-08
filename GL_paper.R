
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(data.table)
library(DescTools)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(devtools)
library(tidyverse)
library(microbiome)
library(metacoder)


#### Alpha diversity
plot_richness <- plot_richness(water_sediment_r, measures=c("Observed"), shape="fraction", color="mountain_range") +
  theme(axis.text.x=element_text(angle=70, hjust=1)) +
  xlab("") + ylab("Observed") + ggtitle("ASV_richness")

alphadt_watersed<- data.table(plot_richness$data)

library(cowplot)
ncolors <- length(levels(alphadt_watersed$fraction))
primary <- 'fraction'
primary_color_list <- brewer.pal(ncolors, 'Dark2') 

asvrichnes_plot <- ggplot(alphadt_watersed, aes(x = fraction, y = value, fill=fraction)) +
  geom_boxplot() + 
  ylab('ASVRichness') +
  xlab(primary) +
  scale_fill_manual(name = primary, values = primary_color_list) +
  labs(fill = primary) +
  facet_grid(~ mountain_range)+
  theme_cowplot() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank()) 

obs.stat <- aov((alphadt_watersed$value) ~ fraction * mountain_range, data=alphadt_watersed)

residual_obs_stat <- resid(obs.stat); F1_obs.stat <- fitted(obs.stat)
plot(F1_obs_stat, residual_obs_stat, xlab = "Fitted Values", ylab = "Normalized residuals"); abline(h = 0, lty = 2) 
summary(obs.stat)
anova(obs.stat)
hist(residual_obs.stat) 
qqnorm(residual_obs_stat)
shapiro.test(residual_obs_stat)
LeveneTest(residual_obs_stat ~ fraction*mountain_range, data=alphadt_watersed)

#### Community composition 

##NMDS

watersed.ord <- ordinate(water_sediment_r, "NMDS", "bray", k=2)
watersed.ord$stress

p1 = plot_ordination(water_sediment_r, watersed.ord, type="Samples", color="gl_name") 
p1 = p1 + geom_point(size=4, alpha=4) 
p1 + geom_text(aes(label=code),hjust=1, vjust=1)
p1 + scale_color_manual(values=primary_color_list)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))         
p1 + geom_text(aes(label=code),hjust=1, vjust=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))         

### manyglm

otu_watersed <- t(otu_table(water_sediment_r, taxa_are_rows=T))
ab <- mvabund(otu_watersed)
mod_bino<- manyglm(ab ~ mountain_range*fraction,
                    data = alphadt_watersed, family = 'negative binomial')
pairwise<- anova(mod_bino, pairwise.comp=~mountain_range*fraction, nBoot = 99)

#### Compositional variability

## betadisper

bdisp <- vegdist(vegan_matrix_sedwater, "bray")
bdisp_mr<- betadisper((bdisp), alphadt_watersed$mountain_range, type=c("centroid"))
bdisp_mr
aov.bdisp_mr <-anova(bdisp_mr)
aov.bdisp_mr
permutest(bdisp_mr)

box_mr <- cbind(alphadt_watersed, dist=bdisp.mr$distances)
box_mr

ggplot(box_mr, aes(x=mountain_range, y=dist))+
  geom_boxplot()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bdisp_mr<- betadisper((bdisp), alphadt_watersed$fraction, type=c("centroid"))
bdisp_fr
aov.bdisp_fr <-anova(bdisp_fr)
aov.bdisp_fr
permutest(bdisp_fr)

box_fr <- cbind(alphadt_watersed, dist=bdisp.fr$distances)
box_fr

ggplot(box_fr, aes(x=fraction, y=dist))+
  geom_boxplot()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tukey_bdisp_mr<- TukeyHSD(bdisp.mr)
tukey_bdisp_fr<- TukeyHSD(bdisp.fr)


#### Upsetplot
library(UpSetR)

fraction <- c('sedup.x', 'waternf.x','waterwf.x')
alphadt_watersed$fraction <- factor(alphadt_watersed$fraction, levels = fraction)

ps2.venn <- merge_samples(water_sediment_r, 'fraction', fun = sum)
ncolors <- length(levels(alphadt_watersed$fraction))
primary_color_list <- brewer.pal(ncolors, 'Dark2') 

venn_obj <- as.data.frame(t(otu_table(ps2.venn)))

## transform into binary
venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, TRUE, FALSE),
                          USE.NAMES = T)
## transform relative abundance
venn_obj.ab <- sapply(venn_obj, function(x) ifelse(x > 0, x/17309, 0),
                      USE.NAMES = T)
## give names to rows
rownames(venn_obj.ab) <- rownames(venn_obj)
rownames(venn_obj.binary) <- rownames(venn_obj)

## always transform as dataframe
venn_obj_binary_df <- as.data.frame(venn_obj.binary)
venn_obj_ab_df<- as.data.frame(venn_obj.ab)

## compute Sum of the taxa
venn_obj_ab_df$Ab_Sum <- rowSums(venn_obj_ab_df)

## merging tables
venn_obj_merge <- merge(venn_obj_binary_df, venn_obj_ab_df, by=0)
venn_obj_merge_df <- as.data.frame(venn_obj_merge)

## compute taxonomy
venn_ra_transform  = transform_sample_counts(ps2.venn, function(x) x / sum(x))

data_glom_venn<- psmelt(venn_ra_transform) # create dataframe from phyloseq object
data_glom_venn$Family <- as.character(data_glom_venn$Family) #convert to character

order_abundance <- data_glom_venn

## subset three vectors
order_abundance_short<-cbind(order_abundance$OTU,order_abundance$Phylum, order_abundance$Family)
order_ab_short_df <- as.data.frame(order_abundance_short)

## rename columns
colnames(order_ab_short_df) <- c("OTU","Phylum","Family")

## Only pick unique rows
order_ab_short_df <- order_ab_short_df [!duplicated(order_ab_short_df ),]

## merge dataframes
merge_abundance_taxo <- merge(venn_obj_merge_df, order_ab_short_df, by.x="Row.names", by.y="OTU")

## Top families
taxglom_ra_transform  = transform_sample_counts(taxglom_sedwater_family_saved, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(taxglom_ra_transform), TRUE)[1:14])
top15_NOMIS_f <- prune_species(TopASV_f, taxglom_ra_transform )
top15_NOMIS_f <- prune_taxa(taxa_sums(top15_NOMIS_f)>0, top15_NOMIS_f)
top_family<-as.data.frame(tax_table(top15_NOMIS_f))

## Other for families that are less abundant
merge_abundance_taxo$Family[!(merge_abundance_taxo$Family %in% top_family$Family)] <- "Other"
merge_abundance_taxo$Family[(merge_abundance_taxo$Family == "")] <- "Other"

## with Rbrewer
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

col=colorRampPalette(brewer.pal(11,"Spectral"))(n)

merge_abundance_taxo_df <- as.data.frame(merge_abundance_taxo)

upset(
  merge_abundance_taxo_df, fraction2
)

fraction2 <- colnames(merge_abundance_taxo_df)
fraction2 <- fraction2[2:4]

## Plot the upset
coucou<-upset(
  merge_abundance_taxo_df,
  fraction,
  annotations = list(
    'Abundance'=(
      ggplot(mapping=aes(y=Ab_Sum,fill=Family))
      + geom_bar(stat="identity", position="fill")
      + scale_y_continuous(labels=scales::percent_format())
      +scale_fill_manual(values=c("#7FC97F","#1B9E77","#D95F02","#7570B3","#E7298A", 
                                  "#66A61E", "#E6AB02","#A6761D","#326CF0","#A6CEE3","#1F78B4","#E2DF8A","#FACD87","#444f92","#559ffd"))
    )
  ),
  width_ratio=0.1)


#### Differential abundance analyses using Metacoder

metadata_sedwater <- plot_richness(taxglom_sedwater_family, measures=c("Observed"), shape="fraction", color="mountain_range") +
  theme(axis.text.x=element_text(angle=70, hjust=1)) +
  xlab("") + ylab("Observed") + ggtitle("ASV richness")

alphadt_sedwater_tax<- data.table(metadata_sedwater$data)

obj <- parse_tax_data(sedwater_tib,
                      class_cols = "lineage", 
                      class_sep = ";")

##Removing low-abundance counts
obj$data$tax_data <- zero_low_counts(obj, dataset = "tax_data", min_count = 5)
no_reads <- rowSums(obj$data$tax_data[, alphadt_sedwater_tax$Sample]) == 0
sum(no_reads)

obj <- filter_obs(obj, target = "tax_data", ! no_reads, drop_taxa = TRUE)
print(obj)

## Taking the proportion
obj$data$tax_data <- calc_obs_props(obj, "tax_data")


## Getting per-taxon information
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = alphadt_sedweater_tax$Sample)

obj$data$tax_abund <- obj$data$tax_abund %>%
  mutate(taxon = taxon_names(obj)) %>%
  select(taxon_id, taxon, everything())

obj

## Compare different treatments
obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund",
                                      cols = alphadt_sedwater_tax$Sample, # What columns of sample data to use
                                      groups = alphadt_sedwater_tax$fraction) # What category each sample is assigned to
range(obj$data$diff_table$wilcox_p_value, finite = TRUE) 

##Include Holm correction
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "holm")

range(obj$data$diff_table$wilcox_p_value, finite = TRUE) 

obj$data$diff_table %>% pull(log2_median_ratio)

obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

diftable<-as.data.frame(obj$data$diff_table)
taxabund <-as.data.frame(obj$data$tax_abund)

dont_print <- c("f__unknown", "f__unidentified")

ncolors <- length(levels(alphadt_sedwater_tax$fraction))
primary <- 'fraction'
primary_color_list <- brewer.pal(ncolors, 'Dark2') 

set.seed(1)
heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs, 
                 node_label = taxon_names,
                 node_color = log2_median_ratio, 
                 node_color_range = c("#D95F02","gray","#1B9E77"), 
                 node_color_trans = "linear", 
                 node_color_interval = c(-3, 3), 
                 edge_color_interval = c(-3, 3), 
                 node_size_axis_label = "Number of ASVs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 node_label_max = 500,
                 layout = "davidson-harel", 
                 initial_layout = "reingold-tilford",
                 output_file = "differential_heat_tree_family.pdf") 












