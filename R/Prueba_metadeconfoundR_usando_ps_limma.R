#Data
```{r MetadeconfoundR heatmaps_limma}
library(devtools)
install_github("TillBirkner/metadeconfoundR")

ps_limma <- readRDS("/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/outputs/ps_limma.Rds")
library(readxl)
metadata_hiper <- data.frame(sample_data(ps_limma))
View(metadata_hiper)
```

```{r Changing names ASVs}
library(tidyr)
library(tibble)
library(purrr)
require(phyloseq)
require(tidyverse)
require(magrittr)


get_latest_annotation <- function(phyloseq_obj) {
  tax_table <- phyloseq_obj@tax_table %>%
    as.data.frame() %>%
    rownames_to_column('ASV') %>%
    as_tibble() %>%
    tidyr::gather('Rank', 'Value', -ASV) %>%
    nest(data = -ASV) %>%
    mutate(TaxaID = map_chr(data, function(x) {
      x %>%
        dplyr::filter(!is.na(Value)) %>%
        magrittr::use_series(Value) %>%
        tail(1)
    })) %>%
    mutate(TaxaUp = map_chr(data, function(x) {
      x %>%
        dplyr::filter(!is.na(Value)) %>%
        magrittr::use_series(Value) %>%
        tail(2) %>% head(1)
    })) %>%
    mutate(TaxaID = paste(ASV, TaxaID)) %>%
    dplyr::select(-c(data, TaxaUp))
  
  return(tax_table)
}

##Usage latest_annotations <- get_latest_annotation(phyloseq_obj)
latest_annotations <- get_latest_annotation(ps_limma)
latest_annotations<-as.data.frame(latest_annotations)
short_names<-as.data.frame(otu_table(ps_limma))
rownames(short_names)<-latest_annotations$TaxaID
```

```{r Changing names ASVs long}
# generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(ps_limma))
                             [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 
                             sep = "__"))  # to distinguish from "_" within tax ranks

# turn the otu_table into a data.frame
otu_export <- as.data.frame(t(otu_table(ps_limma)))
tmp <- names(otu_export)

# paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

# overwrite old names
names(otu_export) <- names(tmp)
otu_export<-t(otu_export)

head(otu_export)[5]
```

##MetadeconfoundR
```{r MetadeconfoundR heatmaps_limma}
library(metadeconfoundR)
library(dplyr)
ps_limma.rel<-ps_limma
latest_annotations <- get_latest_annotation(ps_limma.rel)
latest_annotations<-as.data.frame(latest_annotations)
short_names<-as.data.frame(otu_table(ps_limma.rel))
rownames(short_names)<-latest_annotations$TaxaID
Meta_Species <- t(short_names) 

metadata_hiper_R <- subset(metadata_hiper, select = -c(SampleID, Sample_type))
metadata_hiper_R <- metadata_hiper_R %>% mutate(Group = ifelse(Group=="Hiperammonemic",1,0))
metadata_hiper_R <- as.data.frame(metadata_hiper_R)
rownames(metadata_hiper_R) <- (metadata_hiper$SampleID)

# check correct ordering
all(rownames(metadata_hiper_R) == rownames(Meta_Species))
## [1] TRUE
all(order(rownames(metadata_hiper_R)) == order(rownames(Meta_Species)))
## [1] TRUE
Meta_Species<- as.data.frame(Meta_Species)
Output1 <- MetaDeconfound(featureMat = as.data.frame(Meta_Species),
                          metaMat = as.data.frame(metadata_hiper_R), nnodes = 14)

Output2_batch <- MetaDeconfound(featureMat = as.data.frame(Meta_Species),
                                metaMat = as.data.frame(metadata_hiper_R), nnodes = 14, randomVar = list("+ (1|Batch)",
                                                                                                         c("Batch")))
View(Output1)
View(Output2_batch)

left <- BuildHeatmap(Output1)
right <- BuildHeatmap(Output2_batch)
#cannot plot it because there is only one variale significant (group)
saveRDS(Output1, "/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/output/Output1_metadec.Rds")
saveRDS(Output2_batch, "/fast/AG_Forslund/Lola/Datos_CIPF_2023/Rats_HE_HiperAmmon_model/output/Output2_batch_all.Rds")

```

```{r plotting}
raw_p<- Output2_batch[1]
corr_p<- Output2_batch[2]
effect_size <- Output2_batch[3]
status<- Output2_batch[4]
raw_p_df<- data.frame(raw_p$Ps)
raw_p_df<- raw_p_df %>% subset(select = c(Group)) %>%
  rownames_to_column()%>%
  mutate(p_Group=Group)%>%
  select(!c(Group))

corr_p_df<- data.frame(corr_p$Qs)
corr_p_df<- corr_p_df %>%subset(select = c(Group)) %>%
  rownames_to_column()%>%
  mutate(corr_Group=Group)%>%
  select(!c(Group))


effect_size_df<- data.frame(effect_size$Ds)
effect_size_df<- effect_size_df %>% subset(select = c(Group)) %>%
  rownames_to_column()%>%
  mutate(effect_Group=Group)%>%
  select(!c(Group))


status_df<- data.frame(status$status)
status_df<- status_df %>% subset(select = c(Group)) %>%
  rownames_to_column()%>%
  mutate(status_Group=Group)%>%
  select(!c(Group))

#create two-column-dataframe containing corresponding "human-readable" names to the "machine-readable" feature names used as row.names in metaDeconfOutput.  
taxtable <- latest_annotations
taxtable$rowname <- latest_annotations$ASV
taxtable<- cbind(rowname=taxtable$rowname,subset(taxtable,select = -c(rowname)))

effect_table <- raw_p_df%>%
  full_join(corr_p_df, by="rowname")%>%
  full_join(effect_size_df, by="rowname")%>%
  full_join(status_df, by="rowname")%>%
  full_join(taxtable, by="rowname")

# remove the entries which have NS and AD in status
effect_table_sig <- effect_table%>%
  filter(status_Group=="OK_nc")

#pivot long format

effect_table_sig_long <- effect_table_sig%>%
  pivot_longer(cols = starts_with("status"), names_to = "comparison_status", values_to = "status")%>%
  separate(comparison_status, c("variable" , "Group"), "_")%>%
  mutate(comparison_status=paste(Group, sep="_"))%>%
  select(-c(variable, Group))%>%
  pivot_longer(cols = starts_with("p"), names_to = "comparison_p", values_to = "raw_p")%>%
  separate(comparison_p, c("variable", "Group"), "_")%>%
  mutate(comparison_p=paste(Group, sep="_"))%>%
  select(-c(variable, Group))%>%
  filter(comparison_p==comparison_status)%>%
  pivot_longer(cols = starts_with("effect"), names_to = "comparison_effectSize", values_to = "effectSize")%>%
  separate(comparison_effectSize, c("variable", "Group"), "_")%>%
  mutate(comparison_effectSize=paste(Group, sep="_"))%>%
  select(-c(variable, Group))%>%
  filter(comparison_p==comparison_effectSize)%>%
  pivot_longer(cols = starts_with("corr"), names_to = "comparison_q", values_to = "corr_p")%>%
  separate(comparison_q, c("variable", "Group"), "_")%>%
  mutate(comparison_q=paste(Group, sep="_"))%>%
  select(-c(variable))%>%
  filter(comparison_p==comparison_q)%>%
  select(-c(comparison_q, comparison_effectSize,comparison_status))%>%
  mutate(fdr= as_factor(case_when(corr_p <= 0.05 ~ "*", corr_p <= 0.01 ~ "**", corr_p <= 0.001 ~ "***", corr_p <= 0.1 ~ ".")))

library(gtools)
effect_table_sig_long%>%
  ggplot(aes (x = Group, y = reorder(rowname, effectSize)))+ 
  geom_point (aes (fill = effectSize, shape = as.factor (sign (effectSize)), size = abs (effectSize), color=(fdr))) +
  scale_shape_manual (values = c (25, 24)) + 
  scale_fill_gradient2 (low = "blue", high = "red", mid = "white", midpoint = 0) +
  scale_color_manual(values=c("gray22","gray1","gray85"))+
  geom_text (aes (label = stars.pval(corr_p)))+
  theme_grey() +
  theme(axis.text.x = element_text(), #angle = 0, hjust = 0, vjust = 1
        #To remove legend = legend.position = "none"
        legend.position = "right",axis.text.y = element_text(face = "italic")) + scale_y_discrete(position = "left") + ylab("Functional modules")

#Heatmap
ggplot(effect_table_sig_long, aes(x = Group, y = reorder(rowname, effectSize))) +
  # do the heatmap tile coloring based on effect sizes
  geom_tile(aes(fill = effectSize)) +
  scale_fill_gradient2 (low = "blue", high = "red", mid = "white", midpoint = 0) +
  # add significance stars/circles for deconfounded/confounded associations
  geom_text (aes (label = stars.pval(corr_p)))+
  guides(color = guide_legend(override.aes = list(shape = c(1,8)) ) ) +
  
  # make it pretty
  theme_classic() +
  theme(axis.text.x = element_text(size = 7,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = 0.3),
        axis.text.y = element_text(size = 7,
                                   angle = 0,
                                   hjust = 1,
                                   vjust = 0.35),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0),
        plot.subtitle=element_text(size=8)) +
  labs(title="Groups_limma and modules correlation heatmap",
       subtitle="FDR-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
       x = "Group comparison (hiperammon=1)",
       y = "ASVs")


df <- data.frame(Abundance = t(abundances(ps_limma)),
                 Group = meta(ps_limma)$Group)
p1 <- ggplot(df, aes(x = Group, y = df$Abundance.ASV589)) +
  geom_boxplot() +
  labs(title = "Relative abundances", subtitle = "Add here the taxa name",  y = "Abundance")+ theme(plot.title = element_text(size=18), axis.text.x = element_text(size = 14), axis.title = element_text(size = 16))

p1