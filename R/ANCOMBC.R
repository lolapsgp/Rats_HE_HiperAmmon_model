library(ANCOMBC)
library(tidyr)
library(tibble)
library(purrr)
require(phyloseq)
require(tidyverse)
require(magrittr)

ps_0 <- readRDS("~/Documents/Doctorado/Estudios/CIPF_2022/Analisis_resultados_ratas_hiperammon/outside_git_outputs_and_inputs/outputs/ps_0.Rds")

wh0 <- genefilter_sample(ps_0, filterfun_sample(function(x) x > 5), A=0.1*nsamples(ps_0))
ps_sub<- prune_taxa(wh0, ps_0)

ps_sub <- tax_glom(ps_sub, taxrank = "Genus")

res_bc <- ancombc(ps_sub,                    # phyloseq object
                  formula = "Group + Batch",        # fixed effects only
                  p_adj_method = "BH",
                  lib_cut = 1000,                   # min library size
                  group = "Group",                  # primary contrast
                  struc_zero = TRUE,                # structural zeros
                  neg_lb = TRUE,
                  alpha = 0.05,
                  global = FALSE,
)

res = res_bc$res

library(dplyr)

# Create a tidy table with only GroupHiperammonemic results
group_table <- data.frame(
  ASV = rownames(res$beta),
  Beta_Group = res$beta$GroupHiperammonemic,
  SE_Group = res$se$GroupHiperammonemic,
  W_Group = res$W$GroupHiperammonemic,
  Pval_Group = res$p_va$GroupHiperammonemic,
  Qval_Group = res$q_va$GroupHiperammonemic,
  Diff_Abn = res$diff_abn$GroupHiperammonemic
)

# View the table
group_table

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
latest_annotations <- get_latest_annotation(ps_0)

group_table$TaxaID <- latest_annotations$TaxaID[match(group_table$ASV, latest_annotations$ASV)]



# Significant results (FDR < 0.05)
sig_group <- group_table %>%
  filter(Diff_Abn == TRUE)

sig_group_strong <- sig_group 

library(ggplot2)
library(dplyr)
library(ggtext)  # for significance stars

# Prepare data in long format
effect_table_sig_long <- sig_group_strong %>%
  mutate(
    Group = "Hiperammonemic",  # vs Control reference
    effectSize = Beta_Group,
    se = SE_Group,
    fdr = Qval_Group,
    corr_p = Pval_Group,
    sig_stars = case_when(
      Qval_Group < 0.001 ~ "***",
      Qval_Group < 0.01 ~ "**",
      Qval_Group < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(TaxaID, Group, effectSize, se, fdr, corr_p, sig_stars)

# Create volcano-style effect plot (adapted from your template)
effect_table_sig_long %>%
  ggplot(aes(x = Group, y = reorder(TaxaID, effectSize))) +
  geom_point(aes(fill = effectSize, 
                 shape = as.factor(sign(effectSize)), 
                 size = abs(effectSize)), 
             alpha = 0.8, stroke = 0.5) +
  scale_shape_manual(values = c("1" = 25, "-1" = 24), 
                     name = "Direction", 
                     labels = c("↑ Control", "↑ HA")) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, name = "logFC") +
  scale_size_continuous(range = c(2, 8), name = "logFoldChange") +
  geom_text(aes(label = sig_stars), 
            hjust = -0.3, size = 4, fontface = "bold") +
  theme_grey(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 0),
        axis.text.y = element_text(face = "italic", size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank()) +
  labs(y = "ASVs", x = "Group", 
       title = "Significant ASVs (FDR < 0.05)") +
  scale_y_discrete(position = "right")

# ASVs increased in Hiperammonemic (your disease group of interest)
ha_increased <- group_table %>%
  filter(Beta_Group > 0, Qval_Group < 0.05, Diff_Abn == TRUE) %>%
  arrange(desc(Beta_Group))

# ASVs increased in Control  
control_increased <- group_table %>%
  filter(Beta_Group < 0, Qval_Group < 0.05, Diff_Abn == TRUE) %>%
  arrange(Beta_Group)  # most negative first

# Check counts
nrow(ha_increased)
nrow(control_increased)

# Top 5 from each
head(ha_increased, 5)
head(control_increased, 5)




