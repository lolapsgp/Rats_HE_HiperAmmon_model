# Check density

# Control
density_boot_Control <- numeric(n_boot)
ps_filtered2 <- readRDS("~/Lola/Rats_HE_HiperAmmon_model/outputs/ps_filtered2.Rds")
ps_genus <- tax_glom(ps_filtered2, taxrank = "Genus")
ps_genus_Control <- subset_samples(ps_genus, Group == "Control")
for(i in 1:n_boot){
  
  samp_ids <- sample(sample_names(ps_genus_Control),
                     replace = TRUE)
  
  ps_boot_C <- prune_samples(samp_ids, ps_genus_Control)
  
  otu_boot_C <- as.data.frame(otu_table(ps_boot_C))
  otu_boot_C <- t(otu_boot_C)
  
  net_boot_C <- netConstruct(otu_boot_C,
                           measure = "sparcc",
                           normMethod = "clr",
                           zeroMethod = "none",
                           sparsMethod = "t-test",
                           seed = 123456)
  
  props_boot_C <- netAnalyze(net_boot_C)
  
  density_boot_Control[i] <- props_boot_C$globalProps$density1
}



#HA
density_boot_HA <- numeric(n_boot)
ps_filtered2 <- readRDS("~/Lola/Rats_HE_HiperAmmon_model/outputs/ps_filtered2.Rds")
ps_genus <- tax_glom(ps_filtered2, taxrank = "Genus")
ps_genus_HA <- subset_samples(ps_genus, Group == "Hiperammonemic")
for(i in 1:n_boot){
  
  samp_ids <- sample(sample_names(ps_genus_HA),
     replace = TRUE)
  
  ps_boot <- prune_samples(samp_ids, ps_genus_HA)
  
  otu_boot <- as.data.frame(otu_table(ps_boot))
  otu_boot <- t(otu_boot)
  
  net_boot <- netConstruct(otu_boot,
   measure = "sparcc",
   normMethod = "clr",
   zeroMethod = "none",
   sparsMethod = "t-test",
   seed = 123456)
  
  props_boot <- netAnalyze(net_boot)
  
  density_boot_HA[i] <- props_boot$globalProps$density1
}


#Compare
wilcox.test(density_boot_HA, density_boot_Control)
#Wilcoxon rank sum test with continuity correction

#data:  density_boot_HA and density_boot_Control
#W = 26627, p-value = 8.983e-09
#alternative hypothesis: true location shift is not equal to 0
mean(density_boot_HA)
#[1] 0.005245408
mean(density_boot_Control)
#[1] 0.003965673
median(density_boot_HA)
#[1] 0.005269497
median(density_boot_Control)
#[1] 0.0009033424
quantile(density_boot_HA, c(0.025, 0.975))
#2.5%        97.5% 
#0.0003011141 0.0120445649 
quantile(density_boot_Control, c(0.025, 0.975))
#2.5%        97.5% 
#0.0003011141 0.0225835592

#The higher density of the HA network remained consistent across bootstrap resampling


#Other tests
library(phyloseq)
library(NetCoMi)

# Keep only taxa present in all samples
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 0, ps_genus)

# Fix taxa order once
taxa_fixed <- taxa_names(ps_genus)


set.seed(123)
n_perm <- 500
density_diff_null <- numeric()
successful <- 0

sample_data(ps_genus)$Group_original <- sample_data(ps_genus)$Group

for(i in 1:n_perm){
  
  result <- tryCatch({
    
    # Shuffle labels
    sample_data(ps_genus)$Group_perm <- 
      sample(sample_data(ps_genus)$Group_original)
    
    ps_perm_C  <- subset_samples(ps_genus, Group_perm == "Control")
    ps_perm_HA <- subset_samples(ps_genus, Group_perm == "Hiperammonemic")
    
    otu_perm_C  <- t(as.data.frame(otu_table(ps_perm_C)))
    otu_perm_HA <- t(as.data.frame(otu_table(ps_perm_HA)))
    
    net_perm_C <- netConstruct(otu_perm_C,
                               measure = "sparcc",
                               normMethod = "clr",
                               zeroMethod = "none",
                               sparsMethod = "t-test",
                               seed = 123456)
    
    net_perm_HA <- netConstruct(otu_perm_HA,
                                measure = "sparcc",
                                normMethod = "clr",
                                zeroMethod = "none",
                                sparsMethod = "t-test",
                                seed = 123456)
    
    dens_perm_C  <- netAnalyze(net_perm_C)$globalProps$density1
    dens_perm_HA <- netAnalyze(net_perm_HA)$globalProps$density1
    
    dens_perm_HA - dens_perm_C
    
  }, error = function(e) {
    return(NA)
  })
  
  if(!is.na(result)){
    density_diff_null <- c(density_diff_null, result)
    successful <- successful + 1
  }
  
  if(i %% 50 == 0){
    cat("Completed:", i, " Successful:", successful, "\n")
  }
}

cat("Total successful permutations:", successful, "\n")

obs_diff<-mean(density_boot_HA - density_boot_Control)
# Empirical p-value
p_empirical <- mean(abs(density_diff_null) >= abs(obs_diff))
p_empirical

#[1] 0.7470238

