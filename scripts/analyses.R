
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # #             local genetic sex differences                 # # #  
# # #             analyses in manuscript                        # # #  
# # #             2022                                          # # #        
# # #             Emil Uffelmann                                # # #       
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##### set-up##### 

## load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

## load data
univ     <- fread("results/quantitative_traits_univ.txt")
equal    <- fread("results/quantitative_traits_equal.txt")
bivar    <- fread("results/quantitative_traits_bivar.txt")
bivar    <- fread("results/quantitative_traits_bivar.txt")
gene_loc <- fread("/Users/schneemil/genetics/magma_v1.10_mac/gene_locations_NCBI37.3/NCBI37.3.gene.loc")
domain_mapping <- fread(file = "data/2022_04_22_domain_mapping.txt", na.strings = "NA")
domain_mapping <- domain_mapping[domain_mapping$keep == "yes",]

#
#
#
#
#
#

##### local heritability #####

univ_wide <- univ %>%
  pivot_wider(names_from = phen, values_from = c(h2.obs, p)) %>%
  na.omit

##  loci that are significant in both males and females
univ_wide_sig <- univ_wide[univ_wide$p_female <= 0.05 / {2495 * 157} & univ_wide$p_male <= 0.05 / {2495 * 157}, ]

## correlation estimates with and without outlying lipoprotein A, direct bilirubin, and total bilirubin loci
# without outliers
cor.test(univ_wide_sig$h2.obs_female[univ_wide_sig$h2.obs_female < 0.2], univ_wide_sig$h2.obs_male[univ_wide_sig$h2.obs_male < 0.2], na.rm = T)
plot(univ_wide_sig$h2.obs_female[univ_wide_sig$h2.obs_female < 0.2], univ_wide_sig$h2.obs_male[univ_wide_sig$h2.obs_male < 0.2])
# with outliers
cor.test(univ_wide_sig$h2.obs_female, univ_wide_sig$h2.obs_male, na.rm = T)
plot(univ_wide_sig$h2.obs_female, univ_wide_sig$h2.obs_male)

## percentage of significant loci
n_sig_male    <- sum(univ_wide$p_female > 0.05 / {2495 * 157} & univ_wide$p_male <= 0.05 / {2495 * 157}, na.rm = T)
n_sig_female  <- sum(univ_wide$p_female <= 0.05 / {2495 * 157} & univ_wide$p_male > 0.05 / {2495 * 157}, na.rm = T)
n_sig_both    <- sum(univ_wide$p_female <= 0.05 / {2495 * 157} & univ_wide$p_male <= 0.05 / {2495 * 157}, na.rm = T)
n_sig_neither <- sum(univ_wide$p_female > 0.05 / {2495 * 157} & univ_wide$p_male > 0.05 / {2495 * 157}, na.rm = T)

n_sig_female / (n_sig_male + n_sig_female + n_sig_both)
n_sig_male / (n_sig_male + n_sig_female + n_sig_both)
(n_sig_male + n_sig_female) / (n_sig_male + n_sig_female + n_sig_both)

## loci with strong differences between males and females
# testosterone
print(univ_wide[univ_wide$phenotype_code == "30850_raw" & univ_wide$h2.obs_male > 10 * univ_wide$h2.obs_female,], n = 40)
print(univ_wide[univ_wide$phenotype_code == "30850_raw" & univ_wide$h2.obs_female > 10 * univ_wide$h2.obs_male,], n = 40)
# Microalbumin in urine
print(univ_wide[univ_wide$phenotype_code == "30500_raw" & univ_wide$h2.obs_female > 100 * univ_wide$h2.obs_male,], n = 40)
# Rheumatoid factor
print(univ_wide[univ_wide$phenotype_code == "30820_raw" & univ_wide$h2.obs_female > 10 * univ_wide$h2.obs_male,], n = 40)
# Urate
print(univ_wide[univ_wide$phenotype_code == "30880_raw" & univ_wide$h2.obs_female > 0.005,], n = 40)
gene_loc[gene_loc$V2 == 4 & gene_loc$V3 >= 8882617 & gene_loc$V3 <= 11050119,] # genes within the loci

#
#
#
#
#
#

##### global genetic correlations #####

## remove non-significant phenotypes
ldsc_rg_sig <- ldsc_rg[ldsc_rg$p <= 0.05 / 157,]

## phenotypes with at least 10 correlations
bivar_mean <- data.frame(phenotype = as.character(), mean_rho = as.numeric(), median_rho = as.numeric())
phenos <- unique(bivar$phenotype_code)
count <- 0

for (pheno in phenos) {
  
  temp <- bivar %>%
    drop_na() %>%
    filter(phenotype_code == pheno)
  
  ## only look at phenotypes with at least 10 rho estimates
  if (nrow(temp) >= 10) {
    count <- count + 1
    bivar_mean[count, "phenotype"] <- pheno
    bivar_mean[count, "mean_rho"] <- mean(bivar$rho[bivar$phenotype_code == pheno], na.rm = T)
    bivar_mean[count, "median_rho"] <- median(bivar$rho[bivar$phenotype_code == pheno], na.rm = T) 
    bivar_mean[count, "se_rho"] <- sd(bivar$rho[bivar$phenotype_code == pheno], na.rm = T) / sqrt(nrow(temp))
    
  }
}

## mean absolute difference
mean(abs(bivar_ldsc$mean_rho - bivar_ldsc$rg))

#
#
#
#
#
#

##### local genetic correlations #####

## check traits with most significant loci
l <- list()
for (pheno in unique(rho1$phenotype_code)) {
  l[pheno] <- sum(rho1$p[rho1$phenotype_code == pheno] < {0.05 / {157 * 2495}}) 
}

n_sig <- data.frame(phenotype_code = names(l), n_sig = unlist(l), row.names = NULL)
n_sig <- merge(n_sig, domain_mapping[, c("name", "phenotype")], by.x = "phenotype_code", by.y = "phenotype", all.x = T)
n_sig[order(n_sig$n_sig),]; nrow(n_sig[n_sig$n_sig >= 1,])

## significant locus for leg fat percentage
rho1[rho1$phenotype_code == "23111_raw" & rho1$p < 1.5e-7,]
bivar[bivar$phenotype_code == "23111_raw" & bivar$locus == 1,]
ldsc_rg[ldsc_rg$phenotype_code == "23111_raw",]
gene_loc[gene_loc$V2 == 1 & gene_loc$V3 >= 10539 & gene_loc$V3 <= 1376204,] # genes within the loci

## High light scatter reticulocyte percentage
rho1[rho1$phenotype_code == "30290_raw" & rho1$p < 1.5e-7,]
bivar[bivar$phenotype_code == "30290_raw" & bivar$locus == 1644,]
ldsc_rg[ldsc_rg$phenotype_code == "30290_raw",]

#
#
#
#
#
#

##### equality #####

###### number of significant loci 

l <- list()
for (pheno in unique(equal$phenotype_code)) {
  l[pheno] <- sum(equal$p.both[equal$phenotype_code == pheno] < {0.05 / {157 * 2495}}) 
}

n_sig <- data.frame(phenotype_code = names(l), n_sig = unlist(l), row.names = NULL)
n_sig <- merge(n_sig, domain_mapping[, c("name", "phenotype")], by.x = "phenotype_code", by.y = "phenotype", all.x = T)
n_sig[order(n_sig$n_sig),]; nrow(n_sig[n_sig$n_sig >= 1,])

###### larger effects in males

temp_equal <- equal %>%
  group_by(phenotype_code) %>%
  dplyr::summarise(median_ratio_female_to_male_raw_total = median(median_ratio_female_to_male_raw),
                   median_ratio_female_to_male_std_total = median(median_ratio_female_to_male_std),
                   ratio_varY_female_to_male = mean(varY_female) / mean(varY_male)) %>%
  as.data.frame()
  

## number of phenotypes where median raw effect is larger in males
sum(temp_equal$median_ratio_female_to_male_raw_total < 1)
sum(temp_equal$median_ratio_female_to_male_raw_total < 1) / nrow(temp_equal)

## number of phenotypes where median std effect is larger in males
sum(temp_equal$median_ratio_female_to_male_std_total < 1)
sum(temp_equal$median_ratio_female_to_male_std_total < 1) / nrow(temp_equal)

## number of phenotypes where median phenotypic variance is larger in males
sum(temp_equal$ratio_varY_female_to_male < 1)
sum(temp_equal$ratio_varY_female_to_male < 1) / nrow(temp_equal)

## check if I find something similar for the heritabilities

univ_wide <- univ %>%
  pivot_wider(names_from = phen, values_from = c(h2.obs, p)) %>%
  na.omit %>%
  mutate(median_ratio_female_to_male = h2.obs_female / h2.obs_male)

temp_univ <- univ_wide %>%
  group_by(phenotype_code) %>%
  dplyr::summarise(median_ratio_female_to_male_total = median(median_ratio_female_to_male))

sum(temp_univ$median_ratio_female_to_male_total < 1)
sum(temp_univ$median_ratio_female_to_male_total < 1) / nrow(temp_univ)