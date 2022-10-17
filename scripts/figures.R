
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # #             local genetic sex differences                 # # #  
# # #             figures in manuscript                         # # #  
# # #             2022                                          # # #        
# # #             Emil Uffelmann                                # # #       
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##### set-up #####

rm(list = ls())
library(data.table)
library(tidyverse)
library(scales)
library(ggrepel)
library(ggridges)
library(ggpubr)
univ  <- fread("results/quantitative_traits_univ.txt")
equal <- fread("results/quantitative_traits_equal.txt")
rho1  <- fread("results/quantitative_traits_rho1.txt")
bivar <- fread("results/quantitative_traits_bivar.txt")
ldsc_rg <- fread("results/quantitative_traits_ldsc_rg.txt")
domain_mapping <- fread(file = "data/2022_04_22_domain_mapping.txt", na.strings = "NA")
domain_mapping <- domain_mapping[domain_mapping$keep == "yes",]

#
#
#
#
#
#
#
#

##### univ ####

##### heritability scatter plot 

univ_wide <- univ %>%
  pivot_wider(names_from = phen, values_from = c(h2.obs, p)) %>%
  na.omit %>%
  left_join(domain_mapping, by = c("phenotype_code" = "phenotype"))

n_sig_male    <- sum(univ_wide$p_female > 0.05 / {2495 * 157} & univ_wide$p_male <= 0.05 / {2495 * 157}, na.rm = T)
n_sig_female  <- sum(univ_wide$p_female <= 0.05 / {2495 * 157} & univ_wide$p_male > 0.05 / {2495 * 157}, na.rm = T)
n_sig_both    <- sum(univ_wide$p_female <= 0.05 / {2495 * 157} & univ_wide$p_male <= 0.05 / {2495 * 157}, na.rm = T)
n_sig_neither <- sum(univ_wide$p_female > 0.05 / {2495 * 157} & univ_wide$p_male > 0.05 / {2495 * 157}, na.rm = T)

p_univ_scatter <- univ_wide %>%
  ggplot(aes(x = h2.obs_male, y = h2.obs_female)) +
    geom_point(aes(x = h2.obs_male, y = h2.obs_female, colour = "Neither male nor female"), data=subset(univ_wide, p_female >= 0.05 / {2495 * 157} & p_male >= 0.05 / {2495 * 157})) +
    geom_point(aes(x = h2.obs_male, y = h2.obs_female, colour = "Male (n = 1001)"), data=subset(univ_wide, p_female >= 0.05 / {2495 * 157} & p_male < 0.05 / {2495 * 157})) +
    geom_point(aes(x = h2.obs_male, y = h2.obs_female, colour = "Female (n = 1734)"), data=subset(univ_wide, p_female < 0.05 / {2495 * 157} & p_male >= 0.05 / {2495 * 157})) +
    geom_point(aes(x = h2.obs_male, y = h2.obs_female, colour = "Male & female (n = 2653)"), data=subset(univ_wide, p_female < 0.05 / {2495 * 157} & p_male < 0.05 / {2495 * 157})) +
    scale_color_manual(values = c("Neither male nor female" = "#E3E1E2", 
                                  "Male (n = 1001)" = "#2E4881",
                                  "Female (n = 1734)" = "#639504",
                                  "Male & female (n = 2653)" = "#F3C13F")) +
    geom_abline(slope=1, intercept = 0) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    # label loci where heritabilities are much larger in females
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(h2.obs_female > 200 * h2.obs_male & p_female <= 0.05 / {2495 * 157} & p_male > 0.05 / {2495 * 157}, name, "")),
                    aes(label = label), show.legend = FALSE, max.overlaps = Inf, box.padding = 0.5, size = 2, force = 0.3) +
    # label loci significant for females that have heritability of more than 1%
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(h2.obs_female > 0.01 & p_female <= 0.05 / {2495 * 157} & p_male > 0.05 / {2495 * 157}, name, "")),
                    aes(label = label), show.legend = FALSE, max.overlaps = Inf, box.padding = 0.5, size = 2, force = 1.1) +
    # label loci where heritabilities are much larger in males
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(h2.obs_male > 200 * h2.obs_female & p_male <= 0.05 / {2495 * 157} & p_female > 0.05 / {2495 * 157}, name, "")),
                    aes(label = label), show.legend = FALSE, max.overlaps = Inf, box.padding = 0.5, size = 2) +
    # label loci significant in both males and females that have very high heritabilities
    geom_text_repel(data = . %>% 
                      mutate(label = ifelse(h2.obs_male > 0.1 & p_male <= 0.05 / {2495 * 157} & p_female <= 0.05 / {2495 * 157}, name, "")),
                    aes(label = label), show.legend = FALSE, box.padding = 0.5, size = 2, min.segment.length = 0, force = 5) +
    # label Urate
    geom_text_repel(data = . %>% 
                      filter(locus == 610) %>%
                      mutate(label = ifelse(h2.obs_female > 0.09 & phenotype_code == "30880_raw" & p_male <= 0.05 / {2495 * 157} & p_female <= 0.05 / {2495 * 157}, name, "")),
                    aes(label = label), show.legend = FALSE, box.padding = 0.5, size = 2, min.segment.length = 0, force = 5) +
    xlab(expression("Male "*italic("h")^2)) +
    ylab(expression("Female "*italic("h")^2)) +
    theme_classic() +
    guides(colour=guide_legend(title="Significance")) +
    theme(legend.title = element_text(size = 7), 
          legend.text = element_text(size = 7))

ggsave(plot = p_univ_scatter, filename = "plots/univ_scatter.png",
       width = 6, height = 3.5, dpi = 320)

#
#
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

###### lollipop

bivar_ldsc <- merge(bivar_mean, ldsc_rg_sig, by.x = "phenotype", by.y = "phenotype_code", all = F) %>%
  drop_na()
temp <- merge(bivar_ldsc, domain_mapping[, c("phenotype", "name", "domain")], by = "phenotype", all = F) %>%
  mutate(domain = as.factor(domain),
         order = seq(1:nrow(bivar_ldsc))) %>%
  arrange(domain)

p <- ggplot(temp) +
  geom_segment( aes(y=name, yend=name, x=mean_rho, xend=rg), color="grey") +
  geom_point( aes(y=name, x=mean_rho, color="#2E4881"), size=0.7) +
  geom_point( aes(y=name, x=rg, color="#639504"), size=0.7) +
  theme_minimal() +
  #coord_fixed(ratio=1/20) +
  scale_color_identity(name = "",
                       guide = "legend",
                       labels = c("LAVA", "LDSC")) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size=6),
    legend.position="top",
    legend.text = element_text(size=6),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10)
  ) +
  xlab("Local rho") +
  ylab("")

png(paste0("plots/global_correlations_lollipop.png"), 
    width=8,height=13,units="cm",res=500)
print(p)
dev.off()  

###### pointrange

temp1 <- bivar_mean %>%
  select(!median_rho) %>%
  rename(rg = mean_rho, se = se_rho, phenotype_code = phenotype) %>%
  mutate(method = "LAVA")

temp2 <- ldsc_rg_sig %>%
  select(rg, se, phenotype_code) %>%
  mutate(method = "LDSC")

common_phenos <- merge(temp1, temp2, by = "phenotype_code")[, "phenotype_code"]

bivar_ldsc <- rbind(temp1[temp1$phenotype_code %in% common_phenos,], temp2[temp2$phenotype_code %in% common_phenos,]) %>%
  left_join(domain_mapping[, c("phenotype", "name", "domain")], by = c("phenotype_code" = "phenotype")) %>%
  mutate(domain = as.factor(domain)) %>%
  arrange(domain)

p <- ggplot(bivar_ldsc, aes(x = rg, reorder_within(name, rg, domain), group = method, colour = method)) +
  geom_pointrange(aes(xmin = rg - 1.96 * se, xmax = rg + 1.96 * se), size=0.1, alpha = 0.8) +
  scale_y_reordered() +
  geom_vline(xintercept = 1, color = "#D82148", linetype = 2) +
  scale_color_manual(values=c('#2E4881', '#639504')) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size=6),
    legend.position="top",
    legend.text = element_text(size=6),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    legend.title=element_blank(),
    strip.text.y = element_text(size = 6, angle=0)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  facet_grid(facets = "domain", scales = "free_y", space = "free_y") +
  xlab("Global correlation") +
  ylab("")

png(paste0("plots/global_correlations_pointrange.png"), 
    width=9,height=13,units="cm",res=500)
print(p)
dev.off()  

#
#
#
#
#
#
#

##### rho1 #####

###### bar graph number of significant loci 
l <- list()
for (pheno in unique(rho1$phenotype_code)) {
  l[pheno] <- sum(rho1$p[rho1$phenotype_code == pheno] < {0.05 / {157 * 2495}}) 
}

n_sig <- data.frame(phenotype_code = names(l), n_sig = unlist(l), row.names = NULL)
n_sig <- merge(n_sig, domain_mapping[, c("name", "phenotype")], by.x = "phenotype_code", by.y = "phenotype", all.x = T)
n_sig[order(n_sig$n_sig),]; nrow(n_sig[n_sig$n_sig >= 1,])

bar_plot <- n_sig %>%
  filter(n_sig > 0) %>%
  ggplot(aes(x = reorder(name, n_sig), y = n_sig)) + 
  geom_bar(stat = "identity", colour="black", fill = "#BEBADA")  +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size=10),
    legend.position="top",
    legend.text = element_text(size=10),
    legend.margin=margin(0,0,0,0),
    #   panel.background = element_rect(fill = "transparent"), # bg of the panel
    #   plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #   legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    #   legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  ylab(expression(atop("Number of significant loci", "(rho1 test)"))) +
  xlab("Phenotype")

ggsave(plot = bar_plot, filename = "plots/rho1_number_sig_loci_bar.png",
       width = 12, height = 15, units = "cm", dpi = 500, bg = "white")

#
#
#
#
#
#
#

##### local genetic correlations & equality & univ #####

domain_mapping <- fread(file = "data/2022_04_22_domain_mapping.txt", na.strings = "NA")
domain_mapping <- domain_mapping %>%
  filter(domain_mapping$keep == "yes") %>%
  arrange(domain) %>%
  mutate(order = row_number()) %>%
  as.data.frame()

temp_univ <- left_join(univ, domain_mapping[, c("domain", "name", "phenotype", "order")],by = c("phenotype_code" = "phenotype")) %>%
  pivot_wider(names_from = phen, values_from = c(h2.obs, p)) %>%
  arrange(domain) %>%
  mutate(domain = droplevels(as.factor(domain)),
         h2obs_ratio_female_to_male = h2.obs_female / h2.obs_male) %>%
  filter(!is.na(name), !is.na(h2obs_ratio_female_to_male), phenotype_code %in% traits) %>%
  as.data.frame()

# selecting all traits with at least 1 significant locus in the equality test (and other traits of interest)
traits <- unique(
  c(equal$phenotype_code[equal$phenotype_code %in% equal$phenotype_code[equal$p.both < 0.05 / {2495 * 157}]])) 
traits <- traits[traits %in% rho1$phenotype_code]

## equality data
temp_equal <- left_join(equal, domain_mapping[, c("domain", "name", "phenotype", "order")], by = c("phenotype_code" = "phenotype")) %>%
  arrange(domain) %>%
  mutate(domain = droplevels(as.factor(domain))) %>%
  filter(phenotype_code %in% traits) %>%
  as.data.frame()

legend_ord <- rev(levels(with(droplevels(temp_equal), reorder(domain, order)))) # to get the legend in the right order

p_equal_raw <- temp_equal %>%
  ggplot(aes(x = median_ratio_female_to_male_raw, y = reorder(name, order), fill = domain)) +
  geom_density_ridges(stat = "binline", scale = 0.85, bins = 50) +
  coord_cartesian(clip = "off") + # prevents top part of figure to be cut-off
  geom_point(size = 2, colour = "#686868", shape = "|", show.legend = F) +
  geom_point(data=subset(temp_equal, p.both < 0.05 / {2495 * 157}), fill="#F3C13F", colour = "#DD0000", shape = 21, size=2, show.legend = F) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10), labels = c(0.01, 0.1, 1, 10)) +
  xlab('Ratio of genetic effects (F:M)') +
  ylab('') +
  geom_vline(xintercept = 1, color = "#D82148") +
  scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"),
                    breaks=legend_ord) +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(fill = "Domain")

m <- merge(temp_equal[, c("locus", "chr", "start", "stop", "phenotype_code", "domain", "name", "median_ratio_female_to_male_raw", "median_ratio_female_to_male_std", "varY_female", "varY_male", "order")],
           temp_univ[, c("locus", "chr", "start", "stop", "phenotype_code", "domain", "name", "order", "h2obs_ratio_female_to_male", "p_female", "p_male")],
           by = c("locus", "chr", "start", "stop", "phenotype_code", "domain", "name", "order")) %>%
  mutate(pheno_ratio = varY_female / varY_male)

m_mean <- m %>%
  filter(h2obs_ratio_female_to_male < 100 & h2obs_ratio_female_to_male > 0.01 & {p_female <= 1e-4 | p_male <= 1e-4}) %>% # filter loci with some genetic signal in one of the sexes
  group_by(name) %>% 
  summarise(pheno_ratio = median(pheno_ratio), 
            median_ratio_female_to_male_raw = median(median_ratio_female_to_male_raw),
            median_ratio_female_to_male_std = median(median_ratio_female_to_male_std),
            h2obs_ratio_female_to_male = median(h2obs_ratio_female_to_male)) %>%
  left_join(domain_mapping[, c("domain", "name", "order")], by = "name")

p_equal_raw_median <- m %>%
  filter(h2obs_ratio_female_to_male < 100 & h2obs_ratio_female_to_male > 0.01 & {p_female <= 1e-4 | p_male <= 1e-4}) %>% 
  ggplot(aes(x = median_ratio_female_to_male_raw, y = reorder(name, order), fill = domain)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, scale = 0.95, show.legend = F) +
  geom_text(
    data = m_mean,
    aes(x = 0.1, y = reorder(name, order), label = round(median_ratio_female_to_male_raw, 2)),
    size = 3, nudge_y = 0.3, colour = "black"
  ) +
  geom_vline(xintercept = 1, color = "#D82148", alpha = 0.6) +
  coord_cartesian(clip = "off") + # prevents top part of figure to be cut-off
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10), labels = c(0.01, 0.1, 1, 10)) +
  xlab('Ratio of genetic effects (F:M)') +
  ylab('') +
  geom_vline(xintercept = 1, color = "#D82148") +
  scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"),
                    breaks=legend_ord) +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(fill = "Domain")

p_equal_std_median <- m %>%
  filter(h2obs_ratio_female_to_male < 100 & h2obs_ratio_female_to_male > 0.01 & {p_female <= 1e-4 | p_male <= 1e-4}) %>% 
  ggplot(aes(x = median_ratio_female_to_male_std, y = reorder(name, order), fill = domain)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, scale = 0.95) +
  geom_text(
    data = m_mean,
    aes(x = 0.35, y = reorder(name, order), label = round(median_ratio_female_to_male_std, 2)),
    size = 3, nudge_y = 0.3, colour = "black"
  ) +
  geom_vline(xintercept = 1, color = "#D82148", alpha = 0.6) +
  coord_cartesian(clip = "off") + # prevents top part of figure to be cut-off
  scale_x_log10(breaks = c(0.5, 1, 2), labels = c(0.5, 1, 2)) +
  xlab('Ratio of genetic effects (F:M)') +
  ylab('') +
  geom_vline(xintercept = 1, color = "#D82148") +
  scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"),
                    breaks=legend_ord) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.title = element_text(size = 9)
  ) +
  labs(fill = "Domain")

ggsave(plot = p_equal_std_median, filename = "plots/ratio_std_genetic_effects.png", width =7, height = 10, bg = "white", dpi = 500)


##### univ h2 data 

p_univ <- m %>%
  filter(h2obs_ratio_female_to_male < 100 & h2obs_ratio_female_to_male > 0.01 & {p_female <= 1e-4 | p_male <= 1e-4}) %>%
  ggplot(aes(x = h2obs_ratio_female_to_male, y = reorder(name, order), fill = domain)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, scale = 0.95, show.legend = F) +
  geom_text(
    data = m_mean,
    aes(x = 0.01, y = name, label = round(h2obs_ratio_female_to_male, 2)),
    size = 3, nudge_y = 0.3, colour = "black"
  ) +
  geom_text(
    data = m_mean,
    aes(x = 50, y = reorder(name, order), label = round(pheno_ratio, 2)),
    size = 3, nudge_y = 0.3, colour = "#D82148", alpha = 0.7
  ) +
  geom_vline(xintercept = 1, color = "#D82148", alpha = 0.6) +
  coord_cartesian(clip = "off") + # prevents top part of figure to be cut-off  
  scale_x_log10(breaks=c(.001, 0.1, 1, 10, 1000),labels=c(.001, 0.1, 1, 10, 1000)) +
  xlab('Ratio of heritabilities (F:M)') +
  ylab('') +
  scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"), aesthetics = "fill",
                    breaks=legend_ord) +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.title = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(fill = "Domain")

## rho1 data

temp_rho1 <- left_join(rho1, domain_mapping[, c("domain", "name", "phenotype", "order")], by = c("phenotype_code" = "phenotype")) %>%
  arrange(domain) %>%
  mutate(domain = droplevels(as.factor(domain))) %>%
  filter(phenotype_code %in% traits) %>%
  as.data.frame()

p_rho1 <- temp_rho1 %>%
  ggplot(aes(x = rho, y = reorder(name, order), fill = domain)) +
  geom_density_ridges(stat = "binline", show.legend = F, scale = 0.95, bins = 15) +
  geom_point(size = 2, colour = "#686868", shape = "|", show.legend = F) +
  geom_point(data=subset(temp_rho1, p < 0.05 / {2495 * 157}), fill="#F3C13F", colour = "#DD0000", shape = 21, size=2, show.legend = F) +
  geom_vline(xintercept = 1, color = "#D82148") +
  coord_cartesian(clip = "off") +
  xlab('Local genetic correlation') +
  ylab('') +
  theme_minimal() +
  scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"),
                    breaks=legend_ord) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.title = element_text(size = 9)
  )

## save plot
figure <- ggarrange(p_rho1, p_equal_raw, p_equal_raw_median, p_univ,
                    labels = c("a", "b", "c", "d"), common.legend = T, legend = "right",
                    ncol = 4, nrow = 1, font.label = list(size = 10, family="Helvetica", face = "bold"),
                    align = "h", widths = c(0.4, 0.2, 0.2, 0.2))

ggsave(plot = figure, filename = "plots/equal_rho1_univ.png", width =14, height = 10, bg = "white", dpi = 500)

#
#
#
#
#
#
#

##### equality #####

###### bar graph number of significant loci 
l <- list()
for (pheno in unique(equal$phenotype_code)) {
  l[pheno] <- sum(equal$p.both[equal$phenotype_code == pheno] < {0.05 / {157 * 2495}}) 
}

n_sig <- data.frame(phenotype_code = names(l), n_sig = unlist(l), row.names = NULL)
n_sig <- merge(n_sig, domain_mapping[, c("name", "phenotype")], by.x = "phenotype_code", by.y = "phenotype", all.x = T)
n_sig[order(n_sig$n_sig),]; nrow(n_sig[n_sig$n_sig >= 1,])

bar_plot <- n_sig %>%
  filter(n_sig > 0) %>%
  ggplot(aes(x = reorder(name, n_sig), y = n_sig)) + 
  geom_bar(stat = "identity", colour="black", fill = "#BEBADA")  +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size=10),
    legend.position="top",
    legend.text = element_text(size=10),
    legend.margin=margin(0,0,0,0),
    #   panel.background = element_rect(fill = "transparent"), # bg of the panel
    #   plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #   legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    #   legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  ylab(expression(atop("Number of significant loci", "(Equality test)"))) +
  xlab("Phenotype")

ggsave(plot = bar_plot, filename = "plots/equal_number_sig_loci_bar.png",
       width = 12, height = 15, units = "cm", dpi = 500, bg = "white")


##### scatter plot of mean ratio of genetic effects vs ratio of phenotype variance per trait 

temp_equal <- equal %>%
  group_by(phenotype_code) %>%
  dplyr::summarise(median_ratio_female_to_male_raw_total = median(median_ratio_female_to_male_raw),
                   median_ratio_female_to_male_std_total = median(median_ratio_female_to_male_std),
                   ratio_varY_female_to_male = mean(varY_female) / mean(varY_male))

equal_domain <- merge(temp_equal, domain_mapping[, c("domain", "name", "phenotype")], 
                               by.x = "phenotype_code", 
                               by.y = "phenotype", all.x = T)

p_raw <- equal_domain %>%
  ggplot(aes(x = median_ratio_female_to_male_raw_total, y = ratio_varY_female_to_male)) +
  geom_point(color = ifelse(equal_domain$median_ratio_female_to_male_raw_total < 10^-0.4 | 
                              equal_domain$median_ratio_female_to_male_raw_total > 10^0.4, "red", "grey50")) +
  geom_text_repel(data = . %>% 
                    mutate(label = ifelse(median_ratio_female_to_male_raw_total < 10^-0.4 | median_ratio_female_to_male_raw_total > 10^0.4, name, "")),
                  aes(label = label), show.legend = FALSE, max.overlaps = Inf, nudge_x = 0.15, size = 2.5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ylab("Ratio female to male phenotypic variance") +
  xlab("Ratio of female to male genetic effects") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size=8)
  )

ggsave(plot = p, filename = "plots/equal_ratio_equal_pheno_scatter.png",
       width = 11, height = 9, units = "cm", dpi = 500, bg = "white")


##### density plot of ratio of heritability, phenotypic variance, genetic effects 

temp_equal <- equal %>%
  group_by(phenotype_code) %>%
  dplyr::summarise(median_ratio_female_to_male_raw_total = median(median_ratio_female_to_male_raw),
                   median_ratio_female_to_male_std_total = median(median_ratio_female_to_male_std),
                   ratio_varY_female_to_male = mean(varY_female) / mean(varY_male)) %>%
 # filter(ratio_varY_female_to_male < 10 & ratio_varY_female_to_male > 0.1) %>% # there is one strong outlier
  as.data.frame()

temp <- univ %>%
  pivot_wider(names_from = phen, values_from = c(h2.obs, p)) %>%
  na.omit %>%
  mutate(median_ratio_female_to_male = h2.obs_female / h2.obs_male) %>%
  group_by(phenotype_code) %>%
  dplyr::summarise(median_ratio_female_to_male_total = median(median_ratio_female_to_male)) %>%
  left_join(temp_equal, by = "phenotype_code") %>%
  pivot_longer(cols = c("median_ratio_female_to_male_total", "median_ratio_female_to_male_raw_total", "median_ratio_female_to_male_std_total", "ratio_varY_female_to_male"),
               names_to = "source")

p <- temp %>%
 # filter(source != "median_ratio_female_to_male_std_total") %>%
  ggplot(aes(x=value, group=source, fill=source)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_minimal()
p
