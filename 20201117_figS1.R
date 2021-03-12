samp_comp_ctenos_factor <- samp_comp_ctenos %>%
  # bad actors
  filter(!(eid %in% c("JWL0136", "JWL0150", "JWL0102", "JWL0147"))) %>%
  ungroup() %>%
  mutate(id = factor(id, levels = samp_comp_ctenos %>% arrange(chain, dbond) %>% pull(id) %>% unique()))

### which were detected in significant amounts?
sig_cpds <- samp_comp_ctenos_factor %>%
  group_by(id) %>%
  mutate(
    ttest = t.test(x=frac_molar, alternative="greater") %>% list(),
    pval = ttest[[1]]$p.value,
    min_sig_frac = ttest[[1]]$conf.int %>% min()
  ) %>%
  group_by(id, chain, pval, min_sig_frac) %>%
  summarize() %>%
  ungroup() %>%
  # Bonferroni correction
  mutate(pval_adj = pval*n()) %>%
  rowwise() %>%
  mutate(sig = ifelse(pval_adj<=0.05, max(samp_comp_ctenos$frac_molar),0))

# plot mole fraction distros and significance

pdf(file = "~/Documents/MBARI/Lipids/cteno-lipids-2020/figs/S1/FigS1_var1.pdf", width = 190/25.4, height = 85/25.4)
samp_comp_ctenos_factor %>%
  ggplot(aes(x = id, y = frac_molar)) +
  geom_col(data = sig_cpds, fill="black", alpha=0.4, width=1, aes(x=id, y=sig)) +
  geom_boxplot(outlier.alpha = 0.4) +
  #geom_text(aes(label=eid)) +
  xlab("compound") +
  ylab("mole fraction") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "mole fraction by compound, all samples")
dev.off()
