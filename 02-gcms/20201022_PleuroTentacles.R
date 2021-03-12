# FA profile matchup analysis for MBA

samp_comp_pleuro <- samp_comp %>%
  filter(sp == "Pleu_bach")

# major cpds only
comp_pleuro_major <- samp_comp_pleuro %>%
  # filter for even-chain only
  filter((chain %% 2) == 0) %>%
  group_by(file_data) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
  # average across samples
  group_by(sp, id) %>%
  mutate(avg_frac_molar = mean(frac_molar)) %>%
  group_by(id) %>%
  filter(max(avg_frac_molar) >= 0.025) %>%
  select(-avg_frac_molar) %>%
  # renormalize
  group_by(file_data) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar))

# summarize by major compound
comp_pleuro_summary <- comp_pleuro_major %>%
  group_by(tissue, id, chain, dbond) %>%
  summarize(
    avg_frac_molar = mean(frac_molar),
    sem_frac_molar = sd(frac_molar) / sqrt(n())
  ) %>%
  group_by(tissue) %>%
  mutate(avg_frac_molar = avg_frac_molar/sum(avg_frac_molar))

# plot composition by species
comp_pleuro_summary %>%
  ungroup() %>%
  mutate(
    tissue = factor(tissue, levels = c("body", "whole", "tentacle")),
    id = factor(id, levels = comp_pleuro_summary %>% arrange(chain, dbond) %>% pull(id) %>% unique())
  ) %>%
  ggplot(aes(x=tissue, y=avg_frac_molar, group=id, fill=id)) +
  geom_col() +
  scale_fill_brewer(palette = "Blues") +
  labs(x="sample", y="mole fraction", fill="fatty acid", title="Pleurobrachia bachei\ncomposition by tissue") +
  guides(fill = guide_legend(title.position = "top")) +
  theme_pubr() +
  theme(legend.position = "right")
