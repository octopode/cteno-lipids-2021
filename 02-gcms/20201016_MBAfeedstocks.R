# FA profile matchup analysis for MBA

samp_comp_mba <- samp_comp %>%
  filter(
    sp %in% c(
      "Aure_labi",
      "Holm_cost",
      "Lamp_crue",
      "Leuc_pulc"
    )
  )

# major cpds only
comp_mba_major <- samp_comp_mba %>%
  # filter for even-chain only
  #filter((chain %% 2) == 0) %>%
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
comp_maj_summary <- comp_mba_major %>%
  group_by(sp, id, chain, dbond) %>%
  summarize(
    n = n(),
    avg_frac_molar = mean(frac_molar),
    sem_frac_molar = sd(frac_molar) / sqrt(n())
  ) %>%
  group_by(sp) %>%
  mutate(avg_frac_molar = avg_frac_molar/sum(avg_frac_molar))

# plot composition by species
comp_maj_summary %>%
  ungroup() %>%
  mutate(
    sp = factor(sp, levels = c("Aure_labi", "Holm_cost", "Leuc_pulc", "Lamp_crue")),
    id = factor(id, levels = comp_maj_summary %>% arrange(chain, dbond) %>% pull(id) %>% unique())
  ) %>%
  ggplot(aes(x=sp, y=avg_frac_molar, group=id, fill=id)) +
  geom_col() +
  scale_fill_brewer(palette = "Blues") +
  labs(x="species", y="mole fraction", fill="fatty acid", title="composition by species") +
  guides(fill = guide_legend(title.position = "top")) +
  theme_pubr() +
  theme(legend.position = "right")
