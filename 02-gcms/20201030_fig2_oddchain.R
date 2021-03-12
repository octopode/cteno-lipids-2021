#source("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/20200929_FAME_comp.R")
source("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/20201025_FAME_comp.R")

## FIGURE 2: FA composition and environmental correlations
# main dataframe
fig2_odd_data <- samp_comp_ctenos %>%
  # bad actors
  filter(!(eid %in% c("JWL0136", "JWL0150", "JWL0141"))) %>%
  # major cpds only
  # average across samples
  group_by(id) %>%
  mutate(avg_frac_molar = median(frac_molar)) %>%
  #filter(first(avg_frac_molar) >= 0.005) %>%
  # odd-chain only
  filter(chain %% 2 > 0) %>%
  select(-avg_frac_molar) %>%
  # renormalize
  #group_by(file_data) %>%
  #mutate(frac_molar = frac_molar/sum(frac_molar)) %>%
  # remove odd-chain >50%
  #filter(max(frac_molar) <= 0.5) %>%
  # txome spp only
  filter(sp %in% tree$tip.label) %>%
  # add a column for the total
  group_by(id) %>%
  mutate(tot = mean(c(min(frac_molar), max(frac_molar)))) %>%
  group_by(eid) %>%
  mutate(
    # for color-coding
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other")),
    # label cold and shallow samples
    cold = temp_col <= 7.5,
    shal = depth_col <= 200
  ) %>%
  # log-transform!
  #mutate(depth_col = log10(depth_col)) %>%
  # melt environmental vars
  pivot_longer(c(depth_col, temp_col, do_col, tot), names_to = "env_var", values_to = "env_val") %>%
  mutate(env_var = factor(env_var, levels=c("tot", "temp_col", "depth_col", "do_col"))) %>%
  # do depth and temp filtering
  filter(
    (env_var == "tot") |
      # model supposing UV required for autox
      ((env_var == "do_col") & shal) |
      # so deep stuff doesn't go in temp dataset
      ((env_var == "temp_col") & shal) |
      # and warm stuff doesn't go in depth dataset
      ((env_var == "depth_col") & cold)
  )

# fits dataframe
fig2_odd_fits <- fig2_odd_data %>%
  # to drop empty levels
  mutate(env_var = paste(env_var)) %>%
  group_by(env_var, id) %>%
  filter(env_var != "tot") %>%
  #filter(env_var != "do_col") %>% # temporary; DO models don't converge!
  group_split() %>%
  pblapply(
    cl = 8L,
    .,
    function(xydata){
      xydata %>%
        pgls_opt(., tree_rooted, frac_molar ~ env_val,
                 # exponential series
                 branch_mults = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()
        ) %>%
        # get the relevant params
        mapply(FUN=function(name, mod){get_params(mod) %>% mutate(branch_mult = name)}, names(.), .) %>%
        t() %>% as.tibble() %>% unnest() %>%
        # select the best branch multiplier
        filter(logLik == max(logLik)) %>%
        filter(alpha == min(alpha)) %>%
        # store the covariates
        mutate(
          x = xydata %>% pull(env_var) %>% first(),
          y = xydata %>% pull(id) %>% first(),
        )
    }
  ) %>%
  do.call(rbind, .) %>%
  dplyr::rename(
    env_var = x,
    id = y
  ) %>%
  mutate(env_var = factor(env_var, levels=c("tot", "temp_col", "depth_col", "do_col")))

# calc significance with Bonferroni correction
n_fames <- fig2_odd_data %>%
  pull(id) %>%
  unique() %>%
  length()

#alpha_adj <- 0.05/n_fames
alpha_adj <- 0.05/1

pdf(file = "~/Documents/MBARI/Lipids/cteno-lipids-2020/figs/2/Fig2_oddchain_var1.pdf", width = 110/25.4, height = 150/25.4)
ggplot() +
  facet_grid(
    cols=vars(env_var),
    rows=vars(id),
    scales="free",
    switch="both"
  ) +
  # significance rectangles
  geom_rect(
    # note multiple testing!
    data=fig2_odd_fits %>% mutate(sig = (p<= (0.05/20))),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
    aes(fill = sig),
    alpha=0.2
  ) +
  # raw data
  # total abundance
  geom_boxplot(
    data = fig2_odd_data %>% filter(env_var == "tot"),
    aes(x=frac_molar, y=env_val),
    orientation = "y",
    width = 0.001,
    alpha=0.4,
    size=0.5
  ) +
  # for temp, depth, DO
  geom_point(
    data = fig2_odd_data %>% filter(env_var != "tot"),
    aes(x=env_val, y=frac_molar, color=sprep),
    alpha = 0.4,
    size = 0.5
  ) +
  #geom_text(data = fig2_odd_data, aes(x=env_val, y=frac_molar, label=eid, color=sprep)) + #TEST
  # PGLS trendline
  geom_abline(data = fig2_odd_fits, aes(slope=b, intercept=c)) +
  # scales and theme
  theme_pubr() +
  theme(
    legend.position = "none",
    text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = chroma_sp) +
  #scale_shape_manual(values = shapes_sp) +
  scale_fill_manual(values = c(NA, "black")) +
  labs(title = paste("odd-chain fatty acids\nalpha = ", alpha_adj %>% formatC(., format = "e", digits = 1), sep=''), y="mole fraction")
dev.off()
