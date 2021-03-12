source("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/20200929_FAME_comp.R")

## FIGURE 3: chain length, DBI vs. environmental variables
fig3_data <- samp_comp_fastats %>% 
  # txome spp only
  filter(sp %in% tree$tip.label) %>% 
  # bad actors
  filter(!(eid %in% c("JWL0136"))) %>%
  mutate(
    # for color-coding
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other")),
    # label cold and shallow samples
    cold = temp_col <= 7.5,
    shal = depth_col <= 200
  ) %>% 
  # melt environmental vars
  pivot_longer(c(depth_col, temp_col), names_to = "env_var", values_to = "env_val") %>% 
  mutate(env_var = factor(env_var, levels=c("temp_col", "depth_col", "do_col"))) %>% 
  # melt indices
  pivot_longer(c(chain, dbi), names_to = "idx_var", values_to = "idx_val") %>% 
  # do depth and temp filtering
  filter(
    # so deep stuff doesn't go in temp dataset
    ((env_var == "temp_col") & shal) |
      # and warm stuff doesn't go in depth dataset
      ((env_var == "depth_col") & cold)
  )

# fits dataframe
fig3_fits <- fig3_data %>% 
  # to drop empty levels
  mutate(env_var = paste(env_var)) %>% 
  group_by(env_var, idx_var) %>% 
  group_split() %>% 
  pblapply(
    cl = 4L,
    .,
    function(xydata){
      xydata %>% 
        pgls_opt(., tree_rooted, idx_val ~ env_val, 
                 # exponential series
                 branch_mults = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()
        ) %>% 
        # get the relevant params
        mapply(FUN=function(name, mod){get_params(mod) %>% mutate(branch_mult = name)}, names(.), .) %>% 
        t() %>% as.tibble() %>% unnest() %>% 
        # select the best branch multiplier
        filter(logLik == max(logLik)) %>% 
        #group_by(logLik) %>% 
        filter(alpha == min(alpha)) %>% 
        # store the covariates
        mutate(
          x = xydata %>% pull(env_var) %>% first(),
          y = xydata %>% pull(idx_var) %>% first()
        )
    }
  ) %>% 
  do.call(rbind, .) %>% 
  dplyr::rename(
    env_var = x,
    idx_var = y
  ) %>% 
  mutate(env_var = factor(env_var, levels=c("temp_col", "depth_col", "do_col")))

pdf(file = "~/Documents/MBARI/Lipids/cteno-lipids-2020/figs/3/Fig3_var5.pdf", width = 85/25.4, height = 90/25.4)
ggplot() +
  facet_grid(
    cols=vars(env_var),
    rows=vars(idx_var),
    scales="free"#,
    #switch="both"
  ) +
  # significance rectangles
  geom_rect(
    # note multiple testing!
    data=fig3_fits %>% mutate(sig = (p<= (0.05/6))), 
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
    aes(fill = sig),
    alpha=0.1
  ) +
  # raw data
  # for temp
  geom_point(
    data = fig3_data %>% filter(env_var == "temp_col"),
    aes(x=env_val, y=idx_val, shape=sprep),
    alpha = 0.4,
    size = 0.75
  ) +
  # for depth
  geom_point(
    data = fig3_data %>% filter(env_var == "depth_col"),
    aes(x=env_val, y=idx_val, shape=sprep),
    alpha = 0.4,
    size = 0.75
  ) +
  # PGLS trendline
  geom_abline(data = fig3_fits, aes(slope=b, intercept=c)) +
  # scales and theme
  theme_pubr() +
  theme(
    legend.position = "none",
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    # hide strips
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  ) +
  #scale_y_continuous(position = "right") +
  #scale_color_manual(values = chroma_sp) +
  scale_shape_manual(values = shapes_sp) +
  scale_fill_manual(values = c(NA, "black")) +
  labs(title = "Figure 3: chain length and double bonds")
dev.off()
