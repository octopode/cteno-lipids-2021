source("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/20200929_FAME_comp.R")

# load long-format CTD data
ctd_data <- read.csv("/Users/jwinnikoff/Documents/MBARI/2019_SICB/ICB_manuscript/figs/CtenoDiversity/20190312_CTD_SV-MB-HI.tsv", sep='\t') %>% 
  as_tibble() %>% 
  # drop depths outside specified range
  filter((depth >= 3) & (depth <= 4000)) %>% 
  mutate(temp = as.numeric(paste(temp)))

# FIGURE 1
## SPECIES BY HABITAT
pdf(file = "~/Documents/MBARI/Lipids/cteno-lipids-2020/figs/1/Fig1_var6.pdf", width = 85/25.4, height = 85/25.4)
samp_fastats_phylo %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=depth_col, color=sprep)) +
  # CTD traces
  geom_smooth(
    data=ctd_data, 
    aes(y=depth, x=temp, linetype=source), 
    color="black", 
    size=0.5,
    span=0.7, 
    method='loess', 
    se=FALSE,
    orientation = "y"
  ) +
  geom_jitter(alpha = 0.4, size = 2, width=1, height=10) +
  theme_pubr() +
  scale_color_manual(
    values = chroma_sp,
    labels = Vectorize(function(sprep){
      samp_comp_fastats %>%
        filter(sp == sprep) %>%
        nrow() %>%
        paste(short2full[sprep], " n=", ., sep="")
    }),
    guide = guide_legend(
      direction = "vertical",
      title.position = "top"
    )
  ) +
  scale_y_reverse() +
  labs(x = "Temperature (deg C)", y = "Depth (m)", color = "Species") +
  theme(
    #legend.position = c(0.7, 0.45),
    legend.position = "none",
    text = element_text(size = 8),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(linetype=FALSE)
#scale_y_continuous(trans=reverselog_trans(10))
dev.off()

## PHENOGRAMS!

#samp_phylo_summary <- samp_fastats_phylo %>%
#  group_by(sp) %>%
#  summarise(
#    n = n(),
#    temp_col = mean(temp_col),
#    depth_col = mean(depth_col),
#    dbi = mean(dbi),
#    chain = mean(chain)
#  )

prune_to_one <- function(pheno, phylo){
  # construct branch length table
  branches <- cbind(phylo$edge, phylo$edge.length) %>%
    as_tibble() %>%
    set_names(c("mrca", "node", "dist")) %>%
    mutate(
      node = as.integer(node),
      mrca = as.integer(mrca)
    )
  
  # get duplicated taxa and their intraspecific MRCAs
  leaves <- tibble(sp = phylo$tip.label) %>%
    mutate(node = row_number())
  
  leaves_dup <- leaves %>%
    group_by(sp) %>%
    mutate(count = n()) %>%
    filter(count > 1) %>%
    mutate(mrca = phylo %>% getMRCA(node) %>% as.integer())
  
  # duplicate taxa and their average distance from MRCA
  taxa_dup <- leaves_dup %>%
    rowwise() %>%
    mutate(dist = tibble(node, mrca) %>% left_join(branches, by = c("node", "mrca")) %>% .$dist) %>%
    group_by(sp, mrca) %>%
    summarize(
      node = min(node),
      dist = mean(dist, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # leaves to drop (all but the first, per sp)
  leaves_drop <- leaves_dup %>%
    anti_join(leaves_dup %>% summarize(node = first(node)))
  
  # set first conspecific branch to average length
  branches_adj <- bind_rows(anti_join(branches, taxa_dup, by = c("node", "mrca")), taxa_dup)
  
  # replace conspecific clades with a single tip
  phylo %>%
    compute.brlen(
      phylo$edge %>%
        as_tibble() %>%
        set_names(c("mrca", "node")) %>%
        # put in order
        left_join(branches_adj, by = c("mrca", "node")) %>%
        pull(dist)
    ) %>%
    drop.tip(leaves_drop %>% pull(node)) %>%
    # finally, prune taxa not in phenodata
    drop.tip(phylo$tip.label %>% .[which(!(. %in% pheno$sp))])
}

fig1_data <- samp_fastats_phylo %>% 
  # bad actors
  filter(!(eid %in% c("JWL0136")))

tree_pruned <- prune_to_one(fig1_data, tree_rooted)
fastats_phylo_summary <- fig1_data %>%
  group_by(sp) %>%
  summarize(
    count = n(),
    depth_med = median(depth_col),
    temp_med = median(temp_col),
    depth_mean = mean(depth_col),
    temp_mean = mean(temp_col),
    depth_sd = sd(depth_col),
    temp_sd = sd(temp_col)
  )

# Depth phenogram
phenogram(
  tree_pruned %>% force.ultrametric(method = 'nnls'),
  fastats_phylo_summary$depth_mean %>% desc() %>% setNames(fastats_phylo_summary$sp),
  spread.cost=c(1,0)
)

phenogram(
  tree_pruned %>% force.ultrametric(method = 'nnls'),
  fastats_phylo_summary$temp_mean %>% setNames(fastats_phylo_summary$sp),
  spread.cost=c(1,0)
)

nsim = 1000
anc <- fastBM(tree_pruned)
fancyTree(tree_pruned,type="phenogram95",x=anc)