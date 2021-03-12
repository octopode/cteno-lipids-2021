library(tidyverse)
library(ggpubr)

env_data <- read_tsv(
  "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/metadata/old/20201010_txome_metadata.tsv",
  col_names = c("DateTime24", "RovName", "DiveNumber", "sample", "Depth", "Temperature", "O2", "sal", "transmiss", "lat", "lon", "timecode", "imagelink", "samplenum")
)

# keep in mind that due to MCMC annotation condenser, many OGs have identical annotation!
annot_data <- read_tsv(
  "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/og_annots_v8.tab",
  col_names = c("orthogroup", "annot")
)

# total read counts
spot_data <- read_tsv("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/metadata/txome_spotcounts.tsv")

# desaturase and elongase genes
#dir_tpms <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/fads_elov"
# non-peroxisomal "perox" genes
#dir_tpms <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/counts_grep_-i_perox__grep_-v_oxisom"
# ascorbate
#dir_tpms <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/counts_grep_-i_ascorb"
# glutathione
#dir_tpms <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/counts_grep_-i_glutathione"
# plasmalogenic enzymes
dir_tpms <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/counts_grep_-i_dihydroxyacetonephosphate"

# read TPMs
readdata <- list.files(path = dir_tpms, full.names = T, pattern = ".tsv") %>%
  lapply(
    function(file_data){
      data_this_file <- file_data %>%
        read_tsv(col_names=FALSE) %>%
        mutate(file = file_data %>% basename())
      return(data_this_file)
    }
  ) %>%
  do.call(rbind, .) %>%
  magrittr::set_colnames(c("target_id", "length", "eff_length", "est_counts", "tpm", "file_data")) %>%
  rowwise() %>%
  mutate(sample = strsplit(target_id, "|", fixed = TRUE) %>% unlist() %>% .[[1]]) %>%
  # only well-represented orthogroups %>%
  group_by(file_data) %>%
  mutate(nsamp = sample %>% unique() %>% length()) %>%
  filter(nsamp >= 40)

# scrub incomplete transcripts
readdata_complete <- readdata %>%
  group_by(file_data) %>%
  # length should be above median
  filter(length >= quantile(length, 0.5))

# sum isoforms
readdata_gene <- readdata_complete %>%
  group_by(sample, file_data) %>%
  summarize(
    tpm = sum(tpm),
    n_isoforms = n()
  )

# join to environmental data
readdata_env <- readdata_gene %>%
  left_join(env_data, by="sample")

## check read count distributions across txome samples
#readdata_gene %>%
#  right_join(env_data, by="sample") %>%
#  group_by(sample) %>%
#  ggplot(aes(x=sample, y=tpm)) +
#    geom_boxplot(alpha=0.2) +
#    theme_pubr() +
#    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
#
## does seq depth affect TPM?
#readdata_gene %>%
#  group_by(sample) %>%
#  summarise(tpm = mean(tpm)) %>%
#  right_join(spot_data, by="sample") %>%
#  filter(sample != "CTE_Mertensia_ovum1") %>%
#  ggplot(aes(x=reads, y=tpm)) +
#  geom_point(alpha=0.2) +
#  geom_smooth(method="lm") +
#  theme_pubr()

## plot some stuff
### (all data)
#pdf("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/20201011_fads_elov_lm.pdf", width = 11, height = 11)
#readdata_env %>%
#  group_by(file_data) %>%
#  #filter(tpm >= quantile(tpm, 0.25)) %>%
#  # melt it!
#  pivot_longer(c(Depth, Temperature, O2, sal, transmiss), names_to = "env_var", values_to = "env_val") %>%
#  ggplot(aes(x=env_val, y=tpm)) +
#    facet_grid(rows=vars(file_data), cols=vars(env_var), scales = "free") +
#    geom_point(alpha = 0.2) +
#    geom_smooth(method="lm") +
#    theme_pubr()
#dev.off()
#
### shallow
#readdata_env %>%
#  filter(Depth <= 200) %>%
#  # melt it!
#  pivot_longer(c(Depth, Temperature, O2, sal, transmiss), names_to = "env_var", values_to = "env_val") %>%
#  ggplot(aes(x=env_val, y=tpm)) +
#  facet_grid(rows=vars(file_data), cols=vars(env_var), scales = "free") +
#  geom_point(alpha = 0.2) +
#  geom_smooth(method="lm") +
#  theme_pubr()
#
### cold
#readdata_env %>%
#  filter(Temperature <= 7.5) %>%
#  # melt it!
#  pivot_longer(c(Depth, Temperature, O2, sal, transmiss), names_to = "env_var", values_to = "env_val") %>%
#  ggplot(aes(x=env_val, y=tpm)) +
#  facet_grid(rows=vars(file_data), cols=vars(env_var), scales = "free") +
#  geom_point(alpha = 0.2) +
#  geom_smooth(method="lm") +
#  theme_pubr()
#
# load accessory functions
source("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/20200821_FAME_comp_funcs.R")
library(phytools)
file_tree <- "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200819_iq.tree"
tree_long <- file_tree %>% read.newick()
tree_long_rooted <- reroot(tree_long, 53, 0)

# annotation labeller function for facets, takes OG filename
annot_label <- Vectorize(function(og){
  annot_data %>%
    # remove any filename extension
    filter(strsplit(og, '.', fixed=TRUE) %>% unlist() %>% .[[1]] == orthogroup) %>%
    pull(annot) %>%
    # parse just the annotation from the hit string
    strsplit(., ' ', fixed=TRUE) %>% unlist() %>%
    # cut off SwissProt ID and metadata
    .[2:min((str_which(., '=')-1), length(.))] %>%
    paste(., collapse=' ') %>%
    # wrap to fit in plot!
    str_wrap(width=20)
})

## color/shape palettes for data points
## these don't make so much sense for txomes
chroma_sp <- RColorBrewer::brewer.pal(9, "Set1") %>%
  setNames(c(
    "CTE_Lampocteis_cruentiventer",
    "CTE_Beroe_cucumis",
    "CTE_Bolinopsis_infundibulum",
    "CTE_Bolinopsis_vitrea",
    "CTE_Leucothea_pulchra",
    "CTE_Bathocyroe_fosteri",
    "CTE_Pleurobrachia_bachei",
    "CTE_Tjalfiella_pink",
    "other"))

# PGLSs
# main dataframe
expn_data <- readdata_env %>%
  # don't care about anything else, don't want any NAs
  select(sample, file_data, tpm, Depth, Temperature, O2) %>%
  drop_na() %>%
  mutate(sp = substr(sample, 1, nchar(sample)-1)) %>%
  # txome spp only
  filter(sp %in% tree_long_rooted$tip.label) %>%
  mutate(
    ## for color-coding
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("CTE_Lampocteis_cruentiventer", "CTE_Beroe_cucumis", "CTE_Bolinopsis_infundibulum", "CTE_Bolinopsis_vitrea", "CTE_Leucothea_pulchra", "CTE_Bathocyroe_fosteri", "CTE_Pleurobrachia_bachei", "CTE_Tjalfiella_pink", "other")),
    # label cold and shallow samples
    cold = Temperature <= 7.5,
    shal = Depth <= 200
  ) %>%
  #mutate(Depth = log10(Depth)) %>%
  pivot_longer(c(Depth, Temperature, O2), names_to = "env_var", values_to = "env_val") %>%
  # constrain analyses
  filter(
    ((env_var == "O2") & shal) |
      ((env_var == "Temperature") & shal) |
      ((env_var == "Depth") & cold)
  ) %>%
  mutate(env_var = factor(env_var, levels=c("Depth", "Temperature", "O2"))) %>%
  # point shape mapping for replicated species
  #mutate(sp = factor(sp)) %>%
  group_by(file_data, sp, env_var) %>%
  mutate(ptshape = ifelse(n()==1, "asingleton", sp)) %>%
  ungroup() %>%
  mutate(ptshape = factor(ptshape) %>% as.integer())

# fits dataframe
expn_fits <- expn_data %>%
  # to drop empty levels
  mutate(env_var = paste(env_var)) %>%
  group_by(env_var, file_data) %>%
  group_split() %>%
  pblapply(
    cl = 8L,
    .,
    function(xydata){
      xydata %>%
        pgls_opt(., tree_long, tpm ~ env_val,
                 # exponential series
                 branch_mults = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()
        ) %>%
        # get the relevant params
        mapply(FUN=function(name, mod){get_params(mod) %>% mutate(branch_mult = name)}, names(.), .) %>%
        t() %>% as_tibble() %>% unnest() %>%
        # store the covariates
        mutate(
          x = xydata %>% pull(env_var) %>% first(),
          y = xydata %>% pull(file_data) %>% first()
        )
    }
  ) %>%
  # tidy function robust to missing columns
  bind_rows() %>%
  dplyr::rename(
    env_var = x,
    file_data = y
  ) %>%
  # select the best branch multiplier for each fit
  group_by(env_var, file_data) %>%
  filter(logLik == max(logLik)) %>%
  filter(alpha == min(alpha)) %>%
  # recast
  ungroup() %>%
  mutate(env_var = factor(env_var, levels=c("Depth", "Temperature", "O2"))) %>%
  # dedupe
  distinct(env_var, file_data, c, b, p, logLik, AIC, BIC)

# calc significance with Bonferroni correction
n_ogs <- expn_data %>%
  pull(file_data) %>%
  unique() %>%
  length()

alpha_adj <- 0.05/n_ogs
#alpha_adj <- 0.05

pdf(file = "~/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto/20201213_plasm_pgls.pdf", width = 85/25.4, height = 0.6*n_ogs+1)
ggplot() +
  facet_grid(
    cols=vars(env_var),
    rows=vars(file_data),
    scales="free",
    switch="both",
    labeller=labeller(file_data = annot_label)
  ) +
  # significance rectangles
  geom_rect(
    # note Bonferroni correction!
    data=expn_fits %>%
      mutate(sig = (p<= alpha_adj)),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
    aes(fill = sig),
    alpha=0.1
  ) +
  # raw data
  geom_point(
    data = expn_data,
    aes(x=env_val, y=tpm, shape=paste(ptshape), color=sprep),
    alpha = 0.4,
    size = 1.25
  ) +
  geom_text(data = expn_data, aes(x=env_val, y=tpm, label=sample, color=sprep)) + #TEST
  # PGLS trendline
  geom_abline(data = expn_fits, aes(slope=b, intercept=c)) +
  # scales and theme
  theme_pubr() +
  theme(
    legend.position = "none",
    text = element_text(size = 7.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.y.left = element_text(angle = 0)
  ) +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = chroma_sp) +
  scale_shape_manual(values = c("1"=16,"2"=2,"3"=6,"4"=0,"5"=7)) +
  scale_fill_manual(values = c(NA, "black")) +
  labs(
    title = paste("ADAPS expression level\nalpha = ", alpha_adj %>% formatC(., format = "e", digits = 1), sep=''),
    x="environment",
    y="transcripts per million")
dev.off()
