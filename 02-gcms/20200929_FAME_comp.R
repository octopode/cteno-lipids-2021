library(GGally)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(tidyverse)
library(ggpubr)
library(pbapply)
# load accessory functions
source("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/20200821_FAME_comp_funcs.R")

# convert JRW filename to sample code
filename2samp <- function(filename){
  prefix <- substr(filename, 0, 3)
  samp <- str_split(filename, "[_-]") %>%
    .[[1]] %>% .[1] %>%
    substr(4, str_length(.)) %>%
    str_pad(., 4, "0", side = "left") %>%
    paste(prefix, ., sep="")
  return(samp)
}

# convert lipid name to sat class
colon2fa <- Vectorize(function(cpd){
  if(str_detect(cpd, ":0.*-ME")) {return("SFA")}
  if(str_detect(cpd, ":1.*-ME")) {return("MUFA")}
  if(str_detect(cpd, ":[2-9].*-ME")) {return("PUFA")}
  else {return("other")}
})

# convert lipid name to double bond count
colon2db <- Vectorize(function(cpd){
  str_extract(cpd, ":.") %>%
    str_remove_all("[:]") %>%
    as.integer()
})

# convert lipid name to chain length
colon2cl <- Vectorize(function(cpd){
  str_extract(cpd, "C.*:") %>%
    str_remove_all("[C:]") %>%
    as.integer()
})

# Convert long taxon codes to short
long2short <- Vectorize(function(name){
  name %>%
    str_split("_") %>% unlist() %>%
    .[(length(.)-1): length(.)] %>%
    str_remove("sp") %>%
    str_remove("cf") %>%
    substr(1, 4) %>%
    paste(collapse = "_") %>%
    str_remove("[0-9]")
})

## LOAD PHYLOGENETIC INFORMATION
file_tree <- "/Users/jwinnikoff/Documents/MBARI/Lipids/GCMSData/cdf/20200819_iq.tree"
tree <- file_tree %>% read.newick()
tree$tip.label <- long2short(tree$tip.label) %>% unname()# %>% setNames(., .)
tree_rooted <- reroot(tree, 53, 0)

## LOAD SAMPLE METADATA
data_file <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/metadata/20200816_lipid_metadata.tsv"
envi_data <- read_tsv(data_file)

# This is the main DF to play with!
samp_comp <- comp_unks %>%
  rowwise() %>%
  mutate(eid = filename2samp(file_data %>% basename())) %>%
  filter(str_detect(eid, "JWL") == TRUE) %>% # samples only; no standards
  filter(stds_mix != "BHT") %>% # targets only; no contaminants
  # get the best file for each sample
  group_by(file_data) %>%
  mutate(tot_ions = sum(intb)) %>%
  group_by(eid) %>%
  filter(tot_ions == max(tot_ions)) %>%
  # join to enviro data
  left_join(envi_data, by = "eid") %>%
  ungroup() %>%
  complete(eid, id, fill = list(frac_molar = 0)) %>%
  # samples that never passed QC and that I should just chuck
  filter(!(eid %in% c(
    "JWL0013",
    "JWL0051",
    "JWL0037"
  ))) %>%
  # a little patch-fix
  mutate(id = ifelse(id == "C18:1(E)-ME", "C18:1-ME", id)) %>%
  group_by(file_data, id) %>%
  mutate(frac_molar = sum(frac_molar)) %>%
  filter(roi == last(roi)) %>%
  distinct(eid, id, sp, tissue, depth_col, temp_col, frac_molar) %>%
  # and some compound annotation
  mutate(
    chain = colon2cl(id),
    dbond = colon2db(id),
    satcl = colon2fa(id)
  )

# ctenos only!
samp_comp_ctenos <- samp_comp %>%
  filter(
    !is.na(sp) &
      sp != "Aure_labi" &
      sp != "Holm_cost" &
      tissue == "whole"
  )

# try to focus on eukaryotic lipids
samp_comp_evenfa <- samp_comp_ctenos %>%
  filter(
    ((chain %% 2) == 0) & # even chain
      (chain >= 14) & # medium/long chain
      str_detect(id, "ME") %>% unlist() # FAME
  ) %>%
  # renormalize
  group_by(eid) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar))

# calc samplewise DBI and mean chain length
samp_comp_fastats <- samp_comp_evenfa %>%
  group_by(eid, sp, temp_col, depth_col) %>%
  mutate(
    dbi = dbond * frac_molar,
    chain = chain * frac_molar
  ) %>%
  summarize(
    dbi = sum(dbi),
    chain = sum(chain)
  ) %>%
  # added 20200925 so we can have individual mole fractions!
  left_join(
    samp_comp_evenfa %>%
      pivot_wider(.,
          id_cols = c("eid"),
          names_from = id,
          values_from = frac_molar
      ),
    by = "eid"
  ) %>%
  replace(is.na(.), 0)

## filter to species in the phylogeny
samp_fastats_phylo <- samp_comp_fastats %>%
  filter(sp %in% tree_rooted$tip.label)

# just fatty alcs
samp_fatalcs_phylo <- samp_comp_ctenos %>%
  filter(id %in% c("C20:1-OH", "C22:1-OH"))

## INSPECT RAW DATA

### by compound across all samples
samp_comp_ctenos %>%
  ggplot(aes(x = id, y = frac_molar)) +
    geom_boxplot(alpha = 0.2) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### which were detected in significant amounts?
sig_cpds <- samp_comp_ctenos %>%
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
  mutate(pval_adj = pval*n())

# how many are significant?
sig_cpds %>%
  filter(pval_adj <= 0.05) %>%
  nrow()

# odd-chain\
samp_comp_ctenos %>% group_by(id) %>% filter(chain %in% c(15,17)) %>% filter(frac_molar == max(frac_molar))

areas_unks %>% filter(str_detect(file_data, "JWL") == TRUE) %>% group_by(id) %>% filter(id %in% c("C15:0-ME","C15:1-ME","C17:0-ME","C17:1-ME")) %>% filter(intensity == max(intensity))

### correlation between compounds
# plotting function called by ggpairs()
linregs <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.2) +
    geom_smooth(
      #aes(color = stats::cor.test(data$x, data$y, method = "pearson", use = "pairwise.complete.obs")$cor),
      formula = y~x,
      method=lm,
      ...
      )
}

pairwise_ggmtx_all <- samp_comp_ctenos %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c("file_data", "eid", "temp_col", "depth_col"),
    names_from = id,
    values_from = frac_molar
  ) %>%
  ggpairs(
    columns = colnames(.) %>% .[which(str_detect(.,":") %>% unlist())],
    lower = list(continuous = linregs)
    ) +
    theme_pubr() +
    labs(x = "mole %", y = "mole %") +
    ggtitle("Pairwise correlations")

#ggsave("~/Downloads/20200816_pairwise_allcpds.pdf", width = w, height = h, pairwise_ggmtx_all, width = 40, height = 40, units="in", limitsize=FALSE)

# vs depth
samp_comp_ctenos %>%
  group_by(id) %>%
  mutate(frac_molar_norm = frac_molar / max(frac_molar)) %>%
  ggplot(aes(x = depth_col, y = frac_molar_norm)) +
    facet_wrap(~id) +
    geom_point(alpha = 0.4) +
    geom_smooth(color = "blue", method = "lm") +
    scale_x_log10() +
    theme_pubr() +
    xlab("collection depth (m)") +
    ylab("normalized molar fraction") +
    ggtitle("normalized molar fraction vs. depth")

# vs temp
samp_comp_ctenos %>%
  group_by(id) %>%
  mutate(frac_molar_norm = frac_molar / max(frac_molar)) %>%
  ggplot(aes(x = temp_col, y = frac_molar_norm)) +
  geom_smooth(color = "red", method = "lm") +
  facet_wrap(~id) +
  geom_point(alpha = 0.4) +
  theme_pubr() +
  xlab("collection temp (deg C)") +
  ylab("normalized molar fraction") +
  ggtitle("normalized molar fraction vs. temperature")

## SUBSET DATA



# drop whole samples with outlier contents
samp_comp_evenfa_inliers <- samp_comp_evenfa %>%
  group_by(id) %>%
  mutate(
    summary = summary(frac_molar) %>% list(),
    iqr = IQR(frac_molar),
    outlier = (frac_molar < unlist(first(summary))[["1st Qu."]] - 1.5*iqr) | (frac_molar > unlist(first(summary))[["3rd Qu."]] + 1.5*iqr)
    ) %>%
  # alternative SD-based outlier definition
  #mutate(
  #  outlier = abs(frac_molar - mean(frac_molar)) > 2.5*sd(frac_molar)
  #) %>%
  group_by(eid) %>%
  filter(!any(outlier))

### pairwise of these
pairwise_ggmtx_evenfa <- samp_comp_evenfa %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c("file_data", "eid", "temp_col", "depth_col"),
    names_from = id,
    values_from = frac_molar
  ) %>%
  ggpairs(
    columns = colnames(.) %>% .[which(str_detect(.,":") %>% unlist())] %>% .[1:5],
    lower = list(continuous = linregs)
  ) +
  theme_pubr() +
  labs(x = "mole %", y = "mole %") +
  ggtitle("Pairwise correlations")

#ggsave("~/Downloads/20200816_pairwise_evenfa_inliers.pdf", width = w, height = h, pairwise_ggmtx_evenfa, width = 20, height = 20, units="in", limitsize=FALSE)

# temp and depth
# species added in as a survey
samp_comp_fastats %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(count > 6, sp, "other")
  ) %>%
  ggplot(aes(x=temp_col, y=depth_col, color=dbi)) +
    geom_jitter(alpha = 0.4, width=1, height=50) +
    theme_pubr() +
    scale_color_distiller(palette = "YlGnBu", direction = 1) +
    scale_y_log10() +
    scale_y_reverse()

# vs. temp
samp_comp_fastats %>%
  ggplot(aes(x = temp_col, y = dbi)) +
    geom_point(alpha = 0.2) +
    theme_pubr() +
    ggtitle("DBI vs. temperature")

# vs. depth
samp_comp_fastats %>%
  ggplot(aes(x = depth_col, y = dbi)) +
  geom_point(alpha = 0.2) +
  theme_pubr() +
  ggtitle("DBI vs. depth")

## slices
# vs. temp
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  ggplot(aes(x = temp_col, y = dbi)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method="lm", color="red") +
  theme_pubr() +
  ggtitle("DBI vs. temperature, depth <= shal m")

# vs. depth
samp_comp_fastats %>%
  filter(temp_col <= cold) %>%
  ggplot(aes(x = depth_col, y = dbi)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method="lm", color="blue") +
  theme_pubr() +
  ggtitle("DBI vs. depth, temperature <= cold deg C")

# vs. temp
samp_comp_fastats %>%
  ggplot(aes(x = temp_col, y = chain)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method="lm", color="red") +
  theme_pubr() +
  ggtitle("Chain length vs. temperature")

# vs. depth
samp_comp_fastats %>%
  ggplot(aes(x = depth_col, y = chain)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method="lm", color="blue") +
  theme_pubr() +
  ggtitle("Chain length vs. depth")

## slices
# vs. temp
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  ggplot(aes(x = temp_col, y = chain, shape = sp)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method="lm", color="red") +
  theme_pubr() +
  ggtitle("Chain length vs. temperature, depth <= shal m")

# vs. depth
samp_comp_fastats %>%
  filter(temp_col <= cold) %>%
  ggplot(aes(x = depth_col, y = chain)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method="lm", color="blue") +
  theme_pubr() +
  ggtitle("Chain length vs. depth, temperature <= cold deg C")

### by FA class
#### by species
samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "MUFA", "PUFA"))) %>%
  group_by(eid, sp, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  group_by(sp, satcl) %>%
  mutate(count = n()) %>%
  filter(count > 1) %>%
  ggplot(aes(x = paste(sp, " n=", count, sep=""), y = frac_molar)) +
  facet_grid(rows = vars(satcl)) +
  geom_boxplot(alpha = 0.2) +
  theme_pubr() +
  labs(x = "species", y = "mole fraction", title = "Saturation class by species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#### vs. temp
samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "MUFA", "PUFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = temp_col, y = frac_molar)) +
    facet_grid(cols = vars(satcl)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm") +
    theme_pubr()

#### vs. depth
samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "MUFA", "PUFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = depth_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_pubr()

### Stenothermal and stenobaric slices
#### vs. temp
evenfa_shallow <- samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "MUFA", "PUFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  filter(depth_col <= shal)

evenfa_shallow %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = temp_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", color = "red") +
  theme_pubr() +
  labs(title = paste("FA class vs. temp, n =", evenfa_shallow %>% pull(eid) %>% unique() %>% length()))

#### vs. depth
evenfa_cold <- samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "MUFA", "PUFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  filter(temp_col <= cold)

evenfa_cold %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = depth_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  ggtitle(paste("FA class vs. depth, n =", evenfa_cold %>% pull(eid) %>% unique() %>% length()))

### combined UFAs
#### vs. temp
samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = ifelse(str_detect(satcl, "UFA") %>% unlist(), "UFA", "SFA")) %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "UFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = temp_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_pubr()

#### vs. depth
samp_comp_evenfa %>%
  ungroup() %>%
  mutate(satcl = ifelse(str_detect(satcl, "UFA") %>% unlist(), "UFA", "SFA")) %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "UFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = depth_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_pubr()

### slices
#### vs. temp
evenfa_shallow %>%
  ungroup() %>%
  mutate(satcl = ifelse(str_detect(satcl, "UFA") %>% unlist(), "UFA", "SFA")) %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "UFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = temp_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_pubr()

#### vs. depth
evenfa_cold %>%
  ungroup() %>%
  mutate(satcl = ifelse(str_detect(satcl, "UFA") %>% unlist(), "UFA", "SFA")) %>%
  mutate(satcl = factor(satcl, levels=c("SFA", "UFA"))) %>%
  group_by(eid, temp_col, depth_col, satcl) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  ggplot(aes(x = depth_col, y = frac_molar)) +
  facet_grid(cols = vars(satcl)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_pubr()

# now, only the ones >1% mole fraction
samp_comp_evenfa_major <- samp_comp_evenfa %>%
  ungroup() %>%
  filter(frac_molar >= 0.01) %>%
  # renormalize
  group_by(eid) %>%
  mutate(frac_molar = frac_molar/sum(frac_molar))

# PCA
comp_wide <- samp_comp_evenfa_major %>%
  pivot_wider(
    id_cols = c("file_data", "eid", "temp_col", "depth_col"),
    names_from = id,
    values_from = frac_molar
  ) %>%
  replace(is.na(.), 0) %>%
  ungroup()

# run PCA
fame_pca <- comp_wide %>%
  ungroup() %>%
  select(contains(":")) %>%
  # if scaling is used, minor constituents *MUST* be forced to 0!
  prcomp(center = TRUE, scale = FALSE)

# which PCs correlate best with depth and temperature?
cor(fame_pca$x, comp_wide %>% select(depth_col), use="complete.obs") %>% tibble() %>% # depth corrs
  bind_cols(cor(fame_pca$x, comp_wide %>% select(temp_col), use="complete.obs") %>% tibble()) %>% # temp corrs
  bind_cols(fame_pca$sdev^2 %>% tibble()) %>% # eigenvalues
  setNames(c("depth", "temp", "eigen")) %>%
  mutate(eigen = eigen / sum(eigen)) %>% # properly normalize the eigenvalues
  # add the PC#s as a column
  mutate(PC = seq(length(colnames(fame_pca$x)))) %>%
  gather("param", "weight", -PC) %>%
  mutate(param = factor(.$param, levels = c("eigen", "depth", "temp"))) %>%
  filter(PC <= 4) %>% # only show first n PCs
  ggplot(aes(x = PC, y = weight, fill = param)) +
  geom_bar(position="dodge", stat="identity") +
  scale_x_continuous(breaks = seq(1, length(colnames(fame_pca$x)), by = 1)) +
  scale_fill_brewer(palette = "Dark2") +
  ggtitle("eigenvalues and correlation strengths\nacross principal components") +
  theme_pubr() +
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  theme(legend.position = "right")

# biplot with just vectors
axes = c(2,3) # the PCs to plot
fviz_pca_var(fame_pca, axes = axes, repel = T, col.var = "darkgrey") %>%
  # to log or not to log?
  #fviz_add(., 1*cor(comp_wide %>% select(depth_col),
  fviz_add(., 1*cor(log(comp_wide %>% select(depth_col) + 1),
                    fame_pca$x, use="complete.obs"),
           axes = axes, color ="#6278AB", geom="arrow", labelsize = 3, linetype = "solid") %>%
  fviz_add(., 1*cor(comp_wide %>% select(temp_col),
                    fame_pca$x, use="complete.obs"),
           axes = axes, color ="#BE2D55", geom="arrow", labelsize = 3, linetype = "solid") +
  theme_pubr() +
  ggtitle("fatty acid PC loadings\ndepth and temp supplementary")

# SHOWPLOTS!
chroma_sp <- RColorBrewer::brewer.pal(9, "Set1") %>%
  setNames(c(
    "Lamp_crue",
    "Bero_cucu",
    "Boli_infu",
    "Boli_vitr",
    "Leuc_pulc",
    "Bath_fost",
    "Pleu_bach",
    "Tjal_pink",
    "other"))

shapes_sp <- c(
    "Lamp_crue" = 15, # filled square
    "Bero_cucu" = 18, # filled diamond
    "Boli_infu" = 06, # dn tri open
    "Boli_vitr" = 02, # up tri open
    "Leuc_pulc" = 00, # open square
    "Bath_fost" = 17, # up triangle filled
    "Pleu_bach" = 16, # filled circle
    "Tjal_pink" = 07, # X square
    "other"     = 03  # +
    )

# from https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
library(scales)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

short2full = c(
  "Lamp_crue" = "Lampocteis cruentiventer",
  "Bero_cucu" = "Beroe cucumis",
  "Boli_infu" = "Bolinopsis infundibulum",
  "Boli_vitr" = "Bolinopsis vitrea",
  "Leuc_pulc" = "Leucothea pulchra",
  "Bath_fost" = "Bathocyroe fosteri",
  "Pleu_bach" = "Pleurobrachia bachei",
  "Tjal_pink" = "undescribed platyctene",
  "other" = "other"
)

## species by habitat
pdf(file = "~/Documents/MBARI/2020_eDSBS/sp_by_hab.pdf", width = 6, height = 4.6)
samp_comp_fastats %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=depth_col, color=sprep)) +
  geom_jitter(alpha = 0.4, size = 4, width=1, height=10) +
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
  theme(legend.position = c(0.6, 0.4))
  #scale_y_continuous(trans=reverselog_trans(10))
dev.off()

 ## DBI by habitat
samp_comp_fastats %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=depth_col, color=dbi)) +
  geom_jitter(alpha = 0.4, size = 4, width=1, height=10) +
  geom_text(aes(label=eid)) +
  theme_pubr() +
  scale_color_distiller(
    palette = "YlGnBu", direction = 1,
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  scale_y_reverse() +
  labs(x = "Temperature (deg C)", y = "Depth (m)", color = "Mean # double bonds") +
  theme(legend.position = c(0.6, 0.3))
  #scale_y_continuous(trans=reverselog_trans(10))

## chain length by habitat
samp_comp_fastats %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=depth_col, color=chain)) +
  geom_jitter(alpha = 0.4, size = 4, width=1, height=10) +
  theme_pubr() +
  scale_color_distiller(
    palette = "YlGnBu", direction = 1,
    #breaks = c(16, 18, 20),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  scale_y_reverse() +
  labs(x = "Temperature (deg C)", y = "Depth (m)", color = "Mean chain length") +
  theme(legend.position = c(0.6, 0.3))
#scale_y_continuous(trans=reverselog_trans(10))

## DBI by temp
pdf(file = "~/Documents/MBARI/2020_eDSBS/dbi_vs_temp_shal.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=dbi, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_smooth(method = "lm", color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 2.5) +
  labs(x = "Temperature (deg C)", y = "Mean # double bonds") +
  theme(legend.position = "none")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/dbi_vs_depth_cold.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=depth_col, y=dbi, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  #geom_text(aes(label = eid)) +
  geom_smooth(method = "lm", color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 2.5) +
  labs(x = "Depth (m)", y = "Mean # double bonds") +
  theme(legend.position = "none")
dev.off()

## chain length by temp
pdf(file = "~/Documents/MBARI/2020_eDSBS/chain_vs_temp_shal.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=chain, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  #geom_text(aes(label = eid)) +
  geom_smooth(method = "lm", color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 20) +
  labs(x = "Temperature (deg C)", y = "Mean chain length") +
  theme(legend.position = "none")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/chain_vs_depth_cold.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=depth_col, y=chain, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_smooth(method = "lm", color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 20) +
  labs(x = "Depth (m)", y = "Mean chain length") +
  theme(legend.position = "none")
dev.off()

## INTRASPECIFIC

## DBI by temp
# get p-vals
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  filter(n() >= 6) %>%
  group_split() %>%
  pblapply(
    .,
    function(frame){
      frame %>%
        lm(., formula=dbi ~ temp_col) %>%
        summary() %>%
        .$coefficients %>%
        as_tibble() %>%
        pull(`Pr(>|t|)`)
    }
  )

pdf(file = "~/Documents/MBARI/2020_eDSBS/dbi_vs_temp_shal_intrasp.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  filter(!(sprep %in% c("other"))) %>%
  ggplot(aes(x=temp_col, y=dbi, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 2.5) +
  labs(x = "Temperature (deg C)", y = "Mean # double bonds") +
  theme(legend.position = "none")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/dbi_vs_depth_cold_intrasp.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  filter(!(sprep %in% c("other", "Bero_cucu"))) %>%
  ggplot(aes(x=depth_col, y=dbi, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  #geom_text(aes(label = eid)) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 2.5) +
  labs(x = "Depth (m)", y = "Mean # double bonds") +
  theme(legend.position = "none")
dev.off()

## chain length by temp
# get p-vals
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  filter(n() >= 6) %>%
  group_split() %>%
  pblapply(
    .,
    function(frame){
      frame %>%
        lm(., formula=chain ~ temp_col) %>%
        summary() %>%
        .$coefficients %>%
        as_tibble() %>%
        pull(`Pr(>|t|)`) %>%
        last()
    }
  ) %>% unlist()

pdf(file = "~/Documents/MBARI/2020_eDSBS/chain_vs_temp_shal_intrasp.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  filter(!(sprep %in% c("other"))) %>%
  ggplot(aes(x=temp_col, y=chain, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  #geom_text(aes(label = eid)) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 20) +
  labs(x = "Temperature (deg C)", y = "Mean chain length") +
  theme(legend.position = "none")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/chain_vs_depth_cold_intrasp.pdf", width = w, height = h)
samp_comp_fastats %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  filter(!(sprep %in% c("other", "Bero_cucu"))) %>%
  ggplot(aes(x=depth_col, y=chain, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_smooth(method = "lm") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 20) +
  labs(x = "Depth (m)", y = "Mean chain length") +
  theme(legend.position = "none")
dev.off()

# composition
### by compound across all samples
pdf(file = "~/Documents/MBARI/2020_eDSBS/major_cpds_frac.pdf", width = 3.5, height = 3.5)
samp_comp_evenfa %>%
  ungroup() %>%
  mutate(
    id = factor(id),
    id = factor(id, levels = rev(levels(id)))
  ) %>%
  group_by(id) %>%
  mutate(frac_molar_med = median(frac_molar)) %>%
  filter(frac_molar_med >= 0.025) %>%
  ggplot(aes(y = id, x = frac_molar)) +
  geom_boxplot(alpha = 0.2) +
  theme_pubr() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Fatty acid", x = "Mole fraction")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/major_cpds_db.pdf", width = 4, height = h)
samp_comp_evenfa %>%
  group_by(id) %>%
  mutate(frac_molar_med = median(frac_molar)) %>%
  filter(frac_molar_med >= 0.025) %>%
  group_by(id, chain, dbond) %>%
  summarise() %>%
  ggplot(aes(y = id, x = dbond)) +
  geom_col(alpha = 0.5) +
  theme_pubr() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Fatty acid", x = "# double bonds")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/major_cpds_cl.pdf", width = 4, height = h)
samp_comp_evenfa %>%
  group_by(id) %>%
  mutate(frac_molar_med = median(frac_molar)) %>%
  filter(frac_molar_med >= 0.025) %>%
  group_by(id, chain, dbond) %>%
  summarise() %>%
  ggplot(aes(y = id, x = chain)) +
  geom_col(alpha = 0.4) +
  theme_pubr() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Fatty acid", x = "Chain length")
dev.off()

# temp
pdf(file = "~/Documents/MBARI/2020_eDSBS/major_cpds_temp_shal.pdf", width = 8, height = 2)
samp_comp_evenfa %>%
  group_by(id) %>%
  mutate(frac_molar_med = median(frac_molar)) %>%
  filter(frac_molar_med >= 0.025) %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = length(unique(eid)),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other")
  ) %>%
  ggplot(aes(y = frac_molar, x = temp_col, color = sprep)) +
  facet_wrap(~id, ncol=6) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_manual(values = chroma_sp) +
  theme_pubr() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Mole fraction", x = "Temperature (deg C)") +
  theme(legend.position = "none")
dev.off()

# depth
pdf(file = "~/Documents/MBARI/2020_eDSBS/major_cpds_depth_cold.pdf", width = 8, height = 2.25)
samp_comp_evenfa %>%
  group_by(id) %>%
  mutate(frac_molar_med = median(frac_molar)) %>%
  filter(frac_molar_med >= 0.025) %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = length(unique(eid)),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other")
  ) %>%
  ggplot(aes(y = frac_molar, x = depth_col, color = sprep)) +
  facet_wrap(~id, ncol=6) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  scale_color_manual(values = chroma_sp) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(y = "Mole fraction", x = "Depth (m)") +
  theme(legend.position = "none")
dev.off()

# PHYLOGENETIC CORRECTION!

# problems here: cannot correlate vars with themselves!
# need explicit independents and dependents
cors <- corphylo_indls(
  samp_fastats_phylo %>%
    ungroup() %>%
    #select(-eid, -temp_col, -depth_col) %>%
    select(sp, temp_col, depth_col, `C18:1-ME`, `C18:0-ME`),
  tree
)

samp_comp_top6 <- samp_comp_evenfa %>%
  # major cpds only
  filter(id %in% c("C14:0-ME", "C16:0-ME", "C18:0-ME", "C18:1-ME", "C20:5-ME", "C22:6-ME")) %>%
  arrange(id)

# run pairwise correlations on all samples
pair_models <- pgls_pairwise(samp_comp_top6, tree)

# get the interesting stuff out of those models!
pars_models <- lapply(pair_models, get_params) %>%
  bind_rows() %>%
  dplyr::rename(id.x=x, id.y=y)

# make a grid of scatterplots, GGAlly style, with PGLS trendlines
pair_frac_molar <- samp_comp_top6 %>%
  left_join(samp_comp_top6, by=c("eid", "sp")) %>%
  select(
    eid,
    sp,
    id.x,
    id.y,
    frac_molar.x,
    frac_molar.y,
  ) %>%
  mutate(
    # can switch x and y here to control which margin the distributions go on
    # control whether distributions are on the margin or diagonal
    id.y = ifelse(id.y == id.x, "distribution", id.y),
    id.x = factor(make.names(id.x), levels=samp_comp_top6 %>% pull(id) %>% unique() %>% make.names()),
    # control whether all samples put at top or bottom
    id.y = factor(make.names(id.y), levels=samp_comp_top6 %>% pull(id) %>% unique() %>% make.names() %>% c(., "distribution"))
  ) %>%
  ungroup() %>%
  filter(as.numeric(id.y) > as.numeric(id.x)) # no diagonal
  #filter(as.numeric(id.y) >= as.numeric(id.x)) # with diagonal

# DBI, chain correls

# this routine scans over evolutionary rates to ensure convergence
temp<-samp_fastats_phylo %>%
  filter(temp_col <= cold) %>%
  #filter(depth_col <= shal) %>%
  pgls_opt(., tree_rooted, chain ~ depth_col,
           # exponential series
           branch_mults = lapply(seq(0, 20, 1), function(x){10**x}) %>% unlist()
  ) %>%
  # get the relevant params
  mapply(FUN=function(name, mod){get_params(mod) %>% mutate(branch_mult = name)}, names(.), .) %>%
  t() %>% as.tibble()

# DBI vs. temp, shallow
# significant at p = 3.1E-2
pgls_dbi_temp <- samp_fastats_phylo %>%
  filter(depth_col <= shal) %>%
  pgls_ou(., tree_rooted, dbi ~ temp_col, branch_mult = 1E10)
summary(pgls_dbi_temp)

# DBI vs. depth, cold
# significant at p = 2.8E-4
pgls_dbi_depth <- samp_fastats_phylo %>%
  filter(temp_col <= cold) %>%
  pgls_ou(., tree, dbi ~ depth_col, branch_mult = 1E10)
summary(pgls_dbi_depth)

# chain vs. temp, shallow
# significant at p = 2.9E-7
pgls_chain_temp <- samp_fastats_phylo %>%
  filter(depth_col <= shal) %>%
  pgls_ou(., tree_rooted, chain ~ temp_col, branch_mult = 1E10)
summary(pgls_chain_temp)

# chain vs. depth, cold
# significant at p = 7.7E-3
pgls_chain_depth <- samp_fastats_phylo %>%
  filter(temp_col <= cold) %>%
  pgls_ou(., tree_rooted, chain ~ depth_col, branch_mult = 1E10)
summary(pgls_chain_depth)

# total fatty alcohols vs. temp, shallow
pgls_fatalcs_temp <- samp_fatalcs_phylo %>%
  group_by(eid, sp, depth_col, temp_col) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  filter(depth_col <= shal) %>%
  pgls_ou(., tree_rooted, frac_molar ~ temp_col, branch_mult = 1E10)
summary(pgls_fatalcs_temp)$tTable

# total fatty alcohols vs. depth, cold
pgls_fatalcs_depth <- samp_fatalcs_phylo %>%
  group_by(eid, sp, depth_col, temp_col) %>%
  summarise(frac_molar = sum(frac_molar)) %>%
  filter(temp_col <= cold) %>%
  pgls_ou(., tree_rooted, frac_molar ~ depth_col, branch_mult = 1E10)
summary(pgls_fatalcs_depth)$tTable

pgls_grid <- function(fit){
  ind <- fit$call$model %>% as.character %>% .[[3]]
  dep <- fit$call$model %>% as.character %>% .[[2]]
  rng <- fit$call$data[[ind]] %>% range()
  pre <- predict(fit, tibble(rng) %>% set_names(ind))
  tibble(
    ind = rng,
    dep = pre
  ) %>%
    set_names(c(ind, dep))
}

## FIGURE S2: FA composition and intercorrelations
pdf(file = "~/Documents/MBARI/Lipids/cteno-lipids-2020/figs/S2/Fig2_var15.pdf", width = 185/25.4, height = 135/25.4)
ggplot() +
  facet_grid(
    cols=vars(id.x),
    rows=vars(id.y),
    #scales="free"#,
    #switch="both"
  ) +
  # two-state significance in background
  geom_rect(data=pars_models %>% mutate(sig = (p<= 0.05)),
            xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
            aes(fill = sig),
            alpha=0.1
  ) +
  # continuous significance in background
  #geom_rect(data=pars_models %>% mutate(sig = (p<= 0.05)),
  #  xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
  #  aes(alpha = 1-p),
  #  color = "#E5E5E5"
  #) +
  geom_point(
    data=pair_frac_molar %>%
      mutate(
        sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
        sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
      ),
    aes(
      x=frac_molar.x,
      y=ifelse(paste(id.y) != "distribution", frac_molar.y, NA),
      #y=ifelse(id.x != id.y, frac_molar.y, NA),
      color=sprep
    ),
    alpha=0.2) +
  # regression lines
  geom_abline(
    data=pars_models,
    aes(
      slope=b,
      intercept=c#,
      #alpha=1-p
    )
  ) +
  scale_fill_manual(values = c(NA, "black")) +
  scale_color_manual(
    values = chroma_sp,
    labels = Vectorize(function(sprep){short2full[sprep]}),
    guide = guide_legend(
      direction = "vertical",
      title.position = "top"
    )
  ) +
  # marginal/diagonal boxplots
  geom_boxplot(
    data=pair_frac_molar,
    aes(
      x = mean(c(min(frac_molar.x), max(frac_molar.x))),
      y = ifelse(paste(id.y) == "distribution", frac_molar.y, NA),
    ),
    fill="grey",
    alpha=0.2,
    width=0.4
  ) +
  theme_pubr() +
  theme(
    legend.position = c(0.8, 0.7),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(alpha=FALSE, fill=FALSE) +
  labs(title="Fatty acid intercorrelations and composition", x="mole fraction", y="mole fraction", color="Species")
dev.off()

cold <- 7.5 # how cold is cold?
shal <- 200 # how shallow is shallow?
w <- 4
h <- 4
## DBI by temp
pdf(file = "~/Documents/MBARI/2020_eDSBS/20200925_dbi_vs_temp_pgls.pdf", width = w, height = h)
samp_fastats_phylo %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=dbi, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(data = pgls_grid(pgls_dbi_temp), color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 2.5) +
  labs(x = "Temperature (deg C)", y = "Mean # double bonds") +
  theme(legend.position = "none")
dev.off()

pdf(file = "~/Documents/MBARI/2020_eDSBS/20200925_dbi_vs_depth_pgls.pdf", width = w, height = h)
samp_fastats_phylo %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=depth_col, y=dbi, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(data = pgls_grid(pgls_dbi_depth), color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 2.5) +
  labs(x = "Depth (m)", y = "Mean # double bonds") +
  theme(legend.position = "none")
dev.off()

## chain length by temp
pdf(file = "~/Documents/MBARI/2020_eDSBS/20200925_chain_vs_temp_pgls.pdf", width = w, height = h)
samp_fastats_phylo %>%
  filter(depth_col <= shal) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=temp_col, y=chain, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(data = pgls_grid(pgls_chain_temp), color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 20) +
  labs(x = "Temperature (deg C)", y = "Mean chain length") +
  theme(legend.position = "none")
dev.off()

## by depth
pdf(file = "~/Documents/MBARI/2020_eDSBS/20200925_chain_vs_depth_pgls.pdf", width = w, height = h)
samp_fastats_phylo %>%
  filter(temp_col <= cold) %>%
  group_by(sp) %>%
  mutate(
    count = n(),
    sprep = ifelse(sp %in% names(chroma_sp), sp, "other"),
    sprep = factor(sprep, levels = c("Lamp_crue", "Bero_cucu", "Boli_infu", "Boli_vitr", "Leuc_pulc", "Bath_fost", "Pleu_bach", "Tjal_pink", "other"))
  ) %>%
  ggplot(aes(x=depth_col, y=chain, color=sprep)) +
  geom_point(alpha = 0.4, size = 4) +
  geom_line(data = pgls_grid(pgls_chain_depth), color = "black") +
  theme_pubr() +
  scale_color_manual(values = chroma_sp) +
  ylim(NA, 20) +
  labs(x = "Depth (m)", y = "Mean chain length") +
  theme(legend.position = "none")
dev.off()
