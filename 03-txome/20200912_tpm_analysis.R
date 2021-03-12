library(tidyverse)
library(ggpubr)

env_data <- read_tsv("/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/metadata/20200912_Cteno_depth_temp_EST.tsv")

dir_tpms <- "/Users/jwinnikoff/Documents/MBARI/Lipids/cteno-lipids-2020/kallisto"

key_annot = c(
  "OG0001664.tsv" = "ELOV6",
  "OG0004874.tsv" = "ELOV2",
  "OG0006673.tsv" = "SCD1",
  "OG0009524.tsv" = "ELOV3/6"
)

depth_3state = c(
  "CTE_Aulacoctena_acuminata"	= 2,
  "CTE_Bathyctena_chuni"	= 3,
  "CTE_Bathocyroe_fosteri"	= 2,
  "CTE_Benthic_cteno"	= 3,
  "CTE_Beroe_abyssicola" = 2,
  "CTE_Beroe_cucumis" = 2,
  "CTE_Beroe_forskalii" = 1,
  "CTE_Beroe_ovata" = 1,
  "CTE_Bolinopsis_infundibulum" = 2,
  "CTE_Bolinopsis_vitrea" = 1,
  "CTE_Cestum_veneris" = 1,
  "CTE_Charistephane_fugiens"	= 2,
  "CTE_Ctenoceros_spclear"	= 2,
  "CTE_Cydippid_spblack"	= 3,
  "CTE_Cydippid_sppeach"	= 3,
  "CTE_Cydippid_spredx"	= 2,
  "CTE_Deiopea_kaloktenota" = 1,
  "CTE_Dryodora_glandiformis"	= 1,
  "CTE_Euplokamis_dunlapae"	= 1,
  "CTE_Haeckelia_beehleri" = 1,
  "CTE_Haeckelia_rubra" = 1,
  "CTE_Hormiphora_californensis" = 1,
  "CTE_Kiyohimea_sp"	= 2,
  "CTE_Lampocteis_cruentiventer"	= 2,
  "CTE_Lampea_deep"	= 3,
  "CTE_Lampea_lactea"	= 1,
  "CTE_Leucothea_pulchra" = 1,
  "CTE_Llyria_spbenthic"	= 3,
  "CTE_Llyria_spcopper"	= 3,
  "CTE_Llyria_spdeep"	= 3,
  "CTE_Mertensia_ovum" = 1,
  "CTE_Nepheloctena_red " = 2,
  "CTE_Nepheloctena_whit"	= 2,
  "CTE_Ocyropsis_crystallina" = 1,
  "CTE_Ocyropsis_maculata" = 1,
  "CTE_Tetraphalia_sp"	= 1,
  "CTE_Thalassocalyce_inconstans" = 2,
  "CTE_Velamen_parallelum" = 1,
  "CTE_Vermillion_lobate"	= 2,
  "CTE_Weird_cteno"	= 2
)

# read TPMs
readdata <- list.files(path = dir_tpms, full.names = T) %>%
  lapply(
    function(file_data){
      data_this_file <- file_data %>%
        read_tsv(col_names=FALSE) %>%
        mutate(file = file_data %>% basename())
      return(data_this_file)
    }
  ) %>%
  do.call(rbind, .) %>%
  magrittr::set_colnames(c("target_id", "length", "eff_length", "est_counts", "tpm", "file_data"))

reads_top <- readdata %>%
  rowwise() %>%
  mutate(species = strsplit(target_id, "\\|") %>% unlist() %>% .[1] %>% gsub("[0-9]", "", x=.)) %>%
  group_by(species, file_data) %>%
  filter(tpm == max(tpm)) %>%
  left_join(env_data, by="species")

# max by species
reads_top %>%
  filter(temp_med <= 7.5) %>%
  ungroup() %>%
  ggplot(aes(x=species, y=tpm)) +
    facet_wrap(~file_data, ncol=4) +
    geom_col() +
    #geom_text(aes(label=species)) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

reads_sum <- readdata %>%
  rowwise() %>%
  mutate(
    sample = strsplit(target_id, "\\|") %>% unlist() %>% .[1],
    species = sample %>% gsub("[0-9]", "", x=.)
  ) %>%
  ungroup() %>%
  complete(sample, file_data, fill=list("tpm"=0)) %>%
  group_by(sample, species, file_data) %>%
  summarize(
    tpm_tot = sum(tpm),
    tpm_avg = mean(tpm)
  ) %>%
  rowwise() %>% 
  mutate(depthclass = depth_3state[species])
  #left_join(env_data, by="species") %>%
  #left_join(read_tsv("/Users/jwinnikoff/Documents/MBARI/lab-work/converge/datasets/ctenoPK/trait/Cteno_depth_temp_desc_20180529.txt"), by="sp")

# plot total by species
reads_sum %>%
  filter(file_data == "OG0006673.tsv") %>%
  #filter(temp_med <= 7.5) %>%
  ungroup() %>%
  ggplot(aes(x=sample, y=tpm_tot)) +
  #facet_wrap(~file_data, ncol=4) +
  geom_col() +
  #geom_text(aes(label=species)) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x="sample", y="transcripts per million, sum of all isoforms")

barplots <- reads_sum %>%
  filter(species != "CTE_Euplokamis_dunlapae") %>%
  filter(species != "CTE_Vermillion_lobate") %>%
  filter(file_data != "OG0009524.tsv") %>%
  mutate(annot = key_annot[file_data]) %>%
  group_by(annot) %>%
  group_split() %>%
  lapply(
    .,
    function(data){
      gg <- data %>%
        ggplot(aes(x=sample, y=tpm_tot, fill=desc(depthclass))) +
          geom_col() +
          theme_pubr() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none"
          ) +
          labs(
            title = data %>% pull(annot) %>% first(),
            x = "sample",
            y = "transcripts per million, sum of all isoforms"
          )
      return(gg)
    }
  )

barplots_pretty <- barplots %>%
  .[1:2] %>%
  lapply(
    .,
    function(gg){gg + theme(axis.text.x = element_blank()) + xlab(NULL)}
  ) %>%
  c(., barplots %>% .[3])

barplots_pretty %>%
  grid.arrange(grobs=., ncol=1, heights=c(1,1,2))

