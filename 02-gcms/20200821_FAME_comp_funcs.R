## mappings, not technically functions:

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

# this function takes a phylo object and a table of phenotypic data with column "sp"
# It dereplicates duplicate tip labels, then makes polytomies to match the number of indl's/sp in pheno.
# Finally, it runs PGLS with an OU covariance model, and returns the fit object
pgls_ou <- function(pheno, phylo, form, startval = 1, autoroot = TRUE, branch_mult=1E15){
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
    anti_join(leaves_dup %>% summarize(node = first(node)), by = c("sp", "node"))

  # set first conspecific branch to average length
  branches_adj <- bind_rows(anti_join(branches, taxa_dup, by = c("node", "mrca")), taxa_dup)

  # replace conspecific clades with a single tip
  phylo_pruned <- phylo %>%
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

  # add conspecifics as polytomies
  leaves_add <- pheno %>%
    group_by(sp) %>%
    summarize(add = n()-1) %>%
    filter(add > 0) %>%
    rowwise()

  # prevent crash if there are no duplicate taxa in pheno
  if(nrow(leaves_add) > 0){
    leaves_add <- leaves_add %>%
      summarize(spadd = rep(sp, add) %>% list()) %>%
      pull(spadd) %>%
      unlist()
  }else{
    leaves_add <- list()
  }

  phylo_expanded <- phylo_pruned
  #edge_polytomy <- min(phylo_pruned$edge.length)/100
  edge_polytomy <- 0
  for(sp in leaves_add){
    phylo_expanded <- phylo_expanded %>%
      bind.tip(tip.label = sp, where = match(sp, phylo_expanded$tip.label), edge.length = edge_polytomy)
  }

  # number the individuals in the tree_rooted and the table
  phylo_expanded$tip.label <- tibble(sp = phylo_expanded$tip.label) %>%
    group_by(sp) %>%
    mutate(indl = paste(sp, row_number(), sep="")) %>%
    pull(indl)

  # make the individual names rownames in the data
  pheno_rownames <- pheno %>%
    group_by(sp) %>%
    mutate(indl = paste(sp, row_number(), sep="")) %>%
    column_to_rownames(var = "indl")

  # root the tree_rooted if not rooted
  if(autoroot && !is.rooted(phylo_expanded)){
    message("midpoint-rooting phylogeny")
    phylo_expanded <- phylo_expanded %>% midpoint.root()
  }

  #plot(phylo_expanded) #TEST

  # rescale branches to avoid singular var-cov matrix
  phylo_expanded$edge.length <- phylo_expanded$edge.length * branch_mult

  # do.call is needed to preserve the formula call in the object
  fit <- do.call(
    nlme::gls,
    list(
      data = pheno_rownames,
      model = as.formula(form),
      # Brownian model with lambda estimated concurrently
      #correlation = corPagel(1, phy = phylo_expanded, fixed = FALSE),
      # an OU model, Brownian is garbage
      correlation = corMartins(startval, phy = phylo_expanded, fixed = FALSE),
      method = "ML"#,
      #overrides!
      #control = list(
      #  singular.ok = TRUE,
      #  msTol = 1E-31,
      #  tolerance = 1E-30
      #)
    )
  )

  # store independent variable limits for later extraction
  #fit$limits <- range(pheno[[as.character(form)[[3]]]], na.rm = TRUE)

  return(fit)
}

# run PGLS-OU with a series of different branch multipliers
pgls_opt <- function(pheno, phylo, form, branch_mults, startval = 1, autoroot = TRUE){
  fits <- lapply(
    branch_mults,
    function(mult){
      try(pgls_ou(pheno, phylo, form, startval=startval, autoroot=autoroot, branch_mult = mult), silent=TRUE)
    }
  ) %>%
    setNames(branch_mults)
  # acrobatic, but it works
  fits_tbl <- fits %>%
    # throw out the errors
    tibble() %>%
    rowwise() %>%
    filter(class(`.`) != "try-error")

  fits_names <- fits %>%
    # throw out the errors
    tibble() %>%
    mutate(branch_mult = branch_mults) %>%
    rowwise() %>%
    filter(class(`.`) != "try-error")

  fits_tbl %>%
    pull(`.`) %>%
    setNames(fits_names$branch_mult)
}

# run pairwise PGLSs, return the fits as a named list
pgls_pairwise <- function(pheno, phylo){
  # make id an ordered factor
  pheno <- pheno %>%
    arrange(chain) %>%
    mutate(id = make.names(id, unique=TRUE))
  pheno$id <- factor(pheno$id, levels=pheno %>% pull(id) %>% unique())
  pheno_wide <- pheno %>%
    pivot_wider(.,
                id_cols = c("eid", "sp"),
                names_from = id,
                values_from = frac_molar
    ) %>%
    drop_na()
  # generate fit for each pair of variables
  lapply(
    levels(pheno$id),
    function(x){
      lapply(
        levels(pheno$id)[which(levels(pheno$id) > x)],
        function(y){
          print(paste(y, x, sep='~'))
          pgls_ou(
            pheno_wide,
            phylo,
            as.formula(paste(y, x, sep='~')),
            branch_mult = 1E10
          )
        }
      )
    }
  ) %>%
    # make model spec the name of the fit object
    unlist(recursive = FALSE) %>%
    setNames(lapply(., function(i){i$call$model %>% c() %>% paste()}))
}

# extract slope, intercept, p-val and model spec from fit object
get_params <- function(fit){
  summ <- summary(fit)
  ttab <- summ %>%
    .$tTable %>%
    as_tibble()
  tibble(
    x = fit$call$model[[3]] %>% paste(),
    y = fit$call$model[[2]] %>% paste(),
    alpha = fit$modelStruct$corStruct %>% as.numeric(),
    c = ttab$Value[[1]],
    b = ttab$Value[[2]],
    p = ttab$`p-value`[[2]],
    logLik = fit$logLik,
    AIC = summ$AIC,
    BIC = summ$BIC
  )
}

# Same as above, but accommodates multiple variables using corphylo
corphylo_indls <- function(pheno, phylo, autoroot = TRUE){
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
    anti_join(leaves_dup %>% summarize(node = first(node)), by = c("sp", "node"))

  # set first conspecific branch to average length
  branches_adj <- bind_rows(anti_join(branches, taxa_dup, by = c("node", "mrca")), taxa_dup)

  # replace conspecific clades with a single tip
  phylo_pruned <- phylo %>%
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

  # add conspecifics as polytomies
  leaves_add <- pheno %>%
    group_by(sp) %>%
    summarize(add = n()-1) %>%
    filter(add > 0) %>%
    rowwise() %>%
    summarize(spadd = rep(sp, add) %>% list()) %>%
    pull(spadd) %>%
    unlist()

  phylo_expanded <- phylo_pruned
  # throws warnings!
  for(sp in leaves_add){
    phylo_expanded <- phylo_expanded %>%
      bind.tip(tip.label = sp, where = match(sp, phylo_expanded$tip.label), edge.length = 0)
  }

  # number the individuals in the tree_rooted and the table
  phylo_expanded$tip.label <- tibble(sp = phylo_expanded$tip.label) %>%
    group_by(sp) %>%
    mutate(indl = paste(sp, row_number(), sep="")) %>%
    pull(indl)

  # make the individual names rownames in the data
  pheno_rownames <- pheno %>%
    group_by(sp) %>%
    mutate(indl = paste(sp, row_number(), sep="")) %>%
    ungroup() %>%
    select(-sp) %>%
    column_to_rownames(var = "indl") %>%
    as.matrix()

  # root the tree if not rooted
  if(autoroot && !is.rooted(phylo_expanded)){
    message("midpoint-rooting phylogeny")
    phylo_expanded <- phylo_expanded %>% midpoint.root()
  }

  # rescale branches to avoid singular var-cov matrix
  #phylo_expanded$edge.length <- phylo_expanded$edge.length * 100

  ape::corphylo(
    X = pheno_rownames,
    # matrix format of the elements of U is essential!
    U = pheno_rownames %>% asplit(MARGIN=2) %>% lapply(as.matrix),
    phy = phylo_expanded,
    constrain.d = TRUE
  )
}
