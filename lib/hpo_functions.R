hpoa_path <- "./data/phenotype.20250303.hpoa"
gene_to_phenotype_path <- "./data/genes_to_phenotype.20250303.txt"
pheno_to_gene_path <- "./data/phenotype_to_genes.20250303.txt"

hpo <- ontologyIndex::get_OBO("./data/hp.20250303.obo")

get_all_hpo_terms <- function() {
  message("getting all hpo terms..")
  all_hpos <- get_descendants(hpo, "HP:0000001")
  all_hpos_terms <- unlist(purrr::map(all_hpos, ~ get_term_property(hpo, property_name = "name", term = .)))
  all_hpos_terms <- data.frame(hpo_id = names(all_hpos_terms), hpo_name = as.character(all_hpos_terms)) %>%
    rowwise() %>%
    mutate(hpo_id_name = glue("{hpo_id}:{hpo_name}"))
  return(all_hpos_terms)
}

create_disorder_to_gene_table <- function() {
  # load pheno to gene
  pheno_to_genes_df <- readr::read_delim(pheno_to_gene_path, delim = "\t")
  # add disorder name
  hpoa_df <- readr::read_delim(hpoa_path, delim = "\t", col_types = paste0(rep("c", 12), collapse = ""), skip = 4)
  hpoa_df <- hpoa_df %>% rename("disease_id" = "database_id")
  disorder_to_gene <- pheno_to_genes_df %>%
    select(disease_id, gene_symbol) %>%
    distinct() %>%
    left_join(hpoa_df %>% select(disease_id, disease_name) %>% distinct())
  return(disorder_to_gene)
}

# from https://hpo.jax.org/app/browse/term/HP:0040279 take the mean freq for each value
hpo_frequencies <- list(
  "HP:0040284" = 0.025,
  "HP:0040282" = 0.55,
  "HP:0040283" = 0.17,
  "HP:0040280" = 1,
  "HP:0040281" = 0.9,
  "HP:0040285" = 0
)

# hpo_frequencies<-data.frame(hpo_frequencies)


parse_ratio_column <- function(df, column_name, min_num_o_patients_for_freq = 10) {
  # Extract the specified column
  ratio_strings <- df[[column_name]]

  # Create vectors to store parsed values
  x_values <- c()
  y_values <- c()
  ratio_values <- c()

  for (i in seq_along(ratio_strings)) {
    # print(ratio_strings[i])
    # Check if the string is NA or doesn't match the expected format
    if (is.na(ratio_strings[i]) || !grepl("^\\d+/\\d+$", ratio_strings[i])) {
      x_values[i] <- y_values[i] <- ratio_values[i] <- NA
    } else {
      # Split the string and parse the values
      split_values <- strsplit(ratio_strings[i], "/")[[1]]
      x_values[i] <- as.numeric(split_values[1])
      y_values[i] <- as.numeric(split_values[2])
      # if the numebr of patients is less than min_num_o_patients_for_freq change the frequency to be NA
      ratio_values[i] <- x_values[i] / y_values[i]
      if (y_values[i] <= min_num_o_patients_for_freq) {
        ratio_values[i] <- NA
      }
    }
  }

  # Create a new data frame with the original column and three parsed columns
  df$num_o_patients_with_pheno <- x_values
  df$num_o_patients <- y_values
  df$frequency_numeric <- ratio_values

  return(df)
}

parse_hpo_hpoa_db <- function() {
  # read the pheno to gene hpo table
  pheno_to_genes_df <- readr::read_delim(pheno_to_gene_path, delim = "\t")
  all_diseases <- pheno_to_genes_df %>%
    select(disease_id) %>%
    distinct()
  all_phenos_with_associated_genes <- pheno_to_genes_df %>%
    select(hpo_id, hpo_name) %>%
    mutate(hpo_id_name = glue("{hpo_id}:{hpo_name}")) %>%
    distinct()
  # read the hpoa table (also from hpo)
  hpoa_df <- readr::read_delim(hpoa_path, delim = "\t", col_types = paste0(rep("c", 12), collapse = ""), skip = 4)
  hpoa_df <- hpoa_df %>% rename("disease_id" = "database_id")
  # now only keep disorders that are found in the pheno to gene file
  message("removing hpos and disease without an associated gene..")
  hpoa_df <- hpoa_df %>%
    inner_join(all_diseases) %>%
    inner_join(all_phenos_with_associated_genes)
  # Add numeric frequency term
  # First parse the X/X Frequcny
  message("parsing frequency column..")
  hpoa_df <- parse_ratio_column(hpoa_df, "frequency")
  # Then parse the HPOs based frequencies
  hpoa_df <- hpoa_df %>%
    rowwise() %>%
    mutate(frequency_numeric = ifelse(frequency %in% names(hpo_frequencies),
      hpo_frequencies[[frequency]],
      frequency_numeric
    ))
  # finally, convert the percentage to frequency_numeric
  hpoa_df <- hpoa_df %>%
    rowwise() %>%
    mutate(frequency_numeric = ifelse(grepl("%", frequency),
      as.numeric(str_extract(frequency, "\\d+(.\\d)*")) / 100,
      frequency_numeric
    ))


  # Add categorical frequency term
  message("adding categorical frequency column..")
  hpoa_df <- hpoa_df %>%
    mutate(frequency_cat = case_when(
      !is.na(frequency_numeric) & frequency_numeric <= 0.049 & frequency_numeric >= 0.001 ~ "very_rare",
      !is.na(frequency_numeric) & frequency_numeric == 0 ~ "excluded",
      !is.na(frequency_numeric) & frequency_numeric <= 0.799 & frequency_numeric >= 0.3 ~ "frequent",
      !is.na(frequency_numeric) & frequency_numeric <= 0.999 & frequency_numeric >= 0.8 ~ "very_frequent",
      !is.na(frequency_numeric) & frequency_numeric <= 0.299 & frequency_numeric >= 0.05 ~ "occasional",
      !is.na(frequency_numeric) & frequency_numeric == 1 ~ "obligate",
      is.na(frequency_numeric) ~ "unknown"
    ))
  hpoa_df %>% count(frequency_cat)
  # now change phenotypes that have very low number of observations that are obligate into frequent
  hpoa_df <- hpoa_df %>% mutate(frequency_cat = ifelse(!is.na(num_o_patients) & num_o_patients == 1 & frequency_cat == "obligate",
    "frequent",
    frequency_cat
  ))
  hpoa_df <- hpoa_df %>% mutate(frequency_cat = ifelse(is.na(frequency_cat), "unknown", frequency_cat))
  # Add score for frequency
  hpoa_df <- hpoa_df %>%
    mutate(
      frequency_score = case_when(
        frequency_cat == "very_rare" ~ 0.025,
        frequency_cat == "excluded" ~ -5,
        frequency_cat == "frequent" ~ 0.55,
        frequency_cat == "very_frequent" ~ 0.9,
        frequency_cat == "occasional" ~ 0.17,
        frequency_cat == "obligate" ~ 1,
        frequency_cat == "unknown" ~ 0.55
      ),
      frequency_bayes = case_when(
        frequency_cat == "very_rare" ~ 0.025,
        frequency_cat == "excluded" ~ 0.001,
        frequency_cat == "frequent" ~ 0.55,
        frequency_cat == "very_frequent" ~ 0.9,
        frequency_cat == "occasional" ~ 0.17,
        frequency_cat == "obligate" ~ 0.999,
        frequency_cat == "unknown" ~ 0.2
      )
    )


  # remove duplicates, for each disease (disease_id) take the one with the highest frequency
  # !! TODO !! need to consider if this is the best way
  message("removing duplicate disease_id+hpo_ids - keeping the row with the highest frequency..")
  hpoa_df <- hpoa_df %>%
    group_by(disease_id, hpo_id_name) %>%
    slice_max(frequency_score, n = 1, with_ties = F) %>%
    ungroup()

  # A table containing all the ancestors for each phenotype
  message("getting the ancestors of each HPO..")
  all_hpos_terms <- get_all_hpo_terms()
  all_ancestors <- hpoa_df %>%
    select(hpo_id) %>%
    distinct() %>%
    rowwise() %>%
    mutate(
      ancestors = paste0(get_hpo_ancestors_by_id(hpo_id, with_term = F), collapse = ","),
      num_ancestors = get_num_of_hpo_ancestors_by_id(hpo_id)
    )

  # A table containing all the phenotypes ids and terms
  all_phenos <- hpoa_df %>%
    select(hpo_id, hpo_name, hpo_id_name) %>%
    distinct()
  # join the hpo table with the descendants and then separate the descendants to different rows
  # hpo_specificity_df<-readr::read_delim('./data/hpo_specificity.csv')
  # hpoa_df<-hpoa_df%>%left_join(hpo_specificity_df) # join with specificity before adding descendants
  message("expanding table to include all ancestors..")
  hpoa_df <- hpoa_df %>%
    left_join(all_ancestors) %>%
    rowwise() %>%
    separate_rows(ancestors, sep = ",") %>%
    mutate(
      is_ancestor = ifelse(hpo_id == ancestors, F, T), # set it so that if the ancestor id is the same as the original id - it is not an ancestor
      ancestor_of = hpo_id_name,
      hpo_id = ancestors
    ) %>%
    select(-c(hpo_id_name, hpo_name)) %>% # remove the term and id-term and repopulate them according to the new row HPO ID
    left_join(all_hpos_terms)

  hpo_specificity_df <- generate_hpo_specificity_table(hpoa_df)

  hpoa_df <- hpoa_df %>%
    left_join(hpo_specificity_df) %>%
    mutate(frequency_score_with_specificity = frequency_score / specificity)

  # add specificity score to the hpoa_df - the more disorders a phenotype is associated with the less specific it is and the score is lower
  # you can use this score to modify the unknown frequency score

  hpoa_df <- hpoa_df %>% left_join(
    hpoa_df %>%
      count(hpo_id_name, disease_id) %>%
      count(hpo_id_name) %>%
      rename("num_of_diseases" = "n")
  )
  hpoa_df <- hpoa_df %>% mutate(specificity_score = as.numeric(as.character(cut(num_of_diseases, breaks = c(0, 10, 50, 100, 500, 1000, 9999), labels = c(1, 0.8, 0.6, 0.3, 0.1, 0)))))

  return(hpoa_df)
}


get_disease_list <- function(hpoa_df) {
  disease_list <- hpoa_df %>%
    group_by(disease_id, DiseaseName) %>%
    summarize(hpos = list(hpo_id_name)) %>%
    ungroup()
  return(disease_list)
}

# retrieve all the descendants of a given ID (in oppose to just the nearest children)
get_hpo_descendants_by_id <- function(hpo_id, with_term = T) {
  descendants <- ontologyIndex::get_descendants(hpo, hpo_id)
  if (length(descendants) == 0) {
    message(glue("HPO ID {hpo_id} is not found..skipping.."))
    return(NA)
  }
  get_desc_name <- function(hpo_id) {
    hpo_name <- ontologyIndex::get_term_property(ontology = hpo, property = "name", term = hpo_id)
    return(glue("{hpo_id}:{hpo_name}"))
  }
  if (with_term) {
    descendants_terms <- purrr::map_vec(descendants, get_desc_name)
    return(descendants_terms)
  } else {
    return(descendants)
  }
}


get_hpo_ancestors_by_id <- function(hpo_id, with_term = T) {
  # print(hpo_id)
  ancestors <- ontologyIndex::get_ancestors(hpo, hpo_id)
  # if the hpo_id is not found - return NA
  if (length(ancestors) == 0) {
    message(glue("HPO ID {hpo_id} is not found..skipping.."))
    return(NA)
  }
  ancestors_terms <- ontologyIndex::get_term_property(ontology = hpo, property = "ancestors", term = hpo_id, as_names = TRUE)
  ancestors_terms <- glue("{names(ancestors_terms)}:{ancestors_terms}")
  if (with_term) {
    return(ancestors_terms)
  } else {
    return(ancestors)
  }
}

get_num_of_hpo_ancestors_by_id <- function(hpo_id) {
  hpo_desendants <- ontologyIndex::get_ancestors(hpo, hpo_id)
  return(length(hpo_desendants) - 1)
}

get_num_of_hpo_children_by_id <- function(hpo_id) {
  hpo_desendants <- ontologyIndex::get_descendants(hpo, hpo_id)
  return(length(hpo_desendants) - 1)
}



# !!!! BEFORE you run this script, make sure the updated version of the hpoa_df is loaded, and only then create the specificity
# table, otherwise it will cause discrepancy
generate_hpo_specificity_table <- function(hpoa_df) {
  hpo_specificity_df <- hpoa_df %>%
    # filter(!is_descendant)%>%
    group_by(hpo_id_name, hpo_id, hpo_name) %>%
    filter(!frequency_cat == "excluded") %>%
    summarize(specificity = sum(frequency_score))
  # summarize(specificity=max(1,sum(frequency_score)))# the minimal specificity should be 1 (otherwise if the frequency is rare and it only occurs once it is considered more specific than an obligate)
  write.table(hpo_specificity_df, file = "./data/hpo_specificity.csv", sep = "\t", row.names = F, quote = F)
  return(hpo_specificity_df)
}

depracated_parse_hpo_hpoa_db <- function() {
  genes_to_phenotype_df <- readr::read_delim(gene_to_phenotype_path, delim = "\t")
  hpoa_df <- readr::read_delim(hpoa_path, delim = "\t", col_types = paste0(rep("c", 12), collapse = ""), skip = 4)
  hpoa_df <- hpoa_df %>% rename("disease_id" = "#disease_id")
  # add hpo term
  hpoa_df <- hpoa_df %>%
    left_join(genes_to_phenotype_df %>% select(hpo_id = "HPO-Term-ID", hpo_name = "HPO-Term-Name") %>% distinct()) %>%
    mutate(hpo_id_name = as.character(glue("{hpo_id}:{hpo_name}")))
  terms_to_fix <- hpoa_df %>%
    filter(is.na(hpo_name)) %>%
    pull(hpo_id)
  terms_to_fix <- terms_to_fix[unlist(purrr::map(terms_to_fix, ~ length(ontologyIndex::get_term_frequencies(hpo, .x)) > 0))]
  fixed_hpoa_terms <- hpoa_df %>%
    filter(hpo_id %in% terms_to_fix) %>%
    rowwise() %>%
    mutate(
      hpo_name = ontologyIndex::get_term_property(hpo, "name", hpo_id),
      hpo_id_name = as.character(glue("{hpo_id}:{hpo_name}"))
    )

  # now join the tables
  hpoa_df <- hpoa_df %>%
    filter(!is.na(hpo_name)) %>%
    bind_rows(fixed_hpoa_terms)

  # Add numeric frequency term
  parse_ratio_column <- function(df, column_name) {
    # Extract the specified column
    ratio_strings <- df[[column_name]]

    # Create vectors to store parsed values
    x_values <- c()
    y_values <- c()
    ratio_values <- c()

    for (i in seq_along(ratio_strings)) {
      # print(ratio_strings[i])
      # Check if the string is NA or doesn't match the expected format
      if (is.na(ratio_strings[i]) || !grepl("^\\d+/\\d+$", ratio_strings[i])) {
        x_values[i] <- y_values[i] <- ratio_values[i] <- NA
      } else {
        # Split the string and parse the values
        split_values <- strsplit(ratio_strings[i], "/")[[1]]
        x_values[i] <- as.numeric(split_values[1])
        y_values[i] <- as.numeric(split_values[2])
        ratio_values[i] <- x_values[i] / y_values[i]
      }
    }

    # Create a new data frame with the original column and three parsed columns
    df$num_o_patients_with_pheno <- x_values
    df$num_o_patients <- y_values
    df$frequency_numeric <- ratio_values

    return(df)
  }

  # First parse the X/X Frequcny
  hpoa_df <- parse_ratio_column(hpoa_df, "frequency")
  # Then parse the HPOs based frequencies
  hpoa_df <- hpoa_df %>%
    rowwise() %>%
    mutate(frequency_numeric = ifelse(frequency %in% names(hpo_frequencies),
      hpo_frequencies[[frequency]],
      frequency_numeric
    ))
  # finally, convert the percentage to frequency_numeric
  hpoa_df <- hpoa_df %>%
    rowwise() %>%
    mutate(frequency_numeric = ifelse(grepl("%", frequency),
      as.numeric(str_extract(frequency, "\\d+(.\\d)*")) / 100,
      frequency_numeric
    ))


  # Add categorical frequency term
  hpoa_df <- hpoa_df %>%
    mutate(frequency_cat = case_when(
      !is.na(frequency_numeric) & frequency_numeric <= 0.049 & frequency_numeric >= 0.001 ~ "very_rare",
      !is.na(frequency_numeric) & frequency_numeric == 0 ~ "excluded",
      !is.na(frequency_numeric) & frequency_numeric <= 0.799 & frequency_numeric >= 0.3 ~ "frequent",
      !is.na(frequency_numeric) & frequency_numeric <= 0.999 & frequency_numeric >= 0.8 ~ "very_frequent",
      !is.na(frequency_numeric) & frequency_numeric <= 0.299 & frequency_numeric >= 0.05 ~ "occasional",
      !is.na(frequency_numeric) & frequency_numeric == 1 ~ "obligate",
      is.na(frequency_numeric) ~ "unknown"
    ))
  hpoa_df %>% count(frequency_cat)
  # now change phenotypes that have very low number of observations that are obligate into frequent
  hpoa_df <- hpoa_df %>% mutate(frequency_cat = ifelse(!is.na(num_o_patients) & num_o_patients == 1 & frequency_cat == "obligate",
    "frequent",
    frequency_cat
  ))

  # Add score for frequency
  hpoa_df <- hpoa_df %>%
    mutate(
      frequency_score = case_when(
        frequency_cat == "very_rare" ~ 0.01,
        frequency_cat == "excluded" ~ -5,
        frequency_cat == "frequent" ~ 0.8,
        frequency_cat == "very_frequent" ~ 0.9,
        frequency_cat == "occasional" ~ 0.3,
        frequency_cat == "obligate" ~ 1,
        frequency_cat == "unknown" ~ 0.1
      )
    )
  # remove duplicates, for each disease (disease_id) take the one with the highest frequency
  # !! TODO !! need to consider if this is the best way
  hpoa_df <- hpoa_df %>%
    group_by(disease_id, hpo_id_name) %>%
    slice_max(frequency_score, n = 1, with_ties = F) %>%
    ungroup()

  # A table containing all the ancestors for each phenotype
  all_ancestors <- hpoa_df %>%
    select(hpo_id) %>%
    distinct() %>%
    rowwise() %>%
    mutate(
      ancestors = paste0(get_hpo_ancestors_by_id(hpo_id, with_term = F), collapse = ","),
      num_ancestors = get_num_of_hpo_ancestors_by_id(hpo_id)
    )

  # A table containing all the phenotypes ids and terms
  all_phenos <- hpoa_df %>%
    select(hpo_id, hpo_name, hpo_id_name) %>%
    distinct()
  all_hpos_terms <- get_all_hpo_terms()
  # join the hpo table with the descendants and then separate the descendants to different rows
  # hpo_specificity_df<-readr::read_delim('./data/hpo_specificity.csv')
  # hpoa_df<-hpoa_df%>%left_join(hpo_specificity_df) # join with specificity before adding descendants
  hpoa_df <- hpoa_df %>%
    left_join(all_ancestors) %>%
    rowwise() %>%
    separate_rows(ancestors, sep = ",") %>%
    mutate(
      is_ancestor = ifelse(hpo_id == ancestors, F, T), # set it so that if the ancestor id is the same as the original id - it is not an ancestor
      ancestor_of = hpo_id_name,
      hpo_id = ancestors
    ) %>%
    select(-c(hpo_id_name, hpo_name)) %>% # remove the term and id-term and repopulate them according to the new row HPO ID
    left_join(all_hpos_terms)

  hpo_specificity_df <- generate_hpo_specificity_table(hpoa_df)

  hpoa_df <- hpoa_df %>%
    left_join(hpo_specificity_df) %>%
    mutate(frequency_score_with_specificity = frequency_score / specificity)
  return(hpoa_df)
}
