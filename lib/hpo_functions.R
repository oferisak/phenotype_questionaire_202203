hpoa_path<-'./data/phenotype.20230301.hpoa'
gene_to_phenotype_path<-'./data/genes_to_phenotype.txt'
hpo<-ontologyIndex::get_OBO('./data/hp.20230301.obo')
all_hpos<-get_descendants(hpo,'HP:0000001')
all_hpos_terms<-unlist(purrr::map(all_hpos,~get_term_property(hpo,property_name = 'name',term = .)))
all_hpos_terms<-data.frame(HPO_ID=names(all_hpos_terms),HPO_TERM=as.character(all_hpos_terms))%>%
  rowwise()%>%mutate(HPO_ID_TERM=glue('{HPO_ID}:{HPO_TERM}'))

# from https://hpo.jax.org/app/browse/term/HP:0040279 take the mean freq for each value
hpo_frequencies<-list('HP:0040284'=0.025,
                      'HP:0040282'=0.55,
                      'HP:0040283'=0.17,
                      'HP:0040280'=1,
                      'HP:0040281'=0.9,
                      'HP:0040285'=0)

#hpo_frequencies<-data.frame(hpo_frequencies)

parse_hpo_hpoa_db<-function(){
  genes_to_phenotype_df<-readr::read_delim(gene_to_phenotype_path,delim='\t')  
  hpoa_df<-readr::read_delim(hpoa_path,delim='\t',col_types =paste0(rep('c',12),collapse=''),skip=4)
  hpoa_df<-hpoa_df%>%rename('DatabaseID'='#DatabaseID')
  # add hpo term
  hpoa_df<-hpoa_df%>%left_join(genes_to_phenotype_df%>%select(HPO_ID='HPO-Term-ID',HPO_TERM='HPO-Term-Name')%>%distinct())%>%
    mutate(HPO_ID_TERM=as.character(glue('{HPO_ID}:{HPO_TERM}')))
  terms_to_fix<-hpoa_df%>%filter(is.na(HPO_TERM))%>%pull(HPO_ID)
  terms_to_fix<-terms_to_fix[unlist(purrr::map(terms_to_fix,~length(ontologyIndex::get_term_frequencies(hpo,.x))>0))]
  fixed_hpoa_terms<-hpoa_df%>%
    filter(HPO_ID%in%terms_to_fix)%>%rowwise()%>%
    mutate(HPO_TERM=ontologyIndex::get_term_property(hpo,'name',HPO_ID),
           HPO_ID_TERM=as.character(glue('{HPO_ID}:{HPO_TERM}')))
  
  # now join the tables
  hpoa_df<-hpoa_df%>%filter(!is.na(HPO_TERM))%>%
    bind_rows(fixed_hpoa_terms)
  
  # Add numeric frequency term
  parse_ratio_column <- function(df, column_name) {
    # Extract the specified column
    ratio_strings <- df[[column_name]]
    
    # Create vectors to store parsed values
    x_values <-c() 
    y_values <-c() 
    ratio_values <- c()
    
    for (i in seq_along(ratio_strings)) {
      #print(ratio_strings[i])
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
    df$Frequency_numeric <- ratio_values
    
    return(df)
  }
  
  # First parse the X/X Frequcny
  hpoa_df<-parse_ratio_column(hpoa_df,'Frequency')
  # Then parse the HPOs based frequencies
  hpoa_df<-hpoa_df%>%rowwise()%>%
    mutate(Frequency_numeric=ifelse(Frequency%in%names(hpo_frequencies),
                                                 hpo_frequencies[[Frequency]],
                                                 Frequency_numeric))
  # finally, convert the percentage to Frequency_numeric
  hpoa_df<-hpoa_df%>%rowwise()%>%
    mutate(Frequency_numeric=ifelse(grepl('%',Frequency),
                                    as.numeric(str_extract(Frequency,'\\d+(.\\d)*'))/100,
                                    Frequency_numeric))
  
  
  # Add categorical frequency term
  hpoa_df<-hpoa_df%>%
    mutate(Frequency_cat=case_when(
      !is.na(Frequency_numeric) & Frequency_numeric<=0.049 & Frequency_numeric>=0.001 ~ 'very_rare',
      !is.na(Frequency_numeric) & Frequency_numeric==0 ~ 'excluded',
      !is.na(Frequency_numeric) & Frequency_numeric<=0.799 & Frequency_numeric>=0.3 ~ 'frequent',
      !is.na(Frequency_numeric) & Frequency_numeric<=0.999 & Frequency_numeric>=0.8 ~ 'very_frequent',
      !is.na(Frequency_numeric) & Frequency_numeric<=0.299 & Frequency_numeric>=0.05 ~ 'occasional',
      !is.na(Frequency_numeric) & Frequency_numeric==1 ~ 'obligate',
      is.na(Frequency_numeric) ~ 'unknown'
    ))
  hpoa_df%>%count(Frequency_cat)
  # now change phenotypes that have very low number of observations that are obligate into frequent
  hpoa_df<-hpoa_df%>%mutate(Frequency_cat=ifelse(!is.na(num_o_patients) & num_o_patients==1 & Frequency_cat=='obligate',
                                                 'frequent',
                                                 Frequency_cat))
  
  # Add score for frequency
  hpoa_df<-hpoa_df%>%
    mutate(Frequency_score=case_when(
      Frequency_cat=='very_rare'~0.01,
      Frequency_cat=='excluded'~-5,
      Frequency_cat=='frequent'~0.8,
      Frequency_cat=='very_frequent'~0.9,
      Frequency_cat=='occasional'~0.3,
      Frequency_cat=='obligate'~1,
      Frequency_cat=='unknown'~0.1
    ))
  # remove duplicates, for each disease (DatabaseID) take the one with the highest frequency 
  # !! TODO !! need to consider if this is the best way
  hpoa_df<-hpoa_df%>%
    group_by(DatabaseID,HPO_ID_TERM)%>%
    slice_max(Frequency_score,n=1,with_ties = F)%>%
    ungroup()
  
  # A table containing all the ancestors for each phenotype
  all_ancestors<-hpoa_df%>%select(HPO_ID)%>%distinct()%>%rowwise()%>%
    mutate(ancestors=paste0(get_hpo_ancestors_by_id(HPO_ID,with_term = F),collapse=','),
           num_ancestors=get_num_of_hpo_ancestors_by_id(HPO_ID))
  
  # A table containing all the phenotypes ids and terms
  all_phenos<-hpoa_df%>%select(HPO_ID,HPO_TERM,HPO_ID_TERM)%>%distinct()
  # join the hpo table with the descendants and then separate the descendants to different rows 
  #hpo_specificity_df<-readr::read_delim('./data/hpo_specificity.csv')
  #hpoa_df<-hpoa_df%>%left_join(hpo_specificity_df) # join with specificity before adding descendants
  hpoa_df<-hpoa_df%>%
    left_join(all_ancestors)%>%
    rowwise()%>%
    separate_rows(ancestors,sep=',')%>%
    mutate(is_ancestor=ifelse(HPO_ID==ancestors,F,T), # set it so that if the ancestor id is the same as the original id - it is not an ancestor
           ancestor_of=HPO_ID_TERM,
           HPO_ID=ancestors)%>%
    select(-c(HPO_ID_TERM,HPO_TERM))%>% # remove the term and id-term and repopulate them according to the new row HPO ID
    left_join(all_hpos_terms)
  
  hpo_specificity_df<-generate_hpo_specificity_table(hpoa_df)
  
  hpoa_df<-hpoa_df%>%
    left_join(hpo_specificity_df)%>%
    mutate(Frequency_score_with_specificity=Frequency_score/specificity)
  return(hpoa_df)
}

get_disease_list<-function(hpoa_df){
  disease_list<-hpoa_df%>%group_by(DatabaseID,DiseaseName)%>%
    summarize(hpos=list(HPO_ID_TERM))%>%ungroup()
  return(disease_list)
}

get_hpo_children_by_id<-function(hpo_id){
  ontologyIndex::get_descendants(hpo,hpo_id)
}

get_hpo_ancestors_by_id<-function(hpo_id,with_term=T){
  #print(hpo_id)
  ancestors<-ontologyIndex::get_ancestors(hpo,hpo_id)
  # if the hpo_id is not found - return NA
  if (length(ancestors)==0){
    message(glue('HPO ID {hpo_id} is not found..skipping..'))
    return(NA)
  }
  ancestors_terms<-ontologyIndex::get_term_property(ontology=hpo, property="ancestors", term=hpo_id, as_names=TRUE)
  ancestors_terms<-glue('{names(ancestors_terms)}:{ancestors_terms}')
  if(with_term){
    return(ancestors_terms)
  }else{
    return(ancestors)
  }
}

get_num_of_hpo_ancestors_by_id<-function(hpo_id){
  hpo_desendants<-ontologyIndex::get_ancestors(hpo,hpo_id)
  return(length(hpo_desendants)-1)
}

get_num_of_hpo_children_by_id<-function(hpo_id){
  hpo_desendants<-ontologyIndex::get_descendants(hpo,hpo_id)
  return(length(hpo_desendants)-1)
}

# !!!! BEFORE you run this script, make sure the updated version of the hpoa_df is loaded, and only then create the specificity
# table, otherwise it will cause discrepancy
generate_hpo_specificity_table<-function(hpoa_df){
  hpo_specificity_df<-hpoa_df%>%
    #filter(!is_descendant)%>%
    group_by(HPO_ID_TERM,HPO_ID,HPO_TERM)%>%
    filter(!Frequency_cat=='excluded')%>%
    summarize(specificity=sum(Frequency_score))
    #summarize(specificity=max(1,sum(Frequency_score)))# the minimal specificity should be 1 (otherwise if the frequency is rare and it only occurs once it is considered more specific than an obligate)
  write.table(hpo_specificity_df,file='./data/hpo_specificity.csv',sep = '\t',row.names = F,quote = F)
  return(hpo_specificity_df)
}

