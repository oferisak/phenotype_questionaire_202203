hpoa_path<-'./data/phenotype.hpoa'
gene_to_phenotype_path<-'./data/genes_to_phenotype.txt'

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
  hpoa_df<-readr::read_delim(hpoa_path,delim='\t',col_types =paste0(rep('c',12),collapse=''))
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
  hpoa_df$Frequency_numeric<- as.numeric(
    unlist(
    purrr::map(
      hpoa_df$Frequency,
      ~ ifelse(
        .x %in% names(hpo_frequencies),
        hpo_frequencies[[.x]],
        ifelse(
          grepl('%', .x),
          as.numeric(stringr::str_replace(.x, '%', '')) /
            100,
          ifelse(
            grepl('/', .x),ifelse(as.numeric(stringr::str_split(.x,'/')[[1]][1])<=as.numeric(stringr::str_split(.x, '/')[[1]][2]),
                                  as.numeric(stringr::str_split(.x,'/')[[1]][1]) /
                                    as.numeric(stringr::str_split(.x, '/')[[1]][2]),
                                  as.numeric(stringr::str_split(.x,'/')[[1]][2]) /
                                    as.numeric(stringr::str_split(.x, '/')[[1]][1])),
            .x
          )
        )
      ))))
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
    
  return(hpoa_df)
  
  
}

get_disease_list<-function(hpoa_df){
  disease_list<-hpoa_df%>%group_by(DatabaseID,DiseaseName)%>%
    summarize(hpos=list(HPO_ID_TERM))%>%ungroup()
  return(disease_list)
}

library(ontologyIndex)
data(hpo)
get_hpo_children_by_id<-function(hpo_id){
  ontologyIndex::get_descendants(hpo,hpo_id)
}

# !!!! BEFORE you run this script, make sure the updated version of the hpoa_df is loaded, and only then create the specificity
# table, otherwise it will cause discrepancy
generate_hpo_specificity_table<-function(){
  hpo_specificity_df<-hpoa_df%>%group_by(HPO_ID_TERM,HPO_ID,HPO_TERM)%>%
    filter(!Frequency_cat=='excluded')%>%
    summarize(specificity=sum(Frequency_score))
  write.table(hpo_specificity_df,file='./data/hpo_specificity.csv',sep = '\t',row.names = F,quote = F)
}
