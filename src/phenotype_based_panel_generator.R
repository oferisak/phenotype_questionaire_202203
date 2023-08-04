library(shiny)
library(shinydashboard)
library(DT)
library(tidyr)
library(reactable)
library(dplyr)
library(dashboardthemes)
library(xlsx)
library(stringr)


### add genes to each disorder
### consider using phevaluator to evaluate the method https://github.com/OHDSI/PheValuator

version = 0.1
library(ProjectTemplate)
library(ontologyIndex)
setwd('..')
data(hpo)
load.project()
disease_table<-NULL
# full_hpoa_df<-parse_hpo_hpoa_db()
# disease_list<-get_disease_list(full_hpoa_df)

load(file='./data/preprocessed_data.RData')
#test_table<-full_hpoa_df[1:100000,]
#write.table(test_table,'/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/test_table.csv',row.names = F,sep='\t')

# the main phenotypes defined by the user
main_phenos<-c('Hypertrophic cardiomyopathy')
# from the hpo table filter rows where the HPO_TERM column is found in the main_phenos
original_main_phenos<-full_hpoa_df%>%filter(HPO_TERM %in% main_phenos)
# make a copy of the original table of main phenos
selected_disorders<-original_main_phenos
# from the hpo table, filter only disorders that are found in the selected_disorders (by matching the DatabaseID column) and that the HPO_TERM is not in the main_phenos
selected_disorders_hpos_not_in_main<-full_hpoa_df%>%filter(DatabaseID %in% (selected_disorders%>%pull(DatabaseID)))%>%
  filter(!(HPO_TERM %in% main_phenos))
# a variable to keep all rejected phenotypes
rejected_phenos<-c()
# while there are HPO_TERMs in the selected_disorders_hpos_not_in_main that have a Frequency_cat column equal to "obligate" 
while (sum(selected_disorders_hpos_not_in_main$Frequency_cat=='obligate')>0){
  # for each "obligate" HPO_TERM count the number of times it appears in the table
  hpos_by_count<-selected_disorders_hpos_not_in_main%>%group_by(HPO_TERM,Frequency_cat)%>%
    filter(Frequency_cat=='obligate')%>%
    summarize(n=n())
  message(glue('There are {nrow(hpos_by_count)} phenotypes that i can ask about'))
  # then take the HPO_TERM with the highest number of appearances
  top_obligate_pheno<-hpos_by_count%>%ungroup()%>%filter(!(HPO_TERM%in%main_phenos))%>%slice_max(order_by =n,n=1,with_ties = F)
  # and ask the user whether his/her patient has the phenotype
  is_pheno_present<- readline(prompt = glue('does your patient have {top_obligate_pheno%>%pull(HPO_TERM)}? (y/n) '))
  # if the user responds with yes add that phenotype to the main_phenos
  if (is_pheno_present=='y'){
    main_phenos<-c(main_phenos,top_obligate_pheno$HPO_TERM)
    # selected_disorders_hpos_not_in_main<-full_hpoa_df%>%filter(DatabaseID %in% (selected_disorders%>%pull(DatabaseID)))%>%
    #   filter(!(HPO_TERM %in% main_phenos))
  }
  # if the user responds with no 
  if (is_pheno_present=='n'){
    # add the phenotype to the rejected_phenos variable
    rejected_phenos<-c(rejected_phenos,top_obligate_pheno$HPO_TERM)
    # collect all the disorders that have the HPO as obligate
    rejected_disorders<-selected_disorders_hpos_not_in_main%>%filter(Frequency_cat=='obligate')%>%filter(HPO_TERM %in% rejected_phenos)%>%pull(DatabaseID)
    # remove these disorders from the selected_disorders table
    selected_disorders<-selected_disorders_hpos_not_in_main%>%filter(!(DatabaseID %in% rejected_disorders))
    # regenerate the selected_disorders_hpos_not_in_main by collecting all phenotypes found in disorders that are in the selected_disorders table and remove phenotypes that are already selected as main_phenos
    selected_disorders_hpos_not_in_main<-full_hpoa_df%>%filter(DatabaseID %in% (selected_disorders%>%pull(DatabaseID)))%>%
      filter(!(HPO_TERM %in% main_phenos))
  }
}
# filter all rejected disorders from the original table of main_phenos
final_selected_disorders_ids<-original_main_phenos%>%filter(!(DatabaseID %in% rejected_disorders))%>%pull(DatabaseID)
# now collect all the phenotypes that match the kept disorders
final_selected_disorders<-full_hpoa_df%>%filter(DatabaseID%in%final_selected_disorders)
# now remove disorders that just dont match by collecting their frequet
disorders_to_examine<-final_selected_disorders%>%group_by(DatabaseID,DiseaseName)%>%
  summarize('selected_phenos'=paste0(unique(HPO_TERM[HPO_TERM%in%main_phenos]),collapse = ', '),
            'rejected_phenos'=paste0(unique(HPO_TERM[HPO_TERM%in%rejected_phenos]),collapse = ', '),
            'all_frequent_phenos'=paste0(unique(HPO_TERM[!is_ancestor & (Frequency_cat %in% c('unknown','obligate','very_frequent','frequent'))]),collapse = '\n'))%>%ungroup()

for (i in 1:nrow(disorders_to_examine)){
  disorder_name<-disorders_to_examine%>%slice(i)%>%pull(DiseaseName)
  disorder_id<-disorders_to_examine%>%slice(i)%>%pull(DatabaseID)
  disorder_frequent_phenos<-disorders_to_examine%>%slice(i)%>%pull(all_frequent_phenos)
  is_disorder_relevant<- readline(prompt = glue('Could your patient match {disorder_name} frequent_phenos include:\n {disorder_frequent_phenos}? (y/n)\n '))
  if (is_disorder_relevant=='n'){rejected_disorders<-c(rejected_disorders,disorder_id)}
}

# hpo_specificity_df<-readr::read_delim('./data/hpo_specificity.csv')
# full_hpoa_df<-full_hpoa_df%>%left_join(hpo_specificity_df)
# full_hpoa_df<-full_hpoa_df%>%mutate(Frequency_score_with_specificity=Frequency_score/specificity)

moi_hpos<-c('HP:0000006','HP:0000007','HP:0001425','HP:0001426','HP:0001427','HP:0001428','HP:0001466','HP:0001472','HP:0003743','HP:0003745','HP:0010985')

# Define UI ####
ui <- dashboardPage(dashboardHeader(title = sprintf('Phenotype questionaire'),
                                    titleWidth = 290), 
                    dashboardSidebar(width = 290,
                                     sidebarMenu(menuItem("Phenotype browser",
                                                          tabName = "phenotype_browser", 
                                                          icon = icon('search'))
                                                 
                                     )
                    ),
                    
                    dashboardBody(
                      # stlye definitions from https://stackoverflow.com/questions/52198452/how-to-change-the-background-color-of-the-shiny-dashboard-body
                      
                      
                      tags$head(tags$style(HTML('
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #61c0bd;
                                font-weight: 800;
                                }
                                /* navbar (rest of the header) */
                                .skin-blue .main-header .navbar {
                                background-color: #61c0bd;
                                }
                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #f5fbf8;
                                }
                                /* main sidebar */
                                .main-sidebar {
                                font-size: 20px;
                                }
                                /* sidebar open */
                                .treeview-menu>li>a {
                                font-size: 18px!important;
                                }
                                .dataTables_scrollBody {
                                    transform:rotateX(180deg);
                                }
                                .dataTables_scrollBody table {
                                    transform:rotateX(180deg);
                                }
                                '))),
                      
                      #shinyDashboardThemes(theme = "poor_mans_flatly"),
                      tabItems(
                        # Panel Browser
                        tabItem('phenotype_browser',
                                div(selectizeInput('selected_phenotypes', label='selected phenotypes', multiple=TRUE,choices=unique(full_hpoa_df$HPO_ID_TERM),selected = NULL,width = 800), 
                                    style='font-size:200%;'),
                                div(selectizeInput('ruled_out_phenotypes', label='ruled-out phenotypes', multiple=TRUE,choices=unique(full_hpoa_df$HPO_ID_TERM),selected = NULL,width = 800), 
                                    style='font-size:200%;'),
                                #shinyWidgets::switchInput(inputId = "with_descendants", value = FALSE,onLabel = 'Include descendants',offLabel = 'Without descendants'),
                                checkboxInput(inputId = "with_ancestors", label = 'Add Ancestors',value = F),
                                DT::dataTableOutput(outputId='specific_phenos'),
                                DT::dataTableOutput(outputId='diseases_table')
                        )
                      ))
                    #tags$head(tags$style(HTML('* {font-family: "Kirnberg"};')))
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  
  observeEvent(ignoreInit = TRUE, c(input$selected_phenotypes,input$ruled_out_phenotypes) ,{
    # select only diseases with these phenotypes
    #disease_list%>%filter(mapply(function(x,y) length(setdiff(selected_phenotypes,y))==0, selected_phenotypes, disease_list$hpos))%>%slice(1)%>%pull(hpos)
    original_selected_phenotypes<-input$selected_phenotypes
    if (input$with_ancestors){
      hpoa_df<-full_hpoa_df
      selected_phenos_hpos<-stringr::str_extract(original_selected_phenotypes,'HP:\\d+')
      selected_phenos_ancestors<-unique(unlist(purrr::map(selected_phenos_hpos,function(x) get_hpo_ancestors_by_id(x,with_term = T))))
      selected_phenotypes<- selected_phenos_ancestors
      message(glue('Selected phenotypes: {paste0(original_selected_phenotypes,collapse=", ")}, with ancestors: {paste0(selected_phenotypes,collapse=", ")}'))
    }else{
      hpoa_df<-full_hpoa_df%>%filter(!is_ancestor)
      selected_phenotypes<-original_selected_phenotypes
    }
    ruled_out_phenotypes<-input$ruled_out_phenotypes
    #print(input$with_descendants)
    # if (input$with_descendants){
    #     phenotype_ids<-stringr::str_extract(original_selected_phenotypes,'HP:\\d+')
    #     phenotype_descendants<-get_descendants(hpo,phenotype_ids)
    #     phenotype_terms<-hpoa_df%>%filter(HPO_ID%in%phenotype_descendants)%>%pull(HPO_ID_TERM)%>%unique()
    #     selected_phenotypes<-
    #     message(glue('Selected phenotypes: {paste0(original_selected_phenotypes,collapse=", ")}, with descendants: {paste0(selected_phenotypes,collapse=", ")}'))
    # }
    
    # relevant_diseases<-disease_list%>%
    #     filter(mapply(function(x,y) length(setdiff(selected_phenotypes,y))==0, 
    #                   selected_phenotypes, disease_list$hpos))%>%
    #     pull(DatabaseID)
    
    relevant_diseases<-disease_list%>%
      filter(mapply(function(x,y) length(intersect(selected_phenotypes,y))>0, 
                    selected_phenotypes, disease_list$hpos))%>%
      pull(DatabaseID)
    print(length(relevant_diseases))
    
    #print(relevant_diseases)
    
    # most specific phenotypes
    specific_phenos<-hpoa_df%>%filter(DatabaseID %in% relevant_diseases)%>%
      group_by(HPO_ID_TERM)%>%
      summarize(n=n())
    
    updateSelectInput('selected_phenotypes',
                      session=getDefaultReactiveDomain(),
                      choices = specific_phenos%>%pull(HPO_ID_TERM),
                      selected = original_selected_phenotypes)
    
    output$specific_phenos<-DT::renderDataTable(
      DT::datatable(specific_phenos%>%arrange(desc(n)),
                    # datatable definitions
                    options=list(scrollX=F,pageLength=10),
                    selection='single',
                    rownames= FALSE,
                    filter = list(position = 'top', clear = FALSE))
    )
    
    # disease_table_with_selected_phenotypes<-hpoa_df%>%
    #     filter(DatabaseID %in% relevant_diseases)%>%
    #     group_by(DatabaseID,DiseaseName)%>%
    #     summarize(phenotype_score=sum(ifelse(HPO_ID_TERM%in%input$selected_phenotypes,Frequency_score,0)),
    #               phenotype_freq=paste0(HPO_ID_TERM,"(",Frequency_cat,")",collapse=', '))%>%
    #     arrange(desc(phenotype_score))
    # 
    disease_table_with_selected_phenotypes<-hpoa_df%>%
      filter(DatabaseID %in% relevant_diseases)%>%
      group_by(DatabaseID,DiseaseName)%>%
      # summarize(phenotype_score=sum(ifelse(HPO_ID_TERM%in%selected_phenotypes,Frequency_score,0)),
      #           phenotype_score_with_spec=sum(ifelse(HPO_ID_TERM%in%selected_phenotypes,Frequency_score_with_specificity,0)))%>%
      summarize(phenotype_score=sum(ifelse(HPO_ID_TERM%in%selected_phenotypes,Frequency_score,0))-sum(ifelse(HPO_ID_TERM%in%ruled_out_phenotypes & Frequency_cat!='excluded',Frequency_score,0)),
                phenotype_score_with_spec=sum(ifelse(HPO_ID_TERM%in%selected_phenotypes,Frequency_score_with_specificity,0))-sum(ifelse(HPO_ID_TERM%in%ruled_out_phenotypes & Frequency_cat!='excluded',Frequency_score,0)))%>%
      left_join(
        hpoa_df%>%
          filter(DatabaseID %in% relevant_diseases,HPO_ID%in%moi_hpos)%>%
          group_by(DatabaseID,DiseaseName)%>%
          summarize(moi=paste0(HPO_TERM,collapse=', ')))%>%
      left_join(
        hpoa_df%>%
          filter(!is_ancestor)%>%
          filter(DatabaseID %in% relevant_diseases,HPO_ID_TERM%in%selected_phenotypes)%>%
          group_by(DatabaseID,DiseaseName)%>%
          summarize(selected_phenotype_freq=paste0(HPO_ID_TERM,"(",Frequency_cat,")",collapse=', ')))%>%
      left_join(
        hpoa_df%>%
          filter(!is_ancestor)%>%
          filter(DatabaseID %in% relevant_diseases,!(HPO_ID_TERM%in%selected_phenotypes),Frequency_cat%in%c('obligate','very_frequent','frequent'))%>%
          group_by(DatabaseID,DiseaseName)%>%
          summarize(additional_frequent_phenos=paste0(HPO_ID_TERM,"(",Frequency_cat,")",collapse=', ')))%>%
      left_join(
        hpoa_df%>%
          filter(!is_ancestor)%>%
          filter(DatabaseID %in% relevant_diseases,!(HPO_ID_TERM%in%selected_phenotypes),!(Frequency_cat%in%c('obligate','very_frequent','frequent')))%>%
          group_by(DatabaseID,DiseaseName)%>%
          summarize(additional_rare_phenos=paste0(HPO_ID_TERM,"(",Frequency_cat,")",collapse=', ')))%>%
      arrange(desc(phenotype_score_with_spec))
    
    
    output$diseases_table<-DT::renderDataTable(
      DT::datatable(disease_table_with_selected_phenotypes,
                    # datatable definitions
                    options=list(scrollX=T,pageLength=10),
                    selection='single',
                    rownames= FALSE,
                    filter = list(position = 'top', clear = FALSE))
    )
  })
  
}
#  disable = list(columns = 1:(ncol(final_panel_output)-1)
# Run the application 
shinyApp(ui = ui, server = server)
