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

version = 0.2
library(ProjectTemplate)
setwd('..')
load.project()
disease_table<-NULL
hpoa_df<-parse_hpo_hpoa_db()
disease_list<-get_disease_list(hpoa_df)
hpo_specificity_df<-readr::read_delim('./data/hpo_specificity.csv')
hpoa_df<-hpoa_df%>%left_join(hpo_specificity_df)
hpoa_df<-hpoa_df%>%mutate(Frequency_score_with_specificity=Frequency_score/specificity)

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
                                    div(selectizeInput('selected_phenotypes', label='selected phenotypes', multiple=TRUE,choices=unique(hpoa_df$HPO_ID_TERM),selected = NULL,width = 800), 
                                        style='font-size:200%;'),
                                    div(selectizeInput('ruled_out_phenotypes', label='ruled-out phenotypes', multiple=TRUE,choices=unique(hpoa_df$HPO_ID_TERM),selected = NULL,width = 800), 
                                        style='font-size:200%;'),
                                    #shinyWidgets::switchInput(inputId = "with_descendants", value = FALSE,onLabel = 'Include descendants',offLabel = 'Without descendants'),
                                    checkboxInput(inputId = "with_descendants", label = 'Add descendants',value = F),
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
        original_selected_phenotypes<- input$selected_phenotypes
        selected_phenotypes<-original_selected_phenotypes
        ruled_out_phenotypes<-input$ruled_out_phenotypes
        print(input$with_descendants)
        if (input$with_descendants){
            phenotype_ids<-stringr::str_extract(original_selected_phenotypes,'HP:\\d+')
            phenotype_descendants<-get_descendants(hpo,phenotype_ids)
            phenotype_terms<-hpoa_df%>%filter(HPO_ID%in%phenotype_descendants)%>%pull(HPO_ID_TERM)%>%unique()
            selected_phenotypes<-phenotype_terms
            message(glue('Selected phenotypes: {paste0(original_selected_phenotypes,collapse=", ")}, with descendants: {paste0(selected_phenotypes,collapse=", ")}'))
        }
        
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
                    filter(DatabaseID %in% relevant_diseases,HPO_ID_TERM%in%selected_phenotypes)%>%
                    group_by(DatabaseID,DiseaseName)%>%
                    summarize(selected_phenotype_freq=paste0(HPO_ID_TERM,"(",Frequency_cat,")",collapse=', ')))%>%
            left_join(
                hpoa_df%>%
                    filter(DatabaseID %in% relevant_diseases,!(HPO_ID_TERM%in%selected_phenotypes),Frequency_cat%in%c('obligate','very_frequent','frequent'))%>%
                    group_by(DatabaseID,DiseaseName)%>%
                    summarize(additional_frequent_phenos=paste0(HPO_ID_TERM,"(",Frequency_cat,")",collapse=', ')))%>%
            left_join(
                hpoa_df%>%
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
