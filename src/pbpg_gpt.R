
version = 0.23
library(shiny)
library(shinyWidgets)
library(ProjectTemplate)
library(ontologyIndex)
library(bslib)
library(data.table)
setwd('..')
data(hpo)
load.project()

# prepare the data and save it
# full_hpoa_df<-parse_hpo_hpoa_db()
# save(full_hpoa_df,file='/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/preprocessed_data.RData')

# Tried: 
## modify the questions order to be according to the total frequency in the different disorders (more frequent will be asked before) - 
##  - did not result in a major improvement in the number of questions required

# Read the test table
load('/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/preprocessed_data.RData')
disorder_to_gene<-create_disorder_to_gene_table()
# Read panelapp panels
message('Parsing panelapp panels..')
panelapp<-readr::read_delim('./data/panelapp_db_2023-08-17.csv.gz')

ui <- fluidPage(
  theme = bslib::bs_theme(
    base_font = "Arial",
    heading_font = "Lato",
    bg = "#f7f7f7",
    fg = "#333",
    primary = "#4b3359"
  ),
  
  # ui: tags - head ####
  tags$head(
    tags$style(HTML("
      .main-header {text-align: center;} 
      .subtitle {text-align: center;} 
      .bootstrap-select.btn-group .dropdown-toggle { width: 500px !important; }
      .shiny-input-container { width: 500px !important; }
       h1 { font-weight: bold !important; } /* Bold title header */
      .nav.nav-tabs > li > a { font-weight: bold !important; } /* Bold tab headers */

    "))
  ),
  br(),
  titlePanel(
    tags$h1(tags$i(class = "bi bi-file-earmark-text"), " Phenotype-Questionnaire Gene Panel Generator")
  ),
  br(),
  tabsetPanel(
    # ui: tab - phenotype-questionnaire ####
    tabPanel("Phenotype-questionnaire",
             mainPanel(
               br(),
               fluidRow(
                 column(10,
                        div(style = "max-width: 500px;",  # You can adjust the max-width value as needed
                               sliderInput(inputId = "analysis_depth", 
                                           label = strong("Analysis depth"),
                                           value=-1, max = -0.1, min=-2, step = 0.1),
                               div(style = "display: flex; justify-content: space-between; font-size:10px;",
                                   tags$span("(More questions | More specific)"),
                                   tags$span("(Less questions | More sensitive)")
                               )
                 ))
               ),
               br(),
               fluidRow(
                 column(10, 
                        selectizeInput("main_pheno" ,strong("Select the Main Phenotypes:"), choices = NULL, multiple = T)
                 ),
                 column(10, 
                        selectizeInput("panelapp_panel" ,strong("OR: Select PanelApp Panels:"), choices = unique(panelapp$panel_name), multiple = T)
                 ),
                 column(10, 
                        selectizeInput("additional_pheno" ,"Additional Phenotypes Present:", choices = NULL, multiple = T)
                 ),
                 column(10, 
                        selectizeInput("possible_pheno" ,"Possible Phenotypes:", choices = NULL, multiple = T)
                 ),
                 column(10, 
                        selectizeInput("rejected_pheno" ,"Rejected Phenotypes:", choices = NULL, multiple = T)
                 )
               ),
               br(),
               fluidRow(
                 column(4,
                        actionButton(inputId = "start_analysis", label = tags$span(" Start Analysis"), class = "btn-primary")
                 ),
                 column(6,
                        textOutput("questionsAnsweredText") )
               ),
               br(),
               fluidRow(
                 column(8, dataTableOutput("hpos_by_count_table"))
               ),
               br()
             )
    ),
    # ui: tab - possible disorders ####
    tabPanel("Possible disorders",
             mainPanel(
               br(),
               fluidRow(
                 column(10,downloadButton("download_disorders", "Download Disorders Table"))
               ),
               br(),
               fluidRow(
                 column(10,dataTableOutput("disorders_table"))
               )
             )
    ),
    # ui: tab - panel genes ####
    tabPanel("Panel genes",
             mainPanel(
               br(),
               fluidRow(
                 column(10,sliderInput("likelihood_threshold", "Likelihood Threshold:", -1, min = -5, max = 0, step = 0.1))
               ),
               br(),
               fluidRow(
                 column(10,downloadButton("download_panel_genes", "Download Panel Genes Table"))
               ),
               br(),
               fluidRow(
                 column(10,dataTableOutput('panel_genes_table'))
               )
             )
    )
  )
)


server <- function(input, output, session) {
  # server: define reactive variables ####
  main_phenos_reactive <- reactiveVal()
  additional_phenos_present_reactive <- reactiveVal(c())
  additional_phenos_possible_reactive <- reactiveVal(c())
  panelapp_panels_reactive<-reactiveVal()
  maybe_phenos_reactive<-reactiveVal(c())
  rejected_phenos <- reactiveVal(c())
  rejected_disorders <- reactiveVal(c())
  top_obligate_pheno <- reactiveVal()
  modal_shown <- reactiveVal(FALSE)
  questions_answered <- reactiveVal(0)
  
  # server: update main phenotypes input ####
  observe({
    updateSelectizeInput(session, "main_pheno", choices = unique(full_hpoa_df$hpo_id_name),server=TRUE)
  })
  
  observe({
    updateSelectizeInput(session, "additional_pheno", choices = unique(full_hpoa_df$hpo_id_name),server=TRUE)
  })
  
  observe({
    if (!is.null(input$main_pheno) && length(input$main_pheno) > 0) {
      updateSelectizeInput(session, "panelapp_panel", selected = character(0))
    }
  })
  
  observe({
    if (!is.null(input$panelapp_panel) && length(input$panelapp_panel) > 0) {
      updateSelectizeInput(session, "main_pheno", selected = character(0))
    }
  })
  
  # server: observer: clicking start analysis button ####
  observeEvent(input$start_analysis, {
    updateActionButton(session, "start_analysis", label = "Resume Analysis")
    main_phenos <- input$main_pheno
    main_phenos_reactive(main_phenos)
    additional_phenos_present<-input$additional_pheno
    additional_phenos_present_reactive(additional_phenos_present)
    
    selected_panelapp_panels<-input$panelapp_panel
    panelapp_panels_reactive(selected_panelapp_panels)
    modal_shown(FALSE)
  })
  # server - output - number of questions answered counter ####
  output$questionsAnsweredText <- renderText({
    paste("You have answered", questions_answered(), "questions so far.")
  })
  
  # server: reactive - all_phenos_data ####
  all_phenos_data <- reactive({
    message('Running reactive: all_phenos_data')
    main_phenos <- main_phenos_reactive()
    selected_panelapp_panels<-panelapp_panels_reactive()

    if (is.null(main_phenos) & is.null(selected_panelapp_panels)) return(NULL)
    
    # If generating panel using panelapp panels
    if (!is.null(selected_panelapp_panels)){
      
      panelapp_genes<-panelapp%>%filter(panel_name%in%selected_panelapp_panels)%>%pull(gene_symbol)
      panelapp_disorders<-disorder_to_gene%>%filter(gene_symbol%in%panelapp_genes)
      all_disorders_with_main_phenos <- full_hpoa_df %>%
        filter(disease_id %in% panelapp_disorders$disease_id) %>%
        pull(disease_id) %>%
        unique()
    }
    
    # If generating panel using main HPOs
    if (!is.null(main_phenos)){
      all_disorders_with_main_phenos <- full_hpoa_df %>%
        filter(hpo_id_name %in% main_phenos) %>%
        pull(disease_id) %>%
        unique()
    }
    
    all_phenos_from_disorders_with_main_phenos <- full_hpoa_df %>%
      filter(disease_id %in% all_disorders_with_main_phenos)
    
    list(
      main_phenos = main_phenos,
      all_phenos_from_disorders_with_main_phenos = all_phenos_from_disorders_with_main_phenos#%>%filter(!is_ancestor)
    )
  })
  
  # server: reactive - hpos_by_count_reactive ####
  hpos_by_count_reactive <- reactive({
    message('Running reactive: hpos_by_count_reactive')
    phenos_data <- all_phenos_data()
    if (is.null(phenos_data)) return(NULL)

    # Exclude all ancestors of phenos that the user said are present or maybe present in the patient
    excluded_due_to_addition <- c(phenos_data$main_phenos, 
                                  additional_phenos_present_reactive(),
                                  additional_phenos_possible_reactive(),
                                  maybe_phenos_reactive())
    excluded_due_to_addition_ancestors<-NULL
    for (excluded_pheno_id_term in excluded_due_to_addition){
      excluded_pheno_id<-stringr::str_extract(excluded_pheno_id_term,'HP:\\d+')
      excluded_pheno_ancestors<-get_hpo_ancestors_by_id(excluded_pheno_id,with_term = T)
      excluded_due_to_addition_ancestors<-unique(c(excluded_due_to_addition_ancestors,excluded_pheno_ancestors))
    }
    
    # Exclude all descendants of phenos that the user said are present in the patient
    excluded_due_to_addition <- c(phenos_data$main_phenos, 
                                  additional_phenos_present_reactive(),
                                  additional_phenos_possible_reactive())
    excluded_due_to_addition_descendants<-NULL
    for (excluded_pheno_id_term in excluded_due_to_addition){
      excluded_pheno_id<-stringr::str_extract(excluded_pheno_id_term,'HP:\\d+')
      excluded_pheno_descendants_ids<-get_hpo_descendants_by_id(excluded_pheno_id,with_term = F)
      excluded_pheno_descendants<-phenos_data$all_phenos_from_disorders_with_main_phenos%>%
        filter(hpo_id %in% excluded_pheno_descendants_ids)%>%
        pull(hpo_id_name)
      excluded_due_to_addition_descendants<-unique(c(excluded_due_to_addition_descendants,excluded_pheno_descendants))
    }
    
    # Exclude all descendants of phenos that the user said are not in the patient
    excluded_due_to_rejection <- c(rejected_phenos())
    excluded_due_to_rejection_descendants<-NULL
    for (excluded_pheno_id_term in excluded_due_to_rejection){
      excluded_pheno_id<-stringr::str_extract(excluded_pheno_id_term,'HP:\\d+')
      excluded_pheno_descendants_ids<-get_hpo_descendants_by_id(excluded_pheno_id,with_term = F)
      excluded_pheno_descendants<-phenos_data$all_phenos_from_disorders_with_main_phenos%>%
        filter(hpo_id %in% excluded_pheno_descendants_ids)%>%
        pull(hpo_id_name)
      excluded_due_to_rejection_descendants<-unique(c(excluded_due_to_rejection_descendants,excluded_pheno_descendants))
    }
    
    disorders_with_likelihood<-update_disorders_likelihood()
    disorders_remaining<-disorders_with_likelihood%>%filter(likelihood>input$analysis_depth)%>%pull(disease_id)
    unlikely_disorders<-disorders_with_likelihood%>%filter(likelihood<=input$analysis_depth)%>%pull(disease_id)
    message(glue('There are currently {length(unlikely_disorders)}/{nrow(disorders_with_likelihood)} unlikely disorders..'))
    # Modify to exclude phenos in main and additional_pheno inputs
    #excluded_phenos <- c(main_phenos_ancestors, additional_phenos_present_reactive(),maybe_phenos_reactive(),rejected_phenos())
    hpos_from_disorders_that_are_not_in_main <- phenos_data$all_phenos_from_disorders_with_main_phenos %>%
      #filter(!(disease_id %in% rejected_disorders())) %>%
      filter(!(disease_id %in% unlikely_disorders)) %>%
      filter(!(hpo_id_name %in% excluded_due_to_addition_ancestors))%>%
      filter(!(hpo_id_name %in% excluded_due_to_addition_descendants))%>%
      filter(!(hpo_id_name %in% excluded_due_to_rejection_descendants))
    
    
    final_descendants <- as.data.table(hpos_from_disorders_that_are_not_in_main %>%
                                         filter(is_ancestor == TRUE) %>%
                                         group_by(hpo_id, hpo_id_name) %>%
                                         summarize(final_desc = paste0(unique(ancestor_of)[order(unique(ancestor_of))], collapse='<br>'),
                                                   num_final_desc = length(unique(ancestor_of))))
    
    # final_descendants<-NULL # 
    hpos_by_count <- hpos_from_disorders_that_are_not_in_main %>%
      #filter(!is_ancestor)%>%
      filter(frequency_cat %in% c('obligate','very_frequent','frequent','unknown')) %>%
      #filter(frequency_cat %in% c('obligate','very_frequent')) %>%
      group_by(hpo_id_name) %>%
      summarize(n = length(unique(disease_id)),freq_sum=sum(frequency_numeric,na.rm = T))%>%
      arrange(desc(n))
    
    list(hpos_by_count=hpos_by_count,
         final_descendants=final_descendants,
         disorders_remaining=disorders_remaining)
  })

  # Observer to update the choices of rejected phenotypes selectizeInput
  observe({
    rejected <- rejected_phenos()
    all_phenos <- all_phenos_data()$all_phenos_from_disorders_with_main_phenos$hpo_id_name
    updateSelectizeInput(session, inputId = 'rejected_pheno', choices = unique(all_phenos), selected = rejected)
  })
  
  observe({
    phenos_data <- all_phenos_data()
    if (is.null(phenos_data)) return()
    hpos_by_count_reactive_output<-hpos_by_count_reactive()
    hpos_by_count<-hpos_by_count_reactive_output$hpos_by_count
    disorders_remaining<-hpos_by_count_reactive_output$disorders_remaining
    output$hpos_by_count_table <- renderDataTable({
      hpos_by_count_reactive()$hpos_by_count
    })
    num_of_phenos_that_are_obligated <- nrow(hpos_by_count)
    top_pheno <- hpos_by_count %>%
      ungroup() %>%
      filter(!(hpo_id_name %in% phenos_data$main_phenos)) %>%
      slice_max(order_by = n, n = 1, with_ties = F)
    
    top_obligate_pheno(top_pheno)
    #print(top_obligate_pheno())
    top_pheno_final_descendants<-hpos_by_count_reactive_output$final_descendants
    # if the hpo in question is not an ancestor of anything
    #print(top_pheno_final_descendants)
    if (!is.null(top_pheno_final_descendants)){
      is_final_descendant<-!(top_obligate_pheno()$hpo_id_name %in% (top_pheno_final_descendants%>%pull(hpo_id_name)))
      top_pheno_final_descendants_text<-ifelse(is_final_descendant,
                                               '',
                                               top_pheno_final_descendants%>%
                                                 filter(hpo_id_name==(top_obligate_pheno()$hpo_id_name))%>%
                                                 pull(final_desc))
        
    }
    if (num_of_phenos_that_are_obligated > 0 & !modal_shown()) {
      question_text<-ifelse(top_pheno_final_descendants_text=='' | nchar(top_pheno_final_descendants_text)>800,
                            glue('There are {nrow(hpos_by_count)} more phenotypes from {length(disorders_remaining)} different disorders to ask about.Does your patient have:<br><b>{top_pheno %>% pull(hpo_id_name)}</b><br>'),
                            glue("There are {nrow(hpos_by_count)} more phenotypes from {length(disorders_remaining)} different disorders to ask about.Does your patient have:<br><b>{top_pheno %>% pull(hpo_id_name)}</b><br>Specifically:<br>{top_pheno_final_descendants_text}<br>"))
      print(question_text)
      showModal(modalDialog(
        title = "Confirm Phenotype",
        #glue("There are {nrow(hpos_by_count)} more phenotypes to ask about.\nDoes your patient have {top_pheno %>% pull(hpo_id_name)}\nSpecifically:\n{top_pheno_final_descendants}?"),
        HTML(question_text),
        footer = tagList(
          actionButton("yes_button", "Yes"),
          actionButton("possible_button", "Possible"),
          actionButton("be_more_specific_button", "Be more specific"),
          actionButton("no_button", "No"),
          actionButton("stop_button", "Finish")
        ),
        div(HTML(glue('<br>You have answered {questions_answered()} questions so far.')),style = "font-size:15px;")
      ))
      modal_shown(TRUE)
    } else {
      #output$hpos_by_count_table <-NULL
      # output$hpos_by_count_table <- renderTable({
      #   hpo_combinations_reactive()
      #}, rownames = FALSE)
      
      return()
    }
    
  })
  # server: observer - stop_button ####
  observeEvent(input$stop_button, {
    modal_shown(TRUE)
    removeModal()  

  })
  # server: observer - yes_button ####
  observeEvent(input$yes_button, {
    modal_shown(FALSE)
    questions_answered(questions_answered() + 1)
    removeModal()
    additional_phenos <- c(additional_phenos_present_reactive(), top_obligate_pheno()$hpo_id_name)
    additional_phenos_present_reactive(additional_phenos)
    updateSelectizeInput(session, inputId = 'additional_pheno', choices = unique(full_hpoa_df$hpo_id_name),selected = additional_phenos,server=TRUE) # Update additional_pheno input
  })
  
  # server: observer - possible_button ####
  observeEvent(input$possible_button, {
    modal_shown(FALSE)
    questions_answered(questions_answered() + 1)
    removeModal()
    possible_phenos <- c(additional_phenos_possible_reactive(), top_obligate_pheno()$hpo_id_name)
    additional_phenos_possible_reactive(possible_phenos)
    updateSelectizeInput(session, inputId = 'possible_pheno', choices = unique(full_hpoa_df$hpo_id_name),selected = possible_phenos,server=TRUE) # Update additional_pheno input
  })
  # server: observer - no_button ####
  observeEvent(input$no_button, {
    modal_shown(FALSE)
    questions_answered(questions_answered() + 1)
    removeModal()
    phenos_data <- all_phenos_data()
    rejected_phenos_current <- c(rejected_phenos(), top_obligate_pheno()$hpo_id_name)
    rejected_phenos(rejected_phenos_current)
  })
  
  # server: observer - be_more_specific_button ####
  observeEvent(input$be_more_specific_button, {
    modal_shown(FALSE)
    questions_answered(questions_answered() + 1)
    removeModal()
    maybe_phenos <- c(maybe_phenos_reactive(), top_obligate_pheno()$hpo_id_name)
    maybe_phenos_reactive(maybe_phenos)
    # No modifications made to any lists, simply continue to the next question.
  })

  
  # Function to create the disorders table
  
  create_disorders_table <- function(all_phenos_data, rejected_disorders) {
    # print(all_phenos_data$all_phenos_from_disorders_with_main_phenos %>%
    #   filter(!(disease_id %in% rejected_disorders))%>%count(frequency_cat))
    
    disorders <- all_phenos_data$all_phenos_from_disorders_with_main_phenos %>%
      filter(!is_ancestor)%>%
      filter(frequency_cat %in% c('obligate','very_frequent','frequent','unknown'))%>%
      filter(!(disease_id %in% rejected_disorders)) %>%
      group_by(disease_id, disease_name, frequency_cat) %>%
      summarise(hpo_id_name = paste(unique(hpo_id_name), collapse = ", ")) %>%
      pivot_wider(names_from = frequency_cat, values_from = hpo_id_name, values_fill = "")%>%
      left_join(disorder_to_gene%>%group_by(disease_name)%>%dplyr::summarize(genes=paste0(unique(gene_symbol),collapse=', ')))
    
    return(disorders)
  }
  
  update_disorders_likelihood_v1<-function(){
    all_phenos_data <- all_phenos_data()
    rejected_phen <- rejected_phenos()
    disorders_with_likelihood<-all_phenos_data$all_phenos_from_disorders_with_main_phenos%>%
      mutate(is_rejected=ifelse(hpo_id_name %in% rejected_phen,TRUE,FALSE),
             hpo_likelihood_score=ifelse(is_rejected & frequency_cat!='excluded',-frequency_score,0))%>%
      # if you want only one occurance of each hpo term per disorder group by hpo_id_name use distinct, otherwise remove it from grouping
      select(disease_id,hpo_id_name,hpo_likelihood_score)%>%distinct()%>%
      group_by(disease_id)%>%
      summarize(likelihood=sum(hpo_likelihood_score))%>%
      # re-add the disease name
      left_join(disorder_to_gene%>%group_by(disease_id,disease_name)%>%summarize(genes=paste0(gene_symbol,collapse=','))%>%ungroup())
    return(disorders_with_likelihood)
  }
  
  update_disorders_likelihood<-function(){
    all_phenos_data <- all_phenos_data()
    main_phen<-main_phenos_reactive()
    additional_phen<-additional_phenos_present_reactive()
    rejected_phen <- rejected_phenos()
    disorders_with_likelihood<-all_phenos_data$all_phenos_from_disorders_with_main_phenos%>%
      mutate(is_rejected=ifelse(hpo_id_name %in% rejected_phen,TRUE,FALSE),
             is_main=ifelse(hpo_id_name %in% main_phen,TRUE,FALSE),
             is_additional=ifelse(hpo_id_name %in% additional_phen,TRUE,FALSE),
             hpo_likelihood_score=ifelse(is_rejected & frequency_cat!='excluded',(-frequency_score),0),
             hpo_confidence_score=case_when(is_rejected & frequency_cat!='excluded'~(-frequency_score),
                                            is_main & frequency_cat=='excluded'~-2,
                                            is_main & frequency_cat!='excluded'~frequency_score,
                                            is_additional & frequency_cat!='excluded'~frequency_score,
                                            TRUE~0))%>%
      # if you want only one occurance of each hpo term per disorder group by hpo_id_name use distinct, otherwise remove it from grouping
      select(disease_id,hpo_id_name,hpo_likelihood_score,hpo_confidence_score)%>%distinct()%>%
      group_by(disease_id)%>%
      summarize(likelihood=sum(hpo_likelihood_score),confidence=sum(hpo_confidence_score))%>%
      # re-add the disease name
      left_join(disorder_to_gene%>%group_by(disease_id,disease_name)%>%
                  summarize(genes=paste0(gene_symbol,collapse=','))%>%ungroup())
    return(disorders_with_likelihood)
  }
  
  # Reactive to prepare the data for the disorders table
  disorders_data_reactive <- reactive({
    message('Running reactive: disorders_data_reactive')
    all_phenos_data <- all_phenos_data()
    
    if (is.null(all_phenos_data)) return(NULL)
    disorders_with_likelihood<-update_disorders_likelihood()  
    disorders_with_likelihood<-disorders_with_likelihood%>%
      left_join(
        all_phenos_data$all_phenos_from_disorders_with_main_phenos%>%
          filter(!is_ancestor)%>%
          group_by(disease_id, disease_name, frequency_cat) %>%
          summarise(hpo_id_name = paste(unique(hpo_id_name), collapse = ", ")) %>%
          pivot_wider(names_from = frequency_cat, values_from = hpo_id_name, values_fill = "")%>%
          ungroup()
      )%>%
      arrange(desc(likelihood))
    disorders_with_likelihood

  })
  
  panel_genes_reactive<-reactive({
    message('Running reactive: panel_genes_reactive')
    if (is.null(all_phenos_data)) return(NULL)
    disorders_with_likelihood<-update_disorders_likelihood()
    selected_panelapp_panels<-panelapp_panels_reactive()
    disorder_to_gene_for_table<-disorder_to_gene
    panelapp_genes_not_in_disorder_to_gene_table<-NULL
    if (!is.null(selected_panelapp_panels)){
      # first add the panelapp confidence to all the genes
      disorder_to_gene_for_table<-disorder_to_gene_for_table%>%
        left_join(panelapp%>%
                    filter(panel_name%in%selected_panelapp_panels)%>%
                    select(gene_symbol,panelapp_cat=confidence_level))
      panelapp_genes<-panelapp%>%filter(panel_name%in%selected_panelapp_panels)%>%pull(gene_symbol)
      # add panelapp genes not in disorder to gene genes
      panelapp_genes_not_in_disorder_to_gene_table<-panelapp%>%
        filter(panel_name%in%selected_panelapp_panels)%>%
        filter(!gene_symbol%in%disorder_to_gene_for_table$gene_symbol)%>%
        mutate(disease_id=glue('PANELAPP:{panel_id}'),
               disease_name=glue('{panel_name}'),
               likelihood=as.numeric(confidence_level)-2.5,
               confidence=as.numeric(confidence_level),
               panelapp_cat=as.numeric(confidence_level))%>%
        select(disease_id,gene_symbol,disease_name,likelihood,confidence)
      disorder_to_gene_for_table<-disorder_to_gene_for_table%>%filter(gene_symbol%in%panelapp_genes)
    }
    # filter by user selected threshold
    likely_disorders<-disorders_with_likelihood%>%
      filter(likelihood>input$likelihood_threshold)%>%select(-genes)
    # if no selecte panelapp panel
    if (is.null(panelapp_genes_not_in_disorder_to_gene_table)){
      panel_genes_table<-disorder_to_gene_for_table%>%
        inner_join(likely_disorders)%>%
        mutate(disease_id_name=glue('{disease_id}:{disease_name}'))%>%
        group_by(gene_symbol)%>%
        summarize(disorders=paste0(disease_id_name,collapse=' | '),
                  likelihood=round(max(likelihood),3),
                  confidence=round(max(confidence),3))
    }else{
    # if a panelapp panel was selected, add the panelapp category to each gene
      panel_genes_table<-disorder_to_gene_for_table%>%
        inner_join(likely_disorders)%>%
        bind_rows(panelapp_genes_not_in_disorder_to_gene_table%>%
                    filter(likelihood>input$likelihood_threshold))%>%
        mutate(disease_id_name=glue('{disease_id}:{disease_name}'))%>%
        group_by(gene_symbol)%>%
        summarize(disorders=paste0(disease_id_name,collapse=' | '),
                  likelihood=round(max(likelihood),3),
                  confidence=round(max(confidence),3),
                  panelapp_cat=max(panelapp_cat))
    }
    
    
    panel_genes_table
  })
  
  get_filename<-function(suffix){
    selected_phenos <- input$main_pheno
    selected_panelapp<-input$panelapp_panel
    if (is.null(selected_panelapp)){
      normalized_phenos <- tolower(stringr::str_replace_all(make.names(paste0(selected_phenos,collapse='_')),'HP.\\d+.',''))
    }else{
      normalized_phenos <- tolower(stringr::str_replace_all(make.names(paste0(selected_panelapp,collapse='_')),'HP.\\d+.',''))
    }
    glue('{normalized_phenos}_{Sys.Date()}_{suffix}')
  }
  # server: output - download disorders table button ####
  output$download_disorders <- downloadHandler(
    filename = get_filename('disorders'),
    content = function(file) {
      data <- disorders_data_reactive()
      write.table(data, file,row.names = F,sep = '\t')
    }
  )
  # server: output - download panel genes table button ####
  output$download_panel_genes <- downloadHandler(
    filename = get_filename('genes'),
    content = function(file) {
      data <- panel_genes_reactive()
      write.table(data, file,row.names = F,sep = '\t')
    }
  )
  # output - disorders table ####
  output$disorders_table <- DT::renderDataTable(
    disorders_data_reactive(),
    extensions = 'Buttons',
    options=list(filter = 'top')
  )
  
  # output - panel genes table ####
  output$panel_genes_table<-renderDataTable(
    panel_genes_reactive(),
    extensions = 'Buttons',
    options=list(filter = 'top')
  )
}



shinyApp(ui, server)
