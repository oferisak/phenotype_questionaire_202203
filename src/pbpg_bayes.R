
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

# TODO: Add for each gene - what are the supporting and non supporting phenotypes given

# prepare the data and save it
#full_hpoa_df<-parse_hpo_hpoa_db()
#save(full_hpoa_df,file='/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/preprocessed_data.RData')

# Tried: 
## modify the questions order to be according to the total frequency in the different disorders (more frequent will be asked before) - 
##  - did not result in a major improvement in the number of questions required

# Read the test table
load('/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/preprocessed_data.RData')
disorder_to_gene<-create_disorder_to_gene_table()
# Read panelapp panels
message('Parsing panelapp panels..')
panelapp<-readr::read_delim('./data/panelapp_db_2023-08-17.csv.gz')

bayes_freq_for_unknown<-0.3

message(glue('SETUP: defined frequency for unknown phenotypes: {bayes_freq_for_unknown}'))

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
                                        value=0.1, max = 1, min=0, step = 0.01),
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
                 dataTableOutput("disorders_table")
               )
             )
    ),
    # ui: tab - panel genes ####
    tabPanel("Panel genes",
             mainPanel(
               br(),
               fluidRow(
                 column(10,sliderInput("likelihood_threshold", "Likelihood Threshold:", 0.1, min = 0, max = 1, step = 0.01))
               ),
               br(),
               fluidRow(
                 column(10,downloadButton("download_panel_genes", "Download Panel Genes Table"))
               ),
               br(),
               fluidRow(
                 dataTableOutput('panel_genes_table')
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
  initial_disease_priors<-reactiveVal(c())
  added_pheno<-reactiveVal()
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
  
  # server: observer: start analysis button ####
  observeEvent(input$start_analysis, {
    updateActionButton(session, "start_analysis", label = "Resume Analysis")
    main_phenos <- input$main_pheno
    main_phenos_reactive(main_phenos)
    # initial_disease_priors_vector<-calc_initial_disease_priors()
    # initial_disease_priors(initial_disease_priors_vector)
    #message(glue('there are {length(initial_disease_priors_vector)} intial diseases'))
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
  # server: reactive - intial_disease_priors ####
  calc_initial_disease_priors<-reactive({
    message('Running reactive: calc_initial_disease_priors')
    phenos_data <- all_phenos_data()
    print(phenos_data)
    initial_diseases<-phenos_data$all_phenos_from_disorders_with_main_phenos%>%
      pull(disease_id)%>%unique()
    initial_diseases_priors<-data.frame(disease_id=initial_diseases,prior_prob=rep(1/length(initial_diseases),length(initial_diseases)))
    initial_diseases_priors
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
      filter(disease_id %in% all_disorders_with_main_phenos)%>%
      mutate(frequency_bayes=ifelse(frequency_cat=='unknown',bayes_freq_for_unknown,frequency_bayes)) 
    
    complete_grid <- expand.grid(
      disease_id = unique(all_phenos_from_disorders_with_main_phenos$disease_id),
      hpo_id_name = unique(all_phenos_from_disorders_with_main_phenos$hpo_id_name)
    )
    complete_grid<-complete_grid%>%mutate(disease_id=as.character(disease_id),hpo_id_name=as.character(hpo_id_name))
    all_phenos_from_disorders_with_main_phenos<-complete_grid%>%left_join(all_phenos_from_disorders_with_main_phenos)%>%
      mutate(frequency_bayes=ifelse(is.na(frequency_bayes),0.05,frequency_bayes))
    
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
    
    disorders_with_likelihood<-update_disorders_likelihood_bayes()
    # set the threshold so that all the disorders below the analysis_depth quantile will be considered unlikely and will not be asked about
    likelihood_thresh<-quantile(disorders_with_likelihood$likelihood,probs=input$analysis_depth)
    message(glue('likelihood thresh: {likelihood_thresh}'))
    disorders_remaining<-disorders_with_likelihood%>%filter(likelihood>=likelihood_thresh)%>%pull(disease_id)
    unlikely_disorders<-disorders_with_likelihood%>%filter(likelihood<likelihood_thresh)%>%pull(disease_id)
    message(glue('There are currently {length(unlikely_disorders)}/{nrow(disorders_with_likelihood)} unlikely disorders..'))
    # Modify to exclude phenos in main and additional_pheno inputs
    #excluded_phenos <- c(main_phenos_ancestors, additional_phenos_present_reactive(),maybe_phenos_reactive(),rejected_phenos())
    hpos_from_disorders_that_are_not_in_main <- phenos_data$all_phenos_from_disorders_with_main_phenos %>%
      #filter(!(disease_id %in% rejected_disorders())) %>%
      filter(!(disease_id %in% unlikely_disorders)) %>%
      filter(!(hpo_id_name %in% excluded_due_to_addition_ancestors))%>%
      filter(!(hpo_id_name %in% excluded_due_to_addition_descendants))%>%
      filter(!(hpo_id_name %in% excluded_due_to_rejection_descendants))%>%
      filter(!frequency_cat=='unknown')
    
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
    
    top_pheno_final_descendants<-hpos_by_count_reactive_output$final_descendants
    # if the hpo in question is not an ancestor of anything
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
    added_pheno(top_obligate_pheno()$hpo_id_name)
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
  
  
  # Disorder table functions ####
  
  create_disorders_table <- function(all_phenos_data, rejected_disorders) {

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
  library(dplyr)
  
  # Function to update disease probabilities based on user inputs
  update_disorders_likelihood_bayes <- function() {
    all_phenos_data <- all_phenos_data()
    priors<-calc_initial_disease_priors()
    message(glue('Initial prior probs: {priors%>%pull(prior_prob)%>%unique()}'))
    main_phen<-main_phenos_reactive()
    additional_phen<-additional_phenos_present_reactive()
    selected_phenotypes<-c(main_phen,additional_phen)
    rejected_phenotypes <- rejected_phenos()
    # Add a column to indicate whether a phenotype is selected, rejected, or unspecified
    phenotype_disease_table <- all_phenos_data$all_phenos_from_disorders_with_main_phenos %>%
      mutate(status = case_when(
        hpo_id_name %in% selected_phenotypes ~ 'selected',
        hpo_id_name %in% rejected_phenotypes ~ 'rejected',
        TRUE ~ 'unspecified'
      ))
    
    # Process selected phenotypes
    selected_updates <- phenotype_disease_table %>%
      filter(status == 'selected') %>%
      group_by(disease_id) %>%
      summarize(product_selected = prod(frequency_bayes), .groups = 'drop')
      #summarize(product_selected = prod(ifelse(status=='selected',frequency_bayes,0.05)), .groups = 'drop')
    
    # Process rejected phenotypes
    rejected_updates <- phenotype_disease_table %>%
      filter(status == 'rejected') %>%
      group_by(disease_id) %>%
      summarize(product_rejected = prod(1-frequency_bayes), .groups = 'drop')
    
    # Combine updates and include diseases with unspecified phenotypes (which will not affect the calculation)
    combined_updates <- phenotype_disease_table %>% left_join(priors)%>%
      select(disease_id,prior_prob) %>%
      distinct() %>%
      left_join(selected_updates, by = "disease_id") %>%
      left_join(rejected_updates, by = "disease_id") %>%
      mutate(product_selected = replace_na(product_selected, 1),
             product_rejected = replace_na(product_rejected, 1),
             updated_probability = prior_prob* product_selected * product_rejected)
    
    # Normalize probabilities so they sum to 1
    combined_updates <- combined_updates %>%
      mutate(likelihood = updated_probability / sum(updated_probability),confidence=updated_probability / sum(updated_probability))%>%
      left_join(disorder_to_gene%>%group_by(disease_id,disease_name)%>%
                  summarize(genes=paste0(gene_symbol,collapse=','))%>%ungroup())%>%relocate(disease_name)
    message(glue('TOP DISORDERS:'))
    print(combined_updates%>%slice_max(n=3,order_by=likelihood)%>%select(disease_name,likelihood))
    return(combined_updates)
  }
  
  # Reactive to prepare the data for the disorders table
  disorders_data_reactive <- reactive({
    message('Running reactive: disorders_data_reactive')
    all_phenos_data <- all_phenos_data()
    
    if (is.null(all_phenos_data)) return(NULL)
    disorders_with_likelihood<-update_disorders_likelihood_bayes()  
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
  
  # Generate the panel genes table ####
  panel_genes_reactive<-reactive({
    message('Running reactive: panel_genes_reactive')
    main_phen<-main_phenos_reactive()
    additional_phen<-additional_phenos_present_reactive()
    selected_phenotypes<-c(main_phen,additional_phen)
    rejected_phenotypes <- rejected_phenos()
    
    if (is.null(all_phenos_data)) return(NULL)
    disorders_with_likelihood<-update_disorders_likelihood_bayes()
    all_phenos_data <- all_phenos_data()
    # get the selected phenotypes and their associated disorders
    selected_phenos_for_genes<-
      all_phenos_data$all_phenos_from_disorders_with_main_phenos%>%
      filter(hpo_id_name %in% c(selected_phenotypes,rejected_phenotypes))%>%
      mutate(frequency_cat=ifelse(is.na(frequency_cat),'not_characteristic',frequency_cat),
             hpo_id_name_freq=glue('{hpo_id_name}({frequency_cat})'))
    print(selected_phenos_for_genes)
    selected_panelapp_panels<-panelapp_panels_reactive()
    disorder_to_gene_for_table<-disorder_to_gene
    disorder_to_gene_for_table<-disorder_to_gene_for_table%>%left_join(selected_phenos_for_genes)
    panelapp_genes_not_in_disorder_to_gene_table<-NULL
    if (!is.null(selected_panelapp_panels)){
      # first add the panelapp confidence to all the genes
      disorder_to_gene_for_table<-disorder_to_gene_for_table%>%
        left_join(panelapp%>%
                    filter(panel_name%in%selected_panelapp_panels)%>%
                    select(gene_symbol,panelapp_cat=confidence_level))
      panelapp_genes<-panelapp%>%filter(panel_name%in%selected_panelapp_panels)%>%pull(gene_symbol)
      # add panelapp genes not in disorder to panel genes
      panelapp_genes_not_in_disorder_to_gene_table<-panelapp%>%
        filter(panel_name%in%selected_panelapp_panels)%>%
        filter(!gene_symbol%in%disorder_to_gene_for_table$gene_symbol)%>%
        mutate(disease_id=glue('PANELAPP:{panel_id}'),
               disease_name=glue('{panel_name}'),
               likelihood=as.numeric(confidence_level),
               confidence=as.numeric(confidence_level)-2.5,
               panelapp_cat=as.numeric(confidence_level))%>%
        select(disease_id,gene_symbol,disease_name,likelihood,confidence)
      disorder_to_gene_for_table<-disorder_to_gene_for_table%>%filter(gene_symbol%in%panelapp_genes)
    }
    # filter by user selected threshold
    likely_disorders<-disorders_with_likelihood%>%
      filter(likelihood>quantile(likelihood,probs=input$likelihood_threshold))%>%select(-genes)
    
    # if no selected panelapp panel
    if (is.null(panelapp_genes_not_in_disorder_to_gene_table)){
      panel_genes_table<-disorder_to_gene_for_table%>%
        inner_join(likely_disorders)%>%
        mutate(disease_id_name=glue('{disease_id}:{disease_name}'))%>%
        group_by(gene_symbol,hpo_id_name)%>%
        slice_max(n=1,with_ties = F,order_by = frequency_bayes)%>%ungroup()%>%group_by(gene_symbol)%>%
        summarize(disorders=paste0(unique(disease_id_name),collapse=' | '),
                  likelihood=round(max(likelihood),5),
                  confidence=round(max(confidence),5),
                  phenotypes_answered=paste0(unique(hpo_id_name_freq),collapse=','))
    }else{
      # if a panelapp panel was selected, add the panelapp category to each gene
      panel_genes_table<-disorder_to_gene_for_table%>%
        inner_join(likely_disorders)%>%
        bind_rows(panelapp_genes_not_in_disorder_to_gene_table%>%
                    filter(likelihood>quantile(likelihood,probs=input$likelihood_threshold)))%>%
        mutate(disease_id_name=glue('{disease_id}:{disease_name}'))%>%
        group_by(gene_symbol)%>%
        summarize(disorders=paste0(unique(disease_id_name),collapse=' | '),
                  likelihood=round(max(likelihood),5),
                  confidence=round(max(confidence),5),
                  supporting=paste0(unique(hpo_id_name_freq),collapse=','),
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
