
version = 0.1
library(shiny)
library(shinyWidgets)
library(ProjectTemplate)
library(ontologyIndex)
setwd('..')
data(hpo)
load.project()

# prepare the data and save it
#full_hpoa_df<-parse_hpo_hpoa_db()
#save(full_hpoa_df,file='/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/preprocessed_data.RData')

# Read the test table
load('/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/preprocessed_data.RData')
disorder_to_gene<-create_disorder_to_gene_table()
# full_hpoa_df <- read.delim('/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/test_table.csv')
# all_disorders <- full_hpoa_df %>% pull(disease_id) %>% unique()
ui <- fluidPage(
  theme = bslib::bs_theme(base_font = "Arial", bg = "white", fg = "black"),
  tags$head(
    tags$style(HTML(".main-header {text-align: center;} .subtitle {text-align: center;} .bootstrap-select.btn-group .dropdown-toggle { width: 300px !important; }"))
  ),
  titlePanel("Phenotype-questionnaire Gene Panel Generator"),
  br(),
  fluidRow(
    column(10, offset = 0,
           selectizeInput("main_pheno", "Select the main phenotype:",
                          choices = NULL,
                          multiple = T, width = '500px')
    ),
    # Added new selectizeInput for rejected phenotypes
    column(10, offset = 0,
           selectizeInput("additional_pheno", "Additional phenotypes present:",
                          choices = NULL,
                          multiple = T, width = '500px'),
    ),
    # Added new selectizeInput for rejected phenotypes
    column(10, offset = 0,
           selectizeInput("rejected_pheno", "Select the rejected phenotype:",
                          choices = NULL,
                          multiple = T, width = '500px')
    )
  ),
  br(),
  fluidRow(
    column(4, offset = 4,
           actionButton(inputId = "start_analysis", label = "Start Analysis", class = "btn-primary")
    )
  ),
  br(),
  fluidRow(
    column(8, offset = 0,
           #tags$h4("HPOs by Count"),
           tableOutput("hpos_by_count_table")
    )
  ),
  br(),
  fluidRow(
    column(10, offset = 0,
           dataTableOutput("disorders_table") # Added output for the disorders table
    )
  ),
  br()
)

server <- function(input, output, session) {
  main_phenos_reactive <- reactiveVal()
  rejected_phenos <- reactiveVal(c())
  additional_phenos_reactive <- reactiveVal(c())
  maybe_phenos_reactive<-reactiveVal(c())
  rejected_disorders <- reactiveVal(c())
  top_obligate_pheno <- reactiveVal()
  
  observe({
    updateSelectizeInput(session, "main_pheno", choices = unique(full_hpoa_df$hpo_id_name),server=TRUE)
  })
  
  observeEvent(input$start_analysis, {
    main_phenos <- input$main_pheno
    main_phenos_reactive(main_phenos)
  })
  
  all_phenos_data <- reactive({
    main_phenos <- main_phenos_reactive()
    
    if (is.null(main_phenos)) return(NULL)
    
    all_disorders_with_main_phenos <- full_hpoa_df %>%
      filter(hpo_id_name %in% main_phenos) %>%
      pull(disease_id) %>%
      unique()
    
    all_phenos_from_disorders_with_main_phenos <- full_hpoa_df %>%
      filter(disease_id %in% all_disorders_with_main_phenos)
    
    list(
      main_phenos = main_phenos,
      all_phenos_from_disorders_with_main_phenos = all_phenos_from_disorders_with_main_phenos#%>%filter(!is_ancestor)
    )
  })
  
  # Reactive expression for hpos_by_count
  hpos_by_count_reactive <- reactive({
    phenos_data <- all_phenos_data()
    if (is.null(phenos_data)) return(NULL)
    main_phenos_ancestors<-NULL
    for (main_pheno_id_term in phenos_data$main_phenos){
      main_pheno_id<-stringr::str_extract(main_pheno_id_term,'HP:\\d+')
      main_pheno_ancestors<-get_hpo_ancestors_by_id(main_pheno_id,with_term = T)
      main_phenos_ancestors<-unique(c(main_phenos_ancestors,main_pheno_ancestors))
    }
    # Modify to exclude phenos in main and additional_pheno inputs
    excluded_phenos <- c(main_phenos_ancestors, additional_phenos_reactive(),maybe_phenos_reactive())
    hpos_from_disorders_that_are_not_in_main <- phenos_data$all_phenos_from_disorders_with_main_phenos %>%
      filter(!(disease_id %in% rejected_disorders())) %>%
      filter(!(hpo_id_name %in% excluded_phenos))
    
    # for each ancestor, collect all final descendants
    final_descendants<-hpos_from_disorders_that_are_not_in_main%>%
      filter(is_ancestor==TRUE)%>%
      group_by(hpo_id_name)%>%
      summarize(final_desc=paste0(unique(ancestor_of)[order(unique(ancestor_of))],collapse='<br>'))

    # final_descendants<-NULL # 
    hpos_by_count <- hpos_from_disorders_that_are_not_in_main %>%
      filter(frequency_cat == 'obligate') %>%
      group_by(hpo_id_name) %>%
      summarize(n = n())%>%
      arrange(desc(n))
    
    list(hpos_by_count=hpos_by_count,
         final_descendants=final_descendants)
  })

  # Observer to update the choices of rejected phenotypes selectizeInput
  observe({
    rejected <- rejected_phenos()
    rejected_dis <- rejected_disorders()
    all_phenos <- all_phenos_data()$all_phenos_from_disorders_with_main_phenos$hpo_id_name
    updateSelectizeInput(session, inputId = 'rejected_pheno', choices = unique(all_phenos), selected = rejected)
  })
  
  observe({
    phenos_data <- all_phenos_data()
    if (is.null(phenos_data)) return()
    hpos_by_count_reactive_output<-hpos_by_count_reactive()
    hpos_by_count<-hpos_by_count_reactive_output$hpos_by_count

    num_of_phenos_that_are_obligated <- nrow(hpos_by_count)
    top_pheno <- hpos_by_count %>%
      ungroup() %>%
      filter(!(hpo_id_name %in% phenos_data$main_phenos)) %>%
      slice_max(order_by = n, n = 1, with_ties = F)
    
    top_obligate_pheno(top_pheno)
    print(top_obligate_pheno())
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
    if (num_of_phenos_that_are_obligated > 0) {
      output$hpos_by_count_table <- renderTable({
        hpos_by_count_reactive()$hpos_by_count
      }, rownames = FALSE)
    } else {
      output$hpos_by_count_table <-NULL
      # output$hpos_by_count_table <- renderTable({
      #   hpo_combinations_reactive()
      #}, rownames = FALSE)
      
      return()
    }
    # Render the table for hpos_by_count
    # output$hpos_by_count_table <- renderTable({
    #   hpos_by_count_reactive()
    # }, rownames = FALSE)
    print(top_pheno_final_descendants_text)
    question_text<-ifelse(top_pheno_final_descendants_text=='' | nchar(top_pheno_final_descendants_text)>500,
                          glue('There are {nrow(hpos_by_count)} more phenotypes to ask about.<br>Does your patient have {top_pheno %>% pull(hpo_id_name)}'),
                          glue("There are {nrow(hpos_by_count)} more phenotypes to ask about.<br>Does your patient have {top_pheno %>% pull(hpo_id_name)}<br>Specifically:<br>{top_pheno_final_descendants_text}?"))
    print(question_text)
    showModal(modalDialog(
      title = "Confirm Phenotype",
      #glue("There are {nrow(hpos_by_count)} more phenotypes to ask about.\nDoes your patient have {top_pheno %>% pull(hpo_id_name)}\nSpecifically:\n{top_pheno_final_descendants}?"),
      HTML(question_text),
      footer = tagList(
        actionButton("yes_button", "Yes"),
        actionButton("maybe_button", "Maybe"),
        actionButton("no_button", "No")
      )
    ))
  })
  
  observeEvent(input$yes_button, {
    removeModal()
    additional_phenos <- c(additional_phenos_reactive(), top_obligate_pheno()$hpo_id_name)
    additional_phenos_reactive(additional_phenos)
    updateSelectizeInput(session, inputId = 'additional_pheno', choices = unique(full_hpoa_df$hpo_id_name),selected = additional_phenos,server=TRUE) # Update additional_pheno input
  })
  
  observeEvent(input$no_button, {
    removeModal()
    phenos_data <- all_phenos_data()
    rejected_phenos_current <- c(rejected_phenos(), top_obligate_pheno()$hpo_id_name)
    rejected_phenos(rejected_phenos_current)
    
    rejected_disorders_current <- c(rejected_disorders(), phenos_data$all_phenos_from_disorders_with_main_phenos %>%
                                      filter(frequency_cat == 'obligate') %>%
                                      filter(hpo_id_name %in% rejected_phenos_current) %>%
                                      pull(disease_id))
    
    rejected_disorders(rejected_disorders_current)
  })
  
  observeEvent(input$maybe_button, {
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
      left_join(disorder_to_gene%>%group_by(disease_name)%>%dplyr::summarize(genes=paste0(gene_symbol,collapse=', ')))
    
    return(disorders)
  }
  
  # Reactive to prepare the data for the disorders table
  disorders_data_reactive <- reactive({
    all_phenos_data <- all_phenos_data()
    rejected_dis <- rejected_disorders()
    
    if (is.null(all_phenos_data) || is.null(rejected_dis)) return(NULL)
    
    disorders_data <- create_disorders_table(all_phenos_data, rejected_dis)
    disorders_data
  })
  
  # Render the disorders table
  output$disorders_table <- renderDataTable({
    disorders_data_reactive()
  })
  
}



shinyApp(ui, server)
