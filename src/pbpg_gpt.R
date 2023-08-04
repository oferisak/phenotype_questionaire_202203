library(shiny)
library(shinyWidgets)

# Read the test table
full_hpoa_df <- read.delim('/media/SSD/Bioinformatics/Projects/phenotype_questionaire_202203/data/test_table.csv')
all_disorders <- full_hpoa_df %>% pull(DatabaseID) %>% unique()
ui <- fluidPage(
  theme = bslib::bs_theme(base_font = "Arial", bg = "white", fg = "black"),
  tags$head(
    tags$style(HTML(".main-header {text-align: center;} .subtitle {text-align: center;} .bootstrap-select.btn-group .dropdown-toggle { width: 300px !important; }"))
  ),
  titlePanel("Phenotype-questionnaire Gene Panel Generator"),
  tags$div(class = "subtitle", tags$h3("Select the main Phenotypes:")),
  br(),
  fluidRow(
    column(10, offset = 0,
           selectizeInput("main_pheno", "Select the main phenotype:",
                          choices = unique(test_table$HPO_ID_TERM),
                          multiple = T, width = '500px')
    ),
    # Added new selectizeInput for rejected phenotypes
    column(10, offset = 0,
           selectizeInput("additional_pheno", "Additional phenotypes present:",
                          choices = unique(test_table$HPO_ID_TERM),
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
    column(8, offset = 2,
           tags$h4("HPOs by Count"),
           tableOutput("hpos_by_count_table")
    )
  ),
  br()
)

server <- function(input, output, session) {
  main_phenos_reactive <- reactiveVal()
  rejected_phenos <- reactiveVal(c())
  additional_phenos_reactive <- reactiveVal(c())
  rejected_disorders <- reactiveVal(c())
  top_obligate_pheno <- reactiveVal()
  
  observeEvent(input$start_analysis, {
    main_phenos <- input$main_pheno
    main_phenos_reactive(main_phenos)
  })
  
  all_phenos_data <- reactive({
    main_phenos <- main_phenos_reactive()
    
    if (is.null(main_phenos)) return(NULL)
    
    all_disorders_with_main_phenos <- full_hpoa_df %>%
      filter(HPO_ID_TERM %in% main_phenos) %>%
      pull(DatabaseID) %>%
      unique()
    
    all_phenos_from_disorders_with_main_phenos <- full_hpoa_df %>%
      filter(DatabaseID %in% all_disorders_with_main_phenos)
    
    list(
      main_phenos = main_phenos,
      all_phenos_from_disorders_with_main_phenos = all_phenos_from_disorders_with_main_phenos
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
    excluded_phenos <- c(main_phenos_ancestors, additional_phenos_reactive())
    hpos_from_disorders_that_are_not_in_main <- phenos_data$all_phenos_from_disorders_with_main_phenos %>%
      filter(!(DatabaseID %in% rejected_disorders())) %>%
      filter(!(HPO_ID_TERM %in% excluded_phenos))
    
    hpos_by_count <- hpos_from_disorders_that_are_not_in_main %>%
      filter(Frequency_cat == 'obligate') %>%
      group_by(HPO_ID_TERM) %>%
      summarize(n = n())%>%
      arrange(desc(n))
    
    hpos_by_count
  })
  
  # Observer to update the choices of rejected phenotypes selectizeInput
  observe({
    rejected <- rejected_phenos()
    rejected_dis <- rejected_disorders()
    all_phenos <- all_phenos_data()$all_phenos_from_disorders_with_main_phenos$HPO_ID_TERM
    updateSelectizeInput(session, inputId = 'rejected_pheno', choices = unique(all_phenos), selected = rejected)
  })
  
  observe({
    phenos_data <- all_phenos_data()
    if (is.null(phenos_data)) return()

    hpos_by_count<-hpos_by_count_reactive()
    num_of_phenos_that_are_obligated <- nrow(hpos_by_count)
    top_pheno <- hpos_by_count %>%
      ungroup() %>%
      filter(!(HPO_ID_TERM %in% phenos_data$main_phenos)) %>%
      slice_max(order_by = n, n = 1, with_ties = F)
    
    top_obligate_pheno(top_pheno)
    
    if (num_of_phenos_that_are_obligated == 0 || nrow(top_pheno) == 0) {
      message("Done collecting..")
      return()
    }
    
    # Render the table for hpos_by_count
    output$hpos_by_count_table <- renderTable({
      hpos_by_count_reactive()
    }, rownames = FALSE)
    
    showModal(modalDialog(
      title = "Confirm Phenotype",
      glue("There are {nrow(hpos_by_count)} more phenotypes to ask about.\nDoes your patient have {top_pheno %>% pull(HPO_ID_TERM)}?"),
      footer = tagList(
        actionButton("yes_button", "Yes"),
        actionButton("no_button", "No")
      )
    ))
  })
  
  observeEvent(input$yes_button, {
    removeModal()
    additional_phenos <- c(additional_phenos_reactive(), top_obligate_pheno()$HPO_ID_TERM)
    additional_phenos_reactive(additional_phenos)
    updateSelectizeInput(session, inputId = 'additional_pheno', selected = additional_phenos) # Update additional_pheno input
  })
  
  observeEvent(input$no_button, {
    removeModal()
    phenos_data <- all_phenos_data()
    rejected_phenos_current <- c(rejected_phenos(), top_obligate_pheno()$HPO_ID_TERM)
    rejected_phenos(rejected_phenos_current)
    
    rejected_disorders_current <- c(rejected_disorders(), phenos_data$all_phenos_from_disorders_with_main_phenos %>%
                                      filter(Frequency_cat == 'obligate') %>%
                                      filter(HPO_ID_TERM %in% rejected_phenos_current) %>%
                                      pull(DatabaseID))
    
    rejected_disorders(rejected_disorders_current)
  })
  
  output$disorders_table <- DT::renderDataTable({
    selected_disorders
  }, rownames = FALSE)
  
  output$selected_phenotypes <- renderTable({
    data.frame('Phenotypes' = input$phenotype)
  }, rownames = FALSE)
}



shinyApp(ui, server)
