library(pacman)
## 
p_load(
  here, rio, tidyverse, readr, stringr, colorspace,
  shiny, shinyWidgets, DT, NGLVieweR, plotly, shinyjs, bslib)

source("./Functions/functions_Var_int-nova_v13.R")

##################################################################
## ProtSAVI - Protein Structure-Assisted Variant Interpretation ##
##################################################################

#### Lista de Genes com Dados Pré-Calulados ####
list.genes <- import(here("./Data/List_genes.txt"), header = FALSE)
mane_info <- read.delim("./Data/Mane_links.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------------

addResourcePath(
  "tutorial",
  normalizePath("www", mustWork = TRUE)
)

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "yeti"),
  useShinyjs(),
  
  # ============================================================
  # Importa a fonte Montserrat do Google Fonts
  # ============================================================
  tags$head(
    tags$link(
      rel  = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@400;500;600;700&display=swap"
    )
  ),
  
  # ============================================================
  # Estilo dos SelectInputs (Selectize)
  # ============================================================
  tags$head(
    tags$style(HTML("
    /* Campo principal do selectInput */
    .selectize-input {
      border-radius: 10px !important; /* Bordas arredondadas */
      min-height: 44px;               /* Altura mínima confortável */
      padding: 10px 12px;             /* Espaçamento interno */
      font-size: 14px;                /* Tamanho da fonte */
    }

    /* Dropdown (lista aberta de opções) */
    .selectize-dropdown {
      border-radius: 10px;
    }

    /* Remove sombra padrão do Bootstrap */
    .selectize-control.single .selectize-input {
      box-shadow: none;
    }
  "))
  ),
  
  # ============================================================
  # Estilo da sidebar (.well)
  # ============================================================
  tags$head(
    tags$style(HTML("
    /* Caixa externa da sidebar */
    .well {
      border-radius: 12px;              /* Bordas arredondadas */
      border: 1px solid #E0E0E0;        /* Borda suave */
      box-shadow: none;                /* Remove sombra padrão */
      background-color: #FFFFFF;       /* Fundo branco */
      padding: 25px;                   /* Espaçamento interno */
    }
  "))
  ),
  
  # ============================================================
  # Estilo das abas (Tabs no formato Pills)
  # ============================================================
  tags$head(
    tags$style(HTML("

    /* =========================
       TABS (PILLS)
       ========================= */

    /* Container das tabs */
    .nav-tabs {
      border-bottom: none;
      margin-bottom: 16px;

      display: flex;
      flex-wrap: wrap;
      gap: 12px 16px; /* espaçamento vertical e horizontal */
    }

    /* Tabs padrão (inativas) */
    .nav-tabs > li > a {
      border-radius: 999px !important; /* Formato pill */
      padding: 10px 22px;
      background-color: #F4F6F7;
      color: #4F8F8B;
      border: 1px solid #E0E0E0 !important;
      margin-right: 0px;
      position: relative;
      text-decoration: none;
    }

    /* Remove comportamento clássico do Bootstrap (linha inferior) */
    .nav-tabs > li > a::after,
    .nav-tabs > li.active > a::after {
      display: none !important;
    }

    /* Aba ATIVA */
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:hover,
    .nav-tabs > li.active > a:focus {
      background-color: #FFFFFF;
      color: #1E1E1E;
      border: 1px solid #E0E0E0 !important;
      border-bottom: 1px solid #E0E0E0 !important;
      box-shadow: 0 2px 6px rgba(0,0,0,0.05);
      font-weight: 700 !important;
    }

    /* Hover nas abas inativas */
    .nav-tabs > li:not(.active) > a:hover {
      background-color: #EAF3F3;
      color: #3E6F73;
    }

    /* =========================
       CONTEÚDO DAS ABAS (CARD)
       ========================= */

    .tab-content {
      border: 1px solid #E0E0E0;
      border-radius: 16px;
      padding: 25px;
      background-color: #FFFFFF;
    }

    /* =========================
       CARD LARGO (ABA INFORMATION)
       ========================= */

    .info-card-wide {
      width: 100%;
      max-width: none;
      padding: 30px;
    }

  "))
  ),
  
  # ============================================================
  # Layout das linhas de informação (label + valor)
  # ============================================================
  tags$style(HTML("
  .info-line {
    display: grid;
    grid-template-columns: 200px 1fr; /* coluna fixa + flexível */
    gap: 16px;
    padding: 10px 0;
    border-bottom: 1px solid #eaeaea;
    align-items: start;
  }

  .info-label {
    font-weight: 600;
    color: #333;
    white-space: nowrap;
  }

  .info-value {
    color: #444;
    word-break: break-word;
  }
")),
  
  # ============================================================
  # Ajuste de largura dos inputs
  # ============================================================
  tags$style(HTML("
  /* Ajusta largura dos campos de formulário */
  .form-control {
    display: inline-block;
    width: 90%;
  }
")),
  
  # ============================================================
  # Estilo do rótulo do gene selecionado
  # ============================================================
  tags$style(HTML("
  .gene-label {
    font-family: 'Montserrat', sans-serif;
    font-size: 28px;
    font-weight: 500;
    color: #4F8F8B;
    margin-bottom: 16px;
  }
")),
  
  # ============================================================
  # Fonte global do aplicativo
  # ============================================================
  tags$style(HTML("
  body, 
  h1, h2, h3, h4, h5, h6,
  p, span, label, a,
  input, textarea, select, button,
  .btn {
    font-family: 'Montserrat', sans-serif !important;
  }
")),
  
  ## Alinhas os trem
  
  tags$head(
    tags$style(HTML("
    table.dataTable thead th {
      text-align: center !important;
      vertical-align: middle !important;
    }
  "))
  ),
  
  tags$style(HTML("
  /* Centraliza APENAS o tabsetPanel principal */
  #main_tabs.nav-tabs,
  #main_tabs > .nav-tabs {
    justify-content: center !important;
  }
"))
  ,
  
  titlePanel(
    tagList(
      
      tags$div(
        style = "
        display: flex;
        justify-content: center;
        align-items: center;
        padding: 20px 0 10px 0;
        width: 100%;
      ",
        
        tags$h1(
          "ProtSAVI - Protein Structure-Assisted Variant Interpretation",
          style = "
          font-family: 'Montserrat', sans-serif;
          font-weight: 500;
          color: #4F8F8B;
          font-size: 38px;
          margin: 0;
          text-align: center;
        "
        )
      ),
      
      tags$hr(style = "margin: 10px auto 20px auto; width: 100%;")
    )
  ),
  
  tabsetPanel(
    id = "main_tabs",
    
    ## -------- INTRODUCTION --------
    tabPanel(
      title = "Introduction",
      value = "intro",
      div(
        style = "max-width:900px; margin:auto;",
        includeHTML("www/intro_v3.html")
      )
    ),
    
    ## -------- TUTORIAL --------
    tabPanel(
      title = "Tutorial",
      div(
        style = "max-width:900px; margin:auto;",
        includeHTML("www/VariantGuide_v7.html")
    )),
    
    tabPanel(
      title = "ProtSAVI",
      sidebarLayout(
    sidebarPanel(
      width = 3,
      
      tags$div(
        style = "display: flex; flex-direction: column; align-items: center;",
        
        selectInput(
          inputId = "GENE",
          label = tags$div(class = "gene-label", "Select a Gene"),
          choices = sort(import(here("./Data/list_genes.txt"), header = FALSE)$V1),
          selected = NULL
        )
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'interpret'",
        
        tags$div(
          style = "
    margin-top:10px;
    margin-bottom:15px;
    font-size:16px;
  ",
          
          tags$div(
            style = "display:flex; align-items:center; margin-bottom:6px;",
            HTML("<span style='color:red; font-size:18px; margin-right:8px;'>&#9670;</span>"),
            tags$span("Prioritized criteria:")
          ),
          
          tags$div(
            tags$ul(
              style = "margin-top:4px; margin-bottom:0;",
              tags$li("VUS with ACMG score equal to 3 points, associated with functional information"),
              tags$li("ACMG score = 4 or 5"),
              tags$li("ACMG score = 3 associated with interpretation")
          
        )))),
      
      # Itens só na aba Pathogenicity
      conditionalPanel(
        condition = "input.tabs == 'viewer3d'",
        
        textAreaInput(
          inputId = "Variants",
          label = "Enter Amino Acid variants (comma-separated):",
          placeholder = "Ex:A14C,B15G",
          width = "100%",
          rows = 5
        ),
        actionButton("send_variants", "Send Variants", class = "btn-primary"),

        # Widget de upload recriado via renderUI
        uiOutput("uploadVariantsUI"),
        
        downloadLink(
          outputId = "example_input",
          label = "Example of input format",
          style= "font-size: 14px",
        ),
        
        p("The transcript sequence is in the Information window.",
          style= "font-size: 14px"),
        
        ###
        
        hr(),
        
        # Clinvar Variants
        prettyRadioButtons(
          inputId = "clinvar",
          label = "Clinvar Variants:",
          choices = c(
            "Pathogenic" = 'p',
            "Likely Pathogenic" = 'lp',
            "VUS" = 'vus',
            "Likely Benign" = 'lb',
            "Benign" = 'b',
            "All" = 'all'
          ),
          selected = 'all',
          inline = FALSE,
          fill = TRUE,
          status = "primary",
          animation= "jelly",
          bigger = TRUE
        ),
        tags$style(HTML("
         /* Pathogenic */
          .pretty.p input:checked + label { background-color: #FF0000 !important; color: white; }
           /* Likely Pathogenic */
           .pretty.lp input:checked + label { background-color: #F08080 !important; color: white; }
           /* VUS */
           .pretty.vus input:checked + label { background-color: #00FFFF !important; color: black; }
           /* Likely Benign */
           .pretty.lb input:checked + label { background-color: #B2EDB5 !important; color: black; }
           /* Benign */
           .pretty.b input:checked + label { background-color: #00FF00 !important; color: white; }
           /* All */
           .pretty.all input:checked + label { background-color: #F0F8FF !important; color: black; }
           ")),
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        
        tabPanel(
          title = "Information",
          value = "info",
          
          tags$style(HTML("
          .info-line {
          display: flex;
          align-items: flex-start;
          padding: 10px 0;
          border-bottom: 1px solid #eaeaea;
        }
        
        .info-label {
          width: 180px;
          min-width: 180px;
          font-weight: 600;
          color: #333;
          white-space: normal;
          line-height: 1.3;
        }
        
        .info-value {
          flex: 1;
          color: #444;
        }
        ")),
          
          fluidRow(
            column(6, uiOutput("info_text")),
            column(6, imageOutput("figure"))
          )
        )
        ,
        
        tabPanel(
          title = "Clinvar Statistics",
          value = "stats",
          h3("Statistical Analyses"),
          fluidRow(
            column(6, plotlyOutput("plotb")),
            column(6, plotlyOutput("plotc"))
          ),
          hr(),
          fluidRow(
            column(6, plotlyOutput("plotd")),
            column(6, plotlyOutput("plotf"))
          ),
            hr(),
            fluidRow(
              column(6, plotlyOutput("plotg")),
              column(6, plotlyOutput("ploth"))
            )
        ),
        
        tabPanel(
          title = "In Silico Scores",
          value = "scores",
          h3("In Silico Scores"),
          div(
            style = "width: 100%;",
            DT::dataTableOutput("in_silico_scores")
          ),
          uiOutput("download_table_combined")
        ),
        
        tabPanel(
          title = "Variant Assessment",
          value = "viewer3d",
          conditionalPanel(
            condition = "output.fileUploaded == true",
            fluidRow(
              column(6,
                     h4("Variant Positions", align = 'center'),
                     hr(),
                     NGLVieweROutput("structure", height = "600px"),
                     tags$div(
                       style = "display: flex; align-items: center; font-size: 16px; margin-top: 10px;",
                       tags$div(
                         style = "width: 15px; height: 15px; background-color: yellow; border-radius: 50%; margin-right: 8px;"
                       ),
                       tags$span("Mutated positions")
                     )
              ),
              column(6,
                     h4("ClinVar Variants", align = 'center'),
                     hr(),
                     NGLVieweROutput("structure.2", height = "600px"),
                     # Legenda com 5 círculos
                     tags$div(
                       style = "font-size: 16px; margin-top: 10px;",
                       tags$div(
                         style = "display: flex; flex-direction: column; gap: 4px;",
                         
                         tags$div(
                           style = "display: flex; align-items: center;",
                           tags$div(style = "width: 15px; height: 15px; background-color: #CA0000; border-radius: 50%; margin-right: 8px;"),
                           tags$span("Pathogenic")
                         ),
                         tags$div(
                           style = "display: flex; align-items: center;",
                           tags$div(style = "width: 15px; height: 15px; background-color: #FF9900; border-radius: 50%; margin-right: 8px;"),
                           tags$span("Likely Pathogenic")
                         ),
                         tags$div(
                           style = "display: flex; align-items: center;",
                           tags$div(style = "width: 15px; height: 15px; background-color: #48C9B0; border-radius: 50%; margin-right: 8px;"),
                           tags$span("VUS")
                         ),
                         tags$div(
                           style = "display: flex; align-items: center;",
                           tags$div(style = "width: 15px; height: 15px; background-color: #FCFD9B; border-radius: 50%; margin-right: 8px;"),
                           tags$span("Likely Benign")
                         ),
                         tags$div(
                           style = "display: flex; align-items: center;",
                           tags$div(style = "width: 15px; height: 15px; background-color: #77FC02; border-radius: 50%; margin-right: 8px;"),
                           tags$span("Benign")
                         )
                       )
                     ),
                     style = 'border-left: 1px solid'
              ),
            ),
          ),
          
          conditionalPanel(
            condition = "output.fileUploaded == false",
            #h3("ClinVar variants", align = 'left'),
            p("All ClinVar variants for the selected gene are shown here. Use the left panel to filter by category.", align = 'left'),
            p(tags$b("Waiting for submission"), align = 'left'),
            NGLVieweROutput("structure.3", height = "600px"),
            
            
            # Legenda com 5 círculos
            tags$div(
              style = "font-size: 16px; margin-top: 10px;",
              tags$div(
                style = "display: flex; flex-direction: column; gap: 4px;",
                
                tags$div(
                  style = "display: flex; align-items: center;",
                  tags$div(style = "width: 15px; height: 15px; background-color: #CA0000; border-radius: 50%; margin-right: 8px;"),
                  tags$span("Pathogenic")
                ),
                tags$div(
                  style = "display: flex; align-items: center;",
                  tags$div(style = "width: 15px; height: 15px; background-color: #FF9900; border-radius: 50%; margin-right: 8px;"),
                  tags$span("Likely Pathogenic")
                ),
                tags$div(
                  style = "display: flex; align-items: center;",
                  tags$div(style = "width: 15px; height: 15px; background-color: #48C9B0; border-radius: 50%; margin-right: 8px;"),
                  tags$span("VUS")
                ),
                tags$div(
                  style = "display: flex; align-items: center;",
                  tags$div(style = "width: 15px; height: 15px; background-color: #FCFD9B; border-radius: 50%; margin-right: 8px;"),
                  tags$span("Likely Benign")
                ),
                tags$div(
                  style = "display: flex; align-items: center;",
                  tags$div(style = "width: 15px; height: 15px; background-color: #77FC02; border-radius: 50%; margin-right: 8px;"),
                  tags$span("Benign")
                )
              )
            )
          ),
        ),
        
        tabPanel(
          title = "Variant Interpretation",
          value = "interpret",
          h3("Analysis of Variants"),
          DT::dataTableOutput("resumo"),
          br(),
          uiOutput("download_table")
        ),
        
        tabPanel(
          title = "Structural Descriptors",
          value = "interpret_full",
          #h3("Full Dataset for Variant Interpretation"),
          DT::dataTableOutput("resumo_full"),
          br(),
          uiOutput("download_table_full")
        )
      )
    )
  )),
  
  tabPanel(
    title = "About us",
    value = "aboutus",
    div(
      style = "max-width:900px; margin:auto;",
      includeHTML("www/aboutus_v2.html")
    )
  )),
  
  #tags$hr(style = "margin-top: 20px; margin-bottom: 20px;"),
  
  tags$footer(
    style = "
    width: 100%;
    text-align: center;
    font-size: 0.9em;
    color: #555;
    margin-top: 10px;
  ",
    
    # LOGOS
    tags$div(
      style = "
      display: flex;
      justify-content: center;
      align-items: center;
      gap: 25px;
      max-width: 800px;
      margin: 0 auto;
      padding: 5px 10px;
    ",
      
      #imageOutput("logo_fiocruz", height = "65px")
      tags$div(
        style = "position: relative; top: 3px;",
        imageOutput("logo_fiocruz", height = "60px")
      ),
      
      tags$div(
        style = "position: relative; top: -10px;",
        imageOutput("logo_ufrj", height = "60px")
      ),
      
      tags$div(
        style = "margin-left: -5px;",
        imageOutput("logo_extra", height = "55px")
      ),
      #imageOutput("logo_renomieii", height = "60px")
      tags$div(
        style = "position: relative; top: 16px;margin-left: -15px;",
        imageOutput("logo_renomieii", height = "80px")
    )),
    
    # LINHA
    tags$hr(style = "margin-top: 20px; margin-bottom: 20px;"),
    
    # TEXTO
    tags$div(
      style = "margin-bottom: 30px;",
      
      "Last Update: 04/2026",
      tags$br(),
      
      tags$a(
        href = "mailto:maucosta@gmail.com",
        "📧 Contact Us",
        style = "margin-right: 12px;"
      ),
      
      tags$span("•"),
      
      tags$a(
        href = "https://github.com/paulafranklinn",
        "🐞 Report an issue",
        target = "_blank",
        style = "margin-left: 12px;"
      )
    )
  )
)

server <- function(input, output, session) {

  output$logo_fiocruz <- renderImage({
    list(
      src = "www/fiocruz_.jpg",
      width = "180px",
      alt = "Logo Fiocruz"
    )
  }, deleteFile = FALSE)
  
  output$logo_renomieii <- renderImage({
    list(
      src = "www/logo_renomieii.png",
      width = "130px",
      alt = "Logo Renomieii"
    )
  }, deleteFile = FALSE)
  
  output$logo_ufrj <- renderImage({
    list(
      src = "www/logo_ufrj.png",
      width = "180px",
      alt = "Logo UFRJ"
    )
  }, deleteFile = FALSE)
  
  output$logo_extra <- renderImage({
    list(
      src = "www/logo_iff.jpg",
      width = "120px",
      alt = "Logo Extra"
    )
  }, deleteFile = FALSE)
  
  
  # Reactive values para armazenar os dados do upload e do text area
  uploaded_file <- reactiveVal(NULL)
  text_input_data <- reactiveVal(NULL)
  
  # 1) Recria fileInput sempre que trocar de gene 
  output$uploadVariantsUI <- renderUI({
    req(input$GENE)
    div(id = "fileDiv",
        fileInput("file1", "or Upload Amino Acid Variants", multiple = FALSE, accept = c(".txt", ".csv")),
    )
  })
  
  output$file1 <- renderDataTable({
    # Certifica-se de que all_variants() tem dados antes de criar o data frame
    req(all_variants())
    data.frame(Variante = all_variants(), stringsAsFactors = FALSE)
  })
  
  # 2) Indicador de que dados foram submetidos (arquivo OU texto)
  output$fileUploaded <- reactive({
    !is.null(uploaded_file()) || !is.null(text_input_data())
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  # 3) Esconde abas no início
  hideTab("tabs", "scores")
  hideTab("tabs", "interpret")
  hideTab("tabs", "interpret_full")
  
  # 4) Observador para o UPLOAD DE ARQUIVO
  observeEvent(input$file1, {
    req(input$file1)
    uploaded_file(input$file1)
    text_input_data(NULL) # Limpa o input de texto se um arquivo for enviado
    
    showTab("tabs", "scores")
    showTab("tabs", "viewer3d")
    showTab("tabs", "interpret")
    showTab("tabs", "interpret_full")
    updateTabsetPanel(session, "tabs", selected = "viewer3d")
  })
  
  # Observador para o INPUT DE TEXTO
  observeEvent(input$send_variants, {
    req(input$Variants)
    
    # Processa o texto da caixa de texto
    variants_raw <- unlist(strsplit(input$Variants, ",\\s*|\\s*,", perl = TRUE))
    variants_clean <- trimws(variants_raw)
    variants_clean <- variants_clean[variants_clean != ""]
    
    # Armazena nas reativas
    text_input_data(variants_clean)
    uploaded_file(NULL) # Limpa o arquivo se o texto for enviado
    
    # Limpa o textAreaInput
    updateTextAreaInput(session, "Variants", value = "")
    
    # Mostra abas
    showTab("tabs", "scores")
    showTab("tabs", "viewer3d")
    showTab("tabs", "interpret")
    showTab("tabs", "interpret_full")
    updateTabsetPanel(session, "tabs", selected = "viewer3d")
  })
  
  
  # 5) Ao trocar de gene, esconde abas de novo e zera tudo
  observeEvent(input$GENE, {
    # Reset abas
    hideTab("tabs", "scores")
    hideTab("tabs", "interpret")
    hideTab("tabs", "interpret_full")
    updateTabsetPanel(session, "tabs", selected = "info")
    
    # Reset fileInput e text area
    shinyjs::reset("fileDiv")
    updateTextAreaInput(session, "Variants", value = "")
    
    # Reset reactiveVals
    uploaded_file(NULL)
    text_input_data(NULL)
  })
  
  # Reativo que combina as variantes do arquivo ou do texto
  all_variants <- reactive({
    if (!is.null(uploaded_file())) {
      v <- readLines(uploaded_file()$datapath)
    } else if (!is.null(text_input_data()) && nzchar(text_input_data()[1])) {
      v <- unlist(strsplit(text_input_data(), "[,;\\n\\r\\t ]+"))
    } else {
      return(NULL)
    }
    # Normalização aqui
    v <- trimws(v)
    v <- v[v != ""]
    v <- toupper(v)
    v
  })
  
  # reativas de dados ---------------------------------------------------------------------------
  df_heatmap <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Dataframes", paste0("dataframe_", input$GENE, "_heatmap.rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo de heatmap não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  
  df_REVEL <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Dataframes", paste0("dataframe_", input$GENE, "_REVEL.rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo de REVEL não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  
  df_structural <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Dataframes", paste0("dataframe_", input$GENE, "_structural.rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo estrutural não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  
  df_rosetta <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Dataframes", paste0("dataframe_", input$GENE, "_rosetta.rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo Rosetta não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  
  dfV <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Dataframes", paste0("dataframe_", input$GENE, "_variants.rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo de variantes não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  
  df.info <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Dataframes", paste0("dataframe_", input$GENE, "_information_.rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo de informações não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  
  fig <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Figures", paste0(input$GENE, "_dominios.png"))
    validate(need(file.exists(caminho),
                  paste("Figura não encontrada para o gene:", input$GENE)))
    caminho
  })
  
  var_user_df <- reactive({
    v <- all_variants()
    req(v)
    
    # Normalização
    v <- trimws(v)       # remove espaços no início/fim
    v <- v[v != ""]      # remove entradas vazias
    v <- toupper(v)      # tudo em maiúsculas (se seu dataset de referência usa assim)
    
    data.frame(V1 = v, stringsAsFactors = FALSE)
  })
  
  # Adicione este renderDataTable ao seu server.R para exibir as variantes na tabela
  output$file1 <- renderDataTable({
    req(all_variants())
    data.frame(Variante = all_variants(), stringsAsFactors = FALSE)
  })
  
  add_variants_to_plot <- function(p, clinvar_selection, gene) {
    caminho <- file.path('./Data/Genes', 'CLINVAR', gene, paste0(clinvar_selection, '.txt'))
    if (!file.exists(caminho)) return(p)
    clinvar_data <- tryCatch(read_clinical(path = caminho, df = df_heatmap()),
                             error = function(e) data.frame())
    if (nrow(clinvar_data) > 0) {
      p <- p + geom_point(
        data = clinvar_data,
        aes(x = Resno, y = ResID),
        fill = get_variant_color(clinvar_selection),
        shape = 22, color = 'gray7', size = 3
      )
    }
    p
  }
  
  get_variant_color <- function(tipo) {
    c(
      'p' = '#CA0000', 'lp' = '#FF9900',
      'vus' = '#48C9B0','lb' = '#FCFD9B',
      'b' = '#77FC02', 'all'= '#F0F8FF'
    )[tipo]
  }
  
  uploaded_variants <- reactive({
    v <- all_variants()
    req(v)           # garante que há variantes antes de prosseguir
    as.character(v)
  })
  
  output$variants_table <- renderDataTable({
    req(var_user_df())
    DT::datatable(var_user_df(),
                  options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  heatmap_plot <- reactive({
    req(df_heatmap(), input$GENE)
    variantes <- uploaded_variants()  # agora funciona para arquivo OU textarea
    
    df_plot <- df_heatmap() %>%
      filter(tolower(trimws(mutation)) %in% tolower(trimws(variantes)))
    
    validate(need(nrow(df_plot) > 0, "Nenhuma das variantes enviadas foi encontrada no heatmap."))
    
    p <- plot_heatmap(df_plot, gene = input$GENE, save = FALSE)
    if (input$clinvar == 'all') {
      p <- reduce(c('p','lp','vus','lb','b'),
                  ~ add_variants_to_plot(.x, .y, input$GENE),
                  .init = p)
    } else {
      p <- add_variants_to_plot(p, input$clinvar, input$GENE)
    }
    p
  })
  
  output$plot <- renderPlot({
    heatmap_plot()
  }, width = function() input$width, res = 96)
  
  output$download_heatmap <- downloadHandler(
    filename = function() paste0("heatmap_", input$GENE, ".png"),
    content = function(file) {
      png(file, width = input$width, height = 600, res = 96)
      print(heatmap_plot())
      dev.off()
    }
  )
  
  in_silico_scores_df <- reactive({
    req(uploaded_variants(), df_heatmap(), df_REVEL())
    
    normalize_mut <- function(x) {
      x <- toupper(trimws(as.character(x)))
      x <- gsub("^P\\.", "", x)
      x <- gsub("[^A-Z0-9]", "", x)
      x
    }
    
    vars <- normalize_mut(uploaded_variants())
    
    # --- AlphaMissense ---
    hm <- df_heatmap()
    hm$mutation_norm <- normalize_mut(hm$mutation)
    am_scores <- hm$score[match(vars, hm$mutation_norm)]
    
    am_pred <- ifelse(
      is.na(am_scores), NA_character_,
      ifelse(am_scores <= 0.0853, "Benign strong",
             ifelse(am_scores <= 0.166,  "Benign moderate",
                    ifelse(am_scores <= 0.316,  "Benign supporting",
                           ifelse(am_scores < 0.787,  "Ambiguous",
                                  ifelse(am_scores < 0.956,  "Pathogenic supporting",
                                         ifelse(am_scores < 0.994,  "Pathogenic moderate",
                                                "Pathogenic strong"))))))
    )
    
    # --- REVEL ---
    rv <- df_REVEL()
    rv$mutation_norm <- normalize_mut(rv$mutation)
    revel_scores <- suppressWarnings(as.numeric(rv$REVEL_score))
    revel_scores <- revel_scores[match(vars, rv$mutation_norm)]
    
    revel_pred <- ifelse(
      is.na(revel_scores), NA_character_,
      ifelse(revel_scores <= 0.016, "Benign strong",
             ifelse(revel_scores <= 0.183, "Benign moderate",
                    ifelse(revel_scores <= 0.29,  "Benign supporting",
                           ifelse(revel_scores <= 0.644, "Ambiguous",
                                  ifelse(revel_scores <= 0.773, "Pathogenic supporting",
                                         ifelse(revel_scores <= 0.932, "Pathogenic moderate",
                                                "Pathogenic strong"))))))
    )
    
    # --- MetaRNN ---
    get_metarnn_from_df <- function(df) {
      if (is.null(df) || !"MetaRNN_score" %in% names(df)) {
        return(rep(NA_real_, length(vars)))
      }
      df$mutation_norm <- normalize_mut(df$mutation)
      vals <- suppressWarnings(as.numeric(df$MetaRNN_score))
      vals[match(vars, df$mutation_norm)]
    }
    
    metarnn_scores <- get_metarnn_from_df(rv)
    
    MetaRNN_pred <- ifelse(
      is.na(metarnn_scores), NA_character_,
      ifelse(metarnn_scores <= 0.108, "Benign strong",
             ifelse(metarnn_scores <= 0.267, "Benign moderate",
                    ifelse(metarnn_scores <= 0.43,  "Benign supporting",
                           ifelse(metarnn_scores >= 0.93,  "Pathogenic strong",
                                  ifelse(metarnn_scores >= 0.841, "Pathogenic moderate",
                                         ifelse(metarnn_scores >= 0.748, "Pathogenic supporting",
                                                "Ambiguous"))))))
    )
    
    data.frame(
      Variant = vars,
      AlphaMissense_Score = am_scores,
      AlphaMissense_Prediction = am_pred,
      REVEL_Score = revel_scores,
      REVEL_Prediction = revel_pred,
      MetaRNN_Score = metarnn_scores,
      MetaRNN_Prediction = MetaRNN_pred,
      stringsAsFactors = FALSE
    )
  })
  
  
  output$in_silico_scores <- DT::renderDataTable({
    DT::datatable(
      in_silico_scores_df(),
      options = list(scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$download_in_silico_scores <- downloadHandler(
    filename = function() {
      paste0("in_silico_predictions_", input$GENE, ".csv")
    },
    content = function(file) {
      write.csv(in_silico_scores_df(), file, row.names = FALSE)
    }
  )
  
  output$download_table_combined <- renderUI({
    req(in_silico_scores_df())
    downloadButton(
      "download_in_silico_scores",
      "Download Full Dataset",
      icon = icon("download"),
      class = "btn-primary"
    )
  })
  
  output$example_input <- downloadHandler(
    
    filename = function() {
      "input_example.txt"
    },
    content = function(file) {
      source_file <- file.path("Data/example_txt", "input_example.txt")
      
      validate(
        need(file.exists(source_file),
             paste("File not found:", source_file))
      )
      file.copy(source_file, file)
    }
  )
  
  output$info_text <- renderUI({
    req(df.info())
    info <- df.info()
    
    # converter SOMENTE a coluna Class
    info$Class <- as.character(info$Class)
    
    # renomear label
    info$Class[info$Class == "Information on gene"] <- "Gene name"
    
    # criar link UniProt
    idx <- which(info$Class == "Uniprot code")
    if (length(idx) == 1) {
      code <- as.character(info$Information[idx])
      
      info$Information[idx] <- list(
        tags$a(
          href = paste0("https://www.uniprot.org/uniprotkb/", code, "/entry"),
          target = "_blank",
          code
        )
      )
    }
    
    idx <- which(info$Class == "RefSeq transcript (MANE Select)")
    
    if (length(idx) == 1) {
      nm <- as.character(info$Information[idx])
      
      row <- mane_info[mane_info$NM == nm, ]
      
      if (nrow(row) == 1) {
        info$Information[idx] <- list(
          tags$a(
            href = row$Link,
            target = "_blank",
            nm
          )
        )
      }
    }
    # =========================
    
    # Linha do download
    info <- rbind(
      info,
      data.frame(
        Class = "Transcript sequence",
        Information = I(list(
          downloadLink(
            outputId = "download_fasta",
            label = "Download here"
          )
        )),
        stringsAsFactors = FALSE
      )
    )
    
    tagList(
      lapply(seq_len(nrow(info)), function(i) {
        div(
          class = "info-line",
          span(class = "info-label", info$Class[i]),
          div(class = "info-value", info$Information[[i]])
        )
      })
    )
  })
  

  output$download_fasta <- downloadHandler(
    
    filename = function() {
      
      info <- df.info()
      
      gene <- trimws(info$Information[
        info$Class == "Information on gene"
      ])
      
      nome <- paste0(gene, ".fasta")
      
      print(nome)
      
      nome
    },
    
    content = function(file) {
      
      req(df.info())
      
      info <- df.info()
      
      gene <- trimws(info$Information[
        info$Class == "Information on gene"
      ])
      
      source_file <- file.path("./Data/Fasta", paste0(gene, ".fasta"))
      
      validate(
        need(file.exists(source_file),
             paste("File not found:", source_file))
      )
      
      file.copy(source_file, file)
    }
  )
  
  ### ADICIONADO POR MAURICIO ####
  output$figure <- renderImage({
    width  <- session$clientData$output_figure_width
    height <- session$clientData$output_figure_height
    list(src = fig(), width = width, height = height,
         alt = "This is alternate text")
  }, deleteFile = F)
  
  
  # TABELA INTERPRETAÇÃO DE VARIANTES CLINICAS
  uploaded_variants_table <- reactive({
    req(var_user_df(), df_heatmap(), df_structural(), df_rosetta(), dfV())
    
    heatmap_mut <- toupper(trimws(as.character(df_heatmap()$mutation)))
    
    # encontra variantes que existem no heatmap
    hits_idx <- which(var_user_df()$V1 %in% heatmap_mut)
    validate(need(length(hits_idx) > 0,
                  "Nenhuma das variantes fornecidas foi encontrada no banco de dados do heatmap."))
    
    # chama as funções com user_df_found em vez do var_user_df() bruto
    a <- get_AM_predictions_user_variants(var_user_df(), df_heatmap())
    b <- get_rosetta_score_user_variants(var_user_df(), df_rosetta())
    c <- assign.known.variants.data(var_user_df(), dfV())
    g <- get_REVEL_predictions_user_variants(var_user_df(), df_REVEL())
    
    vars <- var_user_df()$V1
    idx <- parse_number(vars)
    structural.user <- df_structural()[match(idx, df_structural()$resnumber), ]
    
    # Faz consenso dos preditores in silico
    consenso <- join_and_consensus(var_user_df(), a, g)
    consenso <- consenso[match(c$Variants, consenso$V1), ]
    rownames(consenso) <- NULL

    g <- consenso
    
    if (nrow(b) != length(vars)) {
      b <- correct_data_frame_rosetta(user.variants = var_user_df(), processed.rosetta = b, processed.heatmap = a)
    }
    
    df.full <- data.frame(Variant = a$mutation,
                          Solvent = structural.user$sasa,
                          Domain = structural.user$dominios,
                          conserved = structural.user$conserved,
                          Alphamissense_prediction = g$ACMG,
                          Stability = "none",
                          Interpretation = paste0(structural.user$interpretation),
                          stringsAsFactors = FALSE)
    df.full
    
    # #ajustar informações stability
    df.full$Stability[df.full$Solvent == "exposed"] <- NA
    
    #acertar nomes das colunas
    names(df.full) <- c('Variants','Solvent Accessibility','Domain','Conservation',
                        'Alphamissense','Stability', 'Interpretation')
    
    # ### montar data frame para a aba Variant Interpretation ###
    df.temp <- data.frame(Variant = a$mutation,
                          gnomAD = NA,
                          ACMG_criteria = 'NA',
                          Clinvar_variants = structural.user$clinvar.summary,
                          Interpretation = structural.user$interpretation,
                          Literature = c$Consequence,
                          Reference = c$Reference,
                          Clinvar_Classification = c$Classification, stringsAsFactors = FALSE)
    
    
    df.1 <- format.df.final(df.reduzido = df.temp, df.completo = df.full, df.variates.conhecidas = c)
    
    df.final <- add_ACMG_df_final(df.reduzido = df.1, df.estrutural = structural.user, df.heatmap = g, df.variantes.conhecidas = c)
  })
  
  
  clinvar_values <- c(
    "Pathogenic",
    "Likely Pathogenic",
    #"VUS",
    "Likely Benign",
    "Benign"
  )
  
  clinvar_colors <- c(
    "#CA0000",
    "#FF9900",
    #"#48C9B0",
    "#FCFD9B",
    "#77FC02"
  )
  
  df_resumo <- reactive({
    uploaded_variants_table() %>% 
      dplyr::select(
        Variants, `gnomAD Frequency`, `Clinvar Classification`, 
        `ACMG Criteria`, ACMG_points,
        Interpretation, RENOMIEII_Class, Literature, Reference
      ) %>% 
      dplyr::rename(
        `ACMG points` = ACMG_points,
        `SAVI Classification` = RENOMIEII_Class
      ) %>% 
      dplyr::mutate(
        dplyr::across(
          where(is.character),
          ~ gsub("\\bnone\\b", "", ., ignore.case = TRUE)
        )
      ) %>%
      dplyr::mutate(
        `ACMG Criteria` = gsub(",+", ",", `ACMG Criteria`),
        `ACMG Criteria` = gsub("^,|,$", "", `ACMG Criteria`),
        `ACMG Criteria` = trimws(`ACMG Criteria`)
      ) %>%
      dplyr::mutate(
        has_function_term = grepl(
          "Protein-protein interaction site|oligomerization site|binding site|Phosphorylation site",
          Interpretation,
          ignore.case = TRUE
        ),
        acmg_3 = `ACMG points` == 3,
        high_acmg = `ACMG points` >= 4,
        marker_symbol = dplyr::case_when(
          high_acmg ~ "<span style='color: red; font-size:18px;'>&#9670;</span>",
          acmg_3 & has_function_term ~ "<span style='color: red; font-size:18px;'>&#9670;</span>",
          TRUE ~ ""
        ),
        prioridade = marker_symbol != ""
      ) %>%
      dplyr::arrange(desc(prioridade), desc(`ACMG points`)) %>%
      dplyr::select(-has_function_term, -acmg_3, -high_acmg, -prioridade) %>%
      dplyr::rename(`Priority` = marker_symbol) %>%
      dplyr::relocate(`Priority`, .before = 1) %>%
      dplyr::rename("ACMG score" = "ACMG points")
  })
    
  output$resumo <- DT::renderDataTable({
    datatable(
      df_resumo(),
      escape = FALSE,
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(targets = "_all", className = "dt-center")
        )
      )
    ) %>%
      formatStyle(
        columns = "Priority",
        `text-align` = "center",
        `vertical-align` = "middle"
      )
  })
  
  output$download_table <- renderUI({
    req(input$file1)
    downloadButton("download_resumo", "Download dataframe", icon = icon("download"), class = "btn-primary")
  })
  output$download_resumo <- downloadHandler(
    filename = function() paste0("variant_interpretation_", input$GENE, ".csv"),
    content = function(file) write.csv(df_resumo(), file, row.names = FALSE)
  )
  
  uploaded_variants_table_full <- reactive({
    req(var_user_df(), df_heatmap(), df_structural(), df_rosetta(), dfV(), df_REVEL())
    
    heatmap_mut <- toupper(trimws(as.character(df_heatmap()$mutation)))
    hits_idx <- which(var_user_df()$V1 %in% heatmap_mut)
    
    a <- get_AM_predictions_user_variants(var_user_df(), df_heatmap())
    b <- get_rosetta_score_user_variants(var_user_df(), df_rosetta())
    c <- assign.known.variants.data(var_user_df(), dfV())
    d <- assign.mut.type(var_user_df())
    g <- get_REVEL_predictions_user_variants(var_user_df(), df_REVEL())
    
    vars <- var_user_df()$V1
    idx <- parse_number(vars)
    structural.user <- df_structural()[match(idx, df_structural()$resnumber), ]
    
    if (nrow(b) != length(vars)) {
      b.correct <- correct_data_frame_rosetta(user.variants = var_user_df(), processed.rosetta = b, processed.heatmap = a)
      b <- b.correct; rm(b.correct)
    }
    
    df.full <- data.frame(Variant = a$mutation,
                          Solvent = structural.user$sasa,
                          Domain = structural.user$dominios,
                          conserved = structural.user$conserved,
                          #Alphamissense_prediction = a$AM.prediction,
                          #REVEL_prediction = g$REVEL.prediction,
                          Stability = b$Effect,
                          `Mutation type` = d$type.2,
                          Interpretation = paste0(structural.user$interpretation),
                          stringsAsFactors = FALSE)
    df.full
    
    # #ajustar informações stability
    df.full$Stability[df.full$Solvent == "exposed"] <- NA
    
    #acertar nomes das colunas
    names(df.full) <- c('Variants','Solvent Accessibility','Domain','Conservation',
                        'Stability','Physicochemical features','Interpretation')
    df.full
  })
  
  physico_header <- paste0(
    'Physicochemical features ',
    '<i class="fa fa-info-circle" ',
    'data-toggle="tooltip" ',
    'title="',
    'Classification of amino acid substitutions based on physicochemical properties:&#10;',
    '• conserved: similar physicochemical class&#10;',
    '• gly_cys_pro: involves Gly, Cys or Pro&#10;',
    '• polar_charged: polarity or charge change&#10;',
    '• charge_inversion: charge sign inversion',
    '" ',
    'style="cursor:pointer; color:#1f77b4;"></i>'
  )
  
  output$resumo_full <- DT::renderDataTable({
    
    DT::datatable(
      uploaded_variants_table_full(),
      escape = FALSE,   # importante
      colnames = c(
        "Variants",
        "Solvent Accessibility",
        "Domain",
        "Conservation",
        "Stability",
        physico_header,   # importante
        "Interpretation"
      ),
      options = list(
        scrollX = TRUE,
        initComplete = DT::JS(
          "function(settings, json) {",
          "  $('[data-toggle=\"tooltip\"]').tooltip();",
          "}"
        )
      ),
      rownames = FALSE
    )
  })
  
  output$download_table_full <- renderUI({
    req(input$file1)
    downloadButton("download_resumo_full", "Download Full Dataset", icon = icon("download"), class = "btn-primary")
  })
  output$download_resumo_full <- downloadHandler(
    filename = function() paste0("variant_interpretation_full_", input$GENE, ".csv"),
    content = function(file) write.csv(uploaded_variants_table_full(), file, row.names = FALSE)
  )
  
  # legenda para 3D Viewer
  output$legend3d <- renderUI({
    req(input$file1)
    #h2("Mutated positions are shown in yellow")
  })
  
  # 3D Viewer
  output$structure <- renderNGLVieweR({
    req(input$GENE, all_variants(), df_structural())
    pdb <- file.path('Data/Structures', paste0('estrutura_', input$GENE, '.pdb'))
    validate(need(file.exists(pdb),
                  paste("Arquivo PDB não encontrado para o gene selecionado:", input$GENE)))
    
    # Use all_variants() que agora combina arquivo e texto
    indices <- parse_number(all_variants())
    selecoes <- paste(indices, collapse = " or ")
    
    NGLVieweR(pdb) %>%
      stageParameters(backgroundColor = "white") %>%
      setQuality("low") %>%
      addRepresentation("cartoon",
                        param= list(colorScheme = "uniform",
                                    opacity = 1.0)) %>%
      addRepresentation("spacefill",
                        param = list(colorScheme = "element",
                                     colorValue = "yellow",
                                     sele = paste0("(", selecoes, ") and .CA and :A")))
  })
  
  # 3D Viewer - teste 2
  blaalb <- renderNGLVieweR({
    req(input$GENE, df_structural())
    pdb <- file.path('./Data/Structures', paste0('estrutura_', input$GENE, '.pdb'))
    validate(need(file.exists(pdb),
                  paste("Arquivo PDB não encontrado para o gene selecionado:", input$GENE)))
    
    # Define ClinVar categories and mapping values
    clinvar_mapping <- list(
      "p"   = list(values = c(3,8), color = "#CA0000"),        # Pathogenic
      "lp"  = list(values = c(4,9), color = "#FF9900"),     # Likely Pathogenic
      "vus" = list(values = 5,     color = "#48C9B0"),        # VUS
      "lb"  = list(values = c(2,7), color = "#FCFD9B"), # Likely Benign
      "b"   = list(values = c(1,6), color = "#77FC02")   # Benign
    )
    
    selected_clinvar_category <- input$clinvar
    
    # Base viewer
    viewer <- NGLVieweR(pdb) %>%
      stageParameters(backgroundColor = "white") %>%
      setQuality("low") %>%
      addRepresentation("cartoon",
                        param= list(colorScheme = "uniform",
                                    opacity = 0.5))
    
    # Caso seja "all", adiciona cada categoria separadamente
    if (selected_clinvar_category == "all") {
      for (cat in names(clinvar_mapping)) {
        df_cat <- df_structural() %>%
          filter(class.clinvar %in% clinvar_mapping[[cat]]$values)
        
        if (nrow(df_cat) > 0) {
          selecoes <- paste(df_cat$resnumber, collapse = " or ")
          viewer <- viewer %>%
            addRepresentation("cartoon",
                        param= list(colorScheme = "uniform",
                                    opacity = 0.5)) %>%
            addRepresentation("spacefill",
                              param = list(colorScheme = "uniform",
                                           colorValue = clinvar_mapping[[cat]]$color,
                                           sele = paste0("(", selecoes, ") and .CA and :A")))
        }
      }
      return(viewer)
    }
    
    # Caso não seja "all", filtra apenas a categoria selecionada
    clinvar_values_to_filter <- clinvar_mapping[[selected_clinvar_category]]$values
    validate(
      need(!is.null(clinvar_values_to_filter),
           "Categoria ClinVar selecionada não é válida.")
    )
    
    df_structural_filtered <- df_structural() %>%
      filter(class.clinvar %in% clinvar_values_to_filter)
    
    if (nrow(df_structural_filtered) > 0) {
      selecoes.2 <- paste(df_structural_filtered$resnumber, collapse = " or ")
      viewer <- viewer %>%
        addRepresentation("spacefill",
                          param = list(colorScheme = "uniform",
                                       colorValue = clinvar_mapping[[selected_clinvar_category]]$color,
                                       sele = paste0("(", selecoes.2, ") and .CA and :A")))
    }
    
    return(viewer)
  })
  
  output$structure.2 <- blaalb
  output$structure.3 <- blaalb
  
  # Gráficos Plotly na aba Statistics
  plots_list <- reactive({
    req(input$GENE)
    caminho <- file.path("./Data/Graphics", paste0("graficos_", input$GENE, ".rds"))
    validate(need(file.exists(caminho),
                  paste("Arquivo de gráficos não encontrado para o gene:", input$GENE)))
    readRDS(caminho)
  })
  output$plotb <- renderPlotly({
    p <- plots_list()[["b"]]           # no partial matching
    validate(need(!is.null(p), "Gráfico 'b' não encontrado dentro do RDS."))
    if (inherits(p, "plotly")) {
      p
    } else if (inherits(p, "ggplot")) {
      plotly::ggplotly(p)
    } else {
      validate(need(FALSE, paste0(
        "O item 'b' não é ggplot/plotly. Classe: ",
        paste(class(p), collapse = ", "),
        ". Salve o objeto ggplot (não gtable/grob)."
      )))
    }
  })
  
  output$plotc <- renderPlotly({
    p <- plots_list()[["c"]]
    if (inherits(p, "plotly")) p else if (inherits(p, "ggplot")) plotly::ggplotly(p) else
      validate(need(FALSE, paste0("Item 'c' classe: ", paste(class(p), collapse=", "))))
  })
  
  output$plotd <- renderPlotly({
    p <- plots_list()[["d"]]
    if (inherits(p, "plotly")) p else if (inherits(p, "ggplot")) plotly::ggplotly(p) else
      validate(need(FALSE, paste0("Item 'd' classe: ", paste(class(p), collapse=", "))))
  })
  
  output$plotf <- renderPlotly({
    p <- plots_list()[["f"]]
    if (inherits(p, "plotly")) p else if (inherits(p, "ggplot")) plotly::ggplotly(p) else
      validate(need(FALSE, paste0("Item 'f' classe: ", paste(class(p), collapse=", "))))
  })
  
  output$plotg <- renderPlotly({
    p <- plots_list()[["g"]]
    if (inherits(p, "plotly")) p else if (inherits(p, "ggplot")) plotly::ggplotly(p) else
      validate(need(FALSE, paste0("Item 'g' classe: ", paste(class(p), collapse=", "))))
  })
  
  output$ploth <- renderPlotly({
    p <- plots_list()[["h"]]
    if (inherits(p, "plotly")) p else if (inherits(p, "ggplot")) plotly::ggplotly(p) else
      validate(need(FALSE, paste0("Item 'h' classe: ", paste(class(p), collapse=", "))))
  })
  
}

shinyApp(ui, server)