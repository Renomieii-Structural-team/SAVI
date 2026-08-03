df_heatmap <- reactive({
  req(input$GENE)
  caminho <- file.path('dataframes', paste0('dataframe_', input$GENE, '_heatmap.rds'))
  validate(need(file.exists(caminho),
                paste("Arquivo de heatmap não encontrado para o gene:", input$GENE)))
  readRDS(caminho)
})




df_structural <- reactive({
  req(input$GENE)
  caminho <- file.path('dataframes', paste0('dataframe_', input$GENE, '_structural.rds'))
  validate(need(file.exists(caminho),
                paste("Arquivo estrutural não encontrado para o gene:", input$GENE)))
  readRDS(caminho)
})

df_rosetta <- reactive({
  req(input$GENE)
  caminho <- file.path("dataframes", paste0("dataframe_", input$GENE, "_rosetta.rds"))
  validate(need(file.exists(caminho),
                paste("Arquivo Rosetta não encontrado para o gene:", input$GENE)))
  readRDS(caminho)
})

dfV <- reactive({
  req(input$GENE)
  caminho <- file.path("dataframes", paste0("dataframe_", input$GENE, "_variants.rds"))
  validate(need(file.exists(caminho),
                paste("Arquivo de variantes não encontrado para o gene:", input$GENE)))
  readRDS(caminho)
})

#### ADICIONADO POR MAURICIO #####
df.info <- reactive({
  req(input$GENE)
  caminho <- file.path("dataframes", paste0("dataframe_", input$GENE, "_information.rds"))
  validate(need(file.exists(caminho),
                paste("Arquivo de informações não encontrado para o gene:", input$GENE)))
  readRDS(caminho)
})

fig <- reactive({
  req(input$GENE)
  caminho <- file.path("figuras", paste0(input$GENE, "_dominios.png"))
  validate(need(file.exists(caminho),
                paste("Figura não encontrada para o gene:", input$GENE)))
  #(caminho)
  caminho
})
