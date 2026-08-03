### FUNÇÕES PARA ANÁLISE E INTERPRETAÇÃO ESTRUTURAL DE VARIANTES GENOMICAS ### 

####################################################################################################
# PROCESSAMENTO DE DADOS BRUTOS DO ALPHAMISSENSE
# E GERAÇÃO DE OBJETO R A SER USADO EM OUTRAS FUNÇÕES
####################################################################################################
# argumento: 
# pred = arquivo em formato .tsv obtido a partir das predições do alphamissense para proteína de interesse
make_heatmap <- function(pred){
  input.pred <- read.table(pred, sep='\t')
  input.pred$V2 <- as.character(input.pred$V2)
  to_plot <- data.frame(V2=NA, V4=NA, V6=NA, V8=NA)
  tmp <- data.frame(V2=parse_number(input.pred$V2), V3=substr(input.pred$V2, 1, 1), V4=substr(input.pred$V2, nchar(input.pred$V2), nchar(input.pred$V2)), V6=input.pred$V3)
  tmp$V8 = paste0(tmp$V3,tmp$V2,tmp$V4)
  sequence <- unique(sort(tmp$V2))
  for(i in sequence){
    tmp_pred <- tmp[which(tmp$V2 == i),c(1,3,4,5)]
    native <- data.frame(V2=i, V4=unique(sort(tmp$V3[which(tmp$V2 == i)])),V6=0, V8 = "native")
    tmp_pred <- rbind(tmp_pred, native)
    to_plot <- rbind(to_plot, tmp_pred)
  }
  to_plot$V6[which(is.nan(to_plot$V6))] <- 0
  to_plot <- na.omit(to_plot)
  names(to_plot) <- c("Resno","ResID","score","mutation")
  to_plot$ResID <- factor(to_plot$ResID,levels = sort(unique(to_plot$ResID),
                                                      decreasing = T))
  out <- to_plot  
}
####################################################################################################
# função para extrair as predições do Alphamissense # para conjunto 
# de variantes de interesse submetidas pelo usuário
####################################################################################################
# argumentos da função:
# user.var = arquivo texto com variantes em formato de 1 letra (1 por linha)
# dataframe_AM = objeto previamente criado a partir da função make_heatmap
get_AM_predictions_user_variants <- function(user.var,dataframe_AM){
  variants <- user.var
  AM.var <- NULL
  for (i in 1:length(variants$V1)) {
    aux = dataframe_AM[grep(variants$V1[i],dataframe_AM$mutation),]
    AM.var = rbind(AM.var,aux)
  }
  for (i in 1:length(AM.var$mutation)) {
    #AM.var$AM.prediction[i] = 'none'
    #if(AM.var$score[i] <= 0.0853){AM.var$AM.prediction[i] <- "BP4 strong"}
    #if(AM.var$score[i]  > 0.0853 & AM.var$score[i] <= 0.166){AM.var$AM.prediction[i] <- "BP4 moderate"}
    #if(AM.var$score[i]  > 0.166 &AM.var$score[i] <= 0.316){AM.var$AM.prediction[i] <- "BP4 supporting"}
    #if(AM.var$score[i]  > 0.316 &AM.var$score[i] <= 0.787){AM.var$AM.prediction[i] <- 'none'}
    #if(AM.var$score[i] >= 0.787 &AM.var$score[i] < 0.956){AM.var$AM.prediction[i] <- "PP3 supporting"}
    #if(AM.var$score[i] >= 0.956 &AM.var$score[i] < 0.994){AM.var$AM.prediction[i] <- "PP3 moderate"}
    #if(AM.var$score[i] >= 0.994){AM.var$ACMG[i] <- "PP3 strong"}
    
    #AM.var$ACMG[i] = 'VUS'
    AM.var$ACMG[i] = ' '
    if(AM.var$score[i] <= 0.0853){AM.var$ACMG[i] <- "BP4 strong"}
    if(AM.var$score[i]  > 0.0853 & AM.var$score[i] <= 0.166){AM.var$ACMG[i] <- "BP4 moderate"}
    if(AM.var$score[i]  > 0.166 &AM.var$score[i] <= 0.316){AM.var$ACMG[i] <- "BP4 supporting"}
    if(AM.var$score[i]  > 0.316 &AM.var$score[i] < 0.787){AM.var$ACMG[i] <- 'none'}
    if(AM.var$score[i] >= 0.787 &AM.var$score[i] < 0.956){AM.var$ACMG[i] <- "PP3 supporting"}
    if(AM.var$score[i] >= 0.956 &AM.var$score[i] < 0.994){AM.var$ACMG[i] <- "PP3 moderate"}
    if(AM.var$score[i] >= 0.994){AM.var$ACMG[i] <- "PP3 strong"}
    
    # pontuação para critério ACMG 
    AM.var$ACMG.points[i] = 0
    if(AM.var$score[i] <= 0.0853){AM.var$ACMG.points[i] <- -4}
    if(AM.var$score[i]  > 0.0853 & AM.var$score[i] <= 0.166){AM.var$ACMG.points[i] <- -2}
    if(AM.var$score[i]  > 0.166 &AM.var$score[i] <= 0.316){AM.var$ACMG.points[i] <- -1}
    if(AM.var$score[i] >= 0.787 &AM.var$score[i] < 0.956){AM.var$ACMG.points[i] <- 1}
    if(AM.var$score[i] >= 0.956 &AM.var$score[i] < 0.994){AM.var$ACMG.points[i] <- 2}
    if(AM.var$score[i] >= 0.994){AM.var$ACMG.points[i] <- 4}
    
    
  }
  return(AM.var)
}

####################################################################################################
# função para extrair as predições do REVEL # para conjunto  -- PAULA MEXEU AQUI
# de variantes de interesse submetidas pelo usuário
####################################################################################################

get_REVEL_predictions_user_variants <- function(user.var, dataframe) {
  # monta REVEL.var a partir das variantes do usuário
  REVEL.var <- NULL
  for (i in seq_len(nrow(user.var))) {
    aux <- dataframe[grep(user.var$V1[i], dataframe$mutation, fixed = TRUE), , drop = FALSE]
    REVEL.var <- rbind(REVEL.var, aux)
  }
  REVEL.var <- as.data.frame(REVEL.var, stringsAsFactors = FALSE)
  
  if (nrow(REVEL.var) == 0) return(REVEL.var)
  
  # --- crie as colunas ANTES de usar ---
  REVEL.var$ACMG <- NA_character_
  REVEL.var$ACMG.points <- NA_integer_
  REVEL.var$MetaRNN.ACMG <- "none"
  REVEL.var$MetaRNN.ACMG.points <- 0L
  
  # --- REVEL: regras ACMG ---
  for (i in seq_len(nrow(REVEL.var))) {
    sc <- REVEL.var$REVEL_score[i]
    if (is.na(sc)) next
    
    # classe
    if (sc <= 0.016) {
      REVEL.var$ACMG[i] <- "BP4 strong"
    } else if (sc <= 0.183) {
      REVEL.var$ACMG[i] <- "BP4 moderate"
    } else if (sc <= 0.29) {
      REVEL.var$ACMG[i] <- "BP4 supporting"
    } else if (sc <= 0.644) {
      REVEL.var$ACMG[i] <- "none"
    } else if (sc <= 0.773) {
      REVEL.var$ACMG[i] <- "PP3 supporting"
    } else if (sc <= 0.932) {
      REVEL.var$ACMG[i] <- "PP3 moderate"
    } else { sc > 0.932
      REVEL.var$ACMG[i] <- "PP3 strong"
    }
    
    # pontos
    if (sc <= 0.016) {
      REVEL.var$ACMG.points[i] <- -4L
    } else if (sc <= 0.183) {
      REVEL.var$ACMG.points[i] <- -2L
    } else if (sc <= 0.29) {
      REVEL.var$ACMG.points[i] <- -1L
    } else if (sc <= 0.644) {
      REVEL.var$ACMG.points[i] <- 0L
    } else if (sc <= 0.773) {
      REVEL.var$ACMG.points[i] <- 1L
    } else if (sc <= 0.932) {
      REVEL.var$ACMG.points[i] <- 2L
    } else { # sc > 0.932
      REVEL.var$ACMG.points[i] <- 4L
    }
  }
  
  # --- MetaRNN: atenção à ordem das faixas ---
  for (i in seq_len(nrow(REVEL.var))) {
    m <- REVEL.var$MetaRNN_score[i]
    if (is.na(m)) next
    
    if (m >= 0.939) {
      REVEL.var$MetaRNN.ACMG[i] <- "PP3 strong";     REVEL.var$MetaRNN.ACMG.points[i] <- 4L
    } else if (m >= 0.841) {
      REVEL.var$MetaRNN.ACMG[i] <- "PP3 moderate";   REVEL.var$MetaRNN.ACMG.points[i] <- 2L
    } else if (m >= 0.748) {
      REVEL.var$MetaRNN.ACMG[i] <- "PP3 supporting"; REVEL.var$MetaRNN.ACMG.points[i] <- 1L
    } else if (m <= 0.108) {
      REVEL.var$MetaRNN.ACMG[i] <- "BP4 strong";     REVEL.var$MetaRNN.ACMG.points[i] <- -4L
    } else if (m <= 0.267) {
      REVEL.var$MetaRNN.ACMG[i] <- "BP4 moderate";   REVEL.var$MetaRNN.ACMG.points[i] <- -2L
    } else if (m <= 0.43) {
      REVEL.var$MetaRNN.ACMG[i] <- "BP4 supporting"; REVEL.var$MetaRNN.ACMG.points[i] <- -1L
    } else {
      REVEL.var$MetaRNN.ACMG[i] <- "none";           REVEL.var$MetaRNN.ACMG.points[i] <- 0L
    }
  }
  
  return(REVEL.var)
}

####################################################################################################
# função para extrair os resultados do Rosetta
# para conjunto de variantes de interesse submetidas pelo usuário
####################################################################################################
# argumentos da função:
# user.var = arquivo texto com variantes em formato de 1 letra (1 por linha)
# dataframe_rosetta = objeto previamente criado a partir da script dataframe_rosetta.R
get_rosetta_score_user_variants <- function(user.var,dataframe_rosetta){
  dataframe.tapa.buraco = data.frame(Variant = NA, ddG = NA, Effect = NA)
  variants <- user.var
  ROSETTA.var <- NULL
  for (i in 1:length(variants$V1)) {
    aux = dataframe_rosetta[grep(variants$V1[i], dataframe_rosetta$Variant),]
    # if (nrow(aux) == 0) {aux = dataframe.tapa.buraco} # para inserir NA onde não há dados do Rosetta
    ROSETTA.var = rbind(ROSETTA.var,aux)
  }
  return(ROSETTA.var)
}
####################################################################################################
# função para corrigir o data frame gerado a partir dos dados do rosetta
# para o caso de os cálculos estruturais terem sido realizados em proteínas 
# sem todos os resíduos
correct_data_frame_rosetta <- function(user.variants,processed.rosetta,processed.heatmap){
  bla <- user.variants
  bla$V1 <- ifelse(bla$V1 %in% processed.rosetta$Variant, processed.rosetta$Variant, 'none')
  df_none <- as.data.frame(bla[bla$V1 == 'none', ])
  df_none$V2 <- NA
  df_none$V3 <- NA
  names(df_none) <- c( 'Variant', 'ddG', 'Effect')
  teste <- rbind(df_none,processed.rosetta)
  bla <- user.variants
  bla$V1 <- ifelse(bla$V1 %in% processed.rosetta$Variant, 'none', bla$V1) 
  teste$Variant[teste$Variant == 'none'] <- bla$V1[bla$V1 != 'none']
  teste <- teste[match(processed.heatmap$mutation, teste$Variant), ]
  out <- teste
  return(out)
}
####################################################################################################
# função para processar o data frame reduzido (df.final)
format.df.final <- function(df.reduzido,df.completo,df.variates.conhecidas){
  ##### 1 - adiciona coluna stability #####
  # e corrige os NAs da estabilidade para none
  df.reduzido$stability <- df.completo$Stability
  df.reduzido$stability[is.na(df.reduzido$stability)] <- "none"
  
  ##### 2- constrói nova coluna "stability + interpretation" #####
  df.reduzido$teste <- ifelse(
    # Condição: A estabilidade NÃO é 'none'?
    df.reduzido$stability != 'none',
    # Valor se a condição for VERDADEIRA:
    # Concatena os valores
    paste0(df.reduzido$stability, " and ", df.reduzido$Interpretation),
    # Valor se a condição for FALSA:
    # Apenas retorna o valor original da interpretação, sem mudar nada
    df.reduzido$Interpretation
  )
  ##### 3 remove colunas indesejadas ######
  df.reduzido$stability <- NULL
  df.reduzido$Interpretation <- NULL
  
  ##### 4 obtém indices e efeitos conhecidos  ######
  # a partir do dataframe de variantes conhecidas (dfV)
  inds.efeitos.conhecidos <- which(!is.na(df.variates.conhecidas$Consequence))
  efeitos.conhecidos <- df.variates.conhecidas$Consequence[inds.efeitos.conhecidos]
  
  ### 5 loop para adicionar  os dados previamente #####
  # conhecidas (dataframe dfV.user) na coluna literatura e nomad
  n = 0
  for (i in inds.efeitos.conhecidos ) {
    n = n+1
    df.reduzido$Literature[i] <- efeitos.conhecidos[n] 
  }
  
  #### 6 add coluna gnomad #####
  df.reduzido$gnomAD <- df.variates.conhecidas$Frequency
  df.reduzido$gnomAD[is.na(df.reduzido$gnomAD)] <- 'NA'
  ### 7 Acertar nomes das colunas ####
  names(df.reduzido) <- c('Variants', 'gnomAD Frequency', 'ACMG Criteria', 'Clinvar Variants',
                          'Literature', 'Reference','Clinvar Classification','Interpretation')
  
  ### 8 retorna o dataframe modificado ###
  out <- df.reduzido
  return(out)
}

###################################################################################################
# Adequar critério de predição in silico quando houver PM1 ou PM5 naquela posição
# P3 strong/moderate -> supporting / # BP4 strong/moderate -> supporting 
adjust_acmg_rules_in_silico_ <- function(df.heatmap,df.estrutural){
  aux <- df.heatmap
  for (i in 1:length(df.heatmap$ACMG)) {
    # condição para variantes patogenicas com criterios PP3 moderate e strong
    if(df.heatmap$REVEL_score[i] >= 0.773){
      if (df.estrutural$ACMG_pm1[i] != 'none' | df.estrutural$ACMG_PM5[i] != 'none') {
        aux$ACMG[i] <- 'PP3 supporting'
        aux$ACMG.points[i] <- 1}
    }
    # condição para variantes benignas com criterios BP4 moderate e strong
    if(df.heatmap$REVEL_score[i] <= 0.183){
      if (df.estrutural$ACMG_pm1[i] != 'none' | df.estrutural$ACMG_PM5[i] != 'none') {
        aux$ACMG[i] <- 'BP4 supporting'
        aux$ACMG.points[i] <- -1}
    }
  }
  return(aux)
}

###################################################################################################
add_ACMG_df_final <- function(df.reduzido, df.estrutural, df.heatmap, df.variantes.conhecidas) {
  
  # --- 2. Garantir colunas numéricas para soma ---
  safe_num <- function(x) as.numeric(as.character(x))
  
  # substitui NAs por 0 e converte fatores em números
  df.variantes.conhecidas$ACMG.points           <- safe_num(df.variantes.conhecidas$ACMG.points)
  df.variantes.conhecidas$ACMG.clinical.points  <- safe_num(df.variantes.conhecidas$ACMG.clinical.points)
  df.heatmap$ACMG.points                        <- safe_num(df.heatmap$ACMG.points)
  df.estrutural$ACMG.points.pm1                 <- safe_num(df.estrutural$ACMG.points.pm1)
  df.estrutural$ACMG.points.pm5                 <- safe_num(df.estrutural$ACMG.points.pm5)
  df.variantes.conhecidas$ACMG.PP2.points       <- safe_num(df.variantes.conhecidas$ACMG.PP2.points)
  df.variantes.conhecidas$ACMG.BP1.points       <- safe_num(df.variantes.conhecidas$ACMG.BP1.points)
  
  df.reduzido$`ACMG Criteria` <- paste0(
    df.estrutural$ACMG_pm1, ",",
    df.variantes.conhecidas$ACMG_freq, ",",
    df.variantes.conhecidas$ACMG.clinical, ",",
    df.heatmap$ACMG, ",",
    df.estrutural$ACMG_PM5,",",
    df.variantes.conhecidas$ACMG.BP1,",",
    ifelse(df.estrutural$ACMG.points.pm1 == 0,
           df.variantes.conhecidas$ACMG.PP2,
           "none")
  )
  
  df.variantes.conhecidas$ACMG.points[is.na(df.variantes.conhecidas$ACMG.points)] <- 0
  df.variantes.conhecidas$ACMG.clinical.points[is.na(df.variantes.conhecidas$ACMG.clinical.points)] <- 0
  df.heatmap$ACMG.points[is.na(df.heatmap$ACMG.points)] <- 0
  df.estrutural$ACMG.points.pm1[is.na(df.estrutural$ACMG.points.pm1)] <- 0
  df.estrutural$ACMG.points.pm5[is.na(df.estrutural$ACMG.points.pm5)] <- 0
  df.variantes.conhecidas$ACMG.PP2.points[is.na(df.variantes.conhecidas$ACMG.PP2.points)] <- 0
  df.variantes.conhecidas$ACMG.BP1.points[is.na(df.variantes.conhecidas$ACMG.BP1.points)] <- 0
  
  # --- 3. Somar pontos ACMG ---
  df.reduzido$ACMG_points <- 
    df.variantes.conhecidas$ACMG.points +
    df.variantes.conhecidas$ACMG.clinical.points +
    df.heatmap$ACMG.points +
    df.estrutural$ACMG.points.pm1 +
    df.estrutural$ACMG.points.pm5 +
    df.variantes.conhecidas$ACMG.BP1.points +
    ifelse(df.estrutural$ACMG.points.pm1 == 0,
           df.variantes.conhecidas$ACMG.PP2.points,
           0)
  
  # --- 4. Classificação ACMG ---
  df.reduzido$RENOMIEII_Class <- dplyr::case_when(
    df.reduzido$ACMG_points >= 10 ~ 'Pathogenic',
    df.reduzido$ACMG_points >= 6  ~ 'Likely Pathogenic',
    df.reduzido$ACMG_points >= 4  ~ 'VUS Leaning Pathogenic',
    df.reduzido$ACMG_points %in% c(0, 1) ~ "VUS Leaning Benign",
    df.reduzido$ACMG_points <= -7  ~ "Benign",
    df.reduzido$ACMG_points <= -1  ~ "Likely Benign",
    TRUE ~ 'VUS')
  
  return(df.reduzido)
}

####################################################################################################
# função para coletar dados de variantes já caracterizadas na literatura
# ou variantes com alta frequencia em população saudável (gnomAD)
# para conjunto de variantes de interesse submetidas pelo usuário
# caso não haja informações para as variantes submetidas é introduzido NA
##############################################################################################
# argumentos da função:
# user.var = arquivo texto com variantes em formato de 1 letra (1 por linha)
# dataframe_= dataframe com dados de variantes conhecidas e frequencias em população saudável
#assign.known.variants.data <- function(user.var, dataframe){
#  conhecidas <- user.var
#  df.CONHECIDAS <- NULL
#  for (i in 1:length(conhecidas$V1)){
#    aux =dataframe[grep(conhecidas$V1[i], dataframe$Variant),]
#    df.CONHECIDAS = rbind(df.CONHECIDAS, aux)
#  }
#  return(df.CONHECIDAS)
#}

assign.known.variants.data <- function(user.var, dataframe) {
  
  ### Fazendo analise de PP2 e BP1
  
  a <- as.data.frame(table(dataframe$Classification))
  
  den <- sum(a$Freq[a$Var1 != "VUS"])
  
  if (den > 0) {
    
    PP2_value <- round(
      sum(a$Freq[grep("P", a$Var1)]) / den, 3)
    
    BP1_value <- round(
      sum(a$Freq[grep("B", a$Var1)]) / den, 3)
    
  } else {
    
    PP2_value <- 0
    BP1_value <- 0
  }
  ######
  
  conhecidas <- user.var

  # helper: assign a value handling factor/character safely
  set_val <- function(x, idx, val) {
    if (is.factor(x)) {
      lvls <- levels(x)
      if (!val %in% lvls) levels(x) <- c(lvls, val)
      x[idx] <- val
    } else {
      x[idx] <- val
    }
    x
  }
  
  # start with empty df, same columns/types
  df.CONHECIDAS <- dataframe[0, , drop = FALSE]
  
  for (i in seq_len(nrow(conhecidas))) {
    hit <- grepl(conhecidas$V1[i], dataframe$Variant)
    
    if (!any(hit)) {
      # one NA row with same types + mark as unmatched
      aux <- dataframe[1, , drop = FALSE]
      aux[] <- NA
      aux$Variants <- conhecidas$V1[i]
      aux$.__unmatched <- TRUE
    } else {
      aux <- dataframe[hit, , drop = FALSE]
      aux$.__unmatched <- FALSE
    }
    
    df.CONHECIDAS <- rbind(df.CONHECIDAS, aux)
  }
  
  if (PP2_value >= 0.808){
  
  df.CONHECIDAS$ACMG.PP2 <-  "PP2"
  df.CONHECIDAS$ACMG.PP2.points <- 1
    
  } else {
    df.CONHECIDAS$ACMG.PP2 <-  "none"
    df.CONHECIDAS$ACMG.PP2.points <- 0
  }
  
  gene_mechanism <- unique(dataframe$mechanism)
  
  if (gene_mechanism == "LOF"){
    
  if (BP1_value >= 0.569){
    
    df.CONHECIDAS$ACMG.BP1 <-  "BP1"
    df.CONHECIDAS$ACMG.BP1.points <- -1
    
  } else {
    df.CONHECIDAS$ACMG.BP1 <-  "none"
    df.CONHECIDAS$ACMG.BP1.points <- 0
  }
  } else{
    df.CONHECIDAS$ACMG.BP1 <-  "none"
    df.CONHECIDAS$ACMG.BP1.points <- 0
  }
  # fill defaults ONLY on unmatched rows
  um <- df.CONHECIDAS$.__unmatched %in% TRUE
  
  if ("ACMG_freq"   %in% names(df.CONHECIDAS)) 
    df.CONHECIDAS$ACMG_freq   <- set_val(df.CONHECIDAS$ACMG_freq,   um, "PM2")
  
  if ("ACMG.points" %in% names(df.CONHECIDAS)) 
    df.CONHECIDAS$ACMG.points <- set_val(df.CONHECIDAS$ACMG.points, um, 1)
  
  if ("Classification" %in% names(df.CONHECIDAS)) 
    df.CONHECIDAS$Classification <- set_val(df.CONHECIDAS$Classification, um, 'NA')
  
  
  # optional: also set clinical defaults for unmatched
  # if ("ACMG.clinical" %in% names(df.CONHECIDAS))
  #   df.CONHECIDAS$ACMG.clinical <- set_val(df.CONHECIDAS$ACMG.clinical, um, "none")
  # if ("ACMG.clinical.points" %in% names(df.CONHECIDAS))
  #   df.CONHECIDAS$ACMG.clinical.points <- set_val(df.CONHECIDAS$ACMG.clinical.points, um, 0)
  
  # drop helper column
  df.CONHECIDAS$.__unmatched <- NULL
  
  return(df.CONHECIDAS)
}


######################################################################################################
# função para plotar o heatmap 
######################################################################################################
plot_heatmap <- function(df, gene, save = TRUE){
  textcol <- "black" # text color
  size.prot = max(df$Resno)
  
  # sub-rotina pra calcular a separação em x
  calc.x.axis.sep <- function(size){
    if(size < 600){saida = 25}
    if(size > 600 & size < 1200){saida = 50}
    if(size > 1200){saida = 100}
    out <- saida
  }
  sep.scale <- calc.x.axis.sep(size.prot)
  aux <- round(size.prot/sep.scale)
  lab.x <-c(1,seq(sep.scale,aux*sep.scale-1, sep.scale))
  
  # sub-rotina pra calcular o tamanho da figura baseado
  # no n~umero de resíduos da proteína
  calc.plot.size <- function(size){
    if(size < 600){saida = c(6,8)}
    if(size > 600 & size < 1200){saida = c(6,12)}
    if(size > 1200){saida = c(8,12)}
    out <- saida
  }
  size.plot <- calc.plot.size(size.prot)
  
  p <- ggplot(df, aes(Resno, ResID, fill= score)) +
    geom_tile(alpha = 0.5, color = 'gray7', linewidth = 0.1) +
    coord_fixed(ratio = 7) +  # Set the aspect ratio to 1 for squares
    xlab("Amino acid index") +
    ylab("Mutated Residue") +
    theme_grey(base_size=10)+
    theme(legend.position="right", legend.direction="vertical",
          legend.title=element_text(colour=textcol),
          legend.margin=margin(grid::unit(0, "cm")),
          legend.text=element_text(colour=textcol, size=7, face="bold"),
          legend.key.height=grid::unit(0.8, "cm"),
          legend.key.width=grid::unit(0.2, "cm"),
          axis.text.x=element_text(size=10, colour=textcol),
          axis.text.y=element_text(vjust=0.2, colour=textcol),
          axis.ticks=element_line(size=0.4),
          plot.background=element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border=element_blank(),
          plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm"),
          plot.title=element_text(colour=textcol, hjust=0, size=14, face="bold"))+
    scale_x_continuous(expand = c(0, 0),breaks = lab.x, labels = lab.x) +
    scale_y_discrete(expand = c(0, 0))+
    colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", mid = 0.5,l1 = 30, l2 = 100, p1 = .9, p2 = 1.2)
  
  #export object
  if(save){
    saveRDS(p, file = paste0('objetos_ggplot/heatmap_',gene,'.rds'))
    ggsave(p, filename=output, height=size.plot[1], width=size.plot[2], units="in", dpi=200)
  }
  p
}
###################################################################################################
# função para ler dados de mutações misense e converter
# para um data frame no formato do ggplot
read_clinical <- function(path,df){
  data <- read.table(path)
  data$V1 <- as.character(data$V1)
  saida <- data.frame(Resno=parse_number(data$V1), ResID=substr(data$V1, nchar(data$V1), nchar(data$V1)), 
                      score = rep(x = 0, length(data$V1)), mutation=rep(NA, length(data$V1)))
  saida$ResID <- as.character(saida$ResID)
  for(i in 1:length(saida$Resno)){
    saida$score[i] <- df$score[df$ResID == saida$ResID[i] & df$Resno == saida$Resno[i]]
    saida$mutation[i] <- df$mutation[df$ResID == saida$ResID[i] & df$Resno == saida$Resno[i]]
  }
  out <- saida
}
#######################################################################################################
# FUNÇÃO PARA SE DETERMINAR OS TIPOS DE AMINOÁCIDOS PRESENTES
# EM VARIANTES GENÔMICAS E TAMBÉM DA NATUREZA DA MUTAÇÃO
# INPUT DEVE SER UMA TABELA COM A VARIANTE NA PRIMEIRA COLUNA
# E A CLASSIFICAÇÃO DA VARIANTE NO CLINVAR NA SEGUNDA
assign.mut.type <- function(data){
  var <- data$V1
  #clinvar <- data$V2
  # separar tipos de aminoácidos
  hydrophobic = c('A', 'V', 'I', 'L', 'M', 'F', 'W')
  gly_cys_pro = c('G','C','P')
  polar = c('S', 'T', 'Y', 'N', 'Q')
  negative = c('D', 'E')
  positive = c('R', 'H', 'K')
  
  # separar os resíduos WT e mutações em objetos separados
  wt <- substr(var, 1, 1)
  resno <- parse_number(var)
  mut <- substr(var, nchar(var), nchar(var))
  
  # Criar data frame contendo NAs para ser preenchido 
  # com os tipos de resíduos selvage. mutado e a natureza da mutação
  #aux <- data.frame(matrix(NA, nrow = length(var), ncol=5))
  aux <- data.frame(matrix(NA, nrow = length(var), ncol=4))
  #names(aux) <- c('variant','wt.res','mut.res','type','clinvar')
  names(aux) <- c('variant','wt.res','mut.res','type')
  
  # loop para preenchimento do data.frame
  for(i in 1:length(var)){
    # assign type of WT residue 
    if(wt[i] %in% hydrophobic){aux$wt.res[i] <- 'hydrophobic'}
    if(wt[i] %in% gly_cys_pro){aux$wt.res[i] <- 'gly_cys_pro'}
    if(wt[i] %in% polar){aux$wt.res[i] <- 'polar'}
    if(wt[i] %in% positive){aux$wt.res[i] <- 'positive'}
    if(wt[i] %in% negative){aux$wt.res[i] <- 'negative'}
    # assign type of MUTATED residue 
    if(mut[i] %in% hydrophobic){aux$mut.res[i] <- 'hydrophobic'}
    if(mut[i] %in% gly_cys_pro){aux$mut.res[i] <- 'gly_cys_pro'}
    if(mut[i] %in% polar){aux$mut.res[i] <- 'polar'}
    if(mut[i] %in% positive){aux$mut.res[i] <- 'positive'}
    if(mut[i] %in% negative){aux$mut.res[i] <- 'negative'}
    # assign type of substitution residue
    aux$type[i] <- paste0(aux$wt.res[i],"_",aux$mut.res[i])
    if(aux$wt.res[i] == aux$mut.res[i]){aux$type[i] <- 'conserved'}
    if(any(aux[i,] == 'gly_cys_pro', na.rm = T)){aux$type[i] <- 'gly_cys_pro'} # if any residue is gly cys pro - mutation type = 'gly cys pro'
    # assign clinvar status on the data frame
    #aux$clinvar[i] <- clinvar[i]
    # assign varaint on the data frame
    aux$variant[i] <- var[i]
  }
  # remover redundancias para plot exclusivamente - coluna type 2
  aux$type.2 <- aux$type
  aux$type.2[(aux$type == 'polar_hydrophobic'| aux$type == 'hydrophobic_polar')] <-'hydrophobic_polar'
  aux$type.2[(aux$type == 'positive_negative'| aux$type == 'negative_positive')] <-'charge_inversion'
  aux$type.2[(aux$type == 'polar_negative'| aux$type== 'polar_positive')] <-'polar_charged'
  aux$type.2[(aux$type == 'negative_polar'| aux$type == 'positive_polar')] <-'polar_charged'
  aux$type.2[(aux$type == 'hydrophobic_negative'| aux$type == 'hydrophobic_positive')] <-'hydrophobic_charged'
  aux$type.2[(aux$type == 'positive_hydrophobic'| aux$type == 'negative_hydrophobic')] <-'hydrophobic_charged'
  #sort data frame by clinvar status and then by mutation type
  #aux <- aux[order(aux$clinvar, aux$type.2),]
  #aux <- aux[order(aux$type.2),]
  out <- aux
}
#################################################################################################################


## TESTE PAULA

map_label_to_strength_dir_points <- function(lbl) {
  
  if (is.na(lbl) || lbl == "" || tolower(lbl) == "none") {
    return(list(strength="none", dir="none", points=0L))
  }
  
  L <- tolower(lbl)
  
  # ---- direção (procura em qualquer posição)
  if (grepl("\\bbp4\\b", L)) {
    dir <- "benign"
  } else if (grepl("\\bpp3\\b", L)) {
    dir <- "pathogenic"
  } else {
    dir <- "none"
  }
  
  # ---- força
  if (grepl("strong", L)) {
    strength <- "strong"
  } else if (grepl("moderate", L)) {
    strength <- "moderate"
  } else if (grepl("support", L)) {
    strength <- "supporting"
  } else {
    strength <- "none"
  }
  
  # ---- pontos
  mag <- switch(strength,
                "supporting" = 1L,
                "moderate"   = 2L,
                "strong"     = 4L,
                0L)
  
  points <- if (dir == "benign") -mag
  else if (dir == "pathogenic") mag
  else 0L
  
  list(strength=strength, dir=dir, points=points)
}



.strength_rank <- function(x) {
  match(tolower(x),
        c("supporting", "moderate", "strong"),
        nomatch = 0L)
}

.cap_strength <- function(strength, cap_level) {
  
  if (is.infinite(cap_level)) return(strength)
  
  levels <- c("supporting","moderate","strong")
  s_idx  <- .strength_rank(strength)
  c_idx  <- .strength_rank(cap_level)
  
  levels[min(s_idx, c_idx)]
}

.points_from_dir_strength <- function(dir, strength) {
  
  mag <- switch(tolower(strength),
                "supporting" = 1L,
                "moderate"   = 2L,
                "strong"     = 4L,
                0L)
  
  if (tolower(dir) == "benign") -mag
  else if (tolower(dir) == "pathogenic") mag
  else 0L
}

.make_acmg_label <- function(dir, strength) {
  
  code <- if (dir == "benign") "BP4"
  else if (dir == "pathogenic") "PP3"
  else "none"
  
  if (code == "none" || strength == "none") {
    "none"
  } else {
    paste(code, strength)
  }
}

# main function
consensus_acmg <- function(labels, sources = c("AM","REVEL","MetaRNN")) {
  
  info      <- lapply(labels, map_label_to_strength_dir_points)
  strengths <- vapply(info, `[[`, "", "strength")
  dirs      <- vapply(info, `[[`, "", "dir")
  
  # localizar REVEL
  revel_idx <- which(sources == "REVEL")
  if (length(revel_idx) != 1) {
    warning("Fonte REVEL não encontrada; retornando none.")
    return(list(label = "none", points = 0L))
  }
  
  # se REVEL for none → retorna none
  if (dirs[revel_idx] == "none" || strengths[revel_idx] == "none") {
    return(list(label = "none", points = 0L))
  }
  
  # --------------------------------------------------
  # CONSENSO BASEADO EM (dir + strength)
  # --------------------------------------------------
  
  pairs <- paste(dirs, strengths)
  valid <- pairs != "none none"
  
  if (!any(valid)) {
    return(list(label="none", points=0L))
  }
  
  max_freq <- max(table(pairs[valid]))
  
  cap_level <- if (max_freq == 3) {
    Inf
  } else if (max_freq == 2) {
    "moderate"
  } else {
    "supporting"
  }
  
  revel_dir      <- dirs[revel_idx]
  revel_strength <- strengths[revel_idx]
  
  capped_strength <- .cap_strength(revel_strength, cap_level)
  
  out_points <- .points_from_dir_strength(revel_dir, capped_strength)
  out_label  <- .make_acmg_label(revel_dir, capped_strength)
  
  list(label = out_label, points = out_points)
}


# ------------------------------------------------------------
# 4) Junta AM + REVEL (com MetaRNN dentro) e computa consenso
#    - user.var: deve ter 'V1' (ou 'mutation'); usaremos para preservar a ordem/linhas do usuário
#    - AM.var   = output de get_AM_predictions_user_variants()
#    - REVEL.var= output de get_REVEL_predictions_user_variants() (já contém MetaRNN.ACMG*)
join_and_consensus <- function(user.var, AM.var, REVEL.var) {
  # detectar coluna do usuário
  mut_user <- if ("mutation" %in% names(user.var)) "mutation" else if ("V1" %in% names(user.var)) "V1" else NULL
  if (is.null(mut_user)) stop("user.var precisa ter 'mutation' ou 'V1'.")
  
  # normalizar colunas de mutação
  if (!("mutation" %in% names(AM.var))) stop("AM.var precisa ter coluna 'mutation'.")
  if (!("mutation" %in% names(REVEL.var))) stop("REVEL.var precisa ter coluna 'mutation'.")
  
  # selecionar o essencial de cada tabela para evitar colisões de nome
  keep_am    <- intersect(c("mutation","score","AM.prediction","ACMG","ACMG.points"), names(AM.var))
  keep_revel <- intersect(c("mutation","REVEL_score","REVEL.prediction","ACMG","ACMG.points",
                            "MetaRNN_score","MetaRNN.ACMG","MetaRNN.ACMG.points"), names(REVEL.var))
  
  AM.slim    <- AM.var[, keep_am, drop=FALSE]
  names(AM.slim) <- sub("^ACMG$", "AM.ACMG", names(AM.slim))
  names(AM.slim) <- sub("^ACMG\\.points$", "AM.ACMG.points", names(AM.slim))
  
  REVEL.slim <- REVEL.var[, keep_revel, drop=FALSE]
  names(REVEL.slim) <- sub("^ACMG$", "REVEL.ACMG", names(REVEL.slim))
  names(REVEL.slim) <- sub("^ACMG\\.points$", "REVEL.ACMG.points", names(REVEL.slim))
  
  # juntar por mutação, preservando ordem do user.var
  out <- merge(user.var, AM.slim,    by.x = mut_user, by.y = "mutation", all.x = TRUE)
  out <- merge(out,     REVEL.slim,  by.x = mut_user, by.y = "mutation", all.x = TRUE)
  
  # garantir colunas de rótulos existem mesmo se ausentes
  if (!("AM.ACMG" %in% names(out)))           out$AM.ACMG <- NA_character_
  if (!("REVEL.ACMG" %in% names(out)))        out$REVEL.ACMG <- NA_character_
  if (!("MetaRNN.ACMG" %in% names(out)))      out$MetaRNN.ACMG <- NA_character_
  if (!("AM.ACMG.points" %in% names(out)))    out$AM.ACMG.points <- NA_real_
  if (!("REVEL.ACMG.points" %in% names(out))) out$REVEL.ACMG.points <- NA_real_
  if (!("MetaRNN.ACMG.points" %in% names(out))) out$MetaRNN.ACMG.points <- NA_real_
  
  # consenso por linha
  out$ACMG     <- NA_character_
  out$ACMG.points <- NA_real_
  
  for (i in seq_len(nrow(out))) {
    labs <- c(out$AM.ACMG[i], out$REVEL.ACMG[i], out$MetaRNN.ACMG[i])
    cons <- consensus_acmg(labs)
    out$ACMG[i]     <- cons$label
    out$ACMG.points[i] <- cons$points
  }
  
  out
}