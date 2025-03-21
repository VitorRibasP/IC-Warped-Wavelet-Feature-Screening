library(tidyverse)
library(gridExtra)
library(xtable)
library(corrplot)
library(L1pack)
library(segmented)
library("devtools")
devtools::install_github("rodneyfv/waveletFeatureScreening")

source('./wws_functions.R')

Dados <- read_csv("bancoozonio.csv")

Dados <- Dados[,-1]
Dados1 <- Dados[,-12]
Dados$time <- as.Date(Dados$time)

#Criar tabela com Análise Descritiva
Descritiva <- data.frame(Variável = names(Dados1), 
                         Média = sapply(Dados1, mean),
                         Mediana = sapply(Dados1, median),
                         Desvio = sapply(Dados1, sd),
                         Mínimo = sapply(Dados1, min),
                         Máximo = sapply(Dados1, max))

#Plotar Séries Temporais
Grafico_Linhas <- function(y, legenda, titulo) {
  df <- data.frame(x = Dados$time, y = y)
  ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    ggtitle(titulo) +
    ylab(legenda) +
    xlab("Data") +
    theme_minimal()
}

GO3 <- Grafico_Linhas(Dados$o3, "", "O3")
GCO <- Grafico_Linhas(Dados$co, "", "CO")
GNO <- Grafico_Linhas(Dados$no, "", "NO")
GNO2 <- Grafico_Linhas(Dados$no2, "", "NO2")
GTC <- Grafico_Linhas(Dados$tc, "", "TC")
GRH <- Grafico_Linhas(Dados$rh, "", "RH")
GWS <- Grafico_Linhas(Dados$ws, "", "WS")
GCOSWD <- Grafico_Linhas(Dados$cos_wd, "", "COS_WD")
GSINWD <- Grafico_Linhas(Dados$sin_wd, "", "SIN_WD")
GCOSGWD <- Grafico_Linhas(Dados$cos_gwd, "", "COS_GWD")
GSINGWD <- Grafico_Linhas(Dados$sin_gwd, "", "SIN_GWD")

grid.arrange(GO3, GCO, GNO, GNO2, GTC, GRH, GWS, GCOSWD, GSINWD, GCOSGWD, GSINGWD)

#Séries Temporais junto com O3
Grafico_LinhasO3 <- function(y, titulo, legenda) {
  df <- data.frame(x = Dados$time, o3 = Dados$o3, y = y)
  
  # Encontra os limites dos eixos
  lim_o3 <- range(df$o3, na.rm = TRUE)
  lim_y <- range(df$y, na.rm = TRUE)
  
  # Transforma os dados de y para a escala de o3
  scale_factor <- (lim_o3[2] - lim_o3[1]) / (lim_y[2] - lim_y[1])
  df$y_trans <- (df$y - lim_y[1]) * scale_factor + lim_o3[1]
  
  ggplot(df, aes(x = x)) +
    geom_line(aes(y = o3, color = "O3"), show.legend = FALSE) +
    geom_line(aes(y = y_trans, color = "Y"), show.legend = FALSE) +
    scale_y_continuous(
      name = "O3",
      sec.axis = sec_axis(~ (. - lim_o3[1]) / scale_factor + lim_y[1], name = legenda)
    ) +
    labs(title = titulo, x = "Data", y = "O3") +
    scale_color_manual(
      name = "Legenda",
      values = c("O3" = "blue", "Y" = "red"),
      labels = c("O3", legenda)
    ) +
    theme_minimal() +
    theme(
      axis.title.y.right = element_text(color = "red"),
      axis.title.y.left = element_text(color = "blue"),
      legend.position = "bottom"
    )
}

GCOO3 <- Grafico_LinhasO3(Dados$co, "CO", "CO")
GNOO3 <- Grafico_LinhasO3(Dados$no, "NO", "NO")
GNO2O3 <- Grafico_LinhasO3(Dados$no2, "NO2", "NO2")
GTCO3 <- Grafico_LinhasO3(Dados$tc, "TC", "TC")
GRHO3 <- Grafico_LinhasO3(Dados$rh, "RH", "RH")
GWSO3 <- Grafico_LinhasO3(Dados$ws, "WS", "WS")
GCOSWDO3 <- Grafico_LinhasO3(Dados$cos_wd, "COS_WD", "COS_WD")
GSINWDO3 <- Grafico_LinhasO3(Dados$sin_wd, "SIN_WD", "SIN_WD")
GCOSGWDO3 <- Grafico_LinhasO3(Dados$cos_gwd, "COS_GWD", "COS_GWD")
GSINGWDO3 <- Grafico_LinhasO3(Dados$sin_gwd, "SIN_GWD", "SIN_GWD")

grid.arrange(GCOO3, GNOO3, GNO2O3, GTCO3, GRHO3, GWSO3, GCOSWDO3, GSINWDO3, GCOSGWDO3, GSINGWDO3)


#Transformar e tirar dependência temporal
#O objetivo era detalhar os estudos feitos por 
#Link do artigo: https://www.tandfonline.com/doi/full/10.1080/10618600.2024.2342984?scroll=top&needAccess=true#abstract
#Com isso, foi copiado o código que calcula as diferenças e tira e dependência temporal

# response: Month averages of ozone level
O3 = diff(diff(Dados$o3),11)
# sample size of the period without NAs
n <- length(O3)
# number of covariates (minus two because one is date and another is the 
# response o3)
number_covariates <- dim(Dados)[2] - 2

# matrix of features with original values (before taking differences)
X_raw <- Dados[,c(2:11)]

# matrix of features after computing differences
X = matrix(0, nrow = n, ncol = number_covariates)
colnames(X) <- colnames(X_raw)
k = 1
for(i in 1:number_covariates){
  # variables that need differences at lags 1 and 12
  if(colnames(X)[i] %in% c("co","no","tc")){
    colvar <- diff(diff(X_raw[,i]),11)
  }
  # variables that need differences at the lags 1
  if(colnames(X)[i] %in% c("no2","rh","ws","sin_wd","cos_wd")){
    colvar <- diff(X_raw[,i])
  }
  if(colnames(X)[i] %in% c("cos_gwd", "sin_gwd")){
    # remaining variables do not need differences
    colvar <- X_raw[,i]
  }
  
  colvar <- colvar[(length(colvar)-n+1):length(colvar)]
  X[,i] <- (colvar - min(colvar))/(max(colvar) - min(colvar))
}
X <- as.data.frame(X)
#Com isso, os dados estão transformados
# O3 é o vetor do O3 sem temporalidade
# X é a matriz com as covariáveis sem temporalidade

########################################################################
#DADOS ORIGINAIS
########################################################################

#Mudanças de Fase
MFO3 <- Dados$o3
MFX <- Dados[,2:11]

Covariancias_Originais_Fase <- matrix(NA, 10, 9)

for (i in 1:10) {
  Covariancias_Originais_Fase[i, 1] <- cov(MFO3, MFX[,i])
}

for (i in 2:7) {
  #Mudança de fase de 1 mês
  MFO3 <- MFO3[-length(MFO3)]
  MFX <- MFX[-1,]
  
  for (j in 1:10) {
    Covariancias_Originais_Fase[j,i] <- cov(MFO3, MFX[,j])
  }
}

for (i in 1:10) {
  Covariancias_Originais_Fase[i,8] <- max(abs(Covariancias_Originais_Fase[i,1:7]))
}

Covariancias_Originais_Fase[,9] <- c(1,1,2,1,2,0,4,0,1,4)
rownames(Covariancias_Originais_Fase) <- c(colnames(Dados)[2:11])
colnames(Covariancias_Originais_Fase) <- c(
  "Padrão", "1 mês",paste(2:6, "meses"), "Máximo", "Melhor fase")
#################################################################
Correlacoes_Originais_Fase <- Covariancias_Originais_Fase


for (i in 1:10) {
  Correlacoes_Originais_Fase[i, 1] <- cor(MFO3, MFX[,i])
}

for (i in 2:7) {
  #Mudança de fase de 1 mês
  MFO3 <- MFO3[-length(MFO3)]
  MFX <- MFX[-1,]
  
  for (j in 1:10) {
    Correlacoes_Originais_Fase[j,i] <- cor(MFO3, MFX[,j])
  }
}

for (i in 1:10) {
  Correlacoes_Originais_Fase[i,8] <- max(abs(Correlacoes_Originais_Fase[i,1:7]))
}

Correlacoes_Originais_Fase[,9] <- c(1,1,2,1,2,0,4,0,1,4)

#Gráfico de linhas com mudança de fase

#Séries Temporais Defasadas junto com O3
Grafico_Linhas_Mudanca_Fase <- function(y, defasagem,titulo, legenda) {
  df <- data.frame(x = Dados$time, o3 = Dados$o3, 
                   y = c(y[-(1:defasagem)], rep(NA, defasagem)))
  
  # Encontra os limites dos eixos
  lim_o3 <- range(df$o3, na.rm = TRUE)
  lim_y <- range(df$y, na.rm = TRUE)
  
  # Transforma os dados de y para a escala de o3
  scale_factor <- (lim_o3[2] - lim_o3[1]) / (lim_y[2] - lim_y[1])
  df$y_trans <- (df$y - lim_y[1]) * scale_factor + lim_o3[1]
  
  ggplot(df, aes(x = x)) +
    geom_line(aes(y = o3, color = "O3"), show.legend = FALSE) +
    geom_line(aes(y = y_trans, color = "Y"), show.legend = FALSE) +
    scale_y_continuous(
      name = "O3",
      sec.axis = sec_axis(~ (. - lim_o3[1]) / scale_factor + lim_y[1], name = legenda)
    ) +
    labs(title = titulo, x = "Data", y = "O3") +
    scale_color_manual(
      name = "Legenda",
      values = c("O3" = "blue", "Y" = "red"),
      labels = c("O3", legenda)
    ) +
    theme_minimal() +
    theme(
      axis.title.y.right = element_text(color = "red"),
      axis.title.y.left = element_text(color = "blue"),
      legend.position = "bottom"
    )
}

GMFCO <- Grafico_Linhas_Mudanca_Fase(Dados$co, 1, "CO", "CO")
GMFNO <- Grafico_Linhas_Mudanca_Fase(Dados$no, 1, "NO", "NO")
GMFNO2 <- Grafico_Linhas_Mudanca_Fase(Dados$no2, 2, "NO2", "NO2")
GMFTC <- Grafico_Linhas_Mudanca_Fase(Dados$tc, 1, "TC", "TC")
GMFRH <- Grafico_Linhas_Mudanca_Fase(Dados$rh, 2, "RH", "RH")
GMFWS <- Grafico_LinhasO3(Dados$ws, "WS", "WS")
GMFCOSWD <- Grafico_Linhas_Mudanca_Fase(Dados$cos_wd, 4, "COS_WD", "COS_WD")
GMFSINWD <- Grafico_LinhasO3(Dados$sin_wd, "SIN_WD", "SIN_WD")
GMFCOSGWD <- Grafico_Linhas_Mudanca_Fase(Dados$cos_gwd, 1, "COS_GWD", "COS_GWD")
GMFSINGWD <- Grafico_Linhas_Mudanca_Fase(Dados$sin_gwd, 4, "SIN_GWD", "SIN_GWD")

grid.arrange(GMFCO, GMFNO, GMFNO2, GMFTC, GMFRH, GMFWS, GMFCOSWD, GMFSINWD, GMFCOSGWD, GMFSINGWD)

######################################################################
#DADOS TRANSFORMADOS
######################################################################

MFO3 <- O3
MFX <- X

Covariancias_Transformados_Fase <- Covariancias_Originais_Fase

for (i in 1:10) {
  Covariancias_Transformados_Fase[i, 1] <- cov(MFO3, MFX[,i])
}

for (i in 2:7) {
  #Mudança de fase de 1 mês
  MFO3 <- MFO3[-length(MFO3)]
  MFX <- MFX[-1,]
  
  for (j in 1:10) {
    Covariancias_Transformados_Fase[j,i] <- cov(MFO3, MFX[,j])
  }
}

for (i in 1:10) {
  Covariancias_Transformados_Fase[i,8] <- max(abs(Covariancias_Transformados_Fase[i,1:7]))
}

Covariancias_Transformados_Fase[,9] <- c(6,6,3,2,0,6,0,5,6,4)

#######################################################################

Correlacoes_Transformados_Fase <- Covariancias_Originais_Fase

for (i in 1:10) {
  Correlacoes_Transformados_Fase[i, 1] <- cor(MFO3, MFX[,i])
}

for (i in 2:7) {
  #Mudança de fase de 1 mês
  MFO3 <- MFO3[-length(MFO3)]
  MFX <- MFX[-1,]
  
  for (j in 1:10) {
    Correlacoes_Transformados_Fase[j,i] <- cor(MFO3, MFX[,j])
  }
}

for (i in 1:10) {
  Correlacoes_Transformados_Fase[i,8] <- max(abs(Correlacoes_Transformados_Fase[i,1:7]))
}

Correlacoes_Transformados_Fase[,9] <- c(6,6,3,2,0,6,0,5,6,4)

GMFCO <- Grafico_Linhas_Mudanca_Fase(Dados$co, 6, "CO", "CO")
GMFNO <- Grafico_Linhas_Mudanca_Fase(Dados$no, 6, "NO", "NO")
GMFNO2 <- Grafico_Linhas_Mudanca_Fase(Dados$no2, 3, "NO2", "NO2")
GMFTC <- Grafico_Linhas_Mudanca_Fase(Dados$tc, 2, "TC", "TC")
GMFRH <- Grafico_LinhasO3(Dados$rh, "RH", "RH")
GMFWS <- Grafico_Linhas_Mudanca_Fase(Dados$ws, 6, "WS", "WS")
GMFCOSWD <- Grafico_LinhasO3(Dados$cos_wd, "COS_WD", "COS_WD")
GMFSINWD <- Grafico_Linhas_Mudanca_Fase(Dados$sin_wd, 5, "SIN_WD", "SIN_WD")
GMFCOSGWD <- Grafico_Linhas_Mudanca_Fase(Dados$cos_gwd,6, "COS_GWD", "COS_GWD")
GMFSINGWD <- Grafico_Linhas_Mudanca_Fase(Dados$sin_gwd, 4, "SIN_GWD", "SIN_GWD")

grid.arrange(GMFCO, GMFNO, GMFNO2, GMFTC, GMFRH, GMFWS, GMFCOSWD, GMFSINWD, GMFCOSGWD, GMFSINGWD)

#Criando o Data Frame com as 71 Variáveis

Vetor_Aleatorio_Completo_71 <- as.data.frame(matrix(data = 0, nrow = 42, ncol = 71))
Vetor_Aleatorio_Completo_71[,1] <- O3
colnames(Vetor_Aleatorio_Completo_71)[1] <- "O3"
for (i in 1:10) {
  Vetor_Aleatorio_Completo_71[,2+7*(i-1)] <- X[,i]
  colnames(Vetor_Aleatorio_Completo_71)[2+7*(i-1)] <- colnames(X)[i]
}

X_Defasado <- X
for (j in 1:6) {
  X_Defasado <- X_Defasado[-1,]
  X_Defasado <- rbind(X_Defasado, NA)
  for (i in 1:10) {
    Vetor_Aleatorio_Completo_71[,2+7*(i-1)+j] <- X_Defasado[,i]
    colnames(Vetor_Aleatorio_Completo_71)[2+7*(i-1)+j] <- paste0(colnames(X)[i], " ", j)
  }
}

#Matriz de Correlações 71x71
Matriz_Correlacoes_71P <- matrix(nrow = 71, ncol = 71)
for (i in 1:71) {
  for (j in i:71) {
    Matriz_Correlacoes_71P[i,j] <- cor(Vetor_Aleatorio_Completo_71[,i],
                                 Vetor_Aleatorio_Completo_71[,j], 
                                 use = "complete.obs")
    Matriz_Correlacoes_71P[j,i] <- Matriz_Correlacoes_71P[i,j]
  }
}
colnames(Matriz_Correlacoes_71P) <- names(Vetor_Aleatorio_Completo_71)
rownames(Matriz_Correlacoes_71P) <- names(Vetor_Aleatorio_Completo_71)

Matriz_Correlacoes_71S <- matrix(nrow = 71, ncol = 71)

for (i in 1:71) {
  for (j in i:71) {
    Matriz_Correlacoes_71S[i,j] <- cor(Vetor_Aleatorio_Completo_71[,i],
                                       Vetor_Aleatorio_Completo_71[,j], 
                                       use = "complete.obs", method = "spearman")
    Matriz_Correlacoes_71S[j,i] <- Matriz_Correlacoes_71S[i,j]
  }
}
colnames(Matriz_Correlacoes_71S) <- names(Vetor_Aleatorio_Completo_71)
rownames(Matriz_Correlacoes_71S) <- names(Vetor_Aleatorio_Completo_71)

corrplot(Matriz_Correlacoes_71P, type = "upper", diag = FALSE, tl.cex = 0.6)
corrplot(Matriz_Correlacoes_71S, type = "upper", diag = FALSE, tl.cex = 0.6)


CorrelacoesO3 <- data.frame(Pearson = abs(Matriz_Correlacoes_71P[1,]),
                 Spearman = abs(Matriz_Correlacoes_71S[1,]))


CorrelacoesO3 %>% ggplot() +
  geom_point(aes(x = Pearson, y = Spearman)) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(a = 0, b = 1, col = "red") +
  theme_minimal()

################################################################################

#Regressão por partes L2
l2partes <- function(Y, X, corte) {
  data <- data.frame("O3" = Y,
                     "Var" = X)
  if (length(corte) == 1) {
    cortes <- corte[[1]]
    data$X1 <- cortes + (X - cortes)*(X <= cortes)
    data$X2 <- (X - cortes)*(X > cortes)
    modelo <- lm(Y ~X1+X2, data = data)
    coefs <- c(coef(modelo)[1], coef(modelo)[2],
               coef(modelo)[1] + cortes*(coef(modelo)[2] - coef(modelo)[3]),
               coef(modelo)[3])
    resposta <- (coefs[1] + coefs[2]*X)*(X <= cortes) + 
                (coefs[3] + coefs[4]*X)*(X > cortes)
    return(resposta)
  }
  if (length(corte) == 2) {
    cortes <- c(corte[[1]], corte[[2]])
    
    data$X1 <- ifelse(X <= cortes[1], X, cortes[1])
    data$X2 <- ifelse(X > cortes[1] & X <= cortes[2], X - cortes[1], 0)
    data$X3 <- ifelse(X > cortes[2], X - cortes[2], 0)

    modelo <- lm(Y ~ X1 + X2 + X3, data = data)

    coefs <- coef(modelo)
    
    resposta <- (coefs[1] + coefs[2] * X) * (X <= cortes[1]) + 
      (coefs[1] + coefs[2] * cortes[1] + coefs[3] * (X - cortes[1])) * (X > cortes[1] & X <= cortes[2]) + 
      (coefs[1] + coefs[2] * cortes[1] + coefs[3] * (cortes[2] - cortes[1]) + coefs[4] * (X - cortes[2])) * (X > cortes[2])
    
    return(resposta)
  }
  }

l1partes <- function(Y, X, corte) {
  data <- data.frame("O3" = Y,
                     "Var" = X)
  if (length(corte) == 1) {
    cortes <- corte[[1]]
    data$X1 <- cortes + (X - cortes)*(X <= cortes)
    data$X2 <- (X - cortes)*(X > cortes)
    modelo <- lad(Y ~X1+X2, data = data)
    coefs <- c(coef(modelo)[1], coef(modelo)[2],
               coef(modelo)[1] + cortes*(coef(modelo)[2] - coef(modelo)[3]),
               coef(modelo)[3])
    resposta <- (coefs[1] + coefs[2]*X)*(X <= cortes) + 
      (coefs[3] + coefs[4]*X)*(X > cortes)
    return(resposta)
  }
  if (length(corte) == 2) {
    cortes <- c(corte[[1]], corte[[2]])
    
    data$X1 <- ifelse(X <= cortes[1], X, cortes[1])
    data$X2 <- ifelse(X > cortes[1] & X <= cortes[2], X - cortes[1], 0)
    data$X3 <- ifelse(X > cortes[2], X - cortes[2], 0)
    
    modelo <- lad(Y ~ X1 + X2 + X3, data = data)
    
    coefs <- coef(modelo)
    
    resposta <- (coefs[1] + coefs[2] * X) * (X <= cortes[1]) + 
      (coefs[1] + coefs[2] * cortes[1] + coefs[3] * (X - cortes[1])) * (X > cortes[1] & X <= cortes[2]) + 
      (coefs[1] + coefs[2] * cortes[1] + coefs[3] * (cortes[2] - cortes[1]) + coefs[4] * (X - cortes[2])) * (X > cortes[2])
    
    return(resposta)
  }
}

#Regressões Por erro quadrado vs erro absoluto
DadosFiltrados <- Vetor_Aleatorio_Completo_71[c("O3", "no 6", "rh", "ws 5",
                                                "tc 1", "sin_gwd 4", "tc",
                                                "no2 2", "co 4", "cos_gwd 2",
                                                "no 4")]
diffWavelets <- matrix(NA, 10, 41)

for (i in 2:11) {
  data <- data.frame("O3" = DadosFiltrados$O3,
                     "Var" = DadosFiltrados[,i])
  data <- na.omit(data)
  n <- length(data$Var)
  regWavelets <- warped_wavelet_reg(x=data$Var, y=data$O3)
  
  x1_sorted <- sort(data$Var, decreasing = FALSE, index.return=TRUE)
  O3_Wavelet <- rep(NA, n)
  O3_Wavelet[x1_sorted$ix] <- regWavelets$y_fit
  O3_Wavelet <- c(O3_Wavelet, rep(NA, length(DadosFiltrados$O3) - n))
  
  data <- data.frame("Var" = DadosFiltrados[,i],
                     "Wavelet" = O3_Wavelet)
  data <- data %>% arrange(Var)
  
  diffWavelets[i-1,] <- diff(data$Wavelet)/diff(data$Var)
}

rownames(diffWavelets) <- colnames(DadosFiltrados)[2:11]
colnames(diffWavelets) <- paste0(c(1:41), "-", c(2:42))

cortes <- list(
  quantile(DadosFiltrados$`no 6`, 0.6, na.rm = TRUE),
  c(quantile(DadosFiltrados$rh, 0.3, na.rm = TRUE), quantile(DadosFiltrados$rh, 0.6, na.rm = TRUE)),
  quantile(DadosFiltrados$`ws 5`, 0.6, na.rm = TRUE),
  quantile(DadosFiltrados$`tc 1`, 0.4, na.rm = TRUE),
  c(quantile(DadosFiltrados$`sin_gwd 4`, 0.2, na.rm = TRUE), quantile(DadosFiltrados$`sin_gwd 4`, 0.9, na.rm = TRUE)),
  quantile(DadosFiltrados$tc, 0.3, na.rm = TRUE),
  quantile(DadosFiltrados$`no2 2`, 0.6, na.rm = TRUE),
  quantile(DadosFiltrados$`co 4`, 0.6, na.rm = TRUE),
  quantile(DadosFiltrados$`cos_gwd 2`, 0.6, na.rm = TRUE),
  quantile(DadosFiltrados$`no 4`, 0.5, na.rm = TRUE)
)

TabelaRegFiltrados <- data.frame("Variável" = colnames(DadosFiltrados)[2:11],
                                 "Beta0Q" = rep(NA, 10),
                                 "Beta1Q" = rep(NA, 10),
                                 "Beta0A" = rep(NA, 10),
                                 "Beta1A" = rep(NA, 10))
for (i in 2:11) {
  Var <- DadosFiltrados[,i]
  #Regressão Linear com Erro Quadrático
  RegQuadrado <- lm(O3 ~ Var, data = DadosFiltrados)
  TabelaRegFiltrados[i-1,2] <- coef(RegQuadrado)[1]
  TabelaRegFiltrados[i-1,3] <- coef(RegQuadrado)[2]
  
  #Regressão por Partes L2
  RegPartesL2 <- l2partes(DadosFiltrados$O3, Var, cortes[[i-1]])
  
  #Regressão por Partes L1
  RegPartesL1 <- l1partes(DadosFiltrados$O3, Var, cortes[[i-1]])
  
  #Regressão Linear com Erro Absoluto
  RegAbsoluto <- lad(DadosFiltrados$O3 ~ DadosFiltrados[,i])
  TabelaRegFiltrados[i-1,4] <- coef(RegAbsoluto)[1]
  TabelaRegFiltrados[i-1,5] <- coef(RegAbsoluto)[2]  
  
  data <- data.frame("O3" = DadosFiltrados$O3,
                     "Var" = DadosFiltrados[,i])
  data <- na.omit(data)
  
  #Regressão por Wavelets
  n <- length(data$Var)
  regWavelets <- warped_wavelet_reg(x=data$Var, y=data$O3)
  
  x1_sorted <- sort(data$Var, decreasing = FALSE, index.return=TRUE)
  O3_Wavelet <- rep(NA, n)
  O3_Wavelet[x1_sorted$ix] <- regWavelets$y_fit
  
  O3_Wavelet <- c(O3_Wavelet, rep(NA, length(DadosFiltrados$O3) - n))
  
  print(data.frame("O3" = DadosFiltrados$O3,
                   "Var" = DadosFiltrados[,i],
                   "O3_Wavelet" = O3_Wavelet,
                   "L2Partes" = RegPartesL2,
                   "L1Partes" = RegPartesL1) %>% 
    ggplot(aes(x = Var, y = O3)) +
      
    #Regressão Linear Erro Quadrático
    geom_abline(aes(slope = coef(RegQuadrado)[2], intercept = coef(RegQuadrado)[1],
                color = "L2"), size = 1.25) +
    #Regressão Linear Erro Absoluto
    geom_abline(aes(slope = coef(RegAbsoluto)[2] , intercept = coef(RegAbsoluto)[1],
                color = "L1"), size = 1.25) +
    #Regressão por Wavelets
    geom_line(aes(x = Var, y = O3_Wavelet, color = "Warped Wavelets"), size = 1.5) +
    
    #Regressão por Partes Erro Quadrático
    geom_line(aes(x = Var, y = L2Partes, color = "Piecewise L2"), size = 1.25) +
    
    #Regressão por Partes Erro Absoluto
    geom_line(aes(x = Var, y = L1Partes, color = "Piecewise L1"), size = 1.25) +
        
    geom_point(size = 3) +
    scale_color_manual(values = c("L2" = "red", "L1" = "blue",
                                  "Warped Wavelets" = "darkgreen",
                                  "Piecewise L1" = "orange", "Piecewise L2" = "orchid4")) +
    labs(x = TabelaRegFiltrados[i-1,1], color = "Fitting") +
    theme_minimal())
}
