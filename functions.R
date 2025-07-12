surv_exploratoria <- function(df, variavel, tipo = c("numerica", "categorica"),
                              metodo_discretizacao = "jenks", n_classes = 4,
                              tempo = "LENFOL", censura = "FSTAT") {
  
  tipo <- match.arg(tipo)
  v <- df[[variavel]]
  
  if (tipo == "numerica") {
    breaks <- classIntervals(v, n = n_classes, style = metodo_discretizacao, intervalClosure = "left")$brks
    df$grupo <- cut(v, breaks = breaks, include.lowest = TRUE)
    cat("Intervalos criados:\n")
    print(breaks)
  } else {
    df$grupo <- as.factor(v)
    cat("Frequências e proporções:\n")
    freq <- table(df$grupo)
    print(cbind(Freq = freq, Prop = round(prop.table(freq), 3)))
  }
  
  levels(df$grupo) <- gsub("^.*?:\\s*", "", levels(df$grupo))
  
  # Tabela de frequência
  tab <- as.data.frame(table(df$grupo))
  colnames(tab) <- c("Grupo", "Frequência")
  tab$Proporção <- tab$Frequência / sum(tab$Frequência)
  
  # Gráfico de distribuição
  if (tipo == "numerica") {
    p_dist <- ggplot(df, aes(x = .data[[variavel]])) +
      geom_histogram(fill = "#4B9CD3", color = "black", bins = 30, alpha = 0.8) +
      labs(title = paste("Histograma de", variavel),
           x = variavel, y = "Frequência") +
      theme_minimal(base_size = 10)
  } else {
    p_dist <- ggplot(tab, aes(x = Grupo, y = Frequência, fill = Grupo)) +
      geom_bar(stat = "identity", color = "black", alpha = 0.9) +
      geom_text(aes(label = scales::percent(Proporção, accuracy = 0.1)), 
                vjust = -0.4, size = 4.5) +
      scale_fill_npg() +
      labs(title = paste("Distribuição de", variavel),
           y = "Frequência Absoluta", x = variavel, fill = NULL) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none") +
      ylim(0, max(tab$Frequência) * 1.1)
  }
  
  # Gráfico de tempo de sobrevivência em relação as variáveis explicativas
  if (tipo == "numerica") {
    p_time <- ggplot(df, aes(x = .data[[variavel]], y = .data[[tempo]])) +
      geom_point(color = "#E76F51", alpha = 0.6) +
      geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed") +
      labs(title = paste("Dispersão de", tempo, "por", variavel),
           x = variavel, y = "Tempo") +
      theme_minimal(base_size = 10)
  } else {
    p_time <- ggplot(df, aes(x = grupo, y = .data[[tempo]])) +
      geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.color = "black") +
      labs(title = paste("Boxplot do tempo por", variavel),
           x = variavel, y = "Tempo") +
      theme_minimal(base_size = 10)
  }
  
  # Curvas de sobrevivência
  kms <- survfit(as.formula(paste0("Surv(", tempo, ",", censura, ") ~ grupo")), data = df)
  names(kms$strata) <- gsub("^grupo=", "", names(kms$strata))
  
  p_surv <- ggsurvplot(
    kms,
    conf.int = FALSE,
    pval = FALSE,
    risk.table = FALSE,
    legend.title = variavel,
    data = df,
    palette = "npg",
    xlab = "Tempo",
    ylab = expression("S(t)"),
    ggtheme = theme_minimal(base_size = 10),
    title = paste("Curvas de Sobrevivência por", variavel)
  )
  
  # Testes estatísticos
  logrank <- survdiff(as.formula(paste0("Surv(", tempo, ",", censura, ") ~ grupo")), data = df, rho = 0)
  wilcoxon <- survdiff(as.formula(paste0("Surv(", tempo, ",", censura, ") ~ grupo")), data = df, rho = 1)
  
  return(list(
    tabela_frequencia = tab,
    grafico_distribuicao = p_dist,
    grafico_tempo = p_time,
    grafico_sobrevivencia = p_surv$plot,
    logrank_test = logrank,
    wilcoxon_test = wilcoxon
  ))
}

calc_ttt <- function(tempo) {
  x <- sort(tempo)
  n <- length(x)
  if (n < 2) return(NULL)
  total <- sum(x)
  G <- sapply(1:n, function(r) (sum(x[1:r]) + (n - r) * x[r]) / total)
  data.frame(
    rn = (1:n)/n,
    G = G
  )
}


# Função para calcular TTT
calc_ttt_grupo <- function(tempo, grupo) {
  grupo <- as.factor(grupo)
  res <- lapply(levels(grupo), function(g) {
    x <- sort(tempo[grupo == g])
    n <- length(x)
    if (n < 2) return(NULL)
    total <- sum(x)
    G <- sapply(1:n, function(r) (sum(x[1:r]) + (n - r) * x[r]) / total)
    data.frame(
      rn = (1:n)/n,
      G = G,
      grupo = g
    )
  })
  do.call(rbind, res)
}

# Função para gerar um gráfico TTT com ggplot2 para uma variável
plot_ttt_gg <- function(data, tempo_var, grupo_var, status_var = NULL, tipo = c("numerica","categorica")) {
  
  tipo <- match.arg(tipo)
  
  if (tipo == "numerica") {
    
    v <- data[[grupo_var]]
    
    breaks <- classIntervals(v, n = 4, style = 'jenks', intervalClosure = "left")$brks
    data$grupo_var_jenks <- cut(v, breaks = breaks, include.lowest = TRUE)
    
    grupo <- data$grupo_var_jenks
    
    cat("Intervalos criados:\n")
    print(breaks)
  } else {grupo <- data[[grupo_var]]}
  
  tempo <- data[[tempo_var]]
  status <- if (!is.null(status_var)) data[[status_var]] else rep(1, length(tempo))
  
  
  df_ttt <- calc_ttt_grupo(tempo, grupo)
  
  ggplot(df_ttt, aes(x = rn, y = G, color = grupo)) +
    geom_line(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_npg() +
    labs(
      title = paste("Curva TTT por", grupo_var),
      x = expression(r/n),
      y = expression(G(r/n)),
      color = grupo_var
    ) +
    theme_minimal(base_size = 12)
}

# Função principal para gerar e salvar grid de TTTs
plot_ttt_grid <- function(data, tempo_var, status_var = NULL, vars_grupo, ncol = 2, output = "ttt_grid.pdf", tipo = c("numerica","categorica")) {
  
  
  plots <- lapply(vars_grupo, function(v) {
    plot_ttt_gg(data, tempo_var = tempo_var, grupo_var = v, status_var = status_var, tipo = tipo)
  })
  
  grid_plot <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(title = "Curvas TTT para variáveis categóricas")
  
  ggsave(output, grid_plot, width = 10, height = ceiling(length(plots)/ncol) * 4)
}


extrair_resumo <- function(nome, obj) {
  tibble(
    Variável = nome,
    Chisq = round(obj$chisq,3),
    df = length(obj$obs)-1,
    p_valor = round(obj$pvalue,4)
  )
}