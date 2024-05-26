#Ossman Mejía Guzmán
#S25 - Unidad 3. Evidencia de aprendizaje 1. 
#Pronóstico de tratamiento para paciente con cáncer de mama - Datos multivariados: 
#relaciones, comportamiento y predicción.
#24 de mayo de 2024
#Asesora: Maria Claudia Negret Lopez

# 1.Realiza un análisis exploratorio de datos a partir de una muestra aleatoria simple 
# representativa par analizar el comportamiento de las variablesa partir de los grupos 
#M (maligno) y B (benigno) (establezca un nivel de confianza de 99%)

# Instalar y cargar las librerías necesarias
install.packages(c("dplyr", "ggplot2", "GGally", "tidyr", "corrplot"))
library(dplyr)
library(ggplot2)
library(GGally)
library(tidyr)
library(corrplot)

# Asumiendo que los datos ya están cargados en un data frame llamado 'data'

# Función para calcular el tamaño de la muestra
calculate_sample_size <- function(N, confidence_level, margin_of_error, p = 0.5) {
  z <- qnorm((1 + confidence_level) / 2)
  n0 <- (z^2 * p * (1 - p)) / (margin_of_error^2)
  n <- n0 / (1 + ((n0 - 1) / N))
  return(ceiling(n))
}

# Parámetros
N <- nrow(data) # Número total de instancias
confidence_level <- 0.99
margin_of_error <- 0.05

# Calcular tamaño de muestra
sample_size <- calculate_sample_size(N, confidence_level, margin_of_error)

# Establecer la semilla para reproducibilidad
set.seed(123)

# Crear una muestra aleatoria simple
sample_data <- data %>% sample_n(sample_size)

# Validar la representatividad de la muestra comparando proporciones
prop_table_original <- prop.table(table(data$V2))
prop_table_sample <- prop.table(table(sample_data$V2))

prop_table_original
prop_table_sample

# Cargar las librerías necesarias
library(knitr)

# Calcular las estadísticas resumidas
summary_stats <- sample_data %>%
  select(V2, V3:V32) %>%
  group_by(V2) %>%
  summarise(across(everything(), list(
    mean = ~mean(.),
    median = ~median(.),
    sd = ~sd(.),
    IQR = ~IQR(.),
    min = ~min(.),
    max = ~max(.)
  ), .names = "{col}_{fn}"))

# Verificar la estructura de summary_stats
print(summary_stats)

# Reorganizar los datos con pivot_longer
summary_stats_long <- summary_stats %>%
  pivot_longer(-V2, names_to = c("variable", "statistic"), names_sep = "_")

# Verificar la estructura de summary_stats_long
print(summary_stats_long)

# Reorganizar los datos con pivot_wider
summary_stats_wide <- summary_stats_long %>%
  pivot_wider(names_from = statistic, values_from = value)

# Verificar la estructura de summary_stats_wide
print(summary_stats_wide)

# Imprimir la tabla completa
kable(summary_stats_wide, digits = 2, caption = "Estadísticas Resumidas por Grupo")


# Convertir los datos a formato largo para facilitar la visualización
sample_data_long <- sample_data %>%
  pivot_longer(cols = V3:V32, names_to = "variable", values_to = "value")

# Crear histogramas de las variables por grupo
ggplot(sample_data_long, aes(x = value, fill = V2)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Histogramas de variables por grupo", x = "Valor", y = "Frecuencia") +
  theme_minimal() +
  theme(legend.title = element_blank())

# Crear boxplots de las variables por grupo
ggplot(sample_data_long, aes(x = V2, y = value, fill = V2)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Boxplots de variables por grupo", x = "Grupo", y = "Valor") +
  theme_minimal() +
  theme(legend.position = "none")

# Crear gráficos de densidad de las variables por grupo
ggplot(sample_data_long, aes(x = value, fill = V2)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Gráficos de densidad de variables por grupo", x = "Valor", y = "Densidad") +
  theme_minimal() +
  theme(legend.title = element_blank())

# Calcular matriz de correlación
correlation_matrix <- cor(sample_data[,3:32])

# Visualizar la matriz de correlación
corrplot(correlation_matrix, method = "circle")

# Identificar correlaciones con malignidad y benignidad
cor_with_diagnosis <- correlation_matrix[,1]  # Columna de diagnóstico
relevant_variables <- which(abs(cor_with_diagnosis) > 0.5 & names(cor_with_diagnosis) != "V2")

# Mostrar las variables más relevantes
relevant_variables

# Crear gráficos de cajas y bigotes
sample_data_long %>%
  filter(variable %in% names(sample_data)[3:32]) %>%
  ggplot(aes(x = V2, y = value, fill = V2)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Boxplots de variables por diagnóstico", x = "Diagnóstico", y = "Valor") +
  theme_minimal()

# Realizar pruebas t para cada variable
p_values <- sapply(names(sample_data)[3:32], function(var) {
  t_test <- t.test(sample_data_long$value[sample_data_long$V2 == "M"], 
                   sample_data_long$value[sample_data_long$V2 == "B"],
                   var.equal = TRUE)
  return(t_test$p.value)
})

# Corregir los valores p para múltiples comparaciones (por ejemplo, método de Bonferroni)
p_values_corrected <- p.adjust(p_values, method = "bonferroni")

# Imprimir los valores p corregidos
print(p_values_corrected)

# Variables más relevantes en términos de correlación con el diagnóstico
relevant_variables <- which(abs(cor_with_diagnosis) > 0.5 & names(cor_with_diagnosis) != "V2")
relevant_variable_names <- names(sample_data)[relevant_variables + 2]  # +2 para ajustar el índice

# Seleccionar solo las variables relevantes
relevant_data <- sample_data[, c(2, relevant_variables + 2)]  # +2 para ajustar el índice

# Función para calcular las tablas de contingencia y las probabilidades de error tipo I y tipo II
calculate_error_probabilities <- function(data, variable_names) {
  error_probabilities <- data.frame(variable = character(), p_value = numeric(), type_I_error = numeric(), type_II_error = numeric())
  
  for (variable_name in variable_names) {
    contingency_table <- table(data$V2, data[[variable_name]])
    chi_square_test <- chisq.test(contingency_table)
    
    p_value <- chi_square_test$p.value
    type_I_error <- p_value
    type_II_error <- 1 - pchisq(chi_square_test$statistic, df = chi_square_test$parameter)
    
    error_probabilities <- rbind(error_probabilities, data.frame(variable = variable_name, p_value = p_value, type_I_error = type_I_error, type_II_error = type_II_error))
  }
  
  return(error_probabilities)
}

# Calcular las probabilidades de error tipo I y tipo II para las variables relevantes
error_probabilities <- calculate_error_probabilities(relevant_data, relevant_variable_names)

# Mostrar las probabilidades de error
print(error_probabilities)

# Identificar el error más grave
most_severe_error <- error_probabilities[which.max(error_probabilities$type_II_error), ]
print(most_severe_error)

# Función para calcular las tablas de contingencia
calculate_contingency_tables <- function(data, variable_names) {
  contingency_tables <- list()
  
  for (variable_name in variable_names) {
    contingency_table <- table(data$V2, data[[variable_name]])
    contingency_tables[[variable_name]] <- contingency_table
  }
  
  return(contingency_tables)
}

# Calcular las tablas de contingencia para las variables relevantes
contingency_tables <- calculate_contingency_tables(relevant_data, relevant_variable_names)

# Mostrar las tablas de contingencia
for (variable_name in names(contingency_tables)) {
  cat("Tabla de contingencia para", variable_name, ":\n")
  print(contingency_tables[[variable_name]])
  cat("\n")
}

