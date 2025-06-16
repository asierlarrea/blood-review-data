# Librerías necesarias
library(dplyr)
library(ggplot2)

# ==== 1. Función para cargar y normalizar ====
load_and_normalize <- function(file_path, source_name) {
  df <- read.csv(file_path)
  df <- df %>%
    mutate(
      norm_conc = (Concentration - min(Concentration, na.rm = TRUE)) / 
        (max(Concentration, na.rm = TRUE) - min(Concentration, na.rm = TRUE)),
      db_name = source_name
    )
  return(df)
}

# ==== 2. Cargar PeptideAtlas ====
peptideatlas <- read.csv("C:\\Users\\HP zBook 15v\\Desktop\\peptide.csv")

peptideatlas <- peptideatlas %>%
  mutate(
    norm_conc = (Concentration - min(Concentration, na.rm = TRUE)) /
      (max(Concentration, na.rm = TRUE) - min(Concentration, na.rm = TRUE)),
    db_name = "peptideatlas"
  ) %>%
  filter(norm_conc > 0) %>%
  arrange(norm_conc) %>%
  mutate(order = row_number())

# ==== 3. Rutas a otras bases de datos ====
other_files <- list(
  "hpa_ms" = "C:\\Users\\HP zBook 15v\\Desktop\\hpa_ms.csv",
  "pax" = "C:\\Users\\HP zBook 15v\\Desktop\\pax.csv",
  "gpmdb" = "C:\\Users\\HP zBook 15v\\Desktop\\gpmdb.csv",
  "hpa_immuno" = "C:\\Users\\HP zBook 15v\\Desktop\\hpa_immuno.csv",
  "hpa_pea" = "C:\\Users\\HP zBook 15v\\Desktop\\hpa_pea.csv"
)

# ==== 4. Procesar y unir otras bases de datos ====
other_dbs <- lapply(names(other_files), function(name) {
  db <- load_and_normalize(other_files[[name]], name) %>%
    filter(norm_conc > 0) %>%
    inner_join(peptideatlas %>% select(Protein, order), by = "Protein")
}) %>% bind_rows()

# ==== 5. Unir todo para graficar ====
combined_data <- bind_rows(
  peptideatlas %>% select(order, norm_conc, db_name),
  other_dbs %>% select(order, norm_conc, db_name)
)

# ==== 6. Crear paleta de colores ====
# Obtener nombres únicos excepto peptideatlas
db_names <- unique(combined_data$db_name)
other_names <- setdiff(db_names, "peptideatlas")

# Crear colores: negro para peptideatlas, automáticos para el resto
colors <- c("peptideatlas" = "black")
colors[other_names] <- scales::hue_pal()(length(other_names))

# ==== 7. Graficar con colores definidos ====
p <- ggplot(combined_data, aes(x = order, y = norm_conc, color = db_name)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = colors) +
  scale_y_log10() +
  labs(
    x = "Proteins (sorted by PeptideAtlas concentration)",
    y = "Normalized concentration (Log escale)",
    color = "Database"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Mostrar el plot en R
print(p)

# Guardar la figura como TIFF
ggsave("C:\\Users\\HP zBook 15v\\Desktop\\protein_plot.tiff", plot = p, device = "tiff",
       width = 10, height = 6, units = "in", dpi = 300)

