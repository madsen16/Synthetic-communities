setwd("C:/directory")

library(phyloseq)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(GUniFrac)
library(rgl)

load("seqcenter_merge_phyloseq_cluster99.rds")

existing_metadata <- sample_data(p99)
new_column <- read.csv("C:/directory/Metadata_seqcenter_final_processing_v3.csv")$Engraft
existing_metadata$Engraft <- new_column
sample_data(p99) <- existing_metadata

# Extract the sample names from the existing sample data
existing_sample_names <- rownames(sample_data(p99))

# Initialize a vector to store the modified sample names
new_sample_names <- character(length(existing_sample_names))

# Loop over the sample names and modify them individually
for (i in seq_along(existing_sample_names)) {
  # Extract the part before the last underscore
  name_without_suffix <- sub("[A-Z][0-9]+(_[A-Za-z0-9]+)?$", "", existing_sample_names[i])
  
  # Find the corresponding Media.type and Replicate values
  media_type <- sample_data(p99)[existing_sample_names[i], "Media.type"]
  replicate <- sample_data(p99)[existing_sample_names[i], "Replicate"]
  
  # Create the new name
  modified_name <- paste0(name_without_suffix, media_type, "_", replicate)
  
  # Store the modified name in the vector
  new_sample_names[i] <- modified_name
}

# Assign the modified sample names to the phyloseq object
sample_data(p99)$Modified_Sample_Names <- new_sample_names

unifrac_dist_unweight <- UniFrac(p99, weighted = FALSE)

unifrac_dist_weight <- UniFrac(p99, weighted = TRUE)

jaccard_dist <- phyloseq::distance(p99, method = "jaccard")

bray_curtis_dist <- phyloseq::distance(p99, method = "bray")

pca <- prcomp(jaccard_dist)

pca2 <- prcomp(bray_curtis_dist)

pca3 <- prcomp(unifrac_dist_weight)

pca4 <- prcomp(unifrac_dist_unweight)

###PCA for all samples

library(plotly)

sample_data <- sample_data(p99)

unique_samples <- unique(sample_data$Modified_Sample_Names)
n_samples <- length(unique_samples)
colors <- rainbow(n_samples)
sample_colors <- colors[match(sample_data$Modified_Sample_Names, unique_samples)]

###change PCA number for desired metric
pca_data <- data.frame(
  PC1 = pca3$x[, 1],
  PC2 = pca3$x[, 2],
  PC3 = pca3$x[, 3],
  Sample = sample_data$Modified_Sample_Names,
  Color = colors
)

p <- plot_ly(pca_data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Sample, text = ~Sample, type = "scatter3d", mode = "markers")
p <- p %>% layout(title = "Weighted UniFrac Distance")
p


###plotting for controlled visual

library(plotly)
library(shiny)
library(htmlwidgets)
library(MASS)
library(cluster)
library(ggdist)
library(ggplot2)
library(viridis)


sample_data <- sample_data(p99)

unique_samples <- unique(sample_data$Sample)
n_samples <- length(unique_samples)
colors <- rainbow(n_samples)
sample_colors <- colors[match(sample_data$Sample, unique_samples)]


###change PCA number for desired metric
pca_data <- data.frame(
  PC1 = pca3$x[, 1],
  PC2 = pca3$x[, 2],
  PC3 = pca3$x[, 3],
  Sample = sample_data$Sample,
  Color = colors
)

pca_data_with_metadata <- merge(pca_data, sample_data, by = "Sample", all.x = TRUE)

pca_data_with_metadata <- pca_data
pca_data_with_metadata$Soil <- sample_data$Soil[match(pca_data_with_metadata$Sample, sample_data$Sample)]
pca_data_with_metadata$Media.type <- sample_data$Media.type[match(pca_data_with_metadata$Sample, sample_data$Sample)]
pca_data_with_metadata$Engraft <- sample_data$Engraft[match(pca_data_with_metadata$Sample, sample_data$Sample)]
pca_data_with_metadata$Well_plate <- sample_data$Well_plate[match(pca_data_with_metadata$Sample, sample_data$Sample)]
pca_data_with_metadata$Modified_Sample_Names <- sample_data$Modified_Sample_Names[match(pca_data_with_metadata$Sample, sample_data$Sample)]


###all relevant groups

ui <- fluidPage(
  checkboxGroupInput("soil_groups", "Select Soil Groups:",
                     choices = unique(pca_data_with_metadata$Soil), selected = character(0)),
  checkboxGroupInput("media_groups", "Select Media Groups:",
                     choices = unique(pca_data_with_metadata$Media.type), selected = character(0)),
  checkboxGroupInput("engraft_groups", "Select Engraft Groups:",
                     choices = unique(pca_data_with_metadata$Engraft), selected = character(0)),
  checkboxGroupInput("replicate_groups", "Select Well Plate Groups:",
                     choices = unique(pca_data_with_metadata$Well_plate), selected = character(0)),
  plotlyOutput("pca_plot")
)

server <- function(input, output) {
  filtered_data <- reactive({
    selected_soils <- input$soil_groups
    selected_media <- input$media_groups
    selected_engraft <- input$engraft_groups
    selected_replicate <- input$replicate_groups
    
    filtered_pca_data <- pca_data_with_metadata[
      pca_data_with_metadata$Soil %in% selected_soils &
        pca_data_with_metadata$Media.type %in% selected_media &
        pca_data_with_metadata$Engraft %in% selected_engraft & 
        pca_data_with_metadata$Well_plate %in% selected_replicate,]
    
    return(filtered_pca_data)
  })
  
  output$pca_plot <- renderPlotly({
    filtered_data_combined <- filtered_data()
    filtered_data_combined$Combined <- paste(
      filtered_data_combined$Soil, 
      filtered_data_combined$Media.type, 
      filtered_data_combined$Engraft,
      filtered_data_combined$Well_plate,
      sep = "-"
    )
    
    p <- plot_ly(filtered_data_combined, x = ~PC1, y = ~PC2, z = ~PC3,
                 color = ~Combined, text = ~Sample, type = "scatter3d", mode = "markers")
    p <- p %>% layout(
      title = list(
        text = "Weighted UniFrac Distance PCA",
        font = list(size = 16)  # Adjust title font size
      ),
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3"),
        bgcolor = "white",  # Background color
        domain = list(x = c(0, 1), y = c(0, 1)),
        xaxis_showline = TRUE,
        xaxis_linecolor = "black",
        xaxis_linewidth = 3,
        yaxis_showline = TRUE,
        yaxis_linecolor = "black",
        yaxis_linewidth = 3,
        zaxis_showline = TRUE,
        zaxis_linecolor = "black",
        zaxis_linewidth = 3,
        showbackground = FALSE,  # Disable scene background to avoid overlap with paper background
        aspectmode = "manual",  # Ensure consistent aspect ratio
        aspectratio = list(x = 1, y = 1, z = 1)  # Set aspect ratio to 1:1:1
      ),
      margin = list(t = 50, b = 10, l = 10, r = 10),  # Adjust title position
      paper_bgcolor = "white",  # Plot area background color
      font = list(color = "black"),  # Font color
      hoverlabel = list(bgcolor = "black", font = list(color = "black")),  # Tooltip font color
      legend = list(
        bgcolor = "white", 
        bordercolor = "black", 
        borderwidth = 1,
        font = list(size = 12)  # Adjust legend font size
      ) 
    )
    return(p)
  })
}

shinyApp(ui, server)

###all relevant groups focusing on specific well identity

ui <- fluidPage(
  checkboxGroupInput("soil_groups", "Select Soil Groups:",
                     choices = unique(pca_data_with_metadata$Soil), selected = character(0)),
  checkboxGroupInput("media_groups", "Select Media Groups:",
                     choices = unique(pca_data_with_metadata$Media.type), selected = character(0)),
  checkboxGroupInput("engraft_groups", "Select Engraft Groups:",
                     choices = unique(pca_data_with_metadata$Engraft), selected = character(0)),
  checkboxGroupInput("replicate_groups", "Select Well Plate Groups:",
                     choices = unique(pca_data_with_metadata$Well_plate), selected = character(0)),
  plotlyOutput("pca_plot")
)

server <- function(input, output) {
  filtered_data <- reactive({
    selected_soils <- input$soil_groups
    selected_media <- input$media_groups
    selected_engraft <- input$engraft_groups
    selected_replicate <- input$replicate_groups
    
    filtered_pca_data <- pca_data_with_metadata[
      pca_data_with_metadata$Soil %in% selected_soils &
        pca_data_with_metadata$Media.type %in% selected_media &
        pca_data_with_metadata$Engraft %in% selected_engraft & 
        pca_data_with_metadata$Well_plate %in% selected_replicate,]
    
    return(filtered_pca_data)
  })
  
  output$pca_plot <- renderPlotly({
    filtered_data_combined <- filtered_data()
    filtered_data_combined$Combined <- paste(
      filtered_data_combined$Soil, 
      filtered_data_combined$Media.type, 
      filtered_data_combined$Engraft,
      filtered_data_combined$Well_plate,
      sep = "-"
    )
    
    # Create the plot
    p <- plot_ly(filtered_data_combined, x = ~PC1, y = ~PC2, z = ~PC3,
                 color = ~Combined, text = ~Modified_Sample_Names, type = "scatter3d", mode = "markers+text",
                 textposition = "top right", textfont = list(size = 10))
    
    # Layout adjustments
    p <- p %>% layout(
      title = list(
        text = "Weighted UniFrac Distance PCA",
        font = list(size = 16)
      ),
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3"),
        bgcolor = "white",
        domain = list(x = c(0, 1), y = c(0, 1)),
        xaxis_showline = TRUE,
        xaxis_linecolor = "black",
        xaxis_linewidth = 3,
        yaxis_showline = TRUE,
        yaxis_linecolor = "black",
        yaxis_linewidth = 3,
        zaxis_showline = TRUE,
        zaxis_linecolor = "black",
        zaxis_linewidth = 3,
        showbackground = FALSE,
        aspectmode = "manual",
        aspectratio = list(x = 1, y = 1, z = 1)
      ),
      margin = list(t = 50, b = 10, l = 10, r = 10),
      paper_bgcolor = "white",
      font = list(color = "black"),
      hoverlabel = list(bgcolor = "white", font = list(color = "black")),
      legend = list(
        bgcolor = "white", 
        bordercolor = "black", 
        borderwidth = 1,
        font = list(size = 12)
      )
    )
    
    return(p)
  })
}

shinyApp(ui, server)
