setwd("C:/directory")

library(phyloseq)
library (RColorBrewer)
library(ggpubr)
library (ggplot2)
library(dplyr)
library(viridis)
library(plotly)
library(shiny)
library(shinyjs)
library(htmlwidgets)
library(tidyr)
library(reshape2)

# Load your data
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

p99.comp <- microbiome::transform(p99, "compositional")

composition_data <- as.data.frame(otu_table(p99.comp))

tax_data <- as.data.frame(tax_table(p99.comp))

merged_data <- merge(composition_data, tax_data, by = "row.names")

merged_data_longer <- merged_data %>%
  pivot_longer(cols = -c(Row.names, domain, phylum, class, order, family, genus),
               names_to = "Sample",
               values_to = "Abundance")

merged_data_longer$genus[is.na(merged_data_longer$genus)] <- "Unknown"

merged_data_longer$family[is.na(merged_data_longer$family)] <- "Unknown"

compositional_data_with_metadata <- merged_data_longer
compositional_data_with_metadata$Soil <- sample_data$Soil[match(compositional_data_with_metadata$Sample, sample_data$Sample)]
compositional_data_with_metadata$Media.type <- sample_data$Media.type[match(compositional_data_with_metadata$Sample, sample_data$Sample)]
compositional_data_with_metadata$Engraft <- sample_data$Engraft[match(compositional_data_with_metadata$Sample, sample_data$Sample)]
compositional_data_with_metadata$Well_plate <- sample_data$Well_plate[match(compositional_data_with_metadata$Sample, sample_data$Sample)]
compositional_data_with_metadata$Modified_Sample_Names <- sample_data$Modified_Sample_Names[match(compositional_data_with_metadata$Sample, sample_data$Sample)]

# Determine the taxonomic levels available
available_tax_levels <- c("domain", "phylum", "class", "order", "family", "genus")
default_tax_level <- "genus"

# UI
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  checkboxGroupInput("soil_groups", "Select Soil Groups:",
                     choices = unique(compositional_data_with_metadata$Soil), selected = character(0)),
  checkboxGroupInput("media_groups", "Select Media Groups:",
                     choices = unique(compositional_data_with_metadata$Media.type), selected = character(0)),
  checkboxGroupInput("engraft_groups", "Select Engraft Groups:",
                     choices = unique(compositional_data_with_metadata$Engraft), selected = character(0)),
  checkboxGroupInput("replicate_groups", "Select Well Plate Groups:",
                     choices = unique(compositional_data_with_metadata$Well_plate), selected = character(0)),
  selectInput("tax_level", "Select Taxonomic Level:",
              choices = c("domain", "phylum", "class", "order", "family", "genus"), selected = "genus"),
  plotlyOutput("composition")
)

# Update the server function to dynamically adjust the color based on the selected taxonomic level
server <- function(input, output, session) {
  filtered_data <- reactive({
    selected_soils <- input$soil_groups
    selected_media <- input$media_groups
    selected_engraft <- input$engraft_groups
    selected_replicate <- input$replicate_groups
    
    compositional_data_with_metadata[
      compositional_data_with_metadata$Soil %in% selected_soils &
        compositional_data_with_metadata$Media.type %in% selected_media &
        compositional_data_with_metadata$Engraft %in% selected_engraft & 
        compositional_data_with_metadata$Well_plate %in% selected_replicate,]
  })
  
  output$composition <- renderPlotly({
    filtered_data_combined <- filtered_data()
    if (is.null(filtered_data_combined)) {
      return(NULL)
    }
    
    # Ensure Sample is factor
    filtered_data_combined$Sample <- factor(ifelse(is.na(filtered_data_combined$Sample), "Unknown", filtered_data_combined$Sample))
    
    # Check if Abundance is numeric
    if (!is.numeric(filtered_data_combined$Abundance)) {
      return(NULL)
    }
    
    # Ensure Abundance is positive
    filtered_data_combined <- filtered_data_combined[filtered_data_combined$Abundance > 0, ]
    
    # Get unique values of the selected taxonomic level
    unique_tax_values <- unique(filtered_data_combined[[input$tax_level]])
    
    # Generate color palette based on the number of unique values
    tax_palette <- getPalette(length(unique_tax_values))
    names(tax_palette) <- unique_tax_values
    
    # Plotly bar chart
    p <- plot_ly(filtered_data_combined, x = ~Modified_Sample_Names, y = ~Abundance, type = "bar", color = as.factor(filtered_data_combined[[input$tax_level]]),
                 colors = tax_palette) %>%
      layout(title = "Title", barmode = "stack") %>%
      layout(legend = list(orientation = "h", x = 1, y = 1, font = list(size = 10), trace = list(y=1.5), tracegroupgap = 0, xanchor = "left", bgcolor = "rgba(255, 255, 255, 0.5)", bordercolor = "black")) %>%
      layout(
        margin = list(t = 30),
        xaxis = list(
          linecolor = "black",
          tickfont = list(size = 10),
          titlefont = list(color = "black", size = 10),
          tickangle = -90
        ),
        yaxis = list(
          linecolor = "black",
          tickfont = list(size = 10),
          titlefont = list(color = "black", size = 10)
        ),
        height = 600,
        width = 1272# Adjust the height as needed
      )
    
    return(p)
  })
}

# Run the app
shinyApp(ui, server)