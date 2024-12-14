# Libraries
library(shiny)
library(shinythemes) 
# Define the UI
ui <- fluidPage(
  tags$head( # Inserted JavaScript
 #   tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$script(HTML(
      "
      $(document).ready(function() {
        $('#reset_button').click(function() {
          location.reload();
        });
      });
      "
    ))
  ),
  titlePanel(
    h1("PenicillinX", align = "center")
  ), 
  theme = shinytheme("superhero"),
  sidebarLayout(
    sidebarPanel(
      # Dropdown menu for selecting the first antibiotic with placeholder
      selectInput(
        inputId = "antibiotic1",
        label = "Select Antibiotic 1:",
        choices = NULL,      # Populated by the csv
        selected = "",       # Default no selection
        selectize = TRUE     # Enables autocomplete
      ),
      
      # Dropdown menu for selecting the second antibiotic with placeholder
      selectInput(
        inputId = "antibiotic2",
        label = "Select Antibiotic 2:",
        choices = NULL,      # Populated by the csv
        selected = "",       # Default no selection
        selectize = TRUE     # Enables autocomplete
      ),
      
      # Dropdown menu for selecting the third antibiotic with placeholder
      selectInput(
        inputId = "antibiotic3",
        label = "Select Antibiotic 3:",
        choices = NULL,      # Will be populated by the csv 
        selected = "",       # Default no selection 
        selectize = TRUE     # Enables autocomplete
      ),
      actionButton("reset_button", "Reset", style = "color: #4e5d6c; background-color: #fff; border: 1px solid #4e5d6c; border-radius: 3px;"), 
    ),
    mainPanel(
      # Where the antibiotic pictures will be shown
      h3("Selected Antibiotics:"),
      uiOutput("antibiotic_images"),
      
      # Display the matching attributes 
      h3("Matching Features:"),
      verbatimTextOutput(outputId = "matches"),
      
      # If there are common features among 3 antibiotics, then this will show up 
      uiOutput("common_matches_output")
    )
  ),
 
 tags$footer( # Disclaimer and link to github
   p(align = "center",
     "App currently in development. This should absolutely not be used in any clinical setting.",
     a(href = "https://github.com/liamlah/PenicillinX", target = "_blank", 
       "More information can be found on the GitHub page")
   ),
   style = "
      position: fixed;
      bottom: 0; 
      left: 50%;
      transform: translateX(-50%); 
      width: 80%;
      height: 20px;
      color: black;
      padding: 0px;
      background-color: lightgrey;
      z-index: 100;
    "
 ) 
)

# Define the server
server <- function(input, output, session) {
  # Reads from the CSV ---- semi-dummy data at present
  antibiotics_data <- read.csv("antibiotics_data.csv", stringsAsFactors = FALSE)
  
  # Puts a placeholder with the antibiotics data
  antibiotic_choices <- c("Select antibiotic" = "", antibiotics_data$Antibiotic)
  
  # Dropdown with antibiotic choices - consider whether a user might want more than 3 or a button to add more dynamically 
  observe({
    updateSelectInput(session, "antibiotic1", choices = antibiotic_choices, selected = "")
    updateSelectInput(session, "antibiotic2", choices = antibiotic_choices, selected = "")
    updateSelectInput(session, "antibiotic3", choices = antibiotic_choices, selected = "")
  })
  
  # User can select 2 or more antibiotics
  selected_antibiotics <- reactive({
    selections <- c(input$antibiotic1, input$antibiotic2, input$antibiotic3)
    selections <- selections[selections != ""]  # Remove empty selections
    return(selections)
  })
  
  # Finds the common features among two antibiotics
  find_common <- function(antibiotic_a, antibiotic_b) {
    row_a <- antibiotics_data[antibiotics_data$Antibiotic == antibiotic_a, ]
    row_b <- antibiotics_data[antibiotics_data$Antibiotic == antibiotic_b, ]
    
    side_chain_match <- ifelse(row_a$SideChain == row_b$SideChain, row_a$SideChain, NA)
    core_ring_match <- ifelse(row_a$CoreRing == row_b$CoreRing, row_a$CoreRing, NA)
    
    side_chain <- if (!is.na(side_chain_match)) side_chain_match else NULL
    core_ring <- if (!is.na(core_ring_match)) core_ring_match else NULL
    
    return(list(side_chain = side_chain, core_ring = core_ring))
  }
  
  # Finds common features among the selected antibiotics
  find_common_all <- function(selected) {
    selected_data <- antibiotics_data[antibiotics_data$Antibiotic %in% selected, ]
    
    # First antibiotic features
    common_side_chain <- selected_data$SideChain[1]
    common_core_ring <- selected_data$CoreRing[1]
    
    # find common features using intersect
    for (i in 2:nrow(selected_data)) {
      common_side_chain <- intersect(common_side_chain, selected_data$SideChain[i])
      common_core_ring <- intersect(common_core_ring, selected_data$CoreRing[i])
    }
    
    # clean up method - might not be necessary after dummy data taken out, but might be depending on data
    common_side_chain <- common_side_chain[common_side_chain != "" & !is.na(common_side_chain)]
    common_core_ring <- common_core_ring[common_core_ring != "" & !is.na(common_core_ring)]
    
    return(list(side_chain = common_side_chain, core_ring = common_core_ring))
  }
  
  # Formats the feature matches
  format_matches <- function(features) {
    output_strings <- c()
    
    if (!is.null(features$side_chain) && length(features$side_chain) > 0) {
      output_strings <- c(output_strings, paste("Side Chain:", paste(features$side_chain, collapse = ", ")))
    }
    
    if (!is.null(features$core_ring) && length(features$core_ring) > 0) {
      output_strings <- c(output_strings, paste("Core Ring:", paste(features$core_ring, collapse = ", ")))
    }
    
    if (length(output_strings) > 0) {
      return(paste(output_strings, collapse = " | "))
    } else {
      return("No matching features.")
    }
  }
  
  # Text for the output of matching features
  output$matches <- renderText({
    selected <- selected_antibiotics()
    num_selected <- length(selected)
    
    if (num_selected < 2) {
      return("Please select at least two antibiotics for comparison.")
    }
    
    # Ensures someone hasnt picked the same antibiotic twice
    if (length(unique(selected)) < num_selected) {
      return("Please select distinct antibiotics for comparison.")
    }
    
    output_text <- ""
    
    # Compares in pairs 
    pairs <- combn(selected, 2, simplify = FALSE)
    for (pair in pairs) {
      features <- find_common(pair[1], pair[2])
      match_info <- format_matches(features)
      output_text <- paste0(output_text, 
                            "🔹 Matching features between ", 
                            pair[1], " and ", pair[2], ":\n", 
                            match_info, "\n\n")
    }
    
    # output of the result
    return(output_text)
  })
  
  # Common features of all three. in separate box for readability and distinctness
  output$common_matches_output <- renderUI({
    selected <- selected_antibiotics()
    num_selected <- length(selected)
    
    # Ensures three antibiotics are selected before doing this for aesthetics
    if (num_selected == 3) {
      common_features <- find_common_all(selected)
      common_features_formatted <- format_matches(common_features)
      
      # Prepares the output for the common features
      if ((!is.null(common_features$side_chain) && length(common_features$side_chain) > 0) ||
          (!is.null(common_features$core_ring) && length(common_features$core_ring) > 0)) {
        return(
          div(style = "margin-top: 20px; padding: 10px; border: 1px solid #ccc; background-color: #f0f8ff; color: #000;",
              strong(paste("🔹 Matching features among all three antibiotics (", 
                           paste(selected, collapse = ", "), "):")),
              common_features_formatted
          )
        )
      } else {
        return(
          div(style = "margin-top: 20px; padding: 10px; border: 1px solid #ccc; background-color: #f0f8ff; color: #000;",
              strong("❗ No matching features among all three antibiotics."))
        )
      }
    } else {
      return(NULL)  # No output if not all three antibiotics are selected 
    }
  })
  
  # Images and labels for the antibiotics
  output$antibiotic_images <- renderUI({
    selected <- selected_antibiotics()
    
    if (length(selected) == 0) {
      return(tags$p("No antibiotics selected."))
    }
    
    # Create a list to hold image-tag pairs
    image_tags <- lapply(selected, function(antibiotic) {
      img_src <- paste0("ImagesInv/", antibiotic, ".svg")
      img_path <- file.path("www", "Images", paste0(antibiotic, ".svg"))
      if (!file.exists(img_path)) {
        img_src <- "Images/placeholder.png" # Placeholder image while I haven't rendered all of the molecules yet 
      }
      
      div(
        style = "text-align: center; margin: 10px;",
        img(src = img_src, height = "150px", alt = paste("Image of", antibiotic)),
        br(),
        strong(antibiotic)
      )
    })
    
    # Arrange images in a fluidRow with columns
    fluidRow(
      lapply(image_tags, function(tag) {
        column(width = 4, tag)
      })
    )
  })
}

# Run the app
shinyApp(ui, server)


## TO-DO
# 2. Replace dummy data with researched data
# 3. See if it is possible to highlight parts of the images with the corresponding matching features
# 4. Add a button to add more antibiotics dynamically
# 6. improve the look of the app
# 7. Fix issue where images overlap at certain screen sizes

## DONE 
# 1. Add the rest of the missing antibiotic images
# 5. Add a button to reset the selection
