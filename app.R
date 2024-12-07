library(shiny)

# Interface
ui <- fluidPage(
  titlePanel("Antibiotic Cross Reactivity App"),
  
  sidebarLayout(
    sidebarPanel(
      # Dropdown forfirst antibiotic
      selectInput(
        inputId = "antibiotic1",
        label = "First historical reaction:",
        choices = NULL, # Will be populated dynamically
        selected = '',
        selectize = TRUE # Lets you narrow down with a string
      ),
      
      # Dropdown for second antibiotic
      selectInput(
        inputId = "antibiotic2",
        label = "Second historical reaction:",
        choices = NULL, # Will be populated dynamically
        selected = NULL,
        selectize = TRUE # Lets you narrow down with a string
      )
    ),
    mainPanel(
      # Display the matching attributes (Side Chains & Core Rings)
      h3("Matching Features:"),
      verbatimTextOutput(outputId = "matches")
    )
  )
)

# Server
server <- function(input, output, session) {
  # CSV load, consider native dataframe if issues when hosting
  antibiotics_data <- read.csv("antibiotics_data.csv", stringsAsFactors = FALSE)
  
  # Pulls antibiotic names for the dropbdowns
  observe({
    updateSelectInput(session, "antibiotic1", choices = antibiotics_data$Antibiotic)
    updateSelectInput(session, "antibiotic2", choices = antibiotics_data$Antibiotic)
  })
  
  # Compare selectons of antibiotics and find matches
  output$matches <- renderText({
    # Get the selected antibiotics
    anti1 <- input$antibiotic1
    anti2 <- input$antibiotic2
    
    # stops anything ahppenign if nothing selectedd
    if (is.null(anti1) || is.null(anti2) || anti1 == "" || anti2 == "") {
      return("Please select two antibiotics for comparison.")
    }
    
    # searches the rows for the selected abx
    row1 <- antibiotics_data[antibiotics_data$Antibiotic == anti1, ]
    row2 <- antibiotics_data[antibiotics_data$Antibiotic == anti2, ]
    
    # Compares side chains + core rings
    side_chain_match <- ifelse(row1$SideChain == row2$SideChain, row1$SideChain, NA)
    core_ring_match <- ifelse(row1$CoreRing == row2$CoreRing, row1$CoreRing, NA)
    
    # creates OP text
    output_text <- ""
    if (!is.na(side_chain_match) && side_chain_match != "") {
      output_text <- paste0(output_text, "Matching Side Chain: ", side_chain_match, "\n")
    }
    if (!is.na(core_ring_match) && core_ring_match != "") {
      output_text <- paste0(output_text, "Matching Core Ring: ", core_ring_match, "\n")
    }
    
    # text for no matches
    if (output_text == "") {
      output_text <- "No matches found for Side Chain or Core Ring."
    }
    
    # the OP
    return(output_text)
  })
}

# Run
shinyApp(ui, server)