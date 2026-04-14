library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Gene Expression X Normal Distribution"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("numGenes", "Number of Genes:", value = 1, min = 1, max = 200),
      sliderInput("pfreq", "p frequence:", value = 0.5, min = 0, max = 1, step = 0.01)
      #numericInput("meanExpr", "Mean Expression Level:", value = 10, min = 0),
      #numericInput("sdExpr", "Standard Deviation of Expression:", value = 2, min = 0)
    ),
    
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
output$distPlot <- renderPlot({
# Generate gene expression data
n <- input$numGenes
p <- input$pfreq
q <- 1 - p
n_total_alleles <- 2 * n

if (n_total_alleles <= 200) {
  x_vals <- 0:n_total_alleles
} else {
  x_vals <- seq(0, n_total_alleles, length.out = 200)
}


y_vals <- dbinom(round(x_vals), size = n_total_alleles, prob = p)
df <- data.frame(Genotype = x_vals, Frequency = y_vals)

theoretical_mean <- n_total_alleles * p
theoretical_sd <- sqrt(n_total_alleles * p * q)

if(theoretical_sd == 0) {
  theoretical_sd <- 0.001
}

# Plot the distribution of gene expression levels
ggplot(df, aes(x = Genotype, y = Frequency)) +
  geom_col(fill = "steelblue", color = if(n < 50) "white" else NA, alpha = 0.6) +
  stat_function(
    fun = dnorm,
    args = list(mean = theoretical_mean, sd = theoretical_sd),
    color = "red",
    linewidth = 1.2) +
theme_minimal() +
labs (
    title = paste("Simulation for", n, "Loci | p =", p),
    subtitle = "The red line represents the theoretical normal distribution",
    x = "Number of Genotypes",
    y = "Probability Density") +
coord_cartesian(xlim = c(0, n_total_alleles))})
}


server1 <- function(input, output) {
output$distPlot <- renderPlot({
# Generate gene expression data
n <- input$numGenes
p <- input$pfreq
q <- 1 - p
n_total_alleles <- 2 * n

if (n_total_alleles <= 200) {
  x_vals <- 0:n_total_alleles
} else {
  x_vals <- seq(0, n_total_alleles, length.out = 200)
}


y_vals <- dbinom(round(x_vals), size = n_total_alleles, prob = p)
df <- data.frame(Genotype = x_vals, Frequency = y_vals)

theoretical_mean <- n_total_alleles * p
theoretical_sd <- sqrt(n_total_alleles * p * q)

if(theoretical_sd == 0) {
  theoretical_sd <- 0.001
}

if (n_total_alleles <= 50) {
  xlimi <- c(0, n_total_alleles)
} else {
  xlimi <- c(theoretical_mean - 5 * theoretical_sd, theoretical_mean + 5 * theoretical_sd)
}

# Plot the distribution of gene expression levels
ggplot(df, aes(x = Genotype, y = Frequency)) +
  geom_col(fill = "steelblue", color = if(n < 50) "white" else NA, alpha = 0.6) +
  stat_function(
    fun = dnorm,
    args = list(mean = theoretical_mean, sd = theoretical_sd),
    color = "red",
    linewidth = 1.2) +
theme_minimal() +
labs (
    title = paste("Simulation for", n, "Loci | p =", p),
    subtitle = "The red line represents the theoretical normal distribution",
    x = "Number of Genotypes",
    y = "Probability Density") +


coord_cartesian(xlim = xlimi)})
}


shinyApp(ui = ui, server = server1)

