library(shiny)

# Gray-Scottモデルのシミュレーション関数
gray_scott <- function(u, v, Du, Dv, F, k, steps) {
  n <- nrow(u)
  m <- ncol(u)
  for (i in 1:steps) {
    lap_u <- (cbind(u[,2:m], u[,1]) + cbind(u[,m], u[,1:(m-1)]) + 
                rbind(u[2:n,], u[1,]) + rbind(u[n,], u[1:(n-1),]) - 4 * u)
    
    lap_v <- (cbind(v[,2:m], v[,1]) + cbind(v[,m], v[,1:(m-1)]) + 
                rbind(v[2:n,], v[1,]) + rbind(v[n,], v[1:(n-1),]) - 4 * v)
    
    uvv <- u * v * v
    u <- u + Du * lap_u - uvv + F * (1 - u)
    v <- v + Dv * lap_v + uvv - (F + k) * v
  }
  return(list(u = u, v = v))
}

# UI
ui <- fluidPage(
  titlePanel("Gray-Scott Model Simulation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("Du", "Diffusion Rate of U (Du)", min = 0.01, max = 1, value = 0.16, step = 0.01),
      sliderInput("Dv", "Diffusion Rate of V (Dv)", min = 0.01, max = 1, value = 0.08, step = 0.01),
      sliderInput("F", "Feed Rate (F)", min = 0.01, max = 0.1, value = 0.035, step = 0.001),
      sliderInput("k", "Kill Rate (k)", min = 0.01, max = 0.1, value = 0.065, step = 0.001),
      numericInput("steps", "Number of Steps", value = 1000, min = 1)
    ),
    mainPanel(
      plotOutput("plot_u"),
      plotOutput("plot_v")
    )
  )
)

# Server
server <- function(input, output) {
  simulate <- reactive({
    n <- 100
    u <- matrix(1, n, n)
    v <- matrix(0, n, n)
    u[45:55, 45:55] <- 0.50
    v[45:55, 45:55] <- 0.25
    gray_scott(u, v, input$Du, input$Dv, input$F, input$k, input$steps)
  })
  
  output$plot_u <- renderPlot({
    res <- simulate()
    image(res$u, col = terrain.colors(256), main = "Concentration of U")
  })
  
  output$plot_v <- renderPlot({
    res <- simulate()
    image(res$v, col = terrain.colors(256), main = "Concentration of V")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
