leave_one_out <- function(model) {
  
  # Takes an rgasp model as an input and performs leave one out analysis.
  # Output is an object of the custom class 'leave_one_out', which can be
  # plotted and analysed.
  
  # extract information from the rgasp-class
  inputs <- data.frame(slot(model, 'input'))
  output <- slot(model, 'output')
  kernel <- slot(model, 'kernel_type')
  nugget.est <- slot(model, 'nugget.est')
  num_inputs <-  dim(inputs)[1]
  
  # create empty data frame to store results
  results <- data.frame(
    Actual = numeric(0),
    Predicted = numeric(0),
    Lower95 = numeric(0),
    Upper95 = numeric(0),
    Success = numeric(0),
    Input = list()
  )

  loo <- leave_one_out_rgasp(model)  

  # for each point in the inputs
  #for (i in 1:num_inputs) {
    
    #point <- inputs[i,]
    #actual <- output[i,]
    
    # remove it from the inputs and output
    #design <- inputs[-i,]
    #response <- output[-i,]
    
    # define new trends
    #trend <- as.matrix(cbind(1, design))
    #testing_trend <- as.matrix(cbind(1, point))
    
    # create an emulator with the remaining inputs/outputs
    #LOO_model <- rgasp(
    #  design=design, 
    #  response=response,
    #  kernel_type=kernel,
    #  nugget.est=nugget.est,
    #  trend=trend
    #)
    
    # ask the emulator to predict the original point
    #LOO_model.predict <- predict(LOO_model, point, testing_trend=testing_trend)
    #prediction <- LOO_model.predict$mean
    #lower95 <- LOO_model.predict$lower95
    #upper95 <- LOO_model.predict$upper95
    
    # if actual is within the emulator's uncertainty, call it a success
    if (lower95 > actual || upper95 < actual) {
      success <- FALSE
    } else {
      success <- TRUE
    }
    
    # add this point to the results
    new_line <- data.frame(
      Actual = actual,
      Predicted = prediction,
      Lower95 = lower95,
      Upper95 = upper95,
      Success = success,
      Input = list(point)
    )
    results <- rbind(results, new_line)
  }
  class(results) <- "leave_one_out"
  return(results)
}

  # TO DO: calculate RMSE and R^2
  # rmse <- sqrt(mean((results$Actual-results$Predicted)^2))
  # r_squared <- cor(results$Actual, results$Predicted)^2

plot.leave_one_out <- function(loo) {
  
  # Analyses the leave one out test by plotting the actual BISICLES (or whatever
  # numerical model) output against the emulator's predicted value, including
  # error bars. Highlights any cases where the actual value is outside of the
  # emulator's 95% confidence interval.
  
  # leave_one_out class can't be subsetted, so define an equivalent dataframe
  results <- data.frame(
    Actual = loo$Actual,
    Predicted = loo$Predicted,
    Lower95 = loo$Lower95,
    Upper95 = loo$Upper95,
    Success = loo$Success
  )
  
  # separate into points the emulator successfully predicted and points it didn't
  good <- subset(results, results$Success==TRUE)
  bad <- subset(results, results$Success==FALSE)
  
  # find bounds for axes
  max <- max(results$Upper95)
  min <- min(results$Lower95)

  #plot points and error bars for successes (in blue) and failures (in red)
  plot(
    good$Actual, 
    good$Predicted, 
    xlim=c(min, max), 
    ylim=c(min, max),
    xlab='Actual sea level contribution (m)',
    ylab='Predicted sea level contribution (m)',
    col = 'cornflowerblue',
    pch = 1,
    cex = 0.5
  )
  points(
    bad$Actual,
    bad$Predicted,
    col = 'indianred',
    pch = 1,
    cex = 0.5
  )
  segments(
    good$Actual,
    good$Lower95,
    good$Actual,
    good$Upper95,
    col = rgb(0, 0, 1, alpha = 0.2),
    lwd = 2
    )
  segments(
    bad$Actual,
    bad$Lower95,
    bad$Actual,
    bad$Upper95,
    col = rgb(1, 0, 0, alpha = 0.2),
    lwd = 2
    )
  lines(c(min-2, max+2), c(min-2, max+2)) # y=x
}

summary.leave_one_out <- function(loo) {
  
  fem <- loo$Predicted
  f   <- loo$Actual
  std <- (loo$Upper95 - loo$Predicted) / 2
  N   <- length(loo$Predicted)
  
  d   <- sqrt(sum((fem-f)^2/std))
  rmse <- sqrt(sum((fem-f)^2/N))
  
  cat('Normalized Euclidean Distance:', d, '\n')
  cat('RMSE:', rmse)
}
