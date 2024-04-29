library(RobustGaSP)

find_rectangle_dimensions <- function(N) {
  a <- floor(sqrt(N))
  b <- ceiling(N/a)
  return ( c(a, b) )
}

colors <- c('orangered', 'seagreen3', 'royalblue1', 'plum2', 'slateblue1', 'sienna1', 'goldenrod1', 'darkseagreen', 'orchid')

main_effects <- function(model, inputs) {
  
  N <- length(inputs)
  num_testpoints <- 1000
  
  par(mfrow=find_rectangle_dimensions(N))
  for (i in 1:N) {
    
    param_testpoints <- as.matrix(seq(-sqrt(3), sqrt(3), length.out=num_testpoints))
    testing_inputs <- matrix(0, ncol=N, nrow=num_testpoints)
    testing_inputs[, i] <- param_testpoints

    testing_trend <- as.matrix(cbind(1, testing_inputs))
    
    model.predict <- predict(
      model,
      testing_inputs,
      testing_trend=testing_trend
    )
    
    ymax <- max(model.predict$upper95)
    ymin <- min(model.predict$lower95)
    name <- colnames(inputs)[i]
    
    # empty plot
    plot(
      1, 
      1,
      type='l',
      xlim=c(-sqrt(3), sqrt(3)),
      ylim=c(-15, 25), #c(ymin, ymax),
      xlab=name,
      ylab='sea level contribution (m)'
    )
    
    # shade 5-95% uncertainty region
    polygon(
      c(param_testpoints, rev(param_testpoints)), 
      c(model.predict$lower95, rev(model.predict$upper95)), 
      col=adjustcolor(colors[i], alpha.f=0.3),
      border=F
    )
    
    # plot emulator mean
    lines(param_testpoints, model.predict$mean, type='l', col=colors[i])
  }
}
  
  