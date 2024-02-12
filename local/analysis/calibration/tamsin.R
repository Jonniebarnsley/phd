#########################

# Define a function to calibrate with one observation (e.g. a mean or trend)

do_calibration <- function(distances, calib_type) {
  
  if (calib_type == "history_matching") {
  
    # Create NROY index
    nroy_index <- which( abs(distances/obs_err) < 3 )
    print(sprintf("Selected %i (%.1f%%) of %i ensemble members",
                  length(nroy_index), length(nroy_index)*100/length(distances), length(distances)))
    
    return( nroy_index )
  }
  
  if (calib_type == "Bayesian") {
    
    # Calculate weights with Gaussian univariate likelihood
    weightsraw <- sapply( distances, function(x) exp( -0.5 * sum( x * x , na.rm=TRUE) / (obs_err * obs_err) ) )

    # Normalise weights so they sum to 1
    weights <- weightsraw / sum( weightsraw )
    
    return( weights )
  }
}

#########################

# Calculate distances between each ensemble member and the single observation
dist <- sims[ , “y2014” ] – obs

# Call function above to do Bayesian calibration
weights <- do_calibration( dist, "Bayesian" ) 

# Estimate prior and posterior densities – yes it really is this easy
prior <- density( sims[ , “y2100”] )
posterior <- density( sims[ , “y2100”], weights = weights )

# Plot both
plot( prior )
lines( posterior, col ='red')
       