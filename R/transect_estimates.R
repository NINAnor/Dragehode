library(tidyverse)
transektdata <- readRDS("C:/Users/matthew.grainger/Documents/Projects_in_development/Dragehode2023/data/transektdata.RDS")

lokalitetsdata <- readRDS("C:/Users/matthew.grainger/Documents/Projects_in_development/Dragehode2023/data/lokalitetsdata.RDS")

# Function to perform bootstrapping with mean imputation for missing values
perform_bootstrapping <- function(num_iterations, transect_density, transect_lengths) {
  total_population_samples <- numeric(num_iterations)
  
  for (i in 1:num_iterations) {
    # Resample with replacement
    resampled_density <- sample(transect_density, replace = TRUE)
    
    # Calculate total population for the resampled data
    total_population_samples[i] <- sum(resampled_density) * sum(transect_lengths) / length(transect_density)
  }
  
  return(total_population_samples)
}
unique_sites_years <- unique(transektdata[, c("Lokalitet", "År")])

# Initialize an empty data frame to store results
results_df <- data.frame()

# Loop over unique combinations of site and year
for (i in 1:nrow(unique_sites_years)) {
  # Extract current site and year
  current_site <- unique_sites_years$Lokalitet[i]
  current_year <- unique_sites_years$År[i]
  
  # Subset data for the current site and year
  subset_data <- transektdata[transektdata$Lokalitet == current_site & transektdata$År == current_year, ]
  
  num_transects <- dim(subset_data)[1]
  transect_lengths <- subset_data$Lengde
  plant_counts <- subset_data$Dragehode
  transect_density <- plant_counts / transect_lengths
  
  # Number of bootstrap iterations
  num_iterations <- 1000 
  
  # Perform bootstrapping with mean imputation
  bootstrap_samples <- perform_bootstrapping(num_iterations, transect_density, transect_lengths)
  
  # Remove NAs from bootstrap samples (if any)
  bootstrap_samples <- bootstrap_samples[!is.na(bootstrap_samples)]
  
  # Calculate mean and confidence interval
  mean_estimate <- mean(bootstrap_samples)
  mean_estimate <- mean_estimate * mean(transect_lengths)
  confidence_interval <- quantile(bootstrap_samples, c(0.025, 0.975)) * mean(transect_lengths)
  
  # Create a data frame for the current site and year
  current_results <- data.frame(
    mean_estimate = mean_estimate,
    confidence_interval_low = confidence_interval[1],
    confidence_interval_high = confidence_interval[2],
    lokalitet = current_site,
    Year = current_year
  )
  
  # Append the current results to the overall results data frame
  results_df <- rbind(results_df, current_results)
}

# Print or further analyze the results data frame
print(results_df)
#results_df |> view()

# Loop over unique sites
for (site in unique(lokalitetsdata$Lokalitet)) {
  # Filter data for the current site
  site_data <- results_df %>% filter(lokalitet == site)
  
  # Create the plot
  p <- ggplot(site_data) +
    geom_pointrange(aes(Year, mean_estimate, ymin = confidence_interval_low, ymax = confidence_interval_high), colour = "darkgreen") +
    theme_classic() +
    labs(x = "Survey year", y = "Mean estimate from transect data") +
    ggtitle(paste0(site))
  
  # Save the plot to a PDF file in the 'plots' folder
  ggsave(paste0("Figurer/trans_Est/", site, "_plot.png"), plot = p, device = "png"
         , width = 10, height = 6
         )
}

