# Load necessary libraries
library(furrr)
library(fs)
library(progress)

# Define a function to copy a single file
copy_file <- function(file, source_dir, destination_dir) {
  tryCatch({
    # Determine the relative path of the file
    relative_path <- path_rel(file, start = source_dir)
    # Construct the destination path
    dest <- path(destination_dir, relative_path)
    # Create destination directory if it doesn't exist
    dir_create(path_dir(dest))
    # Copy the file
    file_copy(file, dest, overwrite = TRUE)
  }, error = function(e) {
    message(paste("Error copying", file, ":", e$message))
  })
}

# Main function to copy files in parallel
copy_files_in_parallel <- function(source_dir, destination_dir, n_cores = parallel::detectCores() - 1) {
  # Get a list of all files (including in subdirectories) in the source directory
  files <- dir_ls(source_dir, recurse = TRUE, type = "file")
  
  # Create the destination directory if it doesn't exist
  dir_create(destination_dir)
  
  # Initialize a progress bar
  pb <- progress_bar$new(
    format = "Copying files [:bar] :percent :elapsed / :eta",
    total = length(files),
    clear = FALSE,
    width = 60
  )
  
  # Wrapper to update progress bar
  copy_with_progress <- function(file) {
    copy_file(file, source_dir, destination_dir)
    pb$tick() # Update progress bar
  }
  
  # Setup future plan for parallelism
  plan(multisession, workers = n_cores)
  
  # Copy files in parallel using furrr
  future_walk(files, copy_with_progress)
}

# Example usage
source_directory <- "data/glider/pomBODCREQ-5915/unit_345/pyglider_binaries"
destination_directory <- "data/temporary"

# Adjust the number of cores if needed
num_cores <- 4

copy_files_in_parallel(source_directory, destination_directory, num_cores)
