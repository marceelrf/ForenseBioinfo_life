#!/bin/bash

# Path to the parent folder containing subfolders
parent_folder="path/to/fastq/"

# Path to the file you want to copy
file_to_copy="/path/to/run_PECAT.R"

# The command you want to run inside each subfolder
command_to_run="Rscript run_PECAT.R"

log_file="log.txt"

# Function to write to the log
log_message() {
  echo "$(date +'%Y-%m-%d %H:%M:%S') - $1\n" >> "$log_file"
}

# Loop through each subfolder inside the parent folder
for subfolder in "$parent_folder"/*/; do
  if [ -d "$subfolder" ]; then
    # Copy the file into the subfolder
    cp "$file_to_copy" "$subfolder"

    # Change to the subfolder
    cd "$subfolder" || { 
      log_message "Failed to enter $subfolder"; 
      continue; 
    }

    # Ensure the command runs successfully
    if $command_to_run; then
      log_message "Successfully executed command in: $subfolder"
    else
      log_message "Failed to execute command in: $subfolder"
    fi

    # Return to the parent folder after processing the subfolder
    cd "$parent_folder" || { 
      log_message "Failed to return to parent folder after processing $subfolder"; 
      continue; 
    }
  fi
done
