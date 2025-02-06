mrf_samtools_bam2fq <- function(bam = NULL,
                                output_name = NULL,
                                threads = 1,
                                region = NULL,
                                outpath = NULL) {
  
  # Check if BAM file is provided
  if (is.null(bam)) {
    stop("The 'bam' argument must be specified.")
  }
  
  # Check if the BAM file exists
  if (!file.exists(bam)) {
    stop("The specified BAM file does not exist.")
  }
  
  # Determine output file name
  if (is.null(output_name)) {
    # Generate a default output filename based on the BAM file name
    output_name <- paste0(tools::file_path_sans_ext(bam), "_reads.fastq")
  } else {
    # Ensure the output file name has the correct extension
    if (!grepl("\\.fastq$", output_name)) {
      output_name <- paste0(output_name, ".fastq")
    }
  }
  
  # Ensure the outpath exists
  if (!is.null(outpath) && !dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  
  # Construct the output path
  output_name <- file.path(outpath, output_name)
  
  # Create a temporary BAM file name
  tmp_bam <- "tmp_bam.bam"
  
  # Construct the SAMtools command
  samtools_cmd <- paste0("samtools view -h -@", threads, " ", bam)
  
  # Add region specification if provided
  if (!is.null(region)) {
    samtools_cmd <- paste0(samtools_cmd, " ", region)
  }
  
  # Write the output of samtools view to the temporary BAM file
  samtools_cmd <- paste0(samtools_cmd, " > ", tmp_bam)
  
  # Execute the SAMtools command
  system(samtools_cmd, intern = FALSE)
  
  # Construct the bam2fq command
  bam2fq_cmd <- paste0("samtools bam2fq -@", threads, " ", tmp_bam, " > ", output_name)
  
  # Execute the bam2fq command
  system(bam2fq_cmd, intern = FALSE)
  
  # Remove the temporary BAM file
  file.remove(tmp_bam)
  
  # Inform the user
  cat(paste0("BAM to FASTQ conversion complete. Output written to ", output_name, "\n"))
  
  # Return the output file path
  return(output_name)
}


mrf_ont_pecat_cfgfile <- function(project = NULL,
                                  bed = NULL,
                                  reads = NULL,
                                  bam = NULL,
                                  genome_size = NULL,
                                  threads = 4,
                                  outpath = NULL,
                                  mode = c("Single", "Multiple"),
                                  human = c(TRUE, FALSE),
                                  medaka = c("TRUE", "FALSE"),  # Change to character vector
                                  medaka_model = c("r103_min_high_g345", "r103_min_high_g360", 
                                                   "r103_prom_high_g360", "r103_prom_snp_g3210", 
                                                   "r103_prom_variant_g3210", "r10_min_high_g303",
                                                   "r10_min_high_g340", "r941_min_fast_g303", 
                                                   "r941_min_high_g303", "r941_min_high_g330",
                                                   "r941_min_high_g340_rle", "r941_min_high_g344", 
                                                   "r941_min_high_g351", "r941_min_high_g360", 
                                                   "r941_prom_fast_g303", "r941_prom_high_g303", 
                                                   "r941_prom_high_g330", "r941_prom_high_g344", 
                                                   "r941_prom_high_g360", "r941_prom_high_g4011", 
                                                   "r941_prom_snp_g303", "r941_prom_snp_g322", 
                                                   "r941_prom_snp_g360", "r941_prom_variant_g303", 
                                                   "r941_prom_variant_g322", "r941_prom_variant_g360"),
                                  phase_clair = c(0, 1, 2)) {
  
  # Validate and set mode
  mode <- match.arg(mode)
  medaka_model <- match.arg(medaka_model)
  
  # Check if in a Conda environment
  if (is.na(Sys.getenv("CONDA_PREFIX"))) {
    stop("You are not in a CONDA ENV")
  }
  
  # Create a default config file
  if (system("pecat.pl config cfgfile", intern = FALSE) != 0) {
    stop("Failed to create the default config file.")
  }
  
  if (mode == "Single") {
    # Validate inputs for Single mode
    if (is.null(project)) stop("In single mode, 'project' cannot be NULL.")
    if (is.null(reads)) stop("In single mode, 'reads' must be specified.")
    if (is.null(genome_size)) stop("In single mode, 'genome_size' must be specified.")
    
    # Read and update the config file
    cfgfile <- readr::read_lines("cfgfile")
    cfgfile[1] <- paste0("project=", project)
    cfgfile[2] <- paste0("reads=", reads)  # Ensure this path is inside outpath
    cfgfile[3] <- paste0("genome_size=", genome_size)
    cfgfile[4] <- paste0("threads=", threads)
    cfgfile[11] <- "corr_correct_options=--score=weight:lc=10 --aligner edlib --filter1 oh=1000:ohr=0.01 --candidate n=600:f=30"
    
    # Medaka logic
    if (medaka == "TRUE") {  # Check medaka as a character string
      cfgfile[30] <- 'polish_map_options = -x map-ont -w10 -k19 -I 10g'
      cfgfile[33] <- "polish_medaka=1"
      cfgfile[34] <- "polish_medaka_command=singularity exec --containall -B `pwd -P`:`pwd -P` /path/to/medaka_v1.7.2.sif medaka"
      cfgfile[36] <- paste0("polish_medaka_cns_options=--model ", medaka_model)
    }
    
    # Human logic
    if (human) {
      cfgfile[13] <- "corr_rd2rd_options=-x ava-ont -f 0.005"
      cfgfile[16] <- "align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500 -f 0.005"
    }
    
    # Phase Clair logic
    if (phase_clair == 2) {
      cfgfile[19] <- "phase_method=2"
      cfgfile[22] <- "phase_phase_options= --coverage lc=20 --phase_options icr=0.1:icc=8:sc=10"
      cfgfile[23] <- "phase_clair3_command=singularity exec --containall -B `pwd -P`:`pwd -P` /path/to/clair3_v0.1-r12.sif /opt/bin/run_clair3.sh"
      cfgfile[29] <- 'asm2_assemble_options=--reducer0 "best:cmp=2,0.1,0.1|phase:sc=3" --contig_format dual,prialt --min_identity 0.98'
    }
    cfgfile <- append(cfgfile, "phase_filter_options = --threshold 1000")
    cfgfile <- append(cfgfile, "polish_use_reads=0")
    readr::write_lines(x = cfgfile, file = "cfgfile")
    
    # Warning about BED file in Single mode
    if (!is.null(bed) && !file.exists(bed)) {
      warning("In single mode, the BED file argument is ignored.")
    }
    
  } else if (mode == "Multiple") {
    # Validate inputs for Multiple mode
    if (is.null(bed) || !file.exists(bed)) stop("For multiple mode, a valid BED file must be specified.")
    if (is.null(bam) || !file.exists(bam)) stop("In multiple mode, 'bam' file must be specified.")
    
    # Warning about project argument in Multiple mode
    if (!is.null(project)) {
      warning("In multiple mode, the 'project' argument is ignored.")
    }
    
    # Read BED file
    bed_file <- readr::read_tsv(file = bed, col_names = FALSE)
    
    for (i in seq_len(nrow(bed_file))) {
      # Process each row in BED file
      cfgfile <- readr::read_lines("cfgfile")
      
      # Generate a region string
      region <- paste0(bed_file[i, 1], ":", bed_file[i, 2], "-", bed_file[i, 3])
      
      # Provide feedback
      cat(paste0("Extracting the reads for ", crayon::blue(bed_file[i, 4]), "\n"))
      # Call mrf_samtools_bam2fq function
      reads_fastq <- mrf_samtools_bam2fq(bam = bam, output_name = paste0(bed_file[i, 4], "_reads.fastq"), threads = threads, region = region, outpath = outpath)
      
      # Provide feedback
      cat(paste0("Creating the config file for ", crayon::red(bed_file[i, 4]), "\n"))
      
      # Update the config file for the current BED entry
      cfgfile[1] <- paste0("project=", bed_file[i, 4])
      cfgfile[2] <- paste0("reads=", file.path(outpath, paste0(bed_file[i, 4], "_reads.fastq")))  # Ensure the path points to outpath
      cfgfile[3] <- paste0("genome_size=", bed_file[i, 3] - bed_file[i, 2] + 1)
      cfgfile[4] <- paste0("threads=", threads)
      cfgfile[11] <- "corr_correct_options=--score=weight:lc=10 --aligner edlib --filter1 oh=1000:ohr=0.01 --candidate n=600:f=30"
      
      # Medaka logic
      if (medaka == 1) {  # Check medaka as a character string
        cfgfile[30] <- 'polish_map_options = -x map-ont -w10 -k19 -I 10g'
        cfgfile <- append(cfgfile, "polish_use_reads=0")
        cfgfile[33] <- "polish_medaka=1"
        cfgfile[34] <- "polish_medaka_command=singularity exec --containall -B `pwd -P`:`pwd -P` /path/to/medaka_v1.7.2.sif medaka"
        cfgfile[36] <- paste0("polish_medaka_cns_options=--model ", medaka_model)
      }
      
      # Human logic
      if (human) {
        cfgfile[13] <- "corr_rd2rd_options=-x ava-ont -f 0.005"
        cfgfile[16] <- "align_rd2rd_options=-X -g3000 -w30 -k19 -m100 -r500 -f 0.005"
      }
      
      # Phase Clair logic
      if (phase_clair == 2) {
        cfgfile[19] <- "phase_method=2"
        cfgfile[23] <- "phase_clair3_command=singularity exec --containall -B `pwd -P`:`pwd -P` /path/to/clair3_v0.1-r12.sif /opt/bin/run_clair3.sh"
        cfgfile[22] <- "phase_phase_options= --coverage lc=20 --phase_options icr=0.1:icc=8:sc=10"
        cfgfile[29] <- 'asm2_assemble_options=--reducer0 "best:cmp=2,0.1,0.1|phase:sc=3" --contig_format dual,prialt --min_identity 0.98'
        cfgfile <- append(cfgfile, "phase_filter_options = --threshold 1000")
      }
      
      # Write the updated config file
      config_filename <- file.path(outpath, paste0(bed_file[i, 4], "_cfgfile"))
      readr::write_lines(x = cfgfile, file = config_filename)
    }
    
  } else {
    stop("Invalid mode specified. Choose either 'Single' or 'Multiple'.")
  }
}

path <- "/path/to/your/data/sample.bam"
bam_list <- list.files(path = path,
                       pattern = ".bam$",
                       full.names = T)


BED_path <- "/path/to/your/Genes.bed"

for (i in seq_along(bam_list)) {
  
  amostra <- gsub(pattern = paste0(path,"/"),replacement = "",bam_list[i])
  amostra <- gsub(pattern = ".bam", replacement = "",amostra)
  
  # Create the AMOSTRAS directory
  amostra_dir <- paste0(path, "/PECAT/", amostra)
  dir.create(amostra_dir)
  
  # Run mrf_ont_pecat_cfgfile with the AMOSTRAS directory as outpath
  mrf_ont_pecat_cfgfile(project = amostra,
                        bed = BED_path,
                        bam = bam_list[i],
                        threads = 14,
                        mode = "Multiple",
                        outpath = paste0(amostra_dir, "/"),
                        human = T,
                        phase_clair = 2,
                        medaka=1,
                        medaka_model="r103_prom_high_g360")
}
