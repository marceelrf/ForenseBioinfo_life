mrf_ont_pecat_run <- function(){


  library(crayon)
  cfgfile_list <-
    list.files(pattern = "_cfgfile$")

  hg_theme <- red$bold$underline
  number_theme <- silver$bold$underline

  gene <- gsub(pattern = "_cfgfile$",replacement = "",cfgfile_list)

  for(i in seq_along(cfgfile_list)){



    cat(bold("----------------------------------------\n\n"))
    cat(paste0("Running PECAT for ", hg_theme(gene[i])), "\n\n")


    pecat_cmd <- paste0("pecat.pl unzip ",cfgfile_list[i])

    system(pecat_cmd, intern = FALSE)

    cat(paste0("\n",number_theme(i)," of ", number_theme(length(gene))," done!"), "\n\n")

    cat(bold("----------------------------------------\n\n"))

  }
}



mrf_ont_pecat_run()
