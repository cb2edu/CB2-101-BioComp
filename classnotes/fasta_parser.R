#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = T)
filename<-args[1]
fas_file <- file(filename, "r")

pattern <- "^\\>sp\\|(\\S+)\\|"
id <- "Q6GZX3"
inside <- FALSE
buff <- c("")

while (length(line <- readLines(fas_file, n=1)) > 0) {
  m <- regexec(pattern, line, perl = T)
  if (m[[1]][1] != -1) { # id line not sequnce
    s <- regmatches(line, m)
    #cat(s[[1]][2], s[[1]][3], "\n")
    if (s[[1]][2] == id) { #id of our choice
      inside <- TRUE
      buff <- c(buff, line)
    }else {
      if (inside) {
        cat (buff, sep = "\n")
        inside <- FALSE
        break
      }else {
        
      }
      
    }
  }else { # Sequence line
    if (inside) {
      buff <- c(buff, line)
    }else {
      
    }
  }
}

close(fas_file)

if (inside) {
  cat (buff, sep = "\n")
}