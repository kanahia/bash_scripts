#!/usr/bin/env Rscript

# Author: kanahia
# Quality check: Macio
# Date: 21/08/20
description <-
  "Script extract sequences based on PCR primers from the transcriptome,
  generate fasta files which are later on mapped to the genome,
  converted to bed files and ready for visualization in genome browser.


  Usage:
  args1 = path to the template
  args2 = path to the transcriptome
  args3 = path to the output dir
  "

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  cat(description)
  stop()
}

# Packages
library("readxl")
library("dplyr")
library("stringr")

# Arguments
##path to the template file

# path <- 
#   "/home/jason/data/ET31_48_72/in_situ/insitu_48hpf/probes summary_48hpf.xlsx"
path <- args[1] 
   #"/home/jason/data/ET31_48_72/in_situ/insitu_72/probes_summary_72.xlsx"

path_to_transcriptome <- args[2] 
  #"/home/jason/data/references/Danio_rerio-transcriptomes/GRCz11-release-95/Danio_rerio.GRCz11.cdna.all.fa.gz"

output_path <- args[3]
  #"/home/jason/data/ET31_48_72/in_situ/"



# Main functions:

## Read template file
read_template <- function(path) {
  template <- readxl::read_xlsx(path) %>%
    dplyr::select(!c(1, 3, 4, 7, 9)) %>%
    dplyr::filter(!is.na(Sequence)) %>%
    dplyr::mutate(
      Forward = gsub("5' ", "", gsub(" 3'", "", Sequence)))
  
  template$Forward <-
    ifelse(nchar(template$Forward > 35),
           str_remove_all(
             template$Forward, 
             "TAA TAC GAC TCA CTA TAG GGA GA |TAATACGACTCACTATAGGGAGA|TAA TAC GAC TCA CTA TAG GGA GA"), 
           template$Forward) %>% 
    str_replace_all(" ", "")
  
  template <- template %>%
    dplyr::select(!3)
  
  rev_primer <- template$Forward[seq(2, nrow(template), by =2)] #template$Forward[seq(2,30, by =2)]
  
  template <- template %>%
    dplyr::filter(!is.na(Name)) %>%
    dplyr::mutate(Reverse = str_replace_all(rev_primer, " ", "")) %>%
    dplyr::select(!2) %>%
    dplyr::filter(!is.na(`Gene ID`)) %>%
    as.data.frame()
  
  return(template)
}

#output directory
output_dir <- paste0(output_path, "output_dir")

if(dir.exists(output_dir) == TRUE){
  print(paste0("directory already exists: ", output_dir))
}else {
  dir.create(paste0(output_path, "output_dir"))
}
  

template <- read_template(path)

gene <- template$Name
gene_ensembl_no <- template$`Gene ID`

Fwd <- template$Forward
Rev <- template$Reverse

## Get reverse complement
compReverse <- function(x) {
  tmp1 <- unlist(strsplit(as.character(x), NULL))
  paste0(rev(unlist(lapply(tmp1, function(x) {
    if (x == "A") {
      complementary <- "T"
    } else if (x == "T") {
      complementary <- "A"
    } else if (x == "C") {
      complementary <- "G"
    } else if (x == "G") {
      complementary <- "C"
    }
  }))), collapse = "")
}

## Read fasta file
ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs, stringsAsFactors = FALSE)
  # Return the data frame as a result object from the function
  return(DF)
}

## read transcriptome
transcriptome <-
  if (exists("t_GRCz11r95") == TRUE) {
    t_GRCz11r95
  } else {
    ReadFasta(path_to_transcriptome)
  }

## Extract gene details and generate fasta
extractoR <- function(gene_name, ensembl_id, Fwd, Rev, transcriptome) {
  Rev_in_seq <- compReverse(Rev)
  
  string_base <- 
    dplyr::filter(transcriptome, 
                  grepl(gene_name, .data$name) | grepl(paste0(ensembl_id, ".[0-9]"), .data$name))
  
  if (nrow(string_base) == 0) {
    out <- tibble::tibble(
      ensembl_id = ensembl_id,
      ensembl_transcript = "Not found in provided transcriptome",
      gene_name = gene_name,
      sequence = "Not found in provided transcriptome - possible lincRNA",
      length = "Not found in provided transcriptome",
      fasta = "Not found in provided transcriptome"
    )
    return(out)
  }
  
  string <- 
    dplyr::filter(
      transcriptome, 
      grepl(gene_name, .data$name) | grepl(paste0(ensembl_id, ".[0-9]"), .data$name))[, 2]
  
  
  gene_details <- string_base[, 1]
  
  #check within which transcript the primer bind and amplify the region used for probe generation
  gene_in_probe <-
    gene_details[!is.na(nchar(str_sub(string = string, 
                                      start = str_locate(string, Fwd)[, 1], 
                                      end = str_locate(string, Rev_in_seq)[, 2]))) == TRUE]
  
  gene_id <-
    gsub("gene:", "", gsub("gene_symbol:", "",
                           strsplit(gene_in_probe[1], " ")[[1]][c(1, 4, 7)]))
  
  gene_seq <-
    string_base %>% 
    dplyr::filter(grepl(gene_id[1], name)) %>%
    dplyr::select(sequence) %>%
    str_sub(.,
            str_locate(., Fwd)[, 1],
            str_locate(., Rev_in_seq)[, 2]) %>% na.omit()
  
  out <- tibble::tibble(
    ensembl_id = ensembl_id,   #gene_id[2],
    ensembl_transcript = gene_id[1],
    gene_name = gene_id[3],
    sequence = gene_seq,
    length = nchar(as.character(.data$sequence)),
    fasta = paste0(">", 
                   paste(.data$gene_name, .data$ensembl_id, .data$ensembl_transcript, sep = " "), 
                   "\n", 
                   .data$sequence)
  )
  
  return(out)
}

# Output data
out <- list()
for(i in 1:length(gene)) {
  out[[i]] <- extractoR(gene[i], gene_ensembl_no[i], Fwd[i], Rev[i], transcriptome)
}

df <- do.call(rbind, out)

# write output fasta files
for(i in 1:nrow(df)) {
  myfile <- paste0(output_dir, "/", df[i, "gene_name", drop = TRUE], ".fa")
  writeLines(df[i, "fasta", drop = TRUE], con = myfile)
}

