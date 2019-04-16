library(stringr)
library(purrr)
library(stringdist)
library(Biostrings)

countNucleotides <- function(dnaStr) {
  nucleobases <- c("A", "C", "G", "T")
  return(nucleobases %>% 
           map(function(x) str_count("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC", x))
  )
}

dnaToRna <- function(dnaStr) {
  return(gsub("T", "U", dnaStr))
}

reverseComplement <- function(dnaStr) {
  createComplement <- function(dnaStr) {
    return(chartr("ACGT", "TGCA", dnaStr))
  }
  
  reverseString <- function(str) {
    return(paste(rev(unlist(strsplit(str, split=""))), collapse=""))
  }
  
  return(reverseString(createComplement(dnaStr)))
}

gcContent <- function(dnaStr) {
  return(str_count(dnaStr, "G|C") / str_length(dnaStr))
}

countPointMutations <- function(dnaStr1, dnaStr2) {
  return(stringdist(dnaStr1, dnaStr2, method = "hamming"))
}

rnaToProtein <- function(rnaStr) {
  return(translate(RNAString(rnaStr)))
}