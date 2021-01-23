library(readr)

sig <- readLines("./dunning_2018_bacterial.txt")

sig <- strsplit(sig, " ")[[1]]
sig <- sig[!startsWith(sig, "ILMN")]

sig <- sig[sapply(sig, nchar) > 1]

sig <- unique(sig)

sig <- gsub("\\*", "", sig)

write_lines(sig, "dunning_2018_bacterial_cleaned.txt")
