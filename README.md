## INFOF434 

### TP3  

1. Download the hg19 RefSeq genes annotation file from UCSC genome.ucsc.edu / Tools / Table Browser / refGene Table
2. Load the table in R
3. Filter the coding genes (See the RefSeq nomenclature)
4. Generate a table of the transcripts length (from TSS to TES)
5. Plot an histogram of gene length distribution
6. Select all the long non coding transcript (length ≥ 2000bp)
7. Plot an histogram of long non coding gene length distribution
8. Compare the coding / long non coding transcript distribution

### TP4 

1. Download the chr21 and chr22 (hg19) from genome.ucsc.edu Downloads / Genome Data / Data set by chromosome
2. Use the bioconductor package “Biostrings" (readDNAStringSet) to read the fasta file(s)
3. Reload the hg19 gene annotation file (TP2)
4. Use the stringr library to count the number of A, C, G and T
5. Compute the G+C percentage for each chromosome
6. Select the regions corresponding to the first condon (cdsStart)
7. Plot the first codon distribution
