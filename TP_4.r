library(stringr)
library('Biostrings')

setwd("D:\\workspace\\bioinformatics\\BIOLF449\\INFOF434")

fa <- readDNAStringSet('genome.fa')

chr21 <-fa$chr21
chr22 <- fa$chr22

t <- read.table('data', header = T, comment.char = '')[,-1]

codingAnnotationsByChromosomeAndStrand <- function(chromosomeName, strand){
  annotation <- t[t$chrom == chromosomeName, ]
  annotation <- annotation[annotation$strand == strand, ]
  annotation <- annotation[grep('NM', annotation$name),]
  annotation
}


#Returns start codons in chromosome chr using given annotations
startCodonsPlus<-function(chr, annotations){
  startCodon <- DNAStringSet(chr, start= annotations$cdsStart + 1, end = annotations$cdsStart + 3)
  startCodon
}

#Returns reverse compl of start codons (as it point to the - strand) in chromosome chr using given annotations
startCodonsMinus<-function(chr, annotations){
  startCodon <- reverseComplement(DNAStringSet(chr, start= annotations$cdsStart -2, end = annotations$cdsStart))
  startCodon
}

annotation_21_plus <- codingAnnotationsByChromosomeAndStrand('chr21', '+')
annotation_21_minus <- codingAnnotationsByChromosomeAndStrand('chr21', '-')
annotation_22_plus <- codingAnnotationsByChromosomeAndStrand('chr22', '+')
annotation_22_minus <- codingAnnotationsByChromosomeAndStrand('chr22', '-')

#Select the regions corresponding to the first condon (cdsStart)
startCodon21plus <- startCodonsPlus(chr21, annotation_21_plus)
startCodon21minus <- startCodonsMinus(chr22, annotation_21_minus)

startCodon22plus <- startCodonsPlus(chr22, annotation_22_plus)
startCodon22minus <- startCodonsMinus(chr22, annotation_22_minus)

table(startCodon21minus)
table(startCodon21plus)

table(startCodon22minus)
table(startCodon22plus)

#Use the stringr library to count the number of A, C, G and T
counts21 <- str_count(chr21, c('A', 'C', 'G', 'T'))
counts22 <- str_count(chr22, c('A', 'C', 'G', 'T'))

#Compute the G+C percentage for each chromosome
CGpercentage <- function(chr){
  total <- str_count(chr, c('A', 'C', 'G', 'T'))
  c_g_total <- total[2] + total[3]
  c_g_percentage <- (c_g_total * 100) / sum(total)
  c_g_percentage
}

cgPercentage21 <- CGpercentage(chr21)
cgPercentage22 <- CGpercentage(chr22)

#Plot the first codon distribution
#as explained we want to identify the first codons after the start codons


#Returns start codons in chromosome chr using given annotations
firstCodonsPlus<-function(chr, annotations){
  startCodon <- DNAStringSet(chr, start= annotations$cdsStart + 1, end = annotations$cdsStart + 6)
  startCodon <- startCodon[grep('ATG', startCodon),]
  startCodon <- substring(startCodon, first = 4, last = 6)
  startCodon
}

firstCodon21plus <- firstCodonsPlus(chr21, annotation_21_plus)
firstCodon21minus <- firstCodonsPlus(chr21, annotation_21_minus)

firstCodon22plus <- firstCodonsPlus(chr22, annotation_22_plus)
firstCodon22minus <- firstCodonsPlus(chr22, annotation_22_minus)

# output all histograms to pdf
pdf("firstCodonDistribution.pdf")

plot(table(firstCodon21plus), main = 'First codon distribution chr 21 +', xlab = 'First codon', ylab = 'Frequency')
barplot(table(firstCodon21minus), main = 'First codon distribution chr 21 -', xlab = 'First codon', ylab = 'Frequency')

plot(table(firstCodon22plus), main = 'First codon distribution chr 22 +', xlab = 'First codon', ylab = 'Frequency')
barplot(table(firstCodon22minus), main = 'First codon distribution chr 22 -', xlab = 'First codon', ylab = 'Frequency')


dev.off()