setwd("D:\\workspace\\bioinformatics\\BIOLF449\\INFOF434")

#1. Download the hg19 RefSeq genes annotation file from UCSC
#genome.ucsc.edu / Tools / Table Browser / refGene Table

#Data was downloaded in data file next to TP_3.r


#2. Load table data in variable t
t <- read.table('data', header = T, comment.char = '')[,-1]

#3. Filter the coding genes
coding <- t[grep('NM', t$name),]

#4. Generate a table of the transcripts length (from TSS = Transcription start site 
#to TES = transcription end site)

# First approach: simply make the difference between columns
deltaApproach1 <-  coding$txEnd -coding$txStart

# Second approach : use apply function
# Numeric vector values are converted to string so we need to reconvert them to 
# integer in order to perform the operation
deltaLength <- function(v){
  as.integer(v['txEnd']) - as.integer(v['txStart'])
}
deltaApproach2 <- apply(coding, 1 , deltaLength)

# Check both appoaches produce equal results
all (deltaApproach1 == deltaApproach2)

#5. Plot an histogram of gene length distribution

# output all histograms to pdf
pdf("histograms.pdf")

hist(deltaApproach1, nclass =100, col = 'red', main ="Histogram of coding genes length distribution",
     xlab = 'Gene Length', ylab = 'Frequency')

hist(deltaApproach1, nclass =300, col = 'red', main ="Histogram of coding genes length distribution (xlim = 200000)",
     xlab = 'Gene Length', ylab = 'Frequency', xlim = c(0,200000))

#6. Select all the long non coding transcript (length >= 2000bp)
nonCoding <- t[grep('NR_',t$name),]

nonCodingLengths <- nonCoding$txEnd -nonCoding$txStart

#Approach 1
deltaNonCoding <- nonCodingLengths[nonCodingLengths[] > 20000]
# or approach 2
deltaNonCoding2 <- subset(deltaNonCoding, deltaNonCoding[] > 20000)

# Check both appoaches produce equal results
all (deltaNonCoding == deltaNonCoding2)

#7. Plot an histogram of long non coding gene length distribution
hist(deltaNonCoding, nclass =100, col = 'blue', main ="Histogram of long non coding gene length distribution",
     xlab = 'Gene Length', ylab = 'Frequency')

hist(deltaNonCoding, nclass =300, col = 'blue', main ="Histogram of long non coding gene length distribution (xlim=200000)",
     xlab = 'Gene Length', ylab = 'Frequency', xlim = c(0,200000))


#8. Compare the coding / long non coding transcript distribution

hist(deltaApproach1, nclass =100, col = 'red', main ="Coding/long non coding transcript distribution",
     xlab = 'Gene Length', ylab = 'Frequency')
hist(deltaNonCoding, add=T, nclass =100, col = 'blue')
legend("topright", c("Long non coding lengths", "Coding lengths"), fill=c("blue", "red"))


hist(deltaApproach1, nclass =200, col = 'red', main ="Coding / long non coding transcript distribution (xlim=200000)",
     xlab = 'Gene Length', ylab = 'Frequency', xlim = c(0,200000))
hist(deltaNonCoding, add=T, nclass =200, col = 'blue', xlim = c(0,200000))
legend("topright", c("Long non coding lengths", "Coding lengths"), fill=c("blue", "red"))

dev.off()


