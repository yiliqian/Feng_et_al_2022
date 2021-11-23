# Import Packages
library(DESeq2)

# Define folder structure

countDir = './results/counts-summaries'
rlogDir = './results/rlog'
deDir = './results/deseq'

# Create output directories if necessary
dir.create(file.path(rlogDir))
dir.create(file.path(deDir))

# Species
speciesList = c('BU')

# For each sample, read in the count file and perform a rlog transformation

for (species in speciesList)
{
 # Read in counts-summaries
  countDF = read.table(file.path(countDir, paste(species, "-counts.csv", sep = '')),
					  header=TRUE, sep=',', row.names=1)
 # Transform the data
 rlogDF = rlog(data.matrix(countDF), blind=TRUE)
 
 # Convert back to a matrix and add the index
 rlogDF = as.data.frame(rlogDF)
 row.names(rlogDF) = row.names(countDF)
 
 # Write to file
 write.table(rlogDir, file = file.path(rlogDir, paste(species, ".rlog", sep='')), sep=',')
}

# For each sample, read in the count file and perform a rlog transformation

for (species in speciesList)
{
  # Read in counts
  countDF = read.table(file.path(countDir, paste(species, "-counts.csv", sep = '')), 
                       header=TRUE, sep=',', row.names=1)

  # Define sample information
  sampleInfo = data.frame(condition = c('BMM', 'BMM', 'Mucin', 'Mucin'),
						  row.names = colnames(countDF))
  
  # Create dataset for DE analysis
  deSeqObj = DESeqDataSetFromMatrix(countData = as.matrix(countDF),
                                   colData = sampleInfo,
								   design = ~ condition
								   )
  
  # Perform DEanalysis - comparisons to the Mucin
  deSeqObj$condition <- relevel(deSeqObj$condition, ref = 'Mucin')
  deSeqResults = DESeq(deSeqObj)
  
  # Extract comparison (Mucin is control)
  # Shrink fold-change estimates
  MucinToBMM = lfcShrink(deSeqResults, coef='condition_Mucin_vs_BMM', type='apeglm')
  write.csv(as.data.frame(MucinToBMM), file = file.path(deDir, paste(species, "-MucinVsBMM.csv", sep = '')))
}