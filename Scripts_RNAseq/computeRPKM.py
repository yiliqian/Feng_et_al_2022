import csv
import os
import pandas as pd

gbkDir = '../refGenomes'
inputDir = '../05counts'
resultsDir = '../06computeRPKM'

if not os.path.exists(resultsDir):
	os.makedirs(resultsDir)

# List all count files
fileList = []
for file in os.listdir(inputDir):
	if file.endswith('count'):
		fileList.append(file)
fileList = [file.split('.')[0] for file in fileList]

# File to store total # of counted reads
masterDF = pd.DataFrame(index = [], columns = ['Counts'])

#From Genbank files, create a mapping between locus IDs and products

proteinDict = {}

with open(gbkDir+'/BU.csv', mode='r') as infile:
	csvReader = csv.reader(infile)
	proteinDict = {rows[0]:rows[1] for rows in csvReader}

### Convert feature counts to RPKMs
### For pairs, be sure to normalize to reads mapped to that genome, to account for differential abundance within the pair
### Extract product information from the corresponding genbank file

# samples in fileList are monospecies only
for file in fileList:

	dataFrame = pd.read_csv(inputDIr+'/'+file+'.count', index_col=0, skiprows=[0], sep='\t')
	dataFrame = dataFrame.rename(columns={ dataFrame.columns[5]: 'Counts'}) # rename columns
	
	mappedReads = dataFrame.sum(axis=0)['Counts']
	
	# Compute total number of mapped reads
	masterDF.at[file, 'Counts'] = mappedReads
	mappedReads = mappedReads / 1000000 # convert to million reads basis for M in RPKMs
	
	dataFrame['Length'] = (dataFrame['End'] - dataFrame['Start'] + 1) / 1000 # convert to kilobase gene length for K in RPKMs
	
	dataFrame['RPKM'] = dataFrame.iloc[:, 5] / dataFrame['Length'] / mappedReads

# Add protein products 
for locus in dataFrame.index:
	dataFrame.at[locus, 'Product'] = proteinDict[str(locus)]
	
	dataFrame.to_csv(resultsDir+'/'+file+'.rpkm')

masterDF.to_csv(inputDir+'/counts.txt')