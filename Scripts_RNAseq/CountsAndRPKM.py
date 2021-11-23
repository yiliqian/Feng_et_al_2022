import os
import pandas as pd

# Link species to samples
sampleIdDict = {
        'BU': ['Mucin-1']
        }
speciesList = sampleIdDict.keys()

# Each monospecies should have a corresponding color
sampleColorDict = {
    'BU': 'teal'
        }

# For each species read in RPKM values and create a master dataframe
rawDataDir = '/mnt/rdrive/venturelli/General/Sarah/RNAseq/Jun/Jun_RNAseq20210810/Data/06computeRPKM'
resultsDir = '/mnt/rdrive/venturelli/General/Sarah/RNAseq/Jun/Jun_RNAseq20210810/Data/07CountsAndRPKM/RPKM-summaries'

if not os.path.exists(resultsDir):
    os.makedirs(resultsDir)

for species in speciesList:
    sampleList = sampleIdDict[species]
    masterDF = pd.DataFrame(index = [], columns = sampleList)

    ### Read in RPKM of each sample and assign to masterDF
    for sample in sampleList:
        sampleDF = pd.read_csv(rawDataDir+'/'+sample+'.rpkm', index_col=0)
        masterDF[sample] = sampleDF['RPKM']

        masterDF.to_csv(resultsDir+'/'+species+'-rpkm.csv')

# For each species read in count values and create a master dataframe

rawDataDir = '/mnt/rdrive/venturelli/General/Sarah/RNAseq/Jun/Jun_RNAseq20210810/Data/05counts'
resultsDir = '/mnt/rdrive/venturelli/General/Sarah/RNAseq/Jun/Jun_RNAseq20210810/Data/07CountsAndRPKM/counts-summaries'

if not os.path.exists(resultsDir):
    os.makedirs(resultsDir)

for species in speciesList:
    sampleList = sampleIdDict[species]
    masterDF = pd.DataFrame(index = [], columns = sampleList)

    ### Read in counts from each sample and assign to masterDF
    for sample in sampleList:

        sampleDF = pd.read_csv(rawDataDir+'/'+sample+'.count', index_col=0, skiprows=[0], sep='\t')
        sampleDF = sampleDF.rename(columns={ sampleDF.columns[5]: 'Counts'}) # rename column

        masterDF[sample] = sampleDF['Counts']

    masterDF.to_csv(resultsDir+'/'+species+'-counts.csv')

