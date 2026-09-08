from collections import defaultdict 
import glob
import pandas as pd
import pprint
import pysam
import statistics

DEDUP_FILES = sorted(glob.glob('outFilesR2/*_Dedup.bam'))
GFF_FILE = 'reference/gencode_M37/gencode.vM37.annotation.gff3'
OUT_FILE_FREQ = 'ET/PhotoSeq_TechnicalReplicate_GeneBarcodeCounts.csv'
OUT_FILE = 'ET/PhotoSeq_TechnicalReplicate_GeneCounts.csv'

def readGFF():
  genes = {}
  
  print('Reading GFF file...')
  with open(GFF_FILE) as gffHandle:
    for line in gffHandle:
      if line[0] != '#' and line.split('\t')[2] == 'gene':
        genes[line.split('gene_id=')[1].split(';')[0]] = line
  
  return genes

def main():
  filesToAnalyze = DEDUP_FILES

  genes = readGFF()
  allGenes = []
  allConditionBarcodes = []
  frequencies = defaultdict(lambda: defaultdict(int)) # default to 0 count

  # Create dictionary of gene frequencies
  for geneFile in filesToAnalyze:  
    print('Reading %s...' % geneFile)

    condition = geneFile.split('/')[-1].split('_Dedup')[0]
    samFile = pysam.AlignmentFile(geneFile, "rb")
    samIter = samFile.fetch(until_eof=True)

    for read in samIter:
      barcodeSeq = read.query_name.split('_')[1]
      geneName = read.get_tag('XT')

      condBar = '%s/%s' % (condition, barcodeSeq)

      frequencies[geneName][condBar] += 1

      if not condBar in allConditionBarcodes:
        allConditionBarcodes.append(condBar)

  # Now, create CSV
  allData = defaultdict(list)

  allData['Gene'] = frequencies.keys()

  print('Compiling data...')
  for gene in allData['Gene']:
    allData['Gene name'].append(genes[gene].split('gene_name=')[1] \
                                .split(';')[0])
    for condBar in allConditionBarcodes:
      allData[condBar].append(frequencies[gene][condBar])

  photoseq_counts = pd.DataFrame.from_dict(allData)
  photoseq_counts.to_csv(OUT_FILE_FREQ, index=False)

  # Reorder data by PhotoSeq sample, ROI barcode, and replicate
  del photoseq_counts['Gene']


  photoseq_counts.rename(columns={'Gene name': 'Gene',
                             'EL1rm/TATGGA': 'Breast_early_N1', 
                             'E1Trm/TATGGA': 'Breast_early_N1_techrep2'},
                    inplace=True)

  photoseq_counts = photoseq_counts[['Gene', 'Breast_early_N1', 'Breast_early_N1_techrep2']]

  photoseq_counts.set_index('Gene', inplace=True)
  # Add together counts with duplicate indices
  photoseq_counts = photoseq_counts.groupby(photoseq_counts.index).sum()

  photoseq_counts.to_csv(OUT_FILE)

if __name__ == '__main__':
  main()
