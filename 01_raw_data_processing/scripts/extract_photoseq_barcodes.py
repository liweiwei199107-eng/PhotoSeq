import glob
import os

IN_FILES = sorted(glob.glob('inFiles/*_R1.fastq.gz'))
OUT_DIR = 'outFiles/'

CLEAR_OUT_DIR = True

def parseSequences(inFile, filePrefix, trimmedR1File): 
  print('Analyzing %s...' % inFile.split('_R1')[0])

  barcodePattern = '^(?P<discard_1>.{27})' \
                   '(?P<umi_1>.{12})' \
                   '(?P<discard_2>.*)' \
                   '(?P<cell_1>GTTAGG|TATGGA|AGGGTA)' \
                   '(?P<discard_3>TGAGTTAT){s<=2}' \
                   '(?P<discard_4>.{8})' \
                   '.{15,89}+' \
                   '(?P<discard_5>AAAAAAAAAAA.*|AGATCGGAAGA.*|.*$)'

  # Extract UMIs and PhotoSeq ROI barcodes
  extractLogFile = '%s_extract.log' % filePrefix
  
  print('  %s Extracting PhotoSeq UMI and ROI barcode information...' % filePrefix)
  os.system(('umi_tools extract --stdin %s --extract-method=regex ' \
              '--bc-pattern="%s" -L %s --stdout %s') % \
              (inFile, barcodePattern, extractLogFile, trimmedR1File))

def main():
  if CLEAR_OUT_DIR:
    os.system('rm -rf %s' % OUT_DIR)
    os.system('mkdir %s' % OUT_DIR)

  for inFile in IN_FILES:
    filePrefix = '%s%s' % (OUT_DIR, inFile.split('/')[-1].split('_')[0])
    trimmedR1File = '%s_R1_trimmed.fastq.gz' % filePrefix
    parseSequences(inFile, filePrefix, trimmedR1File)

if __name__ == '__main__':
  main()
