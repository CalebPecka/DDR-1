import codecs
import csv
import pandas
import sys, operator

def getRefExpLevel(lowerTh, upperTh, outputfile):
    with codecs.open("final_out.csv", "r", encoding='utf-8', errors='ignore') as f:
        with open(outputfile, 'w') as csvfile:
          reader = csv.DictReader(f)
          writer = csv.writer(csvfile)
          writer.writerow(["genes", "cov", "mean", "std", "MFC", "covMFC"])
          for row in reader:
            #try:

                ens = []
                if lowerTh <= pandas.to_numeric(row["mean"]) <= upperTh:
                  covMFC = float(row["MFC"])*float(row["cov"])
                  Account = [row["genes"], row["cov"], row["mean"], row["std"], row["MFC"], covMFC]
                #print(row["Gene_ID"] + " " + mg.getgene(row["Gene_ID"])['symbol'])
                  writer.writerow(Account)

def findTop(resultFile):
    df  = pandas.read_csv(resultFile)
    sortedlist = df.sort_values(by='covMFC')
    for index, row in sortedlist.iterrows():
        if row[0] in lines:
            print(row[0])
            break

text_file = open("data/housekeeping_symbol.txt", "r")
lines = text_file.read().splitlines()

#READS THE SPLIT QUANTILES FROM STEP 2
refRange = open("data/referenceRanges.csv", "r")
ULrange = refRange.read().split()
#PROCESSING - REMOVING LABELS
rangeSize = int(len(ULrange) / 2)
ULrange = [float(i) for i in ULrange]

lows = []
ups = []

#APPENDING TO VALUES TO THEIR CORRESPONDING RANGE SETS
for i in range(rangeSize):
	lows.append(ULrange[i])
for i in range(rangeSize):
	ups.append(ULrange[rangeSize+i])

outfiles = []
for i in range(rangeSize):
	outfiles.append('ref' + str(i+1) + '.csv')

for i in range(len(outfiles)):
    getRefExpLevel(lows[i],ups[i],outfiles[i])

for outfile in outfiles:
    findTop(outfile)

sys.exit(0)
sys.exit(1)
