###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

PYTHON.PATH <- args[1]

n_Granularity <- as.numeric(args[2]) #THIS VALUE IS USED TO DETERMINE THE NUMBER OF REFERENCE BIOMARKERS YOU WOULD LIKE TO USE.
#THE HIGHER THIS VALUE IS, THE MORE GRANULARITY YOU WILL HAVE IN YOUR REFERENCE LISTS. THE LOWER IT IS, THE MORE STABLE GENES YOU WILL FIND.
#RECOMMENDED DO NOT EXCEED A VALUE OF 25.

print("--------------------------------------------------")
print("Step 2 - Find Reference Genes")
print("--------------------------------------------------")

print("Reading input files.")
norm.table <- read.csv("results/normalized_table.csv", row.names = 1)
out <- read.csv("final_out.csv", row.names = 1)

intermediate <- norm.table[row.names(norm.table) %in% row.names(out), ]

print("Determining reference ranges.")
propData <- as.numeric(c()) #CREATING INPUT DATA TO BE USED IN ANOTHER SCRIPT - REF_SELF_TEST_SYMBOL
percent <- 1 / (n_Granularity + 1)

for (i in c(1:n_Granularity)){ #DETERMINES THE LOW REFERENCE RANGE FOR AN EVEN QUANTILE DISTRIBUTION
  propData <- append(propData, (percent*i) - 0.02)
}
for (i in c(1:n_Granularity)){ #DETERMINES THE HIGH REFERENCE RANGE FOR AN EVEN QUANTILE DISTRIBUTION
  propData <- append(propData, (percent*i) + 0.02)
}

refRange <- quantile(as.numeric(data.matrix(intermediate)), probs = propData)
outPrint <- as.data.frame(refRange)

sink("data/referenceRanges.csv")
cat(outPrint$refRange)
sink()

#DATA FORMATTING
reference_gene = system2(PYTHON.PATH, args="src/ref_sel_test_symbol.py",stdout=T)

print(paste0("Reference genes: ",paste(reference_gene, collapse = ",")))
print(paste0("* ^^^ There should be ", paste0(n_Granularity, " results here. ^^^ *")))
ref_cpm <- norm.table[reference_gene,]

#CREATES THE GUIDELINE COUNTS PER MILLION TO ASSOCIATE WITH 1-5 CATEGORIZATION
write.csv(ref_cpm, "results/ref_cpm.csv")
print("Saving a bin for the reference genes.")
save(reference_gene, file = "ref.gene.bin")
