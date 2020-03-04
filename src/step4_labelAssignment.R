###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

group1.sampleSize <- as.numeric(args[1]) #NUMBER OF ELEMENTS IN MUTATED GROUP (FIRST GROUP). REMAINING ELEMENTS GO TO GROUP 2
n <- as.numeric(args[2]) #NUMBER OF GENES TO BE ANALYZED AT END OF SEQUENCE. INCLUDES TOP 'n' AND BOTTOM 'n'

print("--------------------------------------------------")
print("Step 4: Label Assignment")
print("--------------------------------------------------")

#FIND SIGNIFICANT GENES
print("Reading input files.")
geneList <- read.csv("results/overlap_test_fdr_05_RNASeq.csv", row.names = 1)
geneList <- geneList[order(geneList$ED), ]

print("Finding list of significant genes.")
norm.table <- read.csv("results/normalized_table.csv", row.names = 1)

sigLow <- as.character(head(geneList$row.names.WNT_cpm., n = n))
#sigLow <- geneList$row.names.TN_cpm.[0:n]
sigHigh <- as.character(tail(geneList$row.names.WNT_cpm., n = n))
#sigHigh <- geneList$row.names.TN_cpm.[(length(geneList$ED) - (n-1)) : (length(geneList$ED))]
sigGenes <- c(as.character(sigLow), as.character(sigHigh))

#SUBSETS THE NORMALIZED TABLE INTO FEATURES WE CARE ABOUT
nameSearch <- row.names(norm.table)

ref_cpm <- read.csv("results/ref_cpm.csv", row.names = 1)

feature <- norm.table[row.names(norm.table) %in% sigGenes, ]

#SEARCH THROUGH CATEGORIZATIONS TO FIND PROPER ASSIGNMENT
print("Assigning features to reference gene category.")
for (i in 1:nrow(feature)){
  for (j in 1:ncol(feature)){
    if (as.numeric(feature[i,j]) <= as.numeric(ref_cpm[1,j])){
      feature[i,j] = 0
    }
    for (k in c(1:(nrow(ref_cpm)-1))){
      if (as.numeric(feature[i,j]) > as.numeric(ref_cpm[k,j]) && as.numeric(feature[i,j]) <= as.numeric(ref_cpm[(k+1),j])){
        feature[i,j] = as.numeric(k)
      }
    }
    if (as.numeric(feature[i,j]) > as.numeric(ref_cpm[nrow(ref_cpm),j])) {
      feature[i,j] = nrow(ref_cpm)
    }
  }
}

#ADD LABELS OF GROUPINGS
print("Writing features and labels.")
groupings <- rep("group2",ncol(feature))
groupings[1:group1.sampleSize] <- "group1"

#anno <- as.data.frame(groupings)
#row.names(anno) <- colnames(feature)
#pheatmap(feature, show_colnames = F, cluster_cols = F, annotation_col = anno)

feature <- cbind(t(feature),label=groupings)
write.csv(feature, "results/feature_and_labels.csv",row.names = F)
