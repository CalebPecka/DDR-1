#####EdgeR###############
library(edgeR)

###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

inputFilePath <- args[1]

#PRIMARY INPUT FILE. ASSUMES THAT FIRST N COLUMNS ARE OF 1 GROUP, AND REMAINING ARE STORED AS OTHER
print("--------------------------------------------------")
print("Step 1 - Calculate Stats")
print("--------------------------------------------------")
print("Reading input file.")
complete <- read.csv(inputFilePath, row.names = 1)
print("Input file read successfuly - Dimensions of file:")
print(dim(complete))

#PROCESSING
print("Normalizing input table.")
#grp <- as.factor(rep(1, ncol(complete)))
#y <- DGEList(complete, group = grp)
#y <- calcNormFactors(y, method="TMM")
norm.table <- complete

#GENERATES CORE STATISTICS FOR THE NORMALIZED TABLE
print("Determining core metrics")
cov <- apply(norm.table, 1, sd)/apply(norm.table, 1, mean)
mean <- apply(norm.table, 1, mean)
std <- apply(norm.table, 1, sd)
MFC <- apply(norm.table, 1, max)/apply(norm.table, 1, min)

#WRITES THE NORMALIZED TABLE
write.csv(norm.table,"results/normalized_table.csv")
out <- data.frame(cov, mean, std, MFC)
out <- out[!is.infinite(out$MFC),]
out <- out[!is.na(out$MFC),] #REMOVES NA VALUES WHICH PRODUCE ERRORS IN READING

out <- cbind("genes"= rownames(out), out)

out$MFC <- abs(out$MFC)
out$cov <- abs(out$cov)

#BINDS WITH GENE NAMES
print("Writing finalized normal table - final_out.csv")
write.csv(out, "final_out.csv", row.names=F)
