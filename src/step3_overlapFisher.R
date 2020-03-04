###### INPUT ARGUEMENTS#########
args <- commandArgs(trailingOnly = T)

group1.sampleSize <- as.numeric(args[1])

print("--------------------------------------------------")
print("Step 3 - Fisher Overlap Tests")
print("--------------------------------------------------")

print("Reading gene bin.")
load("ref.gene.bin")
print("Reading normalized table.")
norm.table <- read.csv("results/normalized_table.csv", row.names = 1)
ref_cpm <- read.csv("results/ref_cpm.csv", row.names = 1)
n <- nrow(ref_cpm)

WNT_cpm <- norm.table[, 1:group1.sampleSize]
OT_cpm <- norm.table[, (group1.sampleSize+1):ncol(norm.table)]
  
#################Caterories################################
WNT_ref <- WNT_cpm[reference_gene, ] 
OT_ref <- OT_cpm[reference_gene, ] 

cat_count_total <- data.frame(row.names(WNT_cpm))
  
print("Creating groupings: 0% completed")
WNT_C0 <- c()
#out <- c()
for (i in 1:nrow(WNT_cpm)){
  test <- c()
  for (j in 1:ncol(WNT_cpm)){
    
    if(as.numeric(WNT_cpm[i,j]) < as.numeric(WNT_ref[1,j])){ #&& as.numeric(WNT_cpm[i,j])<=as.numeric(WNT_ref[4,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  WNT_C0[i] <- sum(test)
}
cat_count_total$WNT_C0 <- WNT_C0

for (k in c(1:(n-1))){
  print(k / (2 * n))
  print("Creating groupings: 10% completed")
  #assign(paste("WNT_C", k, sep = ""), c())
  alphaList <- c()
  #WNT_C1 <- c()
  #out <- c()
  
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[k,j]) && as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[k+1,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    alphaList <- append(alphaList, sum(test))
    #assign(paste("WNT_C", k, sep = ""), alphaList)
  }
  cat_count_total[,paste("cat_count_total$WNT_C", k, sep = "")] <- alphaList
}

print("Creating groupings: 50% completed")
#WNT_C4 <- c()
alphaList <- c()
#out <- c()
for (i in 1:nrow(WNT_cpm)){
  test <- c()
  for (j in 1:ncol(WNT_cpm)){
    
    if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[k,j])){#&& as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[5,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  alphaList <- append(alphaList, sum(test))
}

cat_count_total[,paste("cat_count_total$WNT_C", n, sep = "")] <- alphaList

#THIS COMPLETES THE GROUPINGS FOR WILD TYPES. NOW WE MUST MAKE GROUPINGS FOR THE NORMAL SAMPLES.

print("Creating groupings: 0% completed")
OT_C0 <- c()
#out <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) < as.numeric(OT_ref[1,j])){ #&& as.numeric(WNT_cpm[i,j])<=as.numeric(WNT_ref[4,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  OT_C0[i] <- sum(test)
}
cat_count_total$OT_C0 <- OT_C0

for (k in c(1:(n-1))){
  print(k / (2 * n))
  print("Creating groupings: 10% completed")
  #assign(paste("WNT_C", k, sep = ""), c())
  alphaList <- c()
  #WNT_C1 <- c()
  #out <- c()
  
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[k,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[k+1,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    alphaList <- append(alphaList, sum(test))
    #assign(paste("WNT_C", k, sep = ""), alphaList)
  }
  cat_count_total[,paste("cat_count_total$OT_C", k, sep = "")] <- alphaList
}

print("Creating groupings: 90% completed")
#WNT_C4 <- c()
alphaList <- c()
#out <- c()
for (i in 1:nrow(OT_cpm)){
  test <- c()
  for (j in 1:ncol(OT_cpm)){
    
    if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[k,j])){#&& as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[5,j])){
      test[j] <- 1
    }
    else{
      test[j] <- 0
    }
  }
  alphaList <- append(alphaList, sum(test))
}

cat_count_total[,paste("cat_count_total$OT_C", n, sep = "")] <- alphaList






#CODE SEGMENT NO LONGER USED
if (FALSE){
  print("Creating groupings: 20% completed")
  WNT_C2 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[2,j]) && as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[3,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C2[i] <- sum(test)
  }
  
  print("Creating groupings: 30% completed")
  WNT_C3 <- c()
  #out <- c()
  for (i in 1:nrow(WNT_cpm)){
    test <- c()
    for (j in 1:ncol(WNT_cpm)){
      
      if(as.numeric(WNT_cpm[i,j]) >= as.numeric(WNT_ref[3,j]) && as.numeric(WNT_cpm[i,j])< as.numeric(WNT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    WNT_C3[i] <- sum(test)
  }
  
  print("Creating groupings: 50% completed")
  OT_C0 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) < as.numeric(OT_ref[1,j])){ #&& as.numeric(OT_cpm[i,j])<=as.numeric(OT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C0[i] <- sum(test)
  }
  
  print("Creating groupings: 60% completed")
  OT_C1 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[1,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[2,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C1[i] <- sum(test)
  }
  
  print("Creating groupings: 70% completed")
  OT_C2 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[2,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[3,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C2[i] <- sum(test)
  }
  
  print("Creating groupings: 80% completed")
  OT_C3 <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[3,j]) && as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[4,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
    OT_C3[i] <- sum(test)
  }

  print("Creating groupings: 90% completed")
  #OT_C4 <- c()
  #assign(paste("WNT_C", n, sep = ""), c())
  alphaList <- c()
  #out <- c()
  for (i in 1:nrow(OT_cpm)){
    test <- c()
    for (j in 1:ncol(OT_cpm)){
      
      if(as.numeric(OT_cpm[i,j]) >= as.numeric(OT_ref[n,j])){#&& as.numeric(OT_cpm[i,j])< as.numeric(OT_ref[5,j])){
        test[j] <- 1
      }
      else{
        test[j] <- 0
      }
    }
  alphaList <- append(alphaList, sum(test))
  }

  cat_count_total[,paste("cat_count_total$WNT_C", n, sep = "")] <- alphaList
  
}

print("Creating groupings: 100% completed")
#cat_count_total <- data.frame(row.names(WNT_cpm), WNT_C0, WNT_C1, WNT_C2, WNT_C3, WNT_C4, OT_C0, OT_C1, OT_C2, OT_C3, OT_C4)

print("Performing Fisher tests.") 
overlap_p <- c()
for (i in 1: nrow(cat_count_total)){
  A <- as.numeric(cat_count_total[i,2:(2+n)]) #WILD TYPE GROUP
  B <- as.numeric(cat_count_total[i,(3+n):(3+n+n)]) #NORMAL TYPE GROUP
  #D <- as.numeric(cat_count_total[i,12:16])
  #E <- as.numeric(cat_count_total[i,17:21])
  tab=as.table(rbind(A,B))
  row.names(tab)=c('G4','OTHERS')
  c <- fisher.test(tab, workspace=2e+07,hybrid=TRUE)
  overlap_p[i] <- c$p.value
}

#####FOLD CHANGE########
print("Calculating fold change.")
ED <- c()
for (i in 1: nrow(cat_count_total)){
  A <- as.numeric(cat_count_total[i,2:(2+n)]) #WILD TYPE GROUP
  B <- as.numeric(cat_count_total[i,(3+n):(3+n+n)]) #NORMAL TYPE GROUP
  D <- sum(A * seq_along(A))/ncol(WNT_cpm)
  E <- sum(B * seq_along(B))/ncol(OT_cpm)
  #D <- as.numeric(cat_count_total[i,12:16])
  #E <- as.numeric(cat_count_total[i,17:21])
  c <- (D-E)
  ED[i] <- c
}
  
print("Writing overlap tests")
overlap_test <- data.frame(row.names(WNT_cpm),overlap_p, ED )
overlap_test_padjust <- p.adjust(overlap_p, method = "fdr", n = length(overlap_p))
overlap_test_fdr <- data.frame(row.names(WNT_cpm),overlap_test_padjust, ED)
overlap_pvalue_fdr <- overlap_test_fdr[order(ED, decreasing = TRUE), ]
keep <- (overlap_pvalue_fdr$overlap_test_padjust <= 0.25)
write.csv(overlap_pvalue_fdr[keep, ], "results/overlap_test_fdr_05_RNASeq.csv")
