allelicCount <-
  function(fileName) {
    # Calculates the mean among 999 replicates
    # Args:
    # fileName: the file name with extension
    #Returns: A list from containing the mean allelic counts across replicates, from 1 genotype to n-genotypes (where n equals to the number of genotypes in the dataset).
    
    #Seed dataset with no population subdivision, as imported from GenAlEx:
    Seed <-
      read.genalex(
        fileName,
        ploidy = 2,
        geo = FALSE,
        region = FALSE,
        genclone = TRUE,
        sep = "\t",
        recode = FALSE
      )
    
    #Poppr function for multilocus permutation:
    ShuffledSeed <- replicate(999, {
      shufflepop(Seed, method = 4)
    })
    
    spidResult <- list()
    
    #For each reshuffled dataset
    for (i in ShuffledSeed) {
      sumNonEmpty <- 0
      matrix <- i$tab
      colFound <- numeric(length = 0)
      #Loop each element of the dataset row by row
      for (row in 1:nrow(matrix)) {
        for (col in 1:ncol(matrix)) {
          if (!is.element(col, colFound)) {
            #... and checks that each cell of the dataset contains a value that is Not 'NA' or '0'
            if (!is.na(matrix[row, col]) && matrix[row, col] != "0") {
              #When the allele occurrence is found, increases the total number of occurrences by one.
              sumNonEmpty <- sumNonEmpty + 1
              colFound <- c(colFound, col)
            }
          }
        }
        #Concat. the number of non-empty rows of each ShuffledSeed
        rowNameFound <- rownames(matrix)[row]
        spidResult[[rowNameFound]] <-
          c(spidResult[[rowNameFound]], sumNonEmpty)
        
        
      }
    }
    
    #For each element of the result, return the mean
    i <- 0
    meanResult <- list()
    for (res in spidResult) {
      i <- i+1
      meanResult[i] <- mean(res)
    }
    
    return (meanResult)
    
  }
