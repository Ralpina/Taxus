myBeautifulFunc <- function(replicateShuffleDataset, referenceAlleles, referenceAllelesCount) {
  result <- c()
  resMean <- c()
  inutile <- c()
  
  for (element in replicateShuffleDataset) {
    numRow <- 5
    count <- 0
    matrix <- element$tab
    while (numRow <= nrow(matrix)) {
      sumNonEmpty <- 0
      count <- count + 1
      
      #For each column of the dataset up to the maximum number of columns...
      for (col in 1:ncol(matrix)) {
        #...gets a subSet of rows...
        for (row in 1:numRow) {
          #... and checks that each cell of the dataset contains a value that is Not 'NA' or '0'
          if (is.element(colnames(matrix)[col], referenceAlleles) &&
              !is.na(matrix[row, col]) && matrix[row, col] != "0") {
            #When the allele occurrence is found, the function increases the total number of occurrences by one
            sumNonEmpty <- sumNonEmpty + 1
            break
          }
        }
      }
      if (is.null(result[count]) || is.na(result[count])) {
        result[count] <- sumNonEmpty
      } else{
        result[count] <- result[count] + sumNonEmpty
      }
      numRow <- numRow + 5
      if (numRow < nrow(matrix) + 5 &&  numRow  > nrow(matrix)) {
        numRow <- nrow(matrix)
      }
    }
  }
  i <- 0
  for (res in result) {
    i <- i + 1
    resMean[i] <- res / length(replicateShuffleDataset)
    
    inutile[i] <-
      res / length(replicateShuffleDataset) / referenceAllelesCount
  }
  return(list(res=result, mean=resMean, meanInCount=inutile))
}