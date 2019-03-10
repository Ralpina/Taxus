locallyCommonFunction <-
  function(allelicFrequencies,
           maxValue,
           minValue) {
    locallyCommon <- list()

    # for each column of allelicFrequencies
    for (col in 1:ncol(allelicFrequencies)) {
      # Count number of values < minValue
      countSmall <- 0
      # Count number of values >= maxValue
      countBig <- 0
      
      tmpResult <- list()
      # for each row of allelicFrequencies
      for (row in 1:nrow(allelicFrequencies)) {
        # If there is at least one row for the current column that is < maxValue and >= minValue, we can skip this column and analyse the next one
        if (allelicFrequencies[row, col] < maxValue &&
            allelicFrequencies[row, col] >= minValue) {
          break
          # If there is one value >= minValue is candidate to be a valid result, so we save it in a new vector called as the allelicFrequencies's row name and with value = the allelicFrequencies's column name
        } else if (allelicFrequencies[row, col] >= maxValue &&
                   length(tmpResult) == 0) {
          tmpResult[[rownames(allelicFrequencies)[row]]] <-
            colnames(allelicFrequencies)[col]
          countBig <- countBig + 1
        } else if ((allelicFrequencies[row, col] < minValue)) {
          countSmall <- countSmall + 1
        }
      }
      
      # If at the end of the rows we find eight small values plus one big value, this is a valid result
      if (countSmall == nrow(allelicFrequencies) - 1 &&
          countBig == 1) {
        locallyCommon[[names(tmpResult)]] <-
          c(locallyCommon[[names(tmpResult)]], tmpResult[[1]])
      }
    }
    
    return (locallyCommon)
    
  }