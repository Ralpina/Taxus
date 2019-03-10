rareFunction <-
  function(allelicFrequencies,
           maxValue) {
    rare <- list()
    #For each row in allelicFrequencies
    for (row in 1:nrow(allelicFrequencies)) {
      #For each column in allelicFrequencies
      for (col in 1:ncol(allelicFrequencies)) {
        #If the cell value is < maxValue and >0
        if (allelicFrequencies[row, col] < maxValue &&
            allelicFrequencies[row, col] > 0) {
          #Add the column name to a list of vector and use the row name for the list's name
          rare[[rownames(allelicFrequencies)[row]]] <-
            c(rare[[rownames(allelicFrequencies)[row]]], colnames(allelicFrequencies)[col])
        }
      }
    }
    return(rare)
  }