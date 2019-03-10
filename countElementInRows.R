# Given a list with named rows, returns a new list with the length of each row.
countElementInRows <-
  function(namedList) {
    countList <- list()
    # For each row in the givenList count the length
    for (element in names(namedList)) {
      countList[[element]] <- length(namedList[[element]])
    }
    return(countList)
  }