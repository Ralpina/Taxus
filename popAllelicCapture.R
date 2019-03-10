########################################
######Population allelic capture########
########################################

library(poppr)

##population at Wyndcliff, one of our wild reference populations

Wyndcliff <- read.genalex("Wyndcliff.txt", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = "\t", recode = FALSE)
ShuffleWynd <- replicate(999, {shufflepop(Wyndcliff, method = 4)})
##the following code computes the number of alleles occurring in the dataset, as averaged over all the reshuffled replicates of ##the dataset

result <- 0
for (element in ShuffleWynd) {
  colFound <- numeric(length = 0)
  matrix <- element$tab
  for (col in 1:ncol(matrix)) {
    if (is.element(col, colFound)) {
      break
    }
    for (row in 1:nrow(matrix)) {
      if (!is.na(matrix[row, col]) && matrix[row, col] != "0") {
        colFound <- c(colFound, col)
        result <- result + 1
        break
      }
    }
  }
}
result <- result/length(ShuffleWynd)
