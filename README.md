# Taxus
The following document is explanatory, please also refer to .R documents. 
Total Allelic Capture (see “total.R”)
Wild dataset and seed dataset were prepared with no population (or seed accession) subdivision. We used the function shufflepop implemented in the R package POPPR v2.8 (Kamvar et al., 2014), and the multilocus permutation method (i.e., 4). This method reshuffles genotypes among individuals at each locus, maintaining both allelic frequencies and heterozygosity. Each dataset was reshuffled 999 times (R basic function replicate). For each reshuffled dataset, subsets of individuals were evaluated for allelic count, i.e., starting from 1 individual and progressively increasing the subset by 1, up to the maximum number of individuals in the seed dataset (as the seed dataset < wild dataset). Allelic counts were finally averaged among all replicates (shuffled datasets).
##################################
######Total allelic capture#######
##################################

library(poppr)

##FUNCTION:

allelicCount <-
  function(fileName) {
   
   ## This function computes the mean allelic counts across 999 replicates (multilocus-permuted datasets)
   ## Arguments:
   ## fileName: the file name with extension
   ##Returns: A list from containing the mean allelic counts across replicates, from 1 genotype to n-genotypes (where n equals to the number of genotypes in the dataset).
    

##Code: 

    ##Seed dataset with no population subdivision, as imported from GenAlEx:
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
    
    ##Poppr function for multilocus permutations:
    ShuffledSeed <- replicate(999, {
      shufflepop(Seed, method = 4)
    })
    


    spidResult <- list()
    
    ##For each reshuffled dataset
    for (i in ShuffledSeed) {
      sumNonEmpty <- 0
      matrix <- i$tab
      colFound <- numeric(length = 0)
      ##Loop each element of the dataset row by row
      for (row in 1:nrow(matrix)) {
        for (col in 1:ncol(matrix)) {
          if (!is.element(col, colFound)) {
            ##... and checks that each cell of the dataset contains a value that is Not 'NA' or '0'
            if (!is.na(matrix[row, col]) && matrix[row, col] != "0") {
              ##When the allele occurrence is found, increases the total number of occurrences by one.
              sumNonEmpty <- sumNonEmpty + 1
              colFound <- c(colFound, col)
            }
          }
        }
        ##Concat. the number of non-empty rows of each ShuffledSeed
        rowNameFound <- rownames(matrix)[row]
        spidResult[[rowNameFound]] <-
          c(spidResult[[rowNameFound]], sumNonEmpty)
        
        
      }
    }
    
    ##For each element of the result, returns the mean
    i <- 0
    meanResult <- list()
    for (res in spidResult) {
      i <- i+1
      meanResult[i] <- mean(res)
    }
    
    return (meanResult)
    
  }

The same function was applied to the Wild dataset with no population subdivision, as imported from GenAlEx.
 
Population allelic capture (see “popAllelicCapture.R”)
Datasets divided by population and seed accession were prepared. We used the function shufflepop implemented in the R package POPPR and the multilocus permutation method described above, to obtain 999 permuted datasets.
The code is reported below, considering two datasets (Wyndcliff population and Wyndcliff seed accession). Note that the allelic count in the seed accession is conceived considering a number of genotypes which is equal to the number of genotypes occurring in the wild population, to avoid the bias associated to different sample sizes. To satisfy this requirement, datasets were first reshuffled several times and allelic counts were averaged across the reshuffled replicates, to avoid a genotype-sampling bias. This code can be adjusted to allow for different numbers of genotypes. 
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


##seed accession made at Wyndcliff, to compare to the reference wild population


SWyndcliff <- read.genalex("SWyndcliff.txt", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = "\t", recode = FALSE)
ShuffleSWynd <- replicate(999, {shufflepop(SWyndcliff, method = 4)})

##The following code is similar to the previous one, but it computes the number of alleles being represented in the first 24 ##genotypes, which is the total number of genotypes in our reference dataset, Wyndcliff (above). If the number of genotypes in ##the germplasm accession to compare is equal or lower than the number of genotypes in the reference population, the code above ##can be used.

result <- 0
for (element in ShuffleSWynd) {
     colFound <- numeric(length = 0)
     matrix <- element$tab
     for (col in 1:ncol(matrix)) {
         if (is.element(col, colFound)) {
             break
         }
         for (row in 1:24) {
             if (!is.na(matrix[row, col]) && matrix[row, col] != "0") {
                 colFound <- c(colFound, col)
                 result <- result + 1
                 break
             }
         }
     }
 }
result <- result/length(ShuffleSWynd)

 
Locally common alleles and Rare alleles:
The wild dataset divided by population was used to compute allele frequencies in POPPR with the function makefreq. Two allele categories were then evaluated in the wild dataset.
1.	Locally common alleles: occurring in one population with a frequency > 0.20 and in all the others with a frequency < 0.1. 
2.	Rare alleles: occurring with a frequency < 0.05 (and ≠ 0)
We implemented two functions to obtain the list of alleles falling in these categories and a third function to count the number of alleles falling in each category, respectively. We used this information to evaluate allelic capture in seed accessions, as obtained after the permutations described above.
The codes described below can be adjusted to search for alleles falling in other frequency-categories. 

############################################
######Allelic frequency computations########
############################################

library(poppr)
##importing wild dataset divided by pop:
PopWild <- read.genalex("Pops_Wild.txt", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = "\t", recode = FALSE)
PopWildGI <- genclone2genind(PopWild)
PopWildGP <- genind2genpop(PopWildGI)

##computing allele frequencies:
PopWildFreq <- makefreq(PopWildGP)



See also locallyCommon.R

############################################
#############Locally common alleles######### 
############################################

##FUNCTION:

locallyCommonFunction <-
  function(allelicFrequencies,
           maxValue,
           minValue) {
    locallyCommon <- list()
##locallyCommonFunction returns a list of populations and their respective alleles satisfying the conditions maxValue and
##minValue
##in the reference matrix of allele frequencies, i.e. PopWildFreq in our case
##maxValue refers to the lower threshold frequency > 0.20, which one population has to satisfy to be selected; 
##minValue refers to the upper threshold frequency that all other populations need to satisfy (0.1, in our case). 
##Both values can be adjusted.


##Code: 

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



##Example:

WildLCommon <- locallyCommonFunction(PopWildFreq, 0.20, 0.1)


##Returns:

$`Bryn_Pydew`
[1] "Tax36.186"

 
See rare.R

#################################
############Rare alleles#########
#################################

##FUNCTION:

rareFunction <-
  function(allelicFrequencies,
           maxValue)


##rareFunction returns a list of populations and their respective alleles satisfying the condition maxValue (i.e., rows of 
##allelicFrequencies). In our case, maxValue refers to the threshold of 0.05, which can be adjusted to other values. 

##Code:

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

##Example:

WildRare <- rareFunction(PopWildFreq, 0.05)
 

See countElementInRows.R

##Function to count alleles in each category, i.e., locally common or rare:

countElementInRows <-
  function(namedList)

##Given a list with named rows, the function returns a new list with the length of each row.

##Examples:

CountRare <- countElementInRows(WildRare)
CountLCommon <- countElementInRows(WildLCommon)


##code:

countElementInRows <-
  function(namedList) {
    countList <- list()
    # For each row in the givenList count the length
    for (element in names(namedList)) {
      countList[[element]] <- length(namedList[[element]])
    }
    return(countList)
  }

 
See popCapture.R
##############################################################################
############Rare and locally common allele Capture in Seed accessions#########
##############################################################################

##Rare alleles and locally common alleles captured in seed accessions can be computed by the function:

myBeautifulFunc <- function(replicateShuffleDataset, referenceAlleles, referenceAllelesCount)

##Args: 
##replicateShuffleDataset is the group of shuffled seed accessions (999 replicates)
##referenceAlleles is the output of the function rareFunction or locallyCommonFunction computed on the reference dataset
##referenceAllelesCount is the output of the function countElementInRows computed on the output of rareFunction or
##locallyCommonFunction

##The function returns:
##“res”: number of alleles computed 5by5 rows (genotypes) up to the maximum number of genotypes in each permuted dataset,
##summed over the 999 permuted datasets (of replicateShuffleDataset) 
##“mean”: number of alleles computed 5by5 rows (genotypes) up to the maximum number of genotypes in each permuted dataset,
##averaged over the 999 permuted datasets (of replicateShuffleDataset)
##“meanInCount”: each element of the “mean” output divided by the number of alleles belonging to the required category in the
##reference dataset (Wild dataset in our case). This is the actual allelic capture.

##EXAMPLE, given the computations and functions described above:

SWyndcliff <- read.genalex("SWyndcliff.txt", ploidy = 2, geo = FALSE, region = FALSE, genclone = TRUE, sep = "\t", recode = FALSE)
ShuffleSWynd <- replicate(999, {shufflepop(SWyndcliff, method = 4)})

#rare alleles
WildRare <- rareFunction(PopWildFreq, 0.05)
CountRare <- countElementInRows(WildRare)

myBeautifulFunc(ShuffleSWynd, WildRare[[“Wyndcliff”]], CountRare[[“Wyndcliff”]])

$`res`
[1]  4698  7634  9264 10250 10825 10989

$mean
[1]  4.702703  7.641642  9.273273 10.260260 10.835836 11.000000

$meanInCount
[1] 0.2239382 0.3638877 0.4415844 0.4885838 0.5159922 0.5238095 
#locally common alleles

WildLCommon <- rareFunction(PopWildFreq, 0.05)
CountLCommon <- countElementInRows(WildRare)


myBeautifulFunc(ShuffleSWynd, WildLCommon[[“Wyndcliff”]], CountLCommon[[“Wyndcliff”]])

$`res`
[1]  4698  7634  9264 10250 10825 10989

$mean
[1]  4.702703  7.641642  9.273273 10.260260 10.835836 11.000000

$meanInCount
[1] 0.2239382 0.3638877 0.4415844 0.4885838 0.5159922 0.5238095




##CODE:

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
        #...get a subSet of rows...
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
