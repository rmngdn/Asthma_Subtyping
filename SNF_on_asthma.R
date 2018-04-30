#SNF on asthma data
pathToSource <- "/homes/rag4218/Documents/Local_Work/Asthma_Subtyping/rawToMatrix.R"
source(pathToSource)
library(ggplot2)
library("CancerSubtypes")
library(SNFtool)
library(impute)
#Import the data:

##NOTE: CHANGE THIS TO YOUR PERSONNAL PATH ...

wd <- "/data/Datasets_Local/Asthma"
setwd(wd)
path <- "./dataFused.tsv"

#Take a little while, comment it after loading:

data <- read.delim(file = path, header = TRUE, sep = "\t") 

MansoorAssignment <- read.csv("./mansoor-assign.csv")

#Extract wanted informations:

reducedData <- data[data$PATIENT.ID %in% MansoorAssignment$subject_id,]


#Extract the matrix from each datatype:

matrixList <- rawToMatrix(reducedData)

transcriptomics <- matrixList[[1]]
proteomics <- matrixList[[2]]


print(mean(is.na(transcriptomics))) # no NA values in it!
print(mean(is.na(proteomics))) # ~8% of NA values in proteomics

#PROTEOMICS 

# NAdetect <- function(matrix) {
#   n <- dim(matrix)[1]
#   m <- dim(matrix)[2]
#   listRows <- matrix(0, 1, n)
#   listCols <- matrix(0, 1, m)
#   for (i in 1:n) {
#     for (j in 1:m) {
#       if (is.na(matrix[i, j])) {
#         listRows[i] <- 1
#         listCols[j] <- 1
#       }
#     }
#   }
#   res <- list(listRows, listCols)
#   return(res)
# }


# Detect the percentage of NA values per column or per row:

NAperPROBE <- data.frame(colMeans(is.na(proteomics)), row.names = colnames(proteomics))
colnames(NAperPROBE) <- "miss"
NAperPatient <- data.frame(rowMeans(is.na(proteomics)), row.names = rownames(proteomics))
colnames(NAperPatient) <- "miss"


ggplot(NAperPatient) + aes(x=NAperPatient[,1]) +geom_histogram()
ggplot(NAperPROBE) + aes(x=NAperPROBE[,1]) +geom_histogram()

#Test of SNF:
####WARNINGS####
# All the parameters here have been arbitrary chosen. IT'S JUST A TEST.

C = 3 #Number of clusters
K= 20 
alpha = 0.5
Tparam = 10

pMinMiss = 0.2


# We'll delete the probes that have more than pMinMiss of missingness:

listColumnsKept <- unlist(lapply(NAperPROBE, function(p) {p < pMinMiss}))
listRowsKept <- unlist(lapply(NAperPatient, function(p) {p < pMinMiss}))
filtered.proteomics <- proteomics[listRowsKept,listColumnsKept]

#Impute the missing data using knn imputation
imputed.proteomics <- impute.knn(filtered.proteomics)[[1]] # /!| default parameters /!| 

#Compute the affinity matrix
w.proteomics <- affinityMatrix(dist2(imputed.proteomics, imputed.proteomics))

#Assign the TACs identified previously to make a vector with ith element indicating the TAC number of the ith patient 
# (correspondingly to the matrix order)

groupAttribution <- function(matrix, assignment) {
  listPatient <- rownames(matrix)
  tac.patient.df <- assignment[,c("kuo_tac", "subject_id")]
  names(tac.patient.df)<- c("tac", "id")
  tac.patient.df$charID <- lapply(tac.patient.df$id, as.character)
  #print(head(tac.patient.df))
  A <- function(patientID) {
    #Extract
    TAC <- tac.patient.df[tac.patient.df[,3] == patientID, "tac"]
    #Transform to integer
    TAC <- as.character(TAC)
    TAC <- substr(TAC, 4, 4)
    TAC <- as.integer(TAC)
    return(TAC)
  }
  res <- unlist(lapply(listPatient, A))
  return(res)
}
group <- groupAttribution(filtered.proteomics, MansoorAssignment)

# Display the affinity matrix:

displayClusters(w.proteomics, group)

# TRANSCRIPTOMICS:

#we must have the same patients in the two data matrix so let's remove the same patients than proteomics:
filtered.transcriptomics <- transcriptomics[listRowsKept,listColumnsKept]
d.transcriptomics <- dist2(filtered.transcriptomics, reduced.transcriptomics)
w.transcriptomics <- affinityMatrix(d.transcriptomics)     
displayClusters(w.transcriptomics,group = group)
#this one is not so bad

# FUSION

w.fused <- SNF(list(w.proteomics, w.transcriptomics))#default parameters 
displayClusters(w.fused,group = group)

#not really convincing ...
