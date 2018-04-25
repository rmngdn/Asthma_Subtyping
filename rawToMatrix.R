#rawToMatrix.R
#Romain GUEDON

#----PRESETS----

#The goal here is to create one matrix for each data type. The matrix is a classic one: one row per patient and features as columns

#Import the data:
wd <- "/data/Datasets_Local/Asthma"
path <-  paste0(wd,"/pagapow-DataExport-119281/subset1_UBIO_ADULT_FINAL/mrna.tsv")
data <- read.delim(file = path, header = TRUE, sep = "\t")


#----SUBFUNCTIONS---- 


#Add a column where PATIENT.ID takes an integer form, e.g. A_023 becomes 23
facToNum <- function(f) {
  return(as.integer(substring(as.character(f), 3, 5))) 
}

addColumn <- function(rawData) {
  rawData$intPID <- unlist(lapply(rawData$PATIENT.ID, facToNum)) 
}


#Now split data by TISSUE.TYPE into a list of data frames:

split_by_type <- function(data, data_column = "TISSUE.TYPE") {
  ###we let the choice of what column should be used to separate data type, by default, it's TISSUE.TYPE 
  # To avoid growing object we init the resList and then insert each element in the good place
  listOfType <- unique(data[,data_column])
  n <- length(listOfType)
  resList <- vector('list',n)
  for(i in 1:n) {
    temp <- data[data[,data_column] == listOfType[i],]
    resList[[i]] <- temp
  }
  return(resList)
}


split_by_type2 <- function(data,data_column = "TISSUE.TYPE") {
  ###we let the choice of what column should be used to separate data type, by default: it's TISSUE.TYPE 
  # We suppose that lapply is faster:
  listOfType <- unique(data[,data_column])
  A <- function(type) {
    return(data[data[,data_column] == type,])
  }
  resList <- lapply(listOfType, A)
  return(resList)
}


#typeList <- split_by_type2(data)

#We extract oneType 
#oneType <- typeList[[1]]

## The goal here is to create a matrix with named columns containing all the PROBE of the patients

matrixInit <- function(data,feature="PROBE") {
  ## return an empty (NA values) matrix of adjusted to the "size" of the data, i.e. patients in rows and features in columns
  listOfFeatures <- unique(data[,feature])
  listOfPatients <- unique(data[,"PATIENT.ID"])
  m <- length(listOfFeatures) 
  n <- length(listOfPatients)
  res <- matrix(NA, n, m)
  colnames(res) <- as.character(listOfFeatures)
  rownames(res) <- as.character(listOfPatients)
  return(res)
}

#testInit <- matrixInit(data,data$PROBE)
#testInit[1:3,1:3]

selectionValues <- function(oneType, onePatientID) {
  res <- oneType[oneType$PATIENT.ID == onePatientID, "VALUE"]
}
  
  
oneTypeToMatrix <- function(oneType, feature = "PROBE") {
  # sub-functions 
  valuesOnePatient <- function(onePatientID, feature = "PROBE") {
    #We extract the values of the features for the patient.
    temp <- oneType[oneType$PATIENT.ID == onePatientID, c("VALUE", feature)]
    res <- data.frame(temp[, "VALUE"], 
                      row.names = temp[, feature]) 
    colnames(res) <- onePatientID
    return(res)
  }
  
  # function body
  patientIDList <- unique(oneType$PATIENT.ID)
  resMatrix <- matrixInit(data, feature)
  
  
}       
oneTypeToMatrix(data)


#Let's find the duplicates

onePatient <- data[data$PATIENT.ID == "A_000", ]
listLength<-unlist(lapply(1:14, function(p) length(unique(onePatient[,p]))))
setNames(listLength, names(data))

findDuplicates <- function(data) {
  data2 <- data.frame(data)
  data2$geneID
}




#----MAIN FUNCTION----

rawToMatrix <- function(rawData, feature = "PROBE") {
  # add a column to facilitate the selection of patients
  addColumn(rawData) #maybe optional
  dataTypes <- split_by_type2(rawData)
  matrixList <- lapply(dataTypes, function(p) oneTypeToMatrix(p, feature))
  return(matrixList)
  }


testNames <- c("gene1", "gene2","gene3")
empty <- matrix(character(), 3, 1)
testDF <- data.frame(empty)
