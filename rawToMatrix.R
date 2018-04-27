#rawToMatrix.R
#Romain GUEDON

#----What is it for ?-----

# The goal here is to create one matrix per data type of the inputed dataset. 
# the matrix is in the classic format (rows = patients ; cols = features)



#----SUBFUNCTIONS---- 


#split data by TISSUE.TYPE into a list of data frames:

split_by_type2 <- function(rawData, data_column = "TISSUE.TYPE") {
  # we let the choice of what column should be used to separate data type, by default: it's TISSUE.TYPE 
  listOfType <- unique(rawData[, data_column])
  namesOfType <- lapply(listOfType, as.character)
  A <- function(type) {
    return(rawData[rawData[,data_column] == type,])
  }
  resList <- lapply(listOfType, A)
  return(setNames(resList, namesOfType))
}


# extract the data from a given patient and transform it into 
# a named vector
valuesPatient <- function(patientID, type, feature = "PROBE") {
  # extract the values of the features for the patient
  temp <- type[type$PATIENT.ID == patientID, c("VALUE", feature)]
  value <- t(matrix(temp[, "VALUE"]))
  rownames(value) <- patientID
  colnames(value) <- temp[, feature]
  return(value)
}

# creates a matrix with the good size with NA values. The matrix is ready to receive values from each patient
matrixInit <- function(type,feature="PROBE") {
  listOfFeatures <- unique(type[,feature])
  listOfPatients <- unique(type[,"PATIENT.ID"])
  m <- length(listOfFeatures) 
  n <- length(listOfPatients)
  res <- matrix(NA, n, m)
  colnames(res) <- lapply(listOfFeatures, as.character)
  rownames(res) <- lapply(listOfPatients, as.character)
  return(res)
}


typeToMatrix <- function(type, feature = "PROBE") {
  # extract patient IDs
  patientIDList <- unique(type$PATIENT.ID)
  # initiate the matrix
  resMatrix <- matrixInit(type, feature)
  # extract data from each patient and insert it at the right place
  # in the matrix
  test <- TRUE
  for (patientID in patientIDList) {
    row <- valuesPatient(patientID, type, feature)
    resMatrix[rownames(row),colnames(row)] <- row
  }
  return(resMatrix)
}       



#----MAIN FUNCTION----

rawToMatrix <- function(rawData, feature = "PROBE") {
  dataTypes <- split_by_type2(rawData)
  matrixList <- lapply(dataTypes, function(p) typeToMatrix(p, feature))
  return(matrixList)
  }

#----TESTS---- 

#wd <- "/data/Datasets_Local/Asthma"
#path1 <-  paste0(wd,"/mrna1.tsv")
#data1 <- read.delim(file = path1, header = TRUE, sep = "\t")



