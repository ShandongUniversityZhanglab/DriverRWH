#' read_data
#' 
#' read_data uses reads the mutation file
#' it deletes some row(s) of the mutation data whose rowsum is zero, rownames of the remaining rows are considered possibale genes
#' 
#' @param mutation_file the text that contains the mutation data, a gene*sample matrix with rownames being gene names and column names being sample IDs
#' 
#' @return Mutation, the mutation data, rows represent genes and columns represent samples
#' @examples
#' 
#' ###read_data() reads and processes the mutation data
#' 
#' ###it deletes some row(s) of the mutation data whose rowsum is zero,
#' ###rownames of the remaining rows are considered possibale genes
#' 
#' library(DriverRWH)
#' 
#' ###place the mutation file in the current working directory and name it "mutation_file"
#' ###to run the example code, please uncomment them
#' #my_file_dir=getwd()
#' #mutation_file=paste0(my_file_dir,"/mutation_file.txt")
#' 
#' #Mutation=read_data(mutation_file=mutation_file)
#' 
#' ###the mutation data
#' #print(Mutation)
#' 
#' 
#' @export read_data
read_data <-
function(mutation_file){
  
  print(paste("reading mutation data from",mutation_file))
  data=read.delim(mutation_file,header = T,as.is = T,check.names = F)
  
  especial_gene_index=which(is.na(data[,1])|data[,1]=="")
  if (length(especial_gene_index)>0){
    data=data[-especial_gene_index,]
  } else {data=data }
  
  rownames(data)=data[,1]
  data=data[,-1]
  genes=rownames(data)

  Mutation=as.matrix(data[which(apply(data,1,sum)!=0),which(apply(data,2,sum)!=0)])
  
  return(Mutation)
}
