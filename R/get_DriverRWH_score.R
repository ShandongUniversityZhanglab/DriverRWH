#' get_DriverRWH_score
#' 
#' get_DriverRWH_score uses the mutation data, network data to perform DriverRWH.
#' theta is the restart probability at every step of the random walk.
#' result of the DriverRWH will be in out_file_dir
#' 
#' @param Mutation the mutation data, rows represent genes and columns represent samples
#' @param string network data, it should have at least two columns and the first two elements of each row are considered two vertices of an unweighted edge. Other elements in the row will be omitted
#' @param theta the restart probability at every step of the random walk
#' @param out_file_dir the directory where you wish to put result of the DriverRWH
#' @importFrom  igraph graph.data.frame
#' @importFrom utils write.table
#' @examples
#' 
#' ###the example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
#' ###the example network data is HumanNet
#' ###We will show how to get DriverRWH Results
#' 
#' library(DriverRWH)
#' 
#' ########## load the example mutation data ##########
#' #you can also use your own mutation data
#' #rows represent genes and columns represent samples
#' #rownames and colnames of mutation data cannot be omitted
#' data(luscExampleMutation)
#' 
#' ##########     load the network data      ##########
#' #you can use HumanNet or STRINGv10 by loading HumanNet.rda or STRINGv10.rda
#' #you can also use your own network data
#' #the network data should have at least two columns
#' #the first two elements of each row are considered two vertices of an unweighted edge.
#' #Other elements in the row will be omitted
#' data(HumanNet)
#' 
#' #### set the current working directory as the output file directory ####
#' my_file_dir=getwd()
#' 
#' ##########   get the result of DriverRWH   ##########
#' #the output file would be in the out_data folder in the output file directory
#' #the default value of theta is 0.2
#' get_DriverRWH_score(Mutation=luscExampleMutation,string=HumanNet,out_file_dir=my_file_dir)
#' 
#' @export get_DriverRWH_score
get_DriverRWH_score <-
function(Mutation,string,theta=0.2,out_file_dir){
  
  dir.create(paste(out_file_dir,"/out_data/",sep = ""))
  
  graph=graph.data.frame(string)
  
  print("caculating weight of vertex in each hyperedge based on PPI network")
  vertex_W=get_vertex_W(Mutation,graph)
  
  print("caculating transform matrix P")
  P=get_P(Mutation,vertex_W)
  
  print("random walking in hypergraph")
  out=get_hyper_randomwalk(P,Mutation,theta)
  
  Importance=out[[1]]
  Distance=out[[2]]
  
  print(paste("writing Importance to",out_file_dir))
  write.table( Importance,file=paste(out_file_dir,"/out_data/Importance.txt",sep = ""))
  print(paste("writing Distance to",out_file_dir))
  write.table( Distance,file=paste(out_file_dir,"/out_data/Distance.txt",sep = ""))
  
}
