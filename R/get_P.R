#' get_P
#' 
#' get_P uses the mutation data and the vertex weight of the hypergraph to calculate the transition possibility matrix
#' output of get_P is the transition possibility matrix of the hypergraph
#' 
#' @param Mutation the mutation data, rows represent genes and columns represent samples
#' @param vertex_W the vertex weight of the hypergraph
#' @examples
#' 
#' ###the vertex weight of the hypergraph is calculated by get_vertex_W(Mutation,graph)
#' ###and some relevant data should be calculated
#' 
#' ###the example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
#' ###We will show how to get the transition possibility matrix of the hypergraph
#' 
#' library(DriverRWH)
#' 
#' #load the mutation data
#' data(luscExampleMutation)
#' 
#' #load the network data
#' data(HumanNet)
#' 
#' #get the graph from HumanNet
#' graph=graph.data.frame(HumanNet)
#' 
#' #get the vertex weight of the hypergraph
#' W=get_vertex_W(Mutation=luscExampleMutation,graph=graph)
#' 
#' #get the transition possibility matrix P
#' P=get_P(Mutation=luscExampleMutation,vertex_W=W)
#' 
#' @export get_P
get_P <-
function(Mutation,vertex_W){
  
  Degree_v=apply(Mutation,1,sum)
  D_v_inverse=diag(1/Degree_v)
  
  part1=D_v_inverse%*%Mutation
  rownames(part1)=rownames(Mutation)
  
  Degree_ve=apply(vertex_W,2,sum)
  #Degree_ve=colSums(Mutation)
  D_ve_inverse=diag(1/Degree_ve)
  part2=D_ve_inverse%*%t(vertex_W)
  rownames(part2)=colnames(vertex_W)
  
  P=part1%*%part2
  return(P)
}
