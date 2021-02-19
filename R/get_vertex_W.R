#' get_vertex_W
#' 
#' get_vertex_W uses the mutation data and the graph constructed from network data to calculate the weight of the hypergraph
#' output of get_vertex_W is the vertex weight of the hypergraph
#' 
#' @param Mutation the mutation data, rows represent genes and columns represent samples
#' @param graph the graph constructed from network data
#' @importFrom igraph V
#' @importFrom igraph induced_subgraph
#' @importFrom igraph degree
#' @importFrom igraph graph.data.frame
#' @examples
#' 
#' ###the example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
#' ###the example graph is constructed from the HumanNet
#' ###We will show how to get the vertex weight of the hypergraph
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
#' #or from your own network
#' #network_file=""
#' #string=read.table(file=network_file,header=TRUE)
#' #graph=graph.data.frame(string)
#' 
#' #get the vertex weight of the hypergraph
#' W=get_vertex_W(Mutation=luscExampleMutation,graph=graph)
#' 
#' @export get_vertex_W
get_vertex_W <-
function(Mutation,graph){
  
  vertex_W=matrix(0,nrow=nrow(Mutation),ncol = ncol(Mutation))
  vertexes=rownames(Mutation)
  colnames(vertex_W)=colnames(Mutation)
  rownames(vertex_W)=rownames(Mutation)
  
  for(i in 1:ncol(Mutation)){

    v_i=intersect(vertexes[which(Mutation[,i]==1)],V(graph)$name)

    v_diff=setdiff(vertexes[which(Mutation[,i]==1)],V(graph)$name)
    graph_i=induced_subgraph(graph, v_i)
    
    graph_i_degree=degree(graph_i)[v_i]
    graph_i_degree[which(graph_i_degree==0)]=0.01
    
    vertex_W[v_i,i]=graph_i_degree[v_i]
    vertex_W[v_diff,i]=0.01
  }
  return(vertex_W)
  
}
