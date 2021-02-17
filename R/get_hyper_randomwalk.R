#' get_hyper_randomwalk
#' 
#' get_hyper_randomwalk uses the transition possibility matrix, the mutation data to perform random walk with restart on hypergraph
#' theta is the restart probability at every step of the random walk.
#' output of get_hyper_randomwalk is a list which combines the Importance of the genes and the distance of vectors between one iteration
#' 
#' @param P the transition possibility matrix
#' @param Mutation the mutation data, rows represent genes and columns represent samples
#' @param theta the restart probability at every step of the random walk
#' @examples
#' 
#' ###the transition possibility matrix P is calculated by get_P(Mutation,vertex_W) 
#' ###and some relevant data should be calculated
#' 
#' ###the example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
#' ###We will show how to get the result of random walk on hypergraph
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
#' ##########   get the graph from HumanNet  ##########
#' graph=graph.data.frame(HumanNet)
#' 
#' ####  get the vertex weight of the hypergraph   ####
#' W=get_vertex_W(Mutation=luscExampleMutation,graph=graph)
#' 
#' ####  get the transition possibility matrix P   ####
#' P=get_P(Mutation=luscExampleMutation,vertex_W=W)
#' 
#' #### get the result of randomwalk on hypergraph ####
#' #the default value of theta is 0.2
#' out=get_hyper_randomwalk(P=P,Mutation=luscExampleMutation)
#' 
#' @export get_hyper_randomwalk
get_hyper_randomwalk <-
function(P,Mutation,theta=0.2) {
  
  v0=rep(1/nrow(Mutation),nrow(Mutation))    #initial vector
  teleport=rep(1/nrow(Mutation),nrow(Mutation))  
  
  Distance=c()
  vi=v0
  for(k in 1:10){
    if(k==1) print("1st iteration")
    else if(k==2) print("2nd iteration")
    else print(paste0(k,"th iteration"))
    vj=vi
    vi=theta*t(P)%*%vi+(1-theta)*teleport
    dis=sum(abs(vj-vi))
    Distance=append(Distance,dis)
  }
  Importance=data.frame(vi[order(vi,decreasing = T),])
  colnames(Importance)=c("Importance")
  
  out=list(Importance,Distance)
  names(out)=c("Importance","Distance")
  return(out)
}
