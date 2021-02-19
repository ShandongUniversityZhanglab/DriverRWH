\name{DriverRWH-package}
\alias{DriverRWH-package}
\alias{DriverRWH}
\docType{package}
\title{
DriverRWH
}
\description{
DriverRWH prioritizes altered genes by effective integration of somatic mutation data and molecular interaction data. A weighted Hypergraph model is formed by combining the degree of the PPI sub-network composed of each patient's mutated genes. Mutated genes are ranked by the random walk on the weighted Hypergraph model.
}
\details{
\tabular{ll}{
Package: \tab DriverRWH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2020-10-08\cr
License: \tab MIT\cr
Authors: \tab Naiqian Zhang, Chengye Wang, Junhan Shi \cr
Maintainer: \tab Junhan Shi \cr
}
}
\keyword{ package }
\examples{
### the example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
### the example network data is HumanNet
### We will show how to get example DriverRWH Results
### to use your own mutation data or network data, please refer to get_DriverRWH_score
library(DriverRWH)

# load the mutation data
data(luscExampleMutation)

# load the network data
data(HumanNet)

# set the current working directory as the output file directory
my_file_dir=getwd()

# get the result of DriverRWH
# the output file would be in the out_data folder in the output file directory
# the default value of theta is 0.2
get_DriverRWH_score(Mutation=luscExampleMutation,string=HumanNet,out_file_dir=my_file_dir)

}
