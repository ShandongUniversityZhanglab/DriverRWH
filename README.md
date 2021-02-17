# DriverRWH
## Description
DriverRWH prioritizes altered genes by effective integration of somatic mutation data and molecular interaction data. A weighted Hypergraph model is formed by combining the degree of the PPI sub-network composed of each patient's mutated genes. Mutated genes are ranked by the random walk on the weighted Hypergraph model.
## Install
```r
#install necessary build tool
install.packages("devtools")
```
#### Install from github
```r
#the process may take several minutes

library(devtools)
install_github("ShandongUniversityZhanglab/DriverRWH")

#please restart R Console or Rstudio after installation
library(DriverRWH)
```
#### Download and install DriverRWH
```r
#please download zip of this repository
#after that, extract the files to a folder and change the working directory of R to that folder

library(devtools)
setwd("directory_of_folder")
build()
install()

#please restart R Console or Rstudio after installation
library(DriverRWH)
```
## Example
```r
###the example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
###the example network data is HumanNet
###We will show how to get DriverRWH Results

library(DriverRWH)

########## load the example mutation data ##########
#you can also use your own mutation data
#rows represent genes and columns represent samples
#rownames and colnames of mutation data cannot be omitted
data(luscExampleMutation)

##########     load the network data      ##########
#you can use HumanNet or STRINGv10 by loading HumanNet.rda or STRINGv10.rda
#you can also use your own network data
#the network data should have at least two columns
#the first two elements of each row are considered two vertices of an unweighted edge.
#Other elements in the row will be omitted
data(HumanNet)

#### set the current working directory as the output file directory ####
my_file_dir=getwd()

##########   get the result of DriverRWH   ##########
#the output file would be in the out_data folder in the output file directory
#the default value of theta is 0.2
get_DriverRWH_score(Mutation=luscExampleMutation,string=HumanNet,out_file_dir=my_file_dir)
```
## Contact
If you have any questions, please do not hesitate to contact us.
## Last update
Wednesday February 17，2021
