# DriverRWH
## Description
DriverRWH prioritizes altered genes by effective integration of somatic mutation data and molecular interaction data. A weighted Hypergraph model is formed by combining the degree of the PPI sub-network composed of each patient's mutated genes. Mutated genes are ranked by the random walk on the weighted Hypergraph model.
## Install
```r
# Install necessary build tool
install.packages("devtools")
```
#### Install from github
```r
# The process below may take several minutes

library(devtools)
install_github("ShandongUniversityZhanglab/DriverRWH")

# Please restart R Console or Rstudio after installation
library(DriverRWH)
```
#### Download and install DriverRWH
```r
# Please download zip code of this repository
# After that, extract the files to a folder and change the working directory of R to that folder

library(devtools)

directory_of_folder <- "" # root directory of DriverRWH package
setwd("directory_of_folder")
build()
install()

# Please restart R Console or Rstudio after installation
library(DriverRWH)
```
## Example
```r
### The example mutation data is a Lung Squamous Cell Carcinoma Dataset from TCGA
### The example network data is HumanNet
### We will show how to get DriverRWH Results

library(DriverRWH)

########## Load the example mutation data ##########
# You can also use your own mutation data
# Rows represent genes and columns represent samples
# Rownames(name of genes) and colnames(name of samples) of mutation data cannot be omitted
data(luscExampleMutation)

##########     Load the network data      ##########
# You can use HumanNet or STRINGv10 by loading HumanNet.rda or STRINGv10.rda
# You can also use your own network data
# The first two elements of each row are considered two vertices of an unweighted edge.
# Other elements in the row will be omitted
data(HumanNet)

#### Set the current working directory as the output file directory ####
my_file_dir=getwd()

##########   Get the result of DriverRWH   ##########
# The output file would be in the out_data folder in the output file directory
# The default value of theta is 0.2
get_DriverRWH_score(Mutation=luscExampleMutation,string=HumanNet,out_file_dir=my_file_dir)
```
## Contact
If you have any questions, please do not hesitate to contact us.
## Last update
Friday February 19，2021
