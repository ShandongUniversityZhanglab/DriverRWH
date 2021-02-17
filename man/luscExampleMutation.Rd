\name{luscExampleMutation}
\alias{luscExampleMutation}
\docType{data}
\title{
luscExampleMutation
}
\description{
Mutation data for a Lung Squamous Cell Carcinoma Dataset from The Cancer Gene Atlas (TCGA). The data is presented in matrix with genes being the rows and patients being columns. A binary 0 or 1 is used to indicate whether or not a gene is mutated in that patient.
}
\usage{data(luscExampleMutation)}
\details{
Mutation data for an example Lung Squamous Cell Carcinoma Dataset. Rows are genes, columns are patients
}
\source{
https://xenabrowser.net/datapages/?dataset=mc3_gene_level%2FLUSC_mc3_gene_level.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222
}
\references{
The Cancer Genome Atlas. Somatic mutation (SNP and INDEL) - MC3 gene-level non-silent mutation
Nature. 490 (7418):61-70.
Cell Systems. Volume 6 Issue 3: p271-281.e7, 28 March 2018 10.1016/j.cels.2018.03.002
}
\author{
Kyle Ellrott, Matthew H. Bailey, Gordon Saksena, Li Ding, et. al.
}
\examples{
data(luscExampleMutation)
}
\keyword{datasets}
