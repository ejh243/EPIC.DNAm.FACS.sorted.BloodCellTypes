# EPICBloodCompCalc
This package provides reference data from FACs sorted blood cell types and a function to use these data with the Houseman algorithm to estimate cellular proportion
The reference data were generated with the Illumina Infinium MethylationEPIC BeadChip and includes CD4 T cells, CD8 T cells, B cells, Monocytes and Granulocytes. Th output is a matrix with one column per cell type and one row per sample where the values represent the estimated cellular proportion for that sample and that cell type. 


### Install the package

The commands below will install the package directly from github. NB it assumes minfi is already installed. 

```{r,eval=FALSE}

install.packages("devtools")
library("devtools")
install_github("ejh243/EPICBloodCompCalc")
library(EPICBloodCompCalc)
library(minfi)

```

