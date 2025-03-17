# fromParse2onco

Functions to create oncoplots from the results of parse variants. There is an all in one function to do it in a straight way.

## Installation


## Package contents

* **fromparse2table.R**: The funcion `fromparse2table` creates a dataframe of the variants and change the configuration for the library of maftools to create the maf file. It is possible to use more than one excel and tumor only or paired analysis.

* **filteringtable.R**: The function `filteringtable` filters the data frame based on criteria of VAF, reads, flags o the existence of pahtogenicity in CGI or oncokb.

* **modclassvar.R**: The function `modclassvar` classifies the variants in the groups requested by maftools and generate maf objects like the matrix. The matrix is written in the same folder and the object returned is the Tumor Mutational Burden to use it afterwards.

* **importandfilterbysamples.R**: The function `importandfilterbysamples` imports the matrix and filter the genes out by a minimum number of samples mutated.

* **confoplot.R** The function `confoplot` configures the colors for each type of variant and define the shape for building the oncoplot.

* **fromZero2Plot.R** All in one function to create an oncoplot from the parse variants objects

For more information on the parameters and usage of the functions, please check the documentation in R.
