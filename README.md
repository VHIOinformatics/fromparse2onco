# fromparse2onco

Functions to create oncoplots from the results of parseVCF. There is an all in one function to do it in a straight way.

## Installation

```r
library(devtools)
install_github("VHIOinformatics/fromparse2onco")
```


## Package contents

* **fromParse2MAF.R**: The function `fromParse2MAF` reads variant data from one or more Excel files produced by parseVCF program and processes it in order to match MAF format specifications.

* **filterMAF.R**: The function `filterMAF` filters a MAF data frame of genetic variants based on specified criteria (VAF, reads & flags).

* **prepareForOncoplot.R**: The function `prepareForOncoplot` takes a MAF dataframe, filters it and creates an oncomatrix. It also prints summary plots and creates a TMB table. It does so by executing functions in maftools package.

* **makeOncoplot.R** The function `makeOncoplot` imports an onco_matrix.txt file in the working directory and makes an oncoplot using the oncoPrint function in ComplexHeatmap package.

* **fromParse2Onco.R** The function `fromParse2Onco`reads one or more Excel files produced by parseVCF program and makes and oncoplot. It does so by sequentially executing the previous four functions.

For more information on the parameters and usage of the functions, please check the documentation in R.
