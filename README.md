

# NanoStringNormalizeR

Normalize all NanoString *.RCC* files in a directory using [NanoStringNorm](https://CRAN.R-project.org/package=NanoStringNorm)


## Installation

```
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("sgrote/NanoStringNormalizeR")
```


## Example workflow

##### 1. Load this package

```
library(NanoStringNormalizeR)
```

##### 2. Set a path to the directory that contains the *.RCC* files  

For example something like this on UNIX:

```
input_folder = "/home/user/NanoString/Raw_Data"
```

or something like this on Windows:

```
input_folder = "C:\Users\user\NanoString\Raw_Data"
```

##### 3. Define housekeeping genes

For example:

```
house_genes = c("ACTB", "GUSB", "MRPL19", "PSMC4", "PUM1", "RPLP0", "SF3A1", "TFRC")
```

##### 4. Normalize all *.RCC* files

This puts all the *.RCC* files in `input_folder` together and normalizes them at once

```
norm_fc(input_folder, house_genes)
```

_Note that also _.RCC_ in subdirectories of `input_folder` are taken into account_


## Output

This created a folder `results` in the `input_folder` containing the following files:

file | description |
----- | ----- |
`normalized_data.csv` | normalized NanoString data |
`flagged_samples.txt` | flagged samples (unusual normalization factors) |
`ratio.csv` | ratios for target samples |
`ratio_mvp_reference.csv` | ratios for reference samples (mvp-samples) |
`fold_change.csv` | fold-changes for target samples |
`fold_change_mvp_reference.csv`	| fold-changes for reference samples (mvp-samples) |


_Note that existing files will be overwritten_

