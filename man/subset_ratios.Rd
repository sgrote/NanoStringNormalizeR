\name{subset_ratios}
\alias{subset_ratios}
\title{
Save subset of ratios
}
\description{
Save a subset of a .csv file with genes in the first column, like 'ratio.csv' or 'fold_change.csv' created by \code{norm_fc}, given a list of genes to filter for.  
The output will be written to the same directory as the input \code{ratio_csv} file.
}
\usage{
subset_ratios(ratio_csv, genes_csv, outfile="ratio_subset.csv")
}
\arguments{
  \item{ratio_csv}{
	character(), .csv file with genes in first column, like 'ratio.csv' which is created by \code{norm_fc}
}
  \item{genes_csv}{
	character(), .csv file with genes to filter for in first column
}
  \item{outfile}{
	optional character() for the output filename. The output will be written to the same directory as the input \code{ratio_csv} file.
}
}

\value{
	invisible NULL
	Side effects:
    Creates \code{outfile} in the same folder as the \code{ratio_csv}
}
\author{
	Steffi Grote
}
