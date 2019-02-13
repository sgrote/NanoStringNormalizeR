\name{norm_fc}
\alias{norm_fc}
\title{
Normalize all NanoString .RCC files in a directory
}
\description{
Wrapper for NanoStringNorm to be applied to all .RCC files in the input directory
}
\usage{
norm_fc(input_folder, house_genes)
}
\arguments{
  \item{input_folder}{
    character(), folder with .RCC files and/or subfolders with .RCC files for each run
}
  \item{house_genes}{
    character() vector with housekeeping genes to use for the normalization
}
}
\details{
    TODO: explain steps how NanoStringNorm is applied in detail
}
\value{
    invisible NULL
    Side effects: Creates folder 'results' in \code{input_folder}, containing:
    * a table with normalized NanoString data
    * a table with flagged samples (unusual normalization factors)
 	* a table with ratios for target samples
 	* a table with ratios for reference samples (mvp-samples)
 	* a table with fold-changes
 	* a table with fold-changes for reference samples (mvp-samples)
}
\references{
    TODO: NanoString Guidelines with version
}
\author{
    Steffi Grote
}
