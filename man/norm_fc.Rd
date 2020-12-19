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
	Normalization steps:\cr
	1. collect and merge all .RCC files
	2. positive control normalization using \code{NanoStringNorm(..., CodeCount = "geo.mean")}
	3. flag samples with Positive Control normalization factor > 3 or < 0.3
	4. compute background threshold (mean + 2 * sd of negative controls)
	5. houskeeping gene normalization using \code{NanoStringNorm(..., SampleContent = "housekeeping.geo.mean")}
	6. flag samples with sampleContent normalization factor > 10 or < 0.1
	7. mask samples
}
\value{
	invisible NULL\cr
	Side effects: Creates folder 'results' in \code{input_folder}, containing:\cr
	* a table with normalized NanoString data\cr
	* a table with flagged samples (unusual normalization factors)\cr
 	* a table with ratios for target samples\cr
 	* a table with ratios for reference samples (mvp-samples)\cr
 	* a table with fold-changes\cr
 	* a table with fold-changes for reference samples (mvp-samples)\cr
 	* a table with ratios compared to all non-mvp-samples\cr
 	* a table with fold-changes compared all non-mvp-samples\cr
}
\references{
	NanoString Gene Expression Data Analysis Guidelines (MAN-C0011-04) \cr
	https://cran.r-project.org/web/packages/NanoStringNorm/vignettes/NanoStringNorm_Introduction.pdf \cr

}
\author{
	Steffi Grote
}
