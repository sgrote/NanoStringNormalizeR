\name{plot_ratios}
\alias{plot_ratios}
\title{
Plot ratios as heatmap and save as .png
}
\description{
Created png-file with heatmap of log2-ratios, like those created by \code{norm_fc}.
The png-file is in the same folder as the input, with "_heatmap.png" extension.
}
\usage{
plot_ratios(ratio_csv, title="normalized ratios", size_title=1, size_legend=0.8, size_row_names=0.7, size_col_names=0.8, mar_below=18, mar_right=6, png_width=20, png_height=24, dendrogram="both", Rowv=TRUE, Colv=TRUE)
}
\arguments{
  \item{ratio_csv}{
	character(), .csv file with ratios, like ratio.csv which is created by \code{norm_fc}, 
    or a subset of it
}
  \item{title}{
	optional character() for the plot title
}
  \item{size_title}{
	optional numeric() for the size of the plot title
}
  \item{size_legend}{
	optional numeric() for the size of the legend
}
  \item{size_row_names}{
	optional numeric() for the size of the row names
}
  \item{size_col_names}{
	optional numeric() for the size of the column names
}
  \item{mar_below}{
	optional numeric() for the size of the margin below the heatmap (for column names)
}
  \item{mar_right}{
	optional numeric() for the size of the margin on the right side of the heatmap (for row names)
}
  \item{png_width}{
	optional numeric() for the width of the .png output file
}
  \item{png_height}{
	optional numeric() for the height of the .png output file
}
  \item{dendrogram}{
	optional character() defining whether a dendrogram should be added. Possible values: "both" (default),"row","column","none"
}
  \item{Rowv}{
	optional logical() defining whether rows should be reorderd by mean value
}
  \item{Colv}{
	optional logical() defining whether columns should be reorderd by mean value
}

}

\value{
	invisible NULL
	Side effects:
    Creates png-file in the same folder as the \code{ratio_csv}, with "_heatmap.png" extension
}
\author{
	Steffi Grote
}
