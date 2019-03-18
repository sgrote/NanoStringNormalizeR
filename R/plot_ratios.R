
#### Create heatmap of log2-ratios from .csv file with ratios

# ratio_csv = like 'ratio.csv' which is created in norm_fc()
# for subsets of samples/genes generate a ratio-file with only those genes/samples and use that as input 

#   * png-file with heatmap of log2-ratios
#   * the file is in the same folder as the input, with "_heatmap.png" extension


plot_ratios = function(ratio_csv, title="normalized ratios", size_title=1, size_legend=0.8, size_row_names=0.7, size_col_names=0.8, mar_below=18, mar_right=6, png_width=20, png_height=24, dendrogram="both", Rowv=TRUE, Colv=TRUE){
    
    ratio = read.csv(ratio_csv, check.names=F, header=T, as.is=T)
    
    # first column -> row-names
    row.names(ratio) = ratio[,1]
    ratio = ratio[,-1]
    # heatmap needs data as matrix
    ratio = as.matrix(ratio)
    # compute log2
    log2ratio = log2(ratio)
    # set missing values to 0 (TODO: handle differently?)
    log2ratio[is.na(log2ratio)] = 0
    # (0.2 + 1/log10(ncol(log2ratio))) is heatmap2 default but can get pretty large)
    cexRow = (0.2 + 1/log10(ncol(log2ratio)))*size_row_names
    cexCol = (0.2 + 1/log10(ncol(log2ratio)))*size_col_names # take only size_col_names sample names
    
    # create filename
    # remove ".csv" ending (just last 4 characters)
    outfile = substring(ratio_csv, 1, nchar(ratio_csv)-4)
    # add heatmap extension
    outfile = paste0(outfile, "_heatmap.png")
    
    # open png-file
    png(outfile, width=png_width, height=png_height, units='cm', res=400)
    # cex.main cannot be set inside heatmap
    par(cex.main=size_title) 
    # create heatmap
    heatmap.2(log2ratio, scale="none", col=greenred(100), main=title, density.info="none", trace="none",
              margins=c(mar_below, mar_right), keysize=size_legend, key.xlab="log2(ratio)", key.title="",
              cexCol=cexCol, cexRow=cexRow, dendrogram=dendrogram, Rowv=Rowv, Colv=Colv )
    # close png-file
    dev.off()

    message("Created ", outfile)
    
    return(invisible(NULL))
}
