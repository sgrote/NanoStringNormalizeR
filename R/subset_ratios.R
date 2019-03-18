
#### Filter ratio.csv file with gene list and save output
## usefule prior to plotting


subset_ratios = function(ratio_csv, genes_csv, outfile="ratio_subset.csv"){
    
    # avoid adding X to colnames starting with number
    ratio = read.csv(ratio_csv, check.names=F, header=T, as.is=T)
    genes = read.csv(genes_csv, header=F, as.is=T)[,1]
    ratio_subset = ratio[ratio[,1] %in% genes, ]
    
    not_in = genes[! genes %in% ratio[,1]]
    if (length(not_in) > 0){
        warning("Genes not in ratio_csv:\n", paste(not_in, collapse=", "))
    }

    # save to same folder as ratio_csv
    outfile = paste0(dirname(ratio_csv), "/", outfile)

    write.csv(ratio_subset, outfile, quote=F, row.names=F)

    message("Created ", outfile)
    
    return(invisible(NULL))
}
