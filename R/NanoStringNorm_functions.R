

# go through directory and subdirectories and merge all NanoString counts
collect_rccs = function(input_folder){
	message("\nCollecting .RCC files in ", input_folder ," ...")
	# get subdirectories    
	subdirs = list.dirs(path=input_folder, full.names=TRUE, recursive=FALSE)
	# add top-level directory to not require the data to be in subfolders
	subdirs = c(subdirs, input_folder)
	first = TRUE
	for (s in subdirs){
		# check for .RCC files
		n_rccs = length(grep("RCC$", dir(s)))
		if (n_rccs == 0){
			message("Skipping directory ", s," - no .RCC files found.")
			next
		}
		# read RCC
		rccs =  read.markup.RCC(s)
		# extract mRNA
		mrna = rccs$x
		# merge mRNA from different subdirectories
		if (first){
			merged_mrna = mrna
			first = FALSE
		} else {
			merged_mrna = merge(merged_mrna, mrna)
		}
	}
	if (first){
		stop(input_folder, " does not contain .RCC files or subfolders with .RCC files.")
	}
	# There are 3 gene-specific columns, that should be the same for all runs
	n_samples = ncol(merged_mrna) - 3
	message("\n\nFound ", n_samples, " samples")
	if (nrow(mrna) != nrow(merged_mrna)){
		stop("Merging of subfolders failed. Check that Class, Name and Accession for genes are consistent.")
	}
	colnames(merged_mrna) = remove_leading_X(colnames(merged_mrna))
	# just for consistency
	colnames(merged_mrna)[1] = "Code.Class"
	
	return(merged_mrna)
}

# define housekeeping genes explicitly
set_house_genes = function(mrna, house_genes){
	# specify housekeeping genes explicitly
	mrna[mrna[,1]=="Housekeeping", 1] = "Endogenous"
	mrna[mrna[,2] %in% house_genes, 1] = "Housekeeping"
	# sort by class and gene name
	mrna = mrna[order(mrna[,1], mrna[,2]),]
	return(mrna)
}

# remove leading X from char-vector (useful for sample names where X was added automatically)
remove_leading_X = function(text){
	xstart = substring(text,1,1) == "X"
	text[xstart] = substring(text[xstart], 2)
	return(text)
}

# get normalized data and flagged samples
# Positive control normalization -> mean+2sd-background mask (threshold) -> Housekeeping normalization
# this follows the order of normalization steps in NanoStringNorm
nano_norm = function(mrna){
	message("\nNormalize data...")
	
	## A) positive control normalization
	# CodeCount adjusts each sample based on its relative value to all samples
	normalized = NanoStringNorm(mrna, CodeCount = "geo.mean")
	# get flagged samples from positive control normalization
	norm_factors = normalized$sample.summary.stats.norm
	rownames(norm_factors) = remove_leading_X(rownames(norm_factors))
	flag_pos = norm_factors$pos.norm.factor > 3 | norm_factors$pos.norm.factor < 0.3
	flagged_pos = norm_factors[flag_pos,,drop=F]
	if(sum(flag_pos) > 0){
	  message("Samples with Positive Control normalization factor > 3 or < 0.3:")
	  print(flagged_pos)
	}
	
	## B) background thresholding
	# compute background threshold based on negative controls and create mask
	norm = normalized$normalized.data
	anno = norm[,1:3]
	counts = norm[,4:ncol(norm)]
	norm_bg = counts[anno[,1]=="Negative",,drop=F]
	bg_mean = apply(norm_bg, 2, function(x) {mean(x) + 2*sd(x)})
	mask = t(apply(counts, 1, function(x) {x < bg_mean}))

	## C) housekeeping gene normalization
	# SampleContent normalize for sample or RNA content i.e. pipetting fluctuations
	normalized = NanoStringNorm(norm, SampleContent = "housekeeping.geo.mean")
	norm = normalized$normalized.data
	colnames(norm) = remove_leading_X(colnames(norm))
	# get flagged samples from housekeeping normalization
	norm_factors = normalized$sample.summary.stats.norm
	rownames(norm_factors) = remove_leading_X(rownames(norm_factors))
	flag_house = norm_factors$sampleContent.norm.factor > 10 | norm_factors$sampleContent.norm.factor < 0.1
	flagged_house = norm_factors[flag_house,, drop=F]
	if(sum(flag_house) > 0){
	  message("Samples with sampleContent normalization factor > 10 or < 0.1:")
	  print(flagged_house)
	}

	# combine flagged samples and put back into order
	flagged = unique(c(rownames(flagged_pos), rownames(flagged_house)))
	flagged = flagged[order(match(flagged, colnames(norm)))]

	# mask samples that were below background threshold
	norm[,4:ncol(norm)][mask] = NA


	return(list(norm, flagged))
}


# function for geometric mean
geo_mean = function(data){
	log_data = log(data)
	gm = exp(mean(log_data[is.finite(log_data)]))
	return(gm)
}


# get ratios and fold-changes separately for target and reference (MPV) samples
get_fc = function(normalized_mrna){
	
	# get column indices for all reference samples
	ref_indices = grep("mvp", colnames(normalized_mrna), ignore.case=TRUE)
	if (length(ref_indices) == 0){
		stop("No reference samples found. Note that reference samples are assumed
			to carry 'mvp' or 'MVP' in their sample name.")
	}   
	ref_samples = colnames(normalized_mrna)[ref_indices]
	
	message("\nGet ratio and fold-change using reference-samples:\n",
		paste(ref_samples, collapse="\n"))

	# geometric mean across reference samples (vector, one entry per gene)
	ref_mean = apply(normalized_mrna[,ref_indices, drop=F], 1, geo_mean)
	ratio = (normalized_mrna[4:ncol(normalized_mrna)] / ref_mean)
	# fold-change is defined as (see Guidelines page 19)
	# If ratio > 1, then: Fold change = Ratio
	# If ratio < 1, then: Fold change = −1 ∗ 1/Ratio
	fc = ratio
	fc[fc < 1 & !is.na(fc)] = -1 / fc[fc < 1 & !is.na(fc)]
	
	# add sample-column to avoid row-names in output
	ratio = data.frame(Gene=row.names(ratio), ratio)
	fc = data.frame(Gene=row.names(fc), fc)
	
	# remove autoatically added leading X from sample names
	colnames(ratio) = remove_leading_X(colnames(ratio))
	colnames(fc) = remove_leading_X(colnames(fc))
	
	# separate reference from target samples
	ref_indices = grep("mvp", colnames(ratio), ignore.case=TRUE)
	ratio_ref = ratio[,c(1,ref_indices)]
	ratio = ratio[,-ref_indices]
	fc_ref = fc[,c(1,ref_indices)]
	fc = fc[,-ref_indices]
	
	return(list(ratio, fc, ratio_ref, fc_ref))
}


