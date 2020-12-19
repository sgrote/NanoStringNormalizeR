norm_fc = function(input_folder, house_genes){
	
	# collect all .RCC files and define housekeeping genes
	merged_mrna = collect_rccs(input_folder)
	merged_mrna = set_house_genes(merged_mrna, house_genes)
	
	# normalize and flag
	normalized_and_flagged = nano_norm(merged_mrna)
	normalized_mrna = normalized_and_flagged[[1]]
	flagged_samples = normalized_and_flagged[[2]]

	# get ratio and foldchange
	ratio_and_fc = get_fc(normalized_mrna)
	ratio = ratio_and_fc[[1]]
	fold_change = ratio_and_fc[[2]]
	ratio_ref = ratio_and_fc[[3]]
	fold_change_ref = ratio_and_fc[[4]]
	ratio_wo_ref = ratio_and_fc[[5]]
	fold_change_wo_ref = ratio_and_fc[[6]]
	
	# write output
	output_folder = paste0(input_folder, "/results")
	suppressWarnings(dir.create(output_folder)) 
	message("\nWrite output-files to ", output_folder)
	
	write.csv(merged_mrna, paste0(output_folder, "/raw_data.csv"), quote=F, row.names=F)
	write.csv(normalized_mrna, paste0(output_folder, "/normalized_data.csv"), quote=F, row.names=F)
	write.csv(ratio, paste0(output_folder, "/ratio.csv"), quote=F, row.names=F)
	write.csv(fold_change, paste0(output_folder, "/fold_change.csv"), quote=F, row.names=F)
	write.csv(ratio_ref, paste0(output_folder, "/ratio_mvp_reference.csv"), quote=F, row.names=F)
	write.csv(fold_change_ref, paste0(output_folder, "/fold_change_mvp_reference.csv"), quote=F, row.names=F)
	write.csv(ratio_wo_ref, paste0(output_folder, "/ratio_without_mvp_reference.csv"), quote=F, row.names=F)
	write.csv(fold_change_wo_ref, paste0(output_folder, "/fold_change_without_mvp_reference.csv"), quote=F, row.names=F)
	write.table(flagged_samples, paste0(output_folder, "/flagged_samples.txt"), quote=F, row.names=F, col.names=F)
	message("\nDone.\n")
	
	return(invisible(NULL))
}
