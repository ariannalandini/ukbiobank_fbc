
## author: William John Astle
## author email: wja24@cam.ac.uk/will@astle.net
## Description: Look at main FBC by various covariates and produce plots

## TODO
## script should be split into functions for each plot type

#rsq_vec1=int_tech_r.sq
#n_vec=int_n_vector_sorted

generate_rsq_plot<-function(rsq_vec1, rsq_vec2=NULL, n_vec, output_dir=UKBIOBANK_FBC_OUTPUT_DIR, study_name="UK Biobank", type="technical")
{
	phenotype_order=c("plt","plt_f","plt_i","mpv","pdw","pct","rbc","mcv","hct","mch","mchc","hgb","rdw","rdw_sd_raw","ret","ret_p","irf","hlr","hlr_p","mono","neut","eo","baso","neut_eo_sum","eo_baso_sum","baso_neut_sum","gran","neut_p_gran","baso_p_gran","eo_p_gran","myeloid_wbc","gran_p_myeloid_wbc","lymph","wbc","mono_p","neut_p","eo_p","baso_p","lymph_p")
	rsq_data=data.frame(rsq_vec1=rsq_vec1)
	if(hasArg(rsq_vec2))
	{
		rsq_data$rsq_vec2=rsq_vec2
	}
	if(hasArg(n_vec))
	{
		rsq_data$n_vec=n_vec
	}
	print(dim(rsq_data))
	pull_out_elements=is.element(row.names(rsq_data), phenotype_order)
	rsq_data=rsq_data[pull_out_elements,, drop=FALSE]
	
	element_order=order(match(row.names(rsq_data),phenotype_order))
	rsq_data=rsq_data[element_order,, drop=FALSE]

	factor_names=get_trait_short_name(row.names(rsq_data))
	rsq_data$names=factor(factor_names, levels=factor_names)

	if(type=="technical")
	{
		rsq_data$n_effective=rsq_data$rsq_vec1/(1-rsq_data$rsq_vec1)*rsq_data$n_vec
	}
	print("warning - writing temporary data")
	write.table(rsq_data,"~/tmp.tsv", sep="\t", quote=FALSE)
	rsq_plot<-ggplot(rsq_data, aes(x=names, y=100*rsq_vec1))+geom_bar(stat="identity")+theme_bw(base_size=24)+ylab("Variance Explained on Adjustment Scale (%)")+xlab("Variable Name")+theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.border = element_blank(),axis.line.x = element_line(color = 'black'), axis.line.y=element_line(color = 'black'))

	if(type=="technical")
	{ rsq_plot<-rsq_plot+stat_identity(geom="text", aes(x=names, label=as.character(floor(n_effective))),size=10, vjust=-1)#+labs(title=sprintf("%s: Variance Explained by Seasonal and Technical Effects\n(Counts Show Effective Number of Additional Observations)",study_name))

	}
	if(type=="gwas")
	{
		rsq_plot<-rsq_plot#+labs(title=sprintf("%s: Variance of Variable pre-adjusted for Technical Effects Explained by Environmental Covariates and Sex",study_name))	
	}
	rsq_plot<-rsq_plot+scale_y_continuous(limits=c(0,50))
	return(rsq_plot)
}


generate_rsq_plots_for_paper<-function()
{
	dir.create(sprintf("%s/graphics/gwas_var_explained",INTERVAL_FBC_OUTPUT_DIR), showWarnings = FALSE)
	

	int_n_vector=apply(!is.na(int_extract[,grep("gwas_normalised",colnames(int_extract), value=TRUE)]), 2, sum)
	names(int_n_vector)=gsub("(.+)_gwas_normalised","\\1", names(int_n_vector))
	int_n_vector=c(rep(as.list(int_n_vector)$plt,2), NA, int_n_vector)
	names(int_n_vector)[1:3]=c("plt_f","plt_i", "rdw_sd_raw")
	int_n_vector_sorted=int_n_vector[match(names(int_tech_r.sq), names(int_n_vector))]

	pdf(sprintf("%s/graphics/gwas_var_explained/INTERVAL_tech_var_explained.pdf",INTERVAL_FBC_OUTPUT_DIR), width=24, height=8)
	print(generate_rsq_plot(int_tech_r.sq,n_vec=int_n_vector_sorted, study_name="INTERVAL")+labs(title=""))	
	dev.off()

	pdf(sprintf("%s/graphics/gwas_var_explained/INTERVAL_gwas_var_explained.pdf", INTERVAL_FBC_OUTPUT_DIR), width=14, height=10)
	print(generate_rsq_plot(int_gwas_r.sq, study_name="INTERVAL", type="gwas"))
	dev.off()

	
	dir.create(sprintf("%s/graphics/gwas_var_explained",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)
	
	
	ukbb_n_vector=apply(!is.na(bd_extract[,grep("gwas_normalised",colnames(bd_extract), value=TRUE)]), 2, sum)
	names(ukbb_n_vector)=gsub("(.+)_gwas_normalised","\\1", names(ukbb_n_vector))
	ukbb_n_vector=c(NA, ukbb_n_vector)
	names(ukbb_n_vector)[1]=c("rdw_sd_raw")
	ukbb_n_vector_sorted=ukbb_n_vector[match(names(ukbb_tech_r.sq), names(ukbb_n_vector))]
	
	pdf(sprintf("%s/graphics/gwas_var_explained/UKBB_tech_var_explained.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=24,height=8)
	print(generate_rsq_plot(ukbb_tech_r.sq,n_vec=ukbb_n_vector_sorted, study_name="UK Biobank")+labs(title=""))	
	dev.off()

	pdf(sprintf("%s/graphics/gwas_var_explained/UKBB_gwas_var_explained.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=14, height=10)
	print(generate_rsq_plot(ukbb_gwas_r.sq, study_name="UK Biobank", type="gwas"))
	dev.off()
		
}



#expit<-function(x)
#{	
#	return(exp(x)/(1+exp(x)))
#}
power_plot<-function()
{
	N=173039
	C=qchisq(log(8.31)-9*log(10),df=1, log.p=TRUE, lower.tail=FALSE)
	
	n_maf=1000
	n_effect=1000
	power_frame=data.frame(logit_maf=rep(9.5*((1:n_maf)/n_maf-1),n_effect), log_effect=0.69-6+6*ceiling((1:(n_maf*n_effect)/n_maf))/n_effect)
	power_frame$maf=expit(power_frame$logit_maf)
	power_frame$effect=exp(power_frame$log_effect)
	power_frame$rsq=2*power_frame$maf*(1-power_frame$maf)*power_frame$effect^2
	power_frame$power=1-pchisq(C,nc=N*power_frame$rsq,df=1)

	power_plot<-ggplot(power_frame, aes(x=logit_maf, y=log_effect))+geom_raster(aes(fill=power))+theme_bw()
	power_plot<-power_plot+scale_y_continuous(limits=c(log(0.01), log(2)),breaks=c(log(2), log(1), log(0.1), log(0.01), log(0.001)), labels=c("2","1","0.1","0.01", "0.001"))
	power_plot<-power_plot+scale_x_continuous(limits=c(logit(0.0001),0),breaks=c(logit(0.5) , logit(0.4),logit(0.3), logit(0.2), logit(0.1), logit(0.01), logit(0.001), logit(0.0001)), labels=c("0.5","0.4","0.3","0.2", "0.1", "0.01", "0.001", "0.0001"))+scale_fill_gradientn("Power", colours=c("#FFFFFF",rev(rainbow(5))), limits=c(0,1))
	
	
	power_plot<-power_plot+xlab("Minor Allele Frequency")+ylab("True Effect (Standard Deviations)")+theme_bw(base_size=24)	
	
	pdf("~/test_plot.pdf", width=14)
	print(power_plot)
	dev.off()
}

blue_print_TF_plot<-function()
{
	expression_data=read.table("/scratch/wja24/ukbb_tmp/blueprint_platelet_data/blueprint_platelet_gwas_paper.tsv", stringsAsFactors=FALSE, head=TRUE)
	expression_data$cell=factor(expression_data$cell, levels=rev(c("MK","EB","Mono","Neut", "NK","CD4T","CD8T","CD19B")))
	expression_data$gene_name=factor(expression_data$gene_name, levels=c("ENSG00000179588","ENSG00000169946","ENSG00000102145","ENSG00000179348"))
	expression_plot<-ggplot(expression_data, aes(x=gene_name,y=cell, fill=expression))+geom_tile(color = "white")+ scale_fill_gradient2(low="darkblue",mid= "yellow1", high = "red",  name="Expression", breaks=c(-6, 0 ,4), labels=c("None", "Very Low", "High"))+theme_minimal(base_size=14)+theme(axis.text.x = element_text(angle = 45, vjust = 1,   hjust = 1))+ coord_fixed()+ylab("Cell Type")+xlab("Gene")+scale_y_discrete(labels=c("CD19 B", "CD8 T", "CD4 T","Natural Killer","Neutrophil", "Monocyte","Erythroblast", "Megakaryocyte"))+scale_x_discrete(labels=c("ZFPM1","ZFPM2","GATA1","GATA2"))#+labs(title="Pairwise Correlations Between Quantile Normalised GWAS Phenotypes")

	pdf("~/WHO_TF_ZFPM_GATA.pdf")
	print(expression_plot)
	dev.off()
}


power_curve<-function(x)
{
	C=qchisq(log(8.31)-9*log(10),df=1, log.p=TRUE, lower.tail=FALSE)
	return(sqrt(C/(2.0*x*(1-x)*N)))
}

generate_effect_size_plot_for_paper<-function()
{

	library(data.table)
	OUTPUT_DIR="/scratch/wja24/ukbb_tmp/results_tmp"
	
	clump_lead_wide_table=fread(sprintf("%s/meta_analysis_UKBB_UKBiLEVE_INTERVAL_sorted_trait_wide_summary/all_traits_clump_leads.tsv", OUTPUT_DIR))
	clump_lead_wide_table$abs_effect=abs(clump_lead_wide_table$Effect_clump_lead)
	
	accurate_af=fread("/scratch/wja24/ukbb_tmp/results_tmp/clump_leads_af/alt_af.tsv")
	clump_lead_wide_table$MAF=accurate_af$AF[match(clump_lead_wide_table$SNP, accurate_af$SNP)]
	clump_lead_wide_table$MAF[clump_lead_wide_table$MAF>=0.5]=1-clump_lead_wide_table$MAF[clump_lead_wide_table$MAF>=0.5]
#	clump_lead_wide_table$
	clump_lead_wide_table$clump_lead_trait_type="Compound white cell"
	clump_lead_wide_table$clump_lead_trait_type[is.element(clump_lead_wide_table$minP_trait, ukbb_platelet_traits)]="Platelet"
	clump_lead_wide_table$clump_lead_trait_type[is.element(clump_lead_wide_table$minP_trait, ukbb_myeloid_wbc_traits)]="Myeloid white cell"
	clump_lead_wide_table$clump_lead_trait_type[is.element(clump_lead_wide_table$minP_trait, "lymph")]="Lymphoid white cell"
	clump_lead_wide_table$clump_lead_trait_type[is.element(clump_lead_wide_table$minP_trait, ukbb_red_cell_traits)]="Red cell"
	clump_lead_wide_table$clump_lead_trait_type=factor(clump_lead_wide_table$clump_lead_trait_type, levels=c("Platelet", "Red cell", "Myeloid white cell", "Lymphoid white cell", "Compound white cell"))
	ukbb_n_vector=apply(!is.na(bd_extract[,grep("gwas_normalised",colnames(bd_extract), value=TRUE)]), 2, sum)
	int_n_vector=apply(!is.na(int_extract[,grep("gwas_normalised",colnames(int_extract), value=TRUE)]), 2, sum)
	common=intersect(names(ukbb_n_vector), names(int_n_vector))
	N<<-max(ukbb_n_vector[match(common, names(ukbb_n_vector))]+int_n_vector[match(common, names(int_n_vector))])
	
	clump_lead_wide_table$logit_MAF=logit(clump_lead_wide_table$MAF)
	clump_lead_wide_table$log_abs_effect=log(clump_lead_wide_table$abs_effect)


	effect_plot<-ggplot(data=clump_lead_wide_table, aes(x=logit_MAF, y=log_abs_effect))+theme_bw(base_size=24)+geom_point(aes(colour=clump_lead_trait_type), alpha=0.9)+ylab("Estimated Effect (Standard Deviations)")+xlab("Minor Allele Frequency")+labs(colour="Trait type with min p-value")+scale_color_manual(values=c("mediumorchid3","firebrick3", "orange","steelblue","lightskyblue2"))#+scale_y_continuous(limits=c(0,2))#+stat_function(fun=power_curve, col="grey22",size=0.8, n=50000)
	
	effect_plot<-effect_plot+scale_y_continuous(limits=c(log(0.01), log(2)),breaks=c(log(2), log(1), log(0.1), log(0.01), log(0.001)), labels=c("2","1","0.1","0.01", "0.001"))
	effect_plot<-effect_plot+scale_x_continuous(limits=c(logit(0.0001),0),breaks=c(logit(0.5) , logit(0.4), logit(0.3), logit(0.2), logit(0.1), logit(0.01), logit(0.001), logit(0.0001)), labels=c("0.5","0.4","0.3","0.2", "0.1", "0.01", "0.001", "0.0001"))
	
	dir.create(sprintf("%s/graphics/hits",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)	
	pdf(sprintf("%s/graphics/hits/effect_sizes.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=16.5)
	print(effect_plot)
	dev.off()

	density_effect_plot<-ggplot(data=clump_lead_wide_table, aes(x=logit_MAF, y=log_abs_effect))+theme_bw(base_size=24)+stat_density2d(geom="raster", aes(fill = ..density..), contour=FALSE)+scale_fill_continuous(low="white", high="blue")+ylab("Estimated Effect (Standard Deviations)")+xlab("Minor Allele Frequency")+labs(colour="Trait type with min p-value")
	
	density_effect_plot<-density_effect_plot+scale_y_continuous(limits=c(log(0.01), log(2)),breaks=c(log(2), log(1), log(0.1), log(0.01), log(0.001)), labels=c("2","1","0.1","0.01", "0.001"))
	density_effect_plot<-density_effect_plot+scale_x_continuous(limits=c(logit(0.0001),0),breaks=c(logit(0.5) , logit(0.4), logit(0.3), logit(0.2), logit(0.1), logit(0.01), logit(0.001), logit(0.0001)), labels=c("0.5","0.4","0.3","0.2", "0.1", "0.01", "0.001", "0.0001"))
	
	#effect_plot_small=effect_plot+scale_x_continuous(limits=c(-0.001,0.1))+ylab("")+xlab("")+guides(colour=FALSE)
	#pdf(sprintf("%s/graphics/hits/effect_sizes_small.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=7)
	#print(effect_plot_small)
	#dev.off()

	pdf(sprintf("%s/graphics/hits/effect_sizes_density.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=16.5)
	print(density_effect_plot)
	dev.off()

	
	return(1);

	OUTPUT_DIR="/scratch/wja24/ukbb_tmp/results_tmp"

	full_rbc_summary_data=fread(sprintf("%s/meta_analysis_UKBB_UKBiLEVE_INTERVAL_sorted/meta_rbc_gwas_normalised_sorted.tsv", OUTPUT_DIR))

	hist_plot_data=data.frame(SNP=full_rbc_summary_data$SNP, MAF=full_rbc_summary_data$Freq_EAl)
	hist_plot_data$MAF[hist_plot_data$MAF>=0.5]=1-hist_plot_data$MAF[hist_plot_data$MAF>=0.5]
	hist_plot_data$lead=is.element(hist_plot_data$SNP, clump_lead_wide_table$SNP)

	lead_counts=hist(hist_plot_data$MAF[hist_plot_data$lead], breaks=0:25/50,plot=FALSE)$counts
	hist_counts=hist(hist_plot_data$MAF, breaks=0:25/50,plot=FALSE)$counts
	ratio_data=data.frame(MAF=1:25/51, ratio=lead_counts/hist_counts)


	pdf(sprintf("%s/graphics/hits/MAF_hit_ratio.pdf",UKBIOBANK_FBC_OUTPUT_DIR))
	maf_ratio_plot<-ggplot(data=ratio_data, aes(x=MAF, y=100*ratio))+geom_point()+theme_bw(base_size=24)+scale_x_continuous(limits=c(0,0.5))+xlab("Minor Allele Frequency")+ylab("Sentinel Fraction of Tested Variants (%)")
	print(maf_ratio_plot)
	dev.off()
	
	pdf(sprintf("%s/graphics/hits/MAF_distribution.pdf",UKBIOBANK_FBC_OUTPUT_DIR))
	maf_plot<-ggplot(data=hist_plot_data, aes(x=MAF))+geom_histogram(data=subset(hist_plot_data,lead), aes(y=..density..), binwidth=0.005)+theme_bw(base_size=24)+scale_x_continuous(limits=c(0,0.5))+xlab("Minor Allele Frequency")+ylab("Probability Density")
	print(maf_plot)
	dev.off()

		
		
}



generate_heritability_plot_for_paper<-function()
{
	ukbb=read.table("/lustre4/wja24/data/UK_biobank/extracted_heritability_estimates/UKBiobank.tsv", head=FALSE, stringsAsFactors=FALSE)
	ukbileve=read.table("/lustre4/wja24/data/UK_biobank/extracted_heritability_estimates/UKBiLEVE.tsv", head=FALSE, stringsAsFactors=FALSE)
	interval=read.table("/lustre4/wja24/data/UK_biobank/extracted_heritability_estimates/INTERVAL.tsv", head=FALSE, stringsAsFactors=FALSE)

	ukbb_to_merge=cbind(ukbb, rep("UK Biobank (not UK BiLEVE)", dim(ukbb)[1]))
	colnames(ukbb_to_merge)=c("Phenotype","h2","study")
	
	ukbil_to_merge=cbind(ukbileve, rep("UK BiLEVE", dim(ukbileve)[1]))
	colnames(ukbil_to_merge)=c("Phenotype","h2","study")
	
	interval_to_merge=cbind(interval, rep("INTERVAL", dim(interval)[1]))
	colnames(interval_to_merge)=c("Phenotype","h2","study")

	merged_data=rbind(ukbb_to_merge, ukbil_to_merge, interval_to_merge)


	phenotype_order=c("plt","plt_f","plt_i","mpv","pdw","pct","rbc","mcv","hct","mch","mchc","hgb","rdw","rdw_sd_raw","ret","ret_p","irf","hlr","hlr_p","mono","neut","eo","baso","neut_eo_sum","eo_baso_sum","baso_neut_sum","gran","neut_p_gran","baso_p_gran","eo_p_gran","myeloid_wbc","gran_p_myeloid_wbc","lymph","wbc","mono_p","neut_p","eo_p","baso_p","lymph_p")
	
	pull_out_elements=is.element(merged_data$Phenotype, phenotype_order)
	merged_data=merged_data[pull_out_elements,, drop=FALSE]

	factor_level_raw_names=unique(merged_data$Phenotype)
	factor_level_raw_names=factor_level_raw_names[order(match(factor_level_raw_names, phenotype_order))]

	factor_levels=get_trait_short_name(factor_level_raw_names)
	merged_data$names=factor(get_trait_short_name(merged_data$Phenotype), levels=factor_levels)
		


	herit_dot_plot<-ggplot(data=merged_data, aes(x=names, y=100*h2))+theme_bw(base_size=24)+geom_point(aes(color=factor(study), name="Sub-study"),  position=position_jitter(h=0, w=0.13), size=4)+scale_y_continuous(limits=c(0,30))+xlab("Phenotype")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Variance Explained (%)")+labs(colour="Component Study")#labs(colour="Sub-study", title="Genotype Estimated Trait Heritabilities")

	dir.create(sprintf("%s/graphics/heritability",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)
	pdf(sprintf("%s/graphics/heritability/gwas_herit.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=21)
	print(herit_dot_plot)
	dev.off()
	
	
	merged_data$type=get_trait_type(merged_data[,1])
	herit_bar_plot_data=merged_data %>% group_by(type) %>% summarise(h2_median=100*median(h2))
	herit_bar_plot_data$factor_type=factor(herit_bar_plot_data$type, levels=c("Platelet", "Mature red cell", "Immature red cell", "Myeloid white cell", "Lymphoid white cell", "Compound white cell"))

	herit_bar_plot<-ggplot(herit_bar_plot_data, aes(x=factor_type, y=h2_median, fill=factor_type))+theme_bw(base_size=24)+geom_bar(stat="identity")+scale_y_continuous(limits=c(0,30))+xlab("Trait Type")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Median Variance Explained (%)")+labs(colour="Component Study")+scale_fill_manual(values=c("mediumorchid3","firebrick3","firebrick2", "orange","steelblue","lightskyblue2")) + theme(legend.position="none")
	
	pdf(sprintf("%s/graphics/heritability/gwas_short_var_explained.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=10)
	print(herit_bar_plot)
	dev.off()
	
	
	
}

dace_plots_for_paper<-function()
{

	library(data.table)
	library(grid)

	dace_data_reordCoverage<-fread("/scratch/wja24/ukbb_tmp/results_tmp/dace/reordCoverage.txt")
	dace_data_reordTable<-fread("/scratch/wja24/ukbb_tmp/results_tmp/dace/reordTable.txt")

	dace_data_reordTable$Regulatory.element[dace_data_reordTable$Regulatory.element=="enhancers"]="Enhancer"
	dace_data_reordTable$Regulatory.element[dace_data_reordTable$Regulatory.element=="promoters"]="Promoter"
	dace_data_reordTable$Regulatory.element[dace_data_reordTable$Regulatory.element=="repressors"]="Silencer"

	dace_data_reordCoverage$elementAnnot[dace_data_reordCoverage$elementAnnot=="enhancer"]="Enhancer"
	dace_data_reordCoverage$elementAnnot[dace_data_reordCoverage$elementAnnot=="promoter"]="Promoter"
	dace_data_reordCoverage$elementAnnot[dace_data_reordCoverage$elementAnnot=="repressor"]="Silencer"


	dace_data_reordCoverage$match=sprintf("%s:%s",dace_data_reordCoverage$cellAnnot, dace_data_reordCoverage$elementAnnot)
	dace_data_reordTable$match=sprintf("%s:%s",dace_data_reordTable$Cell.type, dace_data_reordTable$Regulatory.element)
	dace_data_reordCoverage$SNP.ratio=dace_data_reordTable$SNP.ratio[match(dace_data_reordCoverage$match, dace_data_reordTable$match)]

	dace_data_reordCoverage$cell_factor=factor(dace_data_reordCoverage$cellAnnot, levels=c("MK","EB","MONO","NEUT", "EO"))
	dace_data_reordCoverage$element_factor=factor(dace_data_reordCoverage$elementAnnot, levels=c("Promoter", "Enhancer", "Silencer"))




	my_plot<-ggplot(dace_data_reordCoverage)+geom_point(aes(x=100*proportions, y=100*SNP.ratio, colour=cell_factor), size=5)+facet_wrap(~element_factor, scales="free")+theme_bw(base_size=22)+theme(legend.key.size=unit(1.5, "cm"))+xlab("Proportion of Genome Covered by Functional Element (%)")+ylab("Proportion of Sentinel Variants (%)")+scale_color_manual("Cell Type", labels=c("Megakaryocytes","Erythroblasts", "Monocytes", "Neutrophils","Eosinophils"), values=c("mediumorchid3","firebrick2", "yellowgreen","orange2","orange4"))
	pdf("~/test_plot.pdf",width=12)
	print(my_plot)
	dev.off()


}


adam_plots_for_paper<-function()
{
	library(data.table)
	OUTPUT_DIR="/scratch/wja24/ukbb_tmp/results_tmp"
	
	clump_lead_wide_table=fread(sprintf("%s/meta_analysis_UKBB_UKBiLEVE_INTERVAL_sorted_trait_wide_summary/all_traits_clump_leads.tsv", OUTPUT_DIR))
	clump_lead_wide_table$abs_effect=abs(clump_lead_wide_table$Effect_clump_lead)
		
}

generate_qc_example_plots_for_paper<-function()
{
	load_gwas_adjusted_UKBB_data()
	unadjusted_wbc_delay<-mean_plot_by_covariate(bd_extract, "wbc", "delay_in_hours_to_quarter_trim", instrument_id="AJ38695", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")
	adjusted_wbc_delay<-mean_plot_by_covariate(bd_extract, "wbc_tech_adj", "delay_in_hours_to_quarter_trim", instrument_id="AJ38695", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")
	unadjusted_wbc_delay<-unadjusted_wbc_delay+scale_y_continuous(limits=c(5,10))
	adjusted_wbc_delay<-adjusted_wbc_delay+scale_y_continuous(limits=c(5,10))

	unadjusted_mono_p_time<-mean_plot_by_covariate(bd_extract, "mono_p", "quarter_acq_in_hours", instrument_id="AJ38695", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")
	adjusted_mono_p_time<-mean_plot_by_covariate(bd_extract, "mono_p_tech_adj", "quarter_acq_in_hours", instrument_id="AJ38695", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")

	unadjusted_mono_p_time<-unadjusted_mono_p_time+scale_x_continuous(limits=c(7,24))+scale_y_continuous(limits=c(1,14))
	adjusted_mono_p_time<-adjusted_mono_p_time+scale_x_continuous(limits=c(7,24))+scale_y_continuous(limits=c(1,14))

	unadjusted_neut_day_of_study<-mean_plot_by_covariate(bd_extract, "neut", "day_of_study_acq", instrument_id="AK30431", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")
	adjusted_neut_day_of_study<-mean_plot_by_covariate(bd_extract, "neut_tech_adj", "day_of_study_acq", instrument_id="AK30431", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")



	unadjusted_neut_day_of_study<-unadjusted_neut_day_of_study+scale_y_continuous(limits=c(0,7))
	adjusted_neut_day_of_study<-adjusted_neut_day_of_study+scale_y_continuous(limits=c(0,7))

	unadjusted_mcv_day_of_study<-mean_plot_by_covariate(bd_extract, "mcv", "day_of_study_acq", instrument_id="AJ38695", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")
	adjusted_mcv_day_of_study<-mean_plot_by_covariate(bd_extract, "mcv_tech_adj", "day_of_study_acq", instrument_id="AJ38695", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")



	unadjusted_mcv_day_of_study<-unadjusted_mcv_day_of_study+scale_y_continuous(limits=c(80,100))
	adjusted_mcv_day_of_study<-adjusted_mcv_day_of_study+scale_y_continuous(limits=c(80,100))

	unadjusted_neut_by_sex_mean_plot<-mean_plot_by_covariate(bd_extract[is.element(bd_extract$meno, c("pre", "post", "male")),], "neut_tech_adj", "age_to_year_acq", "meno", study_name="UK Biobank")+labs(title="UK Biobank: Mean of NEUT#")+theme_bw(base_size=24)+labs(title="")
	# probably should convert this to correct scale.
	print("Warning commiting terrible sin of rewriting raw data")
	bd_extract$neut_gwas_adj=exp(bd_extract$neut_gwas_adj)
	adjusted_neut_by_sex_mean_plot<-mean_plot_by_covariate(bd_extract[is.element(bd_extract$meno, c("pre", "post", "male")),], "neut_gwas_adj", "age_to_year_acq", "meno", study_name="UK Biobank")+theme_bw(base_size=24)+labs(title="")
	unadjusted_neut_by_sex_mean_plot<-unadjusted_neut_by_sex_mean_plot+scale_x_continuous(limits=c(40,70))+scale_y_continuous(limits=c(3,5.5))
	adjusted_neut_by_sex_mean_plot<-adjusted_neut_by_sex_mean_plot+scale_x_continuous(limits=c(40,70))+scale_y_continuous(limits=c(3,5.5))
	
	dir.create(sprintf("%s/graphics/gwas_adjust_medley",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)
	pdf(sprintf("%s/graphics/gwas_adjust_medley/gwas_adjust_medley.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=5*5)
	multiplot(plotlist=list(unadjusted_mcv_day_of_study, adjusted_mcv_day_of_study,unadjusted_neut_day_of_study, adjusted_neut_day_of_study,unadjusted_mono_p_time, adjusted_mono_p_time, unadjusted_wbc_delay, adjusted_wbc_delay,unadjusted_neut_by_sex_mean_plot, adjusted_neut_by_sex_mean_plot), cols=5)
	dev.off()
	
	pdf(sprintf("%s/graphics/gwas_adjust_medley/wbc_unadjusted_delay.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=10, height=7)
	print(unadjusted_wbc_delay)
	dev.off()

	pdf(sprintf("%s/graphics/gwas_adjust_medley/wbc_adjusted_delay.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=10, height=7)
	print(adjusted_wbc_delay)
	dev.off()
	
	pdf(sprintf("%s/graphics/gwas_adjust_medley/mono_p_unadjusted_time.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=8, height=7)
	print(unadjusted_mono_p_time)
	dev.off()

	pdf(sprintf("%s/graphics/gwas_adjust_medley/mono_p_adjusted_time.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=8, height=7)
	print(adjusted_mono_p_time)
	dev.off()
	
	pdf(sprintf("%s/graphics/gwas_adjust_medley/neut_unadjusted_day_of_study.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=8, height=7)
	print(unadjusted_neut_day_of_study)
	dev.off()

	pdf(sprintf("%s/graphics/gwas_adjust_medley/neut_adjusted_day_of_study.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=8, height=7)
	print(adjusted_neut_day_of_study)
	dev.off()
	
	pdf(sprintf("%s/graphics/gwas_adjust_medley/mcv_unadjusted_day_of_study.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=8, height=7)
	print(unadjusted_mcv_day_of_study)
	dev.off()

	pdf(sprintf("%s/graphics/gwas_adjust_medley/mcv_adjusted_day_of_study.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=8, height=7)
	print(adjusted_mcv_day_of_study)
	dev.off()

	pdf(sprintf("%s/graphics/gwas_adjust_medley/neut_unadjusted_by_sex.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=14, height=7)
	print(unadjusted_neut_by_sex_mean_plot)
	dev.off()

	pdf(sprintf("%s/graphics/gwas_adjust_medley/neut_adjusted_by_sex.pdf",UKBIOBANK_FBC_OUTPUT_DIR), width=14, height=7)
	print(adjusted_neut_by_sex_mean_plot)
	dev.off()

}



generate_qc_plots<-function(data_local,data_filters, trait="plt",covariate_list, output_dir=getwd(), study_name="UK Biobank")
{

	inst_ids=unique(data_local$instrument)
	inst_ids=inst_ids[!is.na(inst_ids)]
	
	time_series_plot=list()

	if(!hasArg(covariate_list))
	{	
		covariate_list=c("week_of_study_acq", "quarter_acq_in_hours", "day_of_week_acq", "day_of_week_ext", "day_of_study_acq", "month_of_year_acq", "week_of_year_acq","delay_in_hours_to_quarter","delay_in_hours_to_quarter_trim")
	}

	missing_vars=covariate_list[!is.element(covariate_list, colnames(data_local))]
	if(length(missing_vars)>0)
	{
		print("The variables below are missing from this dataset and will not have plots drawn")
		print(missing_vars)
		covariate_list=covariate_list[is.element(covariate_list, colnames(data_local))]
	}
	for(covariate in covariate_list)
	{
		eval(parse(text=sprintf("%s_summary_plot=list()", covariate)))
		eval(parse(text=sprintf("%s_iq_plot=list()", covariate)))
		eval(parse(text=sprintf("%s_mean_plot=list()", covariate)))
		eval(parse(text=sprintf("%s_sd_plot=list()", covariate)))
		#eval(parse(text=sprintf("%s_log_var_vs_log_mean_plot=list()", covariate)))

	}


	day_cum_plot=list()
	violin_plot=list()

	time_series_ylab=ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))

#	 simple dot plot by time of aquisition
	for(instrument_id in inst_ids)
	{
		eval(parse(text=sprintf("time_series_plot[[instrument_id]]=ggplot(data=subset(data_local, instrument==instrument_id), aes(x=as.numeric(second_of_study_acq/60/60), y=%s))", trait)))

		time_series_plot[[instrument_id]]=time_series_plot[[instrument_id]]+theme_bw()+time_series_ylab+scale_x_continuous(limits=c(0,max(as.numeric(data_local$second_of_study_acq))/60/60))+xlab("Time Between Start of Study and Acquisition (Hours)")+labs(title=sprintf("%s:%s\n Instrument ID: %s",study_name, toupper(trait), instrument_id))+ geom_point(size=0.5) 

		for(covariate in covariate_list)
		{

			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)

			eval(parse(text=sprintf("%s_summary_plot[[instrument_id]]<-summary_plot_by_covariate(data_local, trait, \"%s\", instrument_id=instrument_id, summary_table=summary_table, study_name=study_name)", covariate, covariate)))
			eval(parse(text=sprintf("%s_iq_plot[[instrument_id]]<-iq_plot_by_covariate(data_local, trait, \"%s\",instrument_id=instrument_id, summary_table=summary_table, study_name=study_name)", covariate, covariate)))
			eval(parse(text=sprintf("%s_mean_plot[[instrument_id]]<-mean_plot_by_covariate(data_local, trait, \"%s\", instrument_id=instrument_id, summary_table=summary_table, study_name=study_name)", covariate, covariate)))
			eval(parse(text=sprintf("%s_sd_plot[[instrument_id]]<-sd_plot_by_covariate(data_local, trait, \"%s\",instrument_id=instrument_id, summary_table=summary_table, study_name=study_name)", covariate, covariate)))
		#	eval(parse(text=sprintf("%s_log_var_vs_log_mean_plot[[instrument_id]]<-log_var_vs_log_mean_plot_by_covariate(data_local, trait, \"%s\", instrument_id=instrument_id, summary_table=summary_table, study_name=study_name)", covariate, covariate)))

		}		
	
	}

	dir.create(sprintf("%s/graphics",output_dir, covariate), showWarnings = FALSE)
	dir.create(sprintf("%s/graphics/time_series",output_dir, covariate), showWarnings = FALSE)
	dir.create(sprintf("%s/graphics/day_of_study_cummu_drift",output_dir), showWarnings = FALSE)
	dir.create(sprintf("%s/graphics/time_series",output_dir), showWarnings = FALSE)

	print("still need to draw cumm sum plots")
	for(covariate in covariate_list)
	{
		dir.create(sprintf("%s/graphics/%s_summary",output_dir, covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_iq",output_dir,covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_mean",output_dir,covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_sd",output_dir, covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_log_var_vs_log_mean",output_dir, covariate), showWarnings = FALSE)
			
		pdf(sprintf("%s/graphics/%s_summary/%s_%s_summary_by_instrument.pdf",output_dir,covariate,toupper(trait), covariate), height=7*length(inst_ids))
		eval(parse(text = sprintf("multiplot(plotlist=%s_summary_plot)", covariate)))
		dev.off()
		
		pdf(sprintf("%s/graphics/%s_iq/%s_%s_iq_by_instrument.pdf",output_dir,covariate,toupper(trait), covariate), height=7*length(inst_ids))

		eval(parse(text = sprintf("multiplot(plotlist=%s_iq_plot)", covariate)))
		dev.off()

		pdf(sprintf("%s/graphics/%s_mean/%s_%s_mean_by_instrument.pdf",output_dir,covariate,toupper(trait), covariate), height=7*length(inst_ids))

		eval(parse(text = sprintf("multiplot(plotlist=%s_mean_plot)", covariate)))
		dev.off()

		pdf(sprintf("%s/graphics/%s_sd/%s_%s_sd_by_instrument.pdf",output_dir,covariate,toupper(trait), covariate), height=7*length(inst_ids))

		eval(parse(text = sprintf("multiplot(plotlist=%s_sd_plot)", covariate)))
		dev.off()
				
	#	pdf(sprintf("%s/graphics/%s_log_var_vs_log_mean/%s_%s_log_var_vs_log_mean_by_instrument.pdf",output_dir,covariate,toupper(trait), covariate), height=7*length(inst_ids))

	#	eval(parse(text = sprintf("multiplot(plotlist=%s_log_var_vs_log_mean_plot)", covariate)))
	#	dev.off()
	}
	

	
	png(sprintf("%s/graphics/time_series/%s_time_series_by_instrument.png",output_dir,toupper(trait)), width=2*1280, height=2*2560, res=8*72)	
	multiplot(plotlist=time_series_plot)
	dev.off()
}


generate_covariate_plots<-function(data_local, data_filters, trait="plt", covariate_list, output_dir=getwd(), study_name="UK Biobank")
{

	if(!hasArg(covariate_list))
	{	
		covariate_list=c("age_to_year_acq")
	}
	
	missing_vars=covariate_list[!is.element(covariate_list, colnames(data_local))]
	if(length(missing_vars)>0)
	{
		print("The variables below are missing from this dataset and will not have plots drawn")
		print(missing_vars)
		covariate_list=covariate_list[is.element(covariate_list, colnames(data_local))]
	}

	dir.create(sprintf("%s/graphics",output_dir), showWarnings = FALSE)
	

	for(covariate in covariate_list)
	{	print("here")

		summary_table=summarise_trait_by_covariate(data_local, trait, covariate, "sex")	
			print("there")


		eval(parse(text=sprintf("%s_by_sex_summary_plot<-summary_plot_by_covariate(data_local, trait, \"%s\", \"sex\", summary_table=summary_table,study_name=study_name)", covariate, covariate)))
		eval(parse(text=sprintf("%s_by_sex_iq_plot<-iq_plot_by_covariate(data_local, trait, \"%s\", \"sex\", summary_table=summary_table, study_name=study_name)", covariate, covariate)))
		eval(parse(text=sprintf("%s_by_sex_mean_plot<-mean_plot_by_covariate(data_local, trait, \"%s\", \"sex\", summary_table=summary_table,study_name=study_name)", covariate, covariate)))
		eval(parse(text=sprintf("%s_by_sex_sd_plot<-sd_plot_by_covariate(data_local, trait, \"%s\",  \"sex\", summary_table=summary_table, study_name=study_name)", covariate, covariate)))
		summary_table=summarise_trait_by_covariate(data_local, trait, covariate, "meno_simple")	
		eval(parse(text=sprintf("%s_by_meno_simple_summary_plot<-summary_plot_by_covariate(data_local, trait, \"%s\", \"meno_simple\",  summary_table=summary_table,study_name=study_name)", covariate, covariate)))
		eval(parse(text=sprintf("%s_by_meno_simple_iq_plot<-iq_plot_by_covariate(data_local, trait, \"%s\",  \"meno_simple\",  summary_table=summary_table, study_name=study_name)", covariate, covariate)))
		eval(parse(text=sprintf("%s_by_meno_simple_mean_plot<-mean_plot_by_covariate(data_local, trait, \"%s\",  \"meno_simple\",  summary_table=summary_table,study_name=study_name)", covariate, covariate)))
		eval(parse(text=sprintf("%s_by_meno_simple_sd_plot<-sd_plot_by_covariate(data_local, trait, \"%s\",   \"meno_simple\", summary_table=summary_table,study_name=study_name)", covariate, covariate)))
	
	}
	for(covariate in covariate_list)
	{
		print("got here")
		dir.create(sprintf("%s/graphics/%s_summary_by_sex",output_dir, covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_iq_by_sex",output_dir,covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_mean_by_sex",output_dir,covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_sd_by_sex",output_dir, covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_summary_by_meno_simple",output_dir, covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_iq_by_meno_simple",output_dir,covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_mean_by_meno_simple",output_dir,covariate), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/%s_sd_by_meno_simple",output_dir, covariate), showWarnings = FALSE)

	
		pdf(sprintf("%s/graphics/%s_summary_by_sex/%s_%s_summary_by_sex.pdf",output_dir,covariate,toupper(trait), covariate), height=7)
		eval(parse(text = sprintf("print(%s_by_sex_summary_plot)", covariate)))
		dev.off()
		pdf(sprintf("%s/graphics/%s_iq_by_sex/%s_%s_iq_by_sex.pdf",output_dir,covariate,toupper(trait), covariate), height=7)

		eval(parse(text = sprintf("print(%s_by_sex_iq_plot)", covariate)))
		dev.off()
	
		pdf(sprintf("%s/graphics/%s_mean_by_sex/%s_%s_mean_by_sex.pdf",output_dir,covariate,toupper(trait), covariate), height=7)

		eval(parse(text = sprintf("print(%s_by_sex_mean_plot)", covariate)))
		dev.off()
		
		pdf(sprintf("%s/graphics/%s_sd_by_sex/%s_%s_sd_by_sex.pdf",output_dir,covariate,toupper(trait), covariate), height=7)

		eval(parse(text = sprintf("print(%s_by_sex_sd_plot)", covariate)))
		dev.off()
			
#	comment2<-function()
#	{
			
		pdf(sprintf("%s/graphics/%s_summary_by_meno_simple/%s_%s_summary_by_meno_simple.pdf",output_dir,covariate,toupper(trait), covariate), height=7)
		eval(parse(text = sprintf("print(%s_by_meno_simple_summary_plot)", covariate)))
		dev.off()

		pdf(sprintf("%s/graphics/%s_iq_by_meno_simple/%s_%s_iq_by_meno_simple.pdf",output_dir,covariate,toupper(trait), covariate), height=7)

		eval(parse(text = sprintf("print(%s_by_meno_simple_iq_plot)", covariate)))
		dev.off()
	
		pdf(sprintf("%s/graphics/%s_mean_by_meno_simple/%s_%s_mean_by_meno_simple.pdf",output_dir,covariate,toupper(trait), covariate), height=7)

		eval(parse(text = sprintf("print(%s_by_meno_simple_mean_plot)", covariate)))
		dev.off()
		pdf(sprintf("%s/graphics/%s_sd_by_meno_simple/%s_%s_sd_by_meno_simple.pdf",output_dir,covariate,toupper(trait), covariate), height=7)

		eval(parse(text = sprintf("print(%s_by_meno_simple_sd_plot)", covariate)))
		dev.off()
				
#	}	
	
	}

}


generate_histograms<-function(data_local, data_filters, trait="plt", output_dir=UKBIOBANK_FBC_OUTPUT_DIR, study_name="UK Biobank")
{

	dir.create(sprintf("%s/graphics",output_dir), showWarnings = FALSE)
	
	histogram_by_sex<-hist_plot_by_binary_covariate(data_local, trait, "sex", study_name)
	histogram_by_meno_simple<-hist_plot_by_catagorical(data_local, trait, "meno_simple", study_name)


	dir.create(sprintf("%s/graphics/histogram_by_sex",output_dir), showWarnings = FALSE)
	dir.create(sprintf("%s/graphics/histogram_by_meno_simple",output_dir), showWarnings = FALSE)
	
	pdf(sprintf("%s/graphics/histogram_by_sex/%s_histogram_by_sex.pdf",output_dir,toupper(trait)), height=7)
	print(histogram_by_sex)
	dev.off()

	pdf(sprintf("%s/graphics/histogram_by_meno_simple/%s_histogram_by_meno_simple.pdf",output_dir,toupper(trait)), height=7)
	print(histogram_by_meno_simple)
	dev.off()

	
}




summarise_trait_by_covariate<-function(data_local,trait, covariate, strata, instrument_id)
{

	if(!hasArg(instrument_id)) {data_subset=rep(TRUE, dim(data_local)[1])} else {data_subset=data_local$instrument==instrument_id}
	data_subset[is.na(data_subset)]=FALSE
	if(!hasArg(strata))
	{
		print("1")
		eval(parse(text=sprintf("subset_data=na.omit(data_local[data_subset,c(\"%s\",\"%s\")])",trait, covariate)))

		eval(parse(text=sprintf("summary_table<-subset_data %%>%% group_by(%s) %%>%% summarise(mean=mean(%s,na.rm=TRUE), ci_mean95upper=mean(%s,na.rm=TRUE)+1.96*sd(%s,na.rm=TRUE)/sqrt(length(%s)),  ci_mean95lower=mean(%s,na.rm=TRUE)-1.96*sd(%s,na.rm=TRUE)/sqrt(length(%s)),n=length(%s),st_dev=sd(%s,na.rm=TRUE),  var=var(%s,na.rm=TRUE), ci_var95lower=var(%s,na.rm=TRUE)*(length(%s)-1)/qchisq(0.05/2, length(%s)-1), ci_var95upper=var(%s,na.rm=TRUE)*(length(%s)-1)/qchisq(1-0.05/2, length(%s)-1))", covariate, trait,trait,trait,trait,trait,trait, trait,trait,trait,trait,trait,trait,trait,trait,trait,trait,trait,trait)))
		for (quant in 0:20/20)
		{
       	     		eval(parse(text = sprintf("quantile_table<-subset_data %%>%% group_by(%s) %%>%% summarise(quantile_%0.2f=quantile(%s, %f, na.rm=TRUE, names=FALSE))",  covariate, quant, trait, quant)))
       	     		eval(parse(text = sprintf("summary_table=cbind(summary_table, quantile_%0.2f=quantile_table$quantile_%0.2f)", quant, quant)))
       		}

	
	}
	else
	{
#	print("1")

		eval(parse(text=sprintf("subset_data=na.omit(data_local[data_subset,c(\"%s\",\"%s\",\"%s\")])",trait, covariate, strata)))
#	print("2")


	eval(parse(text=sprintf("summary_table<-subset_data %%>%% group_by(%s, %s) %%>%% summarise(mean=mean(%s,na.rm=TRUE), ci_mean95upper=mean(%s,na.rm=TRUE)+1.96*sd(%s,na.rm=TRUE)/sqrt(length(%s)),  ci_mean95lower=mean(%s,na.rm=TRUE)-1.96*sd(%s,na.rm=TRUE)/sqrt(length(%s)),n=length(%s),st_dev=sd(%s,na.rm=TRUE),  var=var(%s,na.rm=TRUE), ci_var95lower=var(%s,na.rm=TRUE)*(length(%s)-1)/qchisq(0.05/2, length(%s)-1), ci_var95upper=var(%s,na.rm=TRUE)*(length(%s)-1)/qchisq(1-0.05/2, length(%s)-1))", covariate, strata, trait,trait,trait,trait,trait,trait, trait,trait,trait,trait,trait,trait,trait,trait,trait,trait,trait,trait)))
	#	print("3")

	for (quant in 0:20/20)
		{
			#	print("4")

       	     		eval(parse(text = sprintf("quantile_table<-subset_data %%>%% group_by(%s,%s) %%>%% summarise(quantile_%0.2f=quantile(%s, %f, na.rm=TRUE,names=FALSE))",  covariate, strata, quant, trait, quant)))
       	     	#	print("5")
		#	print(summary_table)
			eval(parse(text = sprintf("summary_table$quantile_%0.2f=quantile_table$quantile_%0.2f", quant, quant)))
       		}


	}
	summary_table$ci_mean95upper[is.na(summary_table$ci_mean95upper)]=Inf
	summary_table$ci_mean95lower[is.na(summary_table$ci_mean95lower)]=-Inf

	summary_table
}

cumm_drift_plot_by_covariate<-function(data_local, trait, covariate, instrument_id, summary_table, study_name="UK Biobank")
{

	if(!hasArg(summary_table))
	{	
		if(!hasArg(instrument_id))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate)
		}
		if(hasArg(instrument_id))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)
		}
	}
	
	eval(parse(text=sprintf("summary_table$cumm_sum=cumsum(summary_table$mean[order(summary_table$%s)]-mean(summary_table$mean,na.rm=TRUE))[order(order(summary_table$%s))]", covariate,covariate)))
	eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=cumm_sum))", covariate)))
	plot=plot+geom_line()+theme_bw()
	plot=plot+get_covariate_xlab(covariate)+get_covariate_xscale(covariate)
	plot=plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))
	
	if(!hasArg(instrument_id))
	{ 	
		if(grepl("_adj", trait))
		{
			plot=plot+labs(title=sprintf("%s: Drift plot of %s with adjustment", study_name, get_trait_short_name(get_base_trait(trait))))
		}
		else
		{
			plot=plot+labs(title=sprintf("%s: Drift plot of %s", study_name, get_trait_short_name(get_base_trait(trait))))

		}
	}
	else
	{
		if(grepl("_adj", trait))
		{
			plot=plot+labs(title=sprintf("%s: Drift plot of %s with adjustment\n Instrument ID: %s",study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		else
		{
			plot=plot+labs(title=sprintf("%s: Drift plot of %s\n Instrument ID: %s", study_name, get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		

	}
	plot
	
}

hist_plot_by_binary_covariate<-function(data_local, trait, bin_strata, study_name="UK Biobank", annotate.df=NULL, quantile.limits=c(0.001,0.999))
{


  limits <- as.numeric(quantile(data_local[,trait], na.rm=TRUE, quantile.limits))
  if(!is.null(annotate.df)) {
    if (! any(colnames(annotate.df) == bin_strata) ||
        ! all(levels(annotate.df[[bin_strata]] %in%  levels(factor(data_local[,bin_strata])))) ||
        ! any(colnames(annotate.df) == trait)) {
      stop("annotate.df must have a column corresponding to bin_strata\n")
    }
    if(!any(colnames(annotate.df) == "colour")) {
      message("annoate.df does not have a colour column, defaulting to red")
      annotate.df[["colour"]] = rep("red",nrow(annotate.df))
    }
    limits <- c(min(annotate.df[,trait], limits[1],na.rm=TRUE), max(annotate.df[,trait], limits[2],na.rm=TRUE))
  }
	bin_width=(max(data_local[,trait], na.rm=TRUE)-min(data_local[,trait],na.rm=TRUE))/100
	bin_width=(quantile(data_local[,trait], 0.995, na.rm=TRUE)-quantile(data_local[,trait],1-0.995, na.rm=TRUE))/50

	eval(parse(text=sprintf("plot<-ggplot(data_local, aes(x=%s,fill=as.factor(%s))) + geom_histogram(aes(x=%s, y=-..density..), data=subset(data_local,%s==unique(data_local[,bin_strata])[1]), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=%s, y=..density..), data=subset(data_local,%s==unique(data_local[,bin_strata])[2]) ,alpha = 0.9, binwidth=bin_width)",trait,bin_strata, trait,bin_strata, trait, bin_strata)))

	plot <-plot+scale_x_continuous(limits = limits)
	plot<-plot+coord_flip() + theme_bw()
	plot<-plot+scale_fill_manual(name="", values=get_strata_graph_info(bin_strata)$colours,labels=get_strata_graph_info(bin_strata)$labels)

	plot<-plot+labs(title = sprintf("%s: Density of %s by %s", study_name, get_trait_short_name(get_base_trait(trait)), tolower(get_covariate_name(bin_strata))))
	plot<-plot+xlab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))
	plot<-plot+ylab(eval(parse(text=sprintf("expression(paste(\"Probability density (\",%s ,\")\"))",get_inverse_trait_linear_units(get_base_trait(trait))))))
  if(!is.null(annotate.df)) {
    ranges=list(-range(eval(parse(text=sprintf("density(subset(data_local,%s==unique(data_local[,bin_strata])[1])[[\"%s\"]], na.rm=TRUE)$y", bin_strata, trait)))),
            range(eval(parse(text=sprintf("density(subset(data_local,%s==unique(data_local[,bin_strata])[2])[[\"%s\"]], na.rm=TRUE)$y", bin_strata, trait)))))
    names(ranges)=unique(data_local[,bin_strata])[1:2]
    anno.list=list()
    for(stratum in unique(data_local[,bin_strata])) {
      for(i in which(annotate.df[,bin_strata]==stratum)) {
        anno.list=c(anno.list, annotate("segment", x=annotate.df[[trait]][i], y=ranges[[stratum]][1],
                      xend=annotate.df[[trait]][i], yend=mean(ranges[[stratum]]),
                      colour=sapply(annotate.df[["colour"]][i],as.character),
                      arrow=arrow(angle=30, ends="first",length=unit(8, "points"), type="closed")))
      }
    }
    plot <- plot + anno.list
  }
	plot
}

hist_plot_by_catagorical<-function(data_local, trait, strata, study_name="UK Biobank", quantile_lines=NULL, quantile.limits=c(0.001, 0.999))
{

	limits <- as.numeric(quantile(data_local[,trait], na.rm=TRUE, quantile.limits))
#	bin_width=(max(data_local[,trait], na.rm=TRUE)-min(data_local[,trait],na.rm=TRUE))/100
	bin_width=(quantile(data_local[,trait], 0.995, na.rm=TRUE)-quantile(data_local[,trait],1-0.995, na.rm=TRUE))/50


	eval(parse(text=sprintf("plot<-ggplot(data=subset(data_local,!is.na(data_local$%s)), aes(x=%s,fill=as.factor(%s))) + geom_histogram(aes(x=%s, y=..density..),  binwidth=bin_width)", strata, trait, strata, trait, strata)))
	
	plot<-plot+theme_bw()
	strata_labels=get_strata_graph_info(strata)$labels
	facet_labeller<-function(var){return(strata_labels[names(strata_labels)==var])}

	eval(parse(text=sprintf("plot<-plot+facet_wrap(~%s, drop=TRUE, labeller=as_labeller(strata_labels))",strata)))
	plot<-plot+scale_fill_manual(name="", values=get_strata_graph_info(strata)$colours,labels=get_strata_graph_info(strata)$labels, guide=FALSE)

	plot<-plot+labs(title = sprintf("%s: Density of %s by %s", study_name, get_trait_short_name(get_base_trait(trait)), tolower(get_covariate_name(strata))))
	plot<-plot+xlab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))
	plot<-plot+ylab(eval(parse(text=sprintf("expression(paste(\"Probability density (\",%s, \")\"))",get_inverse_trait_linear_units(get_base_trait(trait))))))
	if(!is.null(quantile_lines))
	{
		for(q in quantile_lines)
		{
			quant_q<-function(x){return(quantile(x,q,na.rm=TRUE))}
			eval(parse(text=sprintf("quant_data=subset(data_local[,c(trait,strata)],!is.na(data_local$%s)) %%>%% group_by(%s) %%>%% summarise(quant=quant_q(%s))",strata,strata,trait)))
			plot<-plot+geom_vline(data=quant_data, aes(xintercept=quant), size=0.5, colour="deeppink")	
		}	
	}
	plot<-plot+scale_x_continuous(limits = limits)

	plot
}




violin_plot_by_discrete_covariate<-function(data_local, trait, strata, study_name="UK Biobank")
{
	levels<-unique(data_local[,strata])
	eval(parse(text=sprintf("plot<-ggplot(subset(data_local, !is.na(data_local[,\"%s\"])), aes(x=%s,y=%s, fill=as.factor(%s))) + geom_violin(alpha = 0.9,showLegend=FALSE)", strata, strata, trait, strata)))
	
	plot<-plot+coord_flip() + theme_bw()
	plot<-plot+scale_fill_manual(name="", values=get_strata_graph_info(strata, levels)$colours,  guide=FALSE)
	plot<-plot+scale_x_discrete(labels=get_strata_graph_info(strata,levels)$labels)
	plot<-plot+labs(title = sprintf("%s: Violin plot of %s by %s", study_name, get_trait_short_name(get_base_trait(trait)), tolower(get_covariate_name(strata))))
	plot<-plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))
	plot<-plot+xlab(get_covariate_name(strata))

	plot
}

generate_bivariate_plot_by_strata<-function(data_local, trait1, trait2, strata, study_name="UK Biobank")
{

	eval(parse(text=sprintf("plot<-ggplot(data_local,aes(x=%s,y=%s))+facet_grid(~%s, labeller=covariate_facet_labeller)", trait1, trait2, strata)))
	plot<-plot+theme_bw()
	plot<-plot+geom_point(size=0.8, colour="tan4")
#	plot<-plot+scale_colour_manual(values=get_strata_graph_info(strata,levels)$colours)
#	rgb_colours=t(col2rgb(get_strata_graph_info("sex")$colours))
#	rgb_colours_shifted=pmin(rgb_colours+100,matrix(255, dim(rgb_colours)[1], dim(rgb_colours)[2]))
#	rgb_colours_hex=rgb(rgb_colours[,1]/255, rgb_colours[,2]/255, rgb_colours[,3]/255, names=rownames(rgb_colours))
#	rgb_colours_shifted_hex=rgb(rgb_colours_shifted[,1]/255, rgb_colours_shifted[,2]/255, rgb_colours_shifted[,3]/255, names=rownames(rgb_colours_shifted))


	eval(parse(text=sprintf("plot<-plot+stat_density2d(aes(fill=..level..), geom=\"polygon\")", strata)))
	plot<-plot+ scale_fill_gradient2(low="tan",high="tan4", guide=FALSE)
	plot<-plot+xlab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait1)),get_trait_linear_units(get_base_trait(trait1))))))
	plot<-plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait2)),get_trait_linear_units(get_base_trait(trait2))))))
	plot<-plot+labs(title = sprintf("%s: Joint density of %s and %s by %s", study_name, get_trait_short_name(get_base_trait(trait1)), get_trait_short_name(get_base_trait(trait2)), tolower(get_covariate_name(strata))))
	plot
}


summary_plot_by_covariate<-function(data_local,trait, covariate, strata, instrument_id, summary_table, study_name="UK Biobank")
{	

	if(!hasArg(summary_table))
	{	
		if((!hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate)
		}
		if((!hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata)
		}
		if((hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)
		}
		if((hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata, instrument_id=instrument_id)
		}
	}
	if(!hasArg(strata))
	{
		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=mean))", covariate)))
	}
	else
	{
		levels<-unique(data_local[,strata])
		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=mean, group=%s, colour=%s))", covariate, strata, strata)))
		plot<-plot+scale_color_manual(name="", values=get_strata_graph_info(strata,levels)$colours,labels=get_strata_graph_info(strata,levels)$labels)
	}
	
	plot=plot+theme_bw()+geom_ribbon(data=summary_table,aes(ymin=quantile_0.00,ymax=quantile_1.00),alpha=0.1)
	plot=plot+geom_line(data=summary_table, aes(y=quantile_0.25), linetype="dashed")
	plot=plot+geom_line(data=summary_table, aes(y=quantile_0.75), linetype="dashed")
	plot=plot+geom_ribbon(data=summary_table,aes(ymin=ci_mean95lower,ymax=ci_mean95upper),alpha=0.3, colour=NA) + geom_point()+geom_line()
	plot=plot+get_covariate_xlab(covariate)+get_covariate_xscale(covariate)
	plot=plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))
	
	if(!hasArg(instrument_id))
	{ 	
		if(grepl("_adj", trait))
		{
			plot=plot+labs(title=sprintf("%s: Summary of %s with adjustment", study_name, get_trait_short_name(get_base_trait(trait))))
		}
		else
		{
			plot=plot+labs(title=sprintf("%s: Summary of %s", study_name, get_trait_short_name(get_base_trait(trait))))

		}
	}
	else
	{
		if(grepl("_adj", trait))
		{
			plot=plot+labs(title=sprintf("%s: Summary of %s with adjustment\n Instrument ID: %s",study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		else
		{
			plot=plot+labs(title=sprintf("%s: Summary of %s\n Instrument ID: %s", study_name, get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		

	}
	plot

}

iq_plot_by_covariate<-function(data_local,trait, covariate, strata, instrument_id, summary_table, study_name="UK Biobank")
{

	if(!hasArg(summary_table))
	{	
		if((!hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate)
		}
		if((!hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata)
		}
		if((hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)
		}
		if((hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata, instrument_id=instrument_id)
		}
	}
	if(!hasArg(strata))
	{
		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=mean))", covariate)))
	}
	else
	{
		levels<-unique(data_local[,strata])

		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=mean, group=%s, colour=%s))", covariate, strata, strata)))
		plot<-plot+scale_color_manual(name="", values=get_strata_graph_info(strata,levels)$colours,labels=get_strata_graph_info(strata,levels)$labels)
	}
	

	plot=plot+theme_bw()
	plot=plot+geom_line(data=summary_table, aes(y=quantile_0.25),colour="grey", linetype="dashed")
	plot=plot+geom_line(data=summary_table, aes(y=quantile_0.75),colour="grey", linetype="dashed")
	plot=plot+geom_ribbon(data=summary_table,aes(ymin=ci_mean95lower,ymax=ci_mean95upper),alpha=0.3, colour=NA) + geom_point()+geom_line()
	plot=plot+get_covariate_xlab(covariate)+get_covariate_xscale(covariate)
	plot=plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))

	if(!hasArg(instrument_id))
	{ 	
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: mean and IQ range of %s with adjustment", study_name, get_trait_short_name(get_base_trait(trait))))
		}
		else{	
			plot=plot+labs(title=sprintf("%s: mean and IQ range of %s", study_name,  get_trait_short_name(get_base_trait(trait))))
		}

	}
	else
	{
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: mean and IQ range of %s with adjustment\n Instrument ID: %s", study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		else
		{
			plot=plot+labs(title=sprintf("%s: mean and IQ range of %s\n Instrument ID: %s", study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id))      
		}
	}



       	plot

}	

mean_plot_by_covariate<-function(data_local,trait, covariate, strata, instrument_id, summary_table, study_name="UK Biobank")
{
	if(!hasArg(summary_table))
	{	
		if((!hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate)
		}
		if((!hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata)
		}
		if((hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)
		}
		if((hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata, instrument_id=instrument_id)
		}
	}
	if(!hasArg(strata))
	{
		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=mean))", covariate)))
	}
	else
	{
		levels<-unique(data_local[,strata])

		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=mean, group=%s, colour=%s))", covariate, strata, strata)))
		plot<-plot+scale_color_manual(name="", values=get_strata_graph_info(strata,levels)$colours,labels=get_strata_graph_info(strata,levels)$labels)
	}
	

	plot=plot+theme_bw()
	plot=plot+geom_ribbon(data=summary_table,aes(ymin=ci_mean95lower,ymax=ci_mean95upper),alpha=0.3, colour=NA) + geom_point()+geom_point()
	plot=plot+get_covariate_xlab(covariate)+get_covariate_xscale(covariate)
	plot=plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))
	if(!hasArg(instrument_id))
	{ 	
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: Mean of %s with adjustment", study_name, get_trait_short_name(get_base_trait(trait))))
		}
		else{	
			plot=plot+labs(title=sprintf("%s: Mean of %s", study_name, get_trait_short_name(get_base_trait(trait))))
		}

	}
	else
	{
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: Mean of %s with adjustment\n Instrument ID: %s", study_name, get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		else
		{
			plot=plot+labs(title=sprintf("%s: Mean of %s\n Instrument ID: %s", study_name, get_trait_short_name(get_base_trait(trait)), instrument_id))      
		}
	}


       	plot

}


sd_plot_by_covariate<-function(data_local,trait, covariate, strata, instrument_id,summary_table, study_name="UK Biobank")
{
	if(!hasArg(summary_table))
	{	
		if((!hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate)
		}
		if((!hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata)
		}
		if((hasArg(instrument_id))&&(!hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)
		}
		if((hasArg(instrument_id))&&(hasArg(strata)))
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, strata=strata, instrument_id=instrument_id)
		}
	}
	if(!hasArg(strata))
	{
		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=st_dev))", covariate)))
	}
	else
	{
		levels<-unique(data_local[,strata])
		eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=as.numeric(%s), y=st_dev, group=%s, colour=%s))", covariate, strata, strata)))
		plot<-plot+scale_color_manual(name="", values=get_strata_graph_info(strata,levels)$colours,labels=get_strata_graph_info(strata,levels)$labels)
	}
	

	plot=plot+theme_bw()
	plot=plot+geom_ribbon(data=summary_table,aes(ymin=sqrt(ci_var95lower),ymax=sqrt(ci_var95upper)),alpha=0.3, colour=NA) + geom_point()
	plot=plot+get_covariate_xlab(covariate)+get_covariate_xscale(covariate)
	plot=plot+ylab(eval(parse(text=sprintf("expression(paste(\"%s (\",%s, \")\"))",get_trait_long_name(get_base_trait(trait)),get_trait_linear_units(get_base_trait(trait))))))

	if(!hasArg(instrument_id))
	{ 	
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: Standard Deviation of %s with adjustment", study_name, get_trait_short_name(get_base_trait(trait))))
		}
		else{	
			plot=plot+labs(title=sprintf("%s: Standard Deviation of %s", study_name, get_trait_short_name(get_base_trait(trait))))
		}

	}
	else
	{
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: Standard Deviation of %s with adjustment\n Instrument ID: %s", study_name, get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		else
		{
			plot=plot+labs(title=sprintf("%s: Standard Deviation of %s\n Instrument ID: %s", study_name, get_trait_short_name(get_base_trait(trait)), instrument_id))      
		}
	}

       	plot

}	

log_var_vs_log_mean_plot_by_covariate<-function(data_local,trait, covariate, instrument_id, summary_table, study_name="UK Biobank")
{
	if(!hasArg(summary_table))
	{	
		if(!hasArg(instrument_id)) 
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate)
		}
		else
		{
			summary_table=summarise_trait_by_covariate(data_local, trait, covariate, instrument_id=instrument_id)

		}	
	}
	summary_table$log_mean=log(summary_table$mean)
	summary_table$log_var=log(summary_table$var)
	eval(parse(text = sprintf("plot<-ggplot(data=summary_table, aes(x=log_mean, y=log_var))", covariate)))
	plot=plot+theme_bw()+ geom_point()
 	plot=plot + geom_smooth(method='lm',formula=y~x)+stat_summary(fun.data=mean_cl_normal) 
	plot=plot+ylab(eval(parse(text=sprintf("expression(paste(\"log-variance %s\"))",get_trait_long_name(get_base_trait(trait))))))
	plot=plot+xlab(eval(parse(text=sprintf("expression(paste(\"log-mean %s\"))",get_trait_long_name(get_base_trait(trait))))))
	subset_summ=(abs(summary_table$mean)>0)&(summary_table$var>0)

	plot=plot+annotate("text",x=0.9*min(summary_table$log_mean[subset_summ], na.rm=TRUE)+0.1*max(summary_table$log_mean[subset_summ], na.rm=TRUE), y=max(summary_table$log_var[subset_summ],na.rm=TRUE),label=sprintf("Gradient=%0.2f",lm(log_var~log_mean, data=summary_table[subset_summ,])$coeff[2]))

	if(!hasArg(instrument_id))
	{ 	
		plot=plot+labs(title=sprintf("%s: log mean vs log variance of %s\n by %s", study_name, get_trait_short_name(get_base_trait(trait)), covariate))
	}
	else
	{
		plot=plot+labs(title=sprintf("%s: log mean vs log variance of %s\n Instrument ID: %s by %s", study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id, covariate))      
 
	}
	if(!hasArg(instrument_id))
	{ 	
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: log mean vs log variance of %s with adjustment", study_name,  get_trait_short_name(get_base_trait(trait))))
		}
		else{	
			plot=plot+labs(title=sprintf("%s: log mean vs log variance of %s", study_name,  get_trait_short_name(get_base_trait(trait))))
		}

	}
	else
	{
		if(grepl("_adj", trait))
		{	
			plot=plot+labs(title=sprintf("%s: log mean vs log variance of %s with adjustment\n Instrument ID: %s", study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id))      
 		}
		else
		{
			plot=plot+labs(title=sprintf("%s: log mean vs log variance of %s\n Instrument ID: %s", study_name,  get_trait_short_name(get_base_trait(trait)), instrument_id))      
		}
	}

       	plot
}	

	


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
	# Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



generate_correlation_plot_wbc_plt<-function(data_local)
{

	print("This needs updating")
	# remove outliers, locally

	quan_wbc=ecdf(data_local$wbc)(data_local$wbc)
	quan_neut=ecdf(data_local$neut)(data_local$neut)
	quan_lymph=ecdf(data_local$lymph)(data_local$lymph)
	quan_mono=ecdf(data_local$mono)(data_local$mono)
#	quan_baso=ecdf(data_local$baso)(data_local$baso)
	quan_eo=ecdf(data_local$eo)(data_local$eo)
	quan_plt=ecdf(data_local$plt)(data_local$plt)

	local_outlier_filter<-	quan_wbc<0.99 & quan_wbc > 0.01 &
			       	quan_neut<0.99 & quan_neut > 0.01 &
			       	quan_mono<0.99 & quan_mono > 0.01 &
			       	quan_lymph<0.99 & quan_lymph > 0.01 &
 			#	quan_baso<0.99 & quan_baso > 0.01 &
			       	quan_plt<0.99 & quan_plt > 0.01 &
 				quan_eo<0.99 & quan_eo > 0.01 
			    

	data_local=data_local[local_outlier_filter,]
	data_local=data_local[apply(!is.na(data_local),1,all),]

	data_local$rounded_bmi=round(data_local$bmi)
	data_local$wbc_age_sex_bmi_adj=residuals(lm(wbc~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	data_local$neut_age_sex_bmi_adj=residuals(lm(neut~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	data_local$mono_age_sex_bmi_adj=residuals(lm(mono~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	data_local$lymph_age_sex_bmi_adj=residuals(lm(lymph~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
#	data_local$baso_age_sex_bmi_adj=residuals(lm(baso~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	data_local$eo_age_sex_bmi_adj=residuals(lm(eo~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	data_local$plt_age_sex_bmi_adj=residuals(lm(plt~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	
	plot_data=data_local[,c("wbc_age_sex_bmi_adj","neut_age_sex_bmi_adj","mono_age_sex_bmi_adj","eo_age_sex_bmi_adj","lymph_age_sex_bmi_adj","plt_age_sex_bmi_adj")]

	colnames(plot_data)=c("WBC", "NEUT","MONO","EO","LYMPH", "PLT")
	pdf(sprintf("%s/output/graphics/bivariate_counts_wbc_plt.pdf",UKBIOBANK_FBC_REPO_DIR))	
	plot_pairwise_kde(plot_data, main="Bivariate Blood Count Distributions Estimated From UKBB data\nAdjusted for Age, Sex and BMI")
	dev.off()

	data_local$wbc_age_sex_adj=residuals(lm(wbc~sex+as.factor(age_acq), data=data_local,na.action=na.exclude))
	data_local$neut_age_sex_adj=residuals(lm(neut~sex+as.factor(age_acq), data=data_local,na.action=na.exclude))
	data_local$mono_age_sex_adj=residuals(lm(mono~sex+as.factor(age_acq), data=data_local,na.action=na.exclude))
	data_local$lymph_age_sex_adj=residuals(lm(lymph~sex+as.factor(age_acq), data=data_local,na.action=na.exclude))
#	data_local$baso_age_sex_bmi_adj=residuals(lm(baso~sex+as.factor(age_acq)+as.factor(rounded_bmi), data=data_local,na.action=na.exclude))
	data_local$eo_age_sex_adj=residuals(lm(eo~sex+as.factor(age_acq), data=data_local,na.action=na.exclude))
	data_local$plt_age_sex_adj=residuals(lm(plt~sex+as.factor(age_acq), data=data_local,na.action=na.exclude))
	plot_data=data_local[data_local$bmi<30,c("wbc_age_sex_adj","neut_age_sex_adj","mono_age_sex_adj","eo_age_sex_adj","lymph_age_sex_adj","plt_age_sex_adj")]
	colnames(plot_data)=c("WBC", "NEUT","MONO","EO","LYMPH", "PLT")

	pdf(sprintf("%s/output/graphics/bivariate_counts_wbc_plt_low_bmi.pdf",UKBIOBANK_FBC_REPO_DIR))	
	plot_pairwise_kde(plot_data, main="Bivariate Blood Count Distributions Estimated From UKBB data\nAdjusted for Age, Sex; BMI<30")
	dev.off()

	plot_data=data_local[data_local$bmi>30,c("wbc_age_sex_adj","neut_age_sex_adj","mono_age_sex_adj","eo_age_sex_adj","lymph_age_sex_adj","plt_age_sex_adj")]
	colnames(plot_data)=c("WBC", "NEUT","MONO","EO","LYMPH", "PLT")

	pdf(sprintf("%s/output/graphics/bivariate_counts_wbc_plt_high_bmi.pdf",UKBIOBANK_FBC_REPO_DIR))	
	plot_pairwise_kde(plot_data, main="Bivariate Blood Count Distributions Estimated From UKBB data\nAdjusted for Age, Sex; BMI>30")
	dev.off()


}

generate_lymph_neut_plot<-function(data_local)
{
	data_local$neut_lymph_ratio<-data_local$neut/data_local$lymph
	
	data_local$sex<-factor(data_local$sex,levels=c("0","1"), labels=c("Female", "Male"))
	log_neut_lymph_ratio_plot<-ggplot(data=data_local) + geom_histogram(aes(x=neut_lymph_ratio, y=..density..),data=data_local, binwidth=0.02)+facet_grid(sex~.)+theme_bw()+ylab("Density")+xlab("Neutrophil Count / Lymphocyte Count")+scale_x_log10()+labs(title="Ratio of Neutrophil to Lymphocyte Count in UK Biobank")+geom_segment(x=0, y=0, xend=0, yend=3, linetype=2)
	print(range(log10(data_local$neut)-log10(data_local$lymph),na.rm=TRUE, finite=TRUE))
	

	pdf(sprintf("%s/output/graphics/NEUT_LYMPH_ratio.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(log_neut_lymph_ratio_plot)	
	dev.off()

}

generate_menopause_age_plots<-function(data_local, output_dir=UKBIOBANK_OUTPUT_DIR)
{

	# age of menopause
	
	bin_width=1

	age_meno_plot<-ggplot(data_local, aes(x=age_acq,fill=as.factor(post_meno))) + geom_histogram(aes(x=age_acq, y=..density..), data=subset(data_local,post_meno==TRUE), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=age_acq, y=-..density..), data=subset(data_local,post_meno==FALSE) ,alpha = 0.9, binwidth=bin_width)+coord_flip() + theme_bw()+labs(title = sprintf("Density of Age in UK Biobank Females Stratified by Menopause"))+xlab("Age (years)")+ylab("Density")+scale_fill_manual(name="", values=c("indianred4","seagreen4"), labels=c("Pre-menopausal", "Post-menopausal"))+scale_y_continuous(limits=c(-0.1,0.1))

	pdf(sprintf("%s/graphics/age_meno.pdf", output_dir))
	print(age_meno_plot)	
	dev.off()

	
}

comment<-function()
{
draw_UKBB_plots<-function()
{

	for(trait in ukbb_traits)
	{
		print(trait)
		print("we should be passing the filter here")
		print("temp bd_extract->bd here")
	#	eval(parse(text=sprintf("tmp_data_filter=(!is.na(bd$sex))&(!is.na(bd$instrument))&(!is.na(bd$%s))&(abs(bd$%s)<Inf)", trait, trait)))
	#	generate_analysis_count_plots(trait, bd[tmp_data_filter,])
		
	#	generate_qc_plots(trait, bd)

		#generate_analysis_count_plots(trait, bd_extract[data_filter,])
	#		generate_qc_plots(trait, bd[tmp_data_filter,], int_qc_vars, output_dir=UKBIOBANK_FBC_OUTPUT_DIR, "UK Biobank")
		    generate_covariate_plots(trait, bd[tmp_data_filter,], int_qc_vars, output_dir=UKBIOBANK_FBC_OUTPUT_DIR, "UK Biobank")
	#	generate_violin_plots_by_hi_low_med_bmi(trait, bd_extract, data_filter)
	}


	summary_table<-ddply(bd_extract[data_filter,c("yob", "age_acq","sex")], .(yob,sex), summarise, mean=mean(age_acq,na.rm=TRUE), upper=quantile(age_acq,0.95,na.rm=TRUE)[[1]], lower=quantile(age_acq,0.05,na.rm=TRUE)[[1]], ci_mean95upper=mean(age_acq,na.rm=TRUE)+1.96*sd(age_acq,na.rm=TRUE)/sqrt(sum(!is.na(age_acq))),  ci_mean95lower=mean(age_acq,na.rm=TRUE)-1.96*sd(age_acq,na.rm=TRUE)/sqrt(sum(!is.na(age_acq))),n=sum(!is.na(age_acq)))
	
	mean_age_acq_sex_yob_summary_plot<-ggplot(data=summary_table, aes(x=yob, y=mean, group=sex, colour=as.factor(sex)))+theme_bw()+scale_color_manual(name="", values=c("orange", "blue"), labels=c("Females", "Males"))+xlab("Year of Birth")+ylab("Age at Acquisition")+labs(title=sprintf("UK Biobank, Age at Aquisition by Year of Birth"))+geom_ribbon(data=summary_table,aes(ymin=lower,ymax=upper),alpha=0.05)+geom_ribbon(data=summary_table,aes(ymin=ci_mean95lower,ymax=ci_mean95upper),alpha=0.3, colour=NA)+ geom_line() + geom_point()

   	pdf(sprintf("%s/output/graphics/age_aq_sex_yob_filt_mean_summary.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(mean_age_acq_sex_yob_summary_plot)
	dev.off()
}
}


draw_plots<-function(traits=tf_adjustable(tf_coulter()), postfix="", study_name="UK Biobank")
{
	if(study_name=="UK Biobank")
	{
		output_dir="UKBIOBANK_FBC_OUTPUT_DIR"
		fbc_dataset="bd_extract"
		filters="ukbb_data_filters"
	}
	if(study_name=="INTERVAL")
	{
		output_dir="INTERVAL_FBC_OUTPUT_DIR"
		fbc_dataset="int_extract"
		filters="int_data_filters"
	}

	generate_plots_wrap<-function(traits, trait_num)
	{	
		eval(parse(text=sprintf("generate_qc_plots(%s, %s, \"%s%s\", output_dir=%s, study_name=\"%s\")", fbc_dataset, filters, traits[trait_num],postfix, output_dir, study_name)))
		eval(parse(text=sprintf("generate_covariate_plots(%s, %s, \"%s%s\", output_dir=%s, study_name=\"%s\")", fbc_dataset, filters, traits[trait_num],postfix, output_dir, study_name)))
		eval(parse(text=sprintf("generate_histograms(%s, %s, \"%s%s\", output_dir=%s, study_name=\"%s\")", fbc_dataset, filters, traits[trait_num],postfix, output_dir, study_name)))

	}
	
	foreach(n = 1:length(traits)) %dopar% generate_plots_wrap(traits, n)						
	
}



draw_UKBB_tail_plots<-function()
{
	generate_tech_plots_wrap<-function(bd_extract, ukbb_data_filters, traits, trait_num)
	{	
		trait=sprintf("%s_tech_adj",traits[trait_num])
		bin_width=(quantile(bd_extract[,trait], 0.995, na.rm=TRUE)-quantile(bd_extract[,trait],1-0.995, na.rm=TRUE))/50


	
		ids=read.table(file=sprintf("%s/output/tail_sequencing/initial_4000.tsv", UKBIOBANK_DATA_DIR),head=TRUE, sep=" ",stringsAsFactors=FALSE)
		lower_tail=ids$MCV/ids$RBC>20
		dir.create(sprintf("%s/graphics",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)
	
		lower_logical=is.element(bd_extract$subject_id, ids$f.eid_13745[lower_tail])
		upper_logical=is.element(bd_extract$subject_id, ids$f.eid_13745[!lower_tail])
		
		upper_data=bd_extract[upper_logical,]
	#	upper_data$meno_colour=upper_data$meno_simple
		upper_data=upper_data[!is.na(upper_data$meno_simple),]
		
		lower_data=bd_extract[lower_logical,]
	#	lower_data$meno_colour=lower_data$meno_simple
		lower_data=lower_data[!is.na(lower_data$meno_simple),]
	
		minlim=min(c(max(c(min(upper_data[,trait],na.rm=TRUE),median(bd_extract[,trait],na.rm=TRUE)-15*mad(bd_extract[,trait],na.rm=TRUE))),quantile(bd_extract[,trait],1-0.995,na.rm=TRUE)))
		maxlim=max(c(min(c(max(upper_data[,trait],na.rm=TRUE),median(bd_extract[,trait],na.rm=TRUE)+15*mad(bd_extract[,trait],na.rm=TRUE))),quantile(bd_extract[,trait],0.995,na.rm=TRUE)))
		histogram_upper_by_meno_simple<-hist_plot_by_catagorical(bd_extract[is.element(bd_extract$meno_simple,c("male", "pre","post")),], trait, "meno_simple", "UK Biobank TRED High RBC# / Low MCV Tail", quantile_lines=c(0.025, 0.975))
		
		minlim=min(c(max(c(min(lower_data[,trait],na.rm=TRUE),median(bd_extract[,trait],na.rm=TRUE)-15*mad(bd_extract[,trait],na.rm=TRUE))),quantile(bd_extract[,trait],1-0.995,na.rm=TRUE)))
		maxlim=max(c(min(c(max(lower_data[,trait],na.rm=TRUE),median(bd_extract[,trait],na.rm=TRUE)+15*mad(bd_extract[,trait],na.rm=TRUE))),quantile(bd_extract[,trait],0.995,na.rm=TRUE)))


		histogram_lower_by_meno_simple<-hist_plot_by_catagorical(bd_extract[is.element(bd_extract$meno_simple,c("male", "pre","post")),], trait, "meno_simple", "UK Biobank TRED Low RBC# / High MCV Tail", quantile_lines=c(0.025, 0.975))



		eval(parse(text=sprintf("histogram_upper_by_meno_simple<-histogram_upper_by_meno_simple+geom_histogram(data=upper_data,aes(x=%s, y=..density..),fill=alpha(\"black\", 0.5), binwidth=bin_width)", trait)))

		eval(parse(text=sprintf("histogram_lower_by_meno_simple<-histogram_lower_by_meno_simple+geom_histogram(data=lower_data,aes(x=%s, y=..density..),fill=alpha(\"black\", 0.5), binwidth=bin_width)", trait)))


		dir.create(sprintf("%s/graphics/TRED_histogram_upper_by_meno_simple",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)
		dir.create(sprintf("%s/graphics/TRED_histogram_lower_by_meno_simple",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)
	

		pdf(sprintf("%s/graphics/TRED_histogram_upper_by_meno_simple/%s_TRED_histogram_upper_by_meno_simple.pdf",UKBIOBANK_FBC_OUTPUT_DIR,toupper(trait)), height=7, width=14)
		print(histogram_upper_by_meno_simple)
		dev.off()
		pdf(sprintf("%s/graphics/TRED_histogram_lower_by_meno_simple/%s_TRED_histogram_lower_by_meno_simple.pdf",UKBIOBANK_FBC_OUTPUT_DIR,toupper(trait)), height=7, width=14)
		print(histogram_lower_by_meno_simple)
		dev.off()

	

	}
	ukbb_plot_traits=c(ukbb_measured_traits, ukbb_derived_traits)
	registerDoMC(cores=length(ukbb_plot_traits))
	foreach(n = 1:length(ukbb_plot_traits)) %dopar% generate_tech_plots_wrap(bd_extract, ukbb_data_filters, ukbb_plot_traits, n)	


}



draw_INTERVAL_tech_plots<-function()
{
	generate_tech_plots_wrap<-function(int_extract, int_data_filters, traits, trait_num)
	{	
		eval(parse(text=sprintf("generate_qc_plots(int_extract, int_data_filters, \"%s_tech_adj\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")", traits[trait_num])))
		eval(parse(text=sprintf("generate_covariate_plots(int_extract, int_data_filters, \"%s_tech_adj\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")",  traits[trait_num])))
		eval(parse(text=sprintf("generate_histograms(int_extract, int_data_filters, \"%s_tech_adj\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")",  traits[trait_num])))

	}

	int_plot_traits=c(int_measured_traits, int_derived_traits)
	registerDoMC(cores=length(int_plot_traits))
	foreach(n = 1:length(int_plot_traits)) %dopar% generate_tech_plots_wrap(int_extract, int_data_filters, int_plot_traits, n)	

}



draw_UKBB_pre_plots<-function()
{
	generate_tech_plots_wrap<-function(bd_extract, ukbb_data_filters, traits, trait_num)
	{	
		eval(parse(text=sprintf("generate_qc_plots(bd_extract, ukbb_data_filters, \"%s\", output_dir=UKBIOBANK_FBC_OUTPUT_DIR, study_name=\"UK Biobank\")", traits[trait_num])))

	}

	ukbb_plot_traits=c(ukbb_measured_traits, ukbb_derived_traits)
	registerDoMC(cores=length(ukbb_plot_traits))
	foreach(n = 1:length(ukbb_plot_traits)) %dopar% generate_tech_plots_wrap(bd_extract, ukbb_data_filters, ukbb_plot_traits, n)	

}


draw_INTERVAL_pre_plots<-function()
{
	generate_tech_plots_wrap<-function(int_extract, int_data_filters, traits, trait_num)
	{	
		eval(parse(text=sprintf("generate_qc_plots(int_extract, int_data_filters, \"%s\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")", traits[trait_num])))

	}

	int_plot_traits=c(int_measured_traits, int_derived_traits)
	registerDoMC(cores=length(int_plot_traits))
	foreach(n = 1:length(int_plot_traits)) %dopar% generate_tech_plots_wrap(int_extract, int_data_filters, int_plot_traits, n)	

}



draw_UKBB_tech_plots<-function()
{
	generate_tech_plots_wrap<-function(bd_extract, ukbb_data_filters, traits, trait_num)
	{	
		eval(parse(text=sprintf("generate_qc_plots(bd_extract, ukbb_data_filters, \"%s_tech_adj\", output_dir=UKBIOBANK_FBC_OUTPUT_DIR, study_name=\"UK Biobank\")", traits[trait_num])))
		eval(parse(text=sprintf("generate_covariate_plots(bd_extract, ukbb_data_filters, \"%s_tech_adj\", output_dir=UKBIOBANK_FBC_OUTPUT_DIR, study_name=\"UK Biobank\")",  traits[trait_num])))
		eval(parse(text=sprintf("generate_histograms(bd_extract, ukbb_data_filters, \"%s_tech_adj\", output_dir=UKBIOBANK_FBC_OUTPUT_DIR, study_name=\"UK Biobank\")",  traits[trait_num])))

	}

	ukbb_plot_traits=c(ukbb_measured_traits, ukbb_derived_traits)
	registerDoMC(cores=length(ukbb_plot_traits))
	foreach(n = 1:length(ukbb_plot_traits)) %dopar% generate_tech_plots_wrap(bd_extract, ukbb_data_filters, ukbb_plot_traits, n)	

}





draw_INTERVAL_gwas_plots<-function()
{
	generate_gwas_plots_wrap<-function(int_extract, int_data_filters, traits, trait_num)
	{
		#eval(parse(text=sprintf("generate_qc_plots(int_extract, int_data_filters, \"%s_gwas_adj\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")",  traits[trait_num])))
		#eval(parse(text=sprintf("generate_covariate_plots(int_extract, int_data_filters, \"%s_gwas_adj\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")",  traits[trait_num])))
		eval(parse(text=sprintf("generate_histograms(int_extract, int_data_filters, \"%s_gwas_adj\", output_dir=INTERVAL_FBC_OUTPUT_DIR, study_name=\"INTERVAL\")",  traits[trait_num])))


	}
	registerDoMC(cores=length(int_gwas_150k_traits))
	
	#foreach(n = 1:length(int_gwas_150k_traits)) %dopar% generate_gwas_plots_wrap(int_extract, int_data_filters, int_gwas_150k_traits, n)	

	foreach(n = 1:1) %dopar% generate_gwas_plots_wrap(int_extract, int_data_filters, int_gwas_150k_traits, n)	
}


draw_data_clean_plots<-function()
{
	day_min_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition"))+geom_point(size=1.2)

	eo_high_res_err_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: EO Hi-Res Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !eo_high_res_filter), colour=2, size=1.2)

	
	baso_high_res_err_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: BASO Hi-Res Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !baso_high_res_filter), colour=2, size=1.2)


	mono_high_res_err_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: MONO Hi-Res Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !mono_high_res_filter), colour=2, size=1.2)


	neut_high_res_err_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: NEUT Hi-Res Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !neut_high_res_filter), colour=2, size=1.2)


	copied_err_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: Copied Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !copied_data_filter), colour=2, size=1.2)

	zero_iretp_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: Zero iretp Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !zero_iretp_filter), colour=2, size=1.2)

	zero_ret_plot<-ggplot(data=bd_extract, aes(x=day_of_study_acq, y=min_acq))+theme_bw()+xlab("Day of Acquisition")+ylab("Minute of Acquisition")+labs(title=sprintf("UK Biobank, Day vs Min of Aquisition: Zero ret Data"))+geom_point(size=1.2)+geom_point(data=subset(bd_extract, !zero_ret_filter), colour=2, size=1.2)




  	   # machine times for filtered outcomes:	

	png(sprintf("%s/output/graphics/day_min_acq.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(day_min_plot)
	dev.off()

	png(sprintf("%s/output/graphics/day_min_acq_eo_high_res.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(eo_high_res_err_plot)
	dev.off()

	png(sprintf("%s/output/graphics/day_min_acq_baso_high_res.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(baso_high_res_err_plot)
	dev.off()

	png(sprintf("%s/output/graphics/day_min_acq_mono_high_res.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(mono_high_res_err_plot)
	dev.off()
	
	png(sprintf("%s/output/graphics/day_min_acq_neut_high_res.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(neut_high_res_err_plot)
	dev.off()


	png(sprintf("%s/output/graphics/day_min_acq_copied.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(copied_err_plot)
	dev.off()

	
	png(sprintf("%s/output/graphics/day_min_acq_zero_iretp.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(zero_iretp_plot)
	dev.off()

	png(sprintf("%s/output/graphics/day_min_acq_zero_ret.png",UKBIOBANK_FBC_REPO_DIR), width=1280, height=1280)	
	print(zero_ret_plot)
	dev.off()

	## focused plots showing discontiuity in histogram:
	
	#hgb	
	sort_diff=diff(sort(bd_extract$hgb))
	bin_width=min(sort_diff[sort_diff>0])/2
	hgb_main_sex_plot<-ggplot(subset(bd_extract, data_filter), aes(x=hgb,fill=as.factor(sex))) + geom_histogram(aes(x=hgb, y=..density..), data=subset(bd_extract,data_filter&(sex==1)), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=hgb, y=-..density..), data=subset(bd_extract,data_filter&(sex==0)) ,alpha = 0.9, binwidth=bin_width)
	hgb_sex_plot<-hgb_main_sex_plot+coord_flip() + theme_bw()+labs(title = sprintf("Density of %s in UK Biobank Stratified by Sex",toupper("hgb")))+ylab("Density")+xlab("HGB")+scale_fill_manual(name="", values=c("orange","blue"), labels=c("Females", "Males"))+scale_x_continuous(limits = c(10, 18))

	#hct	
	sort_diff=diff(sort(bd_extract$hct))
	bin_width=min(sort_diff[sort_diff>0])/2
	hct_main_sex_plot<-ggplot(subset(bd_extract, data_filter), aes(x=hct,fill=as.factor(sex))) + geom_histogram(aes(x=hct, y=..density..), data=subset(bd_extract,data_filter&(sex==1)), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=hct, y=-..density..), data=subset(bd_extract,data_filter&(sex==0)) ,alpha = 0.9, binwidth=bin_width)
	hct_sex_plot<-hct_main_sex_plot+coord_flip() + theme_bw()+labs(title = sprintf("Density of %s in UK Biobank Stratified by Sex",toupper("hct")))+ylab("Density")+xlab("HCT")+scale_fill_manual(name="", values=c("orange","blue"), labels=c("Females", "Males"))+scale_x_continuous(limits = c(0.35,0.45))

	#rbc	
	sort_diff=diff(sort(bd_extract$rbc))
	bin_width=min(sort_diff[sort_diff>0])/2
	rbc_main_sex_plot<-ggplot(subset(bd_extract, data_filter), aes(x=rbc,fill=as.factor(sex))) + geom_histogram(aes(x=rbc, y=..density..), data=subset(bd_extract,data_filter&(sex==1)), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=rbc, y=-..density..), data=subset(bd_extract,data_filter&(sex==0)) ,alpha = 0.9, binwidth=bin_width)
	rbc_sex_plot<-rbc_main_sex_plot+coord_flip() + theme_bw()+labs(title = sprintf("Density of %s in UK Biobank Stratified by Sex",toupper("rbc")))+ylab("Density")+xlab("RBC")+scale_fill_manual(name="", values=c("orange","blue"), labels=c("Females", "Males"))+scale_x_continuous(limits = c(4,5))

	#mpv	
	sort_diff=diff(sort(bd_extract$mpv))
	bin_width=min(sort_diff[sort_diff>0])/2
	mpv_main_sex_plot<-ggplot(subset(bd_extract, data_filter), aes(x=mpv,fill=as.factor(sex))) + geom_histogram(aes(x=mpv, y=..density..), data=subset(bd_extract,data_filter&(sex==1)), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=mpv, y=-..density..), data=subset(bd_extract,data_filter&(sex==0)) ,alpha = 0.9, binwidth=bin_width)
	mpv_sex_plot<-mpv_main_sex_plot+coord_flip() + theme_bw()+labs(title = sprintf("Density of %s in UK Biobank Stratified by Sex",toupper("mpv")))+ylab("Density")+xlab("MPV")+scale_fill_manual(name="", values=c("orange","blue"), labels=c("Females", "Males"))+scale_x_continuous(limits = c(7.5,11))

	#rdw	
	sort_diff=diff(sort(bd_extract$rdw))
	bin_width=min(sort_diff[sort_diff>0])/2
	rdw_main_sex_plot<-ggplot(subset(bd_extract, data_filter), aes(x=rdw,fill=as.factor(sex))) + geom_histogram(aes(x=rdw, y=..density..), data=subset(bd_extract,data_filter&(sex==1)), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=rdw, y=-..density..), data=subset(bd_extract,data_filter&(sex==0)) ,alpha = 0.9, binwidth=bin_width)
	rdw_sex_plot<-rdw_main_sex_plot+coord_flip() + theme_bw()+labs(title = sprintf("Density of %s in UK Biobank Stratified by Sex",toupper("rdw")))+ylab("Density")+xlab("RDW")+scale_fill_manual(name="", values=c("orange","blue"), labels=c("Females", "Males"))+scale_x_continuous(limits = c(12,14))

	pdf(sprintf("%s/output/graphics/hgb_sex_zoom.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(hgb_sex_plot)	
	dev.off()

	pdf(sprintf("%s/output/graphics/hct_sex_zoom.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(hct_sex_plot)	
	dev.off()

	pdf(sprintf("%s/output/graphics/rbc_sex_zoom.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(rbc_sex_plot)	
	dev.off()

	pdf(sprintf("%s/output/graphics/mpv_sex_zoom.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(mpv_sex_plot)	
	dev.off()

	pdf(sprintf("%s/output/graphics/rdw_sex_zoom.pdf",UKBIOBANK_FBC_REPO_DIR))	
	print(rdw_sex_plot)	
	dev.off()


}

generate_pairwise_gwas_correlation<-function(data_local, output_dir=UKBIOBANK_OUTPUT_DIR)
{
	library(reshape2)
	phenotypes=c("plt_gwas_normalised","mpv_gwas_normalised","pdw_gwas_normalised","pct_gwas_normalised","rbc_gwas_normalised","mcv_gwas_normalised","hct_gwas_normalised","mch_gwas_normalised","mchc_gwas_normalised","hgb_gwas_normalised","rdw_gwas_normalised","ret_gwas_normalised","ret_p_gwas_normalised","irf_gwas_normalised","hlr_gwas_normalised","hlr_p_gwas_normalised","mono_gwas_normalised","neut_gwas_normalised","eo_gwas_normalised","baso_gwas_normalised","neut_eo_sum_gwas_normalised","eo_baso_sum_gwas_normalised","baso_neut_sum_gwas_normalised","gran_gwas_normalised","neut_p_gran_gwas_normalised","baso_p_gran_gwas_normalised","eo_p_gran_gwas_normalised","myeloid_wbc_gwas_normalised","gran_p_myeloid_wbc_gwas_normalised","lymph_gwas_normalised","wbc_gwas_normalised","mono_p_gwas_normalised","neut_p_gwas_normalised","eo_p_gwas_normalised","baso_p_gwas_normalised","lymph_p_gwas_normalised")
	corr=cor(data_local[,phenotypes], use="pairwise")
	corr[upper.tri(corr)]=0
	melted_cor=melt(corr)

	ylabs=get_trait_long_name(phenotypes)
	xlabs=get_trait_short_name(phenotypes)
	names(ylabs)=phenotypes
	my_plot<-ggplot(data = melted_cor, aes(Var2, Var1, fill = value))+
		geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "red", mid = "white",  midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation")+theme_minimal(base_size=14)+theme(axis.text.x = element_text(angle = 45, vjust = 1,  size = 8, hjust = 1))+ coord_fixed()+
		ylab("Phenotype Description")+xlab("Phenotype Short Name")+scale_y_discrete(labels=ylabs)+scale_x_discrete(labels=xlabs)#+labs(title="Pairwise Correlations Between Quantile Normalised GWAS Phenotypes")
	for(i in 2:length(phenotypes))
	{
		my_plot=my_plot+geom_segment(x=i, xend=i, y=0.7, yend=i-0.7, linetype=2, size=0.003, alpha=1, colour="darkgrey")
	}
	pdf(sprintf("%s/gwas_corr.pdf", output_dir), width=10, height=10)
	print(my_plot)
	dev.off()

}

write_data_filter_spreadsheet<-function()
{
	require(xlsx)

	workbook=createWorkbook(type="xlsx")
	
	prep_data=bd_extract[!copied_data_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="copied data")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)
	
	prep_data=bd_extract[!zero_ret_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="zero ret")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!zero_iretp_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="zero iret %")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!low_ret_day_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="low ret day")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!eo_high_res_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="eo hi res")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!baso_high_res_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="baso hi res")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!mono_high_res_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="mono hi res")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!lymph_high_res_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="lymph hi res")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)

	prep_data=bd_extract[!neut_high_res_filter,c("f.eid","raw_accu_time")]
	sheet=createSheet(workbook, sheetName="neut hi res")
	addDataFrame(prep_data[!is.na(prep_data[,1]),],sheet, row.names=FALSE)
	saveWorkbook(workbook,sprintf("%s/output/blood_count_feedback.xlsx",UKBIOBANK_FBC_REPO_DIR))



}

# old QC code - generate UKB hist
	# ukbileve histogram
	#	bin_width=(max(data_local[,trait], na.rm=TRUE)-min(data_local[,trait],na.rm=TRUE))/100
	#	eval(parse(text=sprintf("main_ukbileve_plot<-ggplot(data_local, aes(x=%s,fill=as.factor(ukbileve_approx))) + geom_histogram(aes(x=%s, y=..density..), data=subset(data_local,ukbileve_approx), alpha = 0.9, binwidth=bin_width) + geom_histogram(aes(x=%s, y=-..density..), data=subset(data_local,!ukbileve_approx) ,alpha = 0.9, binwidth=bin_width)",trait,trait,trait)))
	#	ukbileve_plot<-main_ukbileve_plot+coord_flip() + theme_bw()+labs(title = sprintf("Density of %s in UK Biobank\n Stratified by Approximately Inferred UK BiLEVE Status",toupper(trait)))+plot_xlab+plot_dens_ylab+scale_fill_manual(name="", values=c("green","red"), labels=c("Not in UK BiLEVE","In UK BiLEVE"))
comment<-function()
{

	if(1==0)
	{
	if(covariate=="day_of_study_acq")
	{
		
		machine_service_data=read.table(sprintf("%s/machine_service.csv", UKBIOBANK_DATA_DIR), head=TRUE, stringsAsFactors=FALSE)
		machine_service_data$start_date_time=dmy(machine_service_data$Date, tz="Europe/London")+hm(machine_service_data$Start_Time)
		machine_service_data$end_date_time=dmy(machine_service_data$Date, tz="Europe/London")+hm(machine_service_data$End_Time)
		
		machine_service_data$day_of_study_start=round(difftime(floor_date(machine_service_data$start_date_time, unit="day"),first_day_of_study, unit="day"))+1
		machine_service_data$sec_of_study_start=difftime(machine_service_data$start_date_time,first_day_of_study, unit="secs")+1

		machine_service_data$day_of_study_end=round(difftime(floor_date(machine_service_data$end_date_time, unit="day"),first_day_of_study, unit="day"))+1
		machine_service_data$sec_of_study_end=difftime(machine_service_data$end_date_time,first_day_of_study, unit="secs")+1


		diff_days_acq=round(as.numeric(difftime(floor_date(bd$date_time_acq, unit="day"), floor_date(bd$date_time_acq[1], unit="day"), unit="days")))
		diff_start_date_time=round(as.numeric(difftime(floor_date(machine_service_data$start_date_time, unit="day"), floor_date(bd$date_time_acq[1], unit="day"), unit="days")))

		machine_service_data$start_day_of_study=machine_service_data$start_date_time-min(diff_days_acq,na.rm=TRUE)+1
		plot=plot+geom_vline(xintercept = as.numeric(machine_service_data$day_of_study_start[gsub("[A-Z]+(.[0-9]+)","\\1",instrument_id)==machine_service_data$machine]))

	}
	}

}
