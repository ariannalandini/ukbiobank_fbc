
get_trait_linear_units<-function(trait)
{
	if(length(trait)>1)
	{	
		return_vec=rep(NA, length(trait))
		for(i in 1:length(trait))
		{
			return_vec[i]=get_trait_linear_units(trait[i])
		}
		return(return_vec)
	}
	
	unit_str="NOUNITS"

	# extract data required so that we don't keep data_filtering at computational expense
	if(is.element(trait,c("irf")))
	{
		unit_str="1"

	}
	if(is.element(trait,c("baso","eo", "lymph", "mono","neut", "wbc","plt", "nrbc","plt_f","plt_i","myeloid_wbc", "gran","ipf_count")))
 	{
		unit_str="10^9, l^-1"
	}
	if(is.element(trait,c("hgb", "mchc")))
 	{
		unit_str="g,(dl)^-1"
	}
	if(is.element(trait,c("mpv", "mcv", "mrv", "mscv", "pdw", "rcvsd", "rdw_sd")))
	{
		unit_str="fl"

	}
	if(is.element(trait,c("rbc", "ret", "hlr")))
	{
		unit_str="10^{12}, l^-1"

	}
	if(is.element(trait,c("hct","pct","baso_p", "mono_p","neut_p","lymph_p","eo_p","ret_p", "nrbc_p", "hlr_p","gran_p_myeloid_wbc","baso_p_gran", "eo_p_gran","nrbc_p", "neut_p_gran", "ipf", "p_lcr", "rdw_cv","rdw_sd")))
	{
		unit_str="'%'"

	}
	if(is.element(trait,c("mch","ret_he")))
	{
		unit_str="pg"

	}
	if(is.element(trait,c("neut_ssc", "neut_sfl", "neut_fsc","mono_ssc", "mono_sfl", "mono_fsc","lymph_ssc", "lymph_sfl", "lymph_fsc")))
	{
		unit_str="ch"

	}

	return(unit_str)
}

ascii_units<-function(expr_string)
{	
	if(length(expr_string)>1)
	{	
		return_vec=rep(NA, length(expr_string))
		for(i in 1:length(expr_string))
		{
			return_vec[i]=ascii_units(expr_string[i])
		}
		return(return_vec)
	}
	return_str=expr_string
	if(expr_string=="10^9, l^-1") {return_str="per nl"}
	if(expr_string=="10^{12}, l^-1") {return_str="per pl"}
	if(expr_string=="g,(dl)^-1") {return_str="g per dl"}
	if(expr_string=="'%'") {return_str="%"}
	if(expr_string=="10^{-9}, l") {return_str="nl"}
	if(expr_string=="10^{-12}, l") {return_str="pl"}
	if(expr_string=="dlg^-1") {return_str="dl per g"}
	if(expr_string=="'%'^{-1}") {return_str="per %"}


	return(return_str)

}

get_inverse_trait_linear_units<-function(trait)
{

	if(length(trait)>1)
	{	
		return_vec=rep(NA, length(trait))
		for(i in 1:length(trait))
		{
			return_vec[i]=get_inverse_trait_linear_units(trait[i])
		}
		return(return_vec)
	}
	
	unit_str="NOUNITS"

	if(is.element(trait,c("irf")))
	{
		unit_str="1"

	}
	if(is.element(trait,c("nrbc","baso","eo", "lymph", "mono","neut", "wbc","plt","plt_f","plt_i", "nrbc", "myeloid_wbc", "gran","ipf_count")))
 	{
		unit_str="10^{-9}, l"
	}
	if(is.element(trait,c("hgb", "mchc")))
 	{
		unit_str="dlg^-1"
	}
	if(is.element(trait,c("mpv", "mcv", "mrv", "mscv", "pdw", "rcvsd", "rdw_sd")))
	{
		unit_str="(fl)^-1"

	}
	if(is.element(trait,c("rbc", "ret", "hlr")))
	{
		unit_str="10^{-12}, l"

	}
	if(is.element(trait,c("nrbc_p","hct","pct","baso_p", "mono_p","neut_p","rdw_cv","lymph_p","eo_p","ret_p", "nrbc_p", "hlr_p","gran_p_myeloid_wbc","baso_p_gran", "eo_p_gran", "neut_p_gran", "ipf", "p_lcr")))
	{
		unit_str="'%'^{-1}"

	}
	if(is.element(trait,c("mch","ret_he")))
	{
		unit_str="(pg)^{-1}"

	}
	if(is.element(trait,c("neut_ssc", "neut_sfl", "neut_fsc","mono_ssc", "mono_sfl", "mono_fsc","lymph_ssc", "lymph_sfl", "lymph_fsc")))
	{
		unit_str="(ch)^{-1}"

	}

	return(unit_str)

}
get_base_trait<-function(trait)
{
	gsub("(.+?)(_tech_adj|_gwas_adj|_tech_ntd_adj|_gwas_ntd_adj|_age_sex_meno_adj|_gwas_normalised)?", "\\1",trait)
	
}

get_trait_short_name<-function(trait)
{

	if(length(trait)>1)
	{	
		return_vec=rep(NA, length(trait))
		for(i in 1:length(trait))
		{
			return_vec[i]=get_trait_short_name(trait[i])
		}
		return(return_vec)
	}

	base_trait=get_base_trait(trait)

	if(base_trait=="plt"){ return("PLT#")}
	if(base_trait=="plt_i"){ return("PLT#-I")}
	if(base_trait=="plt_f"){ return("PLT#-F")}
	
	if(base_trait=="mpv"){ return("MPV")}
	if(base_trait=="pct"){ return("PCT")}
	if(base_trait=="pdw"){ return("PDW")}
	if(base_trait=="ipf"){ return("IPF%")}
	if(base_trait=="ipf_count"){ return("IPF#")}
	if(base_trait=="p_lcr"){ return("P-LCR")}


	if(base_trait=="rbc"){ return("RBC#")}
	if(base_trait=="mcv"){ return("MCV")}
	if(base_trait=="hct"){ return("HCT")}
	if(base_trait=="hgb"){ return("HGB")}
	if(base_trait=="mch"){ return("MCH")}
	if(base_trait=="mchc"){ return("MCHC")}
	if(base_trait=="rdw"){ return("RDW")}
	if(base_trait=="mscv"){ return("MSCV")}
	if(base_trait=="ret_he"){ return("RET-HE")}

	if(base_trait=="nrbc"){ return("NRBC#")}
	if(base_trait=="nrbc_p"){ return("NRBC%")}
	

	if(base_trait=="ret"){ return("RET#")}
	if(base_trait=="mrv"){ return("MRV")}
	if(base_trait=="ret_p"){ return("RET%")}
	if(base_trait=="irf"){ return("IRF")}
	if(base_trait=="hlr"){ return("HLSR#")}
	if(base_trait=="hlr_p"){ return("HLSR%")}
	if(base_trait=="rcvsd"){ return("RCVSD")}
	if(base_trait=="rdw_sd"){ return("RDW-SD")}
	if(base_trait=="rdw_cv"){ return("RDW-CV")}


	if(base_trait=="wbc"){ return("WBC#")}
	if(base_trait=="mono"){ return("MONO#")}
	if(base_trait=="eo"){ return("EO#")}
	if(base_trait=="neut"){ return("NEUT#")}
	if(base_trait=="baso"){ return("BASO#")}
	if(base_trait=="lymph"){ return("LYMPH#")}
	if(base_trait=="nrbc"){ return("NRBC#")}

	if(base_trait=="mono_p"){ return("MONO%")}
	if(base_trait=="eo_p"){ return("EO%")}
	if(base_trait=="neut_p"){ return("NEUT%")}
	if(base_trait=="baso_p"){ return("BASO%")}
	if(base_trait=="lymph_p"){ return("LYMPH%")}
	if(base_trait=="nrbc_p"){ return("NRBC%")}

	if(base_trait=="eo_p_gran"){ return("EO%GRAN")}
	if(base_trait=="neut_p_gran"){ return("NEUT%GRAN")}
	if(base_trait=="baso_p_gran"){ return("BASO%GRAN")}

	if(base_trait=="gran"){ return("GRAN#")}
	if(base_trait=="myeloid_wbc"){ return("MYELOID#")}
	if(base_trait=="gran_p_myeloid_wbc"){ return("GRAN%MYELOID")}

	if(base_trait=="neut_ssc"){ return("NEUT-SSC")}
	if(base_trait=="neut_fsc"){ return("NEUT-FSC")}
	if(base_trait=="neut_sfl"){ return("NEUT-SFL")}

	if(base_trait=="mono_ssc"){ return("MONO-SSC")}
	if(base_trait=="mono_fsc"){ return("MONO-FSC")}
	if(base_trait=="mono_sfl"){ return("MONO-SFL")}

	if(base_trait=="lymph_ssc"){ return("LYMPH-SSC")}
	if(base_trait=="lymph_fsc"){ return("LYMPH-FSC")}
	if(base_trait=="lymph_sfl"){ return("LYMPH-SFL")}

	
	return("No known name")
	
}

get_trait_cell_type<-function(trait)
{
	mature_red_cell_traits=c("hgb","rbc","mcv","mch","mchc","rdw_cv","rdw_sd","rcvsd","hct")
	immature_red_cell_traits=c("ret","hlr","ret_p","hlr_p","irf", "mrv", "ret_he")
	myeloid_white_cell_traits=c("neut","mono","baso","eo","myeloid_wbc","gran","gran_p_myeloid_wbc","eo_p_gran","neut_p_gran","baso_p_gran","neut_ssc", "neut_sfl", "neut_fsc","mono_ssc", "mono_sfl", "mono_fsc")
	lymphoid_white_cell_traits=c("lymph","lymph_ssc", "lymph_sfl", "lymph_fsc")
	compound_white_cell_traits=c("neut_p","eo_p","mono_p","lymph_p","baso_p","wbc")
	platelet_traits=c("mpv","plt","pdw","pct", "p_lcr", "ipf", "ipf_count", "plt_f","plt_i")
	if(length(trait)>1)
	{	
		return_vec=rep(NA, length(trait))
		for(i in 1:length(trait))
		{
			return_vec[i]=get_trait_type(trait[i])
		}
		return(return_vec)
	}
	base_trait=get_base_trait(trait)

	if(is.element(base_trait, mature_red_cell_traits)) {return("Mature red cell")}
	if(is.element(base_trait, immature_red_cell_traits)) {return("Immature red cell")}
	if(is.element(base_trait, myeloid_white_cell_traits)) {return("Myeloid white cell")}
	if(is.element(base_trait, lymphoid_white_cell_traits)) {return("Lymphoid white cell")}
	if(is.element(base_trait, compound_white_cell_traits)) {return("Compound white cell")}
	if(is.element(base_trait, platelet_traits)) {return("Platelet")}
	
	return("TRAIT TYPE UNKNOWN")	
}


get_trait_long_name<-function(trait)
{
	if(length(trait)>1)
	{	
		return_vec=rep(NA, length(trait))
		for(i in 1:length(trait))
		{
			return_vec[i]=get_trait_long_name(trait[i])
		}
		return(return_vec)
	}
	base_trait=get_base_trait(trait)

	if(base_trait=="plt"){ return("Platelet count")}
	if(base_trait=="plt_f"){ return("Platelet count (flow channel F)")}
	if(base_trait=="plt_i"){ return("Platelet count (Impedance channel)")}

	if(base_trait=="mpv"){ return("Mean platelet volume")}
	if(base_trait=="pct"){ return("Plateletcrit")}
	if(base_trait=="pdw"){ return("Platelet distribution width")}
	if(base_trait=="ipf"){ return("Immature percentage of platelets")}
	if(base_trait=="ipf_count"){ return("Immature platelet count")}
	if(base_trait=="p_lcr"){ return("Platelet large cell ratio")}


	if(base_trait=="rbc"){ return("Red blood cell count")}
	if(base_trait=="mcv"){ return("Mean corpuscular volume")}
	if(base_trait=="hct"){ return("Haematocrit")}
	if(base_trait=="hgb"){ return("Haemoglobin concentration")}
	if(base_trait=="mch"){ return("Mean corpuscular haemoglobin")}
	if(base_trait=="mchc"){ return("Mean corpuscular haemoglobin concentration")}
	if(base_trait=="rdw_cv"){ return("Red cell volume distribution coefficient of variation")}
	if(base_trait=="rdw_sd"){ return("Red cell volume distribution width at 20% of distribution height")}
	if(base_trait=="rcvsd"){ return("Red cell volume distribution standard deviation")}

	if(base_trait=="mscv"){ return("Mean spherical cell volume")}


	if(base_trait=="ret"){ return("Reticulocyte count")}
	if(base_trait=="mrv"){ return("Mean reticulocyte volume")}
	if(base_trait=="ret_p"){ return("Reticulocyte percentage of red cells")}
	if(base_trait=="irf"){ return("Immature fraction of reticulocytes")}
	if(base_trait=="hlr"){ return("High light scatter reticulocyte count")}
	if(base_trait=="hlr_p"){ return("High light scatter reticulocyte percentage of red cells")}
	if(base_trait=="ret_he"){ return("Reticulocyte haemoglobin equivalent")}

	if(base_trait=="wbc"){ return("White blood cell count")}
	if(base_trait=="mono"){ return("Monocyte count")}
	if(base_trait=="eo"){ return("Eosinophil count")}
	if(base_trait=="neut"){ return("Neutrophil count")}
	if(base_trait=="baso"){ return("Basophil count")}
	if(base_trait=="lymph"){ return("Lymphocyte count")}
	if(base_trait=="nrbc"){ return("Nucleated red blood cell count")}

	if(base_trait=="mono_p"){ return("Monocyte percentage of white cells")}
	if(base_trait=="eo_p"){ return("Eosinophil percentage of white cells")}
	if(base_trait=="neut_p"){ return("Neutrophil percentage of white cells")}
	if(base_trait=="baso_p"){ return("Basophil percentage of white cells")}
	if(base_trait=="lymph_p"){ return("Lymphocyte percentage of white cells")}
	if(base_trait=="nrbc_p"){return("Nucleated red blood cell count as percentage of white cells")}

	if(base_trait=="eo_p_gran"){ return("Eosinophil percentage of granulocytes")}
	if(base_trait=="neut_p_gran"){ return("Neutrophil percentage of granulocytes")}
	if(base_trait=="baso_p_gran"){ return("Basophil percentage of granulocytes")}

	if(base_trait=="gran"){ return("Granulocyte count")}

	if(base_trait=="myeloid_wbc"){ return("Myeloid white cell count")}

	if(base_trait=="gran_p_myeloid_wbc"){ return("Granulocyte percentage of myeloid white cells")}

	if(base_trait=="neut_ssc"){ return("Mean neutrophil side scatter")}
	if(base_trait=="neut_fsc"){ return("Mean neutrophil forward scatter")}
	if(base_trait=="neut_sfl"){ return("Mean neutrophil fluorescence")}


	if(base_trait=="mono_ssc"){ return("Mean monocyte side scatter")}
	if(base_trait=="mono_fsc"){ return("Mean monocyte forward scatter")}
	if(base_trait=="mono_sfl"){ return("Mean monocyte fluorescence")}

	if(base_trait=="lymph_ssc"){ return("Mean lymphocyte side scatter")}
	if(base_trait=="lymph_fsc"){ return("Mean lymphocyte forward scatter")}
	if(base_trait=="lymph_sfl"){ return("Mean lymphocyte fluorescence")}



	return("No known name")
}



tf_150kgwas<-function(traits=tf_null(), na.rm=TRUE)
{
	filter=c("plt","mpv","pdw","pct","rbc","mcv","hct","mch","mchc","hgb","rdw_cv","ret","ret_p","irf","hlr","hlr_p","mono","neut","eo","baso","neut_eo_sum","eo_baso_sum","baso_neut_sum","gran","neut_p_gran","baso_p_gran","eo_p_gran","myeloid_wbc","gran_p_myeloid_wbc","lymph","wbc","mono_p","neut_p","eo_p","baso_p","lymph_p")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_500kgwas<-function(traits=tf_null(), na.rm=TRUE)
{
	filter=c("plt","mpv","pdw","pct","rbc","mcv","hct","mch","mchc","hgb","rdw_cv","ret","ret_p","irf","mrv","hlr","hlr_p","mscv","mono","neut","eo","baso","lymph","wbc","mono_p","neut_p","eo_p","baso_p","lymph_p")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}
	

}

tf_01<-function(traits=tf_null(),na.rm=TRUE)
{
	filter=c("irf")	
	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_0100<-function(traits=tf_null(),na.rm=TRUE)
{
	filter=c("neut_p","eo_p","mono_p","lymph_p","baso_p","ret_p","hlr_p","hct","pct","gran_p_myeloid_wbc","eo_p_gran","neut_p_gran","baso_p_gran","p_lcr","ipf","nrbc_p")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}


tf_positive<-function(traits=tf_null(),na.rm=TRUE)
{
	filter=c("hgb","rbc","mcv","mch","mchc","rdw_cv","rdw_sd","rcvsd","ret","hlr","nrbc","mpv","plt","pdw","wbc","neut","mono","baso","eo","lymph","myeloid_wbc","gran","eo_baso_sum","neut_eo_sum","baso_neut_sum","plt_i","plt_f","ret_he","ipf_count","neut_ssc","neut_sfl","neut_fsc","mono_ssc","mono_sfl","mono_fsc","lymph_ssc","lymph_sfl","lymph_fsc","mrv","mscv")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}


tf_coulter<-function(traits=tf_null(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), get_base_trait(c(tf_coulter_measured(), tf_coulter_derived())))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_sysmex<-function(traits=tf_null(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), get_base_trait(c(tf_sysmex_measured(), tf_sysmex_derived())))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

	

}

tf_coulter_measured<-function(traits=tf_null(),na.rm=TRUE)
{
	filter=c("rbc", "hgb", "hct", "rcvsd","ret_p",  "mrv",  "mscv","hlr_p","nrbc_p", "plt", "pdw",  "pct","wbc", "neut_p", "eo_p", "mono_p", "lymph_p", "baso_p")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}
	

}

tf_coulter_derived<-function(traits=tf_null(),na.rm=TRUE)
{
	filter=c("mcv","rdw_cv", "mch", "mchc", "ret","irf", "hlr","nrbc","mpv", "neut", "mono", "baso", "eo","lymph", "myeloid_wbc", "gran","eo_baso_sum","neut_eo_sum","baso_neut_sum","gran_p_myeloid_wbc","eo_p_gran","neut_p_gran","baso_p_gran")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}
	

}

tf_sysmex_measured<-function(traits=tf_null(),na.rm=TRUE)
{
	filter=c("plt_f", "plt_i", "pct","pdw", "rbc", "hgb",  "hct", "irf", "ret_p", "rdw_sd","rcvsd", "nrbc_p","wbc",  "neut_p", "eo_p", "mono_p", "lymph_p", "baso_p", "p_lcr", "ipf_count","ret_he", "neut_ssc", "neut_sfl", "neut_fsc","mono_ssc", "mono_sfl", "mono_fsc","lymph_ssc", "lymph_sfl", "lymph_fsc")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}
	
}

tf_sysmex_derived<-function(traits=tf_null(),na.rm=TRUE)
{

	filter=c("plt", "mpv" ,"mcv", "mch", "mchc", "rdw_cv", "ret", "hlr", "hlr_p", "nrbc", "neut", "mono", "baso", "eo", "lymph", "myeloid_wbc", "gran","eo_baso_sum","neut_eo_sum","baso_neut_sum","gran_p_myeloid_wbc","eo_p_gran","neut_p_gran","baso_p_gran", "ipf")

	base_traits=get_base_trait(traits)
	traits[!is.element(base_traits, filter)]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}
	


}

tf_null<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_adjustable<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), setdiff(sort_traits_canonical(),c("nrbc","nrbc_p")))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}


tf_mature_red<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits),c("rbc","mcv","hct","mch","mchc","hgb","rdw_cv","rdw_sd","rcvsd", "mscv"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_immature_red<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c("ret","mrv","irf","hlr","nrbc"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}
tf_compound_red<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c("ret_p", "hlr_p"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}
tf_misc<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c("nrbc_p"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_red<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), get_base_trait(c(tf_mature_red(), tf_immature_red(), tf_compound_red())))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_platelet<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits),c("plt","mpv","pct","pdw"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}
tf_white<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), get_base_trait(c(tf_lymphocyte(), tf_granulocyte(), tf_compound_white(), tf_myeloid_white())))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_compound_white<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits),c("wbc","mono_p","neut_p","eo_p","baso_p","lymph_p"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}


}

tf_lymphocyte<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c("lymph"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_granulocyte<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c("neut","eo","baso","gran","neut_p_gran","baso_p_gran","eo_p_gran","gran"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

tf_myeloid_white<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c(tf_granulocyte(), "myeloid_wbc","gran_p_myeloid_wbc", "mono"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}
tf_bcx<-function(traits=sort_traits_canonical(),na.rm=TRUE)
{
	traits[!is.element(get_base_trait(traits), c("rbc","hgb","hct","mcv","mch","mchc","rdw_cv","wbc","baso", "eo","lymph","mono","neut","plt","mpv"))]=NA
	if(na.rm==TRUE)	{return(traits[!is.na(traits)])} else {return(traits)}

}

sort_traits_canonical<-function(traits)
{
	
	canonical_order=c("plt","mpv","pdw","pct","rbc","mcv","hct","mch","mchc","hgb","rdw_cv","rdw_sd","rcvsd","ret","ret_p","mrv","irf","hlr","hlr_p","mscv","nrbc","nrbc_p","mono","neut","eo","baso","gran","neut_p_gran","baso_p_gran","eo_p_gran","myeloid_wbc","gran_p_myeloid_wbc","lymph","wbc","mono_p","neut_p","eo_p","baso_p","lymph_p")

	if(!hasArg(traits))
	{
		traits=canonical_order
	}

	base_traits=get_base_trait(traits)

	traits[order(match(base_traits, canonical_order))]
}
sort_trait_types_canonical<-function(types)
{
	canonical_order=c("Platelet","Mature red cell","Immature red cell","Myeloid white cell","Lymphoid white cell", "Compound white cell")
	types[order(match(types, canonical_order))]
}





