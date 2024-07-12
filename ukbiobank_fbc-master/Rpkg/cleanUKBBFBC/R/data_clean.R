
## author: William John Astle
## author email: wja24@cam.ac.uk/will@astle.net
## Description: Prepare a cleaned up data set


logit<-function(x)
{

	x[(x<0)&(x>(-10^((-15))))]=0
	x[(x>1)&((1-x)>(-10^((-15))))]=1


	retval=log(x)
	retval=retval-log(1.0-x)
	retval
}

expit<-function(x)
{
	retval=x
	x_low=x<0
	x_high=x>=0
	x_low[is.na(x_low)]=FALSE
	x_high[is.na(x_high)]=FALSE

	retval[x_low]=exp(x[x_low]-log(1.0+exp(x[x_low])))
	retval[x_high]=exp(-log(1.0+exp(-x[x_high])))
	retval
}


quantile_normalise<-function(data, strata)
{
	levels=unique(strata)
	new_data=rep(NA, length(data))
	for(level in levels)	
	{	
		new_data_strata=qnorm((1:sum((!is.na(data[strata==level])))/(sum(!is.na(data[strata==level]))+1)))
		new_data[(!is.na(data))&(strata==level)]=new_data_strata[order(order(data[(!is.na(data))&(strata==level)]))]
	}
	return(new_data)
}



compute_derived_blood_count_variables<-function(data_local, postfix="", instrument="Coulter")
{
	# Red cells

	eval(parse(text=sprintf("return_data<-data.frame(mcv%s=10*data_local$hct%s/data_local$rbc%s)", postfix, postfix,postfix)))
	eval(parse(text=sprintf("return_data$mch%s<-10*data_local$hgb%s/data_local$rbc%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$mchc%s<-100*data_local$hgb%s/data_local$hct%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$ret%s<-data_local$ret_p%s*data_local$rbc%s/100", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$rdw_cv%s<-data_local$rcvsd%s/return_data$mcv%s*100", postfix, postfix, postfix)))

	if(instrument=="Coulter")
	{
		eval(parse(text=sprintf("return_data$hlr%s<-data_local$hlr_p%s*data_local$rbc%s/100", postfix, postfix, postfix)))
		eval(parse(text=sprintf("return_data$irf%s<-return_data$hlr%s/return_data$ret%s", postfix, postfix, postfix)))

	}
	if(instrument=="Sysmex")
	{
		eval(parse(text=sprintf("return_data$hlr%s<-data_local$irf%s*return_data$ret%s", postfix, postfix, postfix)))
		eval(parse(text=sprintf("return_data$hlr_p%s<-data_local$irf%s*return_data$ret%s*100/data_local$rbc%s", postfix, postfix, postfix, postfix)))

	}
	#platelets - replace dodgy f measurements with i measurements
	if((instrument=="Sysmex")&&((postfix=="")||(postfix=="_tech_adj")||(postfix=="_tech_ntd_adj")))
	{
			eval(parse(text=sprintf("return_data$plt%s<-data_local$plt_f%s",postfix,postfix)))
			eval(parse(text=sprintf("return_data$plt%s[is.na(data_local$plt_f%s)]<-data_local$plt_i%s[is.na(data_local$plt_f%s)]",postfix,postfix,postfix,postfix)))
	
			eval(parse(text=sprintf("return_data$ipf%s<-100*data_local$ipf_count%s/return_data$plt%s",postfix,postfix, postfix)))

	}
	else
	{
		eval(parse(text=sprintf("return_data$plt%s<-data_local$plt%s",postfix,postfix)))		
	}
	
	eval(parse(text=sprintf("return_data$mpv%s<-(data_local$pct%s/return_data$plt%s)*10000", postfix, postfix, postfix)))
		
	# white cells 
	eval(parse(text=sprintf("return_data$neut%s<-data_local$neut_p%s*data_local$wbc%s/100", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$baso%s<-data_local$baso_p%s*data_local$wbc%s/100", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$mono%s<-data_local$mono_p%s*data_local$wbc%s/100", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$eo%s<-data_local$eo_p%s*data_local$wbc%s/100", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$lymph%s<-data_local$lymph_p%s*data_local$wbc%s/100", postfix, postfix, postfix)))
	
	# non-clinical white cell
	eval(parse(text=sprintf("return_data$myeloid_wbc%s<-return_data$neut%s+return_data$eo%s+return_data$baso%s+return_data$mono%s", postfix, postfix,postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$gran%s<-return_data$neut%s+return_data$eo%s+return_data$baso%s", postfix,postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$eo_baso_sum%s<-return_data$eo%s+return_data$baso%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$neut_eo_sum%s<-return_data$neut%s+return_data$eo%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$baso_neut_sum%s<-return_data$baso%s+return_data$neut%s", postfix, postfix, postfix)))

	eval(parse(text=sprintf("return_data$gran_p_myeloid_wbc%s<-100*return_data$gran%s/return_data$myeloid_wbc%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$eo_p_gran%s<-100*return_data$eo%s/return_data$gran%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$neut_p_gran%s<-100*return_data$neut%s/return_data$gran%s", postfix, postfix, postfix)))
	eval(parse(text=sprintf("return_data$baso_p_gran%s<-100*return_data$baso%s/return_data$gran%s", postfix, postfix, postfix)))

	return_data			
}




adjust_variable_for_technical<-function(data_local,data_filters, trait, platform="Coulter")
{

	eval(parse(text=sprintf("data_local$trait=data_local$%s", trait)))

	
	if(is.element(get_base_trait(trait), tf_0100()))
	{
		data_local$adjust_scale_trait=logit(data_local$trait/100)
	}else
	{
	
		if(is.element(get_base_trait(trait), tf_01()))
		{
			data_local$adjust_scale_trait=logit(data_local$trait)

		}else
		{
			if(is.element(get_base_trait(trait), tf_positive()))
			{	
				data_local$adjust_scale_trait=log(data_local$trait)
	
			}
		}
	}

	# throw out data without an instrument
        data_local$adjust_scale_trait[is.na(data_local$instrument)]=NA

	# remove infinite values
	data_local$adjust_scale_trait[abs(data_local$adjust_scale_trait)==Inf]=NA
	
	data_local$adjust_scale_trait_drift_adjusted=adjust_variable_for_drift(data_local,data_filters)$adjusted_trait
	
	# identify outlying days
	outlier_day_data=data_local[,c("adjust_scale_trait_drift_adjusted", "day_of_study_acq","instrument")] %>% group_by(day_of_study_acq, instrument) %>% summarise(n=sum(!is.na(adjust_scale_trait_drift_adjusted)),mean=mean(adjust_scale_trait_drift_adjusted,na.rm=TRUE))
	outlier_day_data$z_score=sqrt(outlier_day_data$n)*abs(outlier_day_data$mean-median(data_local$adjust_scale_trait_drift_adjusted,na.rm=TRUE))/mad(data_local$adjust_scale_trait_drift_adjusted,na.rm=TRUE)

	bad_days_data=na.omit(outlier_day_data[(outlier_day_data$n<10)|(outlier_day_data$z_score>6.5),])
	bad_data=is.element(sprintf("%s:%d",data_local$instrument, data_local$day_of_study_acq),sprintf("%s:%d", bad_days_data$instrument, bad_days_data$day_of_study_acq))
	bad_data[is.na(bad_data)]=FALSE
	print(sprintf("%s %f bad data",trait, sum(bad_data,na.rm=TRUE)))
	data_local$adjust_scale_trait[bad_data]=NA

	drift_adjusted_ret_list=adjust_variable_for_drift(data_local,data_filters)
	
	if(is.element(get_base_trait(trait), tf_0100()))
	{	
		drift_adjusted_ret_list$adjusted_trait=100*expit(drift_adjusted_ret_list$adjusted_trait)
	}
	if(is.element(get_base_trait(trait), tf_01()))
	{	
		drift_adjusted_ret_list$adjusted_trait=expit(drift_adjusted_ret_list$adjusted_trait)
	}

	if(is.element(get_base_trait(trait),tf_positive()))
	{	
		drift_adjusted_ret_list$adjusted_trait=exp(drift_adjusted_ret_list$adjusted_trait)		
	}
	drift_adjusted_ret_list
}


adjust_variable_for_drift<-function(data_local, data_filters)
{	
	return_list=list()


	data_local$change_class=data_local$instrument
	prediction_data=data_local[,c("subject_id", "measurement_id", "adjust_scale_trait", "trait","second_acq", "day_of_year_ext","delay_was_imputed")]

	prediction_data$vene_to_midnight_secs=as.numeric(data_local$imputed_delay_mins)*60-as.numeric(data_local$second_acq)
	
	prediction_data$second_of_study_acq=as.numeric(data_local$second_of_study_acq)
	prediction_data$instrument=as.factor(data_local$instrument)
	prediction_data$day_of_week_acq=as.factor(data_local$day_of_week_acq)
	
	# now use all the data
	CENTRAL_SDS=Inf	
	prediction_data$central_data=abs(prediction_data$adjust_scale_trait-median(prediction_data$adjust_scale_trait,na.rm=TRUE))<CENTRAL_SDS*mad(prediction_data$adjust_scale_trait,na.rm=TRUE)

	prediction_data_nona=na.omit(prediction_data[data_filters$delay_good&data_filters$unique_machine_time_stamp,])

	rownames(prediction_data_nona)=prediction_data_nona$measurement_id

	print("starting gam fit")

	# the circular term used to be bs="cp" but this started crashing so switched to cubic spline # first k was 50
	gam_out=gam(adjust_scale_trait~s(second_of_study_acq,by=as.factor(instrument),k=80, bs="ps")+s(day_of_year_ext, k=30, bs="cc")+s(second_acq, vene_to_midnight_secs, by=interaction(instrument,delay_was_imputed),k=30, bs="tp")+as.factor(instrument)+as.factor(day_of_week_acq), data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
	print("ending gam fit")
	
	print(sprintf("Explained: %f percent of the variance",100*summary(gam_out)$r.sq))

	predicted_vals=predict(gam_out, type="response", newdata=prediction_data_nona, na.action=na.pass)
	adjusted_trait=as.vector(data_local$adjust_scale_trait-predicted_vals[match(data_local$measurement_id,rownames(predicted_vals))]+mean(prediction_data_nona$adjust_scale_trait[prediction_data_nona$central_data], na.rm=TRUE))
	
	return_list$adjusted_trait=adjusted_trait
	return_list$r.sq=summary(gam_out)$r.sq
	return_list$r.sq_full_data=1-var(data_local$adjust_scale_trait-predicted_vals[match(data_local$measurement_id,rownames(predicted_vals))], na.rm=TRUE)/var(data_local$adjust_scale_trait,na.rm=TRUE)

	return(return_list)

}


adjust_variable_for_gwas<-function(data_local, data_filters, trait, study="UK Biobank")
{

	
	if(study=="UK Biobank")
	{
		platform="Coulter"
	}
	if(study=="INTERVAL")
	{
		platform="Sysmex"
	}
	eval(parse(text=sprintf("data_local$trait=data_local$%s", trait)))
	
	if(is.element(get_base_trait(trait), tf_0100()))
	{
		data_local$adjust_scale_trait=logit(data_local$trait/100)
	}else
	{
	
		if(is.element(get_base_trait(trait), tf_01()))
		{
			data_local$adjust_scale_trait=logit(data_local$trait)

		}else
		{
			if(is.element(get_base_trait(trait), tf_positive()))
			{	
				data_local$adjust_scale_trait=log(data_local$trait)
	
			}
		}
	}


	print(trait)

	if(study=="UK Biobank")
	{
		prediction_data=data_local[,c("subject_id", "measurement_id", "adjust_scale_trait", "trait","age_acq", "meno", "height", "weight", "pack_years_smoked","smoking_status", "smoking_amount", "drinking_status","alcohol_intake","alcohol_drunk_yesterday","days_since_period")]


		prediction_data$selfr_diabetes=as.numeric(data_filters$selfr_diabetes)
		prediction_data$selfr_typeI_diabetes=as.numeric(data_filters$selfr_typeI_diabetes)
		prediction_data$selfr_hyper_tension=as.numeric(data_filters$selfr_hyper_tension)
		prediction_data$selfr_typeII_diabetes=as.numeric(data_filters$selfr_typeII_diabetes)
		prediction_data$selfr_asthma=as.numeric(data_filters$selfr_asthma)
		prediction_data$dr_diag_diabetes=as.numeric(data_filters$dr_diag_diabetes)
		prediction_data$selfr_emphysema=as.numeric(data_filters$selfr_emphysema)
		prediction_data$adjust_scale_trait[apply(data_filters[,c(ukbb_filter_lists$blood_cancer_filters, ukbb_filter_lists$non_cancer_blood_condition_filters)],1, any)]=NA
		prediction_data$alcohol_drunk_yesterday_missing=as.numeric(is.na(prediction_data$alcohol_drunk_yesterday))
		prediction_data$days_since_period_missing=as.numeric(is.na(prediction_data$days_since_period))
		
		# '1' is an arbitrary value since we are putting dummy adjustments in
		prediction_data$alcohol_drunk_yesterday[is.na(prediction_data$alcohol_drunk_yesterday)]=1
		prediction_data$days_since_period[is.na(prediction_data$days_since_period)]=1
	}
	if(study=="INTERVAL")
	{
		dup_subs=unique(data_local$subject_id[duplicated(data_local$subject_id)])
		for(sub in dup_subs)
		{
			sub_records=(data_local$subject_id==sub)&(!is.na(data_local$adjust_scale_trait))
			if(sum(sub_records)<2){next;}
			data_local$adjust_scale_trait[sub_records&(data_local$date_time_acq!=min(data_local$date_time_acq[sub_records],na.rm=TRUE))]=NA
		}
		
		prediction_data=data_local[,c("subject_id", "measurement_id", "adjust_scale_trait", "trait","age_acq", "meno", "height", "weight", "pack_years_smoked","smoking_status", "smoking_amount", "drinking_status","alcohol_intake")]
		#remove subject specific duplicates
	}
	
	prediction_data$central_data=abs(prediction_data$adjust_scale_trait-median(prediction_data$adjust_scale_trait,na.rm=TRUE))<CENTRAL_SDS*mad(prediction_data$adjust_scale_trait,na.rm=TRUE)

	prediction_data$height_missing=as.numeric(is.na(prediction_data$height))
	prediction_data$weight_missing=as.numeric(is.na(prediction_data$weight))
	prediction_data$age_acq_missing=as.numeric(is.na(prediction_data$age_acq))
	prediction_data$pack_years_smoked_missing=as.numeric(is.na(prediction_data$pack_years_smoked))

	# this is an arbitrary value since we are putting dummy adjustments in
	prediction_data$height[is.na(prediction_data$height)]=1
	prediction_data$weight[is.na(prediction_data$weight)]=1
	prediction_data$log_height=log(prediction_data$height)
	prediction_data$log_weight=log(prediction_data$weight)
	prediction_data$age_acq[is.na(prediction_data$age_acq)]=1
	prediction_data$pack_years_smoked[is.na(prediction_data$pack_years_smoked)]=1

	prediction_data$meno[is.na(prediction_data$meno)]="data_missing"
	prediction_data$smoking_status[is.na(prediction_data$smoking_status)]="data_missing"
	prediction_data$smoking_amount[is.na(prediction_data$smoking_amount)]="data_missing"
	prediction_data$drinking_status[is.na(prediction_data$drinking_status)]="data_missing"
	prediction_data$alcohol_intake[is.na(prediction_data$alcohol_intake)]="data_missing"
	prediction_data$meno_collapse=prediction_data$meno
	prediction_data$meno_collapse[prediction_data$meno_collapse=="unsure"|prediction_data$meno_collapse=="hyst"]="data_missing"

 	prediction_data_nona=na.omit(prediction_data)
	rownames(prediction_data_nona)=prediction_data_nona$measurement_id

	if(study=="UK Biobank")
	{
		gam_out=gam(adjust_scale_trait~s(age_acq,by=as.factor(meno), k=30, bs="ps")+s(log_weight, log_height, by=as.factor(meno), k=30, bs="tp")+
			s(days_since_period,k=30, bs="ps")+s(alcohol_drunk_yesterday, k=30, bs="ps")+as.factor(drinking_status)+as.factor(alcohol_intake)+
			s(pack_years_smoked, k=30, bs="ps")+as.factor(smoking_status)+as.factor(smoking_amount)+
			days_since_period_missing+age_acq_missing+alcohol_drunk_yesterday_missing+pack_years_smoked_missing+height_missing+weight_missing+
			selfr_diabetes+selfr_typeI_diabetes+selfr_hyper_tension+selfr_typeII_diabetes+selfr_asthma+dr_diag_diabetes+selfr_emphysema, 
			data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
	}
	if(study=="INTERVAL")
	{
		gam_out=gam(adjust_scale_trait~s(age_acq,by=as.factor(meno), k=30, bs="ps")+s(log_weight, log_height, by=as.factor(meno), k=30, bs="tp")+
			as.factor(drinking_status)+as.factor(alcohol_intake)+s(pack_years_smoked, k=30, bs="ps")+as.factor(smoking_status)+as.factor(smoking_amount)+
			age_acq_missing+pack_years_smoked_missing+height_missing+weight_missing,
			data=prediction_data_nona[prediction_data_nona$central_data,], optimizer=c("outer", "newton"))
	}
	print(sprintf("Explained: %f percent of the variance",100*summary(gam_out)$r.sq))
	predicted_vals=predict(gam_out, type="response", newdata=prediction_data_nona, na.action=na.pass)

	prediction_data$adjusted_trait=as.vector(data_local$adjust_scale_trait-predicted_vals[match(data_local$measurement_id,rownames(predicted_vals))]+mean(prediction_data_nona$adjust_scale_trait[prediction_data_nona$central_data], na.rm=TRUE))

	return_list=list()
		
	if(is.element(get_base_trait(trait), tf_0100()))
	{	
		return_list$adjusted_trait=100*expit(prediction_data$adjusted_trait)
	}
	if(is.element(get_base_trait(trait), tf_01()))
	{	
		return_list$adjusted_trait=expit(prediction_data$adjusted_trait)
	}

	if(is.element(get_base_trait(trait),tf_positive()))
	{	
		return_list$adjusted_trait=exp(prediction_data$adjusted_trait)		
	}

	return_list$r.sq=summary(gam_out)$r.sq
	return_list$r.sq_full_data=1-var(data_local$adjust_scale_trait-predicted_vals[match(data_local$measurement_id,rownames(predicted_vals))], na.rm=TRUE)/var(data_local$adjust_scale_trait,na.rm=TRUE)

	return(return_list)
}

identify_max_rsq_for_exclusion_variables<-function()
{
	for(name in names(ukbb_filter_lists))
	{
		max_rsq=0
		max_trait=ukbb_traits[1]
		for(trait in ukbb_traits)
		{
			outcome=quantile_normalise(bd_extract[,trait])
			predictors=apply(ukbb_data_filters[,is.element(colnames(ukbb_data_filters),ukbb_filter_lists[[name]])],2, as.numeric)
			lm_obj=lm(outcome~predictors)
			if(max_rsq<summary(lm_obj)$adj.r.squared)
			{
				max_rsq=summary(lm_obj)$adj.r.squared
				max_trait=trait
			}
		}
		print(name)
		print(max_trait)
		print(100*max_rsq)
		print(max_rsq*140000)
	}

	for(name in ukbb_filter_lists[[8]])
	{
		max_rsq=0
		max_trait=ukbb_traits[1]
		for(trait in ukbb_traits)
		{
			outcome=quantile_normalise(bd_extract[,trait])
			predictors=apply(ukbb_data_filters[,name,drop=FALSE],2, as.numeric)
			lm_obj=lm(outcome~predictors)
			if(max_rsq<summary(lm_obj)$adj.r.squared)
			{
				max_rsq=summary(lm_obj)$adj.r.squared
				max_trait=trait
			}
		}
		print(name)
		print(max_trait)
		print(100*max_rsq)
		print(max_rsq*140000)

	}

}

#EOF

