
## author: William John Astle
## author email: wja24@cam.ac.uk/will@astle.net
## Description: Prepare a cleaned up data set


run_all_UKBB_prep<-function()
{
	

load_merged_UKBB_data()
prepare_UKBB_data()
load_prepared_UKBB_data()
load_ancillary_UKBB_data()
adjust_UKBB_variables_for_technical()
load_tech_adjusted_UKBB_data()
adjust_UKBB_variables_for_gwas()
load_gwas_adjusted_UKBB_data()
prepare_UKBB_traits_for_gwas_and_save()
load_gwas_prepared_UKBB_data()

}

prepare_UKBB_data<-function()
{

	## Setup variables
	## Ids: required to make compatible with INTERVAL QC scripts
	bd$measurement_id=1:dim(bd)[1]
	bd$subject_id=bd$f.eid

	# clinic
	bd$clinic=NA
	bd$clinic[bd$f.54.0.0==11012]="Barts"
	bd$clinic[bd$f.54.0.0==11021]="Birmingham"
	bd$clinic[bd$f.54.0.0==11011]="Bristol"
	bd$clinic[bd$f.54.0.0==11008]="Bury"
	bd$clinic[bd$f.54.0.0==11003]="Cardiff"
	bd$clinic[bd$f.54.0.0==11024]="Cheadle (revisit)"
	bd$clinic[bd$f.54.0.0==11020]="Croydon"
	bd$clinic[bd$f.54.0.0==11005]="Edinburgh"
	bd$clinic[bd$f.54.0.0==11004]="Glasgow"
	bd$clinic[bd$f.54.0.0==11018]="Hounslow"
	bd$clinic[bd$f.54.0.0==11010]="Leeds"
	bd$clinic[bd$f.54.0.0==11016]="Liverpool"
	bd$clinic[bd$f.54.0.0==11001]="Manchester"
	bd$clinic[bd$f.54.0.0==11017]="Middlesborough"
	bd$clinic[bd$f.54.0.0==11009]="Newcastle"
	bd$clinic[bd$f.54.0.0==11013]="Nottingham"
	bd$clinic[bd$f.54.0.0==11002]="Oxford"
	bd$clinic[bd$f.54.0.0==11007]="Reading"
	bd$clinic[bd$f.54.0.0==11014]="Sheffield"
	bd$clinic[bd$f.54.0.0==11003]="Stockport (pilot)"
	bd$clinic[bd$f.54.0.0==11006]="Stoke"
	bd$clinic[bd$f.54.0.0==11022]="Swansea"
	bd$clinic[bd$f.54.0.0==11023]="Wrexham"
	bd$clinic[bd$f.54.0.0==11025]="Cheadle (imaging)"




	
	## Blood traits 
	## platelet
	bd$plt=bd$f.30080.0.0
	bd$mpv=bd$f.30100.0.0
 	bd$pct=bd$f.30090.0.0
       	bd$pdw=bd$f.30110.0.0

	## mature red
	bd$hct=bd$f.30030.0.0
	bd$rbc=bd$f.30010.0.0
	bd$mcv=bd$f.30040.0.0
        bd$mch=bd$f.30050.0.0
	bd$mchc=bd$f.30060.0.0
	bd$hgb=bd$f.30020.0.0
  	bd$rdw_cv=bd$f.30070.0.0
	# red cell volume standard deviation for cv adjustment
	bd$rcvsd=bd$rdw_cv*bd$mcv/100


	## immature red
  	bd$ret=bd$f.30250.0.0
  	bd$ret_p=bd$f.30240.0.0
  	bd$irf=bd$f.30280.0.0
  	bd$mrv=bd$f.30260.0.0
  	bd$mscv=bd$f.30270.0.0
	bd$hlr=bd$f.30300.0.0
	bd$hlr_p=bd$f.30290.0.0
	bd$nrbc=bd$f.30170.0.0
	bd$nrbc_p=bd$f.30230.0.0

	## white cell
	bd$baso=bd$f.30160.0.0
  	bd$eo=bd$f.30150.0.0
  	bd$lymph=bd$f.30120.0.0
  	bd$mono=bd$f.30130.0.0
  	bd$neut=bd$f.30140.0.0
 	bd$wbc=bd$f.30000.0.0
	bd$baso_p=bd$f.30220.0.0
  	bd$eo_p=bd$f.30210.0.0
  	bd$lymph_p=bd$f.30180.0.0
  	bd$mono_p=bd$f.30190.0.0
  	bd$neut_p=bd$f.30200.0.0
 
	## route
	#platelet
	bd$plt_route=bd$f.30084.0.0
	bd$mpv_route=bd$f.30104.0.0
 	bd$pct_route=bd$f.30094.0.0
       	bd$pdw_route=bd$f.30114.0.0

	## mature red
	bd$hct_route=bd$f.30034.0.0
	bd$rbc_route=bd$f.30014.0.0
	bd$mcv_route=bd$f.30044.0.0
        bd$mch_route=bd$f.30054.0.0
	bd$mchc_route=bd$f.30064.0.0
	bd$hgb_route=bd$f.30024.0.0
  	bd$rdw_cv_route=bd$f.30074.0.0

	## immature red
  	bd$ret_route=bd$f.30254.0.0
  	bd$ret_p_route=bd$f.30244.0.0
  	bd$irf_route=bd$f.30284.0.0
  	bd$mrv_route=bd$f.30264.0.0
  	bd$mscv_route=bd$f.30274.0.0
	bd$hlr_route=bd$f.30304.0.0
	bd$hlr_p_route=bd$f.30294.0.0
	bd$nrbc_route=bd$f.30174.0.0
	bd$nrbc_p_route=bd$f.30234.0.0

	## white cell
	bd$baso_route=bd$f.30164.0.0
  	bd$eo_route=bd$f.30154.0.0
  	bd$lymph_route=bd$f.30124.0.0
  	bd$mono_route=bd$f.30134.0.0
  	bd$neut_route=bd$f.30144.0.0
 	bd$wbc_route=bd$f.30004.0.0
	bd$baso_p_route=bd$f.30224.0.0
  	bd$eo_p_route=bd$f.30214.0.0
  	bd$lymph_p_route=bd$f.30184.0.0
  	bd$mono_p_route=bd$f.30194.0.0
  	bd$neut_p_route=bd$f.30204.0.0

	## instrument
	#platelet
	bd$plt_inst=bd$f.30083.0.0
	bd$mpv_inst=bd$f.30103.0.0
 	bd$pct_inst=bd$f.30093.0.0
       	bd$pdw_inst=bd$f.30113.0.0

	## mature red
	bd$hct_inst=bd$f.30033.0.0
	bd$rbc_inst=bd$f.30013.0.0
	bd$mcv_inst=bd$f.30043.0.0
        bd$mch_inst=bd$f.30053.0.0
	bd$mchc_inst=bd$f.30063.0.0
	bd$hgb_inst=bd$f.30023.0.0
  	bd$rdw_cv_inst=bd$f.30073.0.0

	## immature red
  	bd$ret_inst=bd$f.30253.0.0
  	bd$ret_p_inst=bd$f.30243.0.0
  	bd$irf_inst=bd$f.30283.0.0
  	bd$mrv_inst=bd$f.30263.0.0
  	bd$mscv_inst=bd$f.30273.0.0
	bd$hlr_inst=bd$f.30303.0.0
	bd$hlr_p_inst=bd$f.30293.0.0
	bd$nrbc_inst=bd$f.30173.0.0
	bd$nrbc_p_inst=bd$f.30233.0.0

	## white cell
	bd$baso_inst=bd$f.30163.0.0
  	bd$eo_inst=bd$f.30153.0.0
  	bd$lymph_inst=bd$f.30123.0.0
  	bd$mono_inst=bd$f.30133.0.0
  	bd$neut_inst=bd$f.30143.0.0
 	bd$wbc_inst=bd$f.30003.0.0
	bd$baso_p_inst=bd$f.30223.0.0
  	bd$eo_p_inst=bd$f.30213.0.0
  	bd$lymph_p_inst=bd$f.30183.0.0
  	bd$mono_p_inst=bd$f.30193.0.0
  	bd$neut_p_inst=bd$f.30203.0.0

	## date time

	print("needs unhashing")
#	bd$date_time_ext=parse_date_time(bd$f.3166.0.0, "dmy HMS",tz="Europe/London")


	#platelet
	bd$plt_date_time_acq=parse_date_time(bd$f.30082.0.0,"ymd HMS", tz="Europe/London")
	bd$mpv_date_time_acq=parse_date_time(bd$f.30102.0.0,"ymd HMS", tz="Europe/London")
 	bd$pct_date_time_acq=parse_date_time(bd$f.30092.0.0,"ymd HMS", tz="Europe/London")
       	bd$pdw_date_time_acq=parse_date_time(bd$f.30112.0.0,"ymd HMS", tz="Europe/London")

	## mature red
	bd$hct_date_time_acq=parse_date_time(bd$f.30032.0.0, "ymd HMS", tz="Europe/London")
	bd$rbc_date_time_acq=parse_date_time(bd$f.30012.0.0, "ymd HMS", tz="Europe/London")
	bd$mcv_date_time_acq=parse_date_time(bd$f.30042.0.0, "ymd HMS", tz="Europe/London")
        bd$mch_date_time_acq=parse_date_time(bd$f.30052.0.0, "ymd HMS", tz="Europe/London")
	bd$mchc_date_time_acq=parse_date_time(bd$f.30062.0.0, "ymd HMS", tz="Europe/London")
	bd$hgb_date_time_acq=parse_date_time(bd$f.30022.0.0, "ymd HMS", tz="Europe/London")
  	bd$rdw_cv_date_time_acq=parse_date_time(bd$f.30072.0.0, "ymd HMS", tz="Europe/London")

	## immature red
  	bd$ret_date_time_acq=parse_date_time(bd$f.30252.0.0, "ymd HMS", tz="Europe/London")
  	bd$ret_p_date_time_acq=parse_date_time(bd$f.30242.0.0, "ymd HMS", tz="Europe/London")
  	bd$irf_date_time_acq=parse_date_time(bd$f.30282.0.0, "ymd HMS", tz="Europe/London")
  	bd$mrv_date_time_acq=parse_date_time(bd$f.30262.0.0, "ymd HMS", tz="Europe/London")
  	bd$mscv_date_time_acq=parse_date_time(bd$f.30272.0.0, "ymd HMS", tz="Europe/London")
	bd$hlr_date_time_acq=parse_date_time(bd$f.30302.0.0, "ymd HMS", tz="Europe/London")
	bd$hlr_p_date_time_acq=parse_date_time(bd$f.30292.0.0, "ymd HMS", tz="Europe/London")
	bd$nrbc_date_time_acq=parse_date_time(bd$f.30172.0.0, "ymd HMS", tz="Europe/London")
	bd$nrbc_p_date_time_acq=parse_date_time(bd$f.30232.0.0, "ymd HMS", tz="Europe/London")

	## white cell
	bd$baso_date_time_acq=parse_date_time(bd$f.30162.0.0, "ymd HMS", tz="Europe/London")
  	bd$eo_date_time_acq=parse_date_time(bd$f.30152.0.0, "ymd HMS", tz="Europe/London")
  	bd$lymph_date_time_acq=parse_date_time(bd$f.30122.0.0, "ymd HMS", tz="Europe/London")
  	bd$mono_date_time_acq=parse_date_time(bd$f.30132.0.0, "ymd HMS", tz="Europe/London")
  	bd$neut_date_time_acq=parse_date_time(bd$f.30142.0.0, "ymd HMS", tz="Europe/London")
 	bd$wbc_date_time_acq=parse_date_time(bd$f.30002.0.0, "ymd HMS", tz="Europe/London")
	bd$baso_p_date_time_acq=parse_date_time(bd$f.30222.0.0, "ymd HMS", tz="Europe/London")
  	bd$eo_p_date_time_acq=parse_date_time(bd$f.30212.0.0, "ymd HMS", tz="Europe/London")
  	bd$lymph_p_date_time_acq=parse_date_time(bd$f.30182.0.0, "ymd HMS", tz="Europe/London")
  	bd$mono_p_date_time_acq=parse_date_time(bd$f.30192.0.0, "ymd HMS", tz="Europe/London")
  	bd$neut_p_date_time_acq=parse_date_time(bd$f.30202.0.0, "ymd HMS", tz="Europe/London")

	# consensus
	bd$route=apply(bd[,grep("_route",colnames(bd))], 1, function(x){ret=unique(x[!is.na(x)]); if(length(ret)==1){return(ret)}; return(NA)})
	bd$instrument=apply(bd[,grep("_inst",colnames(bd))], 1, function(x){ret=unique(x[!is.na(x)]); if(length(ret)==1){return(ret)}; return(NA)})
	bd$date_time_acq=parse_date_time(apply(bd[,grep("_date_time_acq",colnames(bd))], 1, function(x){ret=unique(x[!is.na(x)]); if(length(ret)==1){return(ret)}; return(NA)}),"ymd HMS", tz="Europe/London")
	
	bd$match=sprintf("%s:%s", bd$instrument, bd$date_time_acq)
	bd$match[grepl("NA", bd$match)]=NA

	## covariates

	## anthropometric 
	bd$weight=bd$f.23098.0.0
	bd$height=bd$f.50.0.0
	bd$bmi=bd$f.21001.0.0

	## smoking
	bd$pack_years_smoked=bd$f.20161.0.0
	bd$smoking_status=NA
	bd$smoking_status[bd$f.20116.0.0==0]="never"
	bd$smoking_status[bd$f.20116.0.0==1]="previous"
	bd$smoking_status[bd$f.20116.0.0==2]="current"
	bd$smoking_status[bd$f.20116.0.0==-3]="no_answer"

	bd$smoking_amount=NA
	bd$smoking_amount[bd$f.1239.0.0==0]="never"
	bd$smoking_amount[bd$f.1239.0.0==1]="most_days"
	bd$smoking_amount[bd$f.1239.0.0==2]="occasional"
	bd$smoking_amount[bd$f.1239.0.0==-3]="no_answer"

	## alcohol
	bd$drinking_status=NA
	bd$drinking_status[bd$f.20117.0.0==0]="never"
	bd$drinking_status[bd$f.20117.0.0==1]="previous"
	bd$drinking_status[bd$f.20117.0.0==2]="current"
	bd$drinking_status[bd$f.20117.0.0==-3]="no_answer"

	bd$alcohol_intake=NA
	bd$alcohol_intake[bd$f.1558.0.0==1]="most_days"
	bd$alcohol_intake[bd$f.1558.0.0==2]="three_four_weekly"
	bd$alcohol_intake[bd$f.1558.0.0==3]="one_two_weekly"
	bd$alcohol_intake[bd$f.1558.0.0==4]="one_to_three_monthly"
	bd$alcohol_intake[bd$f.1558.0.0==5]="special_occasions"
	bd$alcohol_intake[bd$f.1558.0.0==6]="never"
	bd$alcohol_intake[bd$f.1558.0.0==-3]="no_answer"

	bd$alcohol_drunk_yesterday=bd$f.100022.0.0

	## genetic

	bd$selfr_sex=NA	
	bd$selfr_sex=NA
	bd$selfr_sex[bd$f.31.0.0==0]="female"
	bd$selfr_sex[bd$f.31.0.0==1]="male"

	bd$sex=NA
	bd$sex[(bd$selfr_sex=="female")&(bd$genetic_sex=="female")]="female"
	bd$sex[(bd$selfr_sex=="male")&(bd$genetic_sex=="male")]="male"

	## female specific 
	bd$meno_age=bd$f.3581.0.0
	bd$menarche_age=bd$f.2714.0.0
	bd$menst_today=bd$f.3720.0.0
	bd$menst_cycle=bd$f.3710.0.0
	bd$days_since_period=bd$f.3700.0.0	
	
	bd$meno=NA
	bd$meno[bd$f.2724.0.0==1]="post"
	bd$meno[bd$f.2724.0.0==0]="pre"
	bd$meno[bd$f.2724.0.0==2]="hyst"
	bd$meno[bd$f.2724.0.0==3]="unsure"
	bd$meno[bd$f.2724.0.0==4]="no_answer"
	bd$meno[bd$sex=="male"]="male"
	bd$meno[is.na(bd$sex)]=NA

	bd$meno_simple=bd$meno
	bd$meno_simple[is.element(bd$meno_simple,c("unsure"))]=NA
	
	## age

	bd$yob=bd$f.34.0.0
	bd$mob=bd$f.52.0.0
	
	## technical

#	# aquisition time info	

	negative_delay=difftime(bd$date_time_acq, bd$date_time_ext)<0
	negative_delay[is.na(negative_delay)]=FALSE
	
	bd$date_time_acq[negative_delay]=NA
	bd$date_time_ext[negative_delay]=NA
	
	same_day=floor_date(bd$date_time_acq,"day")==floor_date(bd$date_time_ext,"day")
	same_day[is.na(same_day)]=FALSE
	

	bd$date_time_acq[same_day]=NA
	bd$date_time_ext[same_day]=NA

	implausably_long=as.numeric(difftime(bd$date_time_acq, bd$date_time_ext))>400
	implausably_long[is.na(implausably_long)]=FALSE
	bd$date_time_acq[implausably_long]=NA
	bd$date_time_ext[implausably_long]=NA	


	ukbb_first_day_of_study=floor_date(min(bd$date_time_acq, na.rm=TRUE), unit="day")

	bd$day_of_study_acq=as.numeric(round(difftime(floor_date(bd$date_time_acq, unit="day"),ukbb_first_day_of_study, unit="day")))+1
	bd$second_of_study_acq=as.numeric(difftime(bd$date_time_acq,ukbb_first_day_of_study, unit="secs")+1)


	bd$week_of_study_acq=floor(as.numeric(bd$day_of_study_acq/7))

	bd$year_acq=year(bd$date_time_acq)

	bd$day_of_week_acq=wday(bd$date_time_acq)
	bd$day_of_year_acq=yday(bd$date_time_acq)

	bd$week_of_year_acq=floor(bd$day_of_year_acq/7)

	bd$month_of_year_acq=month(bd$date_time_acq)
	bd$month_after_jan2007_acq=as.numeric(bd$month_of_year_acq)+12*(bd$year_acq-2007)

	bd$age_acq=(bd$year_acq-bd$yob)+(bd$month_of_year_acq-bd$mob)/12
	bd$age_to_year_acq=floor(bd$age_acq)
	
	bd$second_acq=hour(bd$date_time_acq)*60*60+minute(bd$date_time_acq)*60+second(bd$date_time_acq)
	bd$quarter_acq_in_hours=floor(bd$second_acq/60/15)/4
	bd$min_acq=floor(bd$second_acq/60)

	# venipuncture time
	
	bd$day_of_study_ext=round(difftime(floor_date(bd$date_time_ext, unit="day"),ukbb_first_day_of_study, unit="day"))+1
	bd$second_of_study_ext=difftime(floor_date(bd$date_time_ext),ukbb_first_day_of_study, unit="secs")+1
	
	bd$week_of_study_ext=floor(bd$day_of_study_ext/7)

	bd$year_ext=year(bd$date_time_ext)

	bd$day_of_week_ext=wday(bd$date_time_ext)
	bd$day_of_year_ext=yday(bd$date_time_ext)

	bd$week_of_year_ext=floor(bd$day_of_year_ext/7)

	bd$month_of_year_ext=month(bd$date_time_ext)
	bd$month_after_jan2007_ext=as.numeric(bd$month_of_year_ext)+12*(bd$year_ext-2007)


	bd$second_ext=hour(bd$date_time_ext)*60*60+minute(bd$date_time_ext)*60+second(bd$date_time_ext)
	bd$quarter_ext_in_hours=floor(bd$second_ext/60/15)/4
	bd$min_ext=floor(bd$second_ext/60)


	bd$delay_mins=as.numeric(difftime(bd$date_time_acq, bd$date_time_ext, unit="mins"))
	bd$delay_in_hours_to_quarter=floor(bd$delay_mins/15)/4

	# suspect this for compatibility with INTERVAL but unused
	print("check if this for compatibility")
	bd$imputed_delay_mins=bd$delay_mins
	bd$imputed_delay_in_hours_to_quarter=bd$delay_in_hours_to_quarter
	bd$delay_was_imputed=FALSE

	# fileters and filtered variables.

	ukbb_data_filters<-data.frame(delay_good=(bd$delay_mins>0)&(bd$delay_mins<36*60))

	ukbb_data_filters$delay_good[is.na(ukbb_data_filters$delay_good)]=FALSE
	ukbb_data_filters$delay_bad=!ukbb_data_filters$delay_good
	
	bd$delay_in_hours_to_quarter_trim=bd$delay_in_hours_to_quarter
	bd$delay_in_hours_to_quarter_trim[!ukbb_data_filters$delay_good]=NA

	ukbb_data_filters$unique_machine_time_stamp=!is.element(bd$match,bd$match[duplicated(bd$match)])
	ukbb_data_filters$unique_machine_time_stamp[is.na(bd$match)]=FALSE
	
	# not many like this
	ukbb_data_filters$age_central=(bd$yob>1936)&(bd$yob<1971)
	ukbb_data_filters$age_central[is.na(ukbb_data_filters$age_central)]=FALSE

	# CANCER
	ever_had_cancer<-function(codes)
	{	
		ret_value=rep(FALSE, dim(bd)[1])
		for(code in codes)
		{
			for(col in grep("f.40006\\.[0-9]+\\.0",colnames(bd),value=TRUE))
			{
				ret_value=ret_value|grepl(code,as.character(bd[,col]))	
			}
		}
		ret_value
	}
	had_self_reported<-function(codes)
	{	
		ret_value=rep(FALSE, dim(bd)[1])
		for(code in codes)
		{
			for(col in grep("f.20002\\.0\\.[0-9]+",colnames(bd),value=TRUE))
			{
				ret_value=ret_value|grepl(code,as.character(bd[,col]))	
			}
		}
		ret_value
	}
	had_self_reported_cancer<-function(codes)
	{	
		ret_value=rep(FALSE, dim(bd)[1])
		for(code in codes)
		{
			for(col in grep("f.20001\\.0\\.[0-9]+",colnames(bd),value=TRUE))
			{
				ret_value=ret_value|grepl(code,as.character(bd[,col]))	
			}
		}
		ret_value
	}
	medication<-function(codes)
	{	
		ret_value=rep(FALSE, dim(bd)[1])
		for(code in codes)
		{
			for(col in grep("f.20003\\.0\\.[0-9]+",colnames(bd),value=TRUE))
			{
				ret_value=ret_value|grepl(code,as.character(bd[,col]))	
			}
		}
		ret_value
	}


	# cancer
	cancer_diag_date_char=apply(bd[,grepl("f.40005\\.[0-9]\\.[0-9]",colnames(bd))],2,as.character)
	cancer_diag_date=as.data.frame(matrix(NA, dim(cancer_diag_date_char)[1], dim(cancer_diag_date_char)[2]))
	cancer_diag_days_after_venepuncture=as.data.frame(matrix(NA, dim(cancer_diag_date_char)[1], dim(cancer_diag_date_char)[2]))

	for(d in 1:dim(cancer_diag_date_char)[2])
	{
		cancer_diag_date[,d]=parse_date_time(as.character(cancer_diag_date_char[,d]),orders="ymd", tz="Europe/London")
		cancer_diag_days_after_venepuncture[,d]=as.numeric(difftime(cancer_diag_date[,d], bd$date_time_ext, unit="day"))
	}

	
	cancer_diag_minimum_days_after_vene=apply(cancer_diag_days_after_venepuncture, 1, function(x){min(c(Inf,x[as.numeric(x)>0]), na.rm=TRUE)})
	cancer_diag_minimum_days_before_vene=apply(cancer_diag_days_after_venepuncture, 1, function(x){min(c(Inf,-x[as.numeric(x)<0]), na.rm=TRUE)})

	cancer_diag_previous_two_years=cancer_diag_minimum_days_before_vene<2*365
	cancer_diag_sub_year=cancer_diag_minimum_days_after_vene<365

	ukbb_data_filters$liver_neoplasm=ever_had_cancer("C22[0-9]+")|ever_had_cancer("D022")
	ukbb_data_filters$pancreatic_neoplasm=ever_had_cancer("C25[0-9]+")
	ukbb_data_filters$lung_neoplasm=ever_had_cancer(c("C34[0-9]+","C39[0-9]+"))|ever_had_cancer("D015")
	ukbb_data_filters$kidney_neoplasm=ever_had_cancer(c("C64[0-9]+","C65[0-9]+"))

	ukbb_data_filters$thymus_neoplasm=ever_had_cancer("C37[0-9]+")
	ukbb_data_filters$heart_pleura_neoplasm=ever_had_cancer("C38[0-9]+")
	ukbb_data_filters$bone_or_cartilage_neoplasm=ever_had_cancer(c("C40[0-9]+","C41[0-9]+"))
	ukbb_data_filters$thyroid_endocrine_neoplasm=ever_had_cancer(c("C73[0-9]+","C74[0-9]+","C75[0-9]+"))
	ukbb_data_filters$unspecified_secondary_neoplasm=ever_had_cancer(c("C76[0-9]+","C77[0-9]+","C78[0-9]+","C79[0-9]+","C80[0-9]+"))
	ukbb_data_filters$heam_neoplasm=ever_had_cancer(c("C81[0-9]+","C82[0-9]+","C83[0-9]+","C84[0-9]+","C85[0-9]+","C88[0-9]+","C90[0-9]+","C91[0-9]+","C92[0-9]+","C93[0-9]+","C94[0-9]+","C95[0-9]+","C96[0-9]+"))
	

	ukbb_data_filters$polycythaemia=grepl("D45[0-9]+",as.character(bd$f.40006.0.0))
	ukbb_data_filters$myelodysplastic_syndromes=grepl("D46[0-9]+",as.character(bd$f.40006.0.0))
	ukbb_data_filters$other_haem_lymph_neoplasms=grepl("D47[0-9]+",as.character(bd$f.40006.0.0))


	# self reported cancer
	ukbb_data_filters$selfr_lung_cancer=had_self_reported_cancer(1003)
	ukbb_data_filters$selfr_liver_cancer=had_self_reported_cancer(1024)
	ukbb_data_filters$selfr_pancreatic_cancer=had_self_reported_cancer(1026)
	ukbb_data_filters$selfr_small_cell_cancer=had_self_reported_cancer(1027)
	ukbb_data_filters$selfr_large_cell_cancer=had_self_reported_cancer(1028)
	ukbb_data_filters$selfr_kidney_cancer=had_self_reported_cancer(1034)
	ukbb_data_filters$selfr_lymphoma=had_self_reported_cancer(1047)|had_self_reported_cancer(1052)|had_self_reported_cancer(1053)
	ukbb_data_filters$selfr_leukaemia=had_self_reported_cancer(1048)
	ukbb_data_filters$selfr_mutiple_myeloma=had_self_reported_cancer(1050)
	ukbb_data_filters$selfr_mutiple_myelofibrosis_or_myelodysplasia=had_self_reported_cancer(1051)
	ukbb_data_filters$selfr_chronic_lympocytic=had_self_reported_cancer(1055)
	ukbb_data_filters$selfr_chronic_myeloid=had_self_reported_cancer(1056)
	ukbb_data_filters$selfr_other_haem_malig=had_self_reported_cancer(1058)
	ukbb_data_filters$selfr_acute_myeloid=had_self_reported_cancer(1074)
	ukbb_data_filters$selfr_adrenal_cancer=had_self_reported_cancer(1066)
	ukbb_data_filters$selfr_parathyroid_cancer=had_self_reported_cancer(1067)
	ukbb_data_filters$selfr_primary_bone_cancer=had_self_reported_cancer(1063)
	ukbb_data_filters$selfr_thyroid_cancer=had_self_reported_cancer(1065)
	ukbb_data_filters$selfr_malignant_lymph=had_self_reported_cancer(1070)
	ukbb_data_filters$selfr_metastatic_cancer=had_self_reported_cancer(1071)
	ukbb_data_filters$selfr_respiratory_cancer=had_self_reported_cancer(1084)
	ukbb_data_filters$selfr_bone_secondaries=had_self_reported_cancer(1085)

	ukbb_data_filters$haem_cancer_histology=is.element(bd$f.40011.0.0, as.character(9590:9999))

	ukbb_data_filters$selfr_asthma=had_self_reported(1111)

	ukbb_data_filters$selfr_hiv_aids=had_self_reported(1429)
	ukbb_data_filters$selfr_hyper_tension=had_self_reported(1065)
	ukbb_data_filters$selfr_alcohol_dependency=had_self_reported(1408)
	ukbb_data_filters$selfr_opiod_dependency=had_self_reported(1409)
	ukbb_data_filters$selfr_other_dependency=had_self_reported(1410)
	ukbb_data_filters$selfr_haemophilia=had_self_reported(1328)
	ukbb_data_filters$selfr_sickle_cell=had_self_reported(1339)
	ukbb_data_filters$selfr_thalassaemia=had_self_reported(1340)
	ukbb_data_filters$selfr_polycythaemia_vera=had_self_reported(1438)
	ukbb_data_filters$selfr_alcoholic_liver_disease=had_self_reported(1604)
	ukbb_data_filters$selfr_hepatitis_a=had_self_reported(1578)
	ukbb_data_filters$selfr_hepatitis_b=had_self_reported(1579)
	ukbb_data_filters$selfr_hepatitis_c=had_self_reported(1580)
	ukbb_data_filters$selfr_hepatitis_d=had_self_reported(1581)
	ukbb_data_filters$selfr_hepatitis_e=had_self_reported(1582)
	ukbb_data_filters$selfr_hepatitis=had_self_reported(1155)
	ukbb_data_filters$selfr_viral_hepatitis=had_self_reported(1156)
	ukbb_data_filters$selfr_non_infective_hepatitis=had_self_reported(1157)
	ukbb_data_filters$selfr_liver_failure=had_self_reported(1158)
	ukbb_data_filters$selfr_respiratory_failure=had_self_reported(1124)
	ukbb_data_filters$selfr_liver_biliary_pancreas_problem=had_self_reported(1136)
	ukbb_data_filters$selfr_hypothyroidism=had_self_reported(1126)
	ukbb_data_filters$selfr_coeliac=had_self_reported(1456)
	ukbb_data_filters$selfr_crohns=had_self_reported(1462)
	ukbb_data_filters$selfr_diabetes=had_self_reported(1220)
	ukbb_data_filters$selfr_typeI_diabetes=had_self_reported(1222)
	ukbb_data_filters$selfr_typeII_diabetes=had_self_reported(1223)
	ukbb_data_filters$dr_diag_diabetes=bd$f.2443.0.0==1
	ukbb_data_filters$dr_diag_diabetes[is.na(ukbb_data_filters$dr_diag_diabetes)]=FALSE
	ukbb_data_filters$selfr_ibd=had_self_reported(1461)

	ukbb_data_filters$selfr_gastritis=had_self_reported(1143)
	ukbb_data_filters$selfr_haemochromatosis=had_self_reported(1507)
	ukbb_data_filters$selfr_heart_failure=had_self_reported(1076)

	ukbb_data_filters$selfr_rheumatoid_arthritis=had_self_reported(1464)

	ukbb_data_filters$erythropoietin_med=medication(1140870568)
	ukbb_data_filters$metformin_med=medication(c(1140884600,1141189090))
	ukbb_data_filters$insulin_med=(bd$f.6177.0.0==3)
	ukbb_data_filters$insulin_med[is.na(ukbb_data_filters$insulin_med)]=FALSE
	ukbb_data_filters$insulin_med=ukbb_data_filters$insulin_med|medication(1140883066)
	
	ukbb_data_filters$selfr_pancytopenia=had_self_reported(1447)
	ukbb_data_filters$selfr_neutropenia_or_lymphopenia=had_self_reported(1448)
	ukbb_data_filters$selfr_myeloprolif=had_self_reported(1449)
	ukbb_data_filters$selfr_hered_haem_disord=had_self_reported(1451)
	ukbb_data_filters$selfr_lymphoedema=had_self_reported(1495)
	ukbb_data_filters$selfr_essential_thrombocytosis=had_self_reported(1546)
	ukbb_data_filters$selfr_monoclonal_gammopathy=had_self_reported(1450)

	ukbb_data_filters$selfr_copd=had_self_reported(1112)
	ukbb_data_filters$selfr_emphysema=had_self_reported(1113)

	ukbb_data_filters$selfr_myelofibrosis=had_self_reported(1658)
	ukbb_data_filters$selfr_sleep_apnoea=had_self_reported(1123)

	ukbb_data_filters$selfr_kidney_nephropathy=had_self_reported(1519)
	ukbb_data_filters$selfr_renal_kidney_failure=had_self_reported(1192)
	ukbb_data_filters$selfr_renal_failure_dialysis=had_self_reported(1193)
	ukbb_data_filters$selfr_renal_failure_no_dialysis=had_self_reported(1194)

	
	ukbb_data_filters$pregnant=bd$f.3140.0.0==1
	ukbb_data_filters$pregnant[is.na(ukbb_data_filters$pregnant)]=FALSE
	
	ukbb_data_filters$haem_cancer_histology=is.element(bd$f.40011.0.0, as.character(9590:9999))
	ukbb_data_filters$heptocellular_carcinoma_histology=is.element(bd$f.40011.0.0, as.character(8170:8175))
	ukbb_data_filters$renalcell_carcinoma_histology=is.element(bd$f.40011.0.0, as.character(c(8312,8317,8318)))
	ukbb_data_filters$pheochromocytoma_histology=is.element(bd$f.40011.0.0, as.character(8700))
	ukbb_data_filters$hemangioblastoma_histology=is.element(bd$f.40011.0.0, as.character(9161))


	#genetics
	ukbb_data_filters$genetic_outlier_filter=!bd$CEU_QCed_genetic_eur
	ukbb_data_filters$genetic_outlier_filters[is.na(ukbb_data_filters$genetic_outlier_filter)]=FALSE
	ukbb_data_filters$MZ_twin_filter=bd$is_MZ_twin
	ukbb_data_filters$MZ_twin_filter[is.na(ukbb_data_filters$MZ_twin_filter)]=FALSE

	ukbb_data_filters=as.data.frame(ukbb_data_filters)

	ukbb_filter_lists=list()
	ukbb_filter_lists$blood_cancer_filters=c("haem_cancer_histology","selfr_myelofibrosis","heam_neoplasm","myelodysplastic_syndromes","other_haem_lymph_neoplasms","selfr_lymphoma","selfr_leukaemia","selfr_malignant_lymph","selfr_mutiple_myeloma","selfr_mutiple_myelofibrosis_or_myelodysplasia","selfr_chronic_lympocytic","selfr_chronic_myeloid","selfr_other_haem_malig","selfr_acute_myeloid","selfr_polycythaemia_vera","polycythaemia","selfr_myeloprolif","selfr_essential_thrombocytosis")

	ukbb_filter_lists$non_blood_cancer_filters=c("liver_neoplasm","pancreatic_neoplasm","lung_neoplasm","kidney_neoplasm","thymus_neoplasm","heart_pleura_neoplasm","bone_or_cartilage_neoplasm","thyroid_endocrine_neoplasm","selfr_pancreatic_cancer","selfr_small_cell_cancer","selfr_thyroid_cancer","selfr_metastatic_cancer","selfr_respiratory_cancer","selfr_bone_secondaries","unspecified_secondary_neoplasm","selfr_lung_cancer","selfr_liver_cancer","selfr_large_cell_cancer","selfr_kidney_cancer","selfr_adrenal_cancer","selfr_parathyroid_cancer","selfr_primary_bone_cancer","heptocellular_carcinoma_histology","renalcell_carcinoma_histology","hemangioblastoma","pheochromocytoma_histology","hemangioblastoma_histology")
		
	ukbb_filter_lists$non_cancer_blood_condition_filters=c("selfr_monoclonal_gammopathy","selfr_hered_haem_disord","selfr_haemochromatosis","selfr_thalassaemia","selfr_haemophilia","selfr_sickle_cell", "selfr_neutropenia_or_lymphopenia","selfr_pancytopenia")
 	ukbb_filter_lists$non_cancer_liver_condition_filters=c("selfr_alcoholic_liver_disease","selfr_hepatitis_a","selfr_hepatitis_b","selfr_hepatitis_c","selfr_hepatitis_d","selfr_hepatitis_e","selfr_hepatitis","selfr_viral_hepatitis","selfr_non_infective_hepatitis","selfr_liver_failure","selfr_liver_biliary_pancreas_problem")
	ukbb_filter_lists$non_cancer_renal_condition_filters=c("selfr_kidney_nephropathy","selfr_renal_failure_no_dialysis","selfr_renal_kidney_failure","selfr_renal_failure_dialysis")
	ukbb_filter_lists$non_cancer_respiratory_condition_filters=c("selfr_copd","selfr_emphysema","selfr_sleep_apnoea","selfr_respiratory_failure")

	ukbb_filter_lists$autoimmune_disorder_filters=c("selfr_typeI_diabetes","selfr_typeII_diabetes","selfr_diabetes","dr_diag_diabetes","selfr_coeliac","selfr_crohns","selfr_rheumatoid_arthritis","selfr_ibd")	
	ukbb_filter_lists$residual_systemic_condition_filters=c("pregnant", "selfr_gastritis","selfr_hypothyroidism","selfr_lymphoedema","selfr_hiv_aids", "selfr_heart_failure", "selfr_asthma", "selfr_hyper_tension")
	ukbb_filter_lists$blood_altering_med_filters=c("erythropoietin_med","metformin_med","insulin_med")
	
	ukbb_filter_lists$drug_dependancy_filters=c("selfr_alcohol_dependency","selfr_opiod_dependency","selfr_other_dependency")

	ukbb_filter_lists$genotypic_exclusion_filters=c("MZ_twin_filter","genetic_outlier_filter")


	bd_extract=bd[,(colnames(bd)!="f.eid")&(!grepl("f.[0-9].+", colnames(bd)))]     	
	
	new_variables=compute_derived_blood_count_variables(bd_extract, "", instrument="Coulter")
	bd_extract=cbind(bd_extract[,!is.element(colnames(bd_extract), colnames(new_variables))], new_variables)

	save(list=c("bd_extract"), file=sprintf("%s/Rdata/ret_ukbb_fbc_prepared.Rdata", UKBIOBANK_DATA_DIR))
	save(list=c("ukbb_first_day_of_study","ukbb_data_filters","ukbb_filter_lists"), file=sprintf("%s/Rdata/ret_ukbb_fbc_ancillary.Rdata", UKBIOBANK_DATA_DIR))

	
}


# adjust specific trait and save to disk (e.g. for HPC operation)
adjust_UKBB_variables_for_technical<-function(trait)
{
	ukbb_tech_r.sq=NULL
	ukbb_tech_r.sq_full_data=NULL

	if(hasArg(trait))
	{
		traits=trait
	}else
	{
		traits=tf_adjustable(tf_coulter_measured())
	}
	registerDoMC(cores=length(traits))
	
	adjust_for_tech_wrap_function<-function(bd_extract, ukbb_data_filters, traits, trait_num)
	{
		adjust_variable_for_technical(bd_extract, ukbb_data_filters, traits[trait_num])
	}
		
	result_list=foreach(n = 1:length(traits)) %dopar% adjust_for_tech_wrap_function(bd_extract, ukbb_data_filters, traits, n)

	for(n in 1:length(traits))
	{
		eval(parse(text=sprintf("bd_extract$%s_tech_adj=result_list[[n]]$adjusted_trait", traits[n])))
		eval(parse(text=sprintf("ukbb_tech_r.sq=c(ukbb_tech_r.sq,\"%s\"=result_list[[n]]$r.sq)", traits[n])))
		eval(parse(text=sprintf("ukbb_tech_r.sq_full_data=c(ukbb_tech_r.sq_full_data,\"%s\"=result_list[[n]]$r.sq_full_data)", traits[n])))
	}
	
	if(hasArg(trait))
	{
		dir.create(sprintf("%s/trait_prep",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)	
		eval(parse(text=sprintf("write.table(bd_extract$%s_tech_adj, file=\"%s/trait_prep/%s_tech_adj.tsv\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote=FALSE)", trait, UKBIOBANK_FBC_OUTPUT_DIR, trait)))
		eval(parse(text=sprintf("write.table(ukbb_tech_r.sq, file=\"%s/trait_prep/%s_tech_adj_r.sq.tsv\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote=FALSE)", UKBIOBANK_FBC_OUTPUT_DIR, trait)))
		eval(parse(text=sprintf("write.table(ukbb_tech_r.sq_full_data, file=\"%s/trait_prep/%s_tech_adj_r.sq_full_data.tsv\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote=FALSE)",  UKBIOBANK_FBC_OUTPUT_DIR, trait)))
	
	}
	else
	{	
		new_variables=compute_derived_blood_count_variables(bd_extract, "_tech_adj", "Coulter")
		bd_extract<-cbind(bd_extract[,!is.element(colnames(bd_extract), colnames(new_variables))], new_variables)
		save(list=c("bd_extract","ukbb_tech_r.sq","ukbb_tech_r.sq_full_data"), file=sprintf("%s/Rdata/ret_ukbb_fbc_tech_adjusted.Rdata", UKBIOBANK_DATA_DIR))
	}
	
}

adjust_UKBB_variables_for_gwas<-function(trait,draw_plots=FALSE)
{
	ukbb_gwas_r.sq=NULL
	ukbb_gwas_r.sq_full_data=NULL

	if(hasArg(trait))
	{
		traits=trait
	}
	else
	{
		traits=tf_adjustable(tf_coulter())
	}
	registerDoMC(cores=length(traits))

	adjust_for_gwas_wrap_function<-function(bd_extract, ukbb_data_filters, traits, trait_num)
	{
		adjust_variable_for_gwas(bd_extract, ukbb_data_filters, sprintf("%s_tech_adj",traits[trait_num]), study="UK Biobank")
	}
	
	result_list=foreach(n = 1:length(traits)) %dopar% adjust_for_gwas_wrap_function(bd_extract, ukbb_data_filters, traits, n)
	
	for(n in 1:length(traits))
	{
		eval(parse(text=sprintf("bd_extract$%s_gwas_adj=result_list[[n]]$adjusted_trait", traits[n])))
		if(length(tf_0100(get_base_trait(traits[n])))>0)
		{
			eval(parse(text=sprintf("bd_extract$%s_gwas_adj_supR=logit(result_list[[n]]$adjusted_trait/100)", traits[n])))
	
		}
		if(length(tf_01(get_base_trait(traits[n])))>0)
		{
			eval(parse(text=sprintf("bd_extract$%s_gwas_adj_supR=logit(result_list[[n]]$adjusted_trait)", traits[n])))
	
		}
		if(length(tf_positive(get_base_trait(traits[n])))>0)
		{
			eval(parse(text=sprintf("bd_extract$%s_gwas_adj_supR=log(result_list[[n]]$adjusted_trait)", traits[n])))
	
		}

		eval(parse(text=sprintf("ukbb_gwas_r.sq=c(ukbb_gwas_r.sq,\"%s\"=result_list[[n]]$r.sq)", traits[n])))
		eval(parse(text=sprintf("ukbb_gwas_r.sq_full_data=c(ukbb_gwas_r.sq_full_data,\"%s\"=result_list[[n]]$r.sq_full_data)", traits[n])))
	}
	
	if(hasArg(trait))
	{
		dir.create(sprintf("%s/trait_prep",UKBIOBANK_FBC_OUTPUT_DIR), showWarnings = FALSE)	
		eval(parse(text=sprintf("write.table(bd_extract$%s_gwas_adj, file=\"%s/trait_prep/%s_gwas_adj.tsv\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote=FALSE)", trait, UKBIOBANK_FBC_OUTPUT_DIR, trait)))
		eval(parse(text=sprintf("write.table(ukbb_gwas_r.sq, file=\"%s/trait_prep/%s_gwas_adj_r.sq.tsv\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote=FALSE)", UKBIOBANK_FBC_OUTPUT_DIR, trait)))
		eval(parse(text=sprintf("write.table(ukbb_gwas_r.sq_full_data, file=\"%s/trait_prep/%s_gwas_adj_r.sq_full_data.tsv\", row.names=FALSE, col.names=FALSE, sep=\"\t\", quote=FALSE)",  UKBIOBANK_FBC_OUTPUT_DIR, trait)))

	}
	else
	{
		save(list=c("bd_extract", "ukbb_gwas_r.sq","ukbb_gwas_r.sq_full_data"), file=sprintf("%s/Rdata/ret_ukbb_fbc_gwas_adjusted.Rdata", UKBIOBANK_DATA_DIR))
	}
	
}


prepare_UKBB_traits_for_gwas_and_save<-function()
{


	bd_local=bd_extract
	genetic_exclusions=apply(ukbb_data_filters[,ukbb_filter_lists$genotypic_exclusion_filters],1,any)
	meno_no_na=bd_local$meno
	meno_no_na[is.na(bd_local$meno)]="missing"

	ukbileve_status_no_na=grepl("UKBiLEVE",bd_extract$genotyping_batch)
	ukbileve_status_no_na[is.na(ukbileve_status_no_na)]=FALSE
	instrument_no_na=bd_local$instrument
	instrument_no_na[is.na(bd_local$instrument)]="missing"

	ukbb_meno_strata=interaction(ukbileve_status_no_na,meno_no_na)

	platelet_exclusions=bd_extract$mpv_tech_adj>quantile(bd_extract$mpv_tech_adj,na.rm=TRUE, 1-0.04377202)
	platelet_exclusions[is.na(platelet_exclusions)]=FALSE

	# identify univariate outliers/large adjustments
	for(trait in tf_500kgwas())
	{
		eval(parse(text=sprintf("original_trait_data=bd_local$%s", trait)))
		eval(parse(text=sprintf("adjust_scale_trait_data=bd_local$%s_gwas_adj_supR", trait)))
	

		if(length(tf_0100(get_base_trait(trait)))>0)
				
		{	adjust_scale_original_trait_data=logit(original_trait_data/100)

		}else
		{
			if(length(tf_01(get_base_trait(trait)))>0)
			{
				adjust_scale_original_trait_data=logit(original_trait_data)

	
			}
			else
			{
				if(length(tf_positive(get_base_trait(trait)))>0)
				{	
					adjust_scale_original_trait_data=log(original_trait_data)

				}
			}
		}
	
	
		trait_diff=adjust_scale_original_trait_data-adjust_scale_trait_data
		big_difference=abs(trait_diff)>4*mad(trait_diff, na.rm=TRUE)
		big_difference[is.na(big_difference)]=FALSE	

		outlier_filter=(adjust_scale_trait_data>mad(adjust_scale_trait_data,na.rm=TRUE)*4.5+median(adjust_scale_trait_data,na.rm=TRUE))|(adjust_scale_trait_data<(-mad(adjust_scale_trait_data, na.rm=TRUE)*4.5+median(adjust_scale_trait_data,na.rm=TRUE)))
		outlier_filter[is.na(outlier_filter)]=FALSE
		print("total lost")
		print(sum(big_difference|outlier_filter))	
		eval(parse(text=sprintf("bd_local$%s_gwas_NAs<-big_difference|outlier_filter", trait)))
	
	}

	rownames(bd_local)=1:dim(bd_local)[1]
	
	
	
	# platelet traits: essentially 3 independent  measurements
	princomp_plt_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_platelet(tf_500kgwas())))])
	pca_plt=princomp(princomp_plt_data, cor=TRUE)
	plt_score_sq=apply(scale(pca_plt$scores[,1:3])^2, 1,sum)
	plt_multivariate_outliers=is.element(rownames(bd_local),names(plt_score_sq)[plt_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=3)])

	# all red traits: essentially 8 independent measurements 
	princomp_red_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_red(tf_500kgwas())))])
	pca_red=princomp(princomp_red_data, cor=TRUE)
	red_score_sq=apply(scale(pca_red$scores[,1:8])^2, 1,sum)
	red_cell_multivariate_outliers=is.element(rownames(bd_local),names(red_score_sq)[red_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=8)])

	# immature red traits: essentially 3 independent  measurements 
	princomp_immature_red_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_immature_red(tf_500kgwas())))])
	pca_immature_red=princomp(princomp_immature_red_data, cor=TRUE)
	immature_red_score_sq=apply(scale(pca_immature_red$scores[,1:3])^2, 1,sum)
	immature_red_multivariate_outliers=is.element(rownames(bd_local),names(immature_red_score_sq)[immature_red_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=3)])
	
	# mature red traits: essentially 5 independent  measurements 
	princomp_mature_red_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_mature_red(tf_500kgwas())))])
	pca_mature_red=princomp(princomp_mature_red_data, cor=TRUE)
	mature_red_score_sq=apply(scale(pca_mature_red$scores[,1:5])^2, 1,sum)
	mature_red_multivariate_outliers=is.element(rownames(bd_local),names(mature_red_score_sq)[mature_red_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=5)])


	# all white traits: essentially 5 independent  measurements 
	princomp_white_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_white(tf_500kgwas())))])
	pca_white=princomp(princomp_white_data, cor=TRUE)
	white_score_sq=apply(scale(pca_white$scores[,1:5])^2, 1,sum)
	white_cell_multivariate_outliers=is.element(rownames(bd_local),names(white_score_sq)[white_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=5)])

	# all traits: essentially: 15 independent measurements 
	princomp_all_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_500kgwas()))])
	pca_all=princomp(princomp_all_data, cor=TRUE)
	all_score_sq=apply(scale(pca_all$scores[,1:15])^2, 1,sum)
	all_multivariate_outliers=is.element(rownames(bd_local),names(all_score_sq)[all_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=15)])


	# granulocyte traits: 3 independent measurements
	princomp_gran_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_granulocyte(tf_500kgwas())))])
	pca_gran=princomp(princomp_gran_data, cor=TRUE)
	gran_score_sq=apply(scale(pca_gran$scores[,1:3])^2, 1,sum)
	gran_multivariate_outliers=is.element(rownames(bd_local),names(gran_score_sq)[gran_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=3)])


	# myeloid traits: 4 independent measurements
	princomp_myeloid_data=na.omit(bd_local[,is.element(colnames(bd_local),sprintf("%s_gwas_adj_supR",tf_myeloid_white(tf_500kgwas())))])
	pca_myeloid=princomp(princomp_myeloid_data, cor=TRUE)
	myeloid_score_sq=apply(scale(pca_myeloid$scores[,1:4])^2, 1,sum)
	myeloid_multivariate_outliers=is.element(rownames(bd_local),names(myeloid_score_sq)[myeloid_score_sq>qchisq(-9*log(10),log.p=TRUE, lower.tail=FALSE, df=4)])


	union_multivariate_outliers=plt_multivariate_outliers|red_cell_multivariate_outliers|mature_red_multivariate_outliers| immature_red_multivariate_outliers|white_cell_multivariate_outliers|all_multivariate_outliers|gran_multivariate_outliers|myeloid_multivariate_outliers
	print("all multi lost")
	print(sum(union_multivariate_outliers))

	# do the outlier removal
	for(trait in tf_500kgwas())
	{
		print(trait)
		
		eval(parse(text=sprintf("prep_trait=bd_local$%s_gwas_adj", trait)))

		prep_trait[all_multivariate_outliers]<-NA
		prep_trait[genetic_exclusions]<-NA

		if(is.element(trait,tf_platelet()))
		{
			prep_trait[plt_multivariate_outliers|platelet_exclusions]<-NA
			large_plts=bd_extract$mpv_tech_adj>13
			large_plts[is.na(large_plts)]=FALSE
			prep_trait[large_plts]<-NA

		}
		if(is.element(trait,tf_red()))
		{
			prep_trait[red_cell_multivariate_outliers]<-NA
		}
		if(is.element(trait,tf_immature_red()))
		{
			prep_trait[immature_red_multivariate_outliers]<-NA
		}
		if(is.element(trait,tf_mature_red()))
		{
			prep_trait[mature_red_multivariate_outliers]<-NA
		}

		if(is.element(trait,tf_white()))
		{
			prep_trait[white_cell_multivariate_outliers|gran_multivariate_outliers|myeloid_multivariate_outliers]<-NA
		}

		eval(parse(text=sprintf("prep_trait[bd_local$%s_gwas_NAs]<-NA", trait)))
		# adjust for clinic / calling batch and standardise
		
		#generate clinic variables
		quantile_normalised=quantile_normalise(quantile_normalise(quantile_normalise(prep_trait,strata=instrument_no_na), strata=ukbb_meno_strata),strata=ukbileve_status_no_na)		
		bd_extract$clinic[is.na(bd_extract$clinic)]="Unknown"
		eval(parse(text=sprintf("bd_extract$%s_gwas_normalised=scale(resid(lm(quantile_normalised~as.factor(bd_extract$clinic)+as.factor(bd_extract$genotyping_batch), na.action=na.exclude)))",trait)))		
	
	}
	
	save(list=c("bd_extract"), file=sprintf("%s/Rdata/ret_ukbb_fbc_gwas_prepared.Rdata", UKBIOBANK_DATA_DIR))
	
}


generate_UKBB_150k_sample_files<-function(normalised=TRUE)
{
	model_sample_file=read.csv(sprintf("%s/7439_sample/impv1.sample",UKBIOBANK_FBC_DATA_DIR), skip=2, head=FALSE, stringsAsFactors=FALSE, sep=" ")

	sample_table=data.frame(FID=as.character(model_sample_file[,1]), stringsAsFactors=FALSE)
	sample_table$IID=sample_table$FID
	sample_table$missing=as.character(0)


	#read PCs
	PCs=read.table(sprintf("%s/gwas/150k_PCs/ukbb_inclusion_list_50PCs_eigenvectors.txt",UKBIOBANK_FBC_DATA_DIR), stringsAsFactors=FALSE)
	inclusion_list=read.table(sprintf("%s/gwas/150k_PCs/inclusion_list_plink.txt",UKBIOBANK_FBC_DATA_DIR), stringsAsFactors=FALSE)
	colnames(PCs)=sprintf("PC%d", 1:10)
	sample_table=cbind(sample_table,PCs[match(sample_table$FID, inclusion_list[,1]),1:10])


	#generate clinic variables
	bd_extract$CLINIC[grepl("Bart",bd_extract$CLINIC)]="Barts_London"
	bd_extract$CLINIC[is.na(bd_extract$CLINIC)]="Unknown"

	clinic=bd_extract$CLINIC[match(sample_table$FID,bd_extract$f.eid_7439)]	
	sample_table=cbind(sample_table, clinic)

	#add phenotypes
	if(normalised)
	{
		sample_table=cbind(sample_table, bd_extract[match(sample_table$FID,bd_extract$f.eid_7439),grepl("gwas_normalised",colnames(bd_extract))])
	}
	else
	{			
		sample_table=cbind(sample_table, bd_extract[match(sample_table$FID,bd_extract$f.eid_7439),grepl("true_scale_gwas_adj",colnames(bd_extract))])

	}	
	# generate genotyping batch variables

#	batch_vars=NULL
#	batches=unique(bd_extract$genotyping_batch)
#	batches=batch[batch!=2000]

#	for(batch in batches)
#	{
#		eval(parse(text=sprintf("batch_vars$batch_%s=as.numeric(bd_extract$genotyping_batch==batch)",batch)))
#	}
#	batch_vars=data.frame(batch_vars, stringsAsFactors=FALSE)
#	sample_table=cbind(sample_table, batch_vars[match(sample_table$FID,bd_extract$f.eid_7439),])

	if(normalised)
	{
		header_line=c("0","0","0",rep("C", 10), "D",rep("P", sum(grepl("gwas_normalised",colnames(bd_extract)))))
	}
	else
	{
		header_line=c("0","0","0",rep("C", 10), "D",rep("P", sum(grepl("true_scale_gwas_adj",colnames(bd_extract)))))
	}
	ukbileve_status=bd_extract$ukbileve_status[match(sample_table$FID,bd_extract$f.eid_7439)]
	ukbileve_status[is.na(ukbileve_status)]=FALSE	

	sample_table_ukbileve=sample_table
	sample_table_non_ukbileve=sample_table
	sample_table_ukbileve[!ukbileve_status,4:dim(sample_table_ukbileve)[2]]=NA
	sample_table_non_ukbileve[ukbileve_status,4:dim(sample_table_non_ukbileve)[2]]=NA

	if(normalised)
	{
		post_fix=""
	}
	else
	{
		post_fix="_unquantised"
	}

	write.table(t(colnames(sample_table)),sprintf("%s/gwas/150k_sample_files/ukbb_150k_not_ukbileve%s.sample",UKBIOBANK_FBC_DATA_DIR, post_fix), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table(t(colnames(sample_table)),sprintf("%s/gwas/150k_sample_files/ukbb_150k_ukbileve%s.sample",UKBIOBANK_FBC_DATA_DIR,post_fix), sep=" ", quote=FALSE, row.names=FALSE,col.names=FALSE)

	write.table(t(header_line),sprintf("%s/gwas/150k_sample_files/ukbb_150k_not_ukbileve%s.sample",UKBIOBANK_FBC_DATA_DIR,post_fix), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	write.table(t(header_line),sprintf("%s/gwas/150k_sample_files/ukbb_150k_ukbileve%s.sample",UKBIOBANK_FBC_DATA_DIR,post_fix), sep=" ", quote=FALSE, row.names=FALSE,col.names=FALSE, append=TRUE)

	write.table(sample_table_non_ukbileve,sprintf("%s/gwas/150k_sample_files/ukbb_150k_not_ukbileve%s.sample",UKBIOBANK_FBC_DATA_DIR,post_fix), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	write.table(sample_table_ukbileve,sprintf("%s/gwas/150k_sample_files/ukbb_150k_ukbileve%s.sample",UKBIOBANK_FBC_DATA_DIR,post_fix), sep=" ", quote=FALSE, row.names=FALSE,col.names=FALSE, append=TRUE)

}


generate_UKBB_500k_sample_files<-function(normalised=TRUE, bcx=FALSE)
{
	sample_table=data.frame(FID=1:max(bd_extract$sample_file_index, na.rm=TRUE), stringsAsFactors=FALSE)
	sample_table$IID=sample_table$FID
	sample_table$missing=as.character(0)

	if(bcx==FALSE)
	{
		traits=tf_500kgwas()
	}
	else
	{
		traits=tf_bcx()
	}
	#read PCs
	sample_table=cbind(sample_table,bd_extract[match(sample_table$IID, bd_extract$sample_file_index),sprintf("PC%d",1:10)])


	#add phenotypes
	if(normalised)
	{
		sample_table=cbind(sample_table, bd_extract[match(sample_table$IID, bd_extract$sample_file_index),sprintf("%s_gwas_normalised",traits)])
	}
	else
	{
		
		sample_table=cbind(sample_table, bd_extract[match(sample_table$IID, bd_extract$sample_file_index),sprintf("%s_gwas_adj_supR",traits)])


	}	


	header_line=c("0","0","0",rep("C", 10),rep("P", length(traits)))
	
	if(normalised)
	{
		post_fix=""
	}
	else
	{
		post_fix="_unquantised"
	}
	if(bcx)
	{
		filename=sprintf("ukbb_500k_bcx%s.sample", post_fix)
	}	
	else
	{	
		filename=sprintf("ukbb_500k%s.sample", post_fix)
	
	}
	write.table(t(colnames(sample_table)),sprintf("%s/gwas/500k_sample_files/%s",UKBIOBANK_DATA_DIR, filename), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
	write.table(t(header_line),sprintf("%s/gwas/500k_sample_files/%s",UKBIOBANK_DATA_DIR,filename), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
	write.table(sample_table,sprintf("%s/gwas/500k_sample_files/%s",UKBIOBANK_DATA_DIR,filename), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

}

write_data_for_bcx<-function()
{
	bcx_data=bd_extract[order(bd_extract$sample_file_index),c("sample_file_index", sprintf("%s_tech_adj",tf_bcx()), sprintf("%s_gwas_adj",tf_bcx()))]
	bcx_data_no_na=bcx_data[!is.na(bcx_data$sample_file_index),]
	write.table(bcx_data_no_na,sprintf("%s/gwas/500k_sample_files/ukbb_500k_adjusted.tsv",UKBIOBANK_DATA_DIR), sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)



}

write_data_for_cardio<-function()
{
	cardio_data=bd_extract[order(bd_extract$sample_file_index),c("sample_file_index", sprintf("%s_tech_adj",tf_500kgwas()), sprintf("%s_gwas_adj",tf_500kgwas()), sprintf("%s_gwas_normalised",tf_500kgwas()))]
	cardio_data_no_na=cardio_data[!is.na(cardio_data$sample_file_index),]
	write.table(cardio_data_no_na,sprintf("%s/gwas/500k_sample_files/ukbb_500k_gwas_traits_adjusted.tsv",UKBIOBANK_DATA_DIR), sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)



}


