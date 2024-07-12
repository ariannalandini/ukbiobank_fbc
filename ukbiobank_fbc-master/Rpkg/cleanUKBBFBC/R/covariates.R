get_covariate_xlab<-function(covariate)
{
	if(covariate=="day_of_study_acq"){ return(xlab("Acquisition time (Days into study)")) }
	if(covariate=="week_of_study_acq"){ return(xlab("Acquisition time (Weeks into study)")) }
	if(covariate=="week_of_year_acq"){ return(xlab("Acquisition time (Weeks into year)")) }
	if(covariate=="month_of_year_acq"){ return(xlab("Acquisition time (Month of year)")) }
	if(covariate=="quarter_acq_in_hours"){ return(xlab("Acquisition time (Hours after midnight)")) }
	if(covariate=="day_of_week_acq"){ return(xlab("Aquisition time (Day of week; Sun=1)")) }

	if(covariate=="month_after_jan2007_acq"){ return(xlab("Acquisition time (Months after January 2007)")) }
	if(covariate=="month_after_jan2012_acq"){ return(xlab("Acquisition time (Months after January 2012)")) }

	if(covariate=="day_of_study_ext"){ return(xlab("Venepuncture time (Days into study)")) }
	if(covariate=="week_of_study_ext"){ return(xlab("Venepuncture time (Weeks into study)")) }
	if(covariate=="week_of_year_ext"){ return(xlab("Venepuncture time (Weeks into year)")) }
	if(covariate=="month_of_year_ext"){ return(xlab("Venepuncture time (Month of year)")) }
	if(covariate=="quarter_ext_in_hours"){ return(xlab("Venepunture time after midnight (Hours)")) }
	if(covariate=="month_after_jan2007_ext"){ return(xlab("Venepuncture time after January 2007 (Months)")) }
	if(covariate=="month_after_jan2012_ext"){ return(xlab("Venepuncture time after January 2012 (Months)")) }

	if(covariate=="day_of_week_ext"){ return(xlab("Venepuncture time (Day of week; Sun=1)")) }

	if((covariate=="delay_in_hours_to_quarter")||(covariate=="delay_in_hours_to_quarter_trim")){ return(xlab("Delay between venepuncture and acquisition (Hours)")) }
	if((covariate=="age_to_year_acq")){ return(xlab("Age (Years)")) }
	if((covariate=="instrument")){ return(xlab("Instrument")) }


	return(xlab("NO LABEL"))
}

get_covariate_name<-function(covariate)
{
	if(covariate=="day_of_study_acq"){ return("Acquisition time") }
	if(covariate=="week_of_study_acq"){ return("Acquisition time") }
	if(covariate=="week_of_year_acq"){ return("Acquisition time") }
	if(covariate=="month_of_year_acq"){ return("Acquisition time") }
	if(covariate=="quarter_acq_in_hours"){ return("Acquisition time") }
	if(covariate=="day_of_week_acq"){ return("Aquisition time") }

	if(covariate=="month_after_jan2007_acq"){ return("Acquisition time") }
	if(covariate=="month_after_jan2012_acq"){ return("Acquisition time") }

	if(covariate=="day_of_study_ext"){ return("Venepuncture time") }
	if(covariate=="week_of_study_ext"){ return("Venepuncture time") }
	if(covariate=="week_of_year_ext"){ return("Venepuncture time") }
	if(covariate=="month_of_year_ext"){ return("Venepuncture time") }
	if(covariate=="quarter_ext_in_hours"){ return("Venepunture time") }
	if(covariate=="month_after_jan2007_ext"){ return("Venepuncture time") }
	if(covariate=="month_after_jan2012_ext"){ return("Venepuncture time") }

	if(covariate=="day_of_week_ext"){ return("Venepuncture time") }

	if((covariate=="delay_in_hours_to_quarter")||(covariate=="delay_in_hours_to_quarter_trim")){ return("Delay between venepuncture and acquisition") }
	if((covariate=="age_to_year_acq")){ return("Age") }
	if((covariate=="sex")){ return("Sex") }
	if((covariate=="meno")){ return("Sex & Menopause") }
	if((covariate=="meno_simple")){ return("Sex & Menopause") }

	if((covariate=="instrument")){ return(xlab("Instrument")) }



	return("NO LABEL")
}

get_strata_graph_info<-function(strata, levels=NA)
{
	returnval=list()
	if(strata=="sex")
	{	
		returnval$labels=c("female"="Female", "male"="Male")
		returnval$colours=c("female"="orange", "male"="blue")
		return(returnval)
	}
	if(strata=="meno")
	{	
		returnval$labels=c("male"="Male", "pre"="Pre-menopausal", "post"="Post-menopausal","hyst"="Hysterectomy", "unsure"="Unsure","noanswer"="No-answer")
		returnval$colours=c("pre"="indianred4", "post"="seagreen4", "male"="blue", "hyst"="black","unsure"="black","noanswer"="black")
		return(returnval)
	}
	if(strata=="meno_simple")
	{
		returnval$labels=c("male"="Male", "pre"="Pre-menopausal", "post"="Post-menopausal","hyst"="Hysterectomy")
		returnval$colours=c("pre"="indianred4", "post"="seagreen4", "male"="blue", "hyst"="gold3")
		return(returnval)
	}

	returnval$labels=levels
	returnval$colours=rep("grey", length(levels))
	return(returnval)
}

covariate_facet_labeller<-function(variable, value)
{
	return(get_strata_graph_info(variable)$label[value])
	
}
	



get_covariate_xscale<-function(covariate)
{
	if((covariate=="month_of_year_acq")||(covariate=="month_of_year_ext")){ return(scale_x_continuous(breaks=1:12))}
	if((covariate=="day_of_week_acq")||(covariate=="day_of_meek_ext")){ return(scale_x_continuous(breaks=1:7))}

	return(scale_x_continuous())
}



