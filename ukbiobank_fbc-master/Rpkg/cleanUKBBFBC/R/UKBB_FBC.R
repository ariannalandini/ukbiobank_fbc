
## author: William John Astle
## author email: wja24@cam.ac.uk/will@astle.net
## Description: Prepare the dataset and set environment variables

.onLoad<-function(libname, pkgname)
{
	
	UKBIOBANK_DATA_DIR<<-Sys.getenv("UKBIOBANK_DATA_DIR")
	UKBIOBANK_FBC_OUTPUT_DIR<<-Sys.getenv("UKBIOBANK_FBC_OUTPUT_DIR")
	print("Is central sds used")
	CENTRAL_SDS<<-3.5
}

merge_UKBB_data<-function()
{
	
	print("replace this with a directory PATH")
	
	# read in link tables 
	link13745_987=read.csv(sprintf("%s/link_tables/ukb4759.csv", UKBIOBANK_DATA_DIR))
	link13745_7439=read.csv(sprintf("%s/link_tables/ukb4760.csv", UKBIOBANK_DATA_DIR))


	# datasets for application 13745
	# note that datsets 6835  and 4735 contains no data which are different from 8893
	source("/home/wja24/data/UK_biobank/8893_in_13745/ukb8893.r", local=TRUE)
	bd_8893=bd
	rm(bd)
	source("/home/wja24/data/UK_biobank/4755_in_13745/ukb4755.r",local=TRUE)
	bd_4755=bd[,c("f.eid",colnames(bd)[!is.element(colnames(bd),colnames(bd_8893))])]
	rm(bd)

	# datasets for application 987
	source("/home/wja24/data/UK_biobank/4286_in_987/ukb4286.r",local=TRUE)
	bd_4286=bd[,c("f.eid",colnames(bd)[!is.element(colnames(bd),c(colnames(bd_4755),colnames(bd_8893)))])]
	bd_4286$f.eid=link13745_987$Encoded.ID.ANID13745[match(bd_4286$f.eid,link13745_987$Encoded.ID.ANID987)]
	bd_4286=bd_4286[!is.na(bd_4286$f.eid),]
	rm(bd)

	source("/home/wja24/data/UK_biobank/2204_in_987/ukb2204.r",local=TRUE)
	bd_2204=bd[,c("f.eid",colnames(bd)[!is.element(colnames(bd),c(colnames(bd_4286),colnames(bd_4755),colnames(bd_8893)))])]
	bd_2204$f.eid=link13745_987$Encoded.ID.ANID13745[match(bd_2204$f.eid,link13745_987$Encoded.ID.ANID987)]
	bd_2204=bd_2204[!is.na(bd_2204$f.eid),]

	rm(bd)
	source("/home/wja24/data/UK_biobank/1764_in_987/ukb1764.r",local=TRUE)
	bd_1764=bd[,c("f.eid",colnames(bd)[!is.element(colnames(bd),c(colnames(bd_2204),colnames(bd_4286),colnames(bd_4755),colnames(bd_8893)))])]
	bd_1764$f.eid=link13745_987$Encoded.ID.ANID13745[match(bd_1764$f.eid,link13745_987$Encoded.ID.ANID987)]
	bd_1764=bd_1764[!is.na(bd_1764$f.eid),]
	rm(bd)

	## merge the data
	bd=merge(merge(merge(merge(bd_1764, bd_2204, by="f.eid", all=TRUE),bd_4286, by="f.eid",all=TRUE), bd_4755, by="f.eid", all=TRUE), bd_8893, by="f.eid", all=TRUE)
	bd$f.eid_7439=link13745_7439$Encoded.ID.ANID7439[match(bd$f.eid,link13745_7439$Encoded.ID.ANID13745)]

	## order the columns by label
	bd=bd[,c(1,dim(bd)[2],1+order(as.numeric(gsub("f.(.+)\\.(.+)\\.(.+)","\\1.\\2\\3",colnames(bd)[2:(dim(bd)[2]-1)]))))]

	# read in the link table to application 7439 which corresponds to Cambridge genotype release
	link13745_7439=read.csv(sprintf("%s/link_tables/ukb4760.csv", UKBIOBANK_DATA_DIR))

	# Add id column for application 7439
	bd$f.eid_7439=link13745_7439$Encoded.ID.ANID7439[match(bd$f.eid,link13745_7439$Encoded.ID.ANID13745)]

	# remove withdrawn individuals 
	withdrawn=read.table(sprintf("%s/withdrawn/13745/w1374_20170726.csv", UKBIOBANK_DATA_DIR)) 
	bd=bd[!is.element(bd$f.eid, withdrawn),]

	bd$FID=NULL
	bd$IID=NULL

	save(list=c("bd"),file=sprintf("%s/Rdata/ret_ukbb_fbc_merged.Rdata", UKBIOBANK_DATA_DIR))

}

load_merged_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_merged.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
}

load_prepared_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_prepared.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
	
}

load_tech_adjusted_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_tech_adjusted.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
	
}

load_gwas_adjusted_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_gwas_adjusted.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
	
}

load_gwas_prepared_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_gwas_prepared.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
	
}

load_gwas_prepared_original_scale_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_gwas_prepared_original_scale.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
	
}

load_ancillary_UKBB_data<-function()
{
	load(sprintf("%s/Rdata/ret_ukbb_fbc_ancillary.Rdata", UKBIOBANK_DATA_DIR), envir=.GlobalEnv)
	
}




