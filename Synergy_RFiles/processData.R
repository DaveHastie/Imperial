# TODO: Add comment
# 
# Author: dhastie
###############################################################################


checkCategoryBoundaries<-function(dirPath='/home/dhastie/Imperial/Analyses/Synergy/Data/Raw',sex='all'){

	# Link to only those subjects that have all data
   subjectData<-read.csv(file.path(dirPath,'subject.csv'),header=T)
	subjectIDVec<-as.integer(subjectData$subjctid)
	tobaccoData<-read.csv(file.path(dirPath,'tobacco.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(tobaccoData$subjctid))
	dieselData<-read.csv(file.path(dirPath,'dme.csv'),header=T)	
	subjectIDVec<-intersect(subjectIDVec,as.integer(dieselData$subjctid))
	ppdData<-read.csv(file.path(dirPath,'ppd.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(ppdData$subjctid))
	exposureData<-read.csv(file.path(dirPath,'allexposures.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(exposureData$subjctid))
	etsHomeData<-read.csv(file.path(dirPath,'etsHome.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsHomeData$subjctid))
	etsWorkData<-read.csv(file.path(dirPath,'etsWork.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsWorkData$subjctid))
	
	
	
	subjectIDVec<-sort(subjectIDVec)
	subjectData<-subjectData[match(subjectIDVec,as.integer(subjectData$subjctid)),]
	tobaccoData<-tobaccoData[match(subjectIDVec,as.integer(tobaccoData$subjctid)),]
	dieselData<-dieselData[match(subjectIDVec,as.integer(dieselData$subjctid)),]
	ppdData<-ppdData[match(subjectIDVec,as.integer(ppdData$subjctid)),]
	exposureData<-exposureData[match(subjectIDVec,as.integer(exposureData$subjctid)),]
	etsHomeData<-etsHomeData[match(subjectIDVec,as.integer(etsHomeData$subjctid)),]
	etsWorkData<-etsWorkData[match(subjectIDVec,as.integer(etsWorkData$subjctid)),]

	
	# Cut down to the confounders we have here 
	keepSubjects<-which(subjectData$age>0)
	keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp>-1000))
	keepSubjects<-intersect(keepSubjects,which(subjectData$educ>=0))
	if(sex=='Men'){
   	keepSubjects<-intersect(keepSubjects,which(subjectData$sex==1))
	}else{
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==0))		
	}
   keepSubjects<-sort(keepSubjects)
	
	subjectData<-subjectData[keepSubjects,]
	tobaccoData<-tobaccoData[keepSubjects,]
	dieselData<-dieselData[keepSubjects,]
	ppdData<-ppdData[keepSubjects,]
	exposureData<-exposureData[keepSubjects,]
	etsHomeData<-etsHomeData[keepSubjects,]
	etsWorkData<-etsWorkData[keepSubjects,]
	
	# Tobacco
	browser()
	cat('Intensity\n')
   cat(paste(quantile(na.omit(tobaccoData$intensity[tobaccoData$intensity>0]),0.25),
					quantile(na.omit(tobaccoData$intensity[tobaccoData$intensity>0]),0.5),
					quantile(na.omit(tobaccoData$intensity[tobaccoData$intensity>0]),0.75),sep=' - '))
	cat('\n')

	cat('Duration\n')
   cat(paste(quantile(na.omit(tobaccoData$duration[tobaccoData$duration>0]),0.25),
					quantile(na.omit(tobaccoData$duration[tobaccoData$duration>0]),0.5),
					quantile(na.omit(tobaccoData$duration[tobaccoData$duration>0]),0.75),sep=' - '))
	cat('\n')
	
	cat('Time since cessation\n')
   cat(paste(quantile(na.omit(tobaccoData$smokestop_cig[tobaccoData$smokestop_cig>0]),0.33),
					quantile(na.omit(tobaccoData$smokestop_cig[tobaccoData$smokestop_cig>0]),0.66),sep=' - '))
	cat('\n')
	
	cat('PackYears\n')
	cat(paste(quantile(na.omit(dieselData$packyears[dieselData$packyears>0]),0.25),
					quantile(na.omit(dieselData$packyears[dieselData$packyears>0]),0.5),
					quantile(na.omit(dieselData$packyears[dieselData$packyears>0]),0.75),sep=' - '))
	cat('\n')

	cat('DME\n')
	cat(paste(quantile(na.omit(dieselData$cum_mean[dieselData$cum_mean>0]),0.25),
					quantile(na.omit(dieselData$cum_mean[dieselData$cum_mean>0]),0.5),
					quantile(na.omit(dieselData$cum_mean[dieselData$cum_mean>0]),0.75),sep=' - '))
	cat('\n')
	
	cat('PAH\n')
	cat(paste(quantile(na.omit(exposureData$PAH_cum0[exposureData$PAH_cum0>0]),0.25),
					quantile(na.omit(exposureData$PAH_cum0[exposureData$PAH_cum0>0]),0.5),
					quantile(na.omit(exposureData$PAH_cum0[exposureData$PAH_cum0>0]),0.75),sep=' - '))
	cat('\n')
	
	cat('Silica\n')
	cat(paste(quantile(na.omit(exposureData$silica_cum0[exposureData$silica_cum0>0]),0.25),
					quantile(na.omit(exposureData$silica_cum0[exposureData$silica_cum0>0]),0.5),
					quantile(na.omit(exposureData$silica_cum0[exposureData$silica_cum0>0]),0.75),sep=' - '))
	cat('\n')
	
	cat('Nickel\n')
	cat(paste(quantile(na.omit(exposureData$nickel_cum0[exposureData$nickel_cum0>0]),0.25),
					quantile(na.omit(exposureData$nickel_cum0[exposureData$nickel_cum0>0]),0.5),
					quantile(na.omit(exposureData$nickel_cum0[exposureData$nickel_cum0>0]),0.75),sep=' - '))
	cat('\n')
	
	cat('Asbestos\n')
	cat(paste(quantile(na.omit(exposureData$asbestos_cum0[exposureData$asbestos_cum0>0]),0.25),
					quantile(na.omit(exposureData$asbestos_cum0[exposureData$asbestos_cum0>0]),0.5),
					quantile(na.omit(exposureData$asbestos_cum0[exposureData$asbestos_cum0>0]),0.75),sep=' - '))
	cat('\n')
	
	cat('Chrome\n')
	cat(paste(quantile(na.omit(exposureData$chrom_cum0[exposureData$chrom_cum0>0]),0.25),
					quantile(na.omit(exposureData$chrom_cum0[exposureData$chrom_cum0>0]),0.5),
					quantile(na.omit(exposureData$chrom_cum0[exposureData$chrom_cum0>0]),0.75),sep=' - '))
	cat('\n')
	
	
}

# Discard histology -9999
# study center (study and country)
# Age discard any < 0
# Education - put 0 and 1 together, discard all < 0

# Intensity 0 10 20
# Duration 0 20 40
# Cessation time NS 15 5 0
# Packyears 0 20 40
# DME 10 40
# PAH 0.2 0.75
# Silica 1 4
# Nickel 0.01 0.06
# Asbestos 0.05 0.2
# Chrome 0.02 0.1

createTobaccoCategorical<-function(outFile,dirPath='/home/dhastie/Imperial/Analyses/Synergy/Data',sex='Men'){
	
	# Link to only those subjects that have all data
   subjectData<-read.csv(file.path(dirPath,'Raw','subject.csv'),header=T)
	subjectIDVec<-as.integer(subjectData$subjctid)
	tobaccoData<-read.csv(file.path(dirPath,'Raw','tobacco.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(tobaccoData$subjctid))
	dieselData<-read.csv(file.path(dirPath,'Raw','dme.csv'),header=T)	
	subjectIDVec<-intersect(subjectIDVec,as.integer(dieselData$subjctid))
	ppdData<-read.csv(file.path(dirPath,'Raw','ppd.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(ppdData$subjctid))
	exposureData<-read.csv(file.path(dirPath,'Raw','allexposures.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(exposureData$subjctid))
	etsHomeData<-read.csv(file.path(dirPath,'Raw','etsHome.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsHomeData$subjctid))
	etsWorkData<-read.csv(file.path(dirPath,'Raw','etsWork.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsWorkData$subjctid))
	
	
	
	subjectIDVec<-sort(subjectIDVec)
	subjectData<-subjectData[match(subjectIDVec,as.integer(subjectData$subjctid)),]
	tobaccoData<-tobaccoData[match(subjectIDVec,as.integer(tobaccoData$subjctid)),]
	dieselData<-dieselData[match(subjectIDVec,as.integer(dieselData$subjctid)),]
	ppdData<-ppdData[match(subjectIDVec,as.integer(ppdData$subjctid)),]
	exposureData<-exposureData[match(subjectIDVec,as.integer(exposureData$subjctid)),]
	etsHomeData<-etsHomeData[match(subjectIDVec,as.integer(etsHomeData$subjctid)),]
	etsWorkData<-etsWorkData[match(subjectIDVec,as.integer(etsWorkData$subjctid)),]
	
	
	# Cut down to the confounders we have here 
	keepSubjects<-which(subjectData$age>0)
	keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp>-1000))
	keepSubjects<-intersect(keepSubjects,which(subjectData$educ>=0))
	if(sex=='Men'){
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==1))
	}else{
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==0))		
	}
	keepSubjects<-sort(keepSubjects)
	
	subjectData<-subjectData[keepSubjects,]
	tobaccoData<-tobaccoData[keepSubjects,]
	dieselData<-dieselData[keepSubjects,]
	ppdData<-ppdData[keepSubjects,]
	exposureData<-exposureData[keepSubjects,]
	etsHomeData<-etsHomeData[keepSubjects,]
	etsWorkData<-etsWorkData[keepSubjects,]	

	nSubjects<-0
	
	nCovariates<-4
	covNames<-c('Intensity','Duration','CessationTime','PackYears')
	
	# Response
	outData<-subjectData$status
   outRawX<-NULL

	# Intensity
	tobaccoData$intensity[tobaccoData$intensity<0|is.na(tobaccoData$intensity)]<--999
	intVec<-with(tobaccoData,ifelse(smokecig==0,0,ifelse(intensity<0,-999,ifelse(intensity<=10,1,ifelse(intensity<=20,2,3)))))
	intVecRaw<-round(with(tobaccoData,ifelse(smokecig==0,0,ifelse(intensity<0,-999,intensity))),2)
   outData<-paste(outData,intVec)
   outRawX<-paste(outRawX,intVecRaw)
	
	# Duration
	tobaccoData$duration[tobaccoData$duration<0|is.na(tobaccoData$duration)]<--999
	durVec<-with(tobaccoData,ifelse(smokecig==0,0,ifelse(duration<0,-999,ifelse(duration<=20,1,ifelse(duration<40,2,3)))))
   durVecRaw<-round(with(tobaccoData,ifelse(smokecig==0,0,ifelse(duration<0,-999,duration))),2)
	outData<-paste(outData,durVec)
	outRawX<-paste(outRawX,durVecRaw)
   
	# Time since cessation
	tobaccoData$smokestopcig[tobaccoData$smokestop_cig<0|is.na(tobaccoData$smokestop_cig)]<--999
	quitVec<-with(tobaccoData,ifelse(smokecig==0,0,ifelse(smokecig==2,2,ifelse(smokestop_cig<0,-999,1))))
   quitVecRaw<-round(with(tobaccoData,ifelse(smokecig==0,100,ifelse(smokecig==2,0,ifelse(smokestop_cig<0,-999,smokestop_cig)))),2)
   outData<-paste(outData,quitVec)
   outRawX<-paste(outRawX,quitVecRaw)
	
	# Pack years
	dieselData$packyears[dieselData$packyears<0|is.na(dieselData$packyears)]<--999
	packyearsVec<-with(dieselData,ifelse(tobaccoData$smokecig==0,0,ifelse(packyears<0,-999,ifelse(packyears<=20,1,ifelse(packyears<=40,2,3)))))
	packyearsRaw<-round(with(dieselData,ifelse(tobaccoData$smokecig==0,0,ifelse(packyears<0,-999,packyears))),2)
   outData<-paste(outData,packyearsVec)
   outRawX<-paste(outRawX,packyearsRaw)
		
	# Confounders
	# Add in age
   confNames<-'Log_Age'
	confData<-log(subjectData$age)-mean(log(subjectData$age))
	
	# Add in education
	# If all 0 then have lowest level of educ
   confNames<-c(confNames,'Educ2','Educ3','Educ4')
	educVec<-with(subjectData,ifelse(educ<2,'0 0 0',ifelse(educ==2,'1 0 0',ifelse(educ==3,'0 1 0','0 0 1'))))
	confData<-paste(confData,educVec)
	
	# Finally add in the center
	uCenters<-unique(dieselData$centers)
	confNames<-c(confNames,paste('Center',uCenters[-1],sep='_'))
	nCenters<-length(uCenters)
	nSubjects<-nrow(subjectData)
	centersMat<-matrix(0,nrow=nSubjects,ncol=nCenters-1)
	centerVec<-rep('',nSubjects)
	for(i in 1:nSubjects){
		if(dieselData$centers[i]>uCenters[1]){
			centersMat[i,match(dieselData$centers[i],uCenters[-1])]<-1
		}
		centerVec[i]<-paste(centersMat[i,],collapse=' ')
	}
	confData<-paste(confData,centerVec)

	# Put it all together
	outData<-paste(outData,confData)
	outData<-c(nSubjects,4,covNames,4+nCenters-1,confNames,'4 4 3 4',outData)

	outFile<-file.path(dirPath,'Input',outFile)
	write(outData,outFile,append=F)
   outFileRawX<-gsub('.txt','_RawX.txt',outFile)
   write(c(nSubjects,4),file=outFileRawX,append=F,ncolumn=1)   
   write(covNames,file=outFileRawX,append=T,ncolumn=1)
   write(outRawX,outFileRawX,append=T,ncolumn=1)
   
   
	
}

# Intensity 0 10 20
# Duration 0 20 40
# Cessation time NS ES CS
# Packyears 0 20 40
# DME 10 40
# PAH 0.2 0.75
# Silica 1 4
# Nickel 0.01 0.06
# Asbestos 0.05 0.2
# Chrome 0.02 0.1
createFullCategorical<-function(outFile,dirPath='/home/dhastie/Imperial/Analyses/Synergy/Data',sex='Men',histology='all',categoricalY=F){
	
	# Link to only those subjects that have all data
   subjectData<-read.csv(file.path(dirPath,'Raw','subject.csv'),header=T)
	subjectIDVec<-as.integer(subjectData$subjctid)
	tobaccoData<-read.csv(file.path(dirPath,'Raw','tobacco.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(tobaccoData$subjctid))
	dieselData<-read.csv(file.path(dirPath,'Raw','dme.csv'),header=T)	
	subjectIDVec<-intersect(subjectIDVec,as.integer(dieselData$subjctid))
	ppdData<-read.csv(file.path(dirPath,'Raw','ppd.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(ppdData$subjctid))
	exposureData<-read.csv(file.path(dirPath,'Raw','allexposures.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(exposureData$subjctid))
	etsHomeData<-read.csv(file.path(dirPath,'Raw','etsHome.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsHomeData$subjctid))
	etsWorkData<-read.csv(file.path(dirPath,'Raw','etsWork.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsWorkData$subjctid))
		
	subjectIDVec<-sort(subjectIDVec)
	subjectData<-subjectData[match(subjectIDVec,as.integer(subjectData$subjctid)),]
	tobaccoData<-tobaccoData[match(subjectIDVec,as.integer(tobaccoData$subjctid)),]
	dieselData<-dieselData[match(subjectIDVec,as.integer(dieselData$subjctid)),]
	ppdData<-ppdData[match(subjectIDVec,as.integer(ppdData$subjctid)),]
	exposureData<-exposureData[match(subjectIDVec,as.integer(exposureData$subjctid)),]
	etsHomeData<-etsHomeData[match(subjectIDVec,as.integer(etsHomeData$subjctid)),]
	etsWorkData<-etsWorkData[match(subjectIDVec,as.integer(etsWorkData$subjctid)),]
	
	
	# Cut down to the confounders we have here 
	keepSubjects<-which(subjectData$age>0)
	keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp>-1000))
	keepSubjects<-intersect(keepSubjects,which(subjectData$educ>=0))
	if(sex=='Men'){
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==1))
	}else{
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==0))		
	}
	# Restrict to histology if necessary
	if(histology=='Adenocarcinoma'){
		keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp==3|subjectData$status==0))
	}else if(histology=='Squamous'){
		keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp==1|subjectData$status==0))
		
	}else if(histology=='SmallCell'){
		keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp==2|subjectData$status==0))		
	}

	if(categoricalY){
		# Ignore anything other than squamous, small cell or adenocarcinoma
		keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp<4))
	}
	
	keepSubjects<-sort(keepSubjects)
	
	
	
	subjectData<-subjectData[keepSubjects,]
	tobaccoData<-tobaccoData[keepSubjects,]
	dieselData<-dieselData[keepSubjects,]
	ppdData<-ppdData[keepSubjects,]
	exposureData<-exposureData[keepSubjects,]
	etsHomeData<-etsHomeData[keepSubjects,]
	etsWorkData<-etsWorkData[keepSubjects,]	
	
	nSubjects<-0
	
	nCovariates<-10
	covNames<-c('Intensity','Duration','CessationTime','PackYears','DME','PAH','Silica','NickelOrChromium','Asbestos','PPD')
	
	# Response
   if(categoricalY){
		# 1 = squamous, 2 = small cell, 3 = squamous
   	outData<-with(subjectData,ifelse(status==0,0,ifelse(histotyp<4,histotyp,4)))
	}else{
		outData<-subjectData$status		
	}
	
   # Intensity
	tobaccoData$intensity[tobaccoData$intensity<0|is.na(tobaccoData$intensity)]<--999
	intVec<-with(tobaccoData,ifelse(smokecig==0,0,ifelse(intensity<0,-999,ifelse(intensity<=10,1,ifelse(intensity<=20,2,3)))))
	outData<-paste(outData,intVec)
	
	# Duration
	tobaccoData$duration[tobaccoData$duration<0|is.na(tobaccoData$duration)]<--999
	durVec<-with(tobaccoData,ifelse(smokecig==0,0,ifelse(duration<0,-999,ifelse(duration<=20,1,ifelse(duration<40,2,3)))))
	outData<-paste(outData,durVec)
	
	# Time since cessation
	tobaccoData$smokestopcig[tobaccoData$smokestop_cig<0|is.na(tobaccoData$smokestop_cig)]<--999
	quitVec<-with(tobaccoData,ifelse(smokecig==0,0,ifelse(smokecig==2,2,ifelse(smokestop_cig<0,-999,1))))
	outData<-paste(outData,quitVec)
	
	# Pack years
	dieselData$packyears[dieselData$packyears<0|is.na(dieselData$packyears)]<--999
	packyearsVec<-with(dieselData,ifelse(tobaccoData$smokecig==0,0,ifelse(packyears<0,-999,ifelse(packyears<=20,1,ifelse(packyears<=40,2,3)))))
	outData<-paste(outData,packyearsVec)

	# DME
	dieselData$cum_mean[dieselData$cum_mean<0|is.na(dieselData$cum_mean)]<--999
	dmeVec<-with(dieselData,ifelse(dieselData$intdme==0,0,ifelse(cum_mean<0,-999,1)))
	outData<-paste(outData,dmeVec)
	
	# PAH
	exposureData$PAH_cum0[exposureData$PAH_cum0<0|is.na(exposureData$PAH_cum0)]<--999
	pahVec<-with(exposureData,ifelse(PAH_cum0==0,0,ifelse(PAH_cum0<0,-999,1)))
	outData<-paste(outData,pahVec)
	
	# Silica
	exposureData$silica_cum0[exposureData$silica_cum0<0|is.na(exposureData$silica_cum0)]<--999
	silicaVec<-with(exposureData,ifelse(silica_cum0==0,0,ifelse(silica_cum0<0,-999,1)))
	outData<-paste(outData,silicaVec)
	
	# Nickel / Chromium
	exposureData$nickel_cum0[exposureData$nickel_cum0<0|is.na(exposureData$nickel_cum0)]<--999
	exposureData$chrom_cum0[exposureData$chrom_cum0<0|is.na(exposureData$chrom_cum0)]<--999
   nickelChromVec<-with(exposureData,ifelse(nickel_cum0==0&chrom_cum0==0,0,ifelse(nickel_cum0>0|chrom_cum0>0,1,-999)))
   outData<-paste(outData,nickelChromVec)
	
	# Asbestos
	exposureData$asbestos_cum0[exposureData$asbestos_cum0<0|is.na(exposureData$asbestos_cum0)]<--999
	asbestosVec<-with(exposureData,ifelse(asbestos_cum0==0,0,ifelse(asbestos_cum0<0,-999,1)))
	outData<-paste(outData,asbestosVec)
		
   # OtherPPD
   ppdData$bronchitis[ppdData$bronchitis<0|is.na(ppdData$bronchitis)]<--999
	ppdData$emphysema[ppdData$emphysema<0|is.na(ppdData$emphysema)]<--999
	ppdData$tuberc[ppdData$tuberc<0|is.na(ppdData$tuberc)]<--999
	ppdData$asbestosis[ppdData$asbestosis<0|is.na(ppdData$asbestosis)]<--999
	ppdData$silicosis[ppdData$silicosis<0|is.na(ppdData$silicosis)]<--999
	ppdData$pneumonia[ppdData$pneumonia<0|is.na(ppdData$pneumonia)]<--999
	
	ppdVec<-with(ppdData,ifelse(bronchitis==1|emphysema==1|tuberc==1|asbestosis==1|silicosis==1|pneumonia==1,1,
                                    ifelse(bronchitis<0|emphysema<0|tuberc<0|asbestosis<0|silicosis<0|pneumonia<0,-999,0)))
	outData<-paste(outData,ppdVec)

	
	# Confounders
	# Add in age
   confNames<-'Log_Age'
	confData<-log(subjectData$age)-mean(log(subjectData$age))
	
	# Add in education
	# If all 0 then have lowest level of educ
   confNames<-c(confNames,'Educ2','Educ3','Educ4')
	educVec<-with(subjectData,ifelse(educ<2,'0 0 0',ifelse(educ==2,'1 0 0',ifelse(educ==3,'0 1 0','0 0 1'))))
	confData<-paste(confData,educVec)
		
	# Finally add in the center
	uCenters<-unique(dieselData$centers)
	confNames<-c(confNames,paste('Center',uCenters[-1],sep='_'))
	nCenters<-length(uCenters)
	nSubjects<-nrow(subjectData)
	centersMat<-matrix(0,nrow=nSubjects,ncol=nCenters-1)
	centerVec<-rep('',nSubjects)
	for(i in 1:nSubjects){
		if(dieselData$centers[i]>uCenters[1]){
			centersMat[i,match(dieselData$centers[i],uCenters[-1])]<-1
		}
		centerVec[i]<-paste(centersMat[i,],collapse=' ')
	}
	confData<-paste(confData,centerVec)
	
	# Farmer and list A jobs
	confNames<-c(confNames,'ListA','Farmer')
	confData<-paste(confData,dieselData$lista,dieselData$farmer)
	
	# Put it all together
	outData<-paste(outData,confData)
	outData<-c(nSubjects,nCovariates,covNames,6+nCenters-1,confNames,'4 4 3 4 2 2 2 2 2 2',outData)
	
	outFile<-file.path(dirPath,'Input',outFile)
	write(outData,outFile,append=F)
	
}

createNonSmokeCategorical<-function(outFile,dirPath='/home/dhastie/Imperial/Analyses/Synergy/Data',sex='Men'){
	
	# Link to only those subjects that have all data
   subjectData<-read.csv(file.path(dirPath,'Raw','subject.csv'),header=T)
	subjectIDVec<-as.integer(subjectData$subjctid)
	tobaccoData<-read.csv(file.path(dirPath,'Raw','tobacco.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(tobaccoData$subjctid))
	dieselData<-read.csv(file.path(dirPath,'Raw','dme.csv'),header=T)	
	subjectIDVec<-intersect(subjectIDVec,as.integer(dieselData$subjctid))
	ppdData<-read.csv(file.path(dirPath,'Raw','ppd.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(ppdData$subjctid))
	exposureData<-read.csv(file.path(dirPath,'Raw','allexposures.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(exposureData$subjctid))
	etsHomeData<-read.csv(file.path(dirPath,'Raw','etsHome.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsHomeData$subjctid))
	etsWorkData<-read.csv(file.path(dirPath,'Raw','etsWork.csv'),header=T)
	subjectIDVec<-intersect(subjectIDVec,as.integer(etsWorkData$subjctid))
	
	
	
	subjectIDVec<-sort(subjectIDVec)
	subjectData<-subjectData[match(subjectIDVec,as.integer(subjectData$subjctid)),]
	tobaccoData<-tobaccoData[match(subjectIDVec,as.integer(tobaccoData$subjctid)),]
	dieselData<-dieselData[match(subjectIDVec,as.integer(dieselData$subjctid)),]
	ppdData<-ppdData[match(subjectIDVec,as.integer(ppdData$subjctid)),]
	exposureData<-exposureData[match(subjectIDVec,as.integer(exposureData$subjctid)),]
	etsHomeData<-etsHomeData[match(subjectIDVec,as.integer(etsHomeData$subjctid)),]
	etsWorkData<-etsWorkData[match(subjectIDVec,as.integer(etsWorkData$subjctid)),]
	
	
	# Cut down to the confounders we have here 
	keepSubjects<-which(subjectData$age>0)
	keepSubjects<-intersect(keepSubjects,which(subjectData$histotyp>-1000))
	keepSubjects<-intersect(keepSubjects,which(subjectData$educ>=0))
	if(sex=='Men'){
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==1))
	}else{
		keepSubjects<-intersect(keepSubjects,which(subjectData$sex==0))		
	}
	# And restrict only to non-smokers with ETS data
	keepSubjects<-intersect(keepSubjects,which(tobaccoData$smokecig==0))
   keepSubjects<-intersect(keepSubjects,which(etsHomeData$home==1|etsHomeData$home==2))
   keepSubjects<-intersect(keepSubjects,which(etsWorkData$work==1|etsWorkData$work==2))

	keepSubjects<-sort(keepSubjects)
	
	subjectData<-subjectData[keepSubjects,]
	tobaccoData<-tobaccoData[keepSubjects,]
	dieselData<-dieselData[keepSubjects,]
	ppdData<-ppdData[keepSubjects,]
	exposureData<-exposureData[keepSubjects,]
	etsHomeData<-etsHomeData[keepSubjects,]
	etsWorkData<-etsWorkData[keepSubjects,]	
	
	nSubjects<-0
	
	nCovariates<-6
	covNames<-c('DME','PAH','Silica','NickelOrChromium','Asbestos','PPD')
	
	# Response
	outData<-subjectData$status
	
	# ETS Home
	etsHomeVec<-with(etsHomeData,ifelse(home==1,1,ifelse(home==2,0,-999)))

	# ETS Work
	etsWorkVec<-with(etsWorkData,ifelse(work==1,1,ifelse(work==2,0,-999)))
	
   	# DME
	dieselData$cum_mean[dieselData$cum_mean<0|is.na(dieselData$cum_mean)]<--999
	dmeVec<-with(dieselData,ifelse(dieselData$intdme==0,0,ifelse(cum_mean<0,-999,1)))
	outData<-paste(outData,dmeVec)
	
	# PAH
	exposureData$PAH_cum0[exposureData$PAH_cum0<0|is.na(exposureData$PAH_cum0)]<--999
	pahVec<-with(exposureData,ifelse(PAH_cum0==0,0,ifelse(PAH_cum0<0,-999,1)))
	outData<-paste(outData,pahVec)
	
	# Silica
	exposureData$silica_cum0[exposureData$silica_cum0<0|is.na(exposureData$silica_cum0)]<--999
	silicaVec<-with(exposureData,ifelse(silica_cum0==0,0,ifelse(silica_cum0<0,-999,1)))
	outData<-paste(outData,silicaVec)
	
	# Nickel / Chromium
	exposureData$nickel_cum0[exposureData$nickel_cum0<0|is.na(exposureData$nickel_cum0)]<--999
	exposureData$chrom_cum0[exposureData$chrom_cum0<0|is.na(exposureData$chrom_cum0)]<--999
   nickelChromVec<-with(exposureData,ifelse(nickel_cum0==0&chrom_cum0==0,0,ifelse(nickel_cum0>0|chrom_cum0>0,1,-999)))
   outData<-paste(outData,nickelChromVec)
	
	# Asbestos
	exposureData$asbestos_cum0[exposureData$asbestos_cum0<0|is.na(exposureData$asbestos_cum0)]<--999
	asbestosVec<-with(exposureData,ifelse(asbestos_cum0==0,0,ifelse(asbestos_cum0<0,-999,1)))
	outData<-paste(outData,asbestosVec)
	
	# PPD
   ppdData$bronchitis[ppdData$bronchitis<0|is.na(ppdData$bronchitis)]<--999
	ppdData$emphysema[ppdData$emphysema<0|is.na(ppdData$emphysema)]<--999
	ppdData$tuberc[ppdData$tuberc<0|is.na(ppdData$tuberc)]<--999
	ppdData$asbestosis[ppdData$asbestosis<0|is.na(ppdData$asbestosis)]<--999
	ppdData$silicosis[ppdData$silicosis<0|is.na(ppdData$silicosis)]<--999
	ppdData$pneumonia[ppdData$pneumonia<0|is.na(ppdData$pneumonia)]<--999
	
	ppdVec<-with(ppdData,ifelse(bronchitis==1|emphysema==1|tuberc==1|asbestosis==1|silicosis==1|pneumonia==1,1,
                                    ifelse(bronchitis<0|emphysema<0|tuberc<0|asbestosis<0|silicosis<0|pneumonia<0,-999,0)))
	outData<-paste(outData,ppdVec)

	# Confounders
	# Add in age
   confNames<-'Log_Age'
	confData<-log(subjectData$age)-mean(log(subjectData$age))
	
	# Add in education
	# If all 0 then have lowest level of educ
   confNames<-c(confNames,'Educ2','Educ3','Educ4')
	educVec<-with(subjectData,ifelse(educ<2,'0 0 0',ifelse(educ==2,'1 0 0',ifelse(educ==3,'0 1 0','0 0 1'))))
	confData<-paste(confData,educVec)
	
   confNames<-c(confNames,'ETSHome','ETSWork')
   confData<-paste(confData,etsHomeVec,etsWorkVec)
   
	# Finally add in the center
	uCenters<-unique(dieselData$centers)
	confNames<-c(confNames,paste('Center',uCenters[-1],sep='_'))
	nCenters<-length(uCenters)
	nSubjects<-nrow(subjectData)
	centersMat<-matrix(0,nrow=nSubjects,ncol=nCenters-1)
	centerVec<-rep('',nSubjects)
	for(i in 1:nSubjects){
		if(dieselData$centers[i]>uCenters[1]){
			centersMat[i,match(dieselData$centers[i],uCenters[-1])]<-1
		}
		centerVec[i]<-paste(centersMat[i,],collapse=' ')
	}
	confData<-paste(confData,centerVec)
	
	# Farmer and list A jobs
	confNames<-c(confNames,'ListA','Farmer')
	confData<-paste(confData,dieselData$lista,dieselData$farmer)
	
	# Put it all together
	outData<-paste(outData,confData)
	outData<-c(nSubjects,nCovariates,covNames,8+nCenters-1,confNames,'2 2 2 2 2 2',outData)
	
	outFile<-file.path(dirPath,'Input',outFile)
	write(outData,outFile,append=F)
	
}
