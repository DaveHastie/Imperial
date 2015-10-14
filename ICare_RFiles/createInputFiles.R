# Function to create the categorical data for the 3 main smoking factors
createCategoricalFile<-function(covariateNames=NULL,filePath="/home/dhastie/Imperial/Analyses/ICare/Data/Input",outFile="icare_cat_all.txt",smokersOnly=F,sex="Men",histology="All",responseAsCovariate=F,writeCSIFile=F,writeDMEFile=F){

   if(is.null(covariateNames)){
      covariateNames<-c('Intensity','Duration','CessationTime','StartAge','PackYears','DarkPackYears','NonFilteredPackYears')
   }
   
   ICareData<-read.csv(file.path(filePath,"IcareTabac_janv2011.csv"),quote='"',header=T)
   
   # Because Belfort is all cases, group with Haut-Rhin
   ICareData$depthab[ICareData$depthab==90]<-68
   
   # Men only
   if(sex=="Men"){
      ICareData<-ICareData[ICareData$sexe==1,]
   }else{
      ICareData<-ICareData[ICareData$sexe==2,]      
   }  
   
   # If necessary restrict to smokers only
   if(smokersOnly){
      ICareData<-ICareData[!is.na(ICareData$csi),]
      ICareData<-ICareData[ICareData$csi>0,]
   }
   
	cat("Initial no. rows: ")
	cat(paste(nrow(ICareData),sum(ICareData$catem)))
	cat("\n")
	
   ICareSubjectData<-read.csv(file.path(filePath,"subject_icare.csv"),quote='"',header=T)
   ICareSubjectData$educ[is.na(ICareSubjectData$educ)]<--999
   idVec<-ICareData$numid
   keepRows<-idVec%in%ICareSubjectData$orig_id
   ICareData<-ICareData[keepRows,]

	cat("Rows after no match to subject file removed: ")
	cat(paste(nrow(ICareData),sum(ICareData$catem)))
	cat("\n")
	
	
   # Remove any rows where the main smoking variables are all missing
   removeRows<-with(ICareData,which(is.na(nbcigjour)&is.na(durcig)&is.na(delaiarret)))
   ICareData<-ICareData[-removeRows,]
	cat("Rows after missing smoking removed: ")
	cat(paste(nrow(ICareData),sum(ICareData$catem)))
	cat("\n")
   
   if(histology!='All'){

      nSubjects<-nrow(ICareData)
      histologyVec<-rep(-999,nSubjects)
      for(i in 1:nSubjects){
         currentID<-ICareData$numid[i]
         if(ICareData$catem[i]==0){
            # No cancer
            histologyVec[i]<-0
         }else{
            histTmp<-ICareSubjectData$histotyp[which(ICareSubjectData$orig_id==currentID)]
            # 1=Squamous, 2=Small cell, 3=Adenocarcinoma
            histologyVec[i]<-ifelse(histTmp<0,-999,ifelse(histTmp>3,4,histTmp))
         }
      }
            
      if(histology=='Adenocarcinoma'){
         ICareData<-ICareData[histologyVec==3|histologyVec==0,]
      }else if(histology=='SmallCell'){
         ICareData<-ICareData[histologyVec==2|histologyVec==0,]         
      }else if(histology=='Squamous'){
         ICareData<-ICareData[histologyVec==1|histologyVec==0,]         
      }else if(histology=='SquamousOrSmallCell'){
         ICareData<-ICareData[histologyVec==1|histologyVec==2|histologyVec==0,]         
      }
   }
      
   ## Confounders
   # Age (centred, on log scale, in years)
   # Education
   
   tmpSubjData<-ICareSubjectData[ICareSubjectData$orig_id%in%ICareData$numid,]
   meanAge<-mean(na.omit(tmpSubjData$age[tmpSubjData$age>0]))
 
   ageVec<-rep(meanAge,nrow(ICareData))
   educMat<-matrix(0,nrow=nrow(ICareData),ncol=3)
   removeRows<-c()
   for(i in 1:nrow(ICareData)){
      currentID<-ICareData$numid[i]
      subjectID<-ICareSubjectData$subjctid[which(ICareSubjectData$orig_id==currentID)]
      ageVec[i]<-ICareSubjectData$age[which(ICareSubjectData$orig_id==currentID)]
      if(ICareSubjectData$educ[which(ICareSubjectData$orig_id==currentID)]>1){
         educMat[i,ICareSubjectData$educ[which(ICareSubjectData$orig_id==currentID)]-1]<-1
      }else if(ICareSubjectData$educ[which(ICareSubjectData$orig_id==currentID)]<0){
         removeRows<-c(removeRows,i)
      }
   }
   logAgeYears<-log(ageVec)-mean(log(ageVec))
   logAgeYears<-logAgeYears[-removeRows]
   ICareData<-ICareData[-removeRows,]
   educMat<-educMat[-removeRows,]
   
	cat("Rows after education removed: ")
	cat(paste(nrow(ICareData),sum(ICareData$catem)))
	cat("\n")
	
   # Remove any centers where only cases or controls
   uCenters<-sort(unique(as.integer(ICareData$depthab)))
   nCenters<-length(uCenters)
   
   for(i in 1:nCenters){
      centerRows<-which(ICareData$depthab==uCenters[i])
      if(sum(ICareData$catem[centerRows])==length(centerRows)||
                     sum(ICareData$catem[centerRows])==0){
         ICareData<-ICareData[-centerRows,]
         logAgeData[-centerRows]
         educMat[-centerRows,]
      }  
   }   

   cat("Rows after bad centers removed: ")
	cat(paste(nrow(ICareData),sum(ICareData$catem)))
	cat("\n")

   if(writeCSIFile){
		csiVec<-ICareData$csi
	}

   # Confounders for centre
   uCentres<-sort(unique(as.integer(ICareData$depthab)))
   uCentres<-uCentres[uCentres!=min(uCentres)]
   nCentres<-length(uCentres)
   centreMat<-matrix(0,nrow=nrow(ICareData),ncol=nCentres)
   for(j in 1:nCentres){
      centreMat[ICareData$depthab==uCentres[j],j]<-1
   }	 
   
   outData<-paste(ICareData$catem)
   outRawX<-NULL
   nCategories<-''
   covNames<-c()
   if(responseAsCovariate){
      outData<-paste(outData,ICareData$catem)
      nCategories<-'2'
      isOrdinal<-'0'
      covNames<-'Response'
   }
      
   if('Intensity'%in%covariateNames){
      # Number of cigarettes per day (intensity)
      cutoff1<-10
      cutoff2<-20
      cutoff3<-30
      cat(paste('Intensity:',0,cutoff1,cutoff2,cutoff3,'\n'))
      nCigPerDay<-with(ICareData,ifelse(is.na(nbcigjour),-999,ifelse(nbcigjour==0,0,ifelse(nbcigjour<cutoff1,1,ifelse(nbcigjour<cutoff2,2,ifelse(nbcigjour<cutoff3,3,4))))))
      nCigPerDayRaw<-round(with(ICareData,ifelse(is.na(nbcigjour),-999,nbcigjour)),2)
      if(smokersOnly){
         nCigPerDay[nCigPerDay>0]<-nCigPerDay[nCigPerDay>0]-1
      }
      outData<-paste(outData,nCigPerDay)
      outRawX<-paste(outRawX,nCigPerDayRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'Intensity')
   }
   
   # Duration of smoking
   if('Duration'%in%covariateNames){
      cutoff1<-20
      cutoff2<-30
      cutoff3<-40
      cat(paste('Duration:',0,cutoff1,cutoff2,cutoff3,'\n'))
      durCig<-with(ICareData,ifelse(is.na(durcig),-999,ifelse(durcig==0,0,ifelse(durcig<=cutoff1,1,ifelse(durcig<=cutoff2,2,ifelse(durcig<=cutoff3,3,4))))))
      durCigRaw<-round(with(ICareData,ifelse(is.na(durcig),-999,durcig)),2)
      if(smokersOnly){
         durCig[durCig>0]<-durCig[durCig>0]-1
      }
      outData<-paste(outData,durCig)
      outRawX<-paste(outRawX,durCigRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'Duration')
   }
   
   # Time since cessation of smoking
   if('CessationTime'%in%covariateNames){
      cutoff1<-10
      cutoff2<-20
      cat(paste('Cessation time:',0,cutoff1,cutoff2,'Current','\n'))
         
      ceaseCig<-with(ICareData,ifelse(is.na(delaiarret),-999,ifelse(delaiarret==0&durCig==0,0,ifelse(delaiarret==0,4,ifelse(delaiarret<=cutoff1,3,ifelse(delaiarret<=cutoff2,2,1))))))
      ceaseCigRaw<-round(with(ICareData,ifelse(is.na(delaiarret),-999,delaiarret)),2)
      if(smokersOnly){
         ceaseCig[ceaseCig>0]<-ceaseCig[ceaseCig>0]-1
      }   
      outData<-paste(outData,ceaseCig)
      outRawX<-paste(outRawX,ceaseCigRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'CessationTime')
   }
   
   if('StartAge'%in%covariateNames){
      cutoff1<-15
      cutoff2<-17
      cutoff3<-20
      cat(paste('Start age:',0,cutoff1,cutoff2,cutoff3,'\n'))
      startAge<-with(ICareData,ifelse(is.na(csi),-999,ifelse(is.na(agedebutcig)&csi>0,-999,ifelse(csi==0,0,ifelse(agedebutcig<cutoff1,4,ifelse(agedebutcig>=cutoff1&agedebutcig<cutoff2,3,ifelse(agedebutcig>=cutoff2&agedebutcig<cutoff3,2,1)))))))
      startAgeRaw<-round(with(ICareData,ifelse(is.na(csi),-999,ifelse(is.na(agedebutcig)&csi>0,-999,ifelse(csi==0,0,agedebutcig)))),2)
      if(smokersOnly){
         startAge[startAge>0]<-startAge[startAge>0]-1
      }         
      outData<-paste(outData,startAge)
      outRawX<-paste(outRawX,startAgeRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'StartAge')
   }

   if('PackYears'%in%covariateNames){
      cutoff1<-15
      cutoff2<-30
      cutoff3<-45
      cat(paste('Pack years:','Non',cutoff1,cutoff2,cutoff3,'\n'))
      packYears<-with(ICareData,ifelse(is.na(csi),-999,ifelse(is.na(packyear),-999,ifelse(csi==0,0,ifelse(packyear<cutoff1,1,ifelse(packyear<cutoff2,2,ifelse(packyear<cutoff3,3,4)))))))
      packYearsRaw<-round(with(ICareData,ifelse(is.na(csi),-999,ifelse(is.na(packyear),-999,packyear))),2)
      if(smokersOnly){
         packYears[packYears>0]<-packYears[packYears>0]-1
      }         
      outData<-paste(outData,packYears)
      outRawX<-paste(outRawX,packYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'PackYears')
   }
      
   if('DarkPackYears'%in%covariateNames||
         'NonFilteredPackYears'%in%covariateNames){
      ICareExtraData<-read.csv(file.path(filePath,"IcareTabac2_janv2011.csv"),quote='"',header=T)
      nSubjects<-nrow(ICareData)
      darkPackYears<-rep(NA,nSubjects)
      nonFilteredPackYears<-rep(NA,nSubjects)         
      for(i in 1:nSubjects){
         currentID<-ICareData$numid[i]
         matchingRows<-ICareExtraData[which(ICareExtraData$numid==currentID),]
         # Check if this is a non-smoker
         if(all(is.na(matchingRows$semaine)&is.na(matchingRows$jour))){            
            darkPackYears[i]<-0
            nonFilteredPackYears[i]<-0
            next
         }         
         nPeriods<-nrow(matchingRows)
         # Check if a smoker but with missing data
         if(nPeriods==0){
            darkPackYears[i]<-NA
            nonFilteredPackYears<-NA
            next
         }
         cigsPerYear<-ifelse(is.na(matchingRows$semaine)&is.na(matchingRows$jour),NA,ifelse(!is.na(matchingRows$semaine),as.numeric(matchingRows$semaine)*(365/7),as.numeric(matchingRows$jour)*365))
         numYears<-ifelse(is.na(matchingRows$a),ageVec[i],as.numeric(matchingRows$a))-as.numeric(matchingRows$de)
         numYears<-ifelse(numYears==0,1,numYears)
         cigsPerPeriod<-numYears*cigsPerYear
         matchingRows<-matchingRows[!is.na(cigsPerPeriod),]
         cigsPerPeriod<-cigsPerPeriod[!is.na(cigsPerPeriod)]
         if(nrow(matchingRows)==0){
            darkPackYears[i]<-NA
            nonFilteredPackYears[i]<-NA
            next            
         }
         darkPeriods<-which(matchingRows$tabac==2)
         nonFilteredPeriods<-which(matchingRows$filtre==2)
         darkPackYears[i]<-sum(cigsPerPeriod[darkPeriods])/(20*365)
         nonFilteredPackYears[i]<-sum(cigsPerPeriod[nonFilteredPeriods])/(20*365)
      }
   }
   
      

   if('DarkPackYears'%in%covariateNames){
      cutoff1<-5
      cutoff2<-20
      cutoff3<-40
      cat(paste('DarkPackYears:','Non',cutoff1,cutoff2,cutoff3,'\n'))
      darkPackYearsRaw<-round(ifelse(is.na(ICareData$csi),-999,ifelse(is.na(darkPackYears),-999,darkPackYears)),2)
      darkPackYears<-ifelse(is.na(ICareData$csi),-999,ifelse(is.na(darkPackYears),-999,ifelse(ICareData$csi==0,0,ifelse(darkPackYears<cutoff1,1,ifelse(darkPackYears<cutoff2,2,ifelse(darkPackYears<cutoff3,3,4))))))
      if(smokersOnly){
         darkPackYears[darkPackYears>0]<-darkPackYears[darkPackYears>0]-1
      }          
      outData<-paste(outData,darkPackYears)
      outRawX<-paste(outRawX,darkPackYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'DarkPackYears')
   }
      
   if('NonFilteredPackYears'%in%covariateNames){
      cutoff1<-5
      cutoff2<-15
      cutoff3<-35
      cat(paste('NonFilteredPackYears:','Non',cutoff1,cutoff2,cutoff3,'\n'))
      nonFilteredPackYearsRaw<-round(ifelse(is.na(ICareData$csi),-999,ifelse(is.na(nonFilteredPackYears),-999,nonFilteredPackYears)),2)
      nonFilteredPackYears<-ifelse(is.na(ICareData$csi),-999,ifelse(is.na(nonFilteredPackYears),-999,ifelse(ICareData$csi==0,0,ifelse(nonFilteredPackYears<cutoff1,1,ifelse(nonFilteredPackYears<cutoff2,2,ifelse(nonFilteredPackYears<cutoff3,3,4))))))
      if(smokersOnly){
         nonFilteredPackYears[nonFilteredPackYears>0]<-nonFilteredPackYears[nonFilteredPackYears>0]-1
      }         
      outData<-paste(outData,nonFilteredPackYears)
      outRawX<-paste(outRawX,nonFilteredPackYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'NonFilteredPackYears')
   }
   
   # Write the output
   outFile<-file.path(filePath,outFile)
   outFileRawX<-gsub('.txt','_RawX.txt',outFile)
	if(writeCSIFile){
		outFileCSI<-gsub('.txt','_CSI.txt',outFile)
		write(ifelse(is.na(csiVec),-999,csiVec),file=outFileCSI,ncolumn=1,append=F)
	}
   write(c(nrow(ICareData),length(covNames)),file=outFile,append=F,ncolumn=1)   
   write(covNames,file=outFile,append=T,ncolumn=1)
   write(c(nrow(ICareData),length(covNames)),file=outFileRawX,append=F,ncolumn=1)   
   write(covNames,file=outFileRawX,append=T,ncolumn=1)
   write(4+nCentres,file=outFile,append=T,ncolumn=1)               
   
   write(c("Age","LowEduc","AvgEduc","HighEduc",paste("Dept",uCentres,sep=""),nCategories),
            file=outFile,append=T,ncolumn=1)     
   outData<-paste(outData,logAgeYears)#,listAVec)
   for(j in 1:3){
      outData<-paste(outData,educMat[,j])
   }
   
   for(j in 1:nCentres){
      outData<-paste(outData,centreMat[,j])
   }
   
   write(outData,file=outFile,append=T,ncolumn=1)
   write(outRawX,file=outFileRawX,append=T,ncolumn=1)
   
}

# Function to create the categorical data for the 3 main smoking factors
createOtherStudyFile<-function(filePath="/home/dhastie/iCare/data/input",smokersOnly=F,sex="Men",histology="All",responseAsCovariate=F){
   
   covariateNames<-c('Intensity','Duration','CessationTime','PackYears')
   
   subjectData<-read.csv(file.path(filePath,"synergy","subject.csv"))
   tobaccoData<-read.csv(file.path(filePath,"synergy","tobacco.csv"))
   dmeData<-read.csv(file.path(filePath,"synergy","dme.csv"))
   
   # Restrict only to subjects where we have tobacco info
   tobaccoRowInd<-match(subjectData$subjctid,tobaccoData$subjctid,nomatch=NA)
   dmeRowInd<-match(subjectData$subjctid,dmeData$subjctid,nomatch=NA)
   
   subjectData<-subjectData[!(is.na(tobaccoRowInd)|is.na(dmeRowInd)),]
   # Reorder the tobacco and DME data so it is in the same order as subject data
   tobaccoData<-tobaccoData[match(subjectData$subjctid,tobaccoData$subjctid),]
   dmeData<-dmeData[match(subjectData$subjctid,dmeData$subjctid),]
   
   # Men only
   if(sex=="Men"){
      maleInd<-which(subjectData$sex==1)
      subjectData<-subjectData[maleInd,]
      tobaccoData<-tobaccoData[maleInd,]
      dmeData<-dmeData[maleInd,]
   }else{
      femaleInd<-which(subjectData$sex==0)
      subjectData<-subjectData[femaleInd,]
      tobaccoData<-tobaccoData[femaleInd,]
      dmeData<-dmeData[femaleInd,]
   }  
   
   # If necessary restrict to smokers only
   if(smokersOnly){
      smokeInd<-which(tobaccoData$smoke>0)
      subjectData<-subjectData[smokeInd,]
      tobaccoData<-tobaccoData[smokeInd,]
      dmeData<-dmeData[smokeInd,]
   }
   
   # Remove any missing centernames
   missingCenterInd<-is.na(subjectData$centername)|(nchar(gsub(' ','',subjectData$centername))==0)
   subjectData<-subjectData[!missingCenterInd,]
   tobaccoData<-tobaccoData[!missingCenterInd,]
   dmeData<-dmeData[!missingCenterInd,]
   
   # Make sure all centers have some cases and controls
   uCenters<-unique(subjectData$centername)
   nCenters<-length(uCenters)
   for(i in 1:nCenters){
      centerRows<-which(subjectData$centername==uCenters[i])
      if(sum(subjectData$status[centerRows])==length(centerRows)||
         sum(subjectData$status[centerRows])==0){
         subjectData<-subjectData[-centerRows,]
         tobaccoData<-tobaccoData[-centerRows,]
         dmeData<-dmeData[-centerRows,]
      }  
   }   

   
   if(histology!='All'){
      histologyVec<-ifelse(subjectData$status==0,0,ifelse(subjectData$histotyp<0,-999,ifelse(subjectData$histotyp>3,4,subjectData$histotyp)))
      
      if(histology=='Adenocarcinoma'){
         subjectData<-subjectData[histologyVec==3|histologyVec==0,]
         tobaccoData<-tobaccoData[histologyVec==3|histologyVec==0,]
         dmeData<-dmeData[histologyVec==3|histologyVec==0,]
      }else if(histology=='SmallCell'){
         subjectData<-subjectData[histologyVec==2|histologyVec==0,]         
         tobaccoData<-tobaccoData[histologyVec==2|histologyVec==0,]
         dmeData<-dmeData[histologyVec==2|histologyVec==0,]
      }else if(histology=='Squamous'){
         subjectData<-subjectData[histologyVec==1|histologyVec==0,]         
         tobaccoData<-tobaccoData[histologyVec==1|histologyVec==0,]
         dmeData<-dmeData[histologyVec==1|histologyVec==0,]
      }else if(histology=='SquamousOrSmallCell'){
         subjectData<-subjectData[histologyVec==1|histologyVec==2|histologyVec==0,]         
         tobaccoData<-tobaccoData[histologyVec==1|histologyVec==2|histologyVec==0,]
         dmeData<-dmeData[histologyVec==1|histologyVec==2|histologyVec==0,]
      }
   }
   
   ## Confounders
   # Age (centred, on log scale, in years)
   # Education
   nSubjects<-nrow(subjectData)

   meanAge<-mean(na.omit(subjectData$age[subjectData$age>0]))
   ageVec<-ifelse(is.na(subjectData$age)|subjectData$age<=0,meanAge,subjectData$age)
   logAgeYears<-log(ageVec)-mean(log(ageVec))
   educVec<-ifelse(is.na(subjectData$educ)|subjectData$educ<0|subjectData$educ>4,2,
                  ifelse(subjectData$educ<2,0,subjectData$educ-1))
   educMat<-matrix(0,ncol=3,nrow=length(educVec))
   for(i in 1:nrow(educMat)){
      if(educVec[i]>0){
         educMat[i,educVec[i]]<-1
      }
   }
      
   # Confounders for centre
   centreVec<-subjectData$centername
   centreVec<-gsub(' ','',centreVec)
   centreVec<-gsub('/','-',centreVec)
   uCentres<-sort(unique(centreVec))   
   nCentres<-length(uCentres)
   centreMat<-matrix(0,nrow=nSubjects,ncol=nCentres)
   for(j in 1:nCentres){
      centreMat[centreVec==uCentres[j],j]<-1
   }   
   
   outData<-paste(subjectData$status)
   outRawX<-NULL
   nCategories<-''
   covNames<-c()
   if(responseAsCovariate){
      outData<-paste(outData,subjectData$status)
      nCategories<-'2'
      isOrdinal<-'0'
      covNames<-'Response'
   }
   
   # Number of cigarettes per day (intensity)
   cutoff1<-10
   cutoff2<-20
   cutoff3<-30
   cat(paste('Intensity:',0,cutoff1,cutoff2,cutoff3,'\n'))
   nCigPerDay<-with(tobaccoData,ifelse(is.na(smoke),-999,ifelse(smoke==0,0,ifelse(is.na(intensity)|intensity<0,-999,ifelse(intensity==0,0,ifelse(intensity<cutoff1,1,ifelse(intensity<cutoff2,2,ifelse(intensity<cutoff3,3,4))))))))
   nCigPerDayRaw<-round(with(tobaccoData,ifelse(is.na(smoke),-999,ifelse(smoke==0,0,ifelse(is.na(intensity)|intensity<0,-999,intensity)))),2)
   if(smokersOnly){
      nCigPerDay[nCigPerDay>0]<-nCigPerDay[nCigPerDay>0]-1
   }
   outData<-paste(outData,nCigPerDay)
   outRawX<-paste(outRawX,nCigPerDayRaw,sep=ifelse(is.null(outRawX),'',' '))
   nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
   covNames<-c(covNames,'Intensity')
   
   # Years of smoking (duration)
   cutoff1<-20
   cutoff2<-30
   cutoff3<-40
   cat(paste('Duration:',0,cutoff1,cutoff2,cutoff3,'\n'))
   durCig<-with(tobaccoData,ifelse(is.na(smoke),-999,ifelse(smoke==0,0,ifelse(is.na(duration)|duration<0,-999,ifelse(duration==0,0,ifelse(duration<=cutoff1,1,ifelse(duration<=cutoff2,2,ifelse(duration<=cutoff3,3,4))))))))
   durCigRaw<-round(with(tobaccoData,ifelse(is.na(smoke),-999,ifelse(smoke==0,0,ifelse(is.na(duration)|duration<0,-999,duration)))),2)
   if(smokersOnly){
      durCig[durCig>0]<-durCig[durCig>0]-1
   }
   outData<-paste(outData,durCig)
   outRawX<-paste(outRawX,durCigRaw,sep=ifelse(is.null(outRawX),'',' '))
   nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
   covNames<-c(covNames,'Duration')
   
   
   # Years since quit smoking
   cutoff1<-10
   cutoff2<-20
   cat(paste('Cessation time:',0,cutoff1,cutoff2,'Current','\n'))
      
   ceaseCig<-with(tobaccoData,ifelse(is.na(smoke),-999,ifelse(smoke==0,0,ifelse(is.na(smokestop_cig)|smokestop_cig<0,-999,ifelse(smokestop_cig==0,4,ifelse(smokestop_cig<=cutoff1,3,ifelse(smokestop_cig<=cutoff2,2,1)))))))
   ceaseCigRaw<-round(with(tobaccoData,ifelse(is.na(smoke),-999,ifelse(smoke==0,0,ifelse(is.na(smokestop_cig)|smokestop_cig<0,-999,smokestop_cig)))),2)
   if(smokersOnly){
      ceaseCig[ceaseCig>0]<-ceaseCig[ceaseCig>0]-1
   }   
   outData<-paste(outData,ceaseCig)
   outRawX<-paste(outRawX,ceaseCigRaw,sep=ifelse(is.null(outRawX),'',' '))
   nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
   covNames<-c(covNames,'CessationTime')
   
   
   cutoff1<-15
   cutoff2<-30
   cutoff3<-45
   cat(paste('Pack years:','Non',cutoff1,cutoff2,cutoff3,'\n'))
   packYears<-with(dmeData,ifelse(is.na(tobaccoData$smoke),-999,ifelse(tobaccoData$smoke==0,0,ifelse(is.na(packyears)|packyears<0,-999,ifelse(packyears==0,0,ifelse(packyears<=cutoff1,1,ifelse(packyears<=cutoff2,2,ifelse(packyears<=cutoff3,3,4))))))))
   packYearsRaw<-round(with(dmeData,ifelse(is.na(tobaccoData$smoke),-999,ifelse(tobaccoData$smoke==0,0,ifelse(is.na(packyears)|packyears<0,-999,packyears)))),2)
   
   if(smokersOnly){
      packYears[packYears>0]<-packYears[packYears>0]-1
   }         
   outData<-paste(outData,packYears)
   outRawX<-paste(outRawX,packYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
   nCategories<-paste(nCategories,ifelse(smokersOnly,4,5),sep=ifelse(nchar(nCategories)>0,' ',''))
   covNames<-c(covNames,'PackYears')
   
   # Remove any subjects who are missing all covariates                    
   deleteInd<-(nCigPerDay<0)&(durCig<0)&(ceaseCig<0)&(packYears<0)   
   subjectData<-subjectData[!deleteInd,]
   outData<-outData[!deleteInd]
   outRawX<-outRawX[!deleteInd]
   logAgeYears<-logAgeYears[!deleteInd]
   educMat<-educMat[!deleteInd,]
   centreMat<-centreMat[!deleteInd,]                    
                       
                       
   # Write the output
   studyVec<-gsub(' ','',subjectData$study_name)
   uStudies<-sort(unique(studyVec))
   nStudies<-length(uStudies)
   for(i in 1:nStudies){
      tmpOutData<-outData[studyVec==uStudies[i]]
      tmpOutRawX<-outRawX[studyVec==uStudies[i]]
      outFile<-file.path(filePath,paste('synergy_',uStudies[i],'_main_packyears.txt',sep=''))
      outFileRawX<-gsub('.txt','_RawX.txt',outFile)
      write(c(length(tmpOutData),length(covNames)),file=outFile,append=F,ncolumn=1)   
      write(covNames,file=outFile,append=T,ncolumn=1)
      write(c(length(tmpOutData),length(covNames)),file=outFileRawX,append=F,ncolumn=1)   
      write(covNames,file=outFileRawX,append=T,ncolumn=1)
      tmpCentreMat<-as.matrix(centreMat[studyVec==uStudies[i],])
      keepInd<-which(apply(tmpCentreMat,2,sum)>0)
      keepInd<-keepInd[-1]
      write(4+length(keepInd),file=outFile,append=T,ncolumn=1)               
      if(length(keepInd)>0){   
         write(c("Age","LowEduc","AvgEduc","HighEduc",paste("Centre-",uCentres[keepInd],sep=""),nCategories),
               file=outFile,append=T,ncolumn=1)
      }else{
         write(c("Age","LowEduc","AvgEduc","HighEduc",nCategories),
               file=outFile,append=T,ncolumn=1)         
      }
      tmpOutData<-paste(tmpOutData,logAgeYears[studyVec==uStudies[i]],listAVec[studyVec==uStudies[i]])
      for(j in 1:3){
         tmpOutData<-paste(tmpOutData,educMat[studyVec==uStudies[i],j])
      }
   
      
      for(j in 1:nCentres){
         if(j%in%keepInd){
            tmpOutData<-paste(tmpOutData,tmpCentreMat[,j])
         }
      }
   
      write(tmpOutData,file=outFile,append=T,ncolumn=1)
      write(tmpOutRawX,file=outFileRawX,append=T,ncolumn=1)
   }
}


