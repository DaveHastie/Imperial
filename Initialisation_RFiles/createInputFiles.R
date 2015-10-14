# Function to create the categorical data for the 3 main smoking factors
createCategoricalFile<-function(covariateNames=NULL,nSubjects=NULL,filePath="/home/dhastie/Imperial/Analyses/Initialisation/Data/Input",outFile="real_discrete_bernoulli.txt"){

   if(is.null(covariateNames)){
      covariateNames<-c('Intensity','Duration','CessationTime','PackYears','DarkPackYears','NonFilteredPackYears')
   }
   
   ICareData<-read.csv(file.path(filePath,"IcareTabac_janv2011.csv"),quote='"',header=T)
   
   # Because Belfort is all cases, group with Haut-Rhin
   ICareData$depthab[ICareData$depthab==90]<-68
   
   # Men only
   ICareData<-ICareData[ICareData$sexe==1,]
   
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
	
	if(!is.null(nSubjects)){
		nCase<-as.integer(nSubjects*sum(ICareData$catem)/nrow(ICareData))
		nControl<-nSubjects-nCase
		caseInd<-which(ICareData$catem==1)
		controlInd<-which(ICareData$catem==0)
		keepSubjects<-sort(c(sample(caseInd,nCase,replace=F),sample(controlInd,nControl,replace=F)))
		ICareData<-ICareData[keepSubjects,]
      logAgeYears<-logAgeYears[keepSubjects]
      educMat<-educMat[keepSubjects,]
	}
	
   cat("Rows after meeting nSubjects criterion: ")
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
         logAgeYears<-logAgeYears[-centerRows]
         educMat<-educMat[-centerRows,]
      }  
   }   

   cat("Rows after bad centers removed: ")
	cat(paste(nrow(ICareData),sum(ICareData$catem)))
	cat("\n")

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
      
   if('Intensity'%in%covariateNames){
      # Number of cigarettes per day (intensity)
      cutoff1<-10
      cutoff2<-20
      cutoff3<-30
      cat(paste('Intensity:',0,cutoff1,cutoff2,cutoff3,'\n'))
      nCigPerDay<-with(ICareData,ifelse(is.na(nbcigjour),-999,ifelse(nbcigjour==0,0,ifelse(nbcigjour<cutoff1,1,ifelse(nbcigjour<cutoff2,2,ifelse(nbcigjour<cutoff3,3,4))))))
      nCigPerDayRaw<-round(with(ICareData,ifelse(is.na(nbcigjour),-999,nbcigjour)),2)
      outData<-paste(outData,nCigPerDay)
      outRawX<-paste(outRawX,nCigPerDayRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,5,sep=ifelse(nchar(nCategories)>0,' ',''))
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
      outData<-paste(outData,durCig)
      outRawX<-paste(outRawX,durCigRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,5,sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'Duration')
   }
   
   # Time since cessation of smoking
   if('CessationTime'%in%covariateNames){
      cutoff1<-10
      cutoff2<-20
      cat(paste('Cessation time:',0,cutoff1,cutoff2,'Current','\n'))
         
      ceaseCig<-with(ICareData,ifelse(is.na(delaiarret),-999,ifelse(delaiarret==0&durCig==0,0,ifelse(delaiarret==0,4,ifelse(delaiarret<=cutoff1,3,ifelse(delaiarret<=cutoff2,2,1))))))
      ceaseCigRaw<-round(with(ICareData,ifelse(is.na(delaiarret),-999,delaiarret)),2)
      outData<-paste(outData,ceaseCig)
      outRawX<-paste(outRawX,ceaseCigRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,5,sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'CessationTime')
   }
   
   if('PackYears'%in%covariateNames){
      cutoff1<-15
      cutoff2<-30
      cutoff3<-45
      cat(paste('Pack years:','Non',cutoff1,cutoff2,cutoff3,'\n'))
      packYears<-with(ICareData,ifelse(is.na(csi),-999,ifelse(is.na(packyear),-999,ifelse(csi==0,0,ifelse(packyear<cutoff1,1,ifelse(packyear<cutoff2,2,ifelse(packyear<cutoff3,3,4)))))))
      packYearsRaw<-round(with(ICareData,ifelse(is.na(csi),-999,ifelse(is.na(packyear),-999,packyear))),2)
      outData<-paste(outData,packYears)
      outRawX<-paste(outRawX,packYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,5,sep=ifelse(nchar(nCategories)>0,' ',''))
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
      outData<-paste(outData,darkPackYears)
      outRawX<-paste(outRawX,darkPackYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,5,sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'DarkPackYears')
   }
      
   if('NonFilteredPackYears'%in%covariateNames){
      cutoff1<-5
      cutoff2<-15
      cutoff3<-35
      cat(paste('NonFilteredPackYears:','Non',cutoff1,cutoff2,cutoff3,'\n'))
      nonFilteredPackYearsRaw<-round(ifelse(is.na(ICareData$csi),-999,ifelse(is.na(nonFilteredPackYears),-999,nonFilteredPackYears)),2)
      nonFilteredPackYears<-ifelse(is.na(ICareData$csi),-999,ifelse(is.na(nonFilteredPackYears),-999,ifelse(ICareData$csi==0,0,ifelse(nonFilteredPackYears<cutoff1,1,ifelse(nonFilteredPackYears<cutoff2,2,ifelse(nonFilteredPackYears<cutoff3,3,4))))))
      outData<-paste(outData,nonFilteredPackYears)
      outRawX<-paste(outRawX,nonFilteredPackYearsRaw,sep=ifelse(is.null(outRawX),'',' '))
      nCategories<-paste(nCategories,5,sep=ifelse(nchar(nCategories)>0,' ',''))
      covNames<-c(covNames,'NonFilteredPackYears')
   }
   
   # Write the output
   outFile<-file.path(filePath,outFile)
   outFileRawX<-gsub('.txt','_RawX.txt',outFile)
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

