analyseModels<-function(){
   require(DiPBaC)

   # Predictive runs
   inputFileNames<-c('paper_allMAll_main_packyears',
                     'paper_allMAll_main_packyears_darkpy_nonfilterpy',
                     'paper_allWAll_main_packyears',
                     'paper_AdMAll_main_packyears',
                     'paper_SqMAll_main_packyears',
                     'paper_SCMAll_main_packyears')
   runFileNames<-c('allMAllMainPY',
                   'allMAllMainPYDarkNonFilter',
                   'allWAllMainPY',
                   'AdMAllMainPY',
                   'SqMAllMainPY',
                   'SCMAllMainPY')
	upperLim<-length(runFileNames)
	for(i in 1:upperLim){
      cat(paste(inputFileNames[i],sep=' - '))
      cat('\n')
      runInfoObj<-readRunInfo('/home/dhastie/Imperial/Analyses/ICare/Data/Output',inputFileNames[i])
      disSimObj<-calcDissimilarityMatrix(runInfoObj)
      clusObj<-calcOptimalClustering(disSimObj)
      riskProfObj<-calcAvgRiskAndProfile(clusObj)
      predObj<-calcPredictions(riskProfObj,doRaoBlackwell=T,fullSweepPredictions=T,fullSweepLogOR=T)
      clusterOrder<-NULL
      clusterOrder<-plotRiskProfile(riskProfObj,outFile=paste('/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/',runFileNames[i],'_riskProf.png',sep=''),showRelativeRisk=T)
      save(list=c("runInfoObj","disSimObj","clusObj","riskProfObj","predObj","clusterOrder"),file=paste("~/Imperial/Analyses/ICare/RAnalysis/Paper/",runFileNames[i],".RData",sep=""))
      rm(list=c("runInfoObj","disSimObj","clusObj","riskProfObj","predObj","clusterOrder")) 
   }
}
  


