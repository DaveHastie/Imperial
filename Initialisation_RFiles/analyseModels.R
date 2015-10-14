analyseModels<-function(){
   require(DiPBaC)

   # Predictive runs
   inputFileNames<-c('real_discrete_bernoulli_1000_4_01',
			'real_discrete_bernoulli_1000_4_10',
			'real_discrete_bernoulli_1000_4_20',
			'real_discrete_bernoulli_2000_4_01',
			'real_discrete_bernoulli_2000_4_10',
			'real_discrete_bernoulli_2000_4_20',
			'real_discrete_bernoulli_3000_4_01',
			'real_discrete_bernoulli_3000_4_10',
			'real_discrete_bernoulli_3000_4_20',
			'real_discrete_bernoulli_4000_4_01',
			'real_discrete_bernoulli_4000_4_10',
			'real_discrete_bernoulli_4000_4_20',
			'real_discrete_bernoulli_4658_4_01',
			'real_discrete_bernoulli_4658_4_05',
			'real_discrete_bernoulli_4658_4_10',
			'real_discrete_bernoulli_4658_4_15',
			'real_discrete_bernoulli_4658_4_20',
			'real_discrete_bernoulli_4658_6_01',
			'real_discrete_bernoulli_4658_6_05',
			'real_discrete_bernoulli_4658_6_10',
			'real_discrete_bernoulli_4658_6_15',
			'real_discrete_bernoulli_4658_6_20')
			
   runFileNames<-c('realDiscBerN1000J4C1','realDiscBerN1000J4C10',
						'realDiscBerN1000J4C20',
						'realDiscBerN2000J4C1','realDiscBerN2000J4C10',
						'realDiscBerN2000J4C20',
						'realDiscBerN3000J4C1','realDiscBerN3000J4C10',
						'realDiscBerN3000J4C20',
						'realDiscBerN4000J4C1','realDiscBerN4000J4C10',
						'realDiscBerN4000J4C20',
						'realDiscBerN4658J4C1','realDiscBerN4658J4C5',
						'realDiscBerN4658J4C10','realDiscBerN4658J4C15',
						'realDiscBerN4658J4C20',
                  'realDiscBerN4658J6C1','realDiscBerN4658J6C5',
						'realDiscBerN4658J6C10','realDiscBerN4658J6C15',
						'realDiscBerN4658J6C20')

				upperLim<-length(runFileNames)
	for(i in 1:upperLim){
      cat(paste(inputFileNames[i],sep=' - '))
      cat('\n')
      runInfoObj<-readRunInfo('/home/dhastie/Imperial/Analyses/Initialisation/Data/Output',inputFileNames[i])
      disSimObj<-calcDissimilarityMatrix(runInfoObj)
      clusObj<-calcOptimalClustering(disSimObj)
      riskProfObj<-calcAvgRiskAndProfile(clusObj)
      predObj<-calcPredictions(riskProfObj,doRaoBlackwell=T,fullSweepPredictions=T,fullSweepLogOR=T)
      clusterOrder<-NULL
      clusterOrder<-plotRiskProfile(riskProfObj,outFile=paste('/home/dhastie/Imperial/Analyses/Initialisation/Plots/',runFileNames[i],'_riskProf.png',sep=''),showRelativeRisk=T)
      save(list=c("runInfoObj","disSimObj","clusObj","riskProfObj","predObj","clusterOrder"),file=paste("~/Imperial/Analyses/Initialisation/RAnalysis/",runFileNames[i],".RData",sep=""))
      rm(list=c("runInfoObj","disSimObj","clusObj","riskProfObj","predObj","clusterOrder")) 
   }
}
  


