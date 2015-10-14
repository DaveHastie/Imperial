# Functions to create plots for Epidemiology paper
# 
# Author: David Hastie
###############################################################################

# Plot of log odds ratio for the clusters in the main+packyears analysis
createPlot1<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot1.eps'){
   require(ggplot2)
	
	load(dataFile)
	
	attach(riskProfObj)
   attach(riskProfClusObj)
   attach(clusObjRunInfoObj)
   	
   	
   postscript(outFile,width=3,height=3,paper="special",horizontal=F)
   	   	
   # Order by posterior median theta risk
   # Compute the means
   orderStat<-apply(risk,2,median)
   # Sort into ascending size
   meanSortIndex<-order(orderStat,decreasing=F) 
		
   # Reorder the risk matrix
   risk<-risk[,meanSortIndex]
	# Show the odds ratio
	logOddsRatio<-matrix(0,nrow(risk),ncol(risk))
	for(c in nClusters:1){
		logOddsRatio[,c]<-log(risk[,c]*(1-risk[,1])/((1-risk[,c])*risk[,1]))
	}
		
   clusterSizes<-clusterSizes[meanSortIndex]
   
	# Recompute the means and now also credible intervals
   logOddsRatioMeans<-apply(logOddsRatio,2,mean,trim=0.005)
   logOddsRatioMean<-sum(logOddsRatioMeans*clusterSizes)/sum(clusterSizes)
   logOddsRatioLower<-apply(logOddsRatio,2,quantile,0.05)
   logOddsRatioUpper<-apply(logOddsRatio,2,quantile,0.95)
   
	# The next line is to avoid outliers spoiling plot scales
   plotMax<-2*max(logOddsRatioUpper)-logOddsRatioMean
   		
   # Get the plot colors
   logOddsRatioColor<-ifelse(logOddsRatioLower>rep(logOddsRatioMean,nClusters),"high",
               ifelse(logOddsRatioUpper<rep(logOddsRatioMean,nClusters),"low","avg"))
   logOddsRatioDF<-data.frame("logOddsRatio"=c(),"cluster"=c(),"meanLogOddsRatio"=c(),
               "lowerLogOddsRatio"=c(),"upperLogOddsRatio"=c(),"fillColor"=c())
   
   # Restructure the data for plotting            
   for(c in 1:nClusters){
   	plotLogOddsRatio<-logOddsRatio[,c]
      plotLogOddsRatio<-plotLogOddsRatio[plotLogOddsRatio<plotMax]
      nPoints<-length(plotLogOddsRatio)
      logOddsRatioDF<-rbind(logOddsRatioDF,data.frame("logOddsRatio"=plotLogOddsRatio,"cluster"=rep(c,nPoints),
                        "meanLogOddsRatio"=rep(logOddsRatioMean,nPoints),
                        "lowerLogOddsRatio"=rep(logOddsRatioLower[c],nPoints),
                        "upperLogOddsRatio"=rep(logOddsRatioUpper[c],nPoints),                        
                        "fillColor"=rep(logOddsRatioColor[c],nPoints)))
  	}
   rownames(logOddsRatioDF)<-seq(1,nrow(logOddsRatioDF),1)
      	
   # Create the plot
   plotObj<-ggplot(logOddsRatioDF)
   plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=logOddsRatio,yintercept=meanLogOddsRatio),size=0.1)
   plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=logOddsRatio,fill=as.factor(fillColor)),size=0.1,outlier.size=0.25)
   plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerLogOddsRatio,colour=as.factor(fillColor)),size=1.0)
   plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperLogOddsRatio,colour=as.factor(fillColor)),size=1.0)
   plotObj<-plotObj+scale_fill_manual(values = c(high ="gray100",low ="gray30", avg ="gray65"))+
               scale_colour_manual(values = c(high ="gray100",low ="gray30", avg ="gray65"))+
               opts(legend.position="none")+labs(x="Cluster",y="Log odds ratio")
   plotObj<-plotObj+opts(axis.title.y=theme_text(size=10,angle=90),axis.title.x=theme_text(size=10))
   
   print(plotObj)
   	
	dev.off()
	detach(riskProfObj)
   detach(riskProfClusObj)
   detach(clusObjRunInfoObj)
   
}


# Plot of high risk cluster covariate probabilities
createPlot2<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot2.eps'){
	
	require(ggplot2)
	
	load(dataFile)
	
	attach(riskProfObj)
   attach(riskProfClusObj)
   attach(clusObjRunInfoObj)
   
   
   postscript(outFile,width=3,height=3,paper="special",horizontal=F)
   
   
   # Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 2, nrow = 2)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
   
   # Order by posterior median theta risk
   # Compute the means
   orderStat<-apply(risk,2,median)
   # Sort into ascending size
   meanSortIndex<-order(orderStat,decreasing=F) 
	highRiskInds<-meanSortIndex[(length(meanSortIndex)-2):length(meanSortIndex)]
	
	for(j in 1:nCovariates){
		profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
         	"lowerProb"=c(),"upperProb"=c())
   	clusterNames<-c("high","higher","highest")
		for(c in 1:3){
   		probMat<-profile[,highRiskInds[c],j,]
   		probMean<-apply(probMat,2,mean)
      	probLower<-apply(probMat,2,quantile,0.05)
      	probUpper<-apply(probMat,2,quantile,0.95)
      
		
      	for(k in 1:nCategories[j]){
         	profileDF<-rbind(profileDF,data.frame("prob"=probMean[k],"cluster"=clusterNames[c],"category"=k-1,
                     "lowerProb"=probLower[k],
                     "upperProb"=probUpper[k]))
            rownames(profileDF)<-seq(1,nrow(profileDF),1)
         }
      }
		
   	
   	plotObj<-ggplot(profileDF)
		plotObj<-plotObj+geom_pointrange(aes(x=as.factor(category),y=prob,ymin=lowerProb,ymax=upperProb,shape=cluster),size=0.15,position=position_dodge(width=0.5))
		plotObj<-plotObj+scale_y_continuous(limits=c(0,1))
      plotObj<-plotObj+opts(legend.position="none")+labs(x="Category")+opts(axis.title.x=theme_text(size=7))+
				opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
		
   	if(j==1||j==3){
      	plotObj<-plotObj+labs(y="Probability")+opts(axis.title.y=theme_text(size=7,angle=90))
   	}else{
      	plotObj<-plotObj+opts(axis.title.y=theme_blank())
   	}
   	plotObj<-plotObj+opts(title=covNames[j],plot.title=theme_text(size=9))
   	plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(j==2||j==4,0.25,0),0.25,ifelse(j==1||j==3,0.25,0.1)),'lines'))+
         	opts(print.margin=unit(c(0,0,0,0),'lines'))
   	
   	print(plotObj,vp=viewport(layout.pos.row=ifelse(j<3,1,2),layout.pos.col=ifelse(j==1||j==3,1,2)))
   }
	
	
	dev.off()
	detach(riskProfObj)
   detach(riskProfClusObj)
   detach(clusObjRunInfoObj)
   
}


# Plot of predicted log odds ratio densities - int / dur combinations, missing tsc
createPlot3<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot3.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFile)
	
	attach(predObj)
   
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-matrix(70:85,4,4,byrow=T)
   # Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 4, nrow = 4)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
	
	maxX<-0
	maxY<-0
	denObj<-vector(mode="list")
	for(i in 70:85){
		denObj[[i-69]]<-density(na.omit(logORPerSweep[,i]))
		maxX<-max(maxX,denObj[[i-69]]$x)
		maxY<-max(maxY,denObj[[i-69]]$y)
	}
	
	for(k2 in 4:1){
		for(k1 in 1:4){
			plotDF<-data.frame('logOddsRatio'=denObj[[relScenarios[k1,k2]-69]]$x,'density'=denObj[[relScenarios[k1,k2]-69]]$y)
			plotObj<-ggplot(plotDF)
			plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density),size=0.2)
			plotObj<-plotObj+opts(legend.position="none")
			
			if(k2==1){
				plotObj<-plotObj+labs(x=paste('Log Odds Ratio (Intensity ',k1,')',sep=''))+opts(axis.title.x=theme_text(size=7))						
			}else{
				plotObj<-plotObj+opts(axis.title.x=theme_blank())
			}
			if(k1==1){
				plotObj<-plotObj+labs(y=paste('Density (Duration ',k2,')',sep=''))+opts(axis.title.y=theme_text(size=7,angle=90))						
			}else{
				plotObj<-plotObj+opts(axis.title.y=theme_blank())
			}
			plotObj<-plotObj+scale_x_continuous(limits=c(0,maxX))+scale_y_continuous(limits=c(0,maxY))
			plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   		plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(k1==4,0.25,0),ifelse(k2==1,0.25,0),0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=5-k2,layout.pos.col=k1))
		
		}
	}
	
	dev.off()
	detach(predObj)
	
}

# Plot of predicted log odds ratio densities - int / dur /tsc combinations
createPlot4<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot4.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFile)
	
	attach(predObj)
   
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-array(6:69,dim=c(4,4,4))
   # Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 4, nrow = 4)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
	
	maxX<-0
	maxY<-0
	denObj<-vector(mode="list")
	for(i in 6:69){
		denObj[[i-5]]<-density(na.omit(logORPerSweep[,i]))
		maxX<-max(maxX,denObj[[i-5]]$x)
		maxY<-max(maxY,denObj[[i-5]]$y)
	}
	
	for(k2 in 4:1){
		for(k1 in 1:4){
			plotDF<-NULL
			for(k3 in 1:4){
				plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObj[[relScenarios[k3,k2,k1]-5]]$x,'density'=denObj[[relScenarios[k3,k2,k1]-5]]$y,'tsc'=rep(k3,length(denObj[[relScenarios[k3,k2,k1]-5]]$y))))
			}
			plotObj<-ggplot(plotDF)
			plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density,group=tsc,linetype=ifelse(tsc==1,2,ifelse(tsc==2,3,ifelse(tsc==3,1,6)))),size=0.2)
			plotObj<-plotObj+opts(legend.position="none")
			
			if(k2==1){
				plotObj<-plotObj+labs(x=paste('Log Odds Ratio (Intensity ',k1,')',sep=''))+opts(axis.title.x=theme_text(size=7))						
			}else{
				plotObj<-plotObj+opts(axis.title.x=theme_blank())
			}
			if(k1==1){
				plotObj<-plotObj+labs(y=paste('Density (Duration ',k2,')',sep=''))+opts(axis.title.y=theme_text(size=7,angle=90))						
			}else{
				plotObj<-plotObj+opts(axis.title.y=theme_blank())
			}
			plotObj<-plotObj+scale_x_continuous(limits=c(0,maxX))+scale_y_continuous(limits=c(0,maxY))
			plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   		plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(k1==4,0.25,0),ifelse(k2==1,0.25,0),0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=5-k2,layout.pos.col=k1))
			
		}
	}
	
	dev.off()
	detach(predObj)
	
}

# Plot of predicted log odds ratio densities - 2 fixed csi scenarios
createPlot5<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot5.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFile)
	
	attach(predObj)
   
	postscript(outFile,width=3,height=3,paper="special",horizontal=F)
   
   # Set up the layout for the plot
	
	maxX<-0
	maxY<-0
	denObj<-vector(mode="list")
	for(i in 86:87){
		denObj[[i-85]]<-density(na.omit(logORPerSweep[,i]))
		maxX<-max(maxX,denObj[[i-85]]$x)
		maxY<-max(maxY,denObj[[i-85]]$y)
	}
	
	plotDF<-data.frame('logOddsRatio'=denObj[[1]]$x,'density'=denObj[[1]]$y,'scenario'=rep(1,length(denObj[[1]]$y)))
	plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObj[[2]]$x,'density'=denObj[[2]]$y,'scenario'=rep(2,length(denObj[[2]]$y))))
	plotObj<-ggplot(plotDF)
	plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density,group=scenario,linetype=scenario),size=0.2)
	plotObj<-plotObj+opts(legend.position="none")
	plotObj<-plotObj+labs(x='Log Odds Ratio')+opts(axis.title.x=theme_text(size=7))						
	plotObj<-plotObj+labs(y='Density')+opts(axis.title.y=theme_text(size=7,angle=90))						
	plotObj<-plotObj+scale_x_continuous(limits=c(0,maxX))+scale_y_continuous(limits=c(0,maxY))
	plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
	print(plotObj)
				
	dev.off()
	detach(predObj)
	
}

# Plot of predicted log odds ratio densities - dark / non-filtered combinations, different smokers
createPlot6<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPYDarkNonFilter.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot6.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFile)
	
	attach(predObj)
   
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-array(86:133,dim=c(3,4,4))
	relScenarios[1,,]<-matrix(86:101,nrow=4,byrow=T)
	relScenarios[2,,]<-matrix(102:117,nrow=4,byrow=T)
   relScenarios[3,,]<-matrix(118:133,nrow=4,byrow=T)
	# Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 4, nrow = 4)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
	
	maxX<-0
	maxY<-0
	denObj<-vector(mode="list")
	for(i in 86:133){
		denObj[[i-85]]<-density(na.omit(logORPerSweep[,i]))
		maxX<-max(maxX,denObj[[i-85]]$x)
		maxY<-max(maxY,denObj[[i-85]]$y)
	}
	
	for(k2 in 4:1){
		for(k1 in 1:4){
         if(k1<4&&k2<4){
			   plotDF<-data.frame('logOddsRatio'=denObj[[relScenarios[1,k1,k2]-85]]$x,'density'=denObj[[relScenarios[1,k1,k2]-85]]$y,'smoker'=rep(1,length(denObj[[relScenarios[1,k1,k2]-85]]$x)))
			   plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObj[[relScenarios[2,k1,k2]-85]]$x,'density'=denObj[[relScenarios[2,k1,k2]-85]]$y,'smoker'=rep(2,length(denObj[[relScenarios[2,k1,k2]-85]]$x))))
			   plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObj[[relScenarios[3,k1,k2]-85]]$x,'density'=denObj[[relScenarios[3,k1,k2]-85]]$y,'smoker'=rep(3,length(denObj[[relScenarios[3,k1,k2]-85]]$x))))
			}else{
            plotDF<-data.frame('logOddsRatio'=denObj[[relScenarios[2,k1,k2]-85]]$x,'density'=denObj[[relScenarios[2,k1,k2]-85]]$y,'smoker'=rep(2,length(denObj[[relScenarios[2,k1,k2]-85]]$x)))
            plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObj[[relScenarios[3,k1,k2]-85]]$x,'density'=denObj[[relScenarios[3,k1,k2]-85]]$y,'smoker'=rep(3,length(denObj[[relScenarios[3,k1,k2]-85]]$x))))            
         }
         plotObj<-ggplot(plotDF)
			plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density,group=smoker,linetype=smoker),size=0.2)
			plotObj<-plotObj+opts(legend.position="none")
			
			if(k2==1){
				plotObj<-plotObj+labs(x=paste('Log Odds Ratio (Dark PY ',k1,')',sep=''))+opts(axis.title.x=theme_text(size=7))						
			}else{
				plotObj<-plotObj+opts(axis.title.x=theme_blank())
			}
			if(k1==1){
				plotObj<-plotObj+labs(y=paste('Density (Non-filtered PY ',k2,')',sep=''))+opts(axis.title.y=theme_text(size=7,angle=90))						
			}else{
				plotObj<-plotObj+opts(axis.title.y=theme_blank())
			}
			plotObj<-plotObj+scale_x_continuous(limits=c(0,maxX))+scale_y_continuous(limits=c(0,maxY))
			plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   		plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(k1==4,0.25,0),ifelse(k2==1,0.25,0),0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=5-k2,layout.pos.col=k1))
			
		}
	}
	
	dev.off()
	detach(predObj)
	
}

# Plot of predicted log odds ratio densities - int / dur combinations, missing tsc, different histologies
createPlot7<-function(dataFiles=c('/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
				'/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/AdMAllMainPY.RData',
				'/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/SqMAllMainPY.RData',
				'/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/SCMAllMainPY.RData'),
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/plot7.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFiles[1])
	predObjAll<-predObj
	load(dataFiles[2])
	predObjAd<-predObj
	load(dataFiles[3])
	predObjSq<-predObj
	load(dataFiles[4])
	predObjSC<-predObj
	
		
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-matrix(70:85,4,4,byrow=T)
   # Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 4, nrow = 4)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
	
	maxX<-0
	maxY<-0
	denObjAll<-vector(mode="list")
	denObjAd<-vector(mode="list")
	denObjSq<-vector(mode="list")
	denObjSC<-vector(mode="list")
	for(i in 70:85){
		denObjAll[[i-69]]<-density(na.omit(predObjAll$logORPerSweep[,i]))
		denObjAd[[i-69]]<-density(na.omit(predObjAd$logORPerSweep[,i]))
		denObjSq[[i-69]]<-density(na.omit(predObjSq$logORPerSweep[,i]))
		denObjSC[[i-69]]<-density(na.omit(predObjSC$logORPerSweep[,i]))
		maxX<-max(maxX,denObjAll[[i-69]]$x,denObjAd[[i-69]]$x,denObjSq[[i-69]]$x,denObjSC[[i-69]]$x)
		maxY<-max(maxY,denObjAll[[i-69]]$y,denObjAd[[i-69]]$y,denObjSq[[i-69]]$y,denObjSC[[i-69]]$y)
	}
	
	for(k2 in 4:1){
		for(k1 in 1:4){
			plotDF<-data.frame('logOddsRatio'=denObjAll[[relScenarios[k1,k2]-69]]$x,'density'=denObjAll[[relScenarios[k1,k2]-69]]$y,'histology'=rep(1,length(denObjAll[[relScenarios[k1,k2]-69]]$y)))
			plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObjAd[[relScenarios[k1,k2]-69]]$x,'density'=denObjAd[[relScenarios[k1,k2]-69]]$y,'histology'=rep(2,length(denObjAd[[relScenarios[k1,k2]-69]]$y))))
			plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObjSq[[relScenarios[k1,k2]-69]]$x,'density'=denObjSq[[relScenarios[k1,k2]-69]]$y,'histology'=rep(3,length(denObjSq[[relScenarios[k1,k2]-69]]$y))))
			plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObjSC[[relScenarios[k1,k2]-69]]$x,'density'=denObjSC[[relScenarios[k1,k2]-69]]$y,'histology'=rep(4,length(denObjSC[[relScenarios[k1,k2]-69]]$y))))
			
			plotObj<-ggplot(plotDF)
			plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density,group=histology,linetype=ifelse(histology==1,1,ifelse(histology==2,6,ifelse(histology==3,3,2)))),size=0.2)
			plotObj<-plotObj+opts(legend.position="none")
			
			if(k2==1){
				plotObj<-plotObj+labs(x=paste('Log Odds Ratio (Intensity ',k1,')',sep=''))+opts(axis.title.x=theme_text(size=7))						
			}else{
				plotObj<-plotObj+opts(axis.title.x=theme_blank())
			}
			if(k1==1){
				plotObj<-plotObj+labs(y=paste('Density (Duration ',k2,')',sep=''))+opts(axis.title.y=theme_text(size=7,angle=90))						
			}else{
				plotObj<-plotObj+opts(axis.title.y=theme_blank())
			}
			plotObj<-plotObj+scale_x_continuous(limits=c(-1,maxX))+scale_y_continuous(limits=c(0,maxY))
			plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   		plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(k1==4,0.25,0),ifelse(k2==1,0.25,0),0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=5-k2,layout.pos.col=k1))
			
		}
	}
	
	dev.off()
			
}

# Create a table of the features of the prediction scenarios
createTable1<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
      packYearsFile='/home/dhastie/Imperial/Analyses/ICare/Data/Output/paper_allMAll_main_packyears_predictPackYears.txt',
      outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/table1.tex'){
   
   require(DiPBaC)
   load(dataFile)
   
   attach(predObj)
   
   relScenarios<-matrix(70:85,4,4,byrow=T)
      
   logORMeans<-apply(logORPerSweep,2,mean)   
   logORStdDev<-apply(logORPerSweep,2,sd)

	nSweeps<-runInfoObj$nSweeps
	nBurn<-runInfoObj$nBurn
	nFilter<-runInfoObj$nFilter
	# Restrict to sweeps after burn in
   firstLine<-2+nBurn/nFilter
   lastLine<-1+(nSweeps+nBurn)/nFilter
   packYearsData<-readLines(packYearsFile)
	packYearsData<-packYearsData[firstLine:lastLine]
	packYearsMat<-matrix(as.numeric(unlist(strsplit(packYearsData," "))),nrow=nSweeps/nFilter,byrow=T)
	
	packYearsMeans<-apply(packYearsMat,2,mean)
	packYearsStdDev<-apply(packYearsMat,2,sd)
	
	outText<-c("\\begin{table}[hp]","\t\\begin{center}","\t\\begin{tabular}{cc|cc}",
			"\t\tIntensity & Duration & Log Odds Ratio & Pack Years \\\\",
			"\t\t\\hline")
	for(j in 1:4){
		for(k in 1:4){
			tmpText<-paste(ifelse(k==1,paste("\t\t\\multirow{4}{*}{",j,"}",sep=""),""),k,sep=" & ")
			tmpText<-paste(tmpText,"&",format(round(logORMeans[relScenarios[j,k]],2),nsmall=2),paste("(",format(round(logORStdDev[relScenarios[j,k]],2),nsmall=2),")",sep=""))
			tmpText<-paste(tmpText,"&",format(round(packYearsMeans[relScenarios[j,k]],2),nsmall=2),paste("(",format(round(packYearsStdDev[relScenarios[j,k]],2),nsmall=2),")",sep=""))
			tmpText<-paste(tmpText,"\\\\")
			outText<-c(outText,tmpText)			
		}
		outText<-c(outText,"\t\t\\hline")
	}	
	
   outText<-c(outText,"\t\\end{tabular}")
	outText<-c(outText,"\t\\end{center}")	
	outText<-c(outText,"\\caption{Table posterior distribution means (standard deviations) for 16 combinations of intensity and duration, where time since cessation and pack years missing, corresponding to the distributions in Figure 3. Also shown for each combination is the mean (standard deviation) of expected the pack years category.}")
	outText<-c(outText,"\\end{table}")
	
	write(outText,file=outFile,append=F)
	detach(predObj)
	   
}

# Create a plot of the clustering
createEPlot1<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/eplot1.eps'){
   
   require(ggplot2)
	require(DiPBaC)
	
	load(dataFile)
   attach(clusObj)
   attach(clusObjRunInfoObj)
	
   postscript(outFile,width=6,height=6,paper="special",horizontal=F)      
   
   # Read in the raw X data
   rawXFileName<-gsub('.txt','_RawX.txt',inputFileName)
   rawXData<-readLines(rawXFileName)
   rawXData<-rawXData[(nCovariates+3):length(rawXData)]
   rawXData<-matrix(as.numeric(unlist(strsplit(rawXData," "))),ncol=nCovariates,byrow=T)
   relInd<-match('CessationTime',covNames)
   meanCat1<-mean(rawXData[xMat[,relInd]==1,relInd])
   rawXData[xMat[,relInd]==0,relInd]<-meanCat1+10
	rawXData[rawXData==-999]<-NA
   
	includeVec<-rep(T,nSubjects)
   for(i in 1:nSubjects){
      if(any(is.na(rawXData[i,]))){
         includeVec[i]<-F
      }
   }
   rawXData<-rawXData[includeVec,]
   clustering<-clustering[includeVec]

   d<-dist(rawXData)
   
   principalComp<-cmdscale(d,k=2)
	plotDF<-data.frame('PC1'=principalComp[,1],'PC2'=principalComp[,2],'cluster'=match(clustering,clusterOrder))
	plotObj<-ggplot(plotDF)
	plotObj<-plotObj+geom_point(aes(x=PC1,y=PC2,shape=factor(cluster),colour=factor(cluster)))
	plotObj<-plotObj+geom_point(data=subset(plotDF,cluster==1),aes(x=PC1,y=PC2,shape=factor(cluster),colour=factor(cluster)))	
	plotObj<-plotObj+labs(x='Principal component 1')+opts(axis.title.x=theme_text(size=7))						
	plotObj<-plotObj+labs(y='Principal component 2')+opts(axis.title.y=theme_text(size=7,angle=90))						
	plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   plotObj<-plotObj+scale_colour_manual(name='Cluster',values=c("1"="black","2"="white","3"="green","4"="pink","5"="blue","6"="yellow","7"="red","8"="purple"))
   plotObj<-plotObj+scale_shape_manual(name='Cluster',values=c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5,"6"=6,"7"=1,"8"=2))
	print(plotObj)

   dev.off()
   
   detach(clusObjRunInfoObj)
   detach(clusObj)
	
}

createEPlot2<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/eplot2.eps'){

	require(ggplot2)
	require(DiPBaC)
	
	load(dataFile)
	
	attach(predObj)
   
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-matrix(70:85,4,4,byrow=T)
   # Set up the layout for the plot
	
	maxX<-0
	maxY<-0
	meanObj<-rep(0,16)
	lQuantileObj<-rep(0,16)
	uQuantileObj<-rep(0,16)
	for(i in 70:85){
		meanObj[i-69]<-mean(na.omit(logORPerSweep[,i]))
		lQuantileObj[i-69]<-quantile(na.omit(logORPerSweep[,i]),0.05)
		uQuantileObj[i-69]<-quantile(na.omit(logORPerSweep[,i]),0.95)
		maxY<-max(maxY,uQuantileObj[i-69]+0.5)
	}
	
	plotDF<-NULL
	for(k2 in 4:1){
		for(k1 in 1:4){
			plotDF<-rbind(plotDF,data.frame('logOddsRatio'=meanObj[relScenarios[k1,k2]-69],'intensity'=k1,'duration'=k2,
					'lowLOR'=lQuantileObj[relScenarios[k1,k2]-69],'highLOR'=uQuantileObj[relScenarios[k1,k2]-69]))
		}
	}
	
	plotObj<-ggplot(plotDF)
	plotObj<-plotObj+geom_line(aes(x=as.factor(intensity),y=logOddsRatio,group=as.factor(duration),colour=as.factor(duration)),size=0.5)
	plotObj<-plotObj+geom_pointrange(aes(x=as.factor(intensity),y=logOddsRatio,ymin=lowLOR,ymax=highLOR,colour=as.factor(duration)),size=0.75)
	plotObj<-plotObj+geom_errorbar(aes(x=as.factor(intensity),ymin=lowLOR,ymax=highLOR,colour=as.factor(duration)),width=0.2)
	plotObj<-plotObj+opts(legend.position="none")
	plotObj<-plotObj+labs(x='Intensity')+opts(axis.title.x=theme_text(size=7))						
	plotObj<-plotObj+labs(y='Log Odds Ratio')+opts(axis.title.y=theme_text(size=7,angle=90))						
	plotObj<-plotObj+scale_y_continuous(limits=c(0,maxY))
	plotObj<-plotObj+scale_colour_manual(name='Duration',values=c("1"="black","2"="green","3"="blue","4"="red"))
	plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   plotObj<-plotObj+opts(plot.margin=unit(c(0.25,0.25,0.25,0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
	print(plotObj)
			
	
	dev.off()
	detach(predObj)
	
}

	
	



# Plot of predicted log odds ratio densities - int / dur combinations, missing tsc, men vs women
createExtraPlot1<-function(dataFiles=c('/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
				'/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allWAllMainPY.RData'),
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/extraplot1.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFiles[1])
	predObjMen<-predObj
	load(dataFiles[2])
	predObjWomen<-predObj
	
	
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-matrix(70:85,4,4,byrow=T)
   # Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 4, nrow = 4)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
	
	maxX<-0
	maxY<-0
	denObjMen<-vector(mode="list")
	denObjWomen<-vector(mode="list")

	for(i in 70:85){
		denObjMen[[i-69]]<-density(na.omit(predObjMen$logORPerSweep[,i]))
		denObjWomen[[i-69]]<-density(na.omit(predObjWomen$logORPerSweep[,i]))
		maxX<-max(maxX,denObjMen[[i-69]]$x,denObjWomen[[i-69]]$x)
		maxY<-max(maxY,denObjMen[[i-69]]$y,denObjWomen[[i-69]]$y)
	}
	
	for(k2 in 4:1){
		for(k1 in 1:4){
			plotDF<-data.frame('logOddsRatio'=denObjMen[[relScenarios[k1,k2]-69]]$x,'density'=denObjMen[[relScenarios[k1,k2]-69]]$y,'sex'=rep(1,length(denObjMen[[relScenarios[k1,k2]-69]]$y)))
			plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObjWomen[[relScenarios[k1,k2]-69]]$x,'density'=denObjWomen[[relScenarios[k1,k2]-69]]$y,'sex'=rep(2,length(denObjWomen[[relScenarios[k1,k2]-69]]$y))))
			
			plotObj<-ggplot(plotDF)
			plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density,group=sex,linetype=sex),size=0.2)
			plotObj<-plotObj+opts(legend.position="none")
			
			if(k2==1){
				plotObj<-plotObj+labs(x=paste('Log Odds Ratio (Intensity ',k1,')',sep=''))+opts(axis.title.x=theme_text(size=7))						
			}else{
				plotObj<-plotObj+opts(axis.title.x=theme_blank())
			}
			if(k1==1){
				plotObj<-plotObj+labs(y=paste('Density (Duration ',k2,')',sep=''))+opts(axis.title.y=theme_text(size=7,angle=90))						
			}else{
				plotObj<-plotObj+opts(axis.title.y=theme_blank())
			}
			plotObj<-plotObj+scale_x_continuous(limits=c(-1,7))+scale_y_continuous(limits=c(0,maxY))
			plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   		plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(k1==4,0.25,0),ifelse(k2==1,0.25,0),0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=5-k2,layout.pos.col=k1))
			
		}
	}
	
	dev.off()
	
}

# Plot of predicted log odds ratio densities - int / dur combinations, missing tsc, men vs women
createExtraPlot2<-function(dataFiles=c('/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
				'/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/EAGLEMAllMainPY.RData'),
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/extraplot2.eps'){
	
	require(ggplot2)
	require(DiPBaC)
	
	load(dataFiles[1])
	predObjICare<-predObj
	load(dataFiles[2])
	predObjEAGLE<-predObj
	
	
	postscript(outFile,width=6,height=6,paper="special",horizontal=F)
   
	# Relevant scenarios
	relScenarios<-matrix(70:85,4,4,byrow=T)
   # Set up the layout for the plot
   plotLayout<-grid.layout(ncol = 4, nrow = 4)
   grid.newpage()
   pushViewport(viewport(layout = plotLayout))   
	
	maxX<-0
	maxY<-0
	denObjICare<-vector(mode="list")
	denObjEAGLE<-vector(mode="list")
	
	for(i in 70:85){
		denObjICare[[i-69]]<-density(na.omit(predObjICare$logORPerSweep[,i]))
		denObjEAGLE[[i-69]]<-density(na.omit(predObjEAGLE$logORPerSweep[,i]))
		maxX<-max(maxX,denObjICare[[i-69]]$x,denObjEAGLE[[i-69]]$x)
		maxY<-max(maxY,denObjICare[[i-69]]$y,denObjEAGLE[[i-69]]$y)
	}
	
	for(k2 in 4:1){
		for(k1 in 1:4){
			plotDF<-data.frame('logOddsRatio'=denObjICare[[relScenarios[k1,k2]-69]]$x,'density'=denObjICare[[relScenarios[k1,k2]-69]]$y,'study'=rep(1,length(denObjICare[[relScenarios[k1,k2]-69]]$y)))
			plotDF<-rbind(plotDF,data.frame('logOddsRatio'=denObjEAGLE[[relScenarios[k1,k2]-69]]$x,'density'=denObjEAGLE[[relScenarios[k1,k2]-69]]$y,'study'=rep(2,length(denObjEAGLE[[relScenarios[k1,k2]-69]]$y))))
			
			plotObj<-ggplot(plotDF)
			plotObj<-plotObj+geom_line(aes(x=logOddsRatio,y=density,group=study,linetype=study),size=0.2)
			plotObj<-plotObj+opts(legend.position="none")
			
			if(k2==1){
				plotObj<-plotObj+labs(x=paste('Log Odds Ratio (Intensity ',k1,')',sep=''))+opts(axis.title.x=theme_text(size=7))						
			}else{
				plotObj<-plotObj+opts(axis.title.x=theme_blank())
			}
			if(k1==1){
				plotObj<-plotObj+labs(y=paste('Density (Duration ',k2,')',sep=''))+opts(axis.title.y=theme_text(size=7,angle=90))						
			}else{
				plotObj<-plotObj+opts(axis.title.y=theme_blank())
			}
			plotObj<-plotObj+scale_x_continuous(limits=c(-1,maxX))+scale_y_continuous(limits=c(0,maxY))
			plotObj<-plotObj+opts(axis.text.x=theme_text(size=7))+opts(axis.text.y=theme_text(size=7))
   		plotObj<-plotObj+opts(plot.margin=unit(c(0.25,ifelse(k1==4,0.25,0),ifelse(k2==1,0.25,0),0.25),'lines'))+
         		opts(print.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=5-k2,layout.pos.col=k1))
			
		}
	}
	
	dev.off()
	
}


createETable2<-function(dataFile='/home/dhastie/Imperial/Analyses/ICare/RAnalysis/Paper/allMAllMainPY.RData',
		CSIFile='/home/dhastie/Imperial/Analyses/ICare/Data/Input/allMAll_main_packyears_CSI.txt',
		outFile='/home/dhastie/Imperial/Analyses/ICare/Plots/Paper/etable2.tex'){
	
	load(dataFile)
		
	attach(riskProfObj)
	attach(riskProfClusObj)
   attach(clusObjRunInfoObj)
   	   	   	
   # Order by posterior median theta risk
   # Compute the means
   orderStat<-apply(risk,2,median)
	# Sort into ascending size
   meanSortIndex<-order(orderStat,decreasing=F) 
	
	risk<-risk[,meanSortIndex]
	profile<-profile[,meanSortIndex,,]	
	
	clusterSizes<-clusterSizes[meanSortIndex]

	# Show the odds ratio
	logOddsRatio<-matrix(0,nrow(risk),ncol(risk))
	for(c in nClusters:1){
		logOddsRatio[,c]<-log(risk[,c]*(1-risk[,1])/((1-risk[,c])*risk[,1]))
	}
	
	# Compute the CSI
	CSIVec<-scan(CSIFile,quiet=T)
	CSIVec[CSIVec<0]<-NA
	CSIMeans<-rep(0,nClusters)
	CSIStdDevs<-rep(0,nClusters)
	for(c in 1:nClusters){
		CSIMeans[c]<-mean(na.omit(CSIVec[clustering==c]))
		CSIStdDevs[c]<-ifelse(any(na.omit(CSIVec[clustering==c]>0)),sd(na.omit(CSIVec[clustering==c])),0)
	}
	CSIMeans<-CSIMeans[meanSortIndex]
	CSIStdDevs<-CSIStdDevs[meanSortIndex]
	
	outText<-c("\\begin{landscape}","\\begin{table}[hp]","\t\\begin{center}",paste("\t\\begin{tabular}{ll|",paste(rep("c",nClusters),collapse=""),"}",sep=""),
				paste("\t\t&\\multicolumn{",nClusters,"}","{c}{Cluster} \\\\",sep=""),
				paste("\t\t& &",paste(1:nClusters,collapse=" & "),"\\\\"),
				"\t\t\\hline")
	
	tmpText<-paste("\t\t\\multicolumn{2}{l|}{No. Subjects} & ",paste(clusterSizes,collapse=" & ")," \\\\",sep="")
	outText<-c(outText,tmpText,"\t\t\\hline")
	tmpText<-"\t\t\\multicolumn{2}{l|}{Log OR}"
	for(c in 1:nClusters){
		tmpText<-paste(tmpText,"&",round(mean(na.omit(logOddsRatio[,c])),2),ifelse(c==1,"",paste("(",round(sd(na.omit(logOddsRatio[,c])),2),")",sep="")))   			
	}
	tmpText<-paste(tmpText,"\\\\")
	outText<-c(outText,tmpText,"\t\t\\hline")
	shortCovs<-c('INT','DUR','TSC','PY')
	for(j in 1:nCovariates){
		for(k in 1:nCategories[j]){
			tmpText<-paste(ifelse(k==1,paste("\t\t\\multirow{",nCategories[j],"}{*}{",shortCovs[j],"} &",sep="")," &"),k-1)
			for(c in 1:nClusters){
				tmpText<-paste(tmpText,"&",format(round(mean(na.omit(profile[,c,j,k])),2),nsmall=2),paste("(",format(round(sd(na.omit(profile[,c,j,k])),2),nsmall=2),")",sep=""))
   		}
			tmpText<-paste(tmpText,"\\\\")
			outText<-c(outText,tmpText)
		}
		outText<-c(outText,"\t\t\\hline")
	}
	tmpText<-"CSI &  "
	for(c in 1:nClusters){
		tmpText<-paste(tmpText,"&",format(round(CSIMeans[c],2),nsmall=2),paste("(",format(round(CSIStdDevs[c],2),nsmall=2),")",sep=""))
   }
	tmpText<-paste(tmpText,"\\\\")
	outText<-c(outText,tmpText)
	
	outText<-c(outText,"\t\\end{tabular}")
	outText<-c(outText,"\t\\end{center}")	
	outText<-c(outText,"\\caption{Table of distribution mean (standard deviation) for characteristics of clusters from the representative clustering of the analysis of Intensity, Duration, Time Since Cessation and PackYears. For the covariates, the distribution is of the probability that the covariate is in each category.}")
	outText<-c(outText,"\\end{table}")
	outText<-c(outText,"\\end{landscape}")
	
	write(outText,file=outFile,append=F)
	
	detach(riskProfObj)
	detach(riskProfClusObj)
   detach(clusObjRunInfoObj)
   
}
