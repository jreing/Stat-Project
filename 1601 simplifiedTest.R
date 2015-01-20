
{
iteration <- setClass ("iteration",
		slots=c(xbars="numeric",
				pvals="numeric",
				refVec="logical",
				numOfTrues="numeric",
				truesMu1="numeric",
				falsesMu1="numeric",
				size="numeric",
				proceduresApplied= "list"
		)
)
procedure <- setClass ("procedure",
		slots=c(procedureName="function",
				rejects="list",
				fr="numeric",
				power="numeric",
				fdr="numeric",
				fwer="numeric"
		)
)


setRefVector <- function (size, numOfTrues){
	#initialize reference vector
	refvec<-rep (NA,size)
	refvec[1:numOfTrues]<- TRUE # these will be true discoveries
	refvec[(numOfTrues+1):size]<- FALSE # these will be False Discoveries
	return (refvec)
}

getRandomXBars <-function (size, numOfTrues, truesMu1=0, falsesMu1, sd=1){
	xbar<- rep(NA,size) #init xbar
	xbar[1:numOfTrues]<- (rnorm(n=numOfTrues,mean=truesMu1,sd=sd)) 
	#define the null first families
	xbar[(numOfTrues+1):size]<- (rnorm(n=size-numOfTrues,mean=falsesMu1,sd=sd)) 
	#define the non-null last families   
	return (xbar)
}

calcPVals<-function(xbars){
	return (1-pnorm(xbars)) #calculate p-vals
}

rejectBH <-function (pvals, details=FALSE, alpha=0.05){
	#function that gets a vector of p-values
	#sorts it, and returns the number of rejections
	#using BH
	# details- whether to print out details or not
	# pvals is the vector of pvals
	# returns a list with 2 vector: 
	#length is the # of rejections, ix is the rejections indices
	sortedPVals<-sort(pvals, index.return=TRUE)
	i<-length(pvals)
	#print (sortedPVals)
	while (sortedPVals$x[i]>(alpha*i/length (sortedPVals$x))&& i>0){
		i<- i-1
	}
	
	if (details==TRUE){
		print (c("Method: Benjamini Hochberg, First non-rejected value: ",sortedPVals$x[i]))  
		print (c("q=",alpha))
		print (c("Number of rejected hypothesis:",i))
		print (sortedPVals$ix)
	}
	if (i==0)
		return (list("length"=0,"ix"=NA))
	else 
		return (list("length"=i, "ix"=sortedPVals$ix[1:i]))
}

Preject <- function (method="BH", pvals=c(0), details=FALSE, alpha =0.05){
	a<-p.adjust(p=pvals,method=method)
	b<- which(a<alpha)
	return (list ("length"=length(b), "ix"=b))
}


rejectBon<- function (pvals, details=FALSE, alpha=0.05){
	i=1;
	m=length(pvals)
	if (m==0) return ("array of length 0")
	rejects<-numeric()
	while (i<=m){
		
		if (pvals[[i]]<=(alpha/m)){
			rejects<-append (rejects,i)
		} 
		i=i+1
	}
	#print (i)
	if (length(rejects)==0)
		return (list("length"=0,"ix"=NA))
	else 
		return (list("length"=length(rejects), "ix"=rejects))
	
}

countFalseRejections<-function (rejects , refVec){
	if (rejects$length==0){ #if there were no rejects at all
		return (0)
	}
	
	fr<-0 #false rejections counter
	for (i in 1:rejects$length){ # count num of false rejects
		if (refVec[rejects$ix[i]]==TRUE) { #if the rejection was false
			fr<-fr+1
		}
	}
	return (fr) 
}

calcPower<-function (rejects, falseRejections, size,numOfTrues){
	return ((rejects$length-falseRejections)/(size-numOfTrues))
}

calcFDR<- function (rejects,falseRejections){
	if (rejects$length==0) return (0)
	else 
		return (falseRejections/rejects$length)
}

calcFWER<- function (rejects, falseRejections){
	if (rejects$length==0)
		return (0)
	if (falseRejections>0) 
		return (1)
	if (falseRejections==0) 
		return (0)
}

##the test
iterationsTest<- function (n=10000,
		minMu1=0, maxMu1=4, interval=0.5){
	
	iters<-list()
	muIterFDR=list();
	muIterPOWER=list();
	muIterFRs=list();
	muIterFWERs=list();
	for (j in 1:2){
		muIterFDR[[j]]=numeric();
		muIterPOWER[[j]]=numeric();
		muIterFRs[[j]]=numeric();
		muIterFWERs[[j]]=numeric();
	}
	
	
	for (mu1 in seq(minMu1,maxMu1,interval)){
		#totalFDRs=numeric();
		#totalPOWERs=numeric();
		#totalFRs=numeric();
		#totalFWERs=numeric();
		totalFDRs=list();
		totalPOWERs=list();
		totalFRs=list();
		totalFWERs=list();
		
		for (j in 1:2){
			totalFDRs[[j]]=numeric(n);
			totalPOWERs[[j]]=numeric(n);
			totalFRs[[j]]=numeric(n);
			totalFWERs[[j]]=numeric(n);
		}
		
		
		for (i in 1:n){
			iters[i]<- iteration(
					size=100,
					numOfTrues=90,
					falsesMu1=mu1
			)
			iters[[i]]@refVec=setRefVector(iters[[i]]@size,iters[[i]]@numOfTrues)
			iters[[i]]@xbars=getRandomXBars(iters[[i]]@size,iters[[i]]@numOfTrues,falsesMu1=iters[[i]]@falsesMu1)
			iters[[i]]@pvals=calcPVals(iters[[i]]@xbars)
			iters[[i]]@proceduresApplied=append (iters[[i]]@proceduresApplied, new ("procedure", procedureName=rejectBH))
			iters[[i]]@proceduresApplied=append (iters[[i]]@proceduresApplied, new ("procedure", procedureName=rejectBon))
			
			for (j in 1:length(iters[[i]]@proceduresApplied)){
				iters[[i]]@proceduresApplied[[j]]@rejects=iters[[i]]@proceduresApplied[[j]]@procedureName(iters[[i]]@pvals)
				iters[[i]]@proceduresApplied[[j]]@fr=countFalseRejections(iters[[i]]@proceduresApplied[[j]]@rejects,iters[[i]]@refVec)
				iters[[i]]@proceduresApplied[[j]]@power=calcPower(iters[[i]]@proceduresApplied[[j]]@rejects, iters[[i]]@proceduresApplied[[j]]@fr, iters[[i]]@size, iters[[i]]@numOfTrues)
				iters[[i]]@proceduresApplied[[j]]@fdr=calcFDR(iters[[i]]@proceduresApplied[[j]]@rejects, iters[[i]]@proceduresApplied[[j]]@fr)
				iters[[i]]@proceduresApplied[[j]]@fwer=calcFWER(iters[[i]]@proceduresApplied[[j]]@rejects, iters[[i]]@proceduresApplied[[j]]@fr)
				totalFDRs[[j]][i]=iters[[i]]@proceduresApplied[[j]]@fdr  
				totalPOWERs[[j]][i]=iters[[i]]@proceduresApplied[[j]]@power
				totalFRs[[j]][i]=iters[[i]]@proceduresApplied[[j]]@fr
				totalFWERs[[j]][i]= iters[[i]]@proceduresApplied[[j]]@fwer
			} #end for j
		} #end for i
		
		#FDRs <- unlist(lapply(iters, function(iteration)iteration@fdr))
		#FWERs <-unlist(lapply(iters, function(iteration)iteration@fwer))
		#POWERs <-unlist(lapply(iters, function(iteration)iteration@power))
		#FRs <- unlist(lapply(iters, function(iteration)iteration@falseRejections))
		#print (totalFDRs)
		
		print(mu1)
		for (j in 1:length(iters[[i]]@proceduresApplied)){
			muIterFDR[[j]]=append(muIterFDR[[j]],mean(totalFDRs[[j]]))
			muIterPOWER[[j]]=append(muIterPOWER[[j]],mean(totalPOWERs[[j]]))
			muIterFRs[[j]]=append(muIterFRs[[j]],mean(totalFRs[[j]]))
			muIterFWERs[[j]]=append(muIterFWERs[[j]],mean(totalFWERs[[j]]))
			
		} #end For j
		#print (muIterFRs)
	} #end for mu
	
	
	
	plot(seq(minMu1,maxMu1,interval), muIterPOWER[[1]], type="o", col ="green",
			main =c("POWER/mu1 in", as.character(substitute(procedure)), n, "iterations"))
	
	for (j in 2:length(iters[[i]]@proceduresApplied)){
		lines(seq(minMu1,maxMu1,interval), muIterPOWER[[j]], type="o", col ="red",
				main =c("POWER/mu1 in", as.character(substitute(procedure)), n, "iterations"))
		
	}
	legend("topleft", 
			legend= c("BH", "Bon"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("green","red")
	)
	
	plot(seq(minMu1,maxMu1,interval), muIterFDR[[1]], type="o", col ="green", 
			ylim=c(0,1),
			main =c("FDR/mu1 in", as.character(substitute(procedure)), n, "iterations"))
	
	for (j in 2:length(iters[[i]]@proceduresApplied)){
		lines(seq(minMu1,maxMu1,interval), muIterFDR[[j]], type="o", col ="red",
				ylim=c(0,1),
				main =c("FDR/mu1 in", as.character(substitute(procedure)), n, "iterations"))
		
	}
	legend("topleft", 
			legend= c("BH", "Bon"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("green","red")
	)
	
	plot(seq(minMu1,maxMu1,interval), muIterFRs[[1]], type="o", col ="green",
			ylim=c(0,1),
			main =c("E[V]/mu1 in", as.character(substitute(procedure)), n, "iterations"))
	
	for (j in 2:length(iters[[i]]@proceduresApplied)){
		lines(seq(minMu1,maxMu1,interval), muIterFRs[[j]], type="o", col ="red",
				ylim=c(0,1),
				main =c("E[V]/mu1 in", as.character(substitute(procedure)), n, "iterations"))
		
	}
	legend("topleft", 
			legend= c("BH", "Bon"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("green","red")
	)
	
	plot(seq(minMu1,maxMu1,interval), muIterFWERs[[1]], type="o", col ="green",
			ylim=c(0,1),
			main =c("FWER/mu1 in", as.character(substitute(procedure)), n, "iterations"))
	
	for (j in 2:length(iters[[i]]@proceduresApplied)){
		lines(seq(minMu1,maxMu1,interval), muIterFWERs[[j]], type="o", col ="red",
				ylim=c(0,1),
				main =c("FWER/mu1 in", as.character(substitute(procedure)), n, "iterations"))
		
	}
	legend("topleft", 
			legend= c("BH", "Bon"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("green","red")
	)
	
	return (list(muIterFDR,muIterPOWER,muIterFWERs,muIterFRs)) 
	
} #end function

compareMethods<-function (n=10000, method1,method2,
		minMu1=0, maxMu1=4, interval=0.5){
	method1stats=iterationsTest(n=n, procedure=method1)
	method2stats=iterationsTest(n=n, procedure=method2)
	plot(seq(minMu1,maxMu1,interval), method1stats@fwer, type="o", col ="blue",
			main =c("FWER /mu1 in", as.character(substitute(method1)), as.character(substitute(method2)), n, "iterations"))
	lines(seq(minMu1,maxMu1,interval), method2stats@fwer, type="o", col ="red")
	legend("topleft", 
			legend= c(as.character(substitute(method1)),as.character(substitute(method2))), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("blue","red")
	)
	plot(seq(minMu1,maxMu1,interval), method1stats@fdr, type="o", col ="blue",
			main =c("FDR /mu1 in", as.character(substitute(method1)), as.character(substitute(method2)), n, "iterations"))
	lines(seq(minMu1,maxMu1,interval), method2stats@fdr, type="o", col ="red")
	legend("topleft", 
			legend= c(as.character(substitute(method1)),as.character(substitute(method2))), 
			lty=c(1,1),
			lwd=c(2,2),
			col=c("blue","red")
	)
} 

###simpified test- receives list of families, each family has K pvals

simpleTest<-function (ListOfFamilies, familySelectionAlpha=0.05, insideFamilyAlpha=0.05){
  #setup  of variables
  numOfFamilies=length(ListOfFamilies);
  familyArray= list(numOfFamilies);
  familiesSimesPVals= numeric(numOfFamilies);
  familiesBonPVals= numeric(numOfFamilies);
  SelectedFamilies<-list(0);
  alpha=0.05

  # generate families data
  for (i in 1:numOfFamilies){
    familyArray[[i]]= new ("iteration",
                           
                           falsesMu1=3
    )
#     familyArray[[i]]@refVec=setRefVector(familyArray[[i]]@size,familyArray[[i]]@numOfTrues)
#     familyArray[[i]]@xbars=getRandomXBars(familyArray[[i]]@size,familyArray[[i]]@numOfTrues,falsesMu1=familyArray[[i]]@falsesMu1)
    familyArray[[i]]@pvals=ListOfFamilies[[i]]
    print (familyArray[[i]]@pvals)
    familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                              new ("procedure"))
    familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                              new ("procedure"))
    familiesSimesPVals[i]=calcSimesPVals(familyArray[[i]]@pvals)
#     familiesBonPVals[i]=calcBonPVals(familyArray[[i]]@pvals)
    print(familiesSimesPVals[i])
  }
  
  #method 1 step 1:
  SelectedFamilies[[1]]<-SelectFamiliesBH(familiesSimesPVals)
  #method 2 step 1:
#   SelectedFamilies[[2]]<-SelectFamiliesBH(familiesBonPVals)

  print ( SelectedFamilies[[1]])
#   print (SelectedFamilies[[2]])
  
  
  
  ##method 1/2 step 2: #process BH/Bon on selected families
  for (methodix in 1:1){
    if (SelectedFamilies[[methodix]]$length>0){	
      for (i in 1:SelectedFamilies[[methodix]]$length){
        
        if (methodix==1){
          #print ("1")
          # reject inside the selected families with BH alpha=0.05*r/m
          familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[1]]@rejects<-
            Preject(method="BH", pvals=familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@pvals, 
                    alpha=alpha*SelectedFamilies[[methodix]]$length/numOfFamilies)

        }
#         if (methodix==2){
          #print ("2")
          # reject inside the selected families with bon alpha=0.05*r/m
#           familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[2]]@rejects<-
#             Preject(method="bon", pvals=familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@pvals, 
#                     alpha=alpha*SelectedFamilies[[methodix]]$length/numOfFamilies)
#         }
        print (familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[1]]@rejects)
        #calc frs,power,fdr,fwer for the ith selectedfamily
#         familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr=
#           countFalseRejections(
#             familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
#             familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@refVec)
#         familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@power=
#           calcPower(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
#                     familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr, 
#                     familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@size, 
#                     familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@numOfTrues)
#         familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr=
#           calcFDR(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
#                   familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr)
#         familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer=
#           calcFWER(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
#                    familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr)
#         
        #print (c("Family number:", SelectedFamilies[[methodix]]$ix[i],
        #         "FRs:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr,
        #         "FDR:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr,
        #        "FWER: ", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer,
        #         "Power: ", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@power
        #))
#         FDRSum[methodix]=FDRSum[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr
#         FWERSum[methodix]=FWERSum[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer					
#         OverallFR[methodix]=OverallFR[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr
#         OverallRejectsLength[methodix]=OverallRejectsLength[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@
#         proceduresApplied[[methodix]]@rejects$length
        #######################################
        #if (methodix==2){
        # 	print (c
        #           ("FALSE REJECTIONS M 2,  family #:", i))
        #   print(c("REJECTS",
        #	         familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects$length))
        #   print (c(
        #           "FRS",
        #  	       familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr))
        #           
        #  }
        ###############################################
        
      } #end for loop of i selected families     
      
#       
#       TotalAVGFDR[methodix, iters]<-FDRSum[methodix]/SelectedFamilies[[methodix]]$length
#       TotalAVGFWER[methodix, iters]<-FWERSum[methodix]/SelectedFamilies[[methodix]]$length
#       TotalOvrPower[methodix, iters]<- (OverallRejectsLength[methodix]-OverallFR[methodix])/(40*(familySize-numOfTrues))
#       TotalOvrFR[methodix, iters]<-OverallFR[methodix]
#       TotalAVGFR[methodix, iters]<-OverallFR[methodix]/SelectedFamilies[[methodix]]$length
#       TotalOvrFDR[methodix, iters]<- OverallFR[methodix]/OverallRejectsLength[methodix]
#       TotalOvrFWER[methodix, iters]<-sign(TotalOvrFDR[methodix, iters])
      
      #print (c("Average FDR:", FDRSum/SelectedFamilies[[methodix]]$length))
      #print (c("Overall FDR: ", OverallFR/OverallRejectsLength))
      #print (c("Average FWER: ", FWERSum/SelectedFamilies[[methodix]]$length))
      #print (c("Overall FR: ", OverallFR))
      #print ("Interesting Families (Bonferroni) are")
      #print (SelectedFamilies[[methodix]])
      #return (list(familiesBonPVals,familiesSimesPVals, SelectedFamilies[[methodix]], SelectedFamilies[[methodix]]))
      
#       
#     } #end if selected families exist
#     else {
#       #"no selected families in iteration"
#       TotalAVGFDR[methodix,iters]<-0
#       TotalAVGFWER[methodix,iters]<-0
#       TotalOvrFWER[methodix,iters]<-0
#       TotalOvrPower[methodix,iters]<-0
#       TotalOvrFR[methodix,iters]<-0
#       TotalAVGFR[methodix,iters]<-0
#       TotalOvrFDR[methodix,iters]<-0
#     }
  } #end for methodix
#   
#   # methods 3/4
#   for (methodix in 3:4){
#     if (methodix==3){
#       jointFamily@proceduresApplied[[methodix-2]]@rejects=Preject(method="BH", pvals=jointFamily@pvals, alpha=0.05)
#     }
#     if (methodix==4){
#       jointFamily@proceduresApplied[[methodix-2]]@rejects=Preject(method="bon", pvals=jointFamily@pvals, alpha=0.05)
#     }
#     
#     
#     jointFamily@proceduresApplied[[methodix-2]]@fr=countFalseRejections(jointFamily@proceduresApplied[[methodix-2]]@rejects,jointFamily@refVec)
#     jointFamily@proceduresApplied[[methodix-2]]@power=calcPower(jointFamily@proceduresApplied[[methodix-2]]@rejects,
#                                                                 jointFamily@proceduresApplied[[methodix-2]]@fr, jointFamily@size, jointFamily@numOfTrues)
#     jointFamily@proceduresApplied[[methodix-2]]@fdr=calcFDR(jointFamily@proceduresApplied[[methodix-2]]@rejects,
#                                                             jointFamily@proceduresApplied[[methodix-2]]@fr)
#     jointFamily@proceduresApplied[[methodix-2]]@fwer=calcFWER(jointFamily@proceduresApplied[[methodix-2]]@rejects, 
#                                                               jointFamily@proceduresApplied[[methodix-2]]@fr)
#     
#     TotalOvrFDR[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@fdr
#     TotalOvrFWER[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@fwer					
#     TotalOvrFR[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@fr
#     TotalOvrPower[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@power						
#     
  }
  #			print ("TOTAL OVR FDR")
  #			print (TotalOvrFDR)
  
}


####### second test

secondTest <-function(n=10000, numOfFamilies=40, familySize=100, numOfTrues=90, 
                      minMu1=0, maxMu1=4, interval=0.5, alpha=0.05){
  
  OvrPowerMeans<-matrix(0,4,maxMu1*2+1)
  AVGFDRMeans<-matrix(0,4,maxMu1*2+1)
  OvrFDRMeans<-matrix(0,4,maxMu1*2+1)
  AVGFWERMeans<-matrix(0,4,maxMu1*2+1)
  OvrFWERMeans<-matrix(0,4,maxMu1*2+1)
  OvrFRMeans<-matrix(0,4,maxMu1*2+1)
  AVGFRMeans<-matrix(0,4,maxMu1*2+1)
  
  for (mu in seq(0,maxMu1,0.5)){
    
    TotalAVGFDR<-matrix(0,2,n)
    TotalAVGFWER<-matrix(0,2,n)
    TotalOvrFWER<-matrix(0,4,n)
    TotalOvrPower<-matrix(0,4,n)
    TotalAVGFR<-matrix(0,2,n)
    TotalOvrFR<-matrix(0,4,n)
    TotalOvrFDR<-matrix(0,4,n)
    print (mu)
    
    for (iters in 1:n){
      
      
      #setup  of variables
      familyArray= list(numOfFamilies);
      familiesSimesPVals= numeric(numOfFamilies);
      familiesBonPVals= numeric(numOfFamilies);
      SelectedFamilies<-list(0);
      OverallRejectsLength<-numeric(2)
      OverallFR <- numeric(2)
      FDRSum <-numeric(2)
      FWERSum <- numeric(2)
      jointFamily <- new ("iteration",
                          size=familySize*numOfFamilies,
                          numOfTrues=numOfTrues*numOfFamilies,
                          falsesMu1=mu
      )
      
      # generate families data
      for (i in 1:numOfFamilies){
        familyArray[[i]]= new ("iteration",
                               size=familySize,
                               numOfTrues=numOfTrues,
                               falsesMu1=mu
        )
        familyArray[[i]]@refVec=setRefVector(familyArray[[i]]@size,familyArray[[i]]@numOfTrues)
        familyArray[[i]]@xbars=getRandomXBars(familyArray[[i]]@size,familyArray[[i]]@numOfTrues,falsesMu1=familyArray[[i]]@falsesMu1)
        familyArray[[i]]@pvals=calcPVals(familyArray[[i]]@xbars)
        familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                                  new ("procedure"))
        familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                                  new ("procedure"))
        familiesSimesPVals[i]=calcSimesPVals(familyArray[[i]]@pvals)
        familiesBonPVals[i]=calcBonPVals(familyArray[[i]]@pvals)
        
        #add to jointFamily for methods 3-4
        jointFamily@refVec=c(jointFamily@refVec,familyArray[[i]]@refVec)
        jointFamily@xbars=c(jointFamily@xbars,familyArray[[i]]@xbars)
        jointFamily@pvals=c(jointFamily@pvals,familyArray[[i]]@pvals)
        jointFamily@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                             new ("procedure"))
        jointFamily@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                             new ("procedure"))
      }
      
      #method 1 step 1:
      SelectedFamilies[[1]]<-SelectFamiliesBH(familiesSimesPVals)
      #method 2 step 1:
      SelectedFamilies[[2]]<-SelectFamiliesBH(familiesBonPVals)
      
      
      
      
      ##method 1/2 step 2: #process BH/Bon on selected families
      for (methodix in 1:2){
        if (SelectedFamilies[[methodix]]$length>0){	
          for (i in 1:SelectedFamilies[[methodix]]$length){
            
            if (methodix==1){
              #print ("1")
              # reject inside the selected families with BH alpha=0.05*r/m
              familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[1]]@rejects<-
                Preject(method="BH", pvals=familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@pvals, 
                        alpha=alpha*SelectedFamilies[[methodix]]$length/numOfFamilies)
            }
            if (methodix==2){
              #print ("2")
              # reject inside the selected families with bon alpha=0.05*r/m
              familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[2]]@rejects<-
                Preject(method="bon", pvals=familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@pvals, 
                        alpha=alpha*SelectedFamilies[[methodix]]$length/numOfFamilies)
            }
            
            #calc frs,power,fdr,fwer for the ith selectedfamily
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr=
              countFalseRejections(
                familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@refVec)
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@power=
              calcPower(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                        familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr, 
                        familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@size, 
                        familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@numOfTrues)
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr=
              calcFDR(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                      familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr)
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer=
              calcFWER(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                       familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr)
            
            #print (c("Family number:", SelectedFamilies[[methodix]]$ix[i],
            #         "FRs:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr,
            #         "FDR:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr,
            #        "FWER: ", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer,
            #         "Power: ", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@power
            #))
            FDRSum[methodix]=FDRSum[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr
            FWERSum[methodix]=FWERSum[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer					
            OverallFR[methodix]=OverallFR[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr
            OverallRejectsLength[methodix]=OverallRejectsLength[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@
            proceduresApplied[[methodix]]@rejects$length
            #######################################
            #if (methodix==2){
            # 	print (c
            #           ("FALSE REJECTIONS M 2,  family #:", i))
            #   print(c("REJECTS",
            #	         familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects$length))
            #   print (c(
            #           "FRS",
            #  	       familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr))
            #           
            #  }
            ###############################################
            
          } #end for loop of i selected families     
          
          
          TotalAVGFDR[methodix, iters]<-FDRSum[methodix]/SelectedFamilies[[methodix]]$length
          TotalAVGFWER[methodix, iters]<-FWERSum[methodix]/SelectedFamilies[[methodix]]$length
          TotalOvrPower[methodix, iters]<- (OverallRejectsLength[methodix]-OverallFR[methodix])/(40*(familySize-numOfTrues))
          TotalOvrFR[methodix, iters]<-OverallFR[methodix]
          TotalAVGFR[methodix, iters]<-OverallFR[methodix]/SelectedFamilies[[methodix]]$length
          TotalOvrFDR[methodix, iters]<- OverallFR[methodix]/OverallRejectsLength[methodix]
          TotalOvrFWER[methodix, iters]<-sign(TotalOvrFDR[methodix, iters])
          
          #print (c("Average FDR:", FDRSum/SelectedFamilies[[methodix]]$length))
          #print (c("Overall FDR: ", OverallFR/OverallRejectsLength))
          #print (c("Average FWER: ", FWERSum/SelectedFamilies[[methodix]]$length))
          #print (c("Overall FR: ", OverallFR))
          #print ("Interesting Families (Bonferroni) are")
          #print (SelectedFamilies[[methodix]])
          #return (list(familiesBonPVals,familiesSimesPVals, SelectedFamilies[[methodix]], SelectedFamilies[[methodix]]))
          
          
        } #end if selected families exist
        else {
          #"no selected families in iteration"
          TotalAVGFDR[methodix,iters]<-0
          TotalAVGFWER[methodix,iters]<-0
					TotalOvrFWER[methodix,iters]<-0
					TotalOvrPower[methodix,iters]<-0
					TotalOvrFR[methodix,iters]<-0
					TotalAVGFR[methodix,iters]<-0
					TotalOvrFDR[methodix,iters]<-0
				}
			} #end for methodix
			
			# methods 3/4
			for (methodix in 3:4){
				if (methodix==3){
					jointFamily@proceduresApplied[[methodix-2]]@rejects=Preject(method="BH", pvals=jointFamily@pvals, alpha=0.05)
				}
				if (methodix==4){
					jointFamily@proceduresApplied[[methodix-2]]@rejects=Preject(method="bon", pvals=jointFamily@pvals, alpha=0.05)
				}
				
				
				jointFamily@proceduresApplied[[methodix-2]]@fr=countFalseRejections(jointFamily@proceduresApplied[[methodix-2]]@rejects,jointFamily@refVec)
				jointFamily@proceduresApplied[[methodix-2]]@power=calcPower(jointFamily@proceduresApplied[[methodix-2]]@rejects,
												jointFamily@proceduresApplied[[methodix-2]]@fr, jointFamily@size, jointFamily@numOfTrues)
				jointFamily@proceduresApplied[[methodix-2]]@fdr=calcFDR(jointFamily@proceduresApplied[[methodix-2]]@rejects,
												jointFamily@proceduresApplied[[methodix-2]]@fr)
				jointFamily@proceduresApplied[[methodix-2]]@fwer=calcFWER(jointFamily@proceduresApplied[[methodix-2]]@rejects, 
												jointFamily@proceduresApplied[[methodix-2]]@fr)
										
				TotalOvrFDR[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@fdr
				TotalOvrFWER[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@fwer					
				TotalOvrFR[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@fr
				TotalOvrPower[methodix,iters]=jointFamily@proceduresApplied[[methodix-2]]@power						
				
			}
#			print ("TOTAL OVR FDR")
#			print (TotalOvrFDR)
		} #end n for - iterations
			
		#print ("Total AVG FDR")
		#print (mean(TotalAVGFDR))
		#print ("Total AVG FWER")
		#print (mean(TotalAVGFWER))
		#print ("Total Ovr FR")
		#print (mean(TotalOvrFR))
		#print ("Total OVr FDR")
		#print (mean(TotalOvrFDR))
		
		for (methodix in 1:4){
			if (methodix<3){
				AVGFDRMeans[methodix,mu*2+1]<-mean(TotalAVGFDR[methodix,])
				AVGFWERMeans[methodix,mu*2+1]<-mean(TotalAVGFWER[methodix,])
				AVGFRMeans[methodix,mu*2+1]<-mean(TotalAVGFR[methodix,])
			}
			OvrPowerMeans[methodix,mu*2+1]<-mean(TotalOvrPower[methodix,])
			OvrFDRMeans[methodix,mu*2+1]<-mean(TotalOvrFDR[methodix,])
			OvrFWERMeans[methodix,mu*2+1]<-mean(TotalOvrFWER[methodix,])
			OvrFRMeans[methodix,mu*2+1]<-mean(TotalOvrFR[methodix,])
			
#			print ('Ovr fdrmeans')
#			print (OvrFDRMeans)
		
		}#end methodix for 
	} #end mu for
	
	
	
#	print (OvrFDRMeans[1,])
#	print (OvrFDRMeans[2,])
#	print (AVGFDRMeans[1,])
#	print (AVGFDRMeans[2,])
#	
	plot (seq(0,maxMu1,0.5),OvrPowerMeans[1,], 
			col="blue", ylim=c(0,1),xlab= "mu", ylab= "Overall Power",
		main =c("Overall Power/mu1 in", n, "iterations with methods 1-2"), type="o")
	lines (seq(0,maxMu1,0.5),OvrPowerMeans[2,], col="red", type="o")
	lines (seq(0,maxMu1,0.5),OvrPowerMeans[3,], col="green", type="o")
	lines (seq(0,maxMu1,0.5),OvrPowerMeans[4,], col="pink", type="o")
	
	
legend("topleft", 
		legend= c("Overall Power method 1", "Overall Power method 2", 
				"Overall Power method 3", "Overall Power method 4"), 
		lty=c(1,1,1,1),
		lwd=c(2,2,2,2),
		col=c("blue","red", "green", "pink")
	)
	plot (seq(0,maxMu1,0.5),AVGFDRMeans[1,], 
			col="blue", ylim=c(0,0.15),xlab= "mu", ylab= "FDR",
			main =c("FDRs/mu1 in", n, "iterations"), type="o")
	lines (seq(0,maxMu1,0.5),OvrFDRMeans[1,], col="red", type="o")
	lines (seq(0,maxMu1,0.5),AVGFDRMeans[2,], col="green", type="o")
	lines (seq(0,maxMu1,0.5),OvrFDRMeans[2,], col="pink", type="o")
	lines (seq(0,maxMu1,0.5),OvrFDRMeans[3,], col="brown", type="o")
	lines (seq(0,maxMu1,0.5),OvrFDRMeans[4,], col="cyan", type="o")
	abline (h=0.05, col="black")
	
	legend("topleft", 
			legend= c("AVG FDR over the selected method 1", "Overall FDR method 1", "AVG FDR over the selected method 2", "Overall FDR method 2", 
					"Overall FDR method 3", "Overall FDR method 4"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("blue","red", "green", "pink", "brown", "cyan")
	)
	plot (seq(0,maxMu1,0.5),AVGFWERMeans[1,], 
			col="blue", ylim=c(0,1.5),xlab= "mu", ylab= "FWER",
			main =c("FWERs/mu1 in", n, "iterations"), type="o")
	lines (seq(0,maxMu1,0.5),OvrFWERMeans[1,], col="red", type="o")
	lines (seq(0,maxMu1,0.5),AVGFWERMeans[2,], col="green", type="o")
	lines (seq(0,maxMu1,0.5),OvrFWERMeans[2,], col="pink", type="o")
	lines (seq(0,maxMu1,0.5),OvrFWERMeans[3,], col="brown", type="o")
	lines (seq(0,maxMu1,0.5),OvrFWERMeans[4,], col="cyan", type="o")
	abline (h=0.05, col="black")
	
	legend("topleft", 
			legend= c("AVG FWER over the selected method 1", "Overall FWER method 1", "AVG FWER over the selected method 2", "Overall FWER method 2", 
					"Overall FWER method 3", "Overall FWER method 4"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("blue","red", "green", "pink", "brown", "cyan")
	)
	
	plot (seq(0,maxMu1,0.5),AVGFRMeans[1,], 
			col="blue", ylim=c(0,25),xlab= "mu", ylab= "FWER",
			main =c("E[V]s/mu1 in", n, "iterations"), type="o")
	lines (seq(0,maxMu1,0.5),OvrFRMeans[1,], col="red", type="o")
	lines (seq(0,maxMu1,0.5),AVGFRMeans[2,], col="green", type="o")
	lines (seq(0,maxMu1,0.5),OvrFRMeans[2,], col="pink", type="o")
	lines (seq(0,maxMu1,0.5),OvrFRMeans[3,], col="brown", type="o")
	lines (seq(0,maxMu1,0.5),OvrFRMeans[4,], col="cyan", type="o")
	abline (h=0.05, col="black")
	
	legend("topleft", 
			legend= c("AVG E[V] over the selected method 1", "Overall E[V] method 1", "AVG E[V] over the selected method 2", "Overall E[V] method 2", 
					"Overall E[V] method 3", "Overall E[V] method 4"), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("blue","red", "green", "pink", "brown", "cyan")
	)
	plot (seq(0,maxMu1,0.5),AVGFWERMeans[1,] , col="blue", ylim=c(0,0.7),
			ylab="FWER/E[V]", xlab="mu",
			main =c("Avg FWER and Avg E[V]/mu1 in", n, "iterations with methods 1,2"), type="o")
	lines (seq(0,maxMu1,0.5),AVGFRMeans[1,], col="red", type="o")
	lines (seq(0,maxMu1,0.5),AVGFWERMeans[2,], col="green", type="o")
	lines (seq(0,maxMu1,0.5),AVGFRMeans[2,], col="pink", type="o")

	legend("topleft", 
			legend= c("Average FWER method 1","E[V] over the selected method 1",
                "Average FWER method 2","E[V] over the selected method 2"	), 
			lty=c(1,1,1,1),
			lwd=c(2,2,2,2),
			col=c("blue","red", "green", "pink")
	)
	print (AVGFRMeans)
}


calcSimesPVals<- function (PVals){
	min(p.adjust(PVals, method="BH"))
}
calcBonPVals<- function (PVals){
	min(p.adjust(PVals, method="bon"))
}

SelectFamiliesBH<-function (familiesSimesPVals){
	#b<-p.adjust(familiesSimesPVals,method="BH")
	#return (which(b<0.05))
	#print ("SELECT: ")
	#print (Preject(familiesSimesPVals))
	#print ("BH SELECT")
	return (Preject(method="BH", pvals=familiesSimesPVals))
	
}  
#return Array of Selected Families and old ixs

SelectFamiliesBon<- function (familiesBonPVals){
	#rejects=rejectBon(familiesBonPVals,alpha=0.05/length(familiesBonPVals))
	#return (rejects$ix)
	#b<-p.adjust(familiesBonPVals,method="bon")
	#return (which(b<0.05))
	#print ("BON SELECT")
	return (Preject(method="bon", pvals=familiesBonPVals))
}
#return Array of Selected Families and old ixs


}
# secondTest(n=100)
