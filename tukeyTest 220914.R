iteration <- setClass ("iteration",
                       #Data structure that holds data regarding a certain iteration (=manipulation of random numbers)        
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
                       #a data structure holding the results of certain procedure applied to an iteration
                       
                       slots=c(procedureName="function",
                               rejects="list",
                               fr="numeric",
                               power="numeric",
                               fdr="numeric",
                               fwer="numeric"
                       )
)

calcSimesPVals<- function (PVals){
  min(p.adjust(PVals, method="BH"))
}

calcPower<-function (rejectsLength, falseRejections, size,numOfFalseH0sInRef,details=FALSE){
  #returns num of  true rejections divided by the number of planned rejections
  if (details==TRUE){
    print ("POWER CALC:")
    print (c("rejects:", rejectsLength))
    print (c("fr:", falseRejections))
    print (c("size:", size))
    print (c("numOfFalseH0sInRef:", numOfFalseH0sInRef))
    print (c("Power",(rejectsLength-falseRejections)/(size-numOfFalseH0sInRef)))
  }
  return ((rejectsLength-falseRejections)/(size-numOfFalseH0sInRef))
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

countFalseRejections<-function (rejects , refVec){
  if (rejects$length==0){ #if there were no rejects at all
    return (0)
  }
  
  fr<-0 #false rejections counter
  for (i in 1:rejects$length){ # count num of false rejects
    #     print (refVec)
    if (refVec[rejects$ix[i]]==FALSE) { #if the rejection was false
      fr<-fr+1
    }
  }
  return (fr) 
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



rejectTukey <- function (xbars,S,groupSize,qStar){
  statistics=as.numeric(dist(xbars)/sqrt(S/groupSize))
  b<-which(statistics>=qStar)
  
  #   print (statistics)
  #   print ("Q*")
  #   print (qStar)
  #   print (list("length"=length(b), "ix"=b))
  return (list("length"=length(b), "ix"=b))
}

Preject <- function (method="BH", pvals=c(0), details=FALSE, alpha =0.05){
  a<-p.adjust(p=pvals,method=method)
  #   print (a)
  b<- which(a<alpha)
  return (list ("length"=length(b), "ix"=b))
}

calcTukeyPVal <-function(xbars,S,groupSize, numOfGroups){
  pval<- 1-ptukey(q=max(dist(xbars))/sqrt(S/groupSize), 
                  nmeans=numOfGroups ,df=numOfGroups * (groupSize-1))
  return (pval)
}

calcPairwisePVals <-function (xbars,S,groupSize, numOfGroups,i=1,details=FALSE){
  ##useful for method C phase 2
  pvals=2*(1-pt(q=dist(xbars)/sqrt(2*S[i]/groupSize),
                df= numOfGroups*(groupSize-1)))
  if (details==TRUE){
    print (c("PAIRWISE", pvals))
    print (as.numeric(pvals))
  }
  #   readline()
  return (as.numeric(pvals))
}

getRandomXBars <-function (size, numOfTrues, truesMu1=0, falsesMu1, sd=1){
  xbar<- rep(NA,size) #init xbar
  xbar[1:numOfTrues]<- (rnorm(n=numOfTrues,mean=truesMu1,sd=sd)) 
  #define the null first families
  xbar[(numOfTrues+1):size]<- (rnorm(n=size-numOfTrues,mean=falsesMu1,sd=sd)) 
  #define the non-null last families   
  #   print (xbar)
  #   readline()
  return (xbar)
}
setDesVector <-function (size, numOfZeros,mu=4){
  a<-rep(0,numOfZeros)
  b<-rep(mu,size-numOfZeros)
  return (c(a,b))
}

setRefVectorBig <-function(size, numOfZeros,mu){
  DesVector<-setDesVector(size,numOfZeros,mu)
  #   print ("DesVector:")
  #   print (DesVector)
  return (as.logical(dist(DesVector)))
}

# this method works only for numOfGroups=3, setRefVectorBig does this right for larger sizes
# setRefVector <- function (size, numOfTrues){
#   #initialize reference vector
#   refvec<-rep (NA,size)
#   refvec[1:numOfTrues]<- TRUE # these will be true discoveries
#   refvec[(numOfTrues+1):size]<- FALSE # these will be False Discoveries
#   return (refvec)
# }


#method A(3): Overall BH (one level)
#method B(2): QTukey Statistics (2 levels) 1: TukeyPval , 2: QTukey Stat
#method C(1): BH with pairwise Pval (2 levels) 1: TukeyPVal, 2:Pairwise PVal
#method 4: To be added later - currently empty

#this algorithm builds "numOfFamilies" families, each family of size "numOfGroups",
#with "numOfTrues" trues and the rest of the families are false, with "numOfSignalFamilies"
#families in which there 
#is signal (=numOfTrues is used), and the rest of the families are with no signal (all mus are 0)
#then uses 3 different methods to reject.
tukeyTest2 <-function(n=10000, numOfFamilies=40, numOfGroups=3, numOfTrues=1, numOfSignalFamilies=5,
                      minMu1=0, maxMu1=4, interval=0.5, alpha=0.05, groupSize=16,details=FALSE){
  
  #definition of constant df
  df=numOfGroups*(groupSize-1) #calculate degs of freedom
  
  #declration of matrices to keep means, size is # of methods * the number of mu values 
  OvrPowerMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  AVGFDRMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  OvrFDRMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  AVGFWERMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  OvrFWERMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  OvrFRMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  AVGFRMeans<-matrix(0,4,(maxMu1-minMu1)/interval+1)
  
  for (mu in seq(minMu1,maxMu1,interval)){
    
    #matrices for each mu1
    TotalAVGFDR<-matrix(0,4,n)
    TotalAVGFWER<-matrix(0,4,n)
    TotalOvrFWER<-matrix(0,4,n)
    TotalOvrPower<-matrix(0,4,n)
    TotalAVGFR<-matrix(0,4,n)
    TotalOvrFR<-matrix(0,4,n)
    TotalOvrFDR<-matrix(0,4,n)
    print (mu)
    #     readline()
    
    #each iteration
    for (iters in 1:n){
      #       print (c("iter:" ,iters))
      
      #setup  of variables
      familyArray= list(numOfFamilies);
      familiesTukeyPVals = numeric(numOfFamilies);
      familiesSimesPVals = numeric(numOfFamilies);
      
      SelectedFamilies<-list(0);
      OverallRejectsLength<-numeric(4) # 4 slots for each method B/C/*/4
      OverallFR <- numeric(4)
      FDRSum <-numeric(4)
      FWERSum <- numeric(4)
      jointFamily <- new ("iteration",
                          size=choose(numOfGroups,2)*numOfFamilies,
                          numOfTrues=0+
                            #num of true hypotheses in the no signal zone
                            numOfSignalFamilies*(numOfGroups-1),
                          #num of true hypotheses in the signal zone
                          falsesMu1=mu
      )
      #       print (jointFamily)
      #       readline()
      ##init variance vector
      S=rchisq(n=numOfFamilies,df=df)/df
      
      for (i in 1:numOfFamilies){
        
        familyArray[[i]]= new ("iteration",
                               size=choose (numOfGroups,2),
                               numOfTrues=numOfTrues, #=1
                               falsesMu1=mu
        ) #init family
        
        
        if (i<=numOfSignalFamilies){ #first 5 families with signal
          familyArray[[i]]@refVec=setRefVectorBig(
            size=numOfGroups,
            numOfZeros=numOfGroups-numOfTrues,
            mu=mu)
          #           print ("Reference vector:")
          #           print (familyArray[[i]]@refVec)
          #build reference Vector
          familyArray[[i]]@xbars=getRandomXBars(size=numOfGroups,
                                                numOfTrues=numOfGroups-1,
                                                falsesMu1=mu,
                                                sd=sqrt(1/groupSize))
        }
        else { #35 families no signal
          familyArray[[i]]@refVec=setRefVectorBig(
            size=numOfGroups,
            numOfZeros=numOfGroups,
            mu=mu)
          
          #build reference Vector (we set numOfTrues as 1 and falsesMu1 as 0 so that all xbars are
          #generated around 0, the value of numOfTrues doesnt matter here)
          familyArray[[i]]@xbars=getRandomXBars(size=numOfGroups,
                                                numOfTrues=1,
                                                falsesMu1=0,
                                                sd=sqrt(1/groupSize))
        }
        #                             print ("Reference vector:")
        #                             print (familyArray[[i]]@refVec)
        if (details==TRUE){
          #                    print (i)
          #                    print (familyArray[[i]]@xbars)         
        } #end details print
        
        familyArray[[i]]@pvals=calcPairwisePVals(familyArray[[i]]@xbars,
                                                 S,
                                                 groupSize,
                                                 numOfGroups,
                                                 i=i,
                                                 details=FALSE)
        #        readline()
        #calc Pairwise PVals for Method C(1)
        
        familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                                  new ("procedure"))
        familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                                  new ("procedure"))
        familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                                  new ("procedure"))
        familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                                  new ("procedure"))
        familiesTukeyPVals[i]=calcTukeyPVal(familyArray[[i]]@xbars
                                            ,S[i],groupSize, numOfGroups)
        familiesSimesPVals[i]=calcSimesPVals(familyArray[[i]]@pvals)
        
        #calc Tukey PVals for method B(2)
        
        #         print (c("family tukey pval", i, familiesTukeyPVals[i]))
        
        #add to jointFamily for methods 3-4
        jointFamily@refVec=append(jointFamily@refVec,familyArray[[i]]@refVec)
        jointFamily@xbars=append(jointFamily@xbars,familyArray[[i]]@xbars)
        jointFamily@pvals=append(jointFamily@pvals,familyArray[[i]]@pvals)
        jointFamily@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                             new ("procedure"))
        jointFamily@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                             new ("procedure"))
      } #end for i
      
      if (details==TRUE){
        #         print (jointFamily)
        #         print (familiesSimesPVals)
        #         readline()
      }
      
      #method B,C step 1: select interesting families using BH
      SelectedFamilies[[1]]<-SelectFamiliesBH(familiesTukeyPVals)
      SelectedFamilies[[2]]<-SelectedFamilies[[1]] #for this stage selection for both methods is the same
      qStar=qtukey(p=1-alpha*SelectedFamilies[[1]]$length/numOfFamilies,
                   nmeans=numOfGroups,df=df)
      SelectedFamilies[[4]]<-SelectFamiliesBH(familiesSimesPVals)
      #       print (familyArray[[i]]@pvals, qStar)
      print ("SELECTS FROM METHOD 4")
      print (SelectedFamilies[[4]])
      #             readline()
      
      
      ##method 1/2 step 2: #process BH/QTukey on selected families
      for (methodix in 1:3){ #run through methods
        if (methodix==3){
          methodix<-4
        }
        print (c("METHOD:", methodix))
        if (SelectedFamilies[[methodix]]$length>0){   #if there are any selected families
          for (i in 1:SelectedFamilies[[methodix]]$length){ #run thru selected families
            
            if (methodix==1 || methodix==4){
              #method 1=method C
              #method 4= BH BH (simes)
              #print ("1")
              # reject inside the selected families with BH alpha=0.05*r/m
              familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects<-
                Preject(method="BH", pvals=familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@pvals, 
                        alpha=alpha*SelectedFamilies[[methodix]]$length/numOfFamilies)
            }
            if (methodix==2){
              #methodB
              #print ("2")
              # reject inside the selected families with QTUKEY threshold
              familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[2]]@rejects<-
                rejectTukey( familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@xbars,
                             S[SelectedFamilies[[methodix]]$ix[i]],
                             groupSize,
                             qStar)
              #               print (dist(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@xbars))
              #               print (familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@refVec)
              #               readline()
            }
            
            #calc frs,power,fdr,fwer for the ith selectedfamily
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr=
              countFalseRejections(
                familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@refVec)
            
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@power=
              calcPower(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects$length,
                        familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr, 
                        familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@size, 
                        length(which(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@refVec==FALSE)))
            
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr=
              calcFDR(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                      familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr)
            
            familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer=
              calcFWER(familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects,
                       familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr)
            if (details==TRUE){ 
              print (c("Family number:", SelectedFamilies[[methodix]]$ix[i]))
              print (c("FRs:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr))
              print (c("FDR:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr))
              print (c("FWER: ", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer))
              print (c("Power: ", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@power))
              print (c("Rejects:", familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@rejects$length))
            } #end details print
            FDRSum[methodix]=FDRSum[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fdr
            FWERSum[methodix]=FWERSum[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fwer  				
            OverallFR[methodix]=OverallFR[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[methodix]]@fr
            OverallRejectsLength[methodix]=OverallRejectsLength[methodix]+familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@
              proceduresApplied[[methodix]]@rejects$length
            
            
          } #end for loop of i selected families     
          #           print (c("OVERALL REJECTS LENGTH method :", methodix ,OverallRejectsLength[methodix]))
          
          
          TotalAVGFDR[methodix, iters] <-FDRSum[methodix]/SelectedFamilies[[methodix]]$length
          if (is.nan(TotalAVGFDR[methodix, iters])) {
            TotalAVGFDR[methodix, iters]<-0
          }
          TotalAVGFWER[methodix, iters]<-FWERSum[methodix]/SelectedFamilies[[methodix]]$length
          if (is.nan(TotalAVGFDR[methodix, iters])) {
            TotalAVGFWER[methodix, iters]<-0
          }
          print ("TOTAL OVR POWER CALC")
          print ("")
          TotalOvrPower[methodix, iters]<- calcPower(rejects=OverallRejectsLength[methodix],
                                                     falseRejections=OverallFR[methodix],
                                                     size=jointFamily@size,
                                                     numOfFalseH0sInRef=jointFamily@size-jointFamily@numOfTrues
          ) 
          
          #             (OverallRejectsLength[methodix]-OverallFR[methodix])/(numOfFamilies*(numOfGroups-numOfTrues))
          TotalOvrFR[methodix, iters]<-OverallFR[methodix]
          TotalAVGFR[methodix, iters]<-OverallFR[methodix]/SelectedFamilies[[methodix]]$length
          
          TotalOvrFDR[methodix, iters]<- OverallFR[methodix]/OverallRejectsLength[methodix]
          if (is.nan(TotalOvrFDR[methodix, iters])) {
            TotalOvrFDR[methodix, iters]<-0
          }
          
          TotalOvrFWER[methodix, iters]<-sign(TotalOvrFDR[methodix, iters])
          if (is.nan(TotalOvrFWER[methodix, iters])) {
            TotalOvrFWER[methodix, iters]<-0
          }
          
          print (c("Iteration:", iters, "mu1:", mu, "method:", methodix))
          print (c("Average FDR:", FDRSum[methodix]/SelectedFamilies[[methodix]]$length))
          print (c("Overall FDR: ",  OverallFR[methodix]/OverallRejectsLength[methodix]))
          print (c("Average FWER: ", FWERSum[methodix]/SelectedFamilies[[methodix]]$length))
          print (c("Overall FR: ", OverallFR[methodix]))
          print (c("Interesting Families () are",SelectedFamilies[[methodix]]))
          #readline()
          
          
        } #end if selected families exist
        else {
          #"no selected families in iteration - all stats  are 0"
          print ("no families selected")
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
      for (methodix in 3:3){
        if (methodix==3){
          #method 3 - overall BH - take joint family and reject in it using BH
          jointFamily@proceduresApplied[[methodix-2]]@rejects=Preject(method="BH", pvals=jointFamily@pvals, alpha=0.05)
        }
        
        #calc stats for joint family
        jointFamily@proceduresApplied[[methodix-2]]@fr=countFalseRejections(jointFamily@proceduresApplied[[methodix-2]]@rejects,jointFamily@refVec)
        jointFamily@proceduresApplied[[methodix-2]]@power=calcPower(jointFamily@proceduresApplied[[methodix-2]]@rejects$length,
                                                                    jointFamily@proceduresApplied[[methodix-2]]@fr, jointFamily@size, jointFamily@size-jointFamily@numOfTrues)
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
      if (methodix!=3){
        AVGFDRMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalAVGFDR[methodix,])
        AVGFWERMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalAVGFWER[methodix,])
        AVGFRMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalAVGFR[methodix,])
      }
      OvrPowerMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalOvrPower[methodix,])
      OvrFDRMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalOvrFDR[methodix,])
      OvrFWERMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalOvrFWER[methodix,])
      OvrFRMeans[methodix,(mu-minMu1)/interval+1]<-mean(TotalOvrFR[methodix,])
      
      #			print ('Ovr fdrmeans')
      #			print (OvrFDRMeans)
      
    }#end methodix for 
  } #end mu for
  
  
  
  #	print (OvrFDRMeans[1,])
  #	print (OvrFDRMeans[2,])
  #	print (AVGFDRMeans[1,])
  #	print (AVGFDRMeans[2,])
  #	
  
  plot (seq(minMu1,maxMu1,interval),OvrPowerMeans[1,], 
        col="blue", ylim=c(0,1),xlab= "mu", ylab= "Overall Power",
        main =c("Overall Power/mu1 in", n, "iterations with methods 1-4"), type="o")
  lines (seq(minMu1,maxMu1,interval),OvrPowerMeans[2,], col="red", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrPowerMeans[3,], col="green", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrPowerMeans[4,], col="pink", type="o")
  
  
  legend("topleft", 
         legend= c("Overall Power method 1", "Overall Power method 2", 
                   "Overall Power method 3","Overall Power method 4"), 
         lty=c(1,1,1,1),
         lwd=c(2,2,2,2),
         col=c("blue","red", "green", "pink")
  )
  plot (seq(minMu1,maxMu1,interval),AVGFDRMeans[1,], 
        col="blue", ylim=c(0,0.2),xlab= "mu", ylab= "FDR",
        main =c("FDRs/mu1 in", n, "iterations"), type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFDRMeans[1,], col="red", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFDRMeans[2,], col="green", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFDRMeans[2,], col="pink", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFDRMeans[3,], col="brown", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFDRMeans[4,], col="cyan", type="o")
  abline (h=0.05, col="black")
  
  legend("topleft", 
         legend= c("AVG FDR over the selected method 1", "Overall FDR method 1",
                   "AVG FDR over the selected method 2", "Overall FDR method 2", 
                   "Overall FDR method 3",
                   "AVG FDR over the selected method 4", "Overall FDR method 4" ), 
         lty=c(1,1,1,1),
         lwd=c(2,2,2,2),
         col=c("blue","red", "green", "pink", "brown","cyan")
  )
  plot (seq(minMu1,maxMu1,interval),AVGFWERMeans[1,], 
        col="blue", ylim=c(0,1),xlab= "mu", ylab= "FWER",
        main =c("FWERs/mu1 in", n, "iterations"), type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFWERMeans[1,], col="red", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFWERMeans[2,], col="green", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFWERMeans[2,], col="pink", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFWERMeans[3,], col="brown", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFWERMeans[4,], col="grey", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFWERMeans[4,], col="cyan", type="o")
  abline (h=0.05, col="black")
  
  legend("topleft", 
         legend= c("AVG FWER over the selected method 1", "Overall FWER method 1",
                   "AVG FWER over the selected method 2", "Overall FWER method 2", 
                   "Overall FWER method 3", 
                   "AVG FWER over the selected method 4","Overall FWER method 4"), 
         lty=c(1,1,1,1),
         lwd=c(2,2,2,2),
         col=c("blue","red", "green", "pink", "brown","grey", "cyan")
  )
  
  plot (seq(minMu1,maxMu1,interval),AVGFRMeans[1,], 
        col="blue", ylim=c(0,3),xlab= "mu", ylab= "E[V]",
        main =c("E[V]s/mu1 in", n, "iterations"), type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFRMeans[1,], col="red", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFRMeans[2,], col="green", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFRMeans[2,], col="pink", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFRMeans[3,], col="brown", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFRMeans[4,], col="yellow", type="o")
  lines (seq(minMu1,maxMu1,interval),OvrFRMeans[4,], col="grey", type="o")
  abline (h=0.05, col="black")
  
  legend("topleft", 
         legend= c("AVG E[V] over the selected method 1", "Overall E[V] method 1", 
                   "AVG E[V] over the selected method 2", "Overall E[V] method 2", 
                   "Overall E[V] method 3",
                   "AVG E[V] over the selected method 4", "Overall E[V] method 4"), 
         lty=c(1,1,1,1),
         lwd=c(2,2,2,2),
         col=c("blue","red", "green", "pink", "brown","yellow","grey")
  )
  plot (seq(minMu1,maxMu1,interval),AVGFWERMeans[1,] , col="blue", ylim=c(0,0.5),
        ylab="FWER/E[V]", xlab="mu",
        main =c("Avg FWER and Avg E[V]/mu1 in", n, "iterations with methods 1,2,4"), type="o")
  #    readline()
  lines (seq(minMu1,maxMu1,interval),AVGFRMeans[1,], col="red", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFWERMeans[2,], col="green", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFRMeans[2,], col="pink", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFWERMeans[4,], col="yellow", type="o")
  lines (seq(minMu1,maxMu1,interval),AVGFRMeans[4,], col="grey", type="o")
  legend("topleft", 
         legend= c("Average FWER method 1","E[V] over the selected method 1",
                   "Average FWER method 2","E[V] over the selected method 2",
                   "Average FWER method 4","E[V] over the selected method 4"), 
         lty=c(1,1,1,1),
         lwd=c(2,2,2,2),
         col=c("blue","red", "green", "pink","yellow","grey")
  )
  
  #change column names in table
#   colnames(OvrPowerMeans)<-seq(minMu1,maxMu1,interval)
  colnames(OvrFDRMeans)<-seq(minMu1,maxMu1,interval)
  colnames(AVGFDRMeans)<-seq(minMu1,maxMu1,interval)
  colnames(OvrFWERMeans)<-seq(minMu1,maxMu1,interval)
  colnames(OvrFRMeans)<-seq(minMu1,maxMu1,interval)
  colnames(AVGFRMeans)<-seq(minMu1,maxMu1,interval)
  colnames(AVGFWERMeans)<-seq(minMu1,maxMu1,interval)
  #rownames(OvrPowerMeans)<-c("QTukey Stat", "Pairwise", "Overall BH", "BH BH")
  
  print ("OVERALL POWER:")
  print (OvrPowerMeans)
#   sprintf("%1.7f", OvrPowerMeans)
#   format (OvrPowerMeans,digits=10)
#  formatC( OvrPowerMeans, format='f',width=8, digits=4 )
#   
  print ("AVG FDR:")
  print( AVGFDRMeans)
  print ("OVERALL FDR:")
  print( OvrFDRMeans)
  print ("OVERALL FWER:")
  print(OvrFWERMeans)
  print ("OVERALL E[V]:")
  print(OvrFRMeans)
  print ("AVG E[V]: ")
  print (AVGFRMeans)
  print ("AVG FWER: ")
  print (AVGFWERMeans)
  
  wb<-loadWorkbook(paste("TukeyTest",Sys.Date(),".xls"), create = TRUE)
  createSheet(wb, name = "OVERALL POWER")
  writeWorksheet(wb, OvrPowerMeans, sheet = "OVERALL POWER")
  
  createSheet(wb, name = "AVG FDR")
  writeWorksheet(wb, AVGFDRMeans, sheet = "AVG FDR")
  createSheet(wb, name = "OVERALL FDR")
  writeWorksheet(wb, OvrFDRMeans, sheet = "OVERALL FDR")
  createSheet(wb, name = "OVERALL FWER")
  writeWorksheet(wb, OvrFWERMeans, sheet = "OVERALL FWER")
  createSheet(wb, name = "OVERALL E(V)")
  writeWorksheet(wb, OvrFRMeans, sheet = "OVERALL E(V)")
  createSheet(wb, name = "AVG E(V)")
  writeWorksheet(wb, AVGFRMeans, sheet = "AVG E(V)")
  createSheet(wb, name = "AVG FWER")
  writeWorksheet(wb, AVGFWERMeans, sheet = "AVG FWER")
  saveWorkbook(wb)
}

#tukeyTest2(n=100)
tukeyTest2(n=1, numOfGroups=5,minMu1=0, maxMu1=9,details=FALSE)
