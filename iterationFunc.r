SelectFamiliesBH<-function (familiesSimesPVals){
  return (Preject(method="BH", pvals=familiesSimesPVals))
  
}  
#return Array of Selected Families and old ixs


OverallBHSelectedFamilies<-function (jointRejects,numOfGroups){
  newIXs<-unique(as.integer((jointRejects$ix-1)/(choose(numOfGroups,2))))+1
  return (list ("length"=length(newIXs), "ix"=newIXs))
  
}


rejectTukey <- function (xbars,S,groupSize,qStar){
<<<<<<< HEAD
 
#    print ("*****AS NUMERIC CHECK*******************")
#   print ("xbars=")
#   print (xbars)
#   print ("S=")
#   print (S)
#   print ("Groupsize=")
#   print (groupSize)
#   print ("qstar")
#   print (qStar)
#   print(dist(xbars)/sqrt(S/groupSize))
#   print (as.numeric(dist(xbars)/sqrt(S/groupSize)))
#   
  
  statistics=as.numeric(dist(xbars)/sqrt(S/groupSize))
  b<-which(statistics>=qStar)
  
  # print (list("length"=length(b), "ix"=b))
  
=======
  statistics=as.numeric(dist(xbars)/sqrt(S/groupSize))
  b<-which(statistics>=qStar)
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
  return (list("length"=length(b), "ix"=b))
}

Preject <- function (method="BH", pvals=c(0), details=FALSE, alpha =0.05){
  a<-p.adjust(p=pvals,method=method)
  b<- which(a<alpha)
  return (list ("length"=length(b), "ix"=b))
}


getRandomXBars <-function (size=-1, numOfTrues=-1, truesMu1=0, falsesMu1, sd=1,DesVector=NaN,details=FALSE){
  #if there's a Design vector supplied, normals will be generated according to it
  # otherwise the xbars will be made with the given parameters.
  
  if (is.nan(DesVector)[1]){
    xbar<- rep(NA,size) #init xbar
    xbar[1:numOfTrues]<- (rnorm(n=numOfTrues,mean=truesMu1,sd=sd)) 
    #define the null first families
    xbar[(numOfTrues+1):size]<- (rnorm(n=size-numOfTrues,mean=falsesMu1,sd=sd)) 
    #define the non-null last families  
  }
  else{
    xbar<- rep(NA,length(DesVector)) #init xbar
    xbar<- (rnorm(n=length(DesVector),mean=DesVector,sd=sd))   
  }
  if (details==TRUE){
    print (DesVector)
    print (xbar)  
  }
  
  #   readline()
  return (xbar)
}

setRefVectorBig <-function(size, numOfZeros,mu, DesVector=NaN, details=FALSE){
  
  if (is.nan(DesVector)[1]){
    DesVector<-setDesVector(size,numOfZeros,mu)
  }
  
  if (details==TRUE){
    print ("DesVector:")
    print (DesVector)
    print ("RefVector:")
    print (as.logical(dist(DesVector)))
  }
  return (as.logical(dist(DesVector)))
}
countFalseRejections<-function (rejects , refVec,details=FALSE){
  if (details==TRUE){
    print ("count FR:")
    print(c("Rejects:" , rejects, "Ref vec:"))
    print (rbind(refVec))
    
  }
  if (rejects$length==0){ #if there were no rejects at all
<<<<<<< HEAD
    print ("0 rejects")
    return (0)
  }
  
=======
     print ("0 rejects")
    return (0)
  }

>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
  fr<-0 #false rejections counter
  for (i in 1:rejects$length){ # count num of false rejects
    #         print (c("len", length(refVec)))
    #         print (c("rejects ix:" , rejects$ix[i]))
    if (refVec[rejects$ix[i]]==FALSE) { #if the rejection was false
      fr<-fr+1
    }
  }
  # print(c("Total false rejections:" , fr))
  if (fr>0){
    print(c("Rejects:" , rejects, "Ref vec:"))
    print (refVec)
  }
  return (fr) 
}


countFalseFamilyRejections <- function (rejects, numOfSignalFamilies){
  if (rejects$length==0){ #if there were no rejects at all
    return (0)
  }
  print ("false family rejections:")
  print (which(rejects$ix>numOfSignalFamilies))
  return (length(which(rejects$ix>numOfSignalFamilies)))
}


calcPower<-function (rejectsLength, falseRejections, size,numOfFalseH0sInRef,details=FALSE){
  #returns num of  true rejections divided by the number of planned rejections
  #Prof Benjamini said: if a family has no signal then its' power is NaN.
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

calcTukeyPVals <-function(xbars,S,groupSize, numOfGroups, details=FALSE){
  pvals<- numeric(nrow(xbars))
  for (i in 1:nrow(xbars)){
    pvals[i]<- 1-ptukey(q=max(dist(xbars[i,]))/sqrt(S[i]/groupSize), 
                        nmeans=numOfGroups ,df=numOfGroups * (groupSize-1))
  }
  if (details) print (pvals)
  return (pvals)
}

calcPairwisePVals <-function (xbars,S,groupSize, numOfGroups,details=FALSE){
  ##useful for method C phase 2
  
  pvals=matrix(1, nrow(xbars), choose (ncol(xbars),2))
  #print(pvals)
  for (i in 1:nrow(xbars)){
    pvals[i,]=as.numeric(2*(1-pt(q=dist(xbars[i,])/sqrt(2*S[i]/groupSize),
                                 df= numOfGroups*(groupSize-1))))
  }
  if (details==TRUE){
    print (c("PAIRWISE", pvals))
  }
  #   readline()
  return (pvals)
}

calcSimesPVals<- function (PVals){
  #gets a matrix of pairwise-pvals as input
  #outputs the min of every row in a vector of simes pvals
  
  # print (PVals[i,])
<<<<<<< HEAD
  
=======
 
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
  simes=numeric(nrow(PVals))
  for (i in 1:nrow(PVals)){
    simes[i]=min(p.adjust(PVals[i,], method="BH"))
  }
  return (simes)
  
}


iteration <-function  (xbars, numOfSignalFamilies, numOfGroups, groupSize, delta,methodix, DesVector){
  alpha=0.05
  #definintion of statsMatrix
  numOfFamilies=nrow(xbars);
  statsMatrix=matrix(NA, numOfFamilies+1, 7)
  colnames(statsMatrix)<-c("FAMILY_#","METHOD_#","SELECTED", "FDR", "FWER", "FR", "POWER")
  statsMatrix[1:numOfFamilies,"FAMILY_#"]=1:numOfFamilies
  statsMatrix[numOfFamilies+1,"FAMILY_#"]=0;
  statsMatrix[,"METHOD_#"]=methodix
  statsMatrix[,"SELECTED"]<-FALSE
  
  #make refernce vector
  refVec=setRefVectorBig(
    size=numOfGroups,
    numOfZeros=numOfGroups-numOfTrues,
    mu=delta,
    DesVector=DesVector,
    details=FALSE)
  noSignalRefVec=rep(FALSE,choose(numOfGroups,2))
  
  #definition of constants df and S
  df=numOfGroups*(groupSize-1) #calculate degs of freedom
  S=rchisq(n=numOfFamilies,df=df)/df
  
  #1. multiply xbars by delta
  xbars=delta*xbars
  
  #2. create PVals accordign to the desired method
  if (methodix==1 || methodix == 3 || methodix == 4){
    pairwisePVals=calcPairwisePVals(xbars,
                                    S,
                                    groupSize,
                                    numOfGroups,
                                    details=FALSE)
<<<<<<< HEAD
#     print(pairwisePVals)
#     print (refVec)
    }
=======
  }
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
  if (methodix==1 || methodix == 2){
    familiesTukeyPVals=calcTukeyPVals(xbars,S,groupSize, 
                                      numOfGroups,details=FALSE)
    #         print (familiesTukeyPVals)
  }
  if (methodix==4){
    print(pairwisePVals)
    familiesSimesPVals=calcSimesPVals(pairwisePVals)
    print(familiesSimesPVals)
  }
  
  #3. select families according to desired method
  if (methodix==1 || methodix==2){
    #methods 1/2: select using Tukey pvals
    SelectedFamilies<-SelectFamiliesBH(familiesTukeyPVals)
    print ("Selected Families:")
    print (SelectedFamilies)
  }
  
  if (methodix==3){
    #method 3 - overall BH - take joint family and reject in it using BH
<<<<<<< HEAD
    
    #IMPORTANT: MUST transpose a matrix if we use as.numeric to create a signle vector 
    #out of it (when each row of the matrix represents a family) 
    bigFamilyRejects=Preject(method="BH", pvals=as.numeric(t(pairwisePVals)), alpha=0.05)
=======
    bigFamilyRejects=Preject(method="BH", pvals=as.numeric(pairwisePVals), alpha=0.05)
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
    SelectedFamilies<-OverallBHSelectedFamilies(bigFamilyRejects, numOfGroups)
    print ("method3 bigfamily rejects")
    print (bigFamilyRejects)
    print ("selected families 3:")      
    print (SelectedFamilies)
  }
  if (methodix==4){
    #method 4: select using Simes pvals
    SelectedFamilies<-SelectFamiliesBH(familiesSimesPVals)
  }
  
  #4. calculate big family FR, FDR,FWER and POWER
  
  statsMatrix[numOfFamilies+1,"FR"]<-countFalseFamilyRejections(SelectedFamilies,numOfSignalFamilies)
  statsMatrix[numOfFamilies+1,"POWER"]<- calcPower(rejects=SelectedFamilies$length,
                                                   falseRejections=statsMatrix[numOfFamilies+1,"FR"],
                                                   size=numOfFamilies,
                                                   numOfFalseH0sInRef=numOfFamilies-numOfSignalFamilies,
                                                   details=FALSE)
  statsMatrix[numOfFamilies+1,"FDR"]<-calcFDR(rejects=SelectedFamilies,
                                              falseRejections=statsMatrix[numOfFamilies+1,"FR"])
  statsMatrix[numOfFamilies+1,"FWER"]<-calcFWER(rejects=SelectedFamilies,
                                                falseRejections=statsMatrix[numOfFamilies+1,"FR"])
  
  #5. Select inside selected families according to the desired method
  #save rejects into rejects.
  
  if (methodix==2){
    #prepare qStar for second phase selection
    qStar=qtukey(p=1-alpha*SelectedFamilies$length/numOfFamilies,
                 nmeans=numOfGroups,df=df)
  }
  for (i in 1:(SelectedFamilies$length)){ #run thru selected families
    #update matrix to indicate this family has been rejected
    #print (SelectedFamilies$ix)
    
    #change refVector after finishing iterating thru signal families
    #to an all FALSE ref vec.
    if (i>numOfSignalFamilies){
      refVec=noSignalRefVec;
<<<<<<< HEAD
    }
=======
  }
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
    statsMatrix[SelectedFamilies$ix[i],"SELECTED"]<-TRUE;
    
    if (methodix==1 || methodix==4){
      #method 1=Tukey BH
      #method 4= BH BH (simes)
      
      # reject inside the selected families with BH alpha=0.05*r/m
      
      rejects<- Preject(method="BH", pvals=pairwisePVals[SelectedFamilies$ix[i],], 
                        alpha=alpha*SelectedFamilies$length/numOfFamilies)
    }
    if (methodix==2){
      #method B
      #print ("2")
      # reject inside the selected families with QTUKEY threshold
<<<<<<< HEAD
     
      # print (xbars[SelectedFamilies$ix[i],])
      
      rejects<-
        rejectTukey(xbars[SelectedFamilies$ix[i],],
=======
      rejects<-
        rejectTukey(xbars[SelectedFamilies$ix[i]],
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
                    S[SelectedFamilies$ix[i]],
                    groupSize,
                    qStar)
      
    }
    
    if (methodix==3){
      #find the subset of the jointfamily which is corresponding to the i'th family rejects.
      
      familyNum<-SelectedFamilies$ix[i]
      
      print (c("familynum", familyNum, "numoFfamilies", numOfFamilies))
      m3rejects<-bigFamilyRejects$ix
      #               print (m3rejects)
      #               print (numOfGroups)
      
      print (i)
      m3rejects<-m3rejects[m3rejects> (familyNum-1)*choose(numOfGroups,2) &
                             m3rejects<= (familyNum)*choose(numOfGroups,2)] 
      #                             print (m3rejects)
      m3rejects<- (m3rejects %% choose(numOfGroups,2))
      #                             print (m3rejects)
      #                             print (m3rejects[m3rejects==0]);
      m3rejects[m3rejects==0]<- choose(numOfGroups,2)                         
      #                             print ("after modulus")
      #                             print ("Method 3 rejects")
      #                                           print (m3rejects)
      #               
      
      rejects<-list("length"=length(m3rejects), "ix"=m3rejects)       
    }
    
    
    #calc frs,power,fdr,fwer for the ith selectedfamily
    
    statsMatrix[SelectedFamilies$ix[i], "FR"]=
      countFalseRejections(rejects=rejects, refVec=refVec,details=FALSE)
    # print ("fr finished")
    statsMatrix[SelectedFamilies$ix[i], "POWER"]=
      calcPower(rejects$length,
                statsMatrix[SelectedFamilies$ix[i], "FR"], 
                choose (ncol(xbars),2),
                length(which(refVec==FALSE)), 
                details=FALSE)       
    statsMatrix[SelectedFamilies$ix[i], "FDR"]=
      calcFDR(rejects, statsMatrix[SelectedFamilies$ix[i], "FR"])       
    statsMatrix[SelectedFamilies$ix[i], "FWER"]=
      calcFWER(rejects, statsMatrix[SelectedFamilies$ix[i], "FR"]) 
  }
  return (statsMatrix)
}

tukeyTestSplit<-function (n=10000, numOfFamilies=40, numOfGroups=3, numOfTrues=1, numOfSignalFamilies=5,
                          mindelta=0, maxdelta=4, interval=0.5, alpha=0.05, groupSize=16,details=FALSE,
                          DesVector=NaN,methodix=1){
  print(methodix)
  numOfGroups=length(DesVector)
  print ("Num Of Entered Groups:")
  print(numOfGroups)
  
  #first numOfSignalFamilies with signal
  xbars<-NULL
  #we get random XBars using sd equals to sqrt(1/groupSize) for the signal 
  #families
  
  for (i in 1:numOfSignalFamilies){
    xbars<-rbind(xbars,getRandomXBars(sd=sqrt(1/groupSize),
                                      DesVector=DesVector))
  }
  #no signal families
  for (i in 1:(numOfFamilies-numOfSignalFamilies)){
    xbars<-rbind(xbars,getRandomXBars(sd=sqrt(1/groupSize),
                                      DesVector=rep(0,length(DesVector)) ))
  }
<<<<<<< HEAD
  print ("XBARS")
  print (xbars)
=======
      print ("XBARS")
      print (xbars)
>>>>>>> 2008f6d3c2b875c5fd2901344f1b67eac6e8b4bc
  print (c("Number of Rows of Xbars:", nrow(xbars)))
  
  stats<-iteration(xbars=xbars, numOfSignalFamilies=numOfSignalFamilies, 
                   numOfGroups=numOfGroups, groupSize=groupSize, delta=1,methodix=methodix, DesVector=DesVector)
  #print only rows representing families that were rejected
  print (stats[!is.na(stats[,"FDR"]),])
  
}
