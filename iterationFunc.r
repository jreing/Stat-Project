SelectFamiliesBH<-function (familiesSimesPVals){
    return (Preject(method="BH", pvals=familiesSimesPVals))
    
}  
#return Array of Selected Families and old ixs


OverallBHSelectedFamilies<-function (jointRejects,numOfGroups){
    newIXs<-unique(as.integer((jointRejects$ix-1)/(choose(numOfGroups,2))))+1
    return (list ("length"=length(newIXs), "ix"=newIXs))
    
}

rejectTukey <- function (xbars,S,groupSize,qStar){
    
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
    
    return (list("length"=length(b), "ix"=b))
}

Preject <- function (method="BH", pvals=c(0), details=FALSE, alpha =0.05){
    a<-p.adjust(p=pvals,method=method)
    b<- which(a<alpha)
    return (list ("length"=length(b), "ix"=b))
}


getRandomXBars <-function (size=-1, numOfTrues=-1, truesMu1=0, 
                           falsesMu1, sd=1,DesVector=NaN,details=FALSE, rows=1){
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
        xbar<- (rnorm(n=length(DesVector)*rows,mean=DesVector,sd=sd))   
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
        # print ("0 rejects")
        return (0)
    }
    
    fr<-0 #false rejections counter
    for (i in 1:rejects$length){ # count num of false rejects
        #         print (c("len", length(refVec)))
        #         print (c("rejects ix:" , rejects$ix[i]))
        if (refVec[rejects$ix[i]]==FALSE) { #if the rejection was false
            fr<-fr+1
        }
    }
    # print(c("Total false rejections:" , fr))
#     if (fr>0){
#         print(c("Rejects:" , rejects, "Ref vec:"))
#         print (refVec)
#     }
    return (fr) 
}


countFalseFamilyRejections <- function (rejects, numOfSignalFamilies){
    if (rejects$length==0){ #if there were no rejects at all
        return (0)
    }
    # print ("false family rejections:")
    # print (which(rejects$ix>numOfSignalFamilies))
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
    
    simes=numeric(nrow(PVals))
    for (i in 1:nrow(PVals)){
        # print (PVals[i,])
        # readline()
        simes[i]=min(p.adjust(PVals[i,], method="BH"))
    }
    
    return (simes)
    
}


iteration <-function  (xbars, numOfSignalFamilies, numOfGroups,
                       groupSize, delta,methodix, DesVector){
    alpha=0.05
   
    
    if (is.null(methodix)){
        startix=1
        endix=4
        # print(startix)
        # print(endix)
        # readline()
    }
    else{
        startix=endix=methodix
    }
    #definItion of res
    res=list(NULL,NULL,NULL,NULL);
    numOfFamilies=nrow(xbars);
   
    noSignalRefVec=rep(FALSE,choose(numOfGroups,2))
    
    #definition of constants df and S
    df=numOfGroups*(groupSize-1) #calculate degs of freedom
    S=rchisq(n=numOfFamilies,df=df)/df
    
    #put null in variables
    pairwisePVals<-NULL
    familiesTukeyPVals<-NULL
    SelectedFamilies<-NULL
    
    #1. multiply xbars of signal families by delta
    for(i in 1:numOfSignalFamilies){
        xbars[i,]=(delta*DesVector)+xbars[i,] 
    }
    # print(xbars)
    # readline()
    
    #start methodix loop
    for (methodix in startix:endix){
        
        #make refernce vector
        refVec=setRefVectorBig(
            size=numOfGroups,
            numOfZeros=numOfGroups-numOfTrues,
            mu=delta,
            DesVector=DesVector,
            details=FALSE)
        
        #init 2nd Phase reject counter
        Overall2ndPhaseRejects=0;
        
        #init OverallNumOfFalseH0s (for OVR POWER calculation)
        OverallNumOfFalseH0s=numOfSignalFamilies*length(which(refVec==FALSE))+
            (numOfFamilies-numOfSignalFamilies)*choose (ncol(xbars),2)

        #definItion of statsMatrix
        statsMatrix=matrix(NA, numOfFamilies+1, 13)
        colnames(statsMatrix)<-c("FAMILY_#","METHOD_#","SELECTED", "FAM FDR", "FAM FWER", "FR", "POWER", "AVG FDR", "AVG FWER","OVR FR", "OVR POWER", "OVR FDR", "OVR FWER")
        statsMatrix[1:numOfFamilies,"FAMILY_#"]=1:numOfFamilies
        statsMatrix[numOfFamilies+1,"FAMILY_#"]=NaN;
        statsMatrix[,"METHOD_#"]=methodix
        statsMatrix[,"SELECTED"]<-FALSE
        # print(statsMatrix)
        
        #2. create PVals accordign to the desired method
        if (methodix==1 || methodix == 3 || methodix == 4){
            if (is.null(pairwisePVals)){
                pairwisePVals<-calcPairwisePVals(xbars,
                                                 S,
                                                 groupSize,
                                                 numOfGroups,
                                                 details=FALSE)
            }
            # print(pairwisePVals)
            #     print (refVec)
        }
        if (methodix==1 || methodix == 2){
            if (is.null(familiesTukeyPVals)){
                familiesTukeyPVals=calcTukeyPVals(xbars,S,groupSize, 
                                                  numOfGroups,details=FALSE)
            }
            # print (familiesTukeyPVals)
            # readline()
        }
        if (methodix==4){
            # print(pairwisePVals)
            familiesSimesPVals=calcSimesPVals(pairwisePVals)
            # print(familiesSimesPVals)
            # readline()
            
        }
        
        #3. select families according to desired method
        if (methodix==1 || methodix==2){
            #methods 1/2: select using Tukey pvals
            if (is.null(SelectedFamilies)){
                SelectedFamilies<-SelectFamiliesBH(familiesTukeyPVals)
                
            }
            #       print ("Selected Families:")
            #       print (SelectedFamilies)
            # readline()
        }
        
        if (methodix==3){
            #method 3 - overall BH - take joint family and reject in it using BH
            
            #IMPORTANT: MUST transpose a matrix if we use as.numeric to create a signle vector 
            #out of it (when each row of the matrix represents a family) 
            bigFamilyRejects=Preject(method="BH", pvals=as.numeric(t(pairwisePVals)), alpha=0.05)
            SelectedFamilies<-OverallBHSelectedFamilies(bigFamilyRejects, numOfGroups)
            #       print ("method3 bigfamily rejects")
            #       print (bigFamilyRejects)
            #       print ("selected families 3:")      
            #       print (SelectedFamilies)
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
        statsMatrix[numOfFamilies+1,"FAM FDR"]<-calcFDR(rejects=SelectedFamilies,
                                                        falseRejections=statsMatrix[numOfFamilies+1,"FR"])
        statsMatrix[numOfFamilies+1,"FAM FWER"]<-calcFWER(rejects=SelectedFamilies,
                                                          falseRejections=statsMatrix[numOfFamilies+1,"FR"])
        statsMatrix[numOfFamilies+1,"SELECTED"]<-SelectedFamilies$length
        
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
            }
            statsMatrix[SelectedFamilies$ix[i],"SELECTED"]<-TRUE;
            
            if (methodix==1 || methodix==4){
                #method 1=Tukey BH
                #method 4= BH BH (simes)
                
                # reject inside the selected families with BH alpha=0.05*r/m
                
                Level2Rejects<- Preject(method="BH", pvals=pairwisePVals[SelectedFamilies$ix[i],], 
                                        alpha=alpha*SelectedFamilies$length/numOfFamilies)
            }
            if (methodix==2){
                #method B
                #print ("2")
                # reject inside the selected families with QTUKEY threshold
                
                # print (xbars[SelectedFamilies$ix[i],])
                
                Level2Rejects<-
                    rejectTukey(xbars[SelectedFamilies$ix[i],],
                                S[SelectedFamilies$ix[i]],
                                groupSize,
                                qStar)
                
            }
            
            if (methodix==3){
                #find the subset of the jointfamily which is corresponding to the i'th family rejects.
                
                familyNum<-SelectedFamilies$ix[i]
                
                # print (c("familynum", familyNum, "numoFfamilies", numOfFamilies))
                m3rejects<-bigFamilyRejects$ix
                #               print (m3rejects)
                #               print (numOfGroups)
                
                # print (i)
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
                Level2Rejects<-list("length"=length(m3rejects), "ix"=m3rejects)       
            }
            
            
            #calc frs,power,fdr,fwer for the ith selectedfamily
            
            statsMatrix[SelectedFamilies$ix[i], "FR"]=
                countFalseRejections(rejects=Level2Rejects, refVec=refVec,details=FALSE)
            statsMatrix[SelectedFamilies$ix[i], "POWER"]=
                calcPower(Level2Rejects$length,
                          statsMatrix[SelectedFamilies$ix[i], "FR"], 
                          choose (ncol(xbars),2),
                          length(which(refVec==FALSE)), 
                          details=FALSE)       
            statsMatrix[SelectedFamilies$ix[i], "FAM FDR"]=
                calcFDR(Level2Rejects, statsMatrix[SelectedFamilies$ix[i], "FR"])       
            statsMatrix[SelectedFamilies$ix[i], "FAM FWER"]=
                calcFWER(Level2Rejects, statsMatrix[SelectedFamilies$ix[i], "FR"]) 
            
            Overall2ndPhaseRejects=Overall2ndPhaseRejects+Level2Rejects$length
            
        }
        #calculate AVGs
        statsMatrix[numOfFamilies+1,"AVG FDR"]=mean(statsMatrix[1:numOfFamilies,"FAM FDR"], na.rm=TRUE)
        statsMatrix[numOfFamilies+1,"AVG FWER"]=mean(statsMatrix[1:numOfFamilies,"FAM FWER"], na.rm=TRUE)
        statsMatrix[numOfFamilies+1,"OVR FR"]=sum(statsMatrix[1:numOfFamilies,"FR"], na.rm=TRUE)
        statsMatrix[numOfFamilies+1,"OVR FDR"]=calcFDR(rejects=list("length"=Overall2ndPhaseRejects, "ix"=NULL),
                                                       falseRejections= statsMatrix[numOfFamilies+1,"OVR FR"])
        statsMatrix[numOfFamilies+1,"OVR FWER"]=
            calcFWER(rejects= list ("length"=Overall2ndPhaseRejects, "ix"=NULL),
                     falseRejections = statsMatrix[numOfFamilies+1,"OVR FR"] )
        
        
#         print (c("Overall2ndPhaseRejects",Overall2ndPhaseRejects, "OVR FR", statsMatrix[numOfFamilies+1,"OVR FR"],
#                  "size " , choose (ncol(xbars),2)*numOfFamilies, 
#                  "numOfFalseH0sInRef", OverallNumOfFalseH0s))
        # readline()
        
        statsMatrix[numOfFamilies+1,"OVR POWER"]=
            calcPower(rejectsLength = Overall2ndPhaseRejects, 
                      falseRejections = statsMatrix[numOfFamilies+1,"OVR FR"],
                      size = choose (ncol(xbars),2)*numOfFamilies,
                      numOfFalseH0sInRef = OverallNumOfFalseH0s)
        if (SelectedFamilies$length==0){
            statsMatrix[numOfFamilies+1,"AVG FDR"]=0
            statsMatrix[numOfFamilies+1,"AVG FWER"]=0
        }
#          print("this is statsMatrix:")
#          print(statsMatrix[numOfFamilies+1,])
#          
         
        res[[methodix]]=statsMatrix[numOfFamilies+1,];
        
    } #end methodix loop
    res=rbind(res[[1]], res[[2]], res[[3]], res[[4]])
    
    # print(res)
    # readline()
    return (res)
}

tukeyTestSplit<-function (n=10000, numOfFamilies=40, numOfGroups=3, numOfTrues=1, numOfSignalFamilies=5,
                          mindelta=0, maxdelta=4, interval=0.5, alpha=0.05, groupSize=16,details=FALSE,
                          DesVector=NaN,methodix=1){
    print(methodix)
    numOfGroups=length(DesVector)
    print ("Num Of Entered Groups:")
    print(numOfGroups)
    
    #initilaize timer
    ptm<-proc.time()
    
    #   print ("XBARS")
    #   print (xbars)
    #   print (c("Number of Rows of Xbars:", nrow(xbars)))
    
    function.names = c(
        'calcFDR', 'calcFWER', 'calcSimesPVals', 'calcPairwisePVals',
        'calcPower','getRandomXBars',
        'calcTukeyPVals', 'countFalseFamilyRejections', 'countFalseRejections',
        'Preject', 'rejectTukey', 'OverallBHSelectedFamilies','SelectFamiliesBH',
        'setRefVectorBig', 'iteration'
    )

    library(parallel)
    library(doParallel)
    library(foreach)
    library(doRNG)
    library(XLConnect)
    nr.cores=4
    cl = makeCluster(nr.cores)
    registerDoParallel(cl)
   
    vecString=paste(DesVector, collapse= ",")
    FILENAME=paste(
        "TukeyTest ",
        format(Sys.time(), "%d%m%y%H%M"),
        "DesVector=",  vecString, 
        # "TotalTime:", totalTime,
        "NumOfSignalFamilies",numOfSignalFamilies,
        "NumOfFamilies", numOfFamilies,
        "n=",n,".xls")
    
    ##parallellized loop to gain all xbars to be used in simulation 
    random.res = foreach(j=1:n , .export=function.names , .packages='foreach') %dopar% {
        #first numOfSignalFamilies with signal
        
       
        #we get random XBars using sd equals to sqrt(1/groupSize) for the signal 
        #families
        
        noSigXBars<-matrix( getRandomXBars(
                                    sd=sqrt(1/groupSize),
                                    DesVector=DesVector, 
                                    rows=numOfSignalFamilies
                                    )
                            , nrow=numOfSignalFamilies
                            , ncol=numOfGroups, byrow=TRUE)
        # print(noSigXBars)
        
        #no signal families
        SigXBars<-matrix( getRandomXBars(
                                sd=sqrt(1/groupSize),
                                DesVector=rep(0,length(DesVector)),
                                rows =numOfFamilies-numOfSignalFamilies
                                ),
                          nrow =numOfFamilies-numOfSignalFamilies,
                          ncol =numOfGroups, byrow=TRUE)
        
        # print(SigXBars) 
    
        res=rbind(noSigXBars, SigXBars);
        
    }
    save (random.res, file = paste (FILENAME, "random xbars"))
	print ("Finished: recieved all random Xbars")
#     # print (random.res)
#     
#     load("xbars001a")
#     print ("xbars loaded")
	
    ##parallelized loop to iterate thru deltas and run the iters on each delta 
    ## and each xbar vector
	
    wb<-loadWorkbook(filename=FILENAME,
                     create = TRUE)
    print ("current Design Vector:")
    print (DesVector)
    deltaRes=foreach (delta= seq(mindelta,maxdelta,interval), .export=function.names) %dopar%{ 
        sumX=iteration(xbars=random.res[[1]], numOfSignalFamilies=numOfSignalFamilies, 
                          numOfGroups=numOfGroups, groupSize=groupSize, delta=delta,
                          methodix=NULL, DesVector=DesVector)
        sumX2=sumX^2
        
        #sd= sum(x_i^2)+2*mean*sum(x_i)+n*mean^2
        
        for (j in 2:n){
            print (j)
            #if (j %% 100 ==0) print (j)
            iter=iteration(xbars=random.res[[j]], numOfSignalFamilies=numOfSignalFamilies, 
                              numOfGroups=numOfGroups, groupSize=groupSize, delta=delta,
                              methodix=NULL, DesVector=DesVector)
            sumX=sumX+iter
            sumX2=sumX2+iter^2
            # }
#             save (wrapper.res, file = paste(FILENAME, "DATA PART" ,k ))
#             print (paste("data part ", k , "was saved"))
#             # print (wrapper.res)
#             rm(wrapper.res)
        }
        #wrapping it up together into means and SDs
        
        # l=length(wrapper.res)
        
        print ("means:")
        EofX=sumX/n
        EofX2=sumX2/n
        wrapper.res2=EofX
        wrapper.res2sd= sqrt(EofX2-EofX^2)
        
        deltaRes=list(wrapper.res2, wrapper.res2sd)
    }
    stopCluster(cl)
    
    #GOOD OUTPUT:
    print (deltaRes)
    save (deltaRes,file=paste(FILENAME,".data"))

    #EXPORT CURRENT RESULT TO EXCEL:
        
#     createSheet(wb, name = paste("MEANS", delta))
#     writeWorksheet(wb, round(wrapper.res2[is.na(wrapper.res2[,"FAMILY_#"]),],5), sheet =  paste("MEANS", delta))
#     createSheet(wb, name =  paste("SD", delta))
#     writeWorksheet(wb, round(wrapper.res2sd[is.na(wrapper.res2[,"FAMILY_#"]),],5), sheet =  paste("SD", delta))
#     
    
    # saveWorkbook(wb)
    print (c("# of iterations:" , n))  
    #stop timer
    totalTime= (proc.time()-ptm)[3]
    
    totalTime=paste(totalTime,collapse="")
    totalTime=strtrim(totalTime,6)
    print(strtrim(totalTime,6))
    # write.csv(simplify2array(wrapper.res),"output.csv")
    
}

# tukeyTestSplit(n=3,numOfSignalFamilies=m1, numOfFamilies=m,interval=0.5, mindelta=1, maxdelta=2
#                               ,DesVector=c(0,0,1), details=FALSE)
#                

N=10000
INTERVAL=0.3
MINDELTA=0.1
MAXDELTA=1

#Series Of Tests:
for (m1 in c(10,50,100,1000)){
   for (m in c(1000,10000)){
       # for (m in c(10000)){
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(0,0,1), details=FALSE)
#       gc()
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(-1,0,1), details=FALSE)
#     gc()
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(-1,0,1,2), details=FALSE)
#     gc()
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(-1,0,0,1), details=FALSE)
#     gc()
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(-1,-1,-1,1,1,1), details=FALSE)
#     gc()
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(-1,-1,0,0,1,1), details=FALSE)
#     gc()
#     tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
#                ,DesVector=c(-2,-1,0,0,1,2), details=FALSE)
#     gc()
    tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
               ,DesVector=c(-1,-1,-1,-1,-1,1,1,1,1,1), details=FALSE)
    gc()
    tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
              ,DesVector=c(-2,-2,-1,-1,0,0,1,1,2,2), details=FALSE)
    gc()
    tukeyTestSplit(n=N,numOfSignalFamilies=m1, numOfFamilies=m,interval=INTERVAL, mindelta=MINDELTA, maxdelta=MAXDELTA
               ,DesVector=c(-3,-2,-1,0,0,0,0,1,2,3), details=FALSE)
  }
}
