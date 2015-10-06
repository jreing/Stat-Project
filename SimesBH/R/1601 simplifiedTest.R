
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


###simpified test- receives list of families, each family has K pvals

x=list(c(0.005, 0.002, 0.001), c(0.08, 0.07, 0.06), c(0.02,0.01,0.03))

y=list(c(0, 0.06, 0.07), c(0,0.055, 0.065), c(0,0.075,0.08))


SimesBH<-function (ListOfFamilies, familySelectionAlpha=0.05, insideFamilyAlpha=0.05,echo="OFF"){
  #setup  of variables
  numOfFamilies=length(ListOfFamilies);
  familyArray= list(numOfFamilies);
  familiesSimesPVals= numeric(numOfFamilies);
  familiesBonPVals= numeric(numOfFamilies);
  SelectedFamilies<-list(0);
  
  # generate families data
  for (i in 1:numOfFamilies){
    familyArray[[i]]= new ("iteration")
    familyArray[[i]]@pvals=ListOfFamilies[[i]]
    if (echo=="ON") print (c("Family", i," pvals:", familyArray[[i]]@pvals))
    familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                              new ("procedure"))
    familyArray[[i]]@proceduresApplied=append(familyArray[[i]]@proceduresApplied,
                                              new ("procedure"))
    familiesSimesPVals[i]=calcSimesPVals(familyArray[[i]]@pvals)
    #     familiesBonPVals[i]=calcBonPVals(familyArray[[i]]@pvals)
    if (echo=="ON")  print(c("Family Simes PVal:", familiesSimesPVals[i]))
  }
  
  #method 1 step 1:
  SelectedFamilies[[1]]<-SelectFamiliesBH(familiesSimesPVals,alpha=familySelectionAlpha)
  #method 2 step 1:
  #   SelectedFamilies[[2]]<-SelectFamiliesBH(familiesBonPVals)
  
  print (c("Selected Families:", SelectedFamilies[[1]]$ix))
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
                    #alpha=alpha*SelectedFamilies[[methodix]]$length/numOfFamilies)
                    alpha=insideFamilyAlpha)
          
        }
        print (c("Family:", SelectedFamilies[[1]]$ix[i], "Rejects:" , familyArray[[SelectedFamilies[[methodix]]$ix[i]]]@proceduresApplied[[1]]@rejects$ix))
        
      } #end for loop of i selected families     
    } #end for methodix
  }
  
}

calcSimesPVals<- function (PVals){
  min(p.adjust(PVals, method="BH"))
}
calcBonPVals<- function (PVals){
  min(p.adjust(PVals, method="bon"))
}

SelectFamiliesBH<-function (familiesSimesPVals,alpha=0.05){
  #b<-p.adjust(familiesSimesPVals,method="BH")
  #return (which(b<0.05))
  #print ("SELECT: ")
  #print (Preject(familiesSimesPVals))
  #print ("BH SELECT")
  return (Preject(method="BH", pvals=familiesSimesPVals, alpha=alpha))
  
}  
#return Array of Selected Families and old ixs
