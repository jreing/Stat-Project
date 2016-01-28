TukeyGraphs<- function (name="TukeyTest  280116-14-14 DesVector= 0,0,1,1,1 n= 10000 .xls"){
    wb1<- loadWorkbook(name)
    sheets <- getSheets(wb1)
    print (sheets)
    METHOD_NAMES=rbind("QTukey Stat", "Pairwise", "Overall BH", "BH BH(Simes)")
    mindelta=0
    maxdelta=4
    interval=0.5
    DATA<-NULL
    arrangedData<-NULL
    
    DATA<-readWorksheetFromFile(name,
                  sheet = sheets, 
                  startRow = 1, endRow = 5,
                  startCol = 2, endCol = 13)
    # print(sheets[1])
    COL_NAMES= colnames(DATA[[sheets[1]]])
    # print (DATA[sheets[1]])
    print (length(sheets)/2)
    for (col_name in COL_NAMES){
        arrangedData[[col_name]]=matrix(NaN, 4, length(sheets)/2)
        colnames(arrangedData[[col_name]])<-seq(mindelta,maxdelta,interval)
        rownames(arrangedData[[col_name]])<-METHOD_NAMES
        for (delta in seq(mindelta,maxdelta,interval)){
           arrangedData[[col_name]][,(delta+0.5)*2]<- 
               DATA[[paste("MEANS",delta)]][,col_name]
        }
            
#         
    }
    print (arrangedData)
      
 
}
