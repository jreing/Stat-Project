TukeyGraphs<- function (name="TukeyTest  270116-20-49 DesVector= -1,0,1 n= 10000 .xls")
    wb1<- loadWorkbook(name1)
    sheets <- getSheets(wb1)
    print (sheets)
    METHOD_NAMES=rbind("QTukey Stat", "Pairwise", "Overall BH", "BH BH(Simes)")
    mindelta=0
    maxdelta=4
    interval=0.5
    i=0
    
    for (sheet in sheets){
        DATA[i]<-readWorksheetFromFile(name1,
                                      sheet = sheet, 
                                      startRow = 1, endRow = 5,
                                      startCol = 2, endCol = 13)
    
        DATA[i]= matrix(as.numeric(unlist(DATA_1)),ncol=12, byrow = TRUE)
        i=i+1
    }
    # COL_NAMES= get colnames
    
    for (col_name in COL_NAMES){
        arrangedData[[col_name]]=matrix(NA, 4, length(COL_NAMES))
        colnames(arrangedData[[col_name]])=sheets
    }
    print(arranged)
    
        print (sheet)
        print (name)
        print (DATA)
      
    }
}
