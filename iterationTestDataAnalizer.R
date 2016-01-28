TukeyGraphs<- function (name="TukeyTest  280116-14-14 DesVector= 0,0,1,1,1 n= 10000 .xls",
                        withSD=FALSE, mindelta=0, maxdelta=4, interval=0.5, n=10000){
    #withSD means that the sd will appear in parenthesis
    #in the output file
    
    #load file
    wb1<- loadWorkbook(name)
    #getting names of sheets (deltas)
    sheets <- getSheets(wb1)
    print (sheets)
    METHOD_NAMES=rbind("QTukey Stat", "Pairwise", "Overall BH", "BH BH(Simes)")
    
    #init vars
    DATA<-NULL
    arrangedData<-NULL
    
    DATA<-readWorksheetFromFile(name,
                  sheet = sheets, 
                  startRow = 1, endRow = 5,
                  startCol = 2, endCol = 13)

    COL_NAMES= colnames(DATA[[sheets[1]]])
    # print (DATA[sheets[1]])
    print (length(sheets)/2)
    for (col_name in COL_NAMES){
        arrangedData[[col_name]]=matrix(NaN, 4, length(sheets)/2)
        colnames(arrangedData[[col_name]])<-seq(mindelta,maxdelta,interval)
        rownames(arrangedData[[col_name]])<-METHOD_NAMES
        for (delta in seq(mindelta,maxdelta,interval)){
           if (withSD==TRUE){
               arrangedData[[col_name]][,(delta+0.5)*2]=
               paste(
                   DATA[[paste("MEANS",delta)]][,col_name],
                   "(",
                   DATA[[paste("SD",delta)]][,col_name],
                   ")"
               )
           }
           else{
               arrangedData[[col_name]][,(delta+0.5)*2]<- 
                   DATA[[paste("MEANS",delta)]][,col_name]
           }
        }
            
#         
    }
    # print (arrangedData)

    ##exporting arranged table to a new Excel file
    wb2<-loadWorkbook(filename=paste("ARRANGED", name), create=TRUE)
                     
    for (col_name in COL_NAMES){
        createSheet(wb2, name = col_name)
        writeWorksheet(wb2, cbind(METHOD_NAMES,arrangedData[[col_name]]), sheet=col_name)
    }
    saveWorkbook(wb2)
    
    if (withSD){
        stop()
    }
    
    ##printing out graphs
    
    print ("ARR DATA")
    print (arrangedData[["OVR.POWER"]])
    print (arrangedData[["AVG.FDR"]])
    print (arrangedData[["OVR.FWER"]])
    plot (seq(mindelta,maxdelta,interval),arrangedData[["OVR.POWER"]][1,], 
          col="blue", ylim=c(0,1),xlim=c(mindelta,maxdelta),xlab= "delta", ylab= "Overall Power",
          main =c("Overall Power/delta in", n, "iterations with methods 1-4"), type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.POWER"]][2,], col="red", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.POWER"]][3,], col="green", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.POWER"]][4,], col="pink", type="o")
    
    legend("bottomright", 
           legend= c("Overall Power QTukey Stat ", "Overall Power Pairwise", 
                     "Overall Power Overall BH ","Overall Power BH BH(Simes) "), 
           lty=c(1,1,1,1),
           lwd=c(2,2,2,2),
           col=c("blue","red", "green", "pink")
    )
    dev.copy2pdf(file=paste(name, "POWER.pdf"))
    
    plot (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FDR"]][1,], 
          col="blue", ylim=c(0,max(arrangedData[["AVG.FDR"]][1,])+0.2),xlab= "delta", ylab= "FDR",
          main =c("FDRs/delta in", n, "iterations"), type="o",axes=TRUE)
    # axis(side=1, at=seq(mindelta,maxdelta,interval));
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FDR"]][1,], col="red", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FDR"]][2,], col="green", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FDR"]][2,], col="pink", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FDR"]][3,], col="brown", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FDR"]][4,], col="cyan", type="o")
    abline (h=0.05, col="black")
    
    legend("topleft", 
           legend= c("AVG FDR over the selected QTukey Stat ", "Overall FDR QTukey Stat ",
                     "AVG FDR over the selected Pairwise", "Overall FDR Pairwise", 
                     "Overall FDR Overall BH ",
                     "AVG FDR over the selected BH BH(Simes) ", "Overall FDR BH BH(Simes) " ), 
           lty=c(1,1,1,1),
           lwd=c(2,2,2,2),
           col=c("blue","red", "green", "pink", "brown","cyan")
    )
    dev.copy2pdf(file=paste(name, "FDR"))
    
    plot (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FWER"]][1,], 
          col="blue", ylim=c(0,1.5),xlab= "delta", ylab= "FWER",
          main =c("FWERs/delta in", n, "iterations"), type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FWER"]][1,], col="red", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FWER"]][2,], col="green", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FWER"]][2,], col="pink", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FWER"]][3,], col="brown", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FWER"]][4,], col="grey", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FWER"]][4,], col="cyan", type="o")
    abline (h=0.05, col="black") 
    legend("topleft", 
           legend= c("AVG FWER over the selected QTukey Stat ", "Overall FWER QTukey Stat ",
                     "AVG FWER over the selected Pairwise", "Overall FWER Pairwise", 
                     "Overall FWER Overall BH ", 
                     "AVG FWER over the selected BH BH(Simes) ","Overall FWER BH BH(Simes) "), 
           lty=c(1,1,1,1),
           lwd=c(2,2,2,2),
           col=c("blue","red", "green", "pink", "brown","grey", "cyan")
    )
    dev.copy2pdf(file=paste(name, "FWER.pdf"))
    
    plot (seq(mindelta,maxdelta,interval),arrangedData[["FR"]][1,], 
          col="blue", ylim=c(0,3),xlab= "delta", ylab= "E[V]",
          main =c("E[V]s/delta in", n, "iterations"), type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FR"]][1,], col="red", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["FR"]][2,], col="green", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FR"]][2,], col="pink", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FR"]][3,], col="brown", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["FR"]][4,], col="yellow", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["OVR.FR"]][4,], col="grey", type="o")
    abline (h=0.05, col="black")
    
    legend("topleft", 
           legend= c("AVG E[V] over the selected QTukey Stat ", "Overall E[V] QTukey Stat ", 
                     "AVG E[V] over the selected Pairwise", "Overall E[V] Pairwise", 
                     "Overall E[V] Overall BH ",
                     "AVG E[V] over the selected BH BH(Simes) ", "Overall E[V] BH BH(Simes) "), 
           lty=c(1,1,1,1),
           lwd=c(2,2,2,2),
           col=c("blue","red", "green", "pink", "brown","yellow","grey")
    )
    dev.copy2pdf(file=paste(name, "E[V].pdf"))
    
    plot (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FWER"]][1,] , col="blue", ylim=c(0,0.5),
          ylab="FWER/E[V]", xlab="delta",
          main =c("Avg FWER and Avg E[V]/delta in", n, "iterations with methods 1,2,4"), type="o")
    #    readline()
    lines (seq(mindelta,maxdelta,interval),arrangedData[["FR"]][1,], col="red", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FWER"]][2,], col="green", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["FR"]][2,], col="pink", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["AVG.FWER"]][4,], col="yellow", type="o")
    lines (seq(mindelta,maxdelta,interval),arrangedData[["FR"]][4,], col="grey", type="o")
    legend("topleft", 
           legend= c("Average FWER QTukey Stat ","E[V] over the selected QTukey Stat ",
                     "Average FWER Pairwise","E[V] over the selected Pairwise",
                     "Average FWER BH BH(Simes) ","E[V] over the selected BH BH(Simes) "), 
           lty=c(1,1,1,1),
           lwd=c(2,2,2,2),
           col=c("blue","red", "green", "pink","yellow","grey")
    )
    dev.copy2pdf(file=paste(name, "FWERandE[V].pdf"))  
    
}
