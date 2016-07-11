
setwd('/home/antoniofisk/Desktop/Uni/MasterThesis/Project/Tests')

par(mfrow=c(1,2))
feasibleSolutionCount <- NULL
Xist <- list(NULL)
Yist <- list(NULL)

for (k in c(1,2,3,4)){
#  i <-8
  filename <- NULL
  filename <- "Objective"
  filename <- c(toString(k), filename)
  filename <- paste(filename, collapse="")
  Y <- as.matrix(read.csv(file=filename, head=TRUE, sep=" "))
  

  Yi <- list(NULL)
  Xi <- list(NULL)
  i <- 1
  j <- 1
  while(i < nrow(Y) ){
    Yii <- NULL
    Xii <- NULL
    while(!is.na(Y[i, 2]) && i < nrow(Y) ){
      Yii <- c(Yii, Y[i, 1])
      Xii <- c(Xii, Y[i, 2])
      i <- i + 1
    }
   
    Yi[[j]] <- Yii
    Xi[[j]] <- Xii
    i <- i+1
    j <- j+1
  }



#filename <- c(toString(k), filename)
  filename <- c(filename, ".jpg")
  filename <-paste(filename, collapse="")
  
  jpeg(filename)
  
  cl <- rainbow(15)
  xlimit <- c(as.numeric(  min(unlist(lapply(Xi, min))) ) ,  as.numeric(  max(unlist(lapply(Xi, max))) )  )
  ylimit <- c(as.numeric(  min(unlist(lapply(Yi, min))) ) ,  as.numeric(  max(unlist(lapply(Yi, max))) )  )
  plot(c(0,1), c(0,1), xlim = xlimit, ylim = ylimit, main = filename, xlab = "#Solutions generated", ylab = "Objective value", type = "n")
  for(i in 1:length(Yi)){
    X <- Xi[[i]]#c(Y[,2], solutionCount)
    Y <- Yi[[i]]#Y[,1]
    #  Y <- c(Y, Y[length(Y)])
  
  #  points(Xi[[i]], Yi[[i]], col = cl[i])
    lines(Xi[[i]], Yi[[i]], col = cl[i])
    
  }
  
  dev.off()
  
  filename2 <-paste("Boxplot", filename, collapse="")
  
  maxObj <- lapply(Yi, as.numeric)
  minObj <- lapply(Yi, as.numeric)
  jpeg(filename2)
  boxplot(minObj, main = filename2, xlab = "#Seed", ylab = "Objective value")
  dev.off()


  minObj <- lapply(Yi, min)
  Yist[[k]] <- minObj 

}

filename3 <-paste("All instances Boxplot", filename, collapse="")

minObj <- lapply(Yist, as.numeric)
jpeg(filename3)
boxplot(minObj, main = filename2, xlab = "Instance", ylab = "Objective value")
dev.off()


