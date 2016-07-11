
setwd('/home/antoniofisk/Desktop/Uni/MasterThesis/Project/Tests done/Tests best irace')

par(mfrow=c(1,2))
feasibleSolutionCount <- NULL
Xist <- list(NULL)
Yist <- list(NULL)

for (k in 1:11){
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
    
    Yi[[j]] <- as.numeric(Yii)
    Xi[[j]] <- as.numeric(Xii)
    i <- i+1
    j <- j+1
  }
  
  Yi <- lapply(Yi, fun <- function(Yii){Yii[Yii<1.0]})
  for(t in 1:length(Xi)){
    Xi[[t]] <- Xi[[t]][(length(Xi[[t]])-length(Yi[[t]])+1):length(Xi[[t]])]
  }
  
  maxl <- max(unlist(lapply(Yi, length)))
  for(t in 1:length(Xi)){
    if(length(Xi[[t]])+1 < maxl){
      Xi[[t]] <- c(  Xi[[t]], rep(as.numeric(  max(unlist(lapply(Xi, max))) ) , maxl - (length(Xi[[t]])+1))  )
      Yi[[t]] <- c(Yi[[t]], rep( Yi[[t]][[length(Yi[[t]])]],  length(Xi[[t]]) - length(Yi[[t]]) )   )
    }
  }
  
  
  #####################################################
  Y.rtd <- list(NULL)
  X.rtd <- list(NULL)
  for(i in 1:length(Yi)){
    Xii <- Xi[[i]]
    Yii <- Yi[[i]]
    Y.complete <- NULL
    x.ind <- 2
    ind <- Xii[x.ind]
    for(j in 1:max(Xii)){
      if(j < ind){
        ;
      }else{
        if(x.ind < length(Xii))
          x.ind <- x.ind + 1
        ind <- Xii[[x.ind]]
      }
      Y.complete <- c(Y.complete, Yii[x.ind-1])
    }
    Y.rtd[[i]] <- as.numeric(Y.complete)
    X.rtd[[i]] <- 1:length(Y.rtd[[i]])
  }
  
  maxl <- max(as.numeric(lapply(Y.rtd, length)))
  for(i in 1:length(Y.rtd)){
    if(length(X.rtd[[i]])+1 < maxl){
      X.rtd[[i]] <- c(  X.rtd[[i]], (length(X.rtd[[i]])+1):maxl  )
      Y.rtd[[i]] <- c(Y.rtd[[i]], rep( Y.rtd[[i]][[length(Y.rtd[[i]])]],  length(X.rtd[[i]]) - length(Y.rtd[[i]]) )   )
    }
    #maxl - (length(Y.rtd[[i]]))
  }
  
  
  Y.rtd.dt <- Reduce(cbind, Y.rtd)
  X.rtd.mean <- X.rtd[[1]]
  Y.rtd.mean <- apply(Y.rtd.dt, 1, mean)
  
  
  filename.rtd <- paste(as.character(k), "RDT", gsub(as.character(k), "", filename), collapse="")
  
  jpeg(paste(filename.rtd, ".jpg"))
  
  cl <- rainbow(15)
  xlimit <- c(min(X.rtd.mean) ,  max(X.rtd.mean)  )
  ylimit <- c(min(Y.rtd.mean) ,  max(Y.rtd.mean)  )
  plot(c(0,1), c(0,1), xlim = xlimit, ylim = ylimit, main = paste(filename, ".jpg"), xlab = "#Iteration", ylab = "Objective value", type = "n")
  
  lines(X.rtd.mean, Y.rtd.mean, col = cl[i])
  
  
  dev.off()
  
  write(unlist(Y.rtd.mean), file = paste(filename.rtd, ".txt"),
        ncolumns = if(is.character(unlist(Y.rtd))) 1 else 1,
        append = FALSE, sep = " ")
  
  #########################################################
  
  
  #filename <- c(toString(k), filename)
  filename <- c(filename)
  filename <-paste(filename, collapse="")
  
  jpeg(filename)
  
  #Xi <- lapply(Xi, fun <- function(Xii){lapply(Xii, as.numeric)})

  cl <- rainbow(15)
  xlimit <- c(as.numeric(  min(unlist(lapply(Xi, min))) ) ,  as.numeric(  max(unlist(lapply(Xi, max))) )  )
  ylimit <- c(as.numeric(  min(unlist(lapply(Yi, min))) ) ,  as.numeric(  max(unlist(lapply(Yi, max))) )  )
  plot(c(0,1), c(0,1), xlim = xlimit, ylim = ylimit, main = paste(filename, ".jpg"), xlab = "#Solutions generated", ylab = "Objective value", type = "n")
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
  jpeg(paste(filename2, ".jpg"))
  boxplot(minObj, main = filename2, xlab = "#Seed", ylab = "Objective value")
  dev.off()
  
  
  minObj <- lapply(Yi, min)
  Yist[[k]] <- minObj 
  
}

filename3 <-paste("All instances Boxplot", filename, collapse="")

minObj <- lapply(Yist, as.numeric)
jpeg(paste(filename3, ".jpg"))
boxplot(minObj, main = filename2, xlab = "Instance", ylab = "Objective value")
dev.off()
