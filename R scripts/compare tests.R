if (!require("data.table")) install.packages("data.table")

main.dir <- '/home/antoniofisk/Desktop/Uni/MasterThesis/Project/tests_processor' 

setwd(main.dir)


dir.tests <- list.dirs('.', recursive=FALSE)

par(mfrow=c(1,2))
feasibleSolutionCount <- NULL
Xist <- list(NULL)
Yist <- list(NULL)

instances <- 1:11#c(1,3,6,11)
Y.tests <- list(NULL)
for (dir in dir.tests){
  #  i <-8
  setwd(file.path(main.dir, dir))
  files <- list(NULL)
#   files <- list.files(path = ".", pattern = " RDT .*.txt", all.files = FALSE,
#                       full.names = FALSE, recursive = FALSE,
#                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

  kk <- 1
  for(k in instances){
    files[[kk]] <-  list.files(path = ".", pattern = paste(as.character(k), "RDT .*.txt"), all.files = FALSE,
                                               full.names = FALSE, recursive = FALSE,
                                               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    if(length(files[[kk]]) > 1)
    files[[kk]] <- files[[kk]][[length(files[[kk]])]]
    
    kk <- kk + 1
  }
  

  files <- Reduce(fun <- function(e1, e2) {
    if (length(e2) == 0) 
      {e1} 
    else
      {
        e1[[length(e1)+1]] <- e2
        return (e1)
      }
    }, 
    files)
  Y <- list(NULL)
  Y <- lapply(files, read.csv, head=TRUE, sep=" ")
#   Y <- Reduce(fun <- function(y1, y2) {
#     y1[[length(y1)+1]] <- y2
#     return(y1)
#     }, Y)
  Y.tests[[length(Y.tests)+1]] <- Y
  
}

Y.tests <- Y.tests[2:length(Y.tests)]

setwd(main.dir)

Yt <- list(NULL)
for(i in 1:min(unlist(lapply(Y.tests, length)))){
  Ytt <- list(NULL)
  for(j in 1:length(Y.tests))
    Ytt[[j]] <- Y.tests[[j]][[i]]
  Yt[[i]] <- Ytt
}


#Xi <- lapply(Xi, fun <- function(Xii){lapply(Xii, as.numeric)})
kk <- 1
for(i in 1:length(Yt)){
  
  Yi <- Yt[[i]]
  Xi <- lapply(Yi, fun <- function(y) {1:length(unlist(y))})
  
  jpeg(paste(as.character(instances[kk]), "Comparision.jpg"))
  
  cl <- rainbow(length(Yi))
  xlimit <- c(as.numeric(  min(unlist(lapply(Xi, min))) ) ,  as.numeric(  max(unlist(lapply(Xi, max))) )  )
  ylimit <- c(as.numeric(  min(unlist(lapply(Yi, min))) ) ,  as.numeric(  max(unlist(lapply(Yi, max))) )  )
  plot(c(0,1), c(0,1), xlim = xlimit, ylim = ylimit, main = paste("Instance", as.character(instances[kk])), xlab = "#Solutions generated", ylab = "Objective value", type = "n")
  legend('topright', c("Arbitrary configuration", "5 mIrace selected configuration", "15 min Irace selected configuration") , 
         lty=1, col=cl, bty='n', cex=.75)
  for(j in 1:length(Yi)){
    

    X <- Xi[[j]]#c(Y[,2], solutionCount)
    Y <- unlist(Yi[[j]])#Y[,1]
    #  Y <- c(Y, Y[length(Y)])
    
    #  points(Xi[[i]], Yi[[i]], col = cl[i])
    lines(X, Y, col = cl[j])
    

    
  }
  

  
  dev.off()

  
#   maxObj <- lapply(Yi, as.numeric)
#   minObj <- lapply(Yi, as.numeric)
#   jpeg(paste(as.character(instances[kk]), "Boxplot.jpg"))
#   boxplot(minObj, main = filename2, xlab = "#Seed", ylab = "Objective value")
#   dev.off()
  

  
  kk <- kk + 1
}


setwd(main.dir)

kk <- 0
for (dir in dir.tests){

  setwd(dir)
  files <- list.files(path = ".", pattern = paste("Objective$"), all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  files <- files[c(3:11, 1, 2)]

  for(file in files){
      #  i <-8
      filename <- NULL
      filename <- "Objective"
      filename <- c(toString(k), filename)
      filename <- paste(filename, collapse="")
      filename <- file
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
  }
  
  Yb[[k]] <- Yi
  
  setwd(main.dir)
  
  filename <- paste(paste(file, "boxplot") , ".jpg")
  jpeg(filename)
  Yb <- unlist(lapply(Y.tests[[1]], fun <- function(y){y[nrow(y), ]}))
  boxplot(Yi, main = filename, xlab = "#Seed", ylab = "Objective value")
  dev.off()

}

setwd(main.dir)
for(k in 1:11){
  filename <- paste(k, "Objective", sep = "")
  
  Yb <- list(NULL)
  kk <- 1
  for(dir in dir.tests){
    setwd(dir)
    Y <- as.matrix(read.csv(file=filename, head=TRUE, sep=" "))
#     Y <- Y[, 1]
    
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
      
      Yi[[j]] <- min(as.numeric(Yii))
      Xi[[j]] <- as.numeric(Xii)
      i <- i+1
      j <- j+1
    }
    
    Yb[[kk]] <- Yi
    kk <- kk + 1
    setwd(main.dir)
  }  
  
  filename <- paste(paste(instances[[k]], "boxplot") , ".jpg")
  jpeg(filename)
  boxplot(lapply(Yb, unlist), main = filename, xlab = "#Seed", ylab = "Objective value")
  dev.off()
  
}

