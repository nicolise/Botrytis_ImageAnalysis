library("EBImage")

label4Nicole <- function(orig.label , image.f = "R1_1205PM.JPG"){

  # generate label image
  MakeLabelMarker <- function(){
    png(filename = "lables.png",  width = 1000 , height = 1000 ,  units = "px")
    #par.default <- par()
    
    par(mfrow=c(20,10),oma = c(0,0,0,0),mar = c(0,0,0,0) , bg ="transparent" )
    for( i in 0:199 ){
      plot.new()
      plot.window(xlim=c(0,1), ylim=c(0,1)) 
      #plot(0:1,0:1,type="n")
      text(0.5,0.5,i,col='#ffff00' , cex= 6)
    }
    dev.off()
    #par(par.default) 
    readImage("lables.png" )
  }
  
  labelMarker <- MakeLabelMarker()
  
  getLabelMarker <- function(i){
    x <- i%%10
    xa <- (x * 100 +1) : ((x+1) * 100)
    y <- (i-x)/10
    ya <- (y * 50 +1) : ((y+1) * 50)
    labelMarker[xa,ya,]
  }

  RawImage <-readImage(image.f )
  #dx <- dim(RawImage)[1] 
  #dy <- dim(RawImage)[2]
  
  # start to draw lables
  for( i in 1:dim(orig.label)[1]) {
    row <- orig.label[i,]
    x <- round(row$Lesion.0.m.cx)
    y <- round(row$Lesion.0.m.cy)
    lablelName <- as.numeric(as.character(row$Object))
    print(paste("DEBUG" , Sys.time() , ": row ",i,",  label " , lablelName , ", Object " , row$Object) )
    label <- getLabelMarker(lablelName)
    
    xa <- (x - 50 + 1 ) : (x + 50 )
    ya <- (y - 50 + 1 ) : (y )
    RawImage[xa,ya,] <- (label[,,1:3] + RawImage[xa,ya,]) # C code would be better
    
  }
  print(paste("DEBUG" , Sys.time() , ": output image " , paste0(image.f,".labled.jpg") ) )
  writeImage(RawImage, paste0(image.f,".labled.jpg"))
}


# get lables 
orig.label <- read.csv("R1_1205PM_Results.csv" ,row.names = 1)
orig.label <- data.frame("Object" = rownames(orig.label),orig.label)
label4Nicole(orig.label,"R1_1205PM.JPG")


# test code 
labels <- read.csv("IMG_7758_Results.csv" )
label4Nicole(labels,"IMG_7758_LeafLabel.jpg")




