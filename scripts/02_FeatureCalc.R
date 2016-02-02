rm(list=ls())
setwd("~/Projects/Botrytis_ImageAnalysis/photos/examples")

#install EBImage and CRImage (only need to do this once per R installation)
#source("https://bioconductor.org/biocLite.R")
#biocLite("EBImage")
#biocLite("CRImage")
library("CRImage")
#inputs = .JPG raw file, .tif masks

#--------------------------------------------------------------------
#WRITE FUNCTION FOR LABELS (FROM VIVIAN)
label4Nicole <- function(orig.label , image.f = RawFiles[f]){

  # generate label image
  MakeLabelMarker <- function(){
    png(filename = (gsub(".JPG","_labels.png",RawFiles[f])),  width = 1000 , height = 1000 ,  units = "px")
    #par.default <- par()
    
    par(mfrow=c(20,10),oma = c(0,0,0,0),mar = c(0,0,0,0) , bg ="transparent" )
    for( i in 0:199 ){
      plot.new()
      plot.window(xlim=c(0,1), ylim=c(0,1)) 
      #plot(0:1,0:1,type="n")
      text(0.5,0.5,i,col='#000000' , cex= 6)
    }
    dev.off()
    #par(par.default) 
    readImage(gsub(".JPG","_labels.png",RawFiles[f]))
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
  
  # start to draw lables
  for( i in 1:dim(orig.label)[1]) {
    row <- orig.label[i,]
    x <- round(row$Lesion.0.m.cx)
    y <- round(row$Lesion.0.m.cy)
    labelName <- as.numeric(as.character(row$Object))
    print(paste("DEBUG" , Sys.time() , ": row ",i,",  label " , labelName , ", Object " , row$Object) )
    label <- getLabelMarker(labelName)
    
    xa <- (x - 50 + 1 ) : (x + 50 )
    ya <- (y - 50 + 1 ) : (y )
    RawImage[xa,ya,] <- (label[,,1:3] + RawImage[xa,ya,]) # C code would be better
    
  }
  print(paste("DEBUG" , Sys.time() , ": output image " , paste0(image.f,"_labled.jpg") ) )
  writeImage(RawImage, paste0(image.f,"_labled.jpg"))
}

#----------------------------------------------------------------------


RawFiles <- list.files(pattern = ".JPG")
LeafMaskFiles <- list.files(pattern= "LeafMask.tif")
LesionMaskFiles <- list.files(pattern="LesionMask.tif")

if(length(RawFiles) == length(LeafMaskFiles) & length(RawFiles) == length(LesionMaskFiles)) {
  for(f in 1:length(RawFiles)) {
    
    
    print(Sys.time())
    print(RawFiles[f])
    print("Reading Files....")
    RawImage <- readImage(RawFiles[f])
    LeafMask <- readImage(LeafMaskFiles[f])
    LeafMask[LeafMask > 0.5] <- 1
    LeafMask[LeafMask < 0.5] <- 0
    LesionMask <- readImage(LesionMaskFiles[f])
    LesionMask[LesionMask > 0.5] <- 1
    LesionMask[LesionMask < 0.5] <- 0
    LeafResidual <- LeafMask & !LesionMask
    
    Image.row <- dim(RawImage)[1]
    Image.col <- dim(RawImage)[2]
    
    print("Calculating hue....")
    ImageHSV <- rgb2hsv(as.vector(RawImage[,,1]), as.vector(RawImage[,,2]), as.vector(RawImage[,,3]), maxColorValue = 16383)
    Image.hue <- Image(ImageHSV[1,], dim = c(Image.row, Image.col))
    #Image.sat <- Image(ImageHSV[2,], dim = c(Image.row, Image.col))
    #Image.val <- Image(ImageHSV[3,]*10000, dim = c(Image.row, Image.col))
    rm(ImageHSV)
    
    #Yellow and green pix
    print("Identifying Yellow and Green pixels....")
    Yellow.pix <- Image.hue > 0.145 & Image.hue < 0.188 & LeafMask > 0.5
    Green.pix <- Image.hue > 0.188 & Image.hue < 0.4167 & LeafMask > 0.5
    
    ###Label leaves and cycle through measurements
    print("Calculating standard leaf metrics....")
    LeafMask.lab <- bwlabel(LeafMask)
    Leaf.Results <- computeFeatures(LeafMask.lab,RawImage,methods.noref=c("computeFeatures.moment", "computeFeatures.shape"), methods.ref=c("computeFeatures.moment"), xname = "Leaf",refnames=c("red","grn","blu"))
    LesionMask.lab <- LesionMask * LeafMask.lab
    Lesion.Results <- computeFeatures(LesionMask.lab,RawImage,methods.noref=c("computeFeatures.moment", "computeFeatures.shape"), methods.ref=c("computeFeatures.moment"), xname = "Lesion",refnames=c("red","grn","blu"))
    
    if(dim(Lesion.Results)[1] != dim(Leaf.Results)[1]) {
      Leaf.Results <- Leaf.Results[rownames(Leaf.Results) %in% rownames(Lesion.Results),]  
    }
    
    Results <- cbind(Lesion.Results, Leaf.Results)
    rm(Lesion.Results, Leaf.Results)
    
    print("Calculating non-standard leaf metrics....")
    MoreResults <- data.frame(NA,NA,NA,NA,NA,NA,NA)
    for(i in 1:max(LeafMask.lab)) {
      
      Lesion.Size <- sum(LesionMask.lab == i)
      Lesion.Grn <- sum(LesionMask.lab == i & Green.pix)
      Lesion.Ylw <- sum(LesionMask.lab == i & Yellow.pix)
      
      Leaf.Size <- sum(LeafMask.lab == i)
      Leaf.Grn <- sum(LeafMask.lab == i & Green.pix & !LesionMask)
      Leaf.Ylw <- sum(LeafMask.lab == i & Yellow.pix  & !LesionMask)
      
      Lesion.Prop <- Lesion.Size/Leaf.Size
      
      MoreResults[i,] <- c(Lesion.Size, Leaf.Size, Lesion.Prop, Leaf.Grn, Leaf.Ylw, Lesion.Grn, Lesion.Ylw)
      
    }
    
    colnames(MoreResults) <- c("Lesion.Size", "Leaf.Size", "Lesion.Prop", "Leaf.Grn", "Leaf.Ylw", "Lesion.Grn", "Lesion.Ylw")
    
    if(dim(Results)[1] != dim(MoreResults)[1]) {
      MoreResults <- MoreResults[rownames(MoreResults) %in% rownames(Results),]  
    }
    
    Results <- cbind(Results, MoreResults)
    
    
    
    
    
    ###Reorder Leaf Labels and Results
    print("Reordering leaf results for silly humans....")
    # find positions of objects and calculate distance matrix
    xy = computeFeatures.moment(LeafMask.lab)[, c('m.cx', 'm.cy')]
    Leaf.pos <- cbind(c(1,as.data.frame(xy)$m.cx), c(1,as.data.frame(xy)$m.cy))
    Temp.dist <- as.matrix(dist(Leaf.pos))
    #Gets the leaf closest to the top left and starts the key
    LeafOrder <- match(min(Temp.dist[-1,1]),Temp.dist[-1,1])
    #removes the reference top-left pixel
    Leaf.pos <- cbind(as.data.frame(xy)$m.cx, as.data.frame(xy)$m.cy)
    #recalculates distance matrix for objects
    Temp.dist <- as.matrix(dist(Leaf.pos))
    
    RowKey <- 0
    #loop through objects to reorder
    for(i in 1:dim(Temp.dist)[2]) {
      
      #Finds the ratio of x to y distance
        #if high, then objects are in the same line.  If low, objects are in diff lines.
      x.dif <- Leaf.pos[,1]-Leaf.pos[LeafOrder[i],1]
      y.dif <- Leaf.pos[,2]-Leaf.pos[LeafOrder[i],2]
      xy.ratio <- abs(x.dif/y.dif)
      
        # Tests if it has finished ordering all leaf objects
      if(length(LeafOrder) == dim(Leaf.pos)[1]) {
        break
        # Test if Image is done and missing values
      } else if(is.na(sum(xy.ratio > 1 & x.dif > 0) == 0)) {
        break
        # Test if it is at the end of the row
      } else if(sum(xy.ratio > 1 & x.dif > 0) == 0) {
        
        FirstLeaf <- LeafOrder[length(LeafOrder)-RowKey]
        
        # Find distances and ratio for all Leaf Objects relative to reference leaf 'FirstLeaf'.
        x.dif <- Leaf.pos[,1]-Leaf.pos[FirstLeaf,1]
        y.dif <- Leaf.pos[,2]-Leaf.pos[FirstLeaf,2]
        xy.ratio <- abs(x.dif/y.dif)
        
        # Find the leaf object with the minimum distance from the reference leaf object that is in the same row.
        #Leaf <- match(min((Temp.dist[,FirstLeaf][xy.ratio < 1 & y.dif >0]), na.rm = TRUE), Temp.dist[,FirstLeaf])
        Leaf <- match(min((Temp.dist[,FirstLeaf][xy.ratio < 1]), na.rm = TRUE), Temp.dist[,FirstLeaf])
        LeafOrder <- append(LeafOrder,Leaf)
        
        RowKey <- 0
        
        # If in the same row and 
      } else {
        Leaf <- match(min(Temp.dist[,LeafOrder[i]][xy.ratio > 1 & x.dif>0],na.rm = TRUE), Temp.dist[,LeafOrder[i]])
        LeafOrder <- append(LeafOrder,Leaf)
        
        RowKey <- RowKey+1
      }

    }
    
    LeafOrder <- unique(LeafOrder[!is.na(LeafOrder)])
    NewOrder <- 1:dim(xy)[1]
    MissingVal <- NewOrder[!NewOrder %in% LeafOrder]
    LeafOrder <- c(LeafOrder,MissingVal)
    Key <- cbind(LeafOrder,NewOrder)
    
    print("Reordering leaf mask objects....")
    for(i in 1:dim(Key)[2]){
      LeafMask.lab[LeafMask.lab == Key[i,1]] <- Key[i,2]  
      
    }
 #   rm(Plate.hue, Plate.sat, Plate.value, Yellow.pix, Obj)
    gc()
    
    ###Print labeled image....
    print("Making reference image....")
    xy <- xy[Key[,1],]
    rownames(xy) <- 1:dim(xy)[1]
    Leaf.label <- paintObjects(LeafMask.lab, RawImage, col='#ffff00')
    gc()
#FIX HERE???
    writeImage(Leaf.label, paste(gsub(".JPG","_LeafLabelTEST.jpg",RawFiles[f])))
    #Leaf.label <- bwlabel(Leaf.label) #it appears there is no way to still mark images with labels
    rm(Leaf.label)
    
    ###Print results for image....
    print("Printing results....")
    Filtered <- Key[,1] %in% as.numeric(rownames(Results))
    Filtered2 <- match(LeafOrder[Filtered], rownames(Results))
    Results <- Results[Filtered2,]
    rownames(Results) <- Key[,2][Filtered]
    Results <- cbind(RawFiles[f], Results)
    colnames(Results)[1] <- "Image"
    write.csv(Results, paste(gsub(".JPG","_Results.csv",RawFiles[f])))
    
    
  }
  
}

#---------------------------------------------------------------------
#PRINT IMAGE LABELS FROM VIVIAN SCRIPT

RawFiles <- list.files(pattern = ".JPG")
ResultsFiles <- list.files(pattern="Results.csv")

#RawFiles <- "208b_08b.JPG"
#ResultsFiles <- "208b_08b_Results.csv"

if(length(RawFiles) == length(ResultsFiles)){
  for(f in 1:length(RawFiles)) {
     
    print(Sys.time())
    print(RawFiles[f])


# get labels 
orig.label <- read.csv(ResultsFiles[f] ,row.names = 1)
orig.label <- data.frame("Object" = rownames(orig.label),orig.label)
label4Nicole(orig.label,RawFiles[f])

# test code 
labels <- read.csv(ResultsFiles[f] )
label4Nicole(labels,gsub(".JPG","_LeafLabel.jpg",RawFiles[f]))
	}
}


print("All done meat wad.  This had better go into Science or Nature for all my hard work!")
