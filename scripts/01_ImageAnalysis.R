#Nicole E Soltis
#From Jason A Corwin
#14\06\10
#Leaf Image Analysis protocol
#-------------------------------------------------------------

#images are transformed into hue/saturation/value (hsv) color space
#leaf mask was created by selecting all pixels above Otsu’s saturation threshold 
#and below Otsu’s value threshold
#opening using a 15-pixel, disc-shaped kernel
#A mask marking the lesion within the leaf is created 
#by excluding all pixels that had a green or yellow hue within the leaf mask 
#opening using a 10-pixel, disc-shaped kernel
#masks are visually confirmed by overlaying the masks with the original image 
#using ImageJ
#all leaves and lesions are phenotyped using the masks and the original raw image 
#using computeFeatures.moment() and computeFeatures.shape() functions in ‘EBImage’ 
#which measures ~40 characteristics per object
#including area, perimeter, major and minor axes,  and eccentricity in pixels
#These thresholds and approaches will be tweaked for each different plant species 
#Leaf and lesion objects were also measured for 25 custom and composite phenotypes
#including the number of green and yellow pixels
#the proportion of yellow to green pixels or residual leaf area
#and the proportion of lesion area to leaf area
#These are selected to measure the senescence zone surrounding the lesion 

#-------------------------------------------------------------
rm(list=ls())
#learn current directory, and set it to the correct one
getwd()
setwd("C:/Users/nesoltis/Desktop/UGHP_Bcpics/03_Masked")

#only run this once to install image analysis
#source("http://bioconductor.org/biocLite.R")

#biocLite("EBImage")
#biocLite("CRImage")

library("EBImage")
library("CRImage")

#note: files already in .JPG format. Else change .JPG to .CR2 in script.
Files <- list.files(pattern=".JPG")
Files
for(i in 1:length(Files)) {
#tracks time to run script
print(Sys.time())

#reads image into a file for EBImage
Plate.cr2 <- readImage(Files[i])
#Plate.cr2 <- Plate.cr2[305:4929,497:2905,]

#retrieves (PIXEL?) dimensions of image file
Plate.row <- dim(Plate.cr2)[1]
Plate.col <- dim(Plate.cr2)[2]

# Import Image, get dimensions, and separate channels
# (what are [,,,])
Plate.red <- as.vector(Plate.cr2[,,1])
Plate.grn <- as.vector(Plate.cr2[,,2])
Plate.blu <- as.vector(Plate.cr2[,,3])

# Convert CR2 Image from rgb to hsv 
# The 14-bit pixel value range is 0:16383  ]
Plate.hsv <- rgb2hsv(Plate.red, Plate.grn, Plate.blu, maxColorValue = 16383)
Plate.hue <- Image(Plate.hsv[1,], dim = c(Plate.row, Plate.col))
Plate.sat <- Image(Plate.hsv[2,], dim = c(Plate.row, Plate.col))
Plate.value <- Image(Plate.hsv[3,]*10000, dim = c(Plate.row, Plate.col))
rm(Plate.red, Plate.grn, Plate.blu, Plate.hsv)

# Set histogram threshold and isolate leaves
#method to change grayscale to binary
#calculates optimum threshold for bi-modal histogram
#minimize intra-class variance 


#OtsuThresh <- calculateOtsu(as.vector(Plate.hue))
#Plate.mask <- Plate.hue>OtsuThresh
#display(Plate.mask)
#fills holes in objects
Lf.mask <- fillHull(Plate.hue > 0.15 & Plate.hue < 0.3)
#array-based structuring element
#for image filtering
#sigma sets SD of gaussian shape
kern <- makeBrush(5,shape = 'gaussian',sigma = 0.3)
#opening = erosion followed by dilation
#erode: applies mask to center over each foreground pixel
#pixels not covered by mask set to background
#dilate: applies mask to center over each background pixel
#pixels not covered by mask set to foreground
Lf.mask <- opening(Lf.mask,kern)

#adds labels to objects
Obj <- bwlabel(Lf.mask)
Obj.minRad <- computeFeatures.shape(Obj)[,'s.radius.min']
Obj.lvs <- as.numeric(names(Obj.minRad[Obj.minRad > 10]))
Lf.mask <- Image(Obj %in% Obj.lvs, dim = c(Plate.row, Plate.col))


#OtsuThresh <- calculateOtsu(as.vector(Plate.sat))
#Plate.mask <- fillHull(Plate.sat>OtsuThresh)
#display(Plate.mask)

#OtsuThresh <- calculateOtsu(as.vector(Plate.value))
#Plate.mask1 <- Plate.value>OtsuThresh
#display(Plate.mask1)

#OtsuThresh <- calculateOtsu(as.vector(Plate.value))
#Plate.mask2 <- !Plate.value>OtsuThresh
#Plate.mask <- Plate.mask | Plate.mask2


#### Remove white pixels
#Plate.mask2 <- Plate.sat>0.0000000001
#Plate.mask <- Plate.mask == Plate.mask2
#rm(Plate.mask2)

# Fill holes
#Plate.mask <- fillHull(!Plate.mask)

# Reduce background and isolate leaf areas
#kern <- makeBrush(15,shape='disc')
#Leaf.Area.Mask <- opening(Plate.mask,kern)
writeImage(Lf.mask, gsub(".JPG", "_LeafMask.JPG", Files[i]))

#Write a labeled image
# Leaf.label <- bwlabel(Lf.mask)
# xy = computeFeatures.moment(Leaf.label)[, c('m.cx', 'm.cy')]
# font = drawfont(weight=600, size=46)
# Leaf.label <- paintObjects(Leaf.label, Plate.cr2, col='#ffff00')
# Leaf.label <- drawtext(Leaf.label, xy=xy, labels=as.character(1:nrow(xy)), font=font, col="yellow")
# writeImage(Leaf.label, gsub(".CR2", "LeafLabel.JPG", Files[i]))
# rm(Leaf.label)

# Identify lesions
#Try modifications here!
Lesion.mask <- fillHull((Plate.sat<0.38 * Lf.mask))

### Make high intensity pixels in lesion low, but not 0
#Lesion.mask[Lesion.mask>0.6] <- 0.1
### Initial Threshold
#OtsuThresh <- calculateOtsu(Lesion.mask[Lesion.mask>0])

#Make Yellow Pixel Mask
###This portion finds the high intensity Green and Red pixel and combines to find yellows
###The high intensity blues are filtered out to remove white from the mask
#Otsu.Red <- calculateOtsu(as.vector(Plate.cr2[,,1])[Plate.cr2[,,1]>0.5])  
#Otsu.Grn <- calculateOtsu(as.vector(Plate.cr2[,,2])[Plate.cr2[,,2]>0.5])
#Otsu.Blu <- calculateOtsu(as.vector(Plate.cr2[,,3][Plate.cr2[,,3]>0.5 & Plate.cr2[,,3]<1]))
#Plate.Ylw <- Plate.cr2[,,1]>Otsu.Red & Plate.cr2[,,2]>Otsu.Grn & Plate.cr2[,,3]<Otsu.Blu

### Need intersection of high intensity Green and Red pixels 
#Lesion.mask <- Lesion.mask<OtsuThresh & Lesion.mask>0 & !Plate.Ylw
#Lesion.mask <- Lesion.mask & Leaf.Area.Mask

### Fill Holes, Filter within Leaves and Open 
#Lesion.mask <- fillHull(Lesion.mask)
#Lesion.mask <- Lesion.mask & Leaf.Area.Mask
kern2 <- makeBrush(9,shape='disc') #changed from 8: next odd number
Lesion.mask <- opening(Lesion.mask,kern2)
writeImage(Lesion.mask, gsub(".JPG", "_LesionMask.JPG", Files[i]))

}


