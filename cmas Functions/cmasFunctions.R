#### Functions for polar transformations


# Helper function for fpCenterMat function defined below (please excuse the
# awful naming scheme). Performs the firing pin estimation using the label
# function in imager
fpCenterCalc <- function(x3p,
                         scheme = 3,
                         high_connectivity = FALSE,
                         tolerance = 0,
                         dilate_prop = .1,
                         centerOffset = c(0,0)){

  mat_bfRegion <- x3p$surface.matrix

  mat_bfRegionBinarized <- mat_bfRegion
  mat_bfRegionBinarized[!is.na(mat_bfRegionBinarized)] <- 1
  mat_bfRegionBinarized[is.na(mat_bfRegionBinarized)] <- 0

  #Label the different regions of the scan using the edges as borders
  mat_bfRegionLabeled <- mat_bfRegionBinarized %>%
    imager::as.cimg() %>%
    imager::dilate_square(size = dilate_prop*ncol(mat_bfRegion)) %>%
    imager::imgradient(scheme = 3) %>%
    imager::enorm() %>%
    imager::add() %>%
    imager::label(high_connectivity = high_connectivity,
                  tolerance = tolerance) %>%
    as.matrix()

  #The pixel in the center of the image should be a part of the firing pin
  #impression hole
  mat_bfRegioncenterLabel <- mat_bfRegionLabeled[round(nrow(mat_bfRegionLabeled)/2),round(ncol(mat_bfRegionLabeled)/2)]

  #Identify all pixels that share a label with the center pixel (these are
  #assumed to be only pixels that are a part of the firing pin impression)
  mat_bfRegionfpHoleIndices <- which(mat_bfRegionLabeled == mat_bfRegioncenterLabel,arr.ind = TRUE)

  #The center pixel of this region is assumed to be the center of the firing pin
  #impression hole
  mat_bfRegionfpHoleCenter <- round(colMeans(mat_bfRegionfpHoleIndices)) + centerOffset

  return(mat_bfRegionfpHoleCenter)
}

#Pads a surface matrix with NA rows/cols based on estimated firing pin center
fpCenterMat <- function(x3p,
                        scheme = 3,
                        high_connectivity = FALSE,
                        tolerance = 0,
                        dilate_prop = .1,
                        centerOffset = c(0,0)){

  mat <- x3p$surface.matrix

  fpCenter <- fpCenterCalc(x3p,
                           scheme = scheme,
                           high_connectivity = high_connectivity,
                           tolerance = tolerance,
                           centerOffset = centerOffset)

  dimToPad <- 2*round(fpCenter - dim(mat)/2)

  rowPad <- matrix(NA,nrow = abs(dimToPad[1]),ncol = ncol(mat))

  #if fp center is below (larger in row index, starting from the top left corner),
  #then we want to pad the first few rows of the matrix
  if(sign(dimToPad[1]) == 1){
    mat1_padded <- rbind(mat,rowPad)
  }
  else{
    mat1_padded <- rbind(rowPad,mat)
  }

  colPad <- matrix(NA,nrow = nrow(mat1_padded),ncol = abs(dimToPad[2]))

  #if fp center is to the right of (larger in col index, starting from the top
  #left corner), then we want to pad the first few cols of the matrix
  if(sign(dimToPad[2]) == 1){
    mat1_padded <- cbind(mat1_padded,colPad)
  }
  else{
    mat1_padded <- cbind(colPad,mat1_padded)
  }
  #update metainformation in x3p too
  x3p$surface.matrix <- mat1_padded
  x3p$header.info$sizeX <- nrow(mat1_padded)
  x3p$header.info$sizeY <- ncol(mat1_padded)
  return(x3p)
}

#adapted from wvtools package. rotationResolution is the change in theta value
#between pixels (e.g., rotationResolution = .5 means output matrix will be
#360/.5 = 720 elements wide). Could perhaps do the same with radius resolution?
#cropWhitespace means that the NA rows/columns around the transformed breech
#face surface will be removed.
car2pol <- function (x3p, method = "nearest",rotationResolution = .5,cropWhitespace = TRUE){

  i.img <- x3p$surface.matrix

  if (method == "nearest") {
    p <- nrow(i.img)
    q <- ncol(i.img)
    m <- max(p, q)
    conv.img <- matrix(data = NA,
                       nrow = (trunc(m/2) - 1), #since we're unwrapping relative the center of the matrix, we want the height of the unwrapped image to be equal to have of the matrix dimension
                       ncol = length(seq(0,360,by = rotationResolution) + rotationResolution))

    for (th in seq(0,360,by = rotationResolution) + rotationResolution){
      for (r in 1:(trunc(m/2) - 1)) {
        th.rad <- pi * (th - 1)/180

        if(any(dim(i.img) < c(round(p/2 + (r - 1) * cos(th.rad)) - 1,
                              round(q/2 + (r - 1) * sin(th.rad))))){
          conv.img[r, th*as.integer(1/rotationResolution) - 1] <- NA
        }
        else if(length(i.img[round(p/2 + (r - 1) * cos(th.rad)) - 1,
                             round(q/2 + (r - 1) * sin(th.rad))]) != 1){
          conv.img[r, th*as.integer(1/rotationResolution) - 1] <- NA
        }

        else{
          conv.img[r, th*as.integer(1/rotationResolution) - 1] <- i.img[round(p/2 + (r - 1) * cos(th.rad)) - 1,
                                                                        round(q/2 + (r - 1) * sin(th.rad))]
        }
      }
    }

  }
  else if (method == "bilinear") {
    p <- nrow(i.img)
    q <- ncol(i.img)
    m <- max(p, q)
    xc <- p%/%2
    yc <- q%/%2
    conv.img <- matrix(data = NA,
                       nrow = (trunc(m/2) - 1),
                       ncol = length(seq(0,360,by = rotationResolution) + rotationResolution))

    for (th in seq(0,360,by = rotationResolution) + rotationResolution) {
      for (r in 1:(trunc(m/2) - 1)) {
        th.rad <- pi * (th - 1)/180
        Sxy <- c(xc + (r - 1) * cos(th.rad), yc +
                   (r - 1) * sin(th.rad))
        xy <- trunc(Sxy)
        dxy <- Sxy - xy

        if(any(nrow(i.img) < c(xy[1],xy[1] + 1)) |
           any(ncol(i.img) < c(xy[2],xy[2] + 1))){
          conv.img[r, th*as.integer(1/rotationResolution) - 1] <- NA
        }
        else{
          a <- matrix(c(i.img[xy[1], xy[2]], i.img[xy[1] + 1, xy[2]],
                        i.img[xy[1], xy[2] + 1], i.img[xy[1] + 1, xy[2] + 1]), 2, 2)

          b <- matrix(c(1 - dxy[1], dxy[1]), 1, 2)

          c <- matrix(c(1 - dxy[2], dxy[2]), 2, 1)

          conv.img[r, th*as.integer(1/rotationResolution) - 1] <- b %*% a %*% c
        }
      }
    }
  }

  #Remove rows/columns containing all 0s (which seems to often be the last
  #column of the matrix, perhaps due to some weirdness with interpolating the
  #end of the image.)
  zerosPerRow <- rowSums({conv.img == 0},na.rm = TRUE)

  conv.img <- conv.img[zerosPerRow < ncol(conv.img),]

  zerosPerCol <- colSums({conv.img == 0},na.rm = TRUE)

  conv.img <- conv.img[,zerosPerCol < nrow(conv.img)]

  if(cropWhitespace){
    #Remove rows/columns containing only NAs
    naPerRow <- conv.img %>%
      is.na() %>%
      rowSums()

    conv.img <- conv.img[naPerRow < ncol(conv.img),]

    naPerCol <- conv.img %>%
      is.na() %>%
      colSums()

    conv.img <- conv.img[,naPerCol < nrow(conv.img)]
  }

  x3p$surface.matrix <- conv.img
  x3p$header.info$sizeX <- nrow(conv.img)
  x3p$header.info$sizeY <- ncol(conv.img)

  return(x3p)
}

polarEstimateRotation <- function(reference,
                                  target,
                                  polarInterpolation = "nearest",
                                  rotationResolution = .5,
                                  rotationInterpolation = 0,
                                  scheme = 3,
                                  high_connectivity = FALSE,
                                  tolerance = 0){

  #transform both surface matrices into polar coordinates after centering on
  #firing pin hole by interpolation
  reference_polar <- reference %>%
    fpCenterMat(scheme = scheme,
                high_connectivity = high_connectivity,
                tolerance = tolerance) %>%
    car2pol(method = polarInterpolation,
            rotationResolution = rotationResolution,
            cropWhitespace = FALSE)

  target_polar <- target %>%
    fpCenterMat(scheme = scheme,
                high_connectivity = high_connectivity,
                tolerance = tolerance) %>%
    car2pol(method = polarInterpolation,
            rotationResolution = rotationResolution,
            cropWhitespace = FALSE)

  #Finding optimal rotation via FFT-based CCF requires replacing missing values
  reference_polar$surface.matrix[is.na(reference_polar$surface.matrix)] <- mean(reference_polar$surface.matrix,na.rm = TRUE)
  target_polar$surface.matrix[is.na(target_polar$surface.matrix)] <- mean(target_polar$surface.matrix,na.rm = TRUE)

  ccfMat <- cmcR:::filterViaFFT(reference_polar$surface.matrix,target_polar$surface.matrix)
  # ccfMat <- imager::correlate(im = as.cimg(reference_polar$surface.matrix),
  #                             filter = as.cimg(target_polar$surface.matrix)) %>%
  #   as.matrix()

  #Find index at which maximum CCF occurs. The corresponding column represents
  #the rotation at which the two scans best align. Need to multiple by the
  #rotation resolution to get back the correct rotation value
  alignment <- (which(ccfMat == max(ccfMat),arr.ind = TRUE) - dim(ccfMat)/2)*rotationResolution


  return(alignment[2])
}

#Rotates a target scan to align with a reference scan by estimating the rotation
#in the polar domain
preProcess_rotateScan <- function(reference,
                                  target,
                                  polarInterpolation = "nearest",
                                  rotationResolution = .5,
                                  rotationInterpolation = 0,
                                  scheme = 3,
                                  high_connectivity = FALSE,
                                  tolerance = 0){

  #transform both surface matrices into polar coordinates after centering on
  #firing pin hole by interpolation
  reference_polar <- reference %>%
    fpCenterMat(scheme = scheme,
                high_connectivity = high_connectivity,
                tolerance = tolerance) %>%
    car2pol(method = polarInterpolation,
            rotationResolution = rotationResolution,
            cropWhitespace = FALSE)

  target_polar <- target %>%
    fpCenterMat(scheme = scheme,
                high_connectivity = high_connectivity,
                tolerance = tolerance) %>%
    car2pol(method = polarInterpolation,
            rotationResolution = rotationResolution,cropWhitespace = FALSE)

  #Finding optimal rotation via FFT-based CCF requires replacing missing values
  reference_polar$surface.matrix[is.na(reference_polar$surface.matrix)] <- mean(reference_polar$surface.matrix,na.rm = TRUE)
  target_polar$surface.matrix[is.na(target_polar$surface.matrix)] <- mean(target_polar$surface.matrix,na.rm = TRUE)

  ccfMat <- cmcR:::filterViaFFT(reference_polar$surface.matrix,target_polar$surface.matrix)
  # ccfMat <- imager::correlate(im = as.cimg(reference_polar$surface.matrix),
  #                             filter = as.cimg(target_polar$surface.matrix)) %>%
  #   as.matrix()

  #Find index at which maximum CCF occurs. The corresponding column represents
  #the rotation at which the two scans best align. Need to multiple by the
  #rotation resolution to get back the correct rotation value
  alignment <- (which(ccfMat == max(ccfMat),arr.ind = TRUE) - dim(ccfMat)/2)*rotationResolution

  #Imager by default pads the images with 0s, so we need to play some tricks to
  #make sure we can distinguish the original height values from these fake 0s
  #after rotating. We re-scale the values (to avoid numerical issues with
  #interpolation) and add 1, which we undo later on
  target_rotated <- target
  target_rotated$surface.matrix <- target_rotated$surface.matrix*1e5 + 1

  target_rotated$surface.matrix <- target_rotated$surface.matrix %>%
    imager::as.cimg() %>%
    imager::imrotate(cx = nrow(.)/2,
                     cy = ncol(.)/2,
                     angle = alignment[2],
                     interpolation = rotationInterpolation) %>%
    as.matrix()

  #Since we rescaled and shifted "true" values above, we know that any of the 0s
  #in the scan are added by imager
  target_rotated$surface.matrix[target_rotated$surface.matrix == 0] <- NA

  target_rotated$surface.matrix <- (target_rotated$surface.matrix - 1)/1e5

  return(target_rotated)
}


#Helper function to pad the rows/columns of a matrix
padByDim <- function(reference,
                     target,
                     dimToPad,
                     side = "pre"){

  #This function assumes that dimToPad represents the difference in dimension
  #between the reference scan and target scan. It will pad whichever scan has
  #the smaller dimension (represented by a negative value in dimToPad)

  if(side == "pre"){
    if(dimToPad[1] < 0){
      reference$surface.matrix <- rbind(matrix(NA,
                                               nrow = abs(dimToPad[1]),
                                               ncol = ncol(reference$surface.matrix)),
                                        reference$surface.matrix)
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[1] > 0){
      target$surface.matrix <- rbind(matrix(NA,
                                            nrow = abs(dimToPad[1]),
                                            ncol = ncol(target$surface.matrix)),
                                     target$surface.matrix)
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }

    if(dimToPad[2] < 0){
      reference$surface.matrix <- cbind(matrix(NA,
                                               ncol = abs(dimToPad[2]),
                                               nrow = nrow(reference$surface.matrix)),
                                        reference$surface.matrix)
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[2] > 0){
      target$surface.matrix <- cbind(matrix(NA,
                                            ncol = abs(dimToPad[2]),
                                            nrow = nrow(target$surface.matrix)),
                                     target$surface.matrix)
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }
  }
  if(side == "post"){
    if(dimToPad[1] < 0){
      reference$surface.matrix <- rbind(reference$surface.matrix,
                                        matrix(NA,
                                               nrow = abs(dimToPad[1]),
                                               ncol = ncol(reference$surface.matrix)))
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[1] > 0){
      target$surface.matrix <- rbind(target$surface.matrix,
                                     matrix(NA,
                                            nrow = abs(dimToPad[1]),
                                            ncol = ncol(target$surface.matrix)))
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }

    if(dimToPad[2] < 0){
      reference$surface.matrix <- cbind(reference$surface.matrix,
                                        matrix(NA,
                                               ncol = abs(dimToPad[2]),
                                               nrow = nrow(reference$surface.matrix)))
      reference$header.info$sizeY <- ncol(reference$surface.matrix)
      reference$header.info$sizeX <- nrow(reference$surface.matrix)
    }
    if(dimToPad[2] > 0){
      target$surface.matrix <- cbind(target$surface.matrix,
                                     matrix(NA,
                                            ncol = abs(dimToPad[2]),
                                            nrow = nrow(target$surface.matrix)))
      target$header.info$sizeY <- ncol(target$surface.matrix)
      target$header.info$sizeX <- nrow(target$surface.matrix)
    }
  }

  return(list(reference,target))
}

#Pads two scans such that (1) their central pixels coincide with their
#respective (estimated) center of the firing pin hole and (2) the matrices are
#the same size. If preProcess_rotateScan is used before this function, we can
#assume that the two scans are translationally & rotationally aligned. The
#differences argument then indicates whether only the intersection of the two
#surfaces should be considered.
alignScans <- function(reference,
                       target,
                       differences = "remove",
                       reference_center = NULL,
                       target_center = NULL){

  if(is.null(reference_center)){
    reference_center <- reference %>%
      fpCenterCalc()
  }

  if(is.null(target_center)){
    target_center <- target %>%
      fpCenterCalc()
  }

  #Pad the matrices so that the centers of their respective firing pin impression
  #holes occupy the same index (shoves both matrices into the top-left corner)
  centerPadded <- padByDim(reference,
                           target,
                           dimToPad = reference_center - target_center,
                           side = "pre")

  #Then perform extra padding to make both matrices the same dimension (shoves
  #both matrices into the top-left corner)
  sameSizePadded <- padByDim(centerPadded[[1]],
                             centerPadded[[2]],
                             dimToPad = dim(centerPadded[[1]]$surface.matrix) - dim(centerPadded[[2]]$surface.matrix),
                             side = "post")

  reference_padded <- sameSizePadded[[1]]
  target_padded <- sameSizePadded[[2]]

  #After padding, fp center index should be same for both scans
  newCenterIndex <- fpCenterCalc(reference_padded)

  #Pad both scans so that the center index of the overall matrices are the same
  #as the fp center
  dimToPad <- 2*round(newCenterIndex - dim(reference_padded$surface.matrix)/2)

  #Pad rows:
  if(dimToPad[1] < 0){
    reference_padded$surface.matrix <- rbind(matrix(NA,
                                                    nrow = abs(dimToPad[1]),
                                                    ncol = ncol(reference_padded$surface.matrix)),
                                             reference_padded$surface.matrix)
    reference_padded$header.info$sizeY <- ncol(reference_padded$surface.matrix)
    reference_padded$header.info$sizeX <- nrow(reference_padded$surface.matrix)

    target_padded$surface.matrix <- rbind(matrix(NA,
                                                 nrow = abs(dimToPad[1]),
                                                 ncol = ncol(target_padded$surface.matrix)),
                                          target_padded$surface.matrix)
    target_padded$header.info$sizeY <- ncol(target_padded$surface.matrix)
    target_padded$header.info$sizeX <- nrow(target_padded$surface.matrix)
  }
  if(dimToPad[1] > 0){
    reference_padded$surface.matrix <- rbind(reference_padded$surface.matrix,
                                             matrix(NA,
                                                    nrow = abs(dimToPad[1]),
                                                    ncol = ncol(reference_padded$surface.matrix)))
    reference_padded$header.info$sizeY <- ncol(reference_padded$surface.matrix)
    reference_padded$header.info$sizeX <- nrow(reference_padded$surface.matrix)

    target_padded$surface.matrix <- rbind(target_padded$surface.matrix,
                                          matrix(NA,
                                                 nrow = abs(dimToPad[1]),
                                                 ncol = ncol(target_padded$surface.matrix)))
    target_padded$header.info$sizeY <- ncol(target_padded$surface.matrix)
    target_padded$header.info$sizeX <- nrow(target_padded$surface.matrix)
  }

  #Pad cols:
  if(dimToPad[2] < 0){
    reference_padded$surface.matrix <- cbind(matrix(NA,
                                                    ncol = abs(dimToPad[2]),
                                                    nrow = nrow(reference_padded$surface.matrix)),
                                             reference_padded$surface.matrix)
    reference_padded$header.info$sizeY <- ncol(reference_padded$surface.matrix)
    reference_padded$header.info$sizeX <- nrow(reference_padded$surface.matrix)

    target_padded$surface.matrix <- cbind(matrix(NA,
                                                 ncol = abs(dimToPad[2]),
                                                 nrow = nrow(target_padded$surface.matrix)),
                                          target_padded$surface.matrix)
    target_padded$header.info$sizeY <- ncol(target_padded$surface.matrix)
    target_padded$header.info$sizeX <- nrow(target_padded$surface.matrix)
  }
  if(dimToPad[2] > 0){
    reference_padded$surface.matrix <- cbind(reference_padded$surface.matrix,
                                             matrix(NA,
                                                    ncol = abs(dimToPad[2]),
                                                    nrow = nrow(reference_padded$surface.matrix)))
    reference_padded$header.info$sizeY <- ncol(reference_padded$surface.matrix)
    reference_padded$header.info$sizeX <- nrow(reference_padded$surface.matrix)

    target_padded$surface.matrix <- cbind(target_padded$surface.matrix,
                                          matrix(NA,
                                                 ncol = abs(dimToPad[2]),
                                                 nrow = nrow(target_padded$surface.matrix)))
    target_padded$header.info$sizeY <- ncol(target_padded$surface.matrix)
    target_padded$header.info$sizeX <- nrow(target_padded$surface.matrix)
  }


  #If desired, will remove any non-missing pixels in one surface matrix that are
  #missing in the other. Could extend this function to "impute" the values
  #(e.g., based on MRF or a simpler interpolation scheme)
  if(differences == "remove"){
    reference_padded$surface.matrix[is.na(target_padded$surface.matrix)] <- NA

    target_padded$surface.matrix[is.na(reference_padded$surface.matrix)] <- NA
  }

  #there may be some redundant padding that we can now remove
  # reference_naRows <- reference_padded$surface.matrix %>%
  #   is.na() %>%
  #   rowSums()
  #
  # reference_rowsToKeep <- {reference_naRows < ncol(reference_padded$surface.matrix)}
  #
  # reference_naCols <- reference_padded$surface.matrix %>%
  #   is.na() %>%
  #   colSums()
  #
  # reference_colsToKeep <- {reference_naCols < nrow(reference_padded$surface.matrix) }
  #
  # target_naRows <- target_padded$surface.matrix %>%
  #   is.na() %>%
  #   rowSums()
  #
  # target_rowsToKeep <- {target_naRows < ncol(target_padded$surface.matrix)}
  #
  # target_naCols <- target_padded$surface.matrix %>%
  #   is.na() %>%
  #   colSums()
  #
  # target_colsToKeep <- {target_naCols < nrow(target_padded$surface.matrix)}
  #
  # reference_padded$surface.matrix <- reference_padded$surface.matrix[reference_rowsToKeep & target_rowsToKeep,
  #                                                                    reference_colsToKeep & target_colsToKeep]
  #
  # target_padded$surface.matrix <- target_padded$surface.matrix[reference_rowsToKeep & target_rowsToKeep,
  #                                                              reference_colsToKeep & target_colsToKeep]
  #
  # reference_padded$header.info$sizeY <- ncol(reference_padded$surface.matrix)
  # reference_padded$header.info$sizeX <- nrow(reference_padded$surface.matrix)
  #
  # target_padded$header.info$sizeY <- ncol(target_padded$surface.matrix)
  # target_padded$header.info$sizeX <- nrow(target_padded$surface.matrix)

  return(list(reference_padded,
              target_padded))
}

#Applies the correlate functions from the imager package (which is a C
#implementation of the computationally-intesive definition of the
#cross-correlation function - not using the faster FFT method). Note that this C
#implementation still requires that missing values be replaced in the scan. The
#periodic argument dictates whether the periodic boundary conditions are
#applied. The function returns either the entire CCF matrix or just the maximum
#value/index in the matrix
imagerCCF <- function(reference,
                      target,
                      periodic = TRUE,
                      normalise = TRUE,
                      ccfMatReturn = TRUE){

  reference$surface.matrix[is.na(reference$surface.matrix)] <- mean(reference$surface.matrix,na.rm = TRUE)
  target$surface.matrix[is.na(target$surface.matrix)] <- mean(target$surface.matrix,na.rm = TRUE)

  #Enforce periodic conditions by copying reference 3 times
  if(periodic){
    reference$surface.matrix <- cbind(reference$surface.matrix,
                                      reference$surface.matrix,
                                      reference$surface.matrix)

  }

  #Calculate full CCF
  ccfMat <- imager::correlate(im = as.cimg(target$surface.matrix),
                              filter = as.cimg(reference$surface.matrix),
                              dirichlet = TRUE,
                              normalise = normalise) %>%
    as.matrix()

  if(periodic){
    #Return center matrix representing the CCF when periodic boundary conditions
    #are enforced
    ccfMat <- ccfMat[,(round(ncol((ccfMat))/3) + 1):round(2*ncol((ccfMat))/3)]
  }

  #Return the whole ccf matrix for further processing
  if(ccfMatReturn){
    return(ccfMat)
  }
  #Or just return the maximum value/indices
  else{
    #Assumed theta resolution is .5 here
    maxInd <- (which(ccfMat == max(ccfMat),arr.ind = TRUE) - dim(ccfMat)/2) * .5

    #Figure out what the interpretation of the first element of maxInd is
    #(something to do with scale?). Multiply by -1 to interpret as aligning
    #reference to target than vice versa
    return(c(max(ccfMat),-1*maxInd))
  }
}

calcVarianceRatio <- function(cmcData,similarityCol = "cmcCount"){
  grand_similarityColAverage <- mean(unlist(cmcData[,similarityCol]))

  withinGroup_similarityCol <- cmcData %>%
    group_by(type) %>%
    summarise(similarityColAverage = mean(!!as.name(similarityCol)),
              similarityColVar = var(!!as.name(similarityCol)),
              .groups = "drop")

  betweenGroupVariability <- withinGroup_similarityCol %>%
    mutate(similarityColSS = (similarityColAverage - grand_similarityColAverage)^2) %>%
    pull(similarityColSS) %>%
    sum()

  withinGroupVariability <- withinGroup_similarityCol %>%
    pull(similarityColVar) %>%
    sum()

  cmcData <- cmcData %>%
    mutate(varRatio = betweenGroupVariability/withinGroupVariability)

  return(cmcData)
}

calcAUC <- function(cmcData){
  decisionRuleAUC <- cmcData %>%
    mutate(type = factor(type,levels = c("non-match","match"))) %>%
    pROC::roc(response = type,predictor = cmcCount,
              levels = c("non-match","match"),
              quiet = TRUE)

  cmcData %>%
    mutate(AUC = as.numeric(decisionRuleAUC$auc))
}

preProcess_topMatch <- function(fileName,
                                downsample = 4,
                                minimumHeightCutoff = .3,
                                exteriorRadiusOffset = -20,
                                firingPinRadius = 160,
                                whitespaceCroppingThresh = 1,
                                roughEstimateExterior = FALSE){

  suppressWarnings({
    #Read-in and immediately downsample scan
    dat1_sampled <- fileName %>%
      x3ptools::x3p_read() %>%
      x3ptools::x3p_sample(m = downsample)

    #Get rid of lower-third of values. This makes scans look more like
    #Sensofar-scanned scans
    dat1_filtered <- dat1_sampled$surface.matrix
    dat1_filtered[dat1_filtered < quantile(dat1_filtered,minimumHeightCutoff,na.rm = TRUE)] <- NA

    #Crop out some of the exterior
    dat1_naRows <- dat1_filtered %>%
      is.na() %>%
      rowSums()
    dat1_naCols<- dat1_filtered %>%
      is.na() %>%
      colSums()


    dat1_filtered1 <- dat1_filtered
    dat1_filtered1[c(1:min(which(dat1_naRows > 0)),max(which(dat1_naRows > 0)):nrow(dat1_filtered1)),] <- NA
    dat1_filtered1[,c(1:min(which(dat1_naCols > 0)),max(which(dat1_naCols > 0)):ncol(dat1_filtered1))] <- NA
    dat1_filtered <- dat1_sampled
    dat1_filtered$surface.matrix <- dat1_filtered1

    #Crop-out more of the exterior
    dat1_exterior <- dat1_filtered %>%
      cmcR::preProcess_crop(region = "exterior",
                            croppingThresh = whitespaceCroppingThresh,
                            roughEstimateExterior = roughEstimateExterior,
                            radiusOffset = exteriorRadiusOffset)

    #Identify firing pin hole based on labeling algorithm
    dat1_fp <- dat1_exterior$surface.matrix %>%
      imager::as.cimg() %>%
      as.data.frame() %>%
      dplyr::rename(x = y,
                    y = x) %>%
      dplyr::filter((((y - nrow(dat1_exterior$surface.matrix)/2)^2 + (x - ncol(dat1_exterior$surface.matrix)/2)^2 < (.9*min(dim(dat1_exterior$surface.matrix))/2)^2) & is.na(value))) %>%
      dplyr::select(c(y,x)) %>%
      dplyr::summarise(centerRow = round(mean(y)),
                       centerCol = round(mean(x)),
                       rad = round(sqrt(nrow(.)/pi)))

    #Filter-out observations within some radius of the firing pin hole (160 here)
    dat1_interior <- dat1_exterior
    dat1_interior$surface.matrix <- dat1_exterior$surface.matrix %>%
      imager::as.cimg() %>%
      as.data.frame() %>%
      # rename(x = y,
      #        y = x) %>%
      dplyr::mutate(value = ifelse((y - dat1_fp$centerCol)^2 + (x - dat1_fp$centerRow)^2 < firingPinRadius^2,NA,value)) %>%
      imager::as.cimg() %>%
      as.matrix()

    #Level the scans
    dat1_leveled <- dat1_interior %>%
      cmcR::preProcess_removeTrend(statistic = "quantile",
                                   method = "fn",
                                   tau = .5)
  })

  return(dat1_leveled)
}

#If the x3p has a mask (i.e., manual labelling), then remove labeled observations
filterByMask <- function(x3p){
  mat <- t(x3p$surface.matrix)

  mat[x3p$mask != paste0(toupper(x3p$matrix.info$Mask$Background[[1]]),"FF")] <- NA

  x3p$surface.matrix <- t(mat)

  x3p <- x3p %>%
    x3p_to_df() %>%
    select(c(x,y,value)) %>%
    df_to_x3p()
}

x3p_delete <- function(x3p, mask_vals) {
  idx <- which(t(as.matrix(x3p$mask)) %in% mask_vals)
  x3p$surface.matrix[idx] <- NA
  x3p
}

x3p_trim_na <- function(x3p) {

  rows <- apply(x3p$surface.matrix, MARGIN=1, FUN=function(x) sum(is.na(x)))
  idx <- which(rows == dim(x3p$surface.matrix)[2])
  if (idx[1] != 1) {
    xmin <- 1
  } else xmin <- max(which(idx==1:length(idx))) +1
  if (idx[length(idx)] != dim(x3p$surface.matrix)[1]) {
    xmax <- dim(x3p$surface.matrix)[1]
  } else xmax <- idx[min(which(idx==(dim(x3p$surface.matrix)[1] - rev(seq_along(idx))+1)))]-1

  cols <- apply(x3p$surface.matrix, MARGIN=2, FUN=function(x) sum(is.na(x)))
  idx <- which(cols == dim(x3p$surface.matrix)[1])

  if (idx[1] != 1) {
    ymin <- 1
  } else ymin <- max(which(idx==1:length(idx))) +1
  if (idx[length(idx)] != dim(x3p$surface.matrix)[2]) {
    ymax <- dim(x3p$surface.matrix)[2]
  } else ymax <- idx[min(which(idx==(dim(x3p$surface.matrix)[2] - rev(seq_along(idx))+1)))]-1

  x3p %>% x3p_crop(x = xmin, y = dim(x3p$surface.matrix)[2]-ymax, width = xmax-xmin, height = ymax-ymin)
}

decision_convergence <- function(cellIndex,
                                 x,
                                 y,
                                 theta,
                                 corr,
                                 direction,
                                 translationThresh = 25,
                                 thetaThresh = 3,
                                 corrThresh = .4){

  comparisonFeaturesDF <- data.frame(cellIndex = cellIndex,
                                     x = x,
                                     y = y,
                                     theta = theta,
                                     corr = corr,
                                     direction = direction)

  convergenceIndicators <- comparisonFeaturesDF %>%
    dplyr::group_by(theta,direction) %>%
    dplyr::mutate(distanceToMed = sqrt((x - median(x))^2 + (y - median(y))^2)) %>%
    dplyr::summarise(distanceToMed = median(distanceToMed)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(direction) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(~ {
      data.frame(direction = unique(.$direction),
                 theta = unique(.$theta[which(.$distanceToMed == min(.$distanceToMed))]),
                 distanceToMed = min(.$distanceToMed))
    })

  thetaDistanceInd <- convergenceIndicators %>%
    dplyr::pull(theta) %>%
    abs() %>%
    diff() %>%
    abs() %>%
    magrittr::is_weakly_less_than(thetaThresh)

  distanceToMedInd <- convergenceIndicators %>%
    dplyr::pull(distanceToMed) %>%
    magrittr::is_weakly_less_than(translationThresh)

  convergenceInd <- thetaDistanceInd & all(distanceToMedInd)

  if(!convergenceInd){
    return("non-CMC (failed)")
  }

  thetaRefs <- comparisonFeaturesDF %>%
    group_by(cellIndex,direction) %>%
    filter(corr == max(corr)) %>%
    ungroup() %>%
    group_by(direction) %>%
    summarise(thetaRef = median(theta))


  convergenceCMCs <- comparisonFeaturesDF %>%
    left_join(thetaRefs,
              by = c("direction")) %>%
    group_by(direction)  %>%
    filter(theta >= thetaRef - thetaThresh & theta <= thetaRef + thetaThresh &
             abs(x - median(x)) <= translationThresh &
             abs(y - median(y)) <= translationThresh) %>%
    filter(corr >= corrThresh) %>%
    ungroup() %>%
    group_by(direction,cellIndex) %>%
    filter(corr == max(corr)) %>%
    mutate(convergenceCMCClassif = "CMC")

  comparisonFeaturesDF %>%
    left_join(convergenceCMCs,
              by = c("cellIndex","x","y","corr","theta","direction")) %>%
    mutate(convergenceCMCClassif = ifelse(is.na(convergenceCMCClassif),"non-CMC",convergenceCMCClassif)) %>%
    pull(convergenceCMCClassif) %>%
    return()
}

x3pToDF <- function(x3p,scale_xy = FALSE){

  ret <- expand.grid("y" = 1:nrow(x3p$surface.matrix),
                     "x" = 1:ncol(x3p$surface.matrix)) %>%
    mutate(value = as.vector(x3p$surface.matrix),
           x = x*ifelse(scale_xy,x3p$header.info$incrementX*1e6,1),
           y = y*ifelse(scale_xy,x3p$header.info$incrementX*1e6,1))

  return(ret)
}

estimatedRotationApp <- function(reference,
                                 target,
                                 reference_v_target_comparison = NULL) {

  polarRotation <- polarEstimateRotation(reference = reference,target = target)

  scansTranslationallyAligned <- alignScans(reference = reference,
                                            target = target,
                                            differences = "do nothing")

  reference <- scansTranslationallyAligned[[1]]
  target <- scansTranslationallyAligned[[2]]

  require(shiny)
  app <- shinyApp(
    ui = fluidPage(
      sidebarLayout(
        sidebarPanel(width = 2,
                     numericInput("referenceAlpha", "Reference Alpha Level", min = 0, max = 1, value = .7),
                     numericInput("targetAlpha", "Target Alpha Level", min = 0, max = 1, value = .7),
                     checkboxInput("overlayCheck","Overlay Scans",value = TRUE),
                     numericInput("plotWidth",label = "Plot Width (px)",min = 1,value = 1000),
                     radioButtons("rotationEstimation","Rotate target scan by angle estimated via:",
                                  choices = c("Polar Transformation" = 2,
                                              "Custom" = 1),
                                  selected = 1),
                     conditionalPanel("input.rotationEstimation == '1'",
                                      inputPanel(numericInput(inputId = "customAngle",
                                                              label = "Custom Angle (deg):",
                                                              value = 0,min = -180,max = 180))),
                     numericInput("corrThresh","Correlation Threshold",value = .5,min = 0,max = 1),
                     numericInput("translationThresh","Translation Threshold",value = 20,min = 0)),
        mainPanel(width = 10,
                  fluidRow(column(6,plotOutput("overlayPlot"))),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  fluidRow(column(6,tableOutput("estimatedRotations"))
                           # ,column(6,plotOutput("differencePlot"))
                           ))
      )
    ),
    server = function(input, output) {

      plotWidth <- reactive({input$plotWidth})

      originalMethodRotation <- NA
      highCMCRotation <- NA

      makeReactiveBinding("originalMethodRotation")
      makeReactiveBinding("highCMCRotation")

      if(!is.null(reference_v_target_comparison)){

        updateRadioButtons(inputId = "rotationEstimation",
                           choices = c("Polar Transformation" = 2,
                                       "Original Method" = 3,
                                       "High CMC Method" = 4,
                                       "Custom" = 1),
                           selected = 1)

        toListen <- reactive({
          list(input$translationThresh,input$corrThresh)
        })
        observeEvent(toListen(), {

          originalMethodRotation <<- reference_v_target_comparison %>%
            group_by(cellIndex) %>%
            filter(pairwiseCompCor == max(pairwiseCompCor)) %>%
            ungroup() %>%
            filter(pairwiseCompCor >= input$corrThresh) %>%
            pull(theta) %>%
            median()

          highCMCRotation <<- reference_v_target_comparison %>%
            mutate(cmcThetaDistrib = cmcR::decision_highCMC_cmcThetaDistrib(cellIndex = cellIndex,x = x,y = y,
                                                                            theta = theta,corr = pairwiseCompCor,
                                                                            xThresh = input$translationThresh,
                                                                            corrThresh = input$corrThresh)) %>%
            group_by(theta) %>%
            summarise(cmcCandidateCount = sum(cmcThetaDistrib == "CMC Candidate"),.groups = "drop") %>%
            ungroup() %>%
            top_n(n = 1,wt = cmcCandidateCount) %>%
            pull(theta) %>%
            median()

        })

      }

      targetRotated <- target

      difference <- bind_rows(reference %>%
                                x3pToDF() %>%
                                mutate(x3pName = "reference"),
                              targetRotated %>%
                                x3pToDF() %>%
                                mutate(x3pName = "target"))  %>%
        pivot_wider(id_cols = c(x,y,x3pName),names_from = x3pName,values_from = value)  %>%
        mutate("Overlap.Difference" = reference - target)

      makeReactiveBinding("targetRotated")
      makeReactiveBinding("difference")
# browser()
      toListen1 <- reactive({
        list(input$rotationEstimation,input$customAngle)
      })

      observeEvent(toListen1(),
                   {

                     # if(input$rotationEstimation == 1){
                     #
                     #   targetRotated <<- target
                     #
                     # }
                     # else{

                       target_rotated <- target
                       target_rotated$surface.matrix <- target_rotated$surface.matrix*1e5 + 1
# browser()
                       target_rotated$surface.matrix <- target_rotated$surface.matrix %>%
                         imager::as.cimg() %>%
                         imager::imrotate(cx = nrow(.)/2,
                                          cy = ncol(.)/2,
                                          angle = c(input$customAngle,polarRotation,originalMethodRotation,highCMCRotation)[as.numeric(input$rotationEstimation)],
                                          interpolation = 0) %>%
                         as.matrix()

                       #Since we rescaled and shifted "true" values above, we know that any of the 0s
                       #in the scan are added by imager
                       target_rotated$surface.matrix[target_rotated$surface.matrix == 0] <- NA

                       target_rotated$surface.matrix <- (target_rotated$surface.matrix - 1)/1e5

                       targetRotated <<- target_rotated

                       #update difference scan (in the form of a data frame) too
                       difference <<- bind_rows(reference %>%
                                                  x3pToDF() %>%
                                                  mutate(x3pName = "reference"),
                                                targetRotated %>%
                                                  x3pToDF() %>%
                                                  mutate(x3pName = "target"))  %>%
                         pivot_wider(id_cols = c(x,y,x3pName),
                                     names_from = x3pName,
                                     values_from = value)  %>%
                         mutate("Overlap.Difference" = reference - target)

                     # }



                   })



      observe({

        output$overlayPlot <- renderPlot(width = plotWidth(),height = plotWidth()/2,{

          if(input$overlayCheck){

            if(prod(dim(reference$surface.matrix)) > prod(dim(targetRotated$surface.matrix))){

              pltDat <- reference %>%
                x3pToDF() %>%
                dplyr::mutate(x3pName = "reference") %>%
                filter(!is.na(value)) %>%
                mutate(value = abs(value - median(value)))

              pltDat1 <- targetRotated %>%
                x3pToDF() %>%
                dplyr::mutate(x3pName = "targetRotated") %>%
                filter(!is.na(value)) %>%
                mutate(value = abs(value - median(value)))

              plt <- ggplot() +
                geom_raster(data = pltDat,
                            aes(x=x,y=y,fill = value),
                            alpha = input$referenceAlpha) +
                ggplot2::scale_fill_gradientn(colours = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858'),
                                              values = scales::rescale(quantile(as.vector(pltDat$value),
                                                                                c(0,.6,.75,.85,.9,.95,.99,.995,1),
                                                                                na.rm = TRUE)),
                                              breaks = function(lims){
                                                dat <- quantile(as.vector(pltDat$value),
                                                                c(0,.5,.99,1),
                                                                # c(0,.01,.5,.99,1),
                                                                # c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                na.rm = TRUE)

                                                dat <- dat %>%
                                                  setNames(paste0(names(dat)," [",round(dat*1e6,1),"]"))

                                                return(dat)
                                              },
                                              na.value = "transparent",
                                              guide = guide_colorbar(title = "target",
                                                                     barheight = 10,
                                                                     order = 2,
                                                                     label.theme = ggplot2::element_text(size = 8),
                                                                     title.theme = ggplot2::element_text(size = 10))) +
                ggnewscale::new_scale_fill() +
                ggplot2::geom_raster(data = pltDat1,
                                     ggplot2::aes(x = x,y = y,fill = value),
                                     alpha = input$targetAlpha) +
                ggplot2::scale_fill_gradientn(colours = c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'),
                                              values = scales::rescale(quantile(as.vector(pltDat1$value),
                                                                                c(0,.6,.75,.85,.9,.95,.99,.995,1),
                                                                                na.rm = TRUE)),
                                              breaks = function(lims){
                                                dat <- quantile(as.vector(pltDat1$value),
                                                                c(0,.5,.99,1),
                                                                # c(0,.01,.5,.99,1),
                                                                # c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                na.rm = TRUE)

                                                dat <- dat %>%
                                                  setNames(paste0(names(dat)," [",round(dat*1e6,1),"]"))

                                                return(dat)
                                              },
                                              na.value = "transparent",
                                              guide = ggplot2::guide_colorbar(title = "reference",
                                                                              barheight = 10,
                                                                              order = 1,
                                                                              label.theme = ggplot2::element_text(size = 8),
                                                                              title.theme = ggplot2::element_text(size = 10),)) +
                coord_fixed(expand = FALSE) +
                ggplot2::theme_minimal() +
                ggplot2::scale_y_reverse() +
                ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.background = ggplot2::element_blank())

            }
            else{

              # updateNumericInput(inputId = "referenceAlpha",value = input$targetAlpha)
              # updateNumericInput(inputId = "targetAlpha",value = input$targetAlpha)


              pltDat <- targetRotated %>%
                x3pToDF() %>%
                dplyr::mutate(x3pName = "target") %>%
                filter(!is.na(value)) %>%
                mutate(value = abs(value - median(value)))

              pltDat1 <- reference %>%
                x3pToDF() %>%
                dplyr::mutate(x3pName = "reference") %>%
                filter(!is.na(value)) %>%
                mutate(value = abs(value - median(value)))

              plt <- ggplot() +
                geom_raster(data = pltDat,
                            aes(x=x,y=y,fill = value),
                            alpha = input$targetAlpha) +
                ggplot2::scale_fill_gradientn(colours = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858'),
                                              values = scales::rescale(quantile(as.vector(pltDat$value),
                                                                                c(0,.6,.75,.85,.9,.95,.99,.995,1),
                                                                                na.rm = TRUE)),
                                              breaks = function(lims){
                                                dat <- quantile(as.vector(pltDat$value),
                                                                c(0,.5,.99,1),
                                                                # c(0,.01,.5,.99,1),
                                                                # c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                na.rm = TRUE)

                                                dat <- dat %>%
                                                  setNames(paste0(names(dat)," [",round(dat*1e6,1),"]"))

                                                return(dat)
                                              },
                                              na.value = "transparent",
                                              guide = guide_colorbar(title = "target",
                                                                     barheight = 10,
                                                                     order = 2,
                                                                     label.theme = ggplot2::element_text(size = 8),
                                                                     title.theme = ggplot2::element_text(size = 10))) +
                ggnewscale::new_scale_fill() +
                ggplot2::geom_raster(data = pltDat1,
                                     ggplot2::aes(x = x,y = y,fill = value),
                                     alpha = input$referenceAlpha) +
                ggplot2::scale_fill_gradientn(colours = c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'),
                                              values = scales::rescale(quantile(as.vector(pltDat1$value),
                                                                                c(0,.6,.75,.85,.9,.95,.99,.995,1),
                                                                                na.rm = TRUE)),
                                              breaks = function(lims){
                                                dat <- quantile(as.vector(pltDat1$value),
                                                                c(0,.5,.99,1),
                                                                # c(0,.01,.5,.99,1),
                                                                # c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                na.rm = TRUE)

                                                dat <- dat %>%
                                                  setNames(paste0(names(dat)," [",round(dat*1e6,1),"]"))

                                                return(dat)
                                              },
                                              na.value = "transparent",
                                              guide = ggplot2::guide_colorbar(title = "reference",
                                                                              barheight = 10,
                                                                              order = 1,
                                                                              label.theme = ggplot2::element_text(size = 8),
                                                                              title.theme = ggplot2::element_text(size = 10),)) +
                coord_fixed(expand = FALSE) +
                ggplot2::theme_minimal() +
                ggplot2::scale_y_reverse() +
                ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               panel.background = ggplot2::element_blank())
            }
          }
          else{

            pltDat <- bind_rows(reference %>%
                                  x3pToDF() %>%
                                  dplyr::mutate(x3pName = "Reference"),
                                targetRotated %>%
                                  x3pToDF() %>%
                                  dplyr::mutate(x3pName = "Target")) %>%
              filter(!is.na(value))

            plt <- pltDat %>%
              ggplot2::ggplot() +
              ggplot2::geom_raster(ggplot2::aes(x = x,y = y,fill = value)) +
              ggplot2::scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                                                            '#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                                            values = scales::rescale(quantile(pltDat$value,
                                                                              c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                              na.rm = TRUE)),
                                            breaks = function(lims){
                                              dat <- quantile(as.vector(pltDat$value),
                                                              c(0,.01,.5,.99,1),
                                                              # c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                              na.rm = TRUE)

                                              dat <- dat %>%
                                                setNames(paste0(names(dat)," [",round(dat*1e6,1),"]"))

                                              return(dat)
                                            },
                                            na.value = "transparent",
                                            guide = ggplot2::guide_colorbar(title = expression("Rel. Height ["*mu*"m]"),
                                                                            barheight = 10,
                                                                            order = 1,
                                                                            label.theme = ggplot2::element_text(size = 8),
                                                                            title.theme = ggplot2::element_text(size = 10))) +
              facet_wrap(~ x3pName,nrow = 1) +
              ggplot2::coord_fixed(expand = FALSE) +
              ggplot2::theme_minimal() +
              ggplot2::scale_y_reverse() +
              ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                             axis.ticks.y = ggplot2::element_blank(),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.background = ggplot2::element_blank())
          }

          plt

        })
      })

      observeEvent({input$rotationEstimation},
                   {

                     output$differencePlot <- renderPlot({

                       pltDat <- difference %>%
                         filter(!is.na(Overlap.Difference)) %>%
                         mutate(x3pName = "Overlap Reference minus Target")

                       plt <- pltDat %>%
                         ggplot2::ggplot() +
                         ggplot2::geom_raster(ggplot2::aes(x = x,y = y,fill = Overlap.Difference)) +
                         ggplot2::scale_fill_gradientn(colours = rev(c('#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                                                                       '#d8daeb','#b2abd2','#8073ac','#542788','#2d004b')),
                                                         #c("#67a9cf","#86b7d6","#a2c5dd","#bdd3e3","#d7e2ea",
                                                                  # "#f1f1f1","#f5dcd3","#f6c8b6","#f6b399","#f39f7d","#ef8a62"),
                                                       values = scales::rescale(quantile(pltDat$Overlap.Difference,
                                                                                         c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                                         na.rm = TRUE)),
                                                       breaks = function(lims){
                                                         dat <- quantile(as.vector(pltDat$Overlap.Difference),
                                                                         c(0,.01,.5,.99,1),
                                                                         # c(0,.01,.025,.1,.25,.5,.75,0.9,.975,.99,1),
                                                                         na.rm = TRUE)

                                                         dat <- dat %>%
                                                           setNames(paste0(names(dat)," [",round(dat*1e6,1),"]"))

                                                         return(dat)
                                                       },
                                                       na.value = "transparent",
                                                       guide = ggplot2::guide_colorbar(title = expression("Rel. Height ["*mu*"m]"),
                                                                                       barheight = 10,
                                                                                       order = 1,
                                                                                       label.theme = ggplot2::element_text(size = 8),
                                                                                       title.theme = ggplot2::element_text(size = 10))) +
                         facet_wrap(~ x3pName,nrow = 1) +
                         ggplot2::coord_fixed(expand = FALSE) +
                         ggplot2::theme_minimal() +
                         ggplot2::scale_y_reverse() +
                         ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                                        axis.ticks.y = ggplot2::element_blank(),
                                        panel.grid.major = ggplot2::element_blank(),
                                        panel.grid.minor = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_blank())

                       plt

                     })

                   })

      observe({

        output$estimatedRotations <- renderTable({

          data.frame("Method" = rep(c("Polar Transform","Original Method","High CMC Method"),each = 1),
                     "Direction" = rep(c("Reference vs. Target"),times = 3),
                     # "Direction" = rep(c("Reference vs. Target","Target vs. Reference"),times = 3),
                     "Estimated Angle" = c(polarRotation,originalMethodRotation,highCMCRotation))

        })

      })
    }
  )

  shiny::runApp(app)

}

scanPhaseShift <- function(x3p,rotation,rowTranslation,colTranslation){

  mat <- x3p$surface.matrix

  if(rotation != 0){

    mat <- cmcR:::rotateSurfaceMatrix(surfaceMat = mat,
                                      theta = rotation)

  }

  if(rowTranslation != 0){

    rowPad <- matrix(NA,nrow = abs(rowTranslation),
                     ncol = ncol(mat))

    if(rowTranslation < 0){

      mat <- rbind(mat,rowPad)

    }
    else if(rowTranslation > 0){

      mat <- rbind(rowPad,mat)

    }

  }

  if(colTranslation != 0){

    colPad <- matrix(NA,ncol = abs(colTranslation),
                     nrow = nrow(mat))

    if(colTranslation < 0){

      mat <- cbind(mat,colPad)

    }
    else if(colTranslation > 0){

      mat <- cbind(colPad,mat)

    }

  }

  x3p$surface.matrix <- mat
  x3p$header.info$sizeX <- nrow(mat)
  x3p$header.info$sizeY <- ncol(mat)

  return(x3p)
}

calcRegistration <- function(comparisonData,xThresh = 20,yThresh = xThresh,
                             corrThresh = .5,thetaThresh = 6,tau = 1,plot = TRUE){

  if(plot){

    plt1 <- diagnosticHistograms(comparisonData %>%
                                   filter(direction == "refToTarget")) +
      plot_annotation(title = "Reference vs. Target Similarity Features")

    plt2 <- diagnosticHistograms(comparisonData %>%
                                   filter(direction == "targetToRef")) +
      plot_annotation(title = "Target vs. Reference Similarity Features")

  }

  ret1 <- comparisonData %>%
    group_by(cellIndex,direction) %>%
    filter(pairwiseCompCor == max(pairwiseCompCor) & pairwiseCompCor >= corrThresh) %>%
    ungroup() %>%
    group_by(direction) %>%
    summarise(theta = median(theta),
              x = median(x),
              y = median(y)) %>%
    ungroup() %>%
    mutate(method = "original") %>%
    select(method,direction,theta,x,y)

  ret2 <- comparisonData %>%
    group_by(direction) %>%
    group_split() %>%
    map_dfr(~ {

      mutate(.,cmcThetaDistribClassif = cmcR::decision_highCMC_cmcThetaDistrib(cellIndex=cellIndex,x=x,y=y,theta=theta,
                                                                               corr=pairwiseCompCor,
                                                                               xThresh=xThresh,yThresh = yThresh,corrThresh = corrThresh)) %>%
        cmcR::decision_highCMC_identifyHighCMCThetas(tau = tau)

    }) %>%
    group_by(theta,direction)  %>%
    filter(thetaCMCIdentif == "High") %>%
    mutate(thetaRange = max(theta) - min(theta)) %>%
    group_by(direction) %>%
    summarise(theta = median(theta),
              x = median(x),
              y = median(y),
              thetaRange = unique(thetaRange)) %>%
    mutate(theta = ifelse(thetaRange <= thetaThresh,
                          theta,NA),
           x = ifelse(thetaRange <= thetaThresh,
                      x,NA),
           y = ifelse(thetaRange <= thetaThresh,
                      y,NA)) %>%
    select(-thetaRange) %>%
    ungroup() %>%
    mutate(method = "highCMC") %>%
    select(method,direction,theta,x,y)

  if(plot){

    return(list("Reference vs. Target Features" = plt1,
                "Target vs. Reference Features" = plt2,
                "Estimated Alignment" = bind_rows(ret1,ret2)))

  }
  else{
    return(bind_rows(ret1,ret2))
  }
}

preProcess_anisotropicFilter <- function(x3p,
                                         amplitude,
                                         highPass = TRUE, #subtracts-off anisotropically filtered matrix
                                         sharpness = 0.7,
                                         anisotropy = 0.6,
                                         alpha = 0.6,
                                         sigma = 1.1,
                                         dl = 0.8,
                                         da = 30,
                                         gauss_prec = 2,
                                         interpolation_type = 0L,
                                         fast_approx = TRUE){

  mat <- x3p$surface.matrix

  mat[is.na(mat)] <- 0

  mat <- mat %>%
    imager::as.cimg() %>%
    imager::blur_anisotropic(amplitude = amplitude,
                             sharpness = sharpness,
                             anisotropy = anisotropy,
                             alpha = alpha,
                             sigma = sigma,
                             dl = dl,
                             da = da,
                             gauss_prec = gauss_prec,
                             interpolation_type = interpolation_type,
                             fast_approx = fast_approx) %>%
    as.matrix()

  mat[is.na(x3p$surface.matrix)] <- NA

  if(highPass){

    x3p$surface.matrix <- x3p$surface.matrix - mat

  }
  else{

    x3p$surface.matrix <- mat

  }

  return(x3p)

}
