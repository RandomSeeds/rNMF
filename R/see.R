#' @title See vectorized or single images represented by a matrix
#' 
#' @description The function prints multiple or single images on screen, depends 
#' on the format of the input matrix (see details). Different color schemes can 
#' be selected. Build-in color include greyscale, blue-red palette and heat 
#' color. 
#' 
#' @details If the input is a matrix of vectorized images (input = "multi", 
#' default setting), that is, each column contains pixels of one vectorized 
#' image, then see() restores each column into a matrix and show all images in 
#' one frame. Current version assumes the images are squared images.
#' 
#' If the input is a matrix of one image (input = "single"), see() shows this 
#' image.
#' 
#' @param X a numeric matrix. X is either a matrix where each column contains 
#' the pixels of a vectorized image, or X is simply the pixel matrix of one 
#' single image. The type of X is indicated by the argument 'input'.
#' @param title a charactor string. Title of the graph. 
#' @param col a vector string. Defult = "heat". What color scheme to use? 
#' Currently allows:
#' \itemize{
#'   \item "heat" for heat color
#'   \item "br" for (blue-cyan-green-yellow-red) palette
#'   \item "grey" for grey scale
#' }
#' @param input a charactor string. Defult = "multi". Possible values: 
#' \itemize{
#'   \item "multi" if X contains multiple vectorized square images.
#'   \item "single" if X is the matrix of single image.
#' }
#' @param layout a vector of 2 numbers or a charactor string "auto". If layout 
#' = "auto", multiple images will be arranged in approximatedly 9 by 16 ratio. 
#' If layout = c(a,b), then images will be arranged in a rows and b columns.
#' @param ... further arguments to pass to image().
#' 
#' @return Null
#'
#' @export
#' 
#' @examples
#' ## Load a build-in data set Symbols, a 5625 by 30 matrix containing 30 75x75 
#' ## images.
#' data(Symbols)
#' see(Symbols, title = "Sample images of four symbols")

see <- function(X, title = "Image", col = "heat", 
                input = "multi", layout = "auto", ...){
  if(col == "heat"){
    mycols <- heat.colors(256)
  }else if(col == "br"){
    mycols <- colorRampPalette(c("blue", "cyan","green","yellow","red"))(256)
  }else if(col == "grey"){
    mycols <- grey(0:256/256)
  }else{
    stop("wrong col")
  }
  low = min(X)
  high = max(X)
  ncol = ncol(X)
  nrow = nrow(X)
  if(input == "multi"){
    if(layout == "auto"){
      nROW = ceiling(sqrt(9/16 * ncol))
      nCOL = ceiling(ncol / nROW)
    }else{
      nROW = layout[1]
      nCOL = layout[2]
    }
    par(mfrow = c(nROW, nCOL), mar = c(0.5,0.5,0.5,0.5), oma = c(0,0,3,0))
    for(i in 1:ncol){
      the.image = matrix(X[,i], nrow = sqrt(nrow))
      the.image = t(the.image[ncol(the.image):1,])
      image(the.image,
            col = mycols, zlim = c(low, high),
            xaxt = "n", yaxt = "n", asp = 1, ...)
    }
    mtext(title, outer = TRUE, cex = 1)
  }else if(input == "single"){
    image(t(X[nrow:1,]), col = mycols, main = title,
          xaxt = "n", yaxt = "n", asp = nrow/ncol, ...)
  }else{
    stop("Wrong 'input' argument. Possible values are 'multi' or 'single'.")
  }
  invisible()
}
