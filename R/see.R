#' Visualize Vectorized Images
#' 
#' The function is a wrapper of image(). It arranges and prints multiple or single images. 
#' 
#' If the input is a matrix of vectorized images (input = "multi", 
#' default setting), that is, each column contains pixels of one vectorized 
#' image, then see() restores each column into a matrix and show all images in 
#' one frame. Current version assumes the images are squared images.
#' If the input is a matrix of one image (input = "single"), see() shows this 
#' image. Different color palette can be selected by specify the "col" argument.
#' Build-in color palette includes greyscale, blue-red and heat color. 
#' 
#' @param X a numeric matrix. X is either a matrix where each column contains 
#' the pixels of a vectorized image, or simply the pixel matrix of one 
#' single image. The type of X is indicated by the argument 'input'.
#' @param title a charactor string. Title of the graph. 
#' @param col a character string. Defult = "heat". What color scheme to use? 
#' Currently allows:
#' \itemize{
#'   \item "heat" for heat color
#'   \item "br" for (blue-cyan-green-yellow-red) palette
#'   \item "grey" for grey scale
#' }
#' @param input a charactor string with default = "multi", specifying the type of
#' images in X. Possible options are:
#' \itemize{
#'   \item "multi" if X contains multiple vectorized square images.
#'   \item "single" if X is the matrix of a single image.
#' }
#' @param layout a vector of 2 possible integers or a charactor string "auto" (default). If layout 
#' = "auto", multiple images will be arranged in an approximatedly 9 by 16 ratio. 
#' If layout = c(a,b), then images will be arranged in a rows and b columns.
#' @param ... further arguments to pass to image().
#' 
#' @return NULL
#'
#' @export
#' 
#' @examples
#' ## Load a build-in data set Symbols, a 5625 by 30 matrix containing 30 75x75 
#' ## images.
#' data(Symbols)
#' see(Symbols, title = "Sample images of four symbols")

see <- function(X, title = NULL,
                col = "yellowRed",
                input = "multi",
                layout = "auto",
                mybreaks = NULL,
                myzlim = NULL,
                balance = FALSE,
                legend = FALSE,
                ...){
    if(is.null(myzlim)){
        low <- min(X)
        high <- max(X)
    }else{
        low <- myzlim[1]
        high <- myzlim[2]
    }
    ncol <- ncol(X)
    nrow <- nrow(X)

    if(is.null(mybreaks)) mybreaks <- seq(low, high, length = 257)
    if(col == "heat"){
        mycols <- heat.colors(256)
    }else if(col == "br"){
        mycols <- colorRampPalette(c("blue", "cyan","green","yellow","red"))(256)
    }else if(col == "grey"){
        mycols <- grey(0:256/256)
    }else if(col == "yellowRed"){
        mypalette <- brewer.pal(11, "RdYlGn")
        mycols <- rev(colorRampPalette(mypalette)(512))[257:512]
    }else if(col == "greenRed"){
        mypalette <- brewer.pal(11, "RdYlGn")
        mycols <- rev(colorRampPalette(mypalette)(256))
        if(is.null(mybreaks)){
            if(low < 0){
                if(balance == TRUE){
                    col.low <- -max(c(abs(low), abs(high)))
                    col.high <- -col.low
                }else{
                    col.low <- low
                    col.high <- high
                }
                mybreaks <- seq(col.low, 0, length = 129)
                mybreaks <- c(mybreaks, seq(0, col.high, length = 129)[-1])
            }else{
                mybreaks <- c(seq(-1, 0, length = 129),
                              seq(0, high, length = 129)[-1])
            }
        }
    }else{
        stop("wrong col")
    }

    if(input == "multi"){
        if(layout[1] == "auto"){
            nROW <- ceiling(sqrt(9/16 * ncol))
            nCOL <- ceiling(ncol / nROW)
        }else{
            nROW <- layout[1]
            nCOL <- layout[2]
        }
        if(!is.null(title))
            par(mfrow = c(nROW, nCOL), mar = c(0.5,0.5,0.5,0.5), oma = c(0,0,3,0))
        else
            par(mfrow = c(nROW, nCOL), mar = c(0.5,0.5,0.5,0.5), oma = c(0,0,0,0))
        for(i in 1:ncol){
            myimage <- matrix(X[,i], nrow = sqrt(nrow))
            myimage <- t(myimage[ncol(myimage):1,])
            imgArgs <- list(x = myimage,
                            col = mycols,
                            breaks = mybreaks,
                            xaxt = "n", yaxt = "n",
                            asp = 1, ...)
            do.call(image, imgArgs)
        }
        mtext(title, outer = TRUE, cex = 1)
    }else if(input == "single"){
        par(mar = c(0.5,0.5,3,0.5), oma = c(0,0,2,0))
        image(t(X[nrow:1,]), col = mycols, main = title,
              xaxt = "n", yaxt = "n", asp = nrow/ncol, ...)
    }else{
        stop("Wrong 'input' argument. Possible options are 'multi' or 'single'.")
    }
    if(legend){
        dev.new()
        par(mar = c(0,0,0,2), oma = c(0,0,0,3))
        browser()

        png("legend2.png", res = 200, unit = "in",
            height = 20, width = 2)
        image(x = 1, y = mybreaks,
              z = matrix(mybreaks, 1, length(mybreaks)),
              col = mycols, breaks = mybreaks, xaxt = "n", yaxt = "n",
              xlab = NA, ylab = NA)
        axis(side = 4, at = round(seq(min(mybreaks), max(mybreaks), 0.05),2),  cex.axis = 2)
        dev.off()
    }
    invisible(list(col = mycols, breaks = mybreaks))
}
