#' Generates files containing sequential plots of results of rNMF. 
#' 
#' This function takes the result of the 'rnmf' function (a list), and generate files containing sequential plots of reconstructed matrices with 1 and 2 basis vectors. To be more specific, let k be the number of columns of W, k >= 2. 'seq.plot' will generate 2k graphs in one or more png files. The first k images are W[,u] %*% H[u,], u = 1,..,k. The next k images are W[,(u1,u2)] %*% H[(u1,u2),], (u1,u2) = (1,2),...,(k-1,k), (k,1). The basis vectors in W are ordered by the range of values in them. 
#' 
#' @param res A list, result from the 'rnmf()' function.
#' @param width Width of the page(s).
#' @param height Height of the page(s).
#' @param mylayout Layouts of the graphs on pages. For example mylayout = c(4,5) would put 20 graphs on each page, in 4 by 5 matrices
#' @param force logic. Force to produce all pages.
#' 
#' @return Nothing. However files will be created.
#' 
#' @export
#' 
#' @examples
#' ## Load a clean single simulated tumor image.
#' data("Tumor")
#' ## Add 5% corruptions.
#' Tumor.corrupted = Tumor
#' set.seed(1)
#' Tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] = 1
#' ## Do rnmf with different settings
#' res = rnmf(A = Tumor.corrupted, trim = 0.05, variation = "cell", my.seed = 1)
#' seq.plot(res, mylayout = c(2,5))

seq.plot = function(res, width = 1280, height = 720,
    mylayout = c(2,5), force = FALSE)
{
    W = res$W
    H = res$H
    k = ncol(W)
    fit = res$fit
    myasp = nrow(fit)/ncol(fit)
    ## Reorder according to largest value in W columns
    W.ord = rev(order(apply(W,2,function(x){max(x)-min(x)})))
    W = W[,W.ord]
    H = H[W.ord,]
    k = ncol(W)
    num.g.page = mylayout[1] * mylayout[2] ## Number of graphs per page
    num.pages = (2 * k - 1) %/% num.g.page + 1 ## Number of pages
    if(num.pages > 30 & !force) stop(paste(num.pages,"pages will be created (too many). To run it anyway, set argument force = TRUE"))
    mycol = heat.colors(100)
    mybreaks = seq(0, max(fit), length = 101)
    for(i in 1:num.pages){
        filename = paste("seqplotting", Sys.Date(),"-page",i,".png", sep = "")
        png(file = filename, width = width, height = height)
        par(mfrow = mylayout, mar = c(1,1,2,1))
        for(u in ((i-1) * num.g.page + 1) : min((i * num.g.page), 2*k))
        if(u <= k){
            image(cbind(W[,u]) %*% rbind(H[u,]), xaxt = "n", yaxt = "n", col = mycol, breaks = mybreaks, main = paste("Basis:",u), asp = myasp)
        }else{
            b1 = u - k
            b2 = b1 %% k + 1
            image(cbind(W[,c(b1,b2)]) %*% rbind(H[c(b1,b2),]), xaxt = "n", yaxt = "n", col = mycol, breaks = mybreaks, main = paste("Basis:",b1,",",b2), asp = myasp)
        }
        dev.off()
    }
}


## Old code.
## seq.plot = function(graph, k = 5, maxit = 50, alpha = 0, beta = 0.1,
##     trim = 0.05, tol = 0.001,  my.seed = 3, my.mode = "vecall", file = FALSE, to.plot = "fit", cex.main = 1){
##     p = nrow(graph); n = ncol(graph)
##     showsub = function(a, to.plot){
##         if(to.plot == "fit"){
##             numbers = paste(a, sep = ",", collapse = ",")
##             image(matrix(W[,a],p,length(a)) %*% matrix(H[a,],nrow = length(a), ncol = n), xaxt = "n", yaxt = "n", main = paste("W(",numbers,") X H(",numbers,")", sep = ""), cex.main = cex.main)
##         }else if(to.plot == "w"){
##             image(t(matrix(W[,a],p,length(a))), xaxt = "n", yaxt = "n", main = paste("W: ", paste(a, sep = ",", collapse = ",")), cex.main = cex.main)
##         }else if(to.plot == "h"){
##             image(t(matrix(H[a,], nrow = length(a), ncol = n)), xaxt = "n", yaxt = "n", main = paste("H: ", paste(a, sep = ",", collapse = ",")), cex.main = cex.main)
##         }else{
##             stop("Wrong 'to.plot' value. Try 'fit', 'w' or 'h'.")
##         }
##     }
##     ans = RNMF(graph, k = k, maxit = maxit, alpha = alpha, beta = beta, trim = trim, tol = tol, my.seed = my.seed, my.mode = my.mode)
##     W = ans$W
##     H = ans$H
##     fit = ans$fit
##     W.ord = rev(order(apply(W,2,max)))
##     W = W[,W.ord]
##     H = H[W.ord,]
##     if(file != FALSE){ png(file = file, width = 1280, height = 720) }
##     par(mfrow = c(4,8), mar = c(1,1,2,1))
##     image(graph, xaxt = "n", yaxt = "n", col = heat.colors(32), breaks = c(seq(0,1, length = 33)), main = "Original", cex.main = cex.main)
##     showsub(1:5, to.plot)
##     for(m in 1:4){
##         mat = combn(1:5, m)
##         nn = ncol(mat)
##         for(i in 1:nn){showsub(mat[,i],to.plot)}
##     }
##     if(file != FALSE){ dev.off() }
##     invisible(return(ans))
## }

## a = seq.plot(graph, beta = 0.19, trim = 0.05, my.seed = 8, file = "seq3h.png", to.plot = "h")

## ---------------------------------------
## ------ Single face image
## ---------------------------------------

## a = as.matrix(read.csv("face.corruption2.csv"))
## graph = as.matrix(read.csv("tum2.csv"))
## res = RNMF(a, k = 10,trimp = 0.99, my.seed = 2, trim = TRUE, mode = "smooth", maxit = 10, beta = 0.1)

## dev.off()

## par(mfrow = c(2,5), mar = c(0,1,2,1), oma = c(0,1,1,0))
## for(i in 1:10){
##     showme(matrix(res$W[,i], 192, 1) %*% matrix(res$H[i,], 1, 168), full = FALSE, title = paste("W column", i))
## }

## showme(res$fit)
