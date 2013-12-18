## Last update Yifan Xu 2013/09/15
source("rNMF.r")

seq.plot = function(res, file = FALSE, width = 1280, height = 720){
    showsub = function(a, to.plot){
        if(to.plot == "fit"){
            numbers = paste(a, sep = ",", collapse = ",")
            image(matrix(W[,a],p,length(a)) %*% matrix(H[a,],nrow = length(a), ncol = n), xaxt = "n", yaxt = "n", main = paste("W(",numbers,") X H(",numbers,")", sep = ""), cex.main = cex.main)
        }else if(to.plot == "w"){
            image(t(matrix(W[,a],p,length(a))), xaxt = "n", yaxt = "n", main = paste("W: ", paste(a, sep = ",", collapse = ",")), cex.main = cex.main)
        }else if(to.plot == "h"){
            image(t(matrix(H[a,], nrow = length(a), ncol = n)), xaxt = "n", yaxt = "n", main = paste("H: ", paste(a, sep = ",", collapse = ",")), cex.main = cex.main)
        }else{
            stop("Wrong 'to.plot' value. Try 'fit', 'w' or 'h'.")
        }
    }
    W = res$W
    H = res$H
    k = ncol(W)
    fit = res$fit
    ## Reorder according to largest value in W columns
    W.ord = rev(order(apply(W,2,max)))
    W = W[,W.ord]
    H = H[W.ord,]
    if(file != FALSE){ png(file = file, width = width, height = height) }
    par(mfrow = c(4,8), mar = c(1,1,2,1))
    ## The following line plots the original image. 
    ## image(graph, xaxt = "n", yaxt = "n", col = heat.colors(32), breaks = c(seq(0,1, length = 33)), main = "Original", cex.main = cex.main)
    showsub(5, to.plot)
    for(m in 1:4){
        mat = combn(1:5, m)
        nn = ncol(mat)
        for(i in 1:nn){showsub(mat[,i],to.plot)}
    }
    if(file != FALSE){ dev.off() }
    return(invisible(res))
}




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
