#' Performs robust non-negative matrix factorization.
#' 
#' Performs robust penalized non-negative matrix factorizaion on a non-negative
#' matrix A to obtain W and H, such that A ~= W %*% H,  where A is a p by n matrix, W is a p by k matrix
#' and H is a k by n matrix. Outliers in A are detected and trimmed (by cells, columns or rows, see details).
#' Here we assume each row of a represent a feature/variable, and each column of A is an observation/sample.
#' The objective function is ||A - W\%*\%H||_{trimmed} + alpha * ||W||_2^2 + beta * sum|H[:,j]|^2. The minimization process is a generalized alternating least square method.
#' 
#' @param A A non-negative numerical matrix. The input non-negative matrix to be decomposed into W \%*\% H.
#' @param k A non-negative integer. The number of columns of W. 
#' @param alpha Non-negative number. Controls ||W||_2^2.
#' @param beta Non-negative number. Controls sparsity of H.
#' @param maxit Positive integer. Maximum number of iterations.
#' @param tol A number between 0 and 1. Tolerance for convergence. We suggest a number smaller than 0.01.
#' @param trim A number between 0 and 1. Percentage of information to be trimmed. If trim = 0, the procedure will be .
#' @param ini.H A numeric matrix. Initial right matrix H in A ~= W \%*\% H. Must have the same number of columns as A, and k rows.
#' @param ini.W A numeric matrix. Initial left matrix W in A ~= W \%*\% H. Must have the same number of rows as A, and k columns.
#' @param ini.zeta A logical matrix. Must have the same dimension as A. A "FALSE" in ini.zeta corresponds to an entry marked as to be trimmed. The number of "FALSE" in ini.zeta must be <=  round(trim * p * n).
#' @param my.seed Numeric. The seed for initialization of W or H. Ignored if ini.H or ini.W are not NULL.
#' @param variation Character string. Which variation to use? The options are: 'cell','col', 'smooth', 'row', 'rowsmooth'
#' @param quiet Logical. Display a short report on screen after the function is run?
#' @param nreg number of runs in "cell" variation. Not tested. DON'T change its default value. (and don't ask questions).
#' @param p1 Don't ask questions.
#' @param n1 Don't ask questions.
#' 
#' @return A list of the following items:
#' "W": Left matrix of the decomposition.
#' "H": Right matrix of the decomposition.
#' "fit": W \%*\% H.
#' "trimmed1": A list of locations of trimmed cells after stage one in each iteration. (to be updated)
#' "trimmed2": A list of locations of trimmed cells after stage two in each iteration. (to be updated)
#' "SS": Sum of squares of the difference. (to be updated)
#' "niter": Number of iterations performed.
#'
#' @keywords Non-negative matrix factorization, Robustness, Sparsity, Trimming
#'
#' @export
#' 
#' @examples
#' ## Load a clean single simulated tumor image.
#' data("Tumor")
#' ## Add 5% corruptions.
#' Tumor.corrupted = Tumor
#' Tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] = 1
#' ## Do rnmf with different settings
#' res.rnmf1 = rnmf(Tumor.corrupted, trim = FALSE)
#' res.rnmf2 = rnmf(Tumor.corrupted, tol = 0.001, trim = 0.06)
#' res.rnmf3 = rnmf(Tumor.corrupted, k = 10, beta = 0.1, tol = 0.001, trim = 0.06, my.seed = 123, variation = "smooth")
#' par(mfrow = c(2,2))
#' image(Tumor.corrupted, main = "Corrupted")
#' image(res.rnmf1$fit, main = "rnmf (no trimming) fit")
#' image(res.rnmf2$fit, main = "rnmf (smooth) fit 2")
#' image(res.rnmf3$fit, main = "rnmf (smooth) fit 3")

rnmf = function(A, k = 5, alpha = 0, beta = 0, maxit = 50, tol = 0.01,
    trim = FALSE, ini.H = NULL, ini.W = NULL, ini.zeta = NULL, my.seed = NULL,
    variation = "cell", quiet = FALSE, nreg = 1, p1 = NA, n1 = NA)
{
    ## initM: Randomly generates a matrix with values between "small" and "large". Uniform distribution. 
    initM = function(large, nrow, ncol, small = 0, my.seed = NULL){
        if(!missing(my.seed)) set.seed(my.seed)
        M = matrix(runif(nrow * ncol, small, large), nrow = nrow, ncol = ncol)
        return(M)
    }
    
    ## Nnls.trimH for variation "cell"
    Nnls.trimH = function(W, A, zeta, beta, k, n){
        fun1 = function(j, W, A, zeta){
            if(all(!c(zeta[,j]))){
                stop("OK, an entire column is trimmed. I really don't know what to do.")
            }else{
                row.keep = c(zeta[,j])
                W.trim = W[row.keep,]
                x.trim = A[row.keep,j]
                W.ext = rbind(W.trim, sqrt(beta) * matrix(1,1,k))
                x.ext = c(x.trim, 0)
            }
            return(nnls(W.ext,x.ext)$x)
        }
        return(sapply(1:n, fun1, W, A, zeta))
    }
    
    ## Nnls.trimW for cell.reg
    Nnls.trimW = function(H, A, zeta, alpha, p1, n1, k, p){
        ## Solution 1: if entire row is trimmed, don't update that row of W.
        ## fun1 = function(i, H, A, zeta){
        ##     if(all(!c(zeta[i,]))){
        ##         ##browser()
        ##         ## These don't work.
        ##         H.ext = rbind(t(H), sqrt(alpha) * diag(k))
        ##         x.ext = c(rep(mean(A[i,]), ncol(A)), rep(0,k))
        ##         return(nnls(H.ext, x.ext)$x)
        ##         ## if(is.na(p1) | is.na(n1)) stop("An entire row is trimmed. I don't know how to smooth it")
        ##         ## return(find.neib(i, p1 = p1, n1 = n1, zeta = zeta, A = A, k = k))
        ##     }else{                   
        ##         col.keep = c(zeta[i,])
        ##         H.trim = H[,col.keep]
        ##         x.trim = A[i,col.keep]
        ##         H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
        ##         x.ext = c(x.trim, rep(0,k))
        ##         return(nnls(H.ext, x.ext)$x)
        ##     }
        ## }
        
        ## Solution 2: if entire row is trimmed, randomly generate a row.
        fun1 = function(i, H, A, zeta){
            if(all(!c(zeta[i,]))){
                ##browser()
                warning(paste("Row",i,"is trimmed entirely in one of the iterations. It is temporarily replaced by random numbers from 0 to 1 to fit H."))
                x.sub = runif(ncol(A), 0, 1)
                H.ext = rbind(t(H), sqrt(alpha) * diag(k))
                x.ext = c(x.sub, rep(0,k))
                return(nnls(H.ext, x.ext)$x)
            }else{
                col.keep = c(zeta[i,])
                H.trim = H[,col.keep]
                x.trim = A[i,col.keep]
                H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
                x.ext = c(x.trim, rep(0,k))
                return(nnls(H.ext, x.ext)$x)
            }
        }
        
        ## ## Solution 3: if entire row is trimmed, randomly select 50% of the row to keep.
        ## fun1 = function(i, H, A, zeta){
        ##     if(all(!c(zeta[i,]))){
        ##         ##browser()
        ##         col.keep = c(zeta[i,])
        ##         col.keep[sample(1:length(col.keep), round(0.5 * length(col.keep)), replace = FALSE )] = TRUE
        ##         H.trim = H[,col.keep]
        ##         x.trim = A[i,col.keep]
        ##         H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
        ##         x.ext = c(x.trim, rep(0,k))
        ##         return(nnls(H.ext, x.ext)$x)
        ##     }else{
        ##         ##browser()
        ##         col.keep = c(zeta[i,])
        ##         H.trim = H[,col.keep]
        ##         x.trim = A[i,col.keep]
        ##         H.ext = rbind(t(H.trim), sqrt(alpha) * diag(k))
        ##         x.ext = c(x.trim, rep(0,k))
        ##     return(nnls(H.ext, x.ext)$x)
        ##     }
        ## }
        return(t(sapply(1:p, fun1, H, A, zeta)))
    }
    
    ## Given a line number in a vectorized column i, find neighboring points of it in a non-vectorized matrix.
    ## p1 = the number of rows of the individual image
    ## n1 = the number of columns of the individual image
    ## zeta = the matrix of points to keep
    ## A = the original image
    find.neib = function(i, p1, n1, zeta, A, k){
        step = 1
        x = find.x(i,p1)
        y = find.y(i,p1)
        while(TRUE){
            p.range = (i - min(step, x - 1)) : (i + min(step, p1 - x)) 
            n.range = (-min(step, y - 1)) : (min(step, n1 - y))
            to.check = c(rep(p.range, length(n.range)) + rep(n.range, each = length(p.range)) * p1)
            if(all(!zeta[to.check,])){
                step = step + 1
            }else{
                return(rep(mean(A[to.check,][zeta[to.check,]]), k))
            }
        }
    }
    
    ## l2 norm.
    l2 = function(x){ return(sqrt(sum(x^2)))} 
    
    ## nmlz: Returns the normalization matrix "J" of columns of M. That is: M %*% J has unit length columns. Zero columns correspond to "1" (unchanged).
    nmlz = function(M){
        J = diag(1/apply(M,2,l2))
        diag(J)[which(diag(J) == Inf)] = 1
        return(J)
    }
    
    ## ## This function reorders (permutes) two matrices A and B which have the same shape so that their column structures are as similar as possible. The output value is the permutation order for the second matrix.
    ## ## WARNING: this function is slow!!
    ## ## NOTE: Currently not enabled in the main function "rNMF".
    ## per = function(A,B){
    ##     fun1 = function(x){
    ##         sum(apply(A - B[,x], 2, l2))
    ##     }
    ##     m = ncol(A)
    ##     per = permn(m)
    ##     diff = lapply(per, fun1)
    ##     return(per[[which.min(diff)]])
    ## }
    
    ## "find.row" takes the (x,y) coordinates in A in the format xypair = c(x,y), and returns the corresponding row index of v.
    ## "find.x" takes the index of v, and returns the x coordinator in A.
    ## "find.y" takes the index of v, and returns the y coordinator in A.
    find.row = function(xypair, p){
        xypair = unlist(xypair)
        return(p * (xypair[2] - 1) + xypair[1])
    }
    find.x = function(i, p){ return((i-1) %% p + 1)}
    find.y = function(i, p){ return((i-1) %/% p + 1)}
    
    ## "smooth2" is used in "smoothing2". It takes input i, the index of the row in image.f that needs to be smoothed. It returns the mean of the closest non-outlier pixels. Searching begins with a disk with radius one. The radius increases if the search fails.
    ## image.f is a data frame that contains x coordinates, y coordinates, pixel values, logic flags "outs".
    ## NOTE: The function needs to be improved by changing the search region from a disk to a square.
    smooth2 = function(i, image.f, p, n){
        r = 1
        x0 = (i-1) %% p + 1
        y0 = (i-1) %/% p + 1
        neibs = which((image.f$x - x0)^2 + (image.f$y - y0)^2 <= 1)
        neibs.index = which(image.f[neibs,"outs"] == FALSE)
        while(length(neibs.index) == 0){
            r = r + 1
            neibs = which((image.f$x - x0)^2 + (image.f$y - y0)^2 <= r^2)
            neibs.index = which(image.f[neibs,"outs"] == FALSE)
        }
        neibs = neibs[neibs.index]
        return(mean(image.f[neibs,"value"]))
    }
    
    ## "smoothing2" uses "smooth2", and is used in "RNMF".
    ## "smoothing2" takes in a "graph" (data frame, or vector, or matrix). It returns either a matrix or a data frame of the smoothed data frame.
    ## NOTE: Examples of "smoothing2" can be found in "smoothing2.r"
    smoothing2 = function(graph, outs, p = NULL, n = NULL, frame = FALSE){
        if(is.vector(graph)){
            if(missing(p) & missing(n)) stop("Missing p and n!")
            graph = data.frame(x = rep(1:p,n), y = rep(1:n,each = p), value = graph)
        }else if(is.matrix(graph)){
            p = nrow(graph); n = ncol(graph)
            graph = data.frame(x = rep(1:p,n), y = rep(1:n,each = p), value = as.vector(graph))
        }
        graph$outs = FALSE
        graph$outs[outs] = TRUE
        graph[outs,"value"] = unlist(lapply(outs, smooth2, graph, p, n))
        if(frame == TRUE){
            return(graph)
        }else{
            return(matrix(graph$value, nrow = p))
        }
    }
    
    ## -------------------
    ## The function starts here
    ## -------------------
    
    tic = proc.time()
    if(!is.matrix(A)) {stop("The input is not a matrix. Consider as.matrix(A)?")}
    if(length(k) != 1) {stop("Something is wrong about k. It should be a positive integer.")}
    if(trim < 0 || trim > 1) {stop("Argument 'trim' must be a value in [0,1), or 'FALSE'.")}
    if(alpha < 0 || beta < 0) {stop("Argument 'alpha' and 'beta' must be non-negative.")}
    p = nrow(A)
    n = ncol(A)
    A.f = data.frame(value = as.vector(A), x = rep(1:p, n), y = rep(1:n, each = p), outs = FALSE)
    if(trim > 0){
        to.trim1 = vector("list", maxit)
        to.trim2 = vector("list", maxit)
    }else{
        to.trim1 = NULL
        to.trim2 = NULL
    }
    ## Initialization.
    ## Initialize W
    ## With mode = "cell", W instead of H as the first step. 
    if(missing(ini.W)){
        if(missing(my.seed)){
            W = initM(large = max(A), nrow = p, ncol = k, small = 0)
        }else{
            W = initM(large = max(A), nrow = p, ncol = k, small = 0, my.seed = my.seed)
        }
    }else{
        if(!is.matrix(ini.W)) {stop("ini.W must be a matrix.")}
        if(nrow(ini.W) != p || ncol(ini.W) != k) {stop("ini.W has the wrong dimension!.")}
        W = ini.W
    }
    ## Initialize H.
    if(missing(ini.H)){
        if(missing(my.seed)){
            H = initM(large = max(A), nrow = k, ncol = n, small = 0)
        }else{
            H = initM(large = max(A), nrow = k, ncol = n, small = 0, my.seed = my.seed)
        }
    }else{
        if(!is.matrix(ini.H)) {stop("ini.H must be a matrix.")}
        if(nrow(ini.H) != k || ncol(ini.H) != n) {stop("ini.H has the wrong dimension!.")}
        H = ini.H
    }
    SS = rep(NA,maxit) ## Initialize some of squares
    ## The following line creates a progress bar.
    pb <- txtProgressBar(min = 0, max = maxit, style = 3)
    obj = 1
    if(variation == "cell" & trim > 0){
        ## zeta denotes points to be kept
        zeta.allTRUE <- matrix(TRUE, nrow = p, ncol = n)
        if(!missing(ini.zeta)){
            if(!is.logical(ini.zeta)) {stop("ini.zeta must be a logical matrix.")}
            if(nrow(ini.zeta) != p || ncol(ini.zeta) != n) {stop("ini.zeta has the wrong dimension. dim(ini.zeta) should be the same as dim(X).")}
            if(sum(c(!ini.zeta)) > round(trim * p * n)) {stop("ini.zeta contains too many FALSES (outliers); Increase trimming percentage?")}
            if(sum(c(!ini.zeta)) < round(trim * p * n)) {
                warning("ini.zeta contains too few FALSES (outliers); other outliers are randomly picked.")
                need.to.fill = round(trim * p * n) - sum(c(!ini.zeta))
                current.true = which(ini.zeta)
                to.change = sample(current.true, need.to.fill)
                ini.zeta[to.change] = FALSE
            }
            zeta = ini.zeta
        }else{
            zeta = zeta.allTRUE
            zeta[sample(1:(p * n), round(trim * p * n), replace = FALSE)] = FALSE ## randomize zeta
        }
        for(i in 1:maxit){
            setTxtProgressBar(pb, i) ## update the progress bar
            ##------Stage 1------##
            ## Fit H
            H = Nnls.trimH(W, A, zeta, beta, k, n)
            for(j in 1:nreg){
                ## Find residuals
                R = (A - W %*% H)^2 
                to.trim1[[i]] = which(rank(R) > round((1 - trim) * n * p))
                ## Update zeta
                zeta = zeta.allTRUE
                zeta[to.trim1[[i]]] = FALSE
                ## Refit H
                H = Nnls.trimH(W, A, zeta, beta, k, n)
            }
            ##------Stage 2------##
            ## Fit W
            W = Nnls.trimW(H, A, zeta, alpha, p1, n1, k, p)
            for(j in 1:nreg){
                ## Find residuals
                ## browser()
                R = (A - W %*% H)^2
                to.trim2[[i]] = which(rank(R) > round((1 - trim) * n * p))
                ## Update zeta
                zeta = zeta.allTRUE
                zeta[to.trim2[[i]]] = FALSE
                ## Refit W
                W = Nnls.trimW(H, A, zeta, alpha, p1, n1, k, p)
                J = nmlz(W) # Find the normalizing matrix of W
                W = W %*% J; H = solve(J) %*% H
            }
            obj[i] = l2((A - W %*% H)[zeta]) + l2(W) + sum(colSums(abs(H))^2)
            if(i > 1){
                if(all(to.trim2[[i]] == to.trim2[[i-1]]) & sum((W - W.prev)^2)/sum(W.prev^2) < tol) break
            }
            W.prev = W
        }
    }else{
        for(i in 1:maxit){
            flag = FALSE
            setTxtProgressBar(pb, i) # update the progress bar
            ## Iteration stage 1. Estimate W. ||MW^T - C||.
            C = rbind(t(A), matrix(0, nrow = k, ncol = p))
            M = rbind(t(H), sqrt(alpha) * diag(k))
            W = t(Nnls(M,C))
            if(trim > 0){
                R = (A - W %*% H)^2
                if(variation == "col" | variation == "all"){
                    rsqr = apply(R, 2, l2) ## Find SS of residuals of columns
                    to.trim1[[i]] = which(rank(rsqr) > round((1 - trim) * n))
                    A.trim = A[,-to.trim1[[i]]]
                    H.trim = H[,-to.trim1[[i]]]
                    ## Fit W with trimmed A and H. 
                    M.trim = rbind(t(H.trim), sqrt(alpha) * diag(k))
                    C.trim = rbind(t(A.trim), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M.trim,C.trim))
                }else if(variation == "cellcol" | variation == "cellall"){
                    to.trim1[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    trim.col = unique((to.trim1[[i]] %/% p) + 1)
                    A.trim = A[,-trim.col]
                    H.trim = H[,-trim.col]
                    ## Fit W with trimmed A and H. 
                    M.trim = rbind(t(H.trim), sqrt(alpha) * diag(k))
                    C.trim = rbind(t(A.trim), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M.trim,C.trim))
                }else if(variation == "smooth"){
                    R = (A - W %*% H)^2
                    to.trim1[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    A.s = matrix(smoothing2(A.f, to.trim1[[i]], p, n, frame = TRUE)$value,p,n)
                    ## Fit W with trimmed A and H. 
                    C.s = rbind(t(A.s), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M,C.s))
                }else if(variation == "rowsmooth"){
                    findmid = function(xvalue, datavalue){
                        rows = (0:(n-1)) * p + xvalue
                        ## all data frame row numbers corresponding to row = xvalue.
                        return(mean(datavalue[rows], na.rm = TRUE))
                    }
                    R = (A - W %*% H)^2
                    trims <- to.trim1[[i]] <- which(rank(abs(R)) > round((1 - trim) * n * p))
                    values = A.f$value ## All matrix values
                    outs = A.f$outs ## Outliers?
                    values[trims] = NA 
                    outs[trims] = TRUE  ## Change outliers to missing values
                    xs = A.f[trims,"x"]
                    uxs = unique(xs)
                    moveto = sapply(uxs, findmid, values)
                    values[trims] = mapvalues(xs, from = uxs, to = moveto)
                    A.s = matrix(values,p,n)
                    all.na = which(is.na(A.s[,1]))
                    if(length(all.na) > 0) stop(paste("Entire rows are trimmed:", all.na))
                    ## Fit W with trimmed A and H. 
                    C.s = rbind(t(A.s), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M,C.s))
                }else if(variation == "row" | variation == "cellrow"){
                    ## No change
                }else{
                    stop("No such mode. Try one of the following: 'col', 'colrow', 'cellcol', 'smooth', 'row', 'cellrow' or 'rowsmooth'")
                }
            }
            ##J = nmlz(W) # Find the normalizing matrix of W
            ##W = W %*% J; H = solve(J) %*% H
            ## Convergence?
            if(i > 1){
                if(setequal(to.trim1[[i]], to.trim1[[i-1]])){
                    if(sum((W - W.prev)^2)/sum(W.prev^2) < tol) break
                }
            }
            W.prev <- W
            ## Iteration stage 2. Estimate H. ||MH - C||
            C = rbind(A, matrix(0,1,n))
            M = rbind(W, sqrt(beta) * matrix(1, 1, k))
            H = Nnls(M,C)
            if(trim > 0){
                R = (A - W %*% H)^2
                if(variation == "col" | variation == "cellcol"){
                    ## No change
                }else if(variation == "cellrow" | variation == "cellall"){
                    to.trim2[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    trim.row = unique(to.trim2[[i]] %% p)
                    A.trim = A[-trim.row,]
                    W.trim = W[-trim.row,]
                    ## Fit H with trimmed A and W. 
                    M.trim = rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
                    C.trim = rbind(A.trim, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M.trim,C.trim)
                }else if(variation == "smooth"){
                    to.trim2[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    A.s = matrix(smoothing2(A.f, to.trim1[[i]], p,n, frame = TRUE)$value,p,n)
                    ## Fit H with trimmed A and W. 
                    C.s = rbind(A.s, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M,C.s)
                }else if(variation == "row" | variation == "all"){
                    rsqr = apply(R, 1, l2) ## Find SS of residuals of rows
                    to.trim2[[i]] = which(rank(rsqr) > round((1 - trim) * p))
                    A.trim = A[-to.trim2[[i]],]
                    W.trim = W[-to.trim2[[i]],]
                    ## Fit H with trimmed A and W. 
                    M.trim = rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
                    C.trim = rbind(A.trim, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M.trim,C.trim)
                }else if(variation == "rowsmooth"){
                    findmid = function(xvalue, datavalue){
                        xrows = (0:(n-1)) * p + xvalue
                        ## all data frame row numbers corresponding to row = xvalue.
                        return(mean(datavalue[xrows], na.rm = TRUE))
                    }
                    R = (A - W %*% H)^2
                    trims <- to.trim2[[i]] <- which(rank(abs(R)) > round((1 - trim) * n * p))
                    ## locations of entries to be trimmed
                    values = A.f$value ## All matrix values
                    outs = rep(FALSE, n * p) ## Initialize outlier flags.
                    outs[trims] = TRUE 
                    values[trims] = NA ## Change outliers to missing values
                    xs = A.f[trims,"x"] 
                    uxs = unique(xs)
                    moveto = sapply(uxs, findmid, values)
                    values[trims] = mapvalues(xs, from = uxs, to = moveto)
                    A.s = matrix(values,p,n)
                    all.na = which(is.na(A.s[,1]))
                    if(length(all.na) > 0) stop(paste("The following entire rows are trimmed:", all.na, ". I don't know how to fill the cells."))
                    ## Fit H with trimmed A and W. 
                    C.s = rbind(A.s, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M,C.s)
                }else{
                    stop("Wrong mode. Try one of the following: 'col', 'cellcol', 'rowrow', 'cellrow', 'smooth','all', 'cellall', 'rowsmooth'")
                }
            }
            SS[i] = (sqrt(sum(A - W %*% H)^2))
        }
    }
    ## Close the link to the progress bar. 
    close(pb)
    fit = W %*% H
    if(quiet == TRUE){
    }else{
        if(trim == FALSE){
            cat("Done. Time used:","\n")
            print(proc.time() - tic)
            cat("No trimming.\n",
                "Input matrix dimension: ",p,"by",n, "\n",
                "Left matrix:",p,"by",k,". Right matrix:",k,"by",n,"\n",
                "alpha = ",alpha,". beta = ",beta,"\n",
                "Number of max iterations = ",maxit,"\n",
                "Number of iterations = ",i,"\n"
                )
        }else{
            cat("Done. Time used:","\n")
            print(proc.time() - tic)
            cat("\n Trimming mode =\"", variation, "\". Proportion trimmed:",trim, "\n",
                "Input matrix dimension: ",p,"by",n, "\n",
                "Left matrix:",p,"by",k,". Right matrix: ",k,"by",n,"\n",
                "alpha = ",alpha,". beta = ",beta,"\n",
                "Number of max iterations = ",maxit,"\n",
                "Number of iterations = ",i,"\n"
                )
        }
    }
    return(invisible(list(W = W, H = H, fit = fit, trimmed1 = to.trim1, trimmed2 = to.trim2, SS = SS, niter = i)))
}

