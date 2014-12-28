#' Robust penalized non-negative matrix factorization.
#' 
#' This function performs robust penalized non-negative matrix factorizaion on a 
#' non-negative matrix A to obtain W and H, such that A ~ WH. Here A is a p by n 
#' matrix where columns are observations and rows are variables. W is a p by k, 
#' and H is k by n. Outliers in A are detected and trimmed (by cells, columns or 
#' rows, see "details"). The objective function is 
#' 
#' ||A - WH||_{trimmed} + alpha * ||W||_2^2 + beta * sum(|H_.j|)^2.
#' 
#' The minimization process is a generalized alternating least square method.
#' 
#' @param A a nonnegative numerical matrix. The matrix to be decomposed into WH.
#' @param k a nonnegative integer. The lower dimension (number of columns of W).  
#' @param alpha a nonnegative number. Default = 0. It controls the magnitude of 
#' ||W||_2^2.
#' @param beta a nonnegative number. Default = 0. It controls the sparsity of H.
#' @param maxit a positive integer. Default = 50. Maximum number of iterations.
#' @param tol a number between 0 and 1. Default = 0.001. Tolerance for 
#' convergence. We suggest a number smaller than 0.01.
#' @param trim a number between 0 and 1. Default = 0. Percentage of contents of 
#' X to be trimmed. rnmf() performs the regular NMF if trim = 0. 
#' @param ini.H a nonnegative numeric matrix with k rows and the same number of 
#' columns as A. Default = NULL. Initial right matrix H in A ~ WH. Must have the 
#' same number of columns as A and k rows. ini.H is ignored if ini.W is not NULL. 
#' @param ini.W a nonnegative numeric matrix with k columns and the same number
#' of rows as A. Default = NULL. Initial left matrix W in A ~ WH.  
#' @param ini.zeta a logical matrix the same size as A. Default = NULL. Initial
#' location of outliers (FALSE in ini.zeta). The number of "FALSE" in ini.zeta 
#' must be <=  round(trim * p * n).
#' @param my.seed a nonnegative integer. Default = NULL. The random seed for 
#' initialization of W or H. my.seed is ignored if ini.H or ini.W is not NULL.
#' @param variation a character string. Default = "cell". Triming variation. The options are: 'cell', 'col', 'row' and 'smooth'.
#' @param quiet a logic string. Default = FALSE. If quiet == TRUE the on screen 
#' report after rnmf() is run is suppressed.
#' @param nreg a positive integer. Default = 1. Number of runs in the "cell" 
#' variation. Not tested yet. DON'T change its default value. (and don't ask questions).
#' @param showprogress a logic string. Default = TRUE. If TRUE show progress 
#' bar during iteration
#' 
#' @return An object of class 'rnmf', which is a list of the following items:
#' \itemize{
#' \item "W": Left matrix of the decomposition, columns of which are basis 
#' vectors of the low dimension projection.
#' \item "H": Right matrix of the decomposition, columns of which are low 
#' dimensional encoding of the data.
#' \item "fit": W \%*\% H.
#' \item "trimmed": A list of locations of trimmed cells in each iteration
#' \item "niter": Number of iterations performed.
#' }
#'
#' @export
#' 
#' @examples
#' ## Load a clean single simulated tumor image.
#' data("Tumor")
#' ## Add 5\% corruptions.
#' Tumor.corrupted = Tumor
#' set.seed(1)
#' Tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] = 1
#' ## Do rnmf with different settings
#' res.rnmf1 = rnmf(A = Tumor.corrupted, trim = FALSE, my.seed = 1)
#' res.rnmf2 = rnmf(A = Tumor.corrupted, tol = 0.001, trim = 0.06, my.seed = 1)
#' res.rnmf3 = rnmf(A = Tumor.corrupted, k = 10, beta = 0.1, tol = 0.001, trim = 0.06, my.seed = 1, variation = "smooth")
#' par(mfrow = c(2,2))
#' image(Tumor.corrupted, main = "Corrupted")
#' image(res.rnmf1$fit, main = "rnmf (no trimming) fit")
#' image(res.rnmf2$fit, main = "rnmf (cell) fit 2")
#' image(res.rnmf3$fit, main = "rnmf (smooth) fit 3")

rnmf <- function(A, k = 5, alpha = 0, beta = 0, maxit = 50, tol = 0.001,
    trim = FALSE, ini.W = NULL, ini.H = NULL,  ini.zeta = NULL, my.seed = NULL,
    variation = "cell", quiet = FALSE, nreg = 1, showprogress = TRUE)
{
    tic = proc.time() # Start a clock.
    p = nrow(A)
    n = ncol(A)
    ## Checks if the arguments are valid. Display an error message if not.
    checkargs(A, k, alpha, beta, maxit, tol, trim, ini.W, ini.zeta, my.seed, 
              variation, quiet, nreg, p, n)
    
    ## Create a data frame of 4 columns: value = entries of A; 
    ## (x,y) = coordinates; outs = is it an outlier?
    A.f = data.frame(value = as.vector(A), x = rep(1 : p, n), y = rep(1 : n, each = p), outs = FALSE) 

    ## to.trim (list) stores entries trimmed in each iteration.
    if(trim > 0) to.trim = vector("list")
    else to.trim = NULL
    ## 'obj' (vector) stores values of the objective function in each iteration.
    obj = 1

    ## Initialize W
    if(missing(ini.W)){
        if(missing(my.seed)){
            W = initM(large = max(A), nrow = p, ncol = k, small = 0)
        }else{ 
            W = initM(large = max(A), nrow = p, ncol = k, small = 0, my.seed = my.seed)
        }
    }else{
        W = ini.W
    }

    if(showprogress) pb = txtProgressBar(min = 0, max = maxit, style = 3) # Creates a progress bar.
    
    if(variation == "cell" & trim > 0){
        ## Initialize zeta (logic matrix the same size as A), indicating cells to be kept (non-outliers)
        if(!missing(ini.zeta))
          {
            if(sum(c(!ini.zeta)) < round(trim * p * n)) {
                warning("The percentage of FALSES (outliers) in ini.zeta is smaller than 'trim'; Other outliers are randomly picked.")
                need.to.fill = round(trim * p * n) - sum(c(!ini.zeta))
                current.true = which(ini.zeta)
                to.change = sample(current.true, need.to.fill)
                ini.zeta[to.change] = FALSE
            }
            zeta = ini.zeta
        }else
          {
            zeta = matrix(TRUE, nrow = p, ncol = n)
            zeta[sample(1 : (p * n), round(trim * p * n))] = FALSE ## randomize zeta
        }

        ## Start iterations.
        for(i in 1 : maxit){
            if(showprogress) setTxtProgressBar(pb, i) ## update the progress bar
            ##------Stage 1------##
            ## Given W, fit H
            H = Nnls.trimH(W, A, zeta, beta, k, n)
            for(j in 1:nreg){
                ## Find residuals
                R = abs(A - W %*% H)
                to.trim[[i]] = order(R, decreasing = TRUE)[1 : round(trim * n * p)]
                ## to.trim[[i]] = which(rank(R) > round((1 - trim) * n * p))
                ## Update zeta
                zeta <- matrix(TRUE, nrow = p, ncol = n)
                zeta[to.trim[[i]]] <- FALSE
                ## Refit H
                H = Nnls.trimH(W, A, zeta, beta, k, n)
            }
            
            ##------Stage 2------##
            ## Fit W
            W = Nnls.trimW(H, A, zeta, alpha, k, p)
            for(j in 1:nreg) {
                ## Find residuals
                R = abs(A - W %*% H)
                ##to.trim[[i]] = which(rank(R) > round((1 - trim) * n * p))
                to.trim[[i]] = order(R, decreasing=TRUE)[1 : round(trim * n * p)]
                ## Update zeta
                zeta <- matrix(TRUE, nrow = p, ncol = n)
                zeta[to.trim[[i]]] <- FALSE
                ## Refit W
                W = Nnls.trimW(H, A, zeta, alpha, k, p)
                J = nmlz(W) # Find the normalizing matrix of W
                W = W %*% J; H = solve(J) %*% H
            }
            obj[i] = l2((A - W %*% H)[zeta]) + l2(W) + sum(colSums(abs(H))^2)
            if(i > 1){
                if(all(to.trim[[i]] == to.trim[[i - 1]]) & sum((W - W.prev)^2) / sum(W.prev^2) < tol) break
            }
            W.prev = W
        }
    }else{
        for(i in 1 : maxit){
            flag = FALSE
            if(showprogress) setTxtProgressBar(pb, i)  # update the progress bar
            ## Iteration stage 1. Given W, estimate H. ||MH - C||
            C = rbind(A, matrix(0,1,n))
            M = rbind(W, sqrt(beta) * matrix(1, 1, k))
            H = Nnls(M,C)
            if(trim > 0){
                R = (A - W %*% H)^2
                if(variation == "col" | variation == "cellcol"){
                    # No change here
                }else if(variation == "cellrow" | variation == "cellall"){
                    to.trim[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    trim.row = unique(to.trim[[i]] %% p)
                    A.trim = A[-trim.row,]
                    W.trim = W[-trim.row,]
                    ## Fit H with trimmed A and W. 
                    M.trim = rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
                    C.trim = rbind(A.trim, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M.trim,C.trim)
                }else if(variation == "smooth"){
                    to.trim[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    A.s = matrix(smoothing2(A.f, to.trim[[i]], p,n, frame = TRUE)$value,p,n)
                    ## Fit H with trimmed A and W. 
                    C.s = rbind(A.s, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M,C.s)
                }else if(variation == "row" | variation == "all"){
                    rsqr = apply(R, 1, l2) ## Find SS of residuals of rows
                    to.trim[[i]] = which(rank(rsqr) > round((1 - trim) * p))
                    A.trim = A[-to.trim[[i]],]
                    W.trim = W[-to.trim[[i]],]
                    ## Fit H with trimmed A and W. 
                    M.trim = rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
                    C.trim = rbind(A.trim, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M.trim,C.trim)
                }else{
                    stop("Wrong mode. Try one of the following: 'cell', 'col', 'row', 'smooth'")
                }
            }

            ## Iteration stage 2. Find W in min||MW^T - C||.
            C = rbind(t(A), matrix(0, nrow = k, ncol = p))
            M = rbind(t(H), sqrt(alpha) * diag(k))
            W = t(Nnls(M,C))
            if(trim > 0){
                R = (A - W %*% H)^2
                if(variation == "col" | variation == "all"){
                    rsqr = apply(R, 2, l2) ## Find SS of residuals of columns
                    to.trim[[i]] = which(rank(rsqr) > round((1 - trim) * n))
                    A.trim = A[,-to.trim[[i]]]
                    H.trim = H[,-to.trim[[i]]]
                    ## Fit W with trimmed A and H. 
                    M.trim = rbind(t(H.trim), sqrt(alpha) * diag(k))
                    C.trim = rbind(t(A.trim), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M.trim,C.trim))
                }else if(variation == "cellcol" | variation == "cellall"){
                    to.trim[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    trim.col = unique((to.trim[[i]] %/% p) + 1)
                    A.trim = A[,-trim.col]
                    H.trim = H[,-trim.col]
                    ## Fit W with trimmed A and H. 
                    M.trim = rbind(t(H.trim), sqrt(alpha) * diag(k))
                    C.trim = rbind(t(A.trim), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M.trim,C.trim))
                }else if(variation == "smooth"){
                    R = (A - W %*% H)^2
                    to.trim[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    A.s = matrix(smoothing2(A.f, to.trim[[i]], p, n, frame = TRUE)$value,p,n)
                    ## Fit W with trimmed A and H. 
                    C.s = rbind(t(A.s), matrix(0, nrow = k, ncol = p))
                    W = t(Nnls(M,C.s))
                ## }else if(variation == "rowsmooth"){
                ##     findmid = function(xvalue, datavalue){
                ##         rows = (0:(n-1)) * p + xvalue
                ##         ## all data frame row numbers corresponding to row = xvalue.
                ##         return(mean(datavalue[rows], na.rm = TRUE))
                ##     }
                ##     R = (A - W %*% H)^2
                ##     trims <- to.trim[[i]] <- which(rank(abs(R)) > round((1 - trim) * n * p))
                ##     values = A.f$value ## All matrix values
                ##     outs = A.f$outs ## Outliers?
                ##     values[trims] = NA 
                ##     outs[trims] = TRUE  ## Change outliers to missing values
                ##     xs = A.f[trims,"x"]
                ##     uxs = unique(xs)
                ##     moveto = sapply(uxs, findmid, values)
                ##     values[trims] = mapvalues(xs, from = uxs, to = moveto)
                ##     A.s = matrix(values,p,n)
                ##     all.na = which(is.na(A.s[,1]))
                ##     if(length(all.na) > 0) stop(paste("Entire rows are trimmed:", all.na))
                ##     ## Fit W with trimmed A and H. 
                ##     C.s = rbind(t(A.s), matrix(0, nrow = k, ncol = p))
                ##     W = t(Nnls(M,C.s))
                }else if(variation == "row" | variation == "cellrow"){
                    ## No change
                }else{
                    stop("Wrong mode. Try one of the following: 'cell', 'col', 'row', 'smooth'")
                }
            }
            ##J = nmlz(W) # Find the normalizing matrix of W
            ##W = W %*% J; H = solve(J) %*% H
            
            ## Convergence?
            #if(i > 1){
            #    if(setequal(to.trim[[i]], to.trim[[i-1]])){
            #        if(sum((W - W.prev)^2)/sum(W.prev^2) < tol) break
            #    }
            #}
            if(i > 1){
              
            }
            W.prev <- W
        }
    }
    ## Close the link to the progress bar. 
    if(showprogress) close(pb)
    
    fit = W %*% H
    if(!quiet){
        if(!trim){
            cat("Done. Time used:","\n")
            print(proc.time() - tic)
            cat("No trimming.\n",
                "Input matrix dimension: ",p," by ",n, "\n",
                "Left matrix: ",p," by ",k,". Right matrix: ",k," by ",n,"\n",
                "alpha = ",alpha,". beta = ",beta,"\n",
                "Number of max iterations = ",maxit,"\n",
                "Number of iterations = ",i,"\n", sep = ""
                )
        }else{
            cat("Done. Time used: ","\n")
            print(proc.time() - tic)
            cat("\n Trimming mode = \"", variation, "\". Proportion trimmed: ",trim, "\n",
                "Input matrix dimension: ",p," by ",n, "\n",
                "Left matrix: ",p," by ",k,". Right matrix: ",k," by ",n,"\n",
                "alpha = ",alpha,". beta = ",beta,"\n",
                "Number of max iterations = ",maxit,"\n",
                "Number of iterations = ",i,"\n", sep = ""
                )
        }
    }
    return(invisible(list(W = W, H = H, fit = fit,
                          trimmed = to.trim, niter = i)))
}

