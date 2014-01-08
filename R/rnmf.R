#' Performs robust penalized non-negative matrix factorization.
#' 
#' Performs robust penalized non-negative matrix factorizaion on a non-negative
#' matrix A to obtain W and H, such that A ~= W %*% H,  where A is a p by n matrix, W is a p by k matrix
#' and H is a k by n matrix. Outliers in A are detected and trimmed (by cells, columns or rows, see details).
#' Here we assume each row of A represents a feature/variable, and each column of A is an observation/sample.
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
#' set.seed(1)
#' Tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] = 1
#' ## Do rnmf with different settings
#' res.rnmf1 = rnmf(A = Tumor.corrupted, trim = FALSE, my.seed = 1)
#' res.rnmf2 = rnmf(A = Tumor.corrupted, tol = 0.001, trim = 0.06, my.seed = 1)
#' res.rnmf3 = rnmf(A = Tumor.corrupted, k = 10, beta = 0.1, tol = 0.001, trim = 0.06, my.seed = 1, variation = "smooth")
#' par(mfrow = c(2,2))
#' image(Tumor.corrupted, main = "Corrupted")
#' image(res.rnmf1$fit, main = "rnmf (no trimming) fit")
#' image(res.rnmf2$fit, main = "rnmf (smooth) fit 2")
#' image(res.rnmf3$fit, main = "rnmf (smooth) fit 3")

rnmf = function(A, k = 5, alpha = 0, beta = 0, maxit = 50, tol = 0.005,
    trim = FALSE, ini.W = NULL, ini.H = NULL,  ini.zeta = NULL, my.seed = NULL,
    variation = "cell", quiet = FALSE, nreg = 1, p1 = NA, n1 = NA)
{
    tic = proc.time() ## Start a clock.
    ## This function checks if the arguments are valid. If any of them are not valid, display an error message.
    checkargs(A, k, alpha, beta, maxit, tol, trim, ini.W, ini.zeta, my.seed, variation, quiet, nreg, p1, n1)
    p = nrow(A)
    n = ncol(A)
    ## Create a data frame contains value = cell values of A; (x,y) = coordinates; outs = is it an outlier?
    A.f = data.frame(value = as.vector(A), x = rep(1:p, n), y = rep(1:n, each = p), outs = FALSE) 
    if(trim > 0){
        to.trim1 = vector("list", maxit)
        ##to.trim2 = vector("list", maxit)
    }else{
        to.trim1 = NULL
        ##to.trim2 = NULL
    }
    
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
    ## The following line creates a progress bar.
    pb = txtProgressBar(min = 0, max = maxit, style = 3)
    obj = 1 ## 'obj' is a vector containing values of the objective function used in corresponding "variations".
    if(variation == "cell" & trim > 0){
        ## Initialize zeta
        ## zeta is a logic matrix which has the same size as A, indicating cells to be kept (non-outliers)
        zeta.allTRUE = matrix(TRUE, nrow = p, ncol = n)
        if(!missing(ini.zeta)){
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

        ## Start iterations.
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
                to.trim1[[i]] = which(rank(R) > round((1 - trim) * n * p))
                ## Update zeta
                zeta = zeta.allTRUE
                zeta[to.trim1[[i]]] = FALSE
                ## Refit W
                W = Nnls.trimW(H, A, zeta, alpha, p1, n1, k, p)
                J = nmlz(W) # Find the normalizing matrix of W
                W = W %*% J; H = solve(J) %*% H
            }
            obj[i] = l2((A - W %*% H)[zeta]) + l2(W) + sum(colSums(abs(H))^2)
            if(i > 1){
                if(all(to.trim1[[i]] == to.trim1[[i-1]]) & sum((W - W.prev)^2)/sum(W.prev^2) < tol) break
            }
            W.prev = W
        }
    }else{
        for(i in 1:maxit){
            flag = FALSE
            setTxtProgressBar(pb, i) # update the progress bar
            ## Iteration stage 1. Estimate H. ||MH - C||
            C = rbind(A, matrix(0,1,n))
            M = rbind(W, sqrt(beta) * matrix(1, 1, k))
            H = Nnls(M,C)
            if(trim > 0){
                R = (A - W %*% H)^2
                if(variation == "col" | variation == "cellcol"){
                    ## No change
                }else if(variation == "cellrow" | variation == "cellall"){
                    to.trim1[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    trim.row = unique(to.trim1[[i]] %% p)
                    A.trim = A[-trim.row,]
                    W.trim = W[-trim.row,]
                    ## Fit H with trimmed A and W. 
                    M.trim = rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
                    C.trim = rbind(A.trim, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M.trim,C.trim)
                }else if(variation == "smooth"){
                    to.trim1[[i]] = which(rank(abs(R)) > round((1 - trim) * n * p))
                    A.s = matrix(smoothing2(A.f, to.trim1[[i]], p,n, frame = TRUE)$value,p,n)
                    ## Fit H with trimmed A and W. 
                    C.s = rbind(A.s, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H = Nnls(M,C.s)
                }else if(variation == "row" | variation == "all"){
                    rsqr = apply(R, 1, l2) ## Find SS of residuals of rows
                    to.trim1[[i]] = which(rank(rsqr) > round((1 - trim) * p))
                    A.trim = A[-to.trim1[[i]],]
                    W.trim = W[-to.trim1[[i]],]
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
                    trims <- to.trim1[[i]] <- which(rank(abs(R)) > round((1 - trim) * n * p))
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

            ## Iteration stage 2. Estimate W. ||MW^T - C||.
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
        }
    }
    ## Close the link to the progress bar. 
    close(pb)
    
    fit = W %*% H
    if(!quiet){
        if(!trim){
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
    return(invisible(list(W = W, H = H, fit = fit, trimmed1 = to.trim1, trimmed2 = to.trim1, niter = i)))
}

