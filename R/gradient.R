#' @title Generate a random initial parameter vector
#' @description
#' Uniform sampling for initial values of a parameter
#' @param len Length of vector
#' @param minTheta Minimum of uniform sampling
#' @param maxTheta Maximum of uniform sampling
#' @param seed (default NULL)
#' @return A matrix with 1 row and len columns
#' @export
initTheta <- function(len, minTheta=-1, maxTheta=1, seed=NULL){
        #create static random
        set.seed(seed)
        #random a value
        thetaList <- stats::runif(len, min=minTheta, max=maxTheta)
        #clear static random
        set.seed(seed)
        #transform into matrix
        result <- matrix(unlist(thetaList), ncol=len, nrow=1, byrow=FALSE)
        return(result)
}

#' @title Use the ADAM optimiser to learn a linear regression model from plink data
#' @description
#' For \eqn{y = \beta X + \epsilon}, use the ADAM optimiser to minimise
#' \deqn{|\epsilon|^2_2= \sum{i=1}^N (y_i - \sum_{j=1}^L \beta_{j} X_{ij})^2.}
#' The ADAM optimiser is an ADAptive Maximiser similar in spirit to Stochastic Gradient Descent. This allows linear regression to be applied to very large problems for which exact solutions are impossible.
#'
#' The components of this model are (p below):
#' \itemize{
#' \item theta (if missing, generate randomly using \code{\link{initTheta}}) our estimate of \eqn{\beta} of length L.
#' \item alpha (default 0.1) The "learning rate"
#' \item beta1 (default 0.9) A component of ADAM
#' \item beta1 (default 0.999) A component of ADAM
#' \item meanMoment (default a vector of 0 with length L) The current momentum for each variable
#' \item varianceMoment (default a vector of 0 with length L) The current momentum for each variables variance
#' \item smooth (default 1e-7) A regularisation on the variance estimate
#' }
#' 
#' @param bed A \code{rbed} object as returned by \code{\link[pcapred]{readbed}}, or a \code{mergedrbed} object as returned by \code{\link[pcapred]{mergeref}}.
#' @param y A vector of length bed$no.ind, giving the predicted value
#' @param maxIter (default=1000) maximum number of iterations to run the optimiser for
#' @param p (default =list()) a vector specifying any control parameters to be set manually; these are described above.
#' @param ssize (default=100) Number of 4-individual binary chunks to use per stochastic update; this is the "minibatch" size.
#' @param seed (default=NULL) Optional seed setting
#' @param iterskip (dfefault=10) If verbose, print to screen only on every iterskip iteration.
#' @param verbose (default=TRUE) whether to display iteration progress
#' @return A list containing:
#' \itemize{
#' \item loss: The history of the average loss (length maxIter), on the subset of data it was calculated on
#' \item change: The history of change size (length maxIter), i.e.  the RMS of the step sizes
#' \item theta: The current best estimate of \eqn{\beta}.
#' \item p: The previous parameter estimate.
#' }
#' @export
bedADAM=function (bed,y,maxIter = 1000,p=list(), 
                  ssize=100, seed = NULL,iterskip=10,verbose=TRUE) 
{
    set.seed(seed)

    if(is.null(p$theta)) p$theta <-initTheta(bed$no.snp, seed = seed)
    if(is.null(p$alpha)) p$alpha <-0.1
    if(is.null(p$beta1)) p$beta1  <- 0.9 
    if(is.null(p$beta2)) p$beta2 <- 0.999
    if(is.null(p$meanMoment)) p$meanMoment <- rep(0,bed$no.snp)
    if(is.null(p$varianceMoment)) p$varianceMoment <- rep(0,bed$no.snp)
    if(is.null(p$smooth)) p$smooth <- 1e-07

    temporaryTheta <- matrix(ncol = length(p$theta), nrow = 1)
    updateRule <- matrix(0, ncol = length(p$theta), nrow = 1)
    lossres=rep(NA,maxIter)
    thetares=rep(NA,maxIter)
    gradientres=rep(NA,maxIter)
    
    best=p$theta
    bestloss=Inf
    for (iteration in 1:maxIter) {
        samplei=sample(1:(bed$no.ind/4),ssize,replace=TRUE)
        sampleX=pcapred::get_data(bed,samplei,meanimpute=TRUE,verbose=FALSE)
        sampleY=y[attr(sampleX,"idx")]
        N=dim(sampleX)[1]
        error <- t(sampleY) - t(sampleX %*% t(p$theta)) 
        lossres[iteration]=sqrt(mean(error^2))
        if(verbose && iteration%%iterskip==0) cat("Iteration ",iteration," Loss = ",lossres[iteration]/N,"\n")
        gradient <- -(2/N) %*% error %*% sampleX

        p$meanMoment <- (p$beta1 * p$meanMoment) + (1 - p$beta1) * gradient
        p$varianceMoment <- (p$beta2 * p$varianceMoment) + (1 - p$beta2) * (gradient^2)
        mean.hat <- p$meanMoment/(1 - p$beta1)
        variance.hat <- p$varianceMoment/(1 - p$beta2)
        updateRule[1, ] <- (p$alpha/(sqrt(variance.hat) + p$smooth)) * mean.hat
        temporaryTheta[1, ] = p$theta[1, ] - updateRule[1,]
        
        thetares[iteration]=sqrt(mean((p$theta-temporaryTheta)^2))
        p$theta <- temporaryTheta
        ## Storing of current best estimate
        if(!is.finite(bestloss) || lossres[iteration]<bestloss) {
            best=p$theta
            bestloss=lossres[iteration]
        }
    }
    return(list(loss=lossres,
                change=thetares,
                coeff=best,
                p=p
                ))
}


#' @title Use the SGD optimiser to learn a linear regression model from plink data
#' @description
#' For \eqn{y = \beta X + \epsilon}, use the SGD optimiser to minimise
#' \deqn{|\epsilon|^2_2= \sum{i=1}^N (y_i - \sum_{j=1}^L \beta_{j} X_{ij})^2.}
#' The SGD (Stochastic Gradient Descent) approach allows linear regression to be applied to very large problems for which exact solutions are impossible.
#'
#' It is very good for large data problems but the learning rate needs to be carefully chosen.
#' 
#' The components of this model are (p below):
#' \itemize{
#' \item theta (if missing, generate randomly using \code{\link{initTheta}}) our estimate of \eqn{\beta} of length L.
#' \item alpha (default 0.1/bed$no.snp) The "learning rate". The algorithm will diverge if this is too large.
#' }
#' 
#' @param bed A \code{rbed} object as returned by \code{\link[pcapred]{readbed}}, or a \code{mergedrbed} object as returned by \code{\link[pcapred]{mergeref}}.
#' @param y A vector of length bed$no.ind, giving the predicted value
#' @param maxIter (default=1000) maximum number of iterations to run the optimiser for
#' @param p (default =list()) a vector specifying any control parameters to be set manually; these are described above.
#' @param ssize (default=100) Number of 4-individual binary chunks to use per stochastic update; this is the "minibatch" size.
#' @param seed (default=NULL) Optional seed setting
#' @param iterskip (dfefault=10) If verbose, print to screen only on every iterskip iteration.
#' @param verbose (default=TRUE) whether to display iteration progress
#' @return A list containing:
#' \itemize{
#' \item loss: The history of the average loss (length maxIter), on the subset of data it was calculated on
#' \item change: The history of change size (length maxIter), i.e.  the RMS of the step sizes
#' \item theta: The current best estimate of \eqn{\beta}.
#' \item p: The previous parameter estimate.
#' }
#' @export
bedSGD=function (bed,y, maxIter = 1000, p=list(),
                  ssize=100, seed = NULL,verbose=TRUE,iterskip=10) 
{
    set.seed(seed)
    if(is.null(p$theta)) p$theta <-initTheta(bed$no.snp, seed = seed)
    if(is.null(p$alpha)) p$alpha <- 1/length(p$theta)
    temporaryTheta <- matrix(ncol = length(p$theta), nrow = 1)
    lossres=rep(NA,maxIter)
    thetares=rep(NA,maxIter)
    gradientres=rep(NA,maxIter)
    
    best=p$theta
    bestloss=Inf
    for (iteration in 1:maxIter) {
        samplei=sample(1:(bed$no.ind/4),ssize,replace=TRUE)
        sampleX=pcapred::get_data(bed,samplei,meanimpute=TRUE,verbose=FALSE)
        sampleY=y[attr(sampleX,"idx")]
        N=dim(sampleX)[1]
        error <- t(sampleY) - t(sampleX %*% t(p$theta)) 
        lossres[iteration]=sqrt(mean(error^2))
        if(verbose && iteration%%iterskip==0) cat("Iteration ",iteration," Loss = ",lossres[iteration]/N,"\n")
        gradient <- -(2/N) %*% error %*% sampleX
        temporaryTheta=p$theta - p$alpha * gradient

        thetares[iteration]=sqrt(mean((p$theta-temporaryTheta)^2))
        p$theta <- temporaryTheta
        if(!is.finite(bestloss) || lossres[iteration]<bestloss) {
            best=p$theta
            bestloss=lossres[iteration]
        }
    }
    return(list(loss=lossres,
                change=thetares,
                coeff=best,
                p=p
                ))
}
