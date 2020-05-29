#' @title Solve an SVD for loadings
#' @description Take a data matrix \eqn{X}, its corresponding eigenvectors \eqn{U} and eigenvalues \eqn{D}, and return an estimate of \eqn{V}, where \eqn{X= UDV'}.
#'
#' You might want to do this if: a) you obtained PCA data that did not contain loadings, b) you wish to use a SNP set that is different to the one used to construct the PCA.
#'
#' This can be dangerous when used for extrapolation, so make sure you test it carefully.
#' @param X a data matrix of dimension N by L
#' @param U the Eigenvector matrix of dimension N by K
#' @param D the Eigenvalue vector of length K
#' @return new loadings in the form of an L by K matrix
#' @export
get_loadings <- function(X,U,D) {
  return( 1/D %*% t(U) %*% X )
}

get_predictions <- function(bed,loadings,verbose=FALSE){
    gd=function(bed,i,y,...){
        data=pcapred::get_data(bed,i,meanimpute=TRUE,verbose=verbose)
        ret=as.numeric(data %*% loadings)
        attr(ret, "idx")=attr(data, "idx")
        ret
    }
    i=1
    pred=rep(NA,bed$no.ind)
    while(i<=ceiling(bed$no.ind/4)){
        tmp=gd(bed,i,y,verbose=verbose)
        pred[attr(tmp,"idx")]=tmp
        i=i+1
    }
    pred
}

bedbiglm<-function(bed,y,verbose=TRUE){
    gd=function(bed,i,y,...){
        data=pcapred::get_data(bed,i,verbose=verbose)
        colnames(data)=paste0("SNP",1:dim(data)[2])
        ret=as.data.frame(cbind(y=y[attr(data,"idx")],data))
        attr(ret,"idx")=attr(data,"idx")
        ret
    }
    if (!requireNamespace("biglm", quietly = TRUE)){
        stop("Package \"biglm\" needed for this function to work. Please install it with install.packages(\"biglm\".",
             call. = FALSE)
    }
    data=gd(bed,1,y,verbose=verbose)
    formula=stats::as.formula(paste0("y~",
                              paste0(colnames(data)[-1],collapse="+"),
                              "-1"
                              ))
    mybiglm=biglm::biglm(formula,data=data)
    i=2
    while(i<=ceiling(bed$no.ind/4)){
        data=gd(bed,i,y,verbose=verbose)
        mybiglm=biglm::update(mybiglm,data)
        i=i+1
    }
    list(biglm=mybiglm,
         coeff=summary(mybiglm)$mat[,1])
}

bedbigfastlm<-function(bed,y,tfile=NULL,verbose=TRUE,...){
    if (!requireNamespace("bigmemory", quietly = TRUE)){
        stop("Package \"bigmemory\" needed for this function to work. Please install it with install.packages(\"bigmemory\".",
             call. = FALSE)
    }
    if (!requireNamespace("bigFastlm", quietly = TRUE)){
        stop("Package \"bigFastlm\" needed for this function to work. Please install it with remotes::install.github(\"jaredhuling/bigFastlm\".",
             call. = FALSE)
    }
    if(verbose) cat("Converting to bigmatrix format\n")
    bigmat=pcapred::bedasbigmatrix(bed,tfile,verbose=verbose)
    if(verbose) cat("Performing computation\n")
    mybigfastlm <- bigFastlm::bigLm(bigmat, y,...)
    ret=list(bigfastlm=mybigfastlm,
             bigmat=bigmat,
             coeff=mybigfastlm$coeff
             )
    return(ret)
}

