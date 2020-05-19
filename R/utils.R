
#' @title Make Genome Wide Predictions From a Bed File
#' @description
#' Create a score for each individual from weighted genetic data
#'  
#' @param bed an rbed object, or a mergedrbed object
#' @param weights a vector of length bed$no.snps for an rbed, or length(bed$DATkeep) for a mergedrbed object
#' @param ... extra parameters to \code{\link[pcapred]{get_data}}
#' @return a vector of length bed$no.inds
#' @export
bedpred=function(bed,weights,...){
    weights=matrix(weights,ncol=1)
    pred=lapply(1:ceiling(bed$no.ind/4),function(i){
        tdata=pcapred::get_data(bed,i,meanimpute=TRUE,verbose=FALSE,...)
        pred= tdata %*% weights
    })
    do.call("c",pred)
}

#' @title Mean Squared Error for Genomewide prediction
#' @description
#' Constructed the genomewide prediction for a score function defined by \code{weights} on an \code{rbed} object, then evaluate  the mean squared error to a provided \code{y}.
#'  
#' @param bed an rbed object, or a mergedrbed object
#' @param y vector of length bed$no.ind describing the truth
#' @param weights a vector of length bed$no.snps for an rbed, or length(bed$DATkeep) for a mergedrbed object
#' @param pred (default=NULL) a previous call to bedloss, to save calling bedpred again
#' @param ... extra parameters to \code{\link[pcapred]{get_data}}
#' @return A list containing pred: the predictions and loss: the mse of those predictions
#' @export
bedloss=function(bed,y,weights,pred=NULL,...){
    if(all(is.null(pred))) pred=bedpred(bed,weights,...)
    loss=mean((y-pred)^2)
    list(loss=loss,pred=pred)
}
