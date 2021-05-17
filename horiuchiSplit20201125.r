## ## Fct for Horiuchi Decomposition using 4 different subpopulations

## endY <- dat[,"2000", "f",,]
## staY <- dat[,"1999", "f",,]

## dimnames(new)

horiuchi.decompo.split4 <- function(endY,staY){

    ## Functions 
    lifetable <- function(mx){
        ax <- c(0.1, rep(0.5, length(mx)-1))
        qx <- mx/(1+(1-ax)*mx)
        qx[length(qx)] <- 1
        px <- 1-qx
        lx <- c(100000, cumprod(px)*100000)
        dx <- -diff(lx)
        Lx1 <- lx[c(-1, -length(lx))]+ax[-length(ax)]*dx[-length(dx)]
        Lx2 <- ifelse(mx[length(mx)] == 0, Lx1[length(Lx1)], dx[length(dx)]/mx[length(mx)])   
        Lx <- c(Lx1, Lx2)
        Tx <- rev(cumsum(rev(Lx)))
        ex <- Tx/lx[-length(lx)]
        out <- ex[1]
        return(out)
        
    }

    ## averages as input for decompo
    average.prop <- (endY[,,"pr"]+staY[,,"pr"])/2
    average.mx <- (endY[,,"mx"]+staY[,,"mx"])/2
    averageTot <- average.prop*average.mx
    average.rates <- apply(averageTot, MARGIN = 1, FUN = sum)

    ## ## Preparation of the basic arrays for the analysis
    rep.col<-function(x,n){
        matrix(rep(x,each=n), ncol=n, byrow=TRUE) 
    }

    ## Array to save results

    outArray <- array(NA, dim=dim(endY), dimnames=dimnames(endY))
    allEth <- c("b", "h", "o", "w")

    ## Death rate cont (I)

    for(i in allEth){
        ethDec <- i
        ethOther <- allEth[!allEth %in% i] 

        rates.begin.I <- rates.end.I <- rep.col(average.rates,
                                                n = length(average.rates))
        beginning.rates.I <- rowSums(averageTot[, ethOther])+average.prop[,ethDec]*staY[,ethDec,"mx"]
        end.rates.I <- rowSums(averageTot[, ethOther])+average.prop[,ethDec]*endY[,ethDec,"mx"]

        diag(rates.begin.I) <- beginning.rates.I
        diag(rates.end.I) <- end.rates.I
        
        ## ## Life Expectancy
        e0.begin.I <- apply(rates.begin.I, MARGIN = 2, FUN = lifetable)
        e0.end.I <- apply(rates.end.I, MARGIN = 2, FUN = lifetable)
        outArray[,ethDec,"mx"] <- e0.end.I-e0.begin.I
    }

    ## Proportion (II)

    for(i in allEth){
        ethDec <- i
        ethOther <- allEth[!allEth %in% i] 

        rates.begin.I <- rates.end.I <- rep.col(average.rates,
                                                n = length(average.rates))
        beginning.rates.I <- rowSums(averageTot[, ethOther])+average.mx[,ethDec]*staY[,ethDec,"pr"]
        end.rates.I <- rowSums(averageTot[, ethOther])+average.mx[,ethDec]*endY[,ethDec,"pr"]

        diag(rates.begin.I) <- beginning.rates.I
        diag(rates.end.I) <- end.rates.I
        
        ## ## Life Expectancy
        e0.begin.I <- apply(rates.begin.I, MARGIN = 2, FUN = lifetable)
        e0.end.I <- apply(rates.end.I, MARGIN = 2, FUN = lifetable)
        outArray[,ethDec,"pr"] <- e0.end.I-e0.begin.I
    }

    return(outArray)
}
