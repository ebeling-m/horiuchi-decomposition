library(tidyverse)
# Usinng the Horiuchi et al algorithm to decompose life expectancy differences into age-specific contributions
# Death rates

mx <- read.table("GBRTENW.Mx_1x1.txt", header = TRUE, skip = 2, stringsAsFactors = FALSE)
mx$Age[mx$Age == "110+"] <- 110
mx$Age <- as.numeric(mx$Age)

mxSel <- 
  mx %>% 
  filter(Year > 2010) %>% 
  mutate(Female = as.numeric(Female))

mxAr <- tapply(mxSel$Female, INDEX = list(mxSel$Age, mxSel$Year), FUN = sum)

# Decompose LE difference between 2012 and 2013
## Lifetable Function
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

mx12 <- mxAr[,"2012"]
mx13 <- mxAr[,"2013"]

e12 <- lifetable(mx12)
e13 <- lifetable(mx13)

# Step 1: calculate average death rates
mean_mx <- (mx12+mx13)/2

# Step 2: Calculate before (2012) and after (2013) life expectancies
before <- after <- matrix(mean_mx, nrow = length(mean_mx), ncol = length(mean_mx), byrow = FALSE)
diag(before) <- mx12
diag(after) <- mx13

e_bef <- apply(before, MARGIN = 2, FUN = lifetable)
e_aft <- apply(after, MARGIN = 2, FUN = lifetable)

e_diff <- e_aft-e_bef

e13-e12
sum(e_diff)







