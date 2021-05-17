# Simple example with just few death rates and two time points

qx1 <- c(0.003, 0.005, 0.07, 1)
qx2 <- c(0.001, 0.002, 0.05, 1)

# Life expectancy difference

## Lifetable Function
lifetable <- function(qx){
  ax <- 0.5
  px <- 1-qx
  lx <- c(100000, cumprod(px)*100000)
  dx <- -diff(lx)
  Lx1 <- lx[c(-1, -length(lx))]+ax*dx[-length(dx)]
  Lx2 <- ax*dx[length(dx)]   
  Lx <- c(Lx1, Lx2)
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/lx[-length(lx)]
  out <- ex[1]
  return(out)
  
}

e1 <- lifetable(qx1)
e2 <- lifetable(qx2)

# Step 1: Average rates

mean_qx <- (qx1+qx2)/2
 
# Step 2: Stepwise calculation with replacement
# One vector with the "before" qx and one with the "after" qx at age I

before1 <- mean_qx
before1[1] <- qx1[1]

after1 <- mean_qx
after1[1] <- qx2[1]

# Calculate LE and the difference
ex_bef1 <- lifetable(before1)
ex_aft1 <- lifetable(after1)
diff1 <- ex_aft1-ex_bef1


# Next steps: do this for every age
# Age II
before2 <- mean_qx
before2[2] <- qx1[2]

after2 <- mean_qx
after2[2] <- qx2[2]

ex_bef2 <- lifetable(before2)
ex_aft2 <- lifetable(after2)
diff2 <- ex_aft2-ex_bef2

# Age III
before3 <- mean_qx
before3[3] <- qx1[3]

after3 <- mean_qx
after3[3] <- qx2[3]

ex_bef3 <- lifetable(before3)
ex_aft3 <- lifetable(after3)
diff3 <- ex_aft3-ex_bef3

# Age IV
before4 <- mean_qx
before4[4] <- qx1[4]

after4 <- mean_qx
after4[4] <- qx2[4]

ex_bef4 <- lifetable(before4)
ex_aft4 <- lifetable(after4)
diff4 <- ex_aft4-ex_bef4


# Age-specific contributions
tot_diff <- c(diff1, diff2, diff3, diff4)

sum(tot_diff)
e2-e1

