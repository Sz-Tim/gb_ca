# Differences between buckthorn CA and merow 2011:
#   - compositional habitats
#   - higher resolution (even @100 acre blocks)
#   - more nuanced dispersal: pr(drop) vs pr(bird)
#   - universal K or habitat specific?
#   - separate life history more -- fruit production step
#     - requires explicit mortality instead of lambdas < 1

##########----------------------------------------------------------------------
##### SIMULATION FUNCTIONS
##########----------------------------------------------------------------------

##---
## simulation wrapper
##---
run_sim <- function(
  lc.mx,
  N.init,
  lambda,
  K,
  sdd.probs,
  sdd.rate,
  bird.pref,
  pr.eat,
  n.ldd,
  tmax
) {
  # Runs the simulation, calling subsequent submodules
  # 1. Initialize populations
  N[,1] <- N.init
  for(t in 2:tmax) {
    # 2. Local fruit production
    N.f <- make_fruits(lc.mx, N[,t], N.recruit[,t-1], fec, K)
    # 3. Short distance dispersal
    N.seed <- sdd_disperse(lc.mx, N.f, sdd.probs, bird.pref, pr.eat)
    # 4. Long distance dispersal
    N.seed <- ldd_disperse(lc.mx, N[,t], N.seed, n.ldd)
    # 5. Seedling establishment
    N.recruit[,t] <- new_seedlings(lc.mx, N.seed, pr.est)
    # 6. Carrying capacity enforcement on adults
    N[,t] <- pmin(N[,t], K)
    N[,t+1] <- N[,t] + N.recruit[,t]
  }
  return(N)
}


##---
## short distance dispersal probabilities
##---
  sdd_set_probs <- function(lc.df, sdd.max, sdd.rate, bird.pref) {
    # Assign base dispersal probabilities from each cell
    # Each layer [1:i,1:j,,n] is the SDD neighborhood for cell n
    # k=1 contains pr(SDD | center,i,j)
    # k=2 contains the ID for each cell in the neighborhood
    # Returns array with dim(i:disp.rows, j:disp.cols, k:2, n:ncell)
    
    n.x <- max(lc.df[,1])
    n.y <- max(lc.df[,2])
    nbr <- 2 * sdd.max + 1
    sdd.i <- array(0, dim=c(nbr, nbr, 2, n.x*n.y))
    
    # generate default dispersal probability matrix
    d.pr <- matrix(0, nbr, nbr)
    ctr <- sdd.max + 1  # center index (i=j) for square mx
    for(i in 1:nbr) {
      for(j in i:nbr) {
        d.pr[i,j] <- dexp((i-ctr)^2 + (j-ctr)^2 - 0.5, sdd.rate)
        d.pr[j,i] <- d.pr[i,j]
      }
    }
    d.pr <- d.pr/sum(d.pr)
    
    # pair cell IDs for each neighborhood
    xx <- apply(lc.df, 1, function(x) seq(x[1]-sdd.max, x[1]+sdd.max)) %>% t
    yy <- apply(lc.df, 1, function(x) seq(x[1]-sdd.max, x[1]+sdd.max)) %>% t
    for(n in 1:(n.x*n.y)) {
      for(i in xx[n,][xx[n,]>0 & xx[n,]<=n.x]) {
        for(j in yy[n,][yy[n,]>0 & yy[n,]<=n.y]) {
          sdd.i[xx[n,]==i, yy[n,]==j, 2, n] <- which(lc.df[,1]==i &
                                                       lc.df[,2]==j)
        }
      }
      # weight by bird habitat preference
      ib <- sdd.i[,,2,n] != 0  # inbound neighbors
      sdd.i[,,1,n][ib] <- d.pr[ib] * bird.pref[lc.df[sdd.i[,,2,n][ib], 3]]
      sdd.i[,,1,n] <- sdd.i[,,1,n]/sum(sdd.i[,,1,n])
    }
    
    return(sdd.i)
  }


##---
## local fruit production
##---
make_fruits <- function(lc.mx, N, N.recruit, fec, K) {
  # Calculate (N.fruit | N, fec, K) for each cell
  # Fecundity rates & fruiting probabilities are habitat specific
  # Assumes no fruit production in first year
  # N.fruit = (N-N.recruit)*pr(Fruit)*fec
  # Returns sparse matrix N.f with:
  #   col(cell.ID, N.fruit)
  #   nrow = sum(N.fruit != 0)
  return(N.f)
}


##---
## short distance dispersal
##---
sdd_disperse <- function(lc.mx, N.f, sdd.probs, bird.pref, pr.eat) {
  # Calculate (N.seeds | N.fruit, sdd.probs, bird.pref, pr.eaten)
  # Accounts for distance from source cell, bird habitat preference,
  #   and the proportion of fruits eaten vs dropped
  # Returns sparse matrix N.sdd with:
  #   cols(cell.ID, (N.seed = 2*(N.fruit - N.eaten + N.deposited)))
  #   nrow = sum(N.seed != 0)
  return(N.seed)
}


##---
## long distance dispersal
##---
ldd_disperse <- function(lc.mx, N, N.seed, n.ldd) {
  # Assign n.ldd random LDD events 
  # A single seed is added to n.ldd target cells
  return(N.seed)
}


##--
## seed germination & establishment
##--
new_seedlings <- function(lc.mx, N.seed, pr.est) {
  # Calculate (N.new | N.seed, pr.est)
  # Allows for incorporation of management effects
  return(N.recruit)
}





