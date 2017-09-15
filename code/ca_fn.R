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
sdd_set_probs <- function(lc.mx, sdd.max, sdd.rate) {
  # Assign base dispersal probabilities from each cell
  # Each layer [1:i,1:j,,n] is the SDD neighborhood for cell n
  # k=1 contains pr(SDD | center,i,j)
  # k=2 contains the ID for each cell in the neighborhood
  # Returns array with dim(i:disp.rows, j:disp.cols, k:2, n:ncell)
  return(sdd.probs)
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





