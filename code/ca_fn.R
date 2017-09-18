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
  run_sim <- function(lc.df, N.init, K, fec, pr.f, pr.eat, sdd.pr, sdd.rate, 
                      n.ldd, pr.est, tmax, stoch=FALSE) {
    # Runs the simulation, calling subsequent submodules
    require(tidyverse); require(magrittr)
    
    # 1. Initialize populations
    ncell <- nrow(lc.df)
    N <- matrix(0, ncell, tmax+1)
    N[,1] <- N.init
    N.recruit <- rep(0, ncell)
    
    for(t in 1:tmax) {
      # 2. Local fruit production
      N.f <- make_fruits(lc.df, N[,t], N.recruit, fec, pr.f, stoch=stoch)
      
      # 3. Short distance dispersal
      N.seed <- sdd_disperse(lc.df, N.f, pr.eat, sdd.pr, sdd.rate, stoch)
      
      # 4. Long distance dispersal
      N.seed <- ldd_disperse(lc.df, N.seed, n.ldd)
      
      # 5. Seedling establishment
      N.recruit <- new_seedlings(lc.df, N.seed, pr.est, stoch)
      
      # 6. Carrying capacity enforcement on adults
      N[,t] <- pmin(N[,t], ceiling(as.matrix(lc.df[,3:8]) %*% K))
      N[,t+1] <- N[,t] + N.recruit
      
      # progress
      cat("Finished year", t, "\n")
    }
    return(N)
  }


  
##---
## short distance dispersal probabilities
##---
  sdd_set_probs <- function(lc.df, sdd.max, sdd.rate, bird.pref, trunc.diag=T) {
    # Assign base dispersal probabilities from each cell
    # Each layer [1:i,1:j,,n] is the SDD neighborhood for cell n
    # trunc.diag: if TRUE, the sdd neighborhood is restricted to within sdd.max
    #   including along diagonals; if FALSE, then sdd.max refers to the maximum
    #   allowable horizontal & vertical distance from the source cell
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
        if(trunc.diag) {
          if( sqrt((i-ctr)^2 + (j-ctr)^2) > sdd.max ) {
            d.pr[i,j] <- 0
          }
        }
        d.pr[j,i] <- d.pr[i,j]
      }
    }
    d.pr <- d.pr/sum(d.pr)
    
    # pair cell IDs for each neighborhood
    xx <- apply(lc.df, 1, function(x) seq(x[1]-sdd.max, x[1]+sdd.max)) %>% t
    yy <- apply(lc.df, 1, function(x) seq(x[2]-sdd.max, x[2]+sdd.max)) %>% t
    for(n in 1:(n.x*n.y)) {
      for(i in xx[n,][xx[n,]>0 & xx[n,]<=n.x]) {
        for(j in yy[n,][yy[n,]>0 & yy[n,]<=n.y]) {
          sdd.i[xx[n,]==i, yy[n,]==j, 2, n] <- which(lc.df[,1]==i &
                                                       lc.df[,2]==j)
        }
      }
      # weight by bird habitat preference
      ib <- sdd.i[,,2,n] != 0  # inbound neighbors
      sdd.i[,,1,n][ib] <- d.pr[ib] * 
        (as.matrix(lc.df[sdd.i[,,2,n][ib], 3:8]) %*% bird.pref)
      sdd.i[,,1,n] <- sdd.i[,,1,n]/sum(sdd.i[,,1,n])
    }
    return(sdd.i)
  }


  
##---
## local population growth: simple (a la Merow 2011)
##---
  grow_pops <- function(lc.df, N.t, lambda, stoch=F) {
    # Calculate updated population sizes after local reproduction 
    # Growth rates are habitat specific
    # Returns sparse dataframe N.new with:
    #   col(id, N.pop, N.new=pop.delta, N.pop.upd=N.pop+N.new.nonemigrants)
    #   nrow = sum(N.new != 0)
    
    if(stoch) {
      
    } else {
      N.id <- which(N.t>0)
      K.id <- as.matrix(lc.df[N.id, 3:8]) %*% K
      lam.id <- as.matrix(lc.df[N.id, 3:8]) %*% lambda
      N.new <- tibble(id = which(N.t>0)) %>%
        mutate(N.pop=N.t[id],
               N.new=(N.pop * (lam.id-1)) %>% ceiling,
               N.pop.upd=pmin(K.id,
                              N.pop + 
                                (lam.id>=1)*N.new*pexp(0.5, sdd.rate) +
                                (lam.id<1)*N.new) %>% ceiling)
    }
    return(N.new)
  }
  
  
  
##---
## local fruit production
##---
  make_fruits <- function(lc.df, N.t, N.recruit, fec, pr.f, stoch=F) {
    # Calculate (N.fruit | N, fec, K) for each cell
    # Fecundity rates & fruiting probabilities are habitat specific
    # Assumes no fruit production in first year
    # N.fruit = (N-N.recruit)*pr(Fruit)*fec
    # Returns sparse dataframe N.f with:
    #   col(id, N.rpr=num.reproducing, N.fruit=total.fruit)
    #   nrow = sum(N.frt != 0)
    
    if(stoch) {
      N.f <- tibble(id = which(N.t>0)) %>%
        mutate(N.rpr = rbinom(n(), N[id]-N.recruit[id],
                            prob=as.matrix(lc.df[id,3:8]) %*% pr.f),
               N.fruit = rpois(n(), 
                             lambda=as.matrix(lc.df[id,3:8]) %*% fec)) %>% 
        filter(N.fruit > 0)
    } else {
      N.f <- tibble(id = which(N.t>0)) %>%
        mutate(N.rpr = ((N.t[id]-N.recruit[id]) * 
                          as.matrix(lc.df[id,3:8]) %*% pr.f) %>% ceiling,
               N.fruit = (N.rpr * 
                            as.matrix(lc.df[id,3:8]) %*% fec) %>% ceiling) %>% 
        filter(N.fruit > 0)
    }
    return(N.f)
  }


  
##---
## short distance dispersal
##---
  sdd_disperse <- function(lc.df, N.f, pr.eat, sdd.pr, sdd.rate, stoch=F) {
    # Calculate (N.seeds | N.fruit, sdd.probs, pr.eaten)
    # Accounts for distance from source cell, bird habitat preference,
    #   and the proportion of fruits eaten vs dropped
    # N.emig excludes seeds deposited within the source cell (1-pexp(0.5, rate))
    # Returns sparse matrix N.sdd with:
    #   cols(cell.ID, (N.seed = 2*(N.fruit - N.eaten + N.deposited)))
    #   nrow = sum(N.seed != 0)
    
    if(stoch) {
      
    } else {
      # calculate seeds deposited within source cell vs emigrants
      N.source <- N.f %>%
        mutate(N.produced = 2*N.fruit,
               N.emig = (N.produced * (1-pexp(.5,sdd.rate)) *
                          as.matrix(lc.df[id,3:8]) %*% pr.f) %>% ceiling,
               N.drop = N.produced - N.emig)
      
      # assign emigrants to target cells
      N.seed <- tibble(id=integer(), N.dep=integer())
      for(i in 1:sum(N.source$N.emig > 0)) {
        n <- N.source$id[i]
        N.seed %<>% add_row(id=c(sdd.pr[,,2,n]),
                             N.dep=c(N.source$N.emig[i] * sdd.pr[,,1,n]) %>% 
                               ceiling)
      }
      
      # sum within each target cell
      N.seed %<>% add_row(id=N.source$id, N.dep=c(N.source$N.drop)) %>%
        group_by(id) %>% 
        summarise(N=sum(N.dep)) %>% 
        filter(N > 0)
    }
    return(N.seed)
  }


  
##---
## long distance dispersal
##---
  ldd_disperse <- function(lc.df, N.seed, n.ldd) {
    # Assign n.ldd random LDD events 
    # A single seed is added to n.ldd target cells
    
    ldd.id <- sample(1:nrow(lc.df), n.ldd, replace=TRUE)
    N.seed %<>% add_row(id=ldd.id, N=rep(1, n.ldd)) %>%
      group_by(id) %>% summarise(N=sum(N))
    
    return(N.seed)
  }


  
##---
## seed germination & establishment
##---
  new_seedlings <- function(lc.df, N.seed, pr.est, stoch=F) {
    # Calculate (N.new | N.seed, pr.est)
    # Allows for incorporation of management effects
    
    if(stoch) {
      
    } else {
      N.recruit <- rep(0, nrow(lc.df))
      N.recruit[N.seed$id] <- (N.seed$N * 
                      as.matrix(lc.df[N.seed$id,3:8]) %*% pr.est) %>% ceiling
    }
    return(N.recruit)
  }







##########----------------------------------------------------------------------
##### MUNGING & PLOTTING FUNCTIONS
##########----------------------------------------------------------------------

##---
## add block IDs for aggregating acres
##---
add_blocks <- function(x, cb.i=cb.i) {
  # adds Block IDs to cells for aggregating 1-acre cells to larger blocks
  # cb.i is a reference dataframe identifying which cells belong to which blocks
  require(magrittr)
  x %>%
    mutate(BlockID=cb.i$BlockID[match(.$CellID, cb.i$CellID)]) %>%
    filter(!is.na(BlockID)) %>% 
    group_by(BlockID)
}







