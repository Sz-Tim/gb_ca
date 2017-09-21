# Differences between buckthorn CA and merow 2011:
#   - compositional habitats
#   - higher resolution (even @100 acre blocks)
#   - more nuanced dispersal: pr(drop) vs pr(bird)
#   - habitat specific K
#   - separate life history more -- fruit production step
#     - requires explicit mortality instead of lambdas < 1

##########----------------------------------------------------------------------
##### SIMULATION FUNCTIONS
##########----------------------------------------------------------------------

##---
## simulation wrapper
##---
  run_sim <- function(lc.df, N.init, K, fec, pr.f, pr.eat, sdd.pr, sdd.rate, 
                      n.ldd, pr.est, tmax, stoch=FALSE, simple=TRUE) {
    # Runs the simulation, calling subsequent submodules
    # simple=T runs model with no fruits/seeds/seedlings -- lambda only
    require(tidyverse); require(magrittr)
    
    # 1. Initialize populations
    ncell <- nrow(lc.df)
    N <- matrix(0, ncell, tmax+1)
    N[,1] <- N.init
    N.recruit <- rep(0, ncell)
    
    if(simple) {
      for(t in 1:tmax){
        # 2. Pre-multiply compositional parameters
        K.agg <- as.matrix(lc.df[,3:8]) %*% K
        lambda.agg <- as.matrix(lc.df[,3:8]) %*% lambda
        
        # 3. Local growth
        cat("Year", t, "- Growing...")
        N.new <- grow_pops(N[,t], lambda.agg, K.agg, stoch)
        
        # 4. Short distance dispersal
        cat("Dispersing locally...")
        N.emig <- sdd_simple(N[,t], N.new, sdd.pr, sdd.rate, K.agg, stoch)
        
        # 5. Long distance dispersal
        cat("Dispersing regionally...")
        N.emig <- ldd_disperse(ncell, N.emig, n.ldd, simple=TRUE)
        
        # 6. Update population sizes
        cat("Updating abundances.\n")
        N[,t+1] <- N.emig$N
      }
    } else {
      for(t in 1:tmax) {
        # 2. Pre-multiply compositional parameters
        lc.mx <- as.matrix(lc.df[,3:8])
        K.agg <- lc.mx %*% K
        fec.agg <- lc.mx %*% fec
        pr.f.agg <- lc.mx %*% pr.f
        pr.eat.agg <- lc.mx %*% pr.eat
        pr.est.agg <- lc.mx %*% pr.est
        
        # 3. Local fruit production
        cat("Year", t, "- Fruiting...")
        N.f <- make_fruits(N[,t], N.recruit, fec.agg, pr.f.agg, stoch)
        
        # 4. Short distance dispersal
        cat("Dispersing locally...")
        N.seed <- sdd_fs(N.f, pr.eat.agg, sdd.pr, sdd.rate, stoch)
        
        # 5. Long distance dispersal
        cat("Dispersing regionally...")
        N.seed <- ldd_disperse(ncell, N.seed, n.ldd, simple=FALSE)
        
        # 6. Seedling establishment
        cat("Establishing...")
        N.recruit <- new_seedlings(ncell, N.seed, pr.est.agg, stoch)
        
        # 7. Carrying capacity enforcement on adults
        cat("Hitting capacity.\n")
        N[,t] <- pmin(N[,t], ceiling(as.matrix(lc.df[,3:8]) %*% K))
        N[,t+1] <- N[,t] + N.recruit
      }
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
    bird.pref.agg <- as.matrix(lc.df[,3:8]) %*% bird.pref
    
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
    xx <- apply(lc.df, 1, function(x) seq(x[1]-sdd.max, x[1]+sdd.max))
    yy <- apply(lc.df, 1, function(x) seq(x[2]-sdd.max, x[2]+sdd.max))
    for(n in 1:(n.x*n.y)) {
      
      # vectorize for speed
      n_i <- expand.grid(xx[,n][xx[,n]>0 & xx[,n]<=n.x],
                         yy[,n][yy[,n]>0 & yy[,n]<=n.y])
      n_x <- apply(n_i, 1, function(x) which(xx[,n]==x[1]))
      n_y <- apply(n_i, 1, function(x) which(yy[,n]==x[2]))
      c_i <- apply(n_i, 1, function(x) which(lc.df[,1]==x[1] & lc.df[,2]==x[2]))
      
      # find cell ID for each cell in neighborhood
      for(i in 1:nrow(n_i)) {
        sdd.i[n_x[i],n_y[i],2,n] <- c_i[i]
      }
      
      # weight by bird habitat preference
      ib <- sdd.i[,,2,n] != 0  # inbounds neighbors
      sdd.i[,,1,n][ib] <- d.pr[ib] * bird.pref.agg[sdd.i[,,2,n][ib]]
      
      # set cell ID to 0 if pr(target) == 0 
      sdd.i[,,2,n][sdd.i[,,1,n]==0] <- 0
      if(n %% 100 == 0) {
        cat("finished cell", n, "\n")
      }
    }
    cat("finished cell", n)
    sdd.i[,,1,] <- apply(sdd.i[,,1,], 3, function(x) x/sum(x))
    return(sdd.i)
  }

  
  
  
##---
## local population growth: simple (a la Merow 2011)
##---
  grow_pops <- function(N.t, lambda.agg, K.agg, stoch=F) {
    # Calculate updated population sizes after local reproduction 
    # Growth rates are habitat specific
    # Returns sparse dataframe N.new with:
    #   col(id, N.pop, N.new=pop.delta, N.pop.upd=N.pop+N.new.nonemigrants)
    #   nrow = sum(N.new != 0)
    
    if(stoch) {
      
    } else {
      N.id <- which(N.t>0)
      lam.id <- lambda.agg[N.id]
      N.new <- tibble(id = which(N.t>0)) %>%
        mutate(N.pop=N.t[id],
               N.new=(N.pop * (lam.id-1)),
               N.pop.upd=pmin(K.agg[N.id,],
                              N.pop + 
                                (lam.id>=1)*N.new*pexp(0.5, sdd.rate) +
                                (lam.id<1)*N.new))
    }
    return(N.new)
  }
  
  
  
##---
## short distance dispersal: simple
##---
  sdd_simple <- function(N.t, N.new, sdd.pr, sdd.rate, K.agg, stoch=F) {
    # Calculate (N.arrivals | N.new, sdd.probs)
    # Accounts for distance from source cell & bird habitat preference
    # Returns dataframe with total population sizes.
    
    # assign emigrants to target cells
    N.source <- N.new %>% filter(N.new > 0)
    N.emig <- tibble(id=integer(), N=integer())
    for(i in 1:nrow(N.source)) {
      n <- N.source$id[i]
      N.emig %<>% add_row(id=c(sdd.pr[,,2,n]),
                          N=c(N.source$N.new[i] * sdd.pr[,,1,n]))
    }
    
    # sum within each target cell
    N.emig %<>% add_row(id=1:length(N.t), N=N.t) %>%
      group_by(id) %>% filter(id != 0) %>%
      summarise(N=sum(N))
    # restrict to K
    N.emig %<>% 
      mutate(N=pmin(K.agg[N.emig$id,], N))
    return(N.emig)
  }
  
  
  
##---
## local fruit production
##---
  make_fruits <- function(N.t, N.recruit, fec.agg, pr.f.agg, stoch=F) {
    # Calculate (N.fruit | N, fec, K) for each cell
    # Fecundity rates & fruiting probabilities are habitat specific
    # Assumes no fruit production in first year
    # N.fruit = (N-N.recruit)*pr(Fruit)*fec
    # Returns sparse dataframe N.f with:
    #   col(id, N.rpr=num.reproducing, N.fruit=total.fruit)
    #   nrow = sum(N.frt != 0)
    
    if(stoch) {
      N.f <- tibble(id = which((N.t-N.recruit)>0)) %>%
        mutate(N.rpr = rbinom(n(), N[id]-N.recruit[id],
                            prob=as.matrix(lc.df[id,3:8]) %*% pr.f),
               N.fruit = rpois(n(), 
                             lambda=as.matrix(lc.df[id,3:8]) %*% fec)) %>% 
        filter(N.fruit > 0)
    } else {
      N.f <- tibble(id = which((N.t-N.recruit)>0)) %>%
        mutate(N.rpr=((N.t[id]-N.recruit[id]) * pr.f.agg[id,]),
               N.fruit=(N.rpr * fec.agg[id,])) %>% 
        filter(N.fruit > 0)
    }
    return(N.f)
  }


  
##---
## short distance dispersal: fruits & seeds
##---
  sdd_fs <- function(N.f, pr.eat.agg, sdd.pr, sdd.rate, stoch=F) {
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
               N.emig = N.produced * (1-pexp(.5,sdd.rate)) * pr.eat.agg[id,],
               N.drop = N.produced - N.emig)
      
      # assign emigrants to target cells
      N.seed <- tibble(id=integer(), N.dep=integer())
      for(i in 1:sum(N.source$N.emig > 0)) {
        n <- N.source$id[i]
        N.seed %<>% add_row(id=c(sdd.pr[,,2,n]),
                             N.dep=c(N.source$N.emig[i] * sdd.pr[,,1,n]))
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
  ldd_disperse <- function(ncell, N.df, n.ldd, simple=TRUE) {
    # Assign n.ldd random LDD events 
    # A single seed is added to n.ldd target cells
    
    ldd.id <- sample(1:ncell, n.ldd, replace=TRUE)
    if(simple) {
      N.df$N[N.df$id==ldd.id] <- N.df$N[N.df$id==ldd.id] + 1
    } else {
      N.df %<>% add_row(id=ldd.id, N=rep(1, n.ldd)) %>%
        group_by(id) %>% summarise(N=sum(N))
    }
    
    return(N.df)
  }


  
##---
## seed germination & establishment
##---
  new_seedlings <- function(ncell, N.seed, pr.est.agg, stoch=F) {
    # Calculate (N.new | N.seed, pr.est)
    # Allows for incorporation of management effects
    
    if(stoch) {
      
    } else {
      N.recruit <- rep(0, ncell)
      N.recruit[N.seed$id] <- N.seed$N * pr.est.agg[N.seed$id,]
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







