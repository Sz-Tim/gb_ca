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
  run_sim <- function(g.p, lc.df, sdd.pr, N.init, control.p=NULL) {
    # Runs the simulation, calling subsequent submodules
    # simple=T runs model with no fruits/seeds/seedlings -- lambda only
    # Storage Structures:
    #   N.lam: matrix(ncell, tmax+1) -- lambda model 
    #   N: array(ncell, tmax+1, 8) -- hybrid model; age 1-7, adults
    #   N.ling: vector(ncell) -- hybrid model; number of new seedlings for t+1
    #   N.sb: vector(ncell) -- hybrid model; number of seeds in seedbank
    require(tidyverse); require(magrittr)
    
    # Unpack parameters
    tmax <- g.p$tmax
    simple <- g.p$simple
    dem.st <- g.p$dem.st
    sdd.st <- g.p$sdd.st
    mat.d <- g.p$mat.d
    bank <- g.p$bank
    ncell <- g.p$lc.r * g.p$lc.c
    K <- g.p$K  # carrying capacity
    pr.s <- g.p$pr.s  # pre-adult survival
    pr.f <- g.p$pr.f  # pr(fruit)
    fec <- g.p$fec  # mean(fruit per adult)
    age.mat <- g.p$age.mat  # mean age at first fruiting
    pr.est <- g.p$pr.est  # pr(seedling est)
    pr.sb <- g.p$pr.sb  # pr(ann.surv seed bank)
    lambda <- g.p$lambda  # pop growth rate
    sdd.rate <- g.p$sdd.rate  # 1/mn for dispersal kernel
    pr.eat <- g.p$pr.eat  # pr(birds eat frt)
    n.ldd <- g.p$n.ldd   # num long distance dispersal events per year
    n.class <- max(g.p$age.mat)
    
    # If buckthorn is being actively managed...
    if(!is.null(control.p)) {}
    
    if(simple) {
      # 1. Initialize populations
      N.lam <- matrix(0, ncell, tmax+1)  
      N.lam[,1] <- rowSums(N.init)
      
      for(t in 1:tmax){
        # 2. Pre-multiply compositional parameters
        K.agg <- as.matrix(lc.df[,4:9]) %*% K
        lambda.agg <- as.matrix(lc.df[,4:9]) %*% lambda
        
        # 3. Local growth
        cat("Year", t, "- Growing...")
        N.new <- grow_pops(N.lam[,t], lambda.agg, K.agg, sdd.rate, dem.st)
        
        # 4. Short distance dispersal
        cat("Dispersing locally...")
        N.emig <- sdd_simple(N.lam[,t], N.new, sdd.pr, sdd.rate, K.agg, sdd.st)
        
        # 5. Long distance dispersal
        cat("Dispersing regionally...")
        N.emig <- ldd_disperse(ncell, N.emig, n.ldd, simple=TRUE)
        
        # 6. Update population sizes
        cat("Updating abundances.\n")
        N[,t+1] <- N.emig$N
      }
    } else {
      # 1. Initialize populations
      N <- array(0, dim=c(ncell, tmax+1, n.class))  
      N[,1,] <- N.init
      N.ling <- rep(0, ncell)
      N.sb <- rep(0, ncell)
      
      for(t in 1:tmax) {
        # 2. Pre-multiply compositional parameters
        lc.mx <- as.matrix(lc.df[,4:9])
        K.agg <- lc.mx %*% K
        pr.s.agg <- c(lc.mx %*% pr.s)
        rel.dens <- t(apply(lc.mx, 1, function(x) K*x/c(x%*%K)))
        fec.agg <- lc.mx %*% fec
        pr.f.agg <- lc.mx %*% pr.f
        pr.eat.agg <- lc.mx %*% pr.eat
        pr.est.agg <- lc.mx %*% pr.est
        pr.sb.agg <- lc.mx %*% pr.sb
        
        # 3. Local fruit production
        cat("Year", t, "- Fruiting...")
        N.f <- make_fruits(N[,t,], lc.mx, age.mat, fec.agg, pr.f.agg,
                           n.class, rel.dens, dem.st, mat.d)
        
        # 4. Short distance dispersal
        cat("Dispersing locally...")
        N.seed <- sdd_fs(N.f, pr.eat.agg, sdd.pr, sdd.rate, sdd.st)
        
        # 5. Long distance dispersal
        cat("Dispersing regionally...")
        N.seed <- ldd_disperse(ncell, N.seed, n.ldd, simple=FALSE)
        
        # 6. Seedling establishment
        cat("Establishing...")
        estab.out <- new_seedlings(ncell, N.seed, N.sb, pr.est.agg, pr.sb.agg,
                                   dem.st, bank)
        N.ling <- estab.out$N.rcrt
        N.sb <- estab.out$N.sb
        
        # 7. Update abundances
        cat("Updating abundances.\n")
        N[,t+1,n.class] <- pmin(round(N[,t,n.class] + 
                                        N[,t,n.class-1] * pr.s.agg),
                          round(K.agg))
        N[,t+1,2:(n.class-1)] <- round(N[,t,1:(n.class-2)] * pr.s.agg)
        N[,t+1,1] <- N.ling
      }
    }
    return(list(N=N, N.sb=N.sb))
  }


  
##---
## short distance dispersal probabilities
##---
  sdd_set_probs <- function(lc.df, g.p) {
    # Assign base dispersal probabilities from each cell
    # Each layer [1:i,1:j,,n] is the SDD neighborhood for cell n
    # trunc.diag: if TRUE, the sdd neighborhood is restricted to within sdd.max
    #   in all directions; if FALSE, then sdd.max refers to the maximum
    #   allowable horizontal & vertical distance from the source cell and the
    #   resulting neighborhood is square
    # k=1 contains pr(SDD | center,i,j)
    # k=2 contains the ID for each cell in the neighborhood
    # Returns array with dim(i:disp.rows, j:disp.cols, k:2, n:ncell)
    require(purrr); require(tidyverse); require(pbapply); require(fastmatch)
    
    # unpack parameters
    sdd.max <- g.p$sdd.max
    sdd.rate <- g.p$sdd.rate
    bird.hab <- g.p$bird.hab
    trunc.diag <- g.p$trunc.diag
    
    # initialize landscape & storage objects
    n.x <- max(lc.df[,1])
    n.y <- max(lc.df[,2])
    ncell <- n.x*n.y
    nbr <- 2 * sdd.max + 1
    sdd.i <- array(0, dim=c(nbr, nbr, 2, ncell))
    bird.hab.agg <- as.matrix(lc.df[,4:9]) %*% (bird.hab %>% divide_by(sum(.)))
    
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
    
    # pair cell IDs for each neighborhood; indexes match neighborhood matrix
    xx <- map(lc.df$x, ~seq(.-sdd.max, .+sdd.max))
    yy <- map(lc.df$y, ~seq(.-sdd.max, .+sdd.max))

    # create lists of inbounds xy neighborhood ranges
    n_ix <- map(xx, ~.[.>0 & .<=n.x])
    n_iy <- map(yy, ~.[.>0 & .<=n.y])
    
    # generate all xy combinations & neighborhood matrix indices
    cat("generating neighborhoods...\n")
    n_i <- map2(n_ix, n_iy, expand_v) 
    n_x <- map2(xx, n_ix, `%fin%`) %>% map(which) %>% map(range)
    n_y <- map2(yy, n_iy, `%fin%`) %>% map(which) %>% map(range)
    
    # match xy combinations with cell IDs
    cat("determining neighborhood cell IDs...\n")
    pboptions(type="none")
    if(g.p$n_cores > 1) {
      p.c <- makeCluster(g.p$n_cores)
      c_i <- pblapply(n_i, function(x) fastmatch::fmatch(x, lc.df$x_y), cl=p.c)
      stopCluster(p.c)
    } else {
      c_i <- pblapply(n_i, function(x) fastmatch::fmatch(x, lc.df$x_y))
    }
    
    cat("calculating probabilities...\n")
    for(n in 1:ncell) {
      # find cell ID for each cell in neighborhood
      sdd.i[n_x[[n]][1]:n_x[[n]][2], n_y[[n]][1]:n_y[[n]][2],2,n] <- c_i[[n]]
      # weight by bird habitat preference
      ib <- sdd.i[,,2,n] != 0  # inbounds neighbors
      sdd.i[,,1,n][ib] <- d.pr[ib] * bird.hab.agg[sdd.i[,,2,n][ib]]
      # set cell ID to 0 if pr(target) == 0 
      sdd.i[,,2,n][sdd.i[,,1,n]==0] <- 0
      # progress update
      if(n %% 5000 == 0) {
        cat("finished cell", n, "\n")
      }
    }
    cat("finished:", n, "cells\n")
    sdd.i[,,1,] <- apply(sdd.i[,,1,], 3, function(x) x/sum(x))
    return(sdd.i)
  }

  
  
  
##---
## local population growth: simple (a la Merow 2011)
##---
  grow_pops <- function(N.t, lambda.agg, K.agg, sdd.rate, dem.st=F) {
    # Calculate updated population sizes after local reproduction 
    # Growth rates are habitat specific
    # Returns sparse dataframe N.new with:
    #   col(id, N.pop, N.new=pop.delta, N.pop.upd=N.pop+N.new.nonemigrants)
    #   nrow = sum(N.new != 0)
    
    if(dem.st) {
      
    } else {
      N.id <- which(N.t>0)
      lam.id <- lambda.agg[N.id]
      N.new <- tibble(id = which(N.t>0)) %>%
        mutate(N.pop=N.t[id],
               N.new=(N.pop * (lam.id-1)),
               N.pop.upd=pmin(K.agg[N.id,],
                              N.pop + 
                                (lam.id>=1)*N.new*pexp(0.5, sdd.rate) +
                                (lam.id<1)*N.new) %>% round)
    }
    return(N.new)
  }
  
  
  
##---
## short distance dispersal: simple
##---
  sdd_simple <- function(N.t, N.new, sdd.pr, sdd.rate, K.agg, sdd.st=F) {
    # Calculate (N.arrivals | N.new, sdd.probs)
    # Accounts for distance from source cell & bird habitat preference
    # Returns dataframe with total population sizes.
    
    N.source <- N.new %>% filter(N.new > 0)
    N.emig <- tibble(id=1:length(N.t), N=N.t) %>%
      add_row(id=apply(N.source, 1, function(x) c(sdd.pr[,,2,x[1]])) %>% c, 
              N=apply(N.source, 1, function(x) c(x[3] * (1-pexp(.5,sdd.rate)) * 
                                      sdd.pr[,,1,x[1]])) %>% c) %>%
      filter(id != 0) %>% group_by(id) %>%
      summarise(N=sum(N))
    N.emig %<>% mutate(N=pmin(K.agg[N.emig$id,], N) %>% round)
    return(N.emig)
  }
  
  
  
##---
## local fruit production
##---
  make_fruits <- function(N.t, lc.mx, age.mat, fec.agg, pr.f.agg, 
                          n.class, rel.dens, dem.st=F, mat.d=F) {
    # Calculate (N.fruit | N, fec, age.mat) for each cell
    # fec, pr.f, & age.mat are habitat specific
    # Assumes no fruit production before age.mat
    # Individuals are assumed to be distributed among LC's relative to K:
    # K.tot = K[1]*lc.mx[1] + ... + K[l]*lc.mx[l]
    # N[l] = N.tot * K[l]*lc.mx[l]/K.tot
    # N.mature = sum(N[,age.mat[l]:8])*K[l]*lc.mx[l]/K.tot
    # N.fruit = N.mature*pr(Fruit)*fec
    # Returns sparse dataframe N.f with:
    #   col(id, N.rpr=num.reproducing, N.fruit=total.fruit)
    #   nrow = sum(N.frt != 0)
    
    
    # calculate N.mature in each LC in each cell
    if(mat.d) {  # does age at maturity differ by LC?
      names(age.mat) <- colnames(lc.mx)
      N.mature <- (map_df(age.mat, 
                          ~rowSums(matrix(N.t[,.:n.class]))) * rel.dens) %>% 
        rowSums %>% round
    } else {
      N.mature <- N.t[,n.class]
    }
    if(dem.st) {
      N.f <- tibble(id = which(N.mature>0)) %>%
        mutate(N.rpr = rbinom(n(), N.mature[id],
                              prob=pr.f.agg[id]),
               N.fruit = rpois(n(), 
                               lambda=N.rpr*fec.agg[id])) %>% 
        filter(N.fruit > 0)
    } else {
      N.f <- tibble(id = which(N.mature>0)) %>%
        mutate(N.rpr=(N.mature[id]) * pr.f.agg[id,],
               N.fruit=(N.rpr * fec.agg[id,]) %>% round) %>% 
        filter(N.fruit > 0)
    }
    return(N.f)
  }


  
##---
## short distance dispersal: fruits & seeds
##---
  sdd_fs <- function(N.f, pr.eat.agg, sdd.pr, sdd.rate, sdd.st=F) {
    # Calculate (N.seeds | N.fruit, sdd.probs, pr.eaten)
    # Accounts for distance from source cell, bird habitat preference,
    #   and the proportion of fruits eaten vs dropped
    # N.emig excludes seeds deposited within the source cell (1-pexp(0.5, rate))
    # Returns sparse matrix N.sdd with:
    #   cols(cell.ID, (N.seed = 2.3*(N.fruit - N.eaten + N.deposited)))
    #   nrow = sum(N.seed != 0)
    
    # calculate seeds deposited within source cell vs emigrants
    N.source <- N.f %>%
      mutate(N.produced=(2.3*N.fruit),
             N.emig=N.produced*(1-pexp(.5,sdd.rate))*pr.eat.agg[id,],
             N.drop=N.produced-N.emig)
    N.seed <- N.source %>% select(id, N.drop) %>% rename(N.dep=N.drop)
    
    if(sdd.st) {
      N.seed$N.dep <- round(N.seed$N.dep)
      SDD_sd <- unlist(apply(N.source, 1,
                             function(x) sample(sdd.pr[,,2,x[1]], x[5], 
                                                replace=TRUE,
                                                prob=sdd.pr[,,1,x[1]])))
      SDD_dep <- tabulate(SDD_sd)  # vector of counts for 1:max(SDD_sd)
      SDD_nonzero <- SDD_dep > 0  # cell id's with N_dep > 0
      N.seed <- add_row(N.seed, 
                        id=which(SDD_nonzero), 
                        N.dep=SDD_dep[SDD_nonzero])
    } else {
      # assign emigrants to target cells & sum within each cell
      N.seed %<>% 
        add_row(id=apply(N.source, 1, 
                         function(x) c(sdd.pr[,,2,x[1]])) %>% c, 
                N.dep=apply(N.source, 1, 
                            function(x) c(x[5] * sdd.pr[,,1,x[1]])) %>% c) %>%
        filter(N.dep > 0)
    }
    N.seed %<>%
      group_by(id) %>% 
      summarise(N=sum(N.dep) %>% round)
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
  new_seedlings <- function(ncell, N.seed, N.sb, pr.est.agg, pr.sb.agg,
                            dem.st=F, bank=F) {
    # Calculate (N.new | N.seed, pr.est)
    # Allows for incorporation of management effects & seedbank
    
    N.rcrt <- rep(0, ncell)
    if(dem.st) {
      N.rcrt[N.seed$id] <- rbinom(nrow(N.seed), N.seed$N, pr.est.agg[N.seed$id])
      if(bank) {
        N.sbEst <- rep(0, ncell)
        # N_est_sb
        N.sbEst[N.seed$id] <- rbinom(nrow(N.seed), N.sb[N.seed$id], 
                                     pr.est.agg[N.seed$id])
        # N_est_tot = N_est + N_est_sb
        N.rcrt[N.seed$id] <- N.rcrt[N.seed$id] + N.sbEst[N.seed$id]
        # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
        N.sb[N.seed$id] <- rbinom(nrow(N.seed),
                                  N.sb[N.seed$id] + N.seed$N - N.rcrt[N.seed$id],
                                  pr.sb.agg[N.seed$id])
      } else {
        N.sb <- rep(0, ncell)
      }
    } else {
      # N_est = N_seed * p(est)
      N.rcrt[N.seed$id] <- (N.seed$N * pr.est.agg[N.seed$id,])
      if(bank) {
        # N_est_tot = N_est + N_est_sb
        N.rcrt[N.seed$id] <- (N.rcrt[N.seed$id] + 
          N.sb[N.seed$id] * pr.est.agg[N.seed$id,]) %>% round
        # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
        N.sb[N.seed$id] <- ((N.sb[N.seed$id]*(1-pr.est.agg[N.seed$id,]) + 
                              N.seed$N - N.rcrt[N.seed$id]) * 
          pr.sb.agg[N.seed$id]) %>% round
      } else {
        N.sb <- rep(0, ncell)
      }
    }
    return(list(N.rcrt=N.rcrt, N.sb=N.sb))
  }







##########----------------------------------------------------------------------
##### MUNGING & PLOTTING FUNCTIONS
##########----------------------------------------------------------------------

##---
## make cell-block reference
##---
make_cb_i <- function(blockSize) {
  read_csv(paste0("data/roads_01_1a.csv")) %>% 
    mutate(CellRow=1:n_distinct(top) %>% rep(n_distinct(left)),
           CellCol=1:n_distinct(left) %>% rep(each=n_distinct(top))) %>%
    filter((CellRow <= max((CellRow %/% blockSize) * blockSize)) &
             (CellCol <= max((CellCol %/% blockSize) * blockSize))) %>%
    mutate(BlockRow=((CellRow-1)%/%blockSize)+1, 
           BlockCol=((CellCol-1)%/%blockSize)+1,
           BlockID=paste(str_pad(BlockCol, 7, "left", "0"), 
                         str_pad(BlockRow, 7, "left", "0"), sep=".") %>% 
             as.numeric %>% factor %>% as.numeric) %>%
    select(c(CellID, CellRow, CellCol, BlockID, BlockRow, BlockCol, left, top))
}
  
  
  
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



##---
## expand all pairwise combinations of two vectors into one character vector
##---
expand_v <- function(x,y) {
  paste(rep.int(x, length(y)), 
        rep.int(y, rep.int(length(x),length(y))),
        sep="_")
}





