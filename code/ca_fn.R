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
  run_sim <- function(ngrid, ncell, g.p, lc.df, sdd.pr, 
                      N.init, control.p=NULL) {
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
    n.lc <- g.p$n.lc
    simple <- g.p$simple
    dem.st <- g.p$dem.st
    sdd.st <- g.p$sdd.st
    bank <- g.p$bank
    K <- g.p$K  # carrying capacity
    pr.s <- g.p$pr.s  # pre-adult survival
    pr.f <- g.p$pr.f  # pr(fruit)
    fec <- g.p$fec  # mean(fruit per adult)
    age.f <- g.p$age.f  # mean age at first fruiting
    pr.est <- g.p$pr.est  # pr(seedling est)
    pr.sb <- g.p$pr.sb  # pr(ann.surv seed bank)
    lambda <- g.p$lambda  # pop growth rate
    sdd.rate <- g.p$sdd.rate  # 1/mn for dispersal kernel
    pr.eat <- g.p$pr.eat  # pr(birds eat frt)
    pr.s.bird <- g.p$pr.s.bird  # pr(viable | digestion)
    n.ldd <- g.p$n.ldd   # num long distance dispersal events per year
    y.ad <- max(g.p$age.f)
    age.f.d <- length(age.f) > 1
    id_i <- lc.df %>% select(id, id_inbd)
    
    # If buckthorn is being actively managed...
    pr.est.trt <- NULL
    if(!is.null(control.p)) {
      nTrt_grd <- control.p$nTrt_grd * ncell
      nTrt_man <- control.p$nTrt_man * ncell
      grd.trt <- control.p$grd.trt
      man.trt <- control.p$man.trt
      t.trt <- control.p$t.trt
      est.trt <- tibble(CellID=numeric(), Trt=character())
      N.trt <- tibble(CellID=numeric(), Trt=character())
    }
    
    if(simple) {
      # 1. Initialize populations
      N.lam <- matrix(0, ncell, tmax+1)  
      N.lam[,1] <- rowSums(N.init)
      
      for(t in 1:tmax){
        # 2. Pre-multiply compositional parameters
        K.ag <- as.matrix(lc.df[,4:9]) %*% K
        lambda.ag <- as.matrix(lc.df[,4:9]) %*% lambda
        
        # 3. Local growth
        cat("Year", t, "- Growing...")
        N.new <- grow_pops(N.lam[,t], lambda.ag, K.ag, sdd.rate, dem.st)
        
        # 4. Short distance dispersal
        cat("Dispersing locally...")
        N.emig <- sdd_simple(N.lam[,t], N.new, sdd.pr, sdd.rate, K.ag, sdd.st)
        
        # 5. Long distance dispersal
        cat("Dispersing regionally...")
        N.emig <- ldd_disperse(ncell, N.emig, n.ldd, simple=TRUE)
        
        # 6. Update population sizes
        cat("Updating abundances.\n")
        N[,t+1] <- N.emig$N
      }
    } else {
      # 1. Initialize populations
      N.sb <- matrix(0, nrow=ngrid, ncol=tmax+1)
      if(age.f.d) {
        N <- array(0, dim=c(ngrid, tmax+1, n.lc, y.ad))  
        N[,1,,] <- N.init
        for(l in 1:n.lc) {
          if(age.f[l] < y.ad) {
            N[,,l,age.f[l]:(y.ad-1)] <- NA
          }
        }
      } else {
        N <- array(0, dim=c(ngrid, tmax+1, y.ad))  
        N[,1,] <- N.init
      }
      
      for(t in 1:tmax) {
        if(age.f.d) {
          N.t <- N[,t,,]
        } else {
          N.t <- N[,t,]
        }
        
        # 2. Implement management
        if(!is.null(control.p) && t >= t.trt) {
          # run function to implement management controls
          # controls may be on whole cells or specific LCs of specific cells
          # controls may affect:
          #  - LC % (e.g., extreme timber harvest)
          #  - pr.est (e.g., litter or cover crops)
          #  - N (e.g., cutting and spraying)
          # Inputs need to specify which LCs & cells will be affected, which
          # method or methods will be used, and potentially the intensity of
          # the method manifested in the decrease in N or pr.est.
          # Ultimately, a new set of inputs will be provided in each year 
          # depending on decisions made in the economic model.
          #
          # For now, the inputs (all sparse) will be:
          #  - lc.trt: df(col=c(CellID, Hwd, WP, Evg, Mxd)) % change to OpI
          #  - est.trt: df(col=c(CellID, Trt)) where Trt=Litter|Cover|Compact
          #  - N.trt: df(col=c(CellID, Trt)) where Trt=Mech|Chem|MechChem
          #  - grd.trt: vector(Litter=p.est, Cover=p.est, Compact=p.est)
          #  - man.trt: vector(Mech=%kill, Chem=%kill, MechChem=%kill)
          #
          # For future complexity, manual.trt could also push a proportion of
          # the adults back to a previous age so they don't fruit for a number
          # of years after the treatment rather than being killed explicitly
          
          # A. Adjust LC %
          lc.trt <- NULL # to do
          
          # B. Adjust p.est
          if(nTrt_grd > 0) {
            est.trt %<>% add_row(CellID=sample(1:ncell, nTrt_grd),
                                 Trt=sample(c("Lit", "Cov", "Com"), 
                                            nTrt_grd, replace=TRUE))
            pr.est.trt <- ground_trt(est.trt, grd.trt)
          }
          
          # C. Adjust N
          if(nTrt_man > 0) {
            N.trt %<>% add_row(CellID=sample(1:ncell, nTrt_man),
                               Trt=sample(c("M", "C", "MC"), 
                                          nTrt_man, replace=TRUE))
            if(age.f.d) {
              N[,t,,] <- manual_trt(N.t, y.ad, N.trt, man.trt)
            } else {
              N[,t,] <- manual_trt(N.t, y.ad, N.trt, man.trt)
            }
          }
        }
        
        # 3. Pre-multiply compositional parameters
        pm <- premultiply(lc.df, K, pr.s, fec, pr.f, pr.eat, pr.est, pr.est.trt)
        
        # 4. Local fruit production
        cat("Year", t, "- Fruiting...")
        N.f <- make_fruits(N.t, pm$lc.mx, age.f.d, pm$fec.ag, pm$pr.f.ag,
                           y.ad, dem.st)
        
        # 5. Short distance dispersal
        cat("Dispersing locally...")
        N.seed <- sdd_fs(id_i, N.f, pm$pr.eat.ag, pr.s.bird, 
                         sdd.pr, sdd.rate, sdd.st)
        
        # 6. Long distance dispersal
        cat("Dispersing regionally...")
        N.seed <- ldd_disperse(ncell, id_i, N.seed, n.ldd, simple=FALSE)
        
        # 7. Seedling establishment
        cat("Establishing...")
        estab.out <- new_seedlings(ngrid, N.seed, N.sb[,t], pm$pr.est.ag, 
                                   pr.sb, dem.st, bank)
        N.sb[,t+1] <- estab.out$N.sb
        
        # 8. Update abundances
        cat("Updating abundances.\n")
        if(age.f.d) {
          for(l in 1:n.lc) {
            N[,t+1,l,y.ad] <- pmin(round(N[,t,l,y.ad] + 
                                           N[,t,l,age.f[l]-1] * pr.s[l]),
                                   pm$K.lc[,l])
            N[,t+1,l,2:(age.f[l]-1)] <- round(N[,t,l,1:(age.f[l]-2)] * pr.s[l])
            N[,t+1,l,1] <- round(estab.out$N.rcrt * pm$rel.dens[,l])
          }
        } else {
          N[,t+1,y.ad] <- pmin(round(N[,t,y.ad] + N[,t,y.ad-1] * pm$pr.s.ag),
                               pm$K.ag)
          N[,t+1,2:(y.ad-1)] <- round(N[,t,1:(y.ad-2)] * pm$pr.s.ag)
          N[,t+1,1] <- estab.out$N.rcrt
        }
      }
    }
    if(age.f.d) {
      N <- apply(N, c(1,2,4), sum, na.rm=TRUE)
    }
    return(list(N=N, N.sb=N.sb))
  }


  
##---
## short distance dispersal probabilities
##---
  sdd_set_probs <- function(ncell, lc.df, g.p) {
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
    
    # initialize landscape & storage objects
    n.x <- range(lc.df[,1])
    n.y <- range(lc.df[,2])
    nbr <- 2 * sdd.max + 1
    sdd.i <- array(0, dim=c(nbr, nbr, 2, ncell))
    bird.hab.ag <- as.matrix(lc.df[,4:9]) %*% (bird.hab %>% divide_by(sum(.)))
    
    # generate default dispersal probability matrix
    d.pr <- matrix(0, nbr, nbr)
    ctr <- sdd.max + 1  # center index (i=j) for square mx
    for(i in 1:nbr) {
      for(j in i:nbr) {
        d.pr[i,j] <- dexp((i-ctr)^2 + (j-ctr)^2 - 0.5, sdd.rate)
        if( sqrt((i-ctr)^2 + (j-ctr)^2) > sdd.max ) {
          d.pr[i,j] <- 0
        }
        d.pr[j,i] <- d.pr[i,j]
      }
    }
    d.pr <- d.pr/sum(d.pr)
    
    # pair cell IDs for each neighborhood; indexes match neighborhood matrix
    xx <- map(lc.df$x[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))
    yy <- map(lc.df$y[lc.df$inbd], ~seq(.-sdd.max, .+sdd.max))

    # create lists of on-map xy neighborhood ranges
    n_ix <- map(xx, ~.[.>=n.x[1] & .<=n.x[2]])
    n_iy <- map(yy, ~.[.>=n.y[1] & .<=n.y[2]])
    
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
      sdd.i[n_y[[n]][1]:n_y[[n]][2],
            n_x[[n]][1]:n_x[[n]][2],2,n] <- matrix(c_i[[n]], 
                                                   ncol=diff(n_x[[n]])+1,
                                                   byrow=TRUE)
      # weight by bird habitat preference
      ib <- sdd.i[,,2,n] != 0  # inbounds neighbors
      sdd.i[,,1,n][ib] <- d.pr[ib] * bird.hab.ag[sdd.i[,,2,n][ib]]
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
  grow_pops <- function(N.t, lambda.ag, K.ag, sdd.rate, dem.st=F) {
    # Calculate updated population sizes after local reproduction 
    # Growth rates are habitat specific
    # Returns sparse dataframe N.new with:
    #   col(id, N.pop, N.new=pop.delta, N.pop.upd=N.pop+N.new.nonemigrants)
    #   nrow = sum(N.new != 0)
    
    if(dem.st) {
      
    } else {
      N.id <- which(N.t>0)
      lam.id <- lambda.ag[N.id]
      N.new <- tibble(id = which(N.t>0)) %>%
        mutate(N.pop=N.t[id],
               N.new=(N.pop * (lam.id-1)),
               N.pop.upd=pmin(K.ag[N.id,],
                              N.pop + 
                                (lam.id>=1)*N.new*pexp(0.5, sdd.rate) +
                                (lam.id<1)*N.new) %>% round)
    }
    return(N.new)
  }
  
  
  
##---
## short distance dispersal: simple
##---
  sdd_simple <- function(N.t, N.new, sdd.pr, sdd.rate, K.ag, sdd.st=F) {
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
    N.emig %<>% mutate(N=pmin(K.ag[N.emig$id,], N) %>% round)
    return(N.emig)
  }
  
  
  
##---
## local fruit production
##---
  make_fruits <- function(N.t, lc.mx, age.f.d, fec.ag, pr.f.ag, 
                          y.ad, dem.st=F) {
    # Calculate (N.fruit | N, fec, age.f) for each cell
    # fec, pr.f, & age.f are habitat specific
    # Assumes no fruit production before age.f
    # Individuals are assumed to be distributed among LC's relative to K:
    # K.tot = K[1]*lc.mx[1] + ... + K[l]*lc.mx[l]
    # N[l] = N.tot * K[l]*lc.mx[l]/K.tot
    # N.mature = sum(N[,age.f[l]:8])*K[l]*lc.mx[l]/K.tot
    # N.fruit = N.mature*pr(Fruit)*fec
    # Returns sparse dataframe N.f with:
    #   col(id, N.rpr=num.reproducing, N.fruit=total.fruit)
    #   nrow = sum(N.frt != 0)
    
    
    # calculate N.mature in each LC in each cell
    if(age.f.d) {
      N.mature <- rowSums(N.t[,,y.ad])
    } else {
      N.mature <- N.t[,y.ad]
    }
    if(dem.st) {
      N.f <- tibble(id = which(N.mature>0)) %>%
        mutate(N.rpr = rbinom(n(), N.mature[id],
                              prob=pr.f.ag[id]),
               N.fruit = rpois(n(), 
                               lambda=N.rpr*fec.ag[id])) %>% 
        filter(N.fruit > 0)
    } else {
      N.f <- tibble(id = which(N.mature>0)) %>%
        mutate(N.rpr=(N.mature[id]) * pr.f.ag[id,],
               N.fruit=(N.rpr * fec.ag[id,]) %>% round) %>% 
        filter(N.fruit > 0)
    }
    return(N.f)
  }


  
##---
## short distance dispersal: fruits & seeds
##---
  sdd_fs <- function(id_i, N.f, pr.eat.ag, pr.s.bird, 
                     sdd.pr, sdd.rate, sdd.st=F) {
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
             N.emig=N.produced*(1-pexp(.5,sdd.rate))*pr.eat.ag[id,],
             N.drop=N.produced-N.emig) %>%
      mutate(N.emig=N.emig*pr.s.bird,
             id_inbd=id_i$id_inbd[id])
    N.seed <- N.source %>% select(id, N.drop) %>% rename(N.dep=N.drop)
    
    if(sdd.st) {
      N.seed$N.dep <- round(N.seed$N.dep)
      SDD_sd <- unlist(apply(N.source, 1,
                             function(x) sample(sdd.pr[,,2,x[7]], x[5], 
                                                replace=TRUE,
                                                prob=sdd.pr[,,1,x[7]])))
      SDD_dep <- tabulate(SDD_sd)  # vector of counts for 1:max(SDD_sd)
      SDD_nonzero <- SDD_dep > 0  # cell id's with N_dep > 0
      N.seed <- add_row(N.seed, 
                        id=which(SDD_nonzero), 
                        N.dep=SDD_dep[SDD_nonzero])
    } else {
      # assign emigrants to target cells & sum within each cell
      N.seed %<>% 
        add_row(id=apply(N.source, 1, 
                         function(x) c(sdd.pr[,,2,x[7]])) %>% c, 
                N.dep=apply(N.source, 1, 
                            function(x) c(x[5] * sdd.pr[,,1,x[7]])) %>% c) %>%
        filter(N.dep > 0)
    }
    N.seed %<>%
      group_by(id) %>% 
      summarise(N=sum(N.dep) %>% round) %>%
      filter(!is.na(id_i$id_inbd[id]) & N > 0)
    return(N.seed)
  }


  
##---
## long distance dispersal
##---
  ldd_disperse <- function(ncell, id_i, N.df, n.ldd, simple=TRUE) {
    # Assign n.ldd random LDD events 
    # A single seed is added to n.ldd target cells
    
    ldd.id <- id_i$id[which(id_i$id_inbd == 
                              sample(1:ncell, n.ldd, replace=TRUE))]
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
  new_seedlings <- function(ngrid, N.seed, N.sb, pr.est.ag, pr.sb,
                            dem.st=F, bank=F) {
    # Calculate (N.new | N.seed, pr.est)
    # Allows for incorporation of management effects & seedbank
    
    N.rcrt <- rep(0, ngrid)
    if(dem.st) {
      N.rcrt[N.seed$id] <- rbinom(nrow(N.seed), N.seed$N, pr.est.ag[N.seed$id])
      if(bank) {
        N.sbEst <- rep(0, ngrid)
        # N_est_sb
        N.sbEst[N.seed$id] <- rbinom(nrow(N.seed), N.sb[N.seed$id], 
                                     pr.est.ag[N.seed$id])
        # N_est_tot = N_est + N_est_sb
        N.rcrt[N.seed$id] <- N.rcrt[N.seed$id] + N.sbEst[N.seed$id]
        # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
        N.sb[N.seed$id] <- rbinom(nrow(N.seed),
                                  N.sb[N.seed$id] + N.seed$N - N.rcrt[N.seed$id],
                                  pr.sb)
      } else {
        N.sb <- rep(0, ngrid)
      }
    } else {
      # N_est = N_seed * p(est)
      N.rcrt[N.seed$id] <- N.seed$N * pr.est.ag[N.seed$id,]
      if(bank) {
        # N_est_tot = N_est + N_est_sb
        N.rcrt[N.seed$id] <- (N.rcrt[N.seed$id] + 
          N.sb[N.seed$id] * pr.est.ag[N.seed$id,]) %>% round
        # N_to_sb = (N_sb_notEst + N_addedToSB) * p(SB)
        N.sb[N.seed$id] <- ((N.sb[N.seed$id]*(1-pr.est.ag[N.seed$id,]) + 
                              N.seed$N - N.rcrt[N.seed$id]) * pr.sb) %>% round
      } else {
        N.sb <- rep(0, ngrid)
      }
    }
    return(list(N.rcrt=N.rcrt, N.sb=N.sb))
  }


  
##---
## initialize populations randomly
##---
pop_init <- function(ngrid, g.p, lc.df) {
  # Initialize populations randomly
  # Adults: 50% K; Subadults: 10% K
  # Accounts for age.f constant vs varying by LC
  # Returns N.init: matrix (row=cell, col=age)
  # Or array (row=cell, col=LC, layer=age)
  
  p.0 <- sample(lc.df$id[lc.df$inbd], g.p$N.p.t0)
  y.ad <- max(g.p$age.f)  # adult age bin
  if(length(g.p$age.f) == 1) {
    N.init <- matrix(0, ngrid, y.ad)  # column for each age class
    N.init[p.0,y.ad] <- round(as.matrix(lc.df[lc.df$id %in% p.0,4:9]) %*% 
                                (g.p$K/2))
    N.init[p.0,-y.ad] <- round(N.init[p.0,y.ad]/5)
  } else {
    N.init <- array(0, dim=c(ncell, g.p$n.lc, y.ad))
    N.init[p.0,,y.ad] <- round(t(t(as.matrix(lc.df[lc.df$id %in% p.0,4:9])) * 
                                   g.p$K/2))
    N.init[p.0,,-y.ad] <- round(N.init[p.0,,y.ad]/5)
  }
  return(N.init)
}



##---
## implement ground cover treatments
##---
ground_trt <- function(est.trt, grd.trt) {
  # given cells, treatments, and treatment effects: adjust pr.est
  pr.est.trt <- tibble(id=est.trt$CellID,
                       pr.est=grd.trt[match(est.trt$Trt, names(grd.trt))])
  return(pr.est.trt)
}



##---
## implement cutting & spraying treatments
##---
manual_trt <- function(N.t, y.ad, N.trt, man.trt) {
  # given cells, treatments, and treatment effects: adjust N.adults
  
  # spraying probably has a %kill for all N
  # some % die, some % get bumped back to juvenile
  # effects on germination rates?
  
  prop.kill <- tibble(id=N.trt$CellID,
                      prop=1-man.trt[match(N.trt$Trt, names(man.trt))])
  if(length(dim(N.t)) == 2) {
    N.t[prop.kill$id,] <- round(N.t[prop.kill$id,] * prop.kill$prop)
  } else {
    N.t[prop.kill$id,,] <- round(N.t[prop.kill$id,,] * prop.kill$prop)
  }
  return(N.t)
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



##---
## premultiply compositional data
##---
premultiply <- function(lc.df, K, pr.s, fec, pr.f, pr.eat, 
                        pr.est, pr.est.trt=NULL) {
  # reformats and calculates cell-means based on land cover composition
  # for relevant parameters. Specifically:
  #  lc.mx: matrix(col=LC, row=cell) with LC proportions
  #  K.ag: vector(cell) with total K
  #  K.lc: matrix(col=LC, row=cell) with K per LC
  #  pr.s.ag: vector(cell) with pr(surv)
  #  rel.dens: matrix(col=LC, row=cell) with relative density among LC
  #  fec.ag: vector(cell) with mn(fruit per adult)
  #  pr.f.ag: vector(cell) with pr(fruit)
  #  pr.eat.ag: vector(cell) with pr(eaten by bird)
  #  pr.est.ag: vector(cell) with pr(establish)
  # If !is.na(pr.s.trt), then the associated pr.s.ag values are substituted in
  # the cells that received a relevant management treatments
  
  lc.mx <- as.matrix(lc.df[,4:9])
  K.ag <- round(lc.mx %*% K)
  K.lc <- round(t(t(lc.mx) * K))
  rel.dens <- t(apply(lc.mx, 1, function(x) K*x/c(x%*%K)))
  pr.s.ag <- c(lc.mx %*% pr.s)
  fec.ag <- lc.mx %*% fec
  pr.f.ag <- lc.mx %*% pr.f
  pr.eat.ag <- lc.mx %*% pr.eat
  pr.est.ag <- lc.mx %*% pr.est
  
  if(!is.null(pr.est.trt)) {
    pr.est.ag[pr.est.trt$id,] <- pr.est.trt$pr.est
  }
  
  return(list(lc.mx=lc.mx, K.ag=K.ag, K.lc=K.lc, rel.dens=rel.dens,
              pr.s.ag=pr.s.ag, fec.ag=fec.ag, pr.f.ag=pr.f.ag,
              pr.eat.ag=pr.eat.ag, pr.est.ag=pr.est.ag))
}





