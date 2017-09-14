########################################################
########################################################
########       PART 1: Load Data            ############
########################################################
########################################################
library(fields)
#The following file contains a matrix called 'nelandscape'. the first two columns give coordinates of each cell on the 82 x 79 grid overlayed on the landscape with (1,1) being the southwest corner. The third column gives the LULC classification (1= water, 2=developed, 3=agriculture, 4=deciduous, 5=coniferous).
print(load('nelandscape-appendix.RData')) #included data file from supplementary materials

########################################################
########################################################
########       PART 2: Functions            ############
########################################################
########################################################
	#FUNCTION FOR SIMULATIONS
cmstarling3dd=function(birddisp,landscape,lat0,long0, generations ,rate,plantprefs, carrying.cap,local.growth=TRUE, long.distance=TRUE,birds=TRUE,num.long.dist=1, recordtimes){
 	#birddisp:an array giving the probability of bird dispersal from each source cell to its respective bird dispersal neighborhood
 	#landscape: matrix with columns giving x coords, y coords and LULC for each cell
 	#lat0, long0: initial populations
 	#generations: # of time steps to interate model
 	#rate: 1/mean of the exponential dispersal kernel
 	#plantprefs: a vector of population growth rates where the first element corresponds to LULC class 1, etc.
 	#carrying.cap: carrying capacity in a single cell
 	#local.growth: logical - include local growth?
 	#long.distance: logical - include LDD?
 	#birds: logical - include local bird dispersal?
 	#num.long.dist: # LDD events/year
 	#recordtimes: steps where a snapshot is take of population distribution
		#########   GENERAL     ###################################################
	init=cbind(long0, lat0)
	lat=max(landscape[,2])
	long=max(landscape[,1]) 
		#record keeping variables
	p=cbind(landscape,rep(0,lat*long))    #state variable for each cell 
	timeseries=matrix(0,lat*long,length(recordtimes) )#to store snapshots of population
	ts.counter=1 
	ocean.sites=which(landscape[,3]==1) #
	good.sites=(1:length(landscape[,1]))[-ocean.sites]
	for(i in 1:nrow(init)){ #plants 1st seed(s)
		p[which(landscape[,1]==init[i,1] & landscape[,2]==init[i,2]),4]=carrying.cap/2 
	}

		########## DISPERSE SEEDS ############################     
	for (t in 1:generations){   
			#put a seed at Durham, NH in 1938 
		if(t==14){p[which(landscape[,1]==35 & landscape[,2]==28),4]=carrying.cap/2}
			#LOCAL GROWTH ______________________________________________
		populated.sites=which(p[,4]>0)
		new.individuals=p[populated.sites,4]*(plantprefs[p[populated.sites,3]]-1)
		p[which(p[,3]==1) ,4]=0 #kill at unsuitable (ocean) sites
		if(local.growth){
			p[populated.sites,4]=ceiling(pmin(carrying.cap,p[populated.sites,4]+ as.numeric(plantprefs[p[populated.sites,3]]>=1)*new.individuals*pexp(.5,rate)+as.numeric(plantprefs[p[populated.sites,3]]<1)*new.individuals))
	    }
    		#RANDOM LONG DISTANCE DISPERSAL _________________________________
    	if (long.distance){
		    for( i in 1:num.long.dist){
		    	start.site=sample(which(p[,4]>=1),1)
		    	target.site=sample(which(p[,3]>1),1) 
		 		while(length(target.site)<1){ #to ensure suitable sites found
		 			start.site=sample(which(p[,4]>=1),1)
		    		target.site=sample(which(p[,3]>1),1)
		    	}
				newsite=sample(target.site,1)
		    	p[newsite,4]=min(carrying.cap,p[newsite,4]+1)
    		}
    	}
    		#BIRD DISPERSAL ______________________________________________________
	    	#set up the # of offspring from each site
	    emigrate=rep(0,length(p[,4]))
	    emigrate[populated.sites]=new.individuals
	    emigrate=pmax(0,emigrate) #get rid of sink pops
    	if( birds ){
    		num.bird.seeds=emigrate*(1-pexp(.5,rate))	
    		for( k in which(round(num.bird.seeds,0)>0) ) {       
       			newsite=try(sample(subset(c(birddisp[,,2,k]),c(birddisp[,,2,k]>0)), round(num.bird.seeds[k],0),prob=subset(c(birddisp[,,1,k]),c(birddisp[,,2,k]>0)) ,replace=T),TRUE) #get the new sites
				p[newsite,4]=pmin(carrying.cap, p[newsite,4]+1) #place the seeds in the new sites
    		}
    	}
    		#_____________________________________________________________
				#for specified time intervals, store the population matrix
    	if ( sum(t==recordtimes)==1 ) {
    		timeseries[ ,ts.counter]=p[1:6478 ,4]
    		ts.counter=ts.counter+1
    	}
    } # time loop
	list(timeseries=timeseries) 
}  

###########################################################
	#FUNCTION FOR DISPERSAL PROBABILITIES
	# determines the probability of being dispersed to each site from site i,j
dispersal.probs=function(landscape,maxdist,rate){
		#this function is always called outside of the cmstarling program
	lat=max(landscape[,2])
	long=max(landscape[,1])
		#avoids calculation at unsuitable sites
	ocean.sites=which(landscape[,3]==1) 
	good.sites=(1:length(landscape[,1]))[-ocean.sites]
		#make a matrix of distance weights
	weights=matrix(0,2*maxdist+1,2*maxdist+1)
	center=c(maxdist+1,maxdist+1)
	for(i in seq(1,2*maxdist+1)){
		for(j in seq(1,2*maxdist+1)){
			weights[i,j]=dexp( (i-center[1])^2 +(j-center[2])^2-.5,rate )
		}
	}
	weights[center[1],center[2]]=0
	weights=weights/sum(weights)
	disp.matrix=array(0,c(nrow(weights),ncol(weights),2,lat*long))
	disp.matrix[,,1,good.sites]=weights
	for( k in good.sites){
		xx=(landscape[k,1]-maxdist):(landscape[k,1]+maxdist)
		yy=(landscape[k,2]-maxdist):(landscape[k,2]+maxdist)
	    for(i in xx[xx>0 & xx<=long]){
	    	for(j in yy[yy>0 & yy<=lat]){
	    		disp.matrix[which(xx==i),which(yy==j),2,k]=which(landscape[,1]==i & landscape[,2]==j)
	    	}
	    }
    }
	return(disp.matrix)
} 

########################################################
########################################################
########       PART 3: Sample Model Run     ############
########################################################
########################################################
landscape=nelandscape
birdtime=c(0,.39,.44,.06,.11)#these correspond to the LULC types (ocean, developed, agriculture, deciduous, coniferous)
	#MAKE DISPERSAL PROBABILITIES ---------------------------------------------
maxdist=3	#sets size of local bird dispersal neighborhood
rate=3.5	#1/mean 
	#this calculates the distance weighted probabilities of dispersal
dist.mat=dispersal.probs(landscape,maxdist,rate)
	#this weights distances by bird habitat use
a=0*dist.mat[,,1,]
for( i in which(!landscape[,3]==1 & !landscape[,3]==0)) {
	for( j in 1:(2*maxdist+1)){
		for(k in 1:(2*maxdist+1)){
			if(!dist.mat[j,k,2,i]==0){
				a[j,k,i]=dist.mat[j,k,1,i]*birdtime[landscape[dist.mat[j,k,2,i],3]]
			}
		}
	}
	a[,,i]=a[,,i]/sum(a[,,i])
}
birddisp=dist.mat
birddisp[,,1,]=a

	#RUN MULTIPLE SIMULATIONS AND PRODUCE AVERAGE SURFACE----------------------
reps=10 # number of replicate model runs. 100 runs were used in the manuscript
lat=max(landscape[,2]);	long=max(landscape[,1]) 
d.reps=array(NA,c(lat*long,5,reps)) # to store output
	#PARAMETERS
plantprefs=c(0, 2.1, 1.5, 1.4, .5) #these correspond to the LULC types (water, developed, agriculture, deciduous, coniferous)
recordtimes=c(21,41,61,81,90) # specify the time steps at which to record the results 
generations=90
num.long.dist=1 # # of seeds per year for LDD
local.growth=TRUE
long.distance=TRUE
birds=TRUE
long0=c(9,39);	lat0=c(5,9)  #initial conditions
carrying.cap=200  
for(i in 1:reps){
	d=cmstarling3dd(birddisp,landscape,lat0,long0,generations,rate,plantprefs, carrying.cap,local.growth,long.distance,birds,num.long.dist,recordtimes=recordtimes)
	d.reps[,,i]=d$timeseries
	print(i)
}
d.reps=ifelse(d.reps>0,1,0) #turn abundance predicitons in to presence/absence
avg.occupancy=matrix(NA,lat*long,5)
for(i in 1:5){ avg.occupancy[,i]=apply(d.reps[,i,],1,mean) }
	
########################################################
########################################################
########       PART 4: Plot Results    #################
########################################################
########################################################
years=c(1940,1960,1980,2000,2009) #labels for plotting
threshold=.50
avg.occupancy[landscape[,3]==1,]=2.2  # a placeholder for water so the plot looks pretty
zc=array(0,c(82,79,5)) #create a grid for use in 'image' function
for(i in 1:nrow(landscape)){ zc[landscape[i,1],landscape[i,2],]=avg.occupancy[i,] }
par(oma=c( 1,1,3,5),mfrow=c(2,3),mar=c(0,0,0,0)) 
for(i in 1:5){
	image(1:82,1:79,zc[,,i],col=colorRampPalette(c('grey90','lightsteelblue2', 'steelblue4','darkslateblue','violetred3', 'white'),bias=3)(100), col.axis='white',xaxt='n',yaxt='n',xlim=c(-1,84),ylim=c(-1,81),bty='n')
	text(15,70,years[i],cex=1.7)
	points(c(long0,35),c(lat0,28),cex=1.2,lwd=2,col='grey90') #introduction points
}
par(oma=c(2,0,8,18))
image.plot(legend.only=T,zlim=c(0,1),col=colorRampPalette(c('grey90','lightsteelblue2','steelblue4','darkslateblue','violetred3'),bias=3)(100),legend.width=8)



