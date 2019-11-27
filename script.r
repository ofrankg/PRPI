# Oswaldo Franco
# CIC-B19
# PRPI
#================Install conda environment===============
#conda install -c conda-forge/label/cf201901 r-msm
#conda install -c conda-forge/label/gcc7 r-msm
#conda install -c conda-forge  r-msm
#conda install -c conda-forge r-markovchain
#conda install -c conda-forge/label/cf201901 r-markovchain
#conda install -c conda-forge/label/gcc7 r-markovchain
#conda install -c conda-forge r-etm
#conda install -c conda-forge/label/cf201901 r-etm
#conda install -c conda-forge/label/gcc7 r-etm
#=================execution==============================
#in R command windows source(file="path/script.r")

library("markovchain")
weatherStates <- c("sunny", "cloudy", "rain")
byRow <- TRUE
weatherMatrix <- matrix(data = c(0.70, 0.2, 0.1,
				 0.3, 0.4, 0.3,
				 0.2, 0.45, 0.35), byrow = byRow, nrow = 3,
				dimnames = list(weatherStates, weatherStates))
#Forma larga de crear la cadena de markov para clima del día
mcWeatherL <- new ("markovchain", states = weatherStates, byrow = byRow, transitionMatrix = weatherMatrix, name = "Weather")

#Forma corta de crear la cadena de markov para clima del día
mcWeatherS <- new("markovchain", states = c("sunny", "cloudy", "rain"),
		transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
					0.3, 0.4, 0.3,
					0.2, 0.45, 0.35), byrow = byRow, nrow = 3),name = "Weather")

#cadena de Markov por defecto
defaultMc <- new("markovchain")

#Mostrar resultados en pantalla
print(mcWeatherL)
print(mcWeatherS)
print(defaultMc)

#============Handling Markovchain objects===============

#row vector results

initialState <- c(0, 1, 0)
after2Days <- initialState * (mcWeatherL * mcWeatherL)
after7Days <- initialState * (mcWeatherL ^ 7)
print(after2Days)
print(round(after7Days,3))

#column vector result
initialState <- c(0, 1, 0)
after2Days <- (t(mcWeatherL) * t(mcWeatherL)) * initialState
after7Days <- (t(mcWeatherL) ^ 7) * initialState
print(after2Days)
print(round(after7Days,3))

#=======================================================
fvals<-function(mchain,initialstate,n) {
	out<-data.frame()
	names(initialstate)<-names(mchain)
	for (i in 0:n)
	{
		iteration<-initialstate*mchain^(i)
		out<-rbind(out,iteration)
	}
	out<-cbind(out, i=seq(0,n))
	out<-out[,c(4,1:3)]
	return(out)
}
print(fvals(mchain=mcWeatherL,initialstate=c(90,5,5),n=4))

#=======================================================
print(states(mcWeatherL))
print(names(mcWeatherL))
print(dim(mcWeatherL))
print(name(mcWeatherL))
name(mcWeatherL) <- "New Name"
print(name(mcWeatherL))

print(markovchain:::sort(mcWeatherL))
print(transitionProbability(mcWeatherL, "cloudy", "rain"))
print(mcWeatherL[2,3])
print(mcWeatherL)

show(mcWeatherL)
#plot(mcWeatherL)
#plot(mcWeatherL, package="diagram", box.size = 0.1)

mcDf <- as(mcWeatherL, "data.frame")
mcNew <- as(mcDf, "markovchain")
print(mcDf)

#mcIgraph <- as(mcWeatherL,"igraph")
#plot(mcIgraph)

require(msm)
Q <- rbind( c(0,0.25,0,0.25),
	   c(0.166,0,0.166,0.166),
	   c(0,0.25,0,0.25),
	   c(0,0,0,0) )
cavmsm <- msm(state~years, subject = PTNUM, data = cav, qmatrix = Q, death = 4)
msmMc <- as(cavmsm, "markovchain")
print(msmMc)

#=====================etm library=======================

library(etm)
data(sir.cont)
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time),]
for (i in 2:nrow(sir.cont)){
	if (sir.cont$id[i]==sir.cont$id[i-1]){
		if (sir.cont$time[i]==sir.cont$time[i-1]){
			sir.cont$time[i-1] <- sir.cont$time[i-1]-0.5
		}
	}
}
tra <- matrix(ncol=3, nrow=3, FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE
tr.prob <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)
print(tr.prob)
etm2mc<-as(tr.prob, "markovchain")
print(etm2mc)

myMatr<-matrix(c(.1,.8,.1,.2,.6,.2,.3,.4,.3), byrow=TRUE, ncol=3)
myMc<-as(myMatr, "markovchain")
print(myMc)

stateNames = c("H", "I", "D")
Q0 <- new("markovchain", states = stateNames, 
	  transitionMatrix =matrix(c(0.7, 0.2, 0.1,0.1, 0.6, 0.3,0, 0, 1), 
				   byrow = TRUE, nrow = 3), name = "state t0")
Q1 <- new("markovchain", states = stateNames, 
	  transitionMatrix = matrix(c(0.5, 0.3, 0.2,0, 0.4, 0.6,0, 0, 1), 
				    byrow = TRUE, nrow = 3), name = "state t1")
Q2 <- new("markovchain", states = stateNames, 
	  transitionMatrix = matrix(c(0.3, 0.2, 0.5,0, 0.2, 0.8,0, 0, 1), 
				    byrow = TRUE,nrow = 3), name = "state t2")
Q3 <- new("markovchain", states = stateNames, 
	  transitionMatrix = matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1), 
				    byrow = TRUE, nrow = 3), name = "state t3")
mcCCRC <- new("markovchainList",markovchains = list(Q0,Q1,Q2,Q3),
name = "Continuous Care Health Community")
print(mcCCRC)
print(mcCCRC[[1]])
print(dim(mcCCRC))
#===========Probability with markovchain objects=======
print(conditionalDistribution(mcWeatherL, "sunny"))
print(steadyStates(mcWeatherL))
gamblerRuinMarkovChain <- function(moneyMax, prob = 0.5) {
	require(matlab)
	matr <- zeros(moneyMax + 1)
	states <- as.character(seq(from = 0, to = moneyMax, by = 1))
	rownames(matr) = states; colnames(matr) = states
	matr[1,1] = 1; matr[moneyMax + 1,moneyMax + 1] = 1
	for(i in 2:moneyMax){
		matr[i,i-1] = 1 - prob; matr[i, i + 1] = prob
	}
	out <- new("markovchain",
		   transitionMatrix = matr,
		   name = paste("Gambler ruin", moneyMax, "dim", sep = " ") )
	return(out)
}
mcGR4 <- gamblerRuinMarkovChain(moneyMax = 4, prob = 0.5)
print( steadyStates(mcGR4) )
print( absorbingStates(mcGR4) )
print( absorbingStates(mcWeatherL) )
.commclassesKernel <- function(P){
	m <- ncol(P)
		stateNames <- rownames(P)
		T <- zeros(m)
		i <- 1
		while (i <= m) {
			a <- i
			b <- zeros(1,m)
			b[1,i] <- 1
			old <- 1
			new <- 0
			while (old != new) {
				old <- sum(find(b > 0))
				n <- size(a)[2]
				matr <- matrix(as.numeric(P[a,]), ncol = m,nrow = n)
				c <- colSums(matr)
				d <- find(c)
				n <- size(d)[2]
				b[1,d] <- ones(1,n)
				new <- sum(find(b>0))
				a <- d
			}
			T[i,] <- b
			i <- i+1 }
		F <- t(T)
		C <- (T > 0)&(F > 0)
		v <- (apply(t(C) == t(T), 2, sum) == m)
		colnames(C) <- stateNames
		rownames(C) <- stateNames
		names(v) <- stateNames
		out <- list(C = C, v = v)
		return(out)
}

P <- matlab::zeros(10)
P[1, c(1, 3)] <- 1/2;
P[2, 2] <- 1/3; P[2,7] <- 2/3;
P[3, 1] <- 1;
P[4, 5] <- 1;
P[5, c(4, 5, 9)] <- 1/3;
P[6, 6] <- 1;
P[7, 7] <- 1/4; P[7,9] <- 3/4;
P[8, c(3, 4, 8, 10)] <- 1/4;
P[9, 2] <- 1;
P[10, c(2, 5, 10)] <- 1/3;
rownames(P) <- letters[1:10]
colnames(P) <- letters[1:10]
probMc <- new("markovchain", transitionMatrix = P,
	      name = "Probability MC")
.commclassesKernel(P)
print(.commclassesKernel(P))
summary(probMc)
print( transientStates(probMc) )
probMcCanonic <- canonicForm(probMc) 
print( probMc )
print( probMcCanonic )
is.accessible(object = probMc, from = "a", to = "c")
is.accessible(object = probMc, from = "g", to = "c")
E <- matrix(0, nrow = 4, ncol = 4)
E[1, 2] <- 1
E[2, 1] <- 1/3; E[2, 3] <- 2/3
E[3,2] <- 1/4; E[3, 4] <- 3/4
E[4, 3] <- 1
mcE <- new("markovchain", states = c("a", "b", "c", "d"),
	   transitionMatrix = E,
	   name = "E")
print( is.irreducible(mcE) )
print( period(mcE) )
require(matlab)
mathematicaMatr <- zeros(5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
		     name = "Mathematica MC", states = statesNames)
summary( mathematicaMc )
#plot( mathematicaMc )
.firstpassageKernel <- function(P, i, n){
	G <- P
	H <- P[i,]
	E <- 1 - diag(size(P)[2])
	for (m in 2:n) {
		G <- P %*% (G * E)
		H <- rbind(H, G[i,])
	}
	return(H)
}
firstPassagePdF <- firstPassage(object = mcWeatherL, state = "sunny",
				n = 10)
print( firstPassagePdF[3, 3] )

#=================Statistical Analysis=======================
weathersOfDays <- rmarkovchain(n = 365, object = mcWeatherL, t0 = "sunny")
print( weathersOfDays[1:30] )
patientStates <- rmarkovchain(n = 5, object = mcCCRC, t0 = "H",
			      include.t0 = TRUE)
print( patientStates[1:10,] )

weatherFittedMLE <- markovchainFit(data = weathersOfDays, method = "mle",name = "Weather MLE")
print(weatherFittedMLE$estimate)
print(weatherFittedMLE$standardError)

weatherFittedLAPLACE <- markovchainFit(data = weathersOfDays,method = "laplace", laplacian = 0.01,name = "Weather LAPLACE")
print(weatherFittedLAPLACE$estimate)
print(createSequenceMatrix(stringchar = weathersOfDays))

myMatr<-matrix(c("a","b","b","a","a","b","b","b","b","a","a","a","b","a"),ncol=2)
print(createSequenceMatrix(stringchar = myMatr,toRowProbs = TRUE))

weatherFittedBOOT <- markovchainFit(data = weathersOfDays, method = "bootstrap", nboot = 20)
print(weatherFittedBOOT$estimate)
print(weatherFittedBOOT$standardError)
weatherFittedBOOTParallel <- markovchainFit(data = weathersOfDays,method = "bootstrap", nboot = 200,parallel = TRUE)
print(weatherFittedBOOTParallel$estimate)
print(weatherFittedBOOTParallel$standardError)



RcppParallel::setThreadOptions(numThreads = 2)

print(weatherFittedMLE$logLikelihood)
print(weatherFittedBOOT$logLikelihood)

print(weatherFittedMLE$confidenceInterval)
print(weatherFittedBOOT$confidenceInterval)

###fix this pararell
#RcppParallel::setNumThreads(2)
#print(weatherFittedMLE$logLikelihood)
#print(weatherFittedBOOT$logLikelihood)

multinomialConfidenceIntervals(transitionMatrix = weatherFittedMLE$estimate@transitionMatrix, countsTransitionMatrix = createSequenceMatrix(weathersOfDays))

data(holson)
singleMc<-markovchainFit(data=holson[,2:12],name="holson")
mcListFit<-markovchainListFit(data=holson[,2:6],name="holson")
print(mcListFit)


c1 <- c("a", "b", "a", "a", "c", "c", "a")
c2 <- c("b")
c3 <- c("c", "a", "a", "c")
c4 <- c("b", "a", "b", "a", "a", "c", "b")
c5 <- c("a", "a", "c", NA)
c6 <- c("b", "c", "b", "c", "a")
myList <- list(c1, c2, c3, c4, c5, c6)
myListMc <- markovchainFit(data = myList)
print(myListMc)


# lo anterior también trabaja 
print(markovchainListFit(data = myList))


## PREDICCIÓN
# predecir desde un objeto markovchain
predict(object = weatherFittedMLE$estimate, newdata = c("cloudy", "sunny"),
        n.ahead = 3)

#predecir desde markovchainList
predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5)


predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5, continue = TRUE)



