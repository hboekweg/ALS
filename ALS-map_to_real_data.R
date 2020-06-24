#categorize data
library(cellAlign)

#### Loading the data ####
data(expGlobalLPS)
data(expGlobalPAM)
data(expLocalLPS)
data(expLocalPAM)
data(trajLPS)
data(trajPAM)

trajLPS = as.numeric(t(trajLPS))
names(trajLPS) = colnames(expGlobalLPS)
trajPAM = as.numeric(t(trajPAM))
names(trajPAM) = colnames(expGlobalPAM)

numPts = 200
interGlobalLPS = cellAlign::interWeights(expDataBatch = expGlobalLPS, trajCond = trajLPS,
                                         winSz = 0.1, numPts = numPts)
interGlobalPAM = cellAlign::interWeights(expDataBatch = expGlobalPAM, trajCond = trajPAM,
                                         winSz = 0.1, numPts = numPts)

require(ggplot2)
require(reshape2)
require(pheatmap)
sharedMarkers = intersect(rownames(expGlobalLPS), rownames(expGlobalPAM))
whichgene=sharedMarkers[1]
selectedLPS<-interGlobalLPS$interpolatedVals[whichgene,]
selectedPAM<-interGlobalPAM$interpolatedVals[whichgene,]

dfLPSi = data.frame(traj = interGlobalLPS$traj, value=(selectedLPS), error=interGlobalLPS$error[whichgene,])
dfLPS = data.frame(traj = trajLPS, t(expGlobalLPS[whichgene,]))
dfPAMi = data.frame(traj = interGlobalPAM$traj, value=(selectedPAM), error=interGlobalPAM$error[whichgene,])
dfPAM = data.frame(traj = trajPAM, t(expGlobalPAM[whichgene,]))
dfLPSM = melt(dfLPS, id.vars = 'traj')
dfPAMM = melt(dfPAM, id.vars = 'traj')
#plot of an example gene and its interpolation with error bars
ggplot(dfLPSi, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=dfLPSM, aes(x=traj,y=value)) + ggtitle(whichgene) 


####Figure out slope function ####
#figure out what the slope is (I'm pretty sure that everything in this dataset has a fold chage of 2--slope of 1)
Y = dfLPSi$value
X = dfLPSi$traj
#For finding the slope: y2-y1/x2-x1
#get the end points of the line. The way below assumes that the arrays are sorted. 
y1 = Y[1]
y2 = Y[length(Y)]
x1 = X[1]
x2 = X[length(X)]
m = (y2-y1)/(x2-x1)


####Figure out standard deviation function ####
#now find the standard deviation
#let look at the error bars
errors = dfLPSi$error
hist(errors, breaks = 20)
median(errors)
sd(errors)

#find residuals and put them in a distribution
#Right now x and Y are the points of the line. We just want the regular points
pointX = dfLPSM$traj
pointY = dfLPSM$value

data = data.frame(pointX, pointY)
ggplot(data, aes(x = pointX, y = pointY)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

plot(pointX, pointY)
m1 <- lm(pointY~pointX)  #Create a linear model
resid(m1) #List of residuals
plot(density(resid(m1))) #A density plot
hist(resid(m1), breaks = 20)
plot(resid(m1))

summary(m1)

data$predicted <- predict(m1) 
data$residuals <- residuals(m1)

####Make Functions ####
#x and y coordinates of the interpolated line
findSlope = function(pointX, pointY){
  data = data.frame(pointX, pointY)
  m1 <- lm(pointY~pointX)
  m = m1$coefficients[2]
  return(m)
}

#the paramters are the data points (not interplated line)
findSD = function(pointX, pointY){
  data = data.frame(pointX, pointY)
  m1 <- lm(pointY~pointX)
  data$predicted <- predict(m1) 
  data$residuals <- residuals(m1)
  sd = sd(data$residuals)
  return(sd)
  
}


####Loope through all the genes and find slope and standard deviation ####
#loop through every gene in the data table 
#get all the gene names

LPSgenes <- vector(mode="list", length=length(sharedMarkers))
names(LPSgenes) <- sharedMarkers
position_in_list = 1
for(gene in sharedMarkers){
  whichgene=sharedMarkers[position_in_list]
  selectedLPS<-interGlobalLPS$interpolatedVals[whichgene,]
  # selectedPAM<-interGlobalPAM$interpolatedVals[whichgene,]
  
  dfLPSi = data.frame(traj = interGlobalLPS$traj, value=(selectedLPS), error=interGlobalLPS$error[whichgene,])
  dfLPS = data.frame(traj = trajLPS, t(expGlobalLPS[whichgene,]))
  # dfPAMi = data.frame(traj = interGlobalPAM$traj, value=(selectedPAM), error=interGlobalPAM$error[whichgene,])
  # dfPAM = data.frame(traj = trajPAM, t(expGlobalPAM[whichgene,]))
  dfLPSM = melt(dfLPS, id.vars = 'traj')
  # dfPAMM = melt(dfPAM, id.vars = 'traj')
  
  #find slope
  Y = dfLPSi$value
  X = dfLPSi$traj
  slope = findSlope(X, Y)
  
  
  #find sd
  pointX = dfLPSM$traj
  pointY = dfLPSM$value
  sd = findSD(pointX, pointY)
  
  LPSgenes[[position_in_list]] = c(slope, sd)
  position_in_list = position_in_list + 1
  
}


PAMgenes <- vector(mode="list", length=length(sharedMarkers))
names(PAMgenes) <- sharedMarkers
position_in_list = 1
for(gene in sharedMarkers){
  whichgene=sharedMarkers[position_in_list]
  # selectedLPS<-interGlobalLPS$interpolatedVals[whichgene,]
  selectedPAM<-interGlobalPAM$interpolatedVals[whichgene,]
  
  # dfLPSi = data.frame(traj = interGlobalLPS$traj, value=(selectedLPS), error=interGlobalLPS$error[whichgene,])
  # dfLPS = data.frame(traj = trajLPS, t(expGlobalLPS[whichgene,]))
  dfPAMi = data.frame(traj = interGlobalPAM$traj, value=(selectedPAM), error=interGlobalPAM$error[whichgene,])
  dfPAM = data.frame(traj = trajPAM, t(expGlobalPAM[whichgene,]))
  # dfLPSM = melt(dfLPS, id.vars = 'traj')
  dfPAMM = melt(dfPAM, id.vars = 'traj')
  
  #find slope
  Y = dfPAMi$value
  X = dfPAMi$traj
  slope = findSlope(X, Y)
  
  
  #find sd
  pointX = dfPAMM$traj
  pointY = dfPAMM$value
  sd = findSD(pointX, pointY)
  
  PAMgenes[[position_in_list]] = c(slope, sd)
  position_in_list = position_in_list + 1
  
}

LPSslopes = c()
LPSdeviations = c()
i = 1
for(gene in LPSgenes){
  # append(myslopes, gene[1])
  LPSslopes[i] = gene[1]
  LPSdeviations[i] = gene[2]
  i = i + 1
  }

PAMslopes = c()
PAMdeviations = c()
i = 1
for(gene in PAMgenes){
  # append(myslopes, gene[1])
  PAMslopes[i] = gene[1]
  PAMdeviations[i] = gene[2]
  i = i + 1
}

#### save for plotting purposes ####
slope_and_sd = data.frame("LPSslope"=LPSslopes, "LPSdeviation"=LPSdeviations, "PAMslope"=PAMslopes, "PAMdeviation"=PAMdeviations)
dir.create("mapping")
write.csv(slope_and_sd, "mapping/cellAlign")

####Generate fake data to test functions####
generate_data = function(slope, sd){
  X = runif(1000)
  X = as.data.frame(sort(X))
  X = t(X)
  
  Ynorm <-  slope * X +1+ rnorm(1000,sd = sd)
  Yflat = rep(mean(Ynorm), 1000)
  Yflat = t(Yflat)
  colnames(Ynorm) = paste0("t", 1:1000)
  rownames(Ynorm) = "gene"
  colnames(Yflat) = paste0("t", 1:1000)
  rownames(Yflat) = "gene"
  
  X = t(X)
  rownames(X) = paste0("t", 1:1000)
  colnames(X) = NULL
  
  data = list("time_traj" = X, "Y" = Ynorm, "flat_Y"=Yflat)
  
  
  
}
interpolate = function(time_traj, expression, subsample=TRUE){
  Y = expression
  X = time_traj
  interGlobalLPS = cellAlign::interWeights(expDataBatch = Y, trajCond =X,
                                           winSz = 0.1, numPts = 200)
  require(reshape2)
  whichgene = "gene"
  selectedLPS<-interGlobalLPS$interpolatedVals[whichgene,]
  
  dfLPSi = data.frame(traj = interGlobalLPS$traj, value=(selectedLPS), error=interGlobalLPS$error[whichgene,])
  
  ###FIX THIS###
  if(subsample==TRUE){
    dfLPS = data.frame(traj = X, gene=t(Y[whichgene,]))
  }else{
    dfLPS = data.frame(traj = X, gene=Y[whichgene,])
  }
  dfLPSM = melt(dfLPS, id.vars = 'traj')
  
  result = list("interpolated" = dfLPSi, "points"=dfLPSM)
  return(result)
}

test_data = generate_data(slope=1, sd=0)
pointX = test_data$time_traj
pointY = t(test_data$Y)
findSD(pointX, pointY)

#need to give it the end points of the interpolated line
inter = interpolate(test_data$time_traj, test_data$Y, subsample=FALSE) 
X = inter$interpolated$traj
Y = inter$interpolated$value
findSlope(X, Y)




#### PLOTTING####
test_data = generate_data(slope=1, sd=.75)
pointX = test_data$time_traj
pointY = t(test_data$Y)
m1 <- lm(pointY~pointX)  #Create a linear model
test_data$predicted <- predict(m1) 
test_data$residuals <- residuals(m1)
test_data = as.data.frame(test_data)

#plot residuals
# Changing alpha of actual values based on absolute value of residuals
library(ggplot2)
ggplot(test_data, aes(x = pointX, y = pointY)) +
  geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +
  geom_segment(aes(xend = pointX, yend = predicted), alpha = .2) +
  
  # > Alpha adjustments made here...
  geom_point(aes(alpha = abs(residuals))) +  # Alpha mapped to abs(residuals)
  guides(alpha = FALSE) +  # Alpha legend removed
  # <
  
  geom_point(aes(y = predicted), shape = 1) +
  theme_bw()

hist(test_data$residuals, breaks = 20)
sd(test_data$residuals)


