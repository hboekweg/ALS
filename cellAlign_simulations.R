library(devtools)
library(cellAlign)


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


subsample = function(oracle, number_of_cells){
  traj = oracle$time_traj
  exp = oracle$Y

  
  trajdf = as.data.frame(traj)
  
  #divide time into quarters
  q = split(as.vector(trajdf), factor(sort(rank(as.vector(trajdf))%%4)))
  first = q$`0`
  second = q$`1`
  third = q$`2`
  fourth = q$`3`
  
  #get an even number of sample from each
  remainder = number_of_cells %% 4
  s_size = number_of_cells %/% 4
  s_size_plus1 = s_size 
  s_size_plus2 = s_size 
  s_size_plus3 = s_size 
  
  if (remainder == 1){
    s_size_plus1 = s_size_plus1 + 1 
  }
  if (remainder == 2){
    s_size_plus1 = s_size_plus1 + 1 
    s_size_plus2 =  s_size_plus2+1 
  }
  if (remainder == 3){
    s_size_plus1 = s_size_plus1 + 1 
    s_size_plus2 =  s_size_plus2+1 
    s_size_plus3 = s_size_plus3 +1 
  }
  
  first_index = rownames(first)
  first_sample = sample(first_index, s_size_plus1)
  
  second_index = rownames(second)
  second_sample = sample(second_index, s_size_plus2)
  
  third_index = rownames(third)
  third_sample = sample(third_index, s_size_plus3)
  
  fourth_index = rownames(fourth)
  fourth_sample = sample(fourth_index, s_size)
  
  all_samples = c(first_sample, second_sample, third_sample, fourth_sample)
  
  #subset traj
  trajdf <- cbind(index = rownames(trajdf), trajdf)
  rownames(trajdf) <- 1:nrow(trajdf)
  traj_subset = trajdf[which((trajdf$index %in% all_samples)==TRUE),]
  rownames(traj_subset) <- traj_subset$index
  traj_subset$index=NULL
  traj_subset = as.matrix(traj_subset)
  colnames(traj_subset) = NULL
  
  
  #subset exp
  exp_subset = exp[ ,which((colnames(exp) %in% all_samples)==TRUE)]
  exp_subset = t(exp_subset)
  rownames(exp_subset)="gene"
  exp_subset = as.data.frame(exp_subset)
  
  subsample = list("traj"=traj_subset, "exp"=exp_subset)
  return(subsample)
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


find_area = function(allx, curve1Y, curve2Y, curve1X, curve2X, my_sample){
  library(geiger)
  library(sp)
  library(rgeos)
  
  #For finding interstections
  # x = paper_results$LPSi$traj
  f1 = curve1Y #paper_results$LPSi$value
  f2 = curve2Y #sample15$LPSi$value
  
  y1 = curve1Y #paper_results$LPSi$value
  y2 = curve2Y #sample15$LPSi$value
  
  # y data 
  x1 <- curve1X #paper_results$LPSi$traj
  x2 <- curve2X #sample15$LPSi$traj
  
  # convert to a sp object (spatial lines)
  l1 <- Line(matrix(c(x1, y1), nc = 2, byrow = F))
  l2 <- Line(matrix(c(x2, y2), nc = 2, byrow = F))
  ll1 <- Lines(list(l1), ID = "1")
  ll2 <- Lines(list(l2), ID = "1")
  sl1 <- SpatialLines(list(ll1), proj4string = CRS("+init=epsg:4269"))
  sl2 <- SpatialLines(list(ll2), proj4string = CRS("+init=epsg:4269"))
  # Calculate locations where spatial lines intersect
  if(is.null(gIntersection(sl1, sl2, byid = TRUE))){
    smallest = min(my_sample$interpolated$traj)
    largest = max(my_sample$interpolated$traj)
    x = c()
    x = append(x, largest, after = length(x))
    x = append(x, smallest, after = length(x))
    x = sort(x, decreasing = FALSE)
    
  }
  else{
    int.pts <- gIntersection(sl1, sl2, byid = TRUE)
    int.coords <- int.pts@coords
    
    # Add end points
    smallest = min(my_sample$interpolated$traj)
    largest = max(my_sample$interpolated$traj)
    x = int.coords[,1]
    x = append(x, largest, after = length(x))
    x = append(x, smallest, after = length(x))
    x = sort(x, decreasing = FALSE)
    
    
  }
  
  total_area = 0
  for (i in c(1:length(x))){
    first = x[i]
    second = x[i+1]
    if(is.na(second)){
      break
    }
    else{
      area = geiger:::.area.between.curves(allx, f1, f2, xrange = c(first, second))
      area = abs(area)
      total_area = total_area + area
    }
  }
  total_area
  
  return(total_area)
}



#CREATE LINES

#flatline (average of all Y values)
data = generate_data(1, 0)


inter = interpolate(data$time_traj, data$Y, subsample=FALSE)
ggplot(inter$interpolated, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=inter$points, aes(x=traj,y=value))

inter_flat = interpolate(data$time_traj, data$flat_Y, subsample=FALSE)
ggplot(inter_flat$interpolated, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=inter_flat$points, aes(x=traj,y=value))


# data = generate_data(1, 0)
mysample = subsample(data, 16)
inter_S = interpolate(mysample$traj, mysample$exp, subsample=TRUE)
ggplot(inter_S$interpolated, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=inter_S$points, aes(x=traj,y=value))

#plot them together
ggplot(inter_S$interpolated, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + 
  geom_line(size = 2) + geom_point(data=inter_S$points, aes(x=traj,y=value)) + 
  geom_line(data=inter_flat$interpolated, color='red') +
  geom_line(data=inter$interpolated, color='blue')


#FIND AREA
find_area(inter$interpolated$traj, inter$interpolated$value, inter_S$interpolated$value, inter$interpolated$traj, inter_S$interpolated$traj, inter_S)
#we're going to have to check ABC for flat vs sample
find_area(inter_flat$interpolated$traj, inter_flat$interpolated$value, inter_S$interpolated$value, inter_flat$interpolated$traj, inter_S$interpolated$traj, inter_S)


accuracy = function(slope, sd, number_of_cells, number_of_subsamples){
  data = generate_data(slope, sd)
  # inter = interpolate(data$time_traj, data$Y, subsample=FALSE) This line is needed only for plotting
  #get flat line
  inter_flat = interpolate(data$time_traj, data$flat_Y, subsample=FALSE)
  
  oracle_area = c()
  flat_area = c()
  for(i in 1:number_of_subsamples){
    mysample = subsample(data, number_of_cells)
    inter_S = interpolate(mysample$traj, mysample$exp, subsample=TRUE)
    
    
    oracle = find_area(inter$interpolated$traj, inter$interpolated$value, inter_S$interpolated$value, inter$interpolated$traj, inter_S$interpolated$traj, inter_S)
    oracle_area[i] = oracle
    flat = find_area(inter_flat$interpolated$traj, inter_flat$interpolated$value, inter_S$interpolated$value, inter_flat$interpolated$traj, inter_S$interpolated$traj, inter_S)
    flat_area[i] = flat
  }

  # areas = list("oracle"= oracle_area, "flat" = flat_line_area)
  areas = list("oracle_area"=oracle_area, "flat_area"=flat_area)
  return(areas)
}

sd_1 = accuracy(slope=1, sd = 1, number_of_cells = 16, 3)


write.csv(areas, "sd_1")






