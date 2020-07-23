library(cellAlign)
require(ggplot2)
generate_data = function(slope, sd){
  X = runif(100000)
  X = as.data.frame(sort(X))
  X = t(X)
  
  Ynorm <-  slope * X +1+ rnorm(100000,sd = sd)
  Yflat = rep(mean(Ynorm), 100000)
  Yflat = t(Yflat)
  colnames(Ynorm) = paste0("t", 1:100000)
  rownames(Ynorm) = "gene"
  colnames(Yflat) = paste0("t", 1:100000)
  rownames(Yflat) = "gene"
  
  X = t(X)
  rownames(X) = paste0("t", 1:100000)
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
  #if there is no area to be caculated (ex when we are calculating slope of stright line with sd of 0)
  if(all.equal(y1,y2)==TRUE){
    total_area = 0
    return(total_area)
  }
  #if there are no interecting points
  else if(is.null(gIntersection(sl1, sl2, byid = TRUE))){
    smallest = min(my_sample$interpolated$traj)
    largest = max(my_sample$interpolated$traj)
    x = c()
    x = append(x, largest, after = length(x))
    x = append(x, smallest, after = length(x))
    x = sort(x, decreasing = FALSE)
    
  }
  #if there are intersecting points
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

accuracy = function(slope, sd, number_of_cells, number_of_subsamples){
  data = generate_data(slope, sd)
  inter = interpolate(data$time_traj, data$Y, subsample=FALSE) 
  #get flat line
  inter_flat = interpolate(data$time_traj, data$flat_Y, subsample=FALSE)
  
  oracle_area = c()
  flat_area = c()
  # number_of_subsamples = number_of_subsamples
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

accuracy_perc = function(oracle_area, flat_area){
  difference = flat_area - oracle_area
  neg = sum(difference < 0) #flat wins
  pos = sum(difference > 0) #truth wins
  accuracy = pos/length(difference)
  return(accuracy)
}


simulate = function(number_of_cells, file_path){
   #### 7CELLS ####
  #cell=7, slope=0
  #cell=7, slope=.5
  sd0 = accuracy(slope=.5, sd = 0, number_of_cells, number_of_subsamples=1000)
  sd.25 = accuracy(slope=.5, sd = .25, number_of_cells, number_of_subsamples=1000)
  sd.5 = accuracy(slope=.5, sd =.5, number_of_cells, number_of_subsamples=1000)
  sd.7 = accuracy(slope=.5, sd =.7, number_of_cells, number_of_subsamples=1000)
  sd1 = accuracy(slope=.5, sd = 1, number_of_cells, number_of_subsamples=1000)
  accuracy_0 = accuracy_perc(sd0$oracle_area, sd0$flat_area)
  accuracy_.25 = accuracy_perc(sd.25$oracle_area, sd.25$flat_area)
  accuracy_.5 = accuracy_perc(sd.5$oracle_area, sd.5$flat_area)
  accuracy_.7 = accuracy_perc(sd.7$oracle_area, sd.7$flat_area)
  accuracy_1 = accuracy_perc(sd1$oracle_area, sd1$flat_area)
  df_slope.5 = data.frame("deviation"=c(0, .25, .5,.7,1), "slope"=c(.5,.5,.5,.5,.5), "accuracy"=c(accuracy_0,accuracy_.25,accuracy_.5, accuracy_.7,accuracy_1))

  #cells=7, slope=1
  sd0 = accuracy(slope=1, sd = 0, number_of_cells, number_of_subsamples=1000)
  sd.25 = accuracy(slope=1, sd = .25, number_of_cells, number_of_subsamples=1000)
  sd.5 = accuracy(slope=1, sd =.5, number_of_cells, number_of_subsamples=1000)
  sd.7 = accuracy(slope=1, sd =.7, number_of_cells, number_of_subsamples=1000)
  sd1 = accuracy(slope=1, sd = 1, number_of_cells, number_of_subsamples=1000)
  accuracy_0 = accuracy_perc(sd0$oracle_area, sd0$flat_area)
  accuracy_.25 = accuracy_perc(sd.25$oracle_area, sd.25$flat_area)
  accuracy_.5 = accuracy_perc(sd.5$oracle_area, sd.5$flat_area)
  accuracy_.7 = accuracy_perc(sd.7$oracle_area, sd.7$flat_area)
  accuracy_1 = accuracy_perc(sd1$oracle_area, sd1$flat_area)
  df_slope1 = data.frame("deviation"=c(0, .25, .5,.7,1), "slope"=c(1,1,1,1,1), "accuracy"=c(accuracy_0,accuracy_.25,accuracy_.5, accuracy_.7,accuracy_1))

  #cell=7, slope=2
  sd0 = accuracy(slope=2, sd = 0, number_of_cells, number_of_subsamples=1000)
  sd.25 = accuracy(slope=2, sd = .25, number_of_cells, number_of_subsamples=1000)
  sd.5 = accuracy(slope=2, sd =.5, number_of_cells, number_of_subsamples=1000)
  sd.7 = accuracy(slope=2, sd =.7, number_of_cells, number_of_subsamples=1000)
  sd1 = accuracy(slope=2, sd = 1, number_of_cells, number_of_subsamples=1000)
  accuracy_0 = accuracy_perc(sd0$oracle_area, sd0$flat_area)
  accuracy_.25 = accuracy_perc(sd.25$oracle_area, sd.25$flat_area)
  accuracy_.5 = accuracy_perc(sd.5$oracle_area, sd.5$flat_area)
  accuracy_.7 = accuracy_perc(sd.7$oracle_area, sd.7$flat_area)
  accuracy_1 = accuracy_perc(sd1$oracle_area, sd1$flat_area)
  df_slope2 = data.frame("deviation"=c(0, .25, .5,.7,1), "slope"=c(2,2,2,2,2), "accuracy"=c(accuracy_0,accuracy_.25,accuracy_.5, accuracy_.7,accuracy_1))


  #cell=7, slope=4
  sd0 = accuracy(slope=4, sd = 0, number_of_cells, number_of_subsamples=1000)
  sd.25 = accuracy(slope=4, sd = .25, number_of_cells, number_of_subsamples=1000)
  sd.5 = accuracy(slope=4, sd =.5, number_of_cells, number_of_subsamples=1000)
  sd.7 = accuracy(slope=4, sd =.7, number_of_cells, number_of_subsamples=1000)
  sd1 = accuracy(slope=4, sd = 1, number_of_cells, number_of_subsamples=1000)
  accuracy_0 = accuracy_perc(sd0$oracle_area, sd0$flat_area)
  accuracy_.25 = accuracy_perc(sd.25$oracle_area, sd.25$flat_area)
  accuracy_.5 = accuracy_perc(sd.5$oracle_area, sd.5$flat_area)
  accuracy_.7 = accuracy_perc(sd.7$oracle_area, sd.7$flat_area)
  accuracy_1 = accuracy_perc(sd1$oracle_area, sd1$flat_area)
  df_slope4 = data.frame("deviation"=c(0, .25, .5,.7,1), "slope"=c(4,4,4,4,4), "accuracy"=c(accuracy_0,accuracy_.25,accuracy_.5, accuracy_.7,accuracy_1))


  # cells_7 = do.call("rbind", list(df_slope.5, df_slope1, df_slope2,df_slope4))
  cells_df = do.call("rbind", list(df_slope.5, df_slope1, df_slope2,df_slope4))
  write.csv(cells_df, file_path)
  return(cells_df)

}

dir.create("ten_reps_fig_1")
dir.create("ten_reps_fig_1/cells_7")
rep1 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep1")
rep2 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep2")
rep3 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep3")
rep4 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep4")
rep5 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep5")
rep6 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep6")
rep7 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep7")
rep8 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep8")
rep9 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep9")
rep10 = simulate(number_of_cells=7, "ten_reps_fig_1/cells_7/rep10")

dir.create("ten_reps_fig_1/cells_16")
rep1 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep1")
rep2 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep2")
rep3 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep3")
rep4 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep4")
rep5 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep5")
rep6 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep6")
rep7 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep7")
rep8 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep8")
rep9 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep9")
rep10 = simulate(number_of_cells=16, "ten_reps_fig_1/cells_16/rep10")

dir.create("ten_reps_fig_1/cells_20")
rep1 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep1")
rep2 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep2")
rep3 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep3")
rep4 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep4")
rep5 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep5")
rep6 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep6")
rep7 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep7")
rep8 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep8")
rep9 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep9")
rep10 = simulate(number_of_cells=20, "ten_reps_fig_1/cells_20/rep10")

dir.create("ten_reps_fig_1/cells_30")
rep1 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep1")
rep2 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep2")
rep3 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep3")
rep4 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep4")
rep5 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep5")
rep6 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep6")
rep7 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep7")
rep8 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep8")
rep9 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep9")
rep10 = simulate(number_of_cells=30, "ten_reps_fig_1/cells_30/rep10")


dir.create("ten_reps_fig_1/cells_100")
rep1 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep1")
rep2 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep2")
rep3 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep3")
rep4 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep4")
rep5 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep5")
rep6 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep6")
rep7 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep7")
rep8 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep8")
rep9 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep9")
rep10 = simulate(number_of_cells=100, "ten_reps_fig_1/cells_100/rep10")


dir.create("ten_reps_fig_1/cells_1000")
rep1 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep1")
rep2 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep2")
rep3 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep3")
rep4 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep4")
rep5 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep5")
rep6 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep6")
rep7 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep7")
rep8 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep8")
rep9 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep9")
rep10 = simulate(number_of_cells=1000, "ten_reps_fig_1/cells_1000/rep10")

#3 dim array: array(1:23, dim(rows, columns, table))
#you can assign dimensions to exist arrays by: dim(my.vector) <- c(3,4,2)






