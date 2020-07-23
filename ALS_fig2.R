library(cellAlign)
require(ggplot2)
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

simulate = function(slope, sd, number_of_cells){
  
  combo_1 = accuracy(slope, sd, number_of_cells, number_of_subsamples=1000)
  combo_2 = accuracy(slope, sd, number_of_cells, number_of_subsamples=1000)
  combo_3 = accuracy(slope, sd, number_of_cells, number_of_subsamples=1000)
  combo_4 = accuracy(slope, sd, number_of_cells, number_of_subsamples=1000)
  combo_5 = accuracy(slope, sd, number_of_cells, number_of_subsamples=1000)
  accuracy_1 = accuracy_perc(combo_1$oracle_area, combo_1$flat_area)
  accuracy_2 = accuracy_perc(combo_2$oracle_area, combo_2$flat_area)
  accuracy_3 = accuracy_perc(combo_3$oracle_area, combo_3$flat_area)
  accuracy_4 = accuracy_perc(combo_4$oracle_area, combo_4$flat_area)
  accuracy_5 = accuracy_perc(combo_5$oracle_area, combo_5$flat_area)
  
  combo_df = data.frame("deviation"=c(sd), "slope"=c(slope), 
                          "accuracy"=c(accuracy_1, accuracy_2, accuracy_3, accuracy_4, accuracy_5))
  
  return(combo_df)
  
}

dir.create("ten_reps_fig_2_second_run")

####Slope/var=.5####
#slope .5 ::  .5/1   .75/1.5  1/2  1.5/3 2/4   3/6
combo_1 = simulate(slope=.5, sd=1, number_of_cells=7)
combo_2 = simulate(slope=.75, sd=1.5, number_of_cells=7)
combo_3 = simulate(slope=1, sd=2, number_of_cells=7)
combo_4 = simulate(slope=2, sd=4, number_of_cells=7)
combo_5 = simulate(slope=3, sd=6, number_of_cells=7)
cells_7 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_7$slope_over_var = ".5"
cells_7$cells = "cells_7"

combo_1 = simulate(slope=.5, sd=1, number_of_cells=16)
combo_2 = simulate(slope=.75, sd=1.5, number_of_cells=16)
combo_3 = simulate(slope=1, sd=2, number_of_cells=16)
combo_4 = simulate(slope=2, sd=4, number_of_cells=16)
combo_5 = simulate(slope=3, sd=6, number_of_cells=16)
cells_16 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_16$slope_over_var = ".5"
cells_16$cells = "cells_16"

combo_1 = simulate(slope=.5, sd=1, number_of_cells=20)
combo_2 = simulate(slope=.75, sd=1.5, number_of_cells=20)
combo_3 = simulate(slope=1, sd=2, number_of_cells=20)
combo_4 = simulate(slope=2, sd=4, number_of_cells=20)
combo_5 = simulate(slope=3, sd=6, number_of_cells=20)
cells_20 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_20$slope_over_var = ".5"
cells_20$cells = "cells_20"

combo_1 = simulate(slope=.5, sd=1, number_of_cells=30)
combo_2 = simulate(slope=.75, sd=1.5, number_of_cells=30)
combo_3 = simulate(slope=1, sd=2, number_of_cells=30)
combo_4 = simulate(slope=2, sd=4, number_of_cells=30)
combo_5 = simulate(slope=3, sd=6, number_of_cells=30)
cells_30 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_30$slope_over_var = ".5"
cells_30$cells = "cells_30"

slope_var_.5 = do.call("rbind", list(cells_7, cells_16, cells_20 ,cells_30))

write.csv(slope_var_.5, "ten_reps_fig_2_second_run/slope_var_.5")

####Slope/var=1####
#.5/.5  1/1  1.75/1.75  2/2   3/3
combo_1 = simulate(slope=.5, sd=.5, number_of_cells=7)
combo_2 = simulate(slope=1, sd=1, number_of_cells=7)
combo_3 = simulate(slope=1.75, sd=1.75, number_of_cells=7)
combo_4 = simulate(slope=2, sd=2, number_of_cells=7)
combo_5 = simulate(slope=3, sd=3, number_of_cells=7)
cells_7 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_7$slope_over_var = "1"
cells_7$cells = "cells_7"

combo_1 = simulate(slope=.5, sd=.5, number_of_cells=16)
combo_2 = simulate(slope=1, sd=1, number_of_cells=16)
combo_3 = simulate(slope=1.75, sd=1.75, number_of_cells=16)
combo_4 = simulate(slope=2, sd=2, number_of_cells=16)
combo_5 = simulate(slope=3, sd=3, number_of_cells=16)
cells_16 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_16$slope_over_var = "1"
cells_16$cells = "cells_16"

combo_1 = simulate(slope=.5, sd=.5, number_of_cells=20)
combo_2 = simulate(slope=1, sd=1, number_of_cells=20)
combo_3 = simulate(slope=1.75, sd=1.75, number_of_cells=20)
combo_4 = simulate(slope=2, sd=2, number_of_cells=20)
combo_5 = simulate(slope=3, sd=3, number_of_cells=20)
cells_20 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_20$slope_over_var = "1"
cells_20$cells = "cells_20"

combo_1 = simulate(slope=.5, sd=.5, number_of_cells=30)
combo_2 = simulate(slope=1, sd=1, number_of_cells=30)
combo_3 = simulate(slope=1.75, sd=1.75, number_of_cells=30)
combo_4 = simulate(slope=2, sd=2, number_of_cells=30)
combo_5 = simulate(slope=3, sd=3, number_of_cells=30)
cells_30 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_30$slope_over_var = "1"
cells_30$cells = "cells_30"

slope_var_1 = do.call("rbind", list(cells_7, cells_16, cells_20 ,cells_30))

write.csv(slope_var_1, "ten_reps_fig_2_second_run/slope_var_1")

####Slope/var=1.5####
# .75/.5  2.25/1.5  1.5/1   3/2   4.5/3  
combo_1 = simulate(slope=.75, sd=.5, number_of_cells=7)
combo_2 = simulate(slope=2.25, sd=1.5, number_of_cells=7)
combo_3 = simulate(slope=1.5, sd=1, number_of_cells=7)
combo_4 = simulate(slope=3, sd=2, number_of_cells=7)
combo_5 = simulate(slope=4.5, sd=3, number_of_cells=7)
cells_7 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_7$slope_over_var = "1.5"
cells_7$cells = "cells_7"

combo_1 = simulate(slope=.75, sd=.5, number_of_cells=16)
combo_2 = simulate(slope=2.25, sd=1.5, number_of_cells=16)
combo_3 = simulate(slope=1.5, sd=1, number_of_cells=16)
combo_4 = simulate(slope=3, sd=2, number_of_cells=16)
combo_5 = simulate(slope=4.5, sd=3, number_of_cells=16)
cells_16 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_16$slope_over_var = "1.5"
cells_16$cells = "cells_16"

combo_1 = simulate(slope=.75, sd=.5, number_of_cells=20)
combo_2 = simulate(slope=2.25, sd=1.5, number_of_cells=20)
combo_3 = simulate(slope=1.5, sd=1, number_of_cells=20)
combo_4 = simulate(slope=3, sd=2, number_of_cells=20)
combo_5 = simulate(slope=4.5, sd=3, number_of_cells=20)
cells_20 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_20$slope_over_var = "1.5"
cells_20$cells = "cells_20"

combo_1 = simulate(slope=.75, sd=.5, number_of_cells=30)
combo_2 = simulate(slope=2.25, sd=1.5, number_of_cells=30)
combo_3 = simulate(slope=1.5, sd=1, number_of_cells=30)
combo_4 = simulate(slope=3, sd=2, number_of_cells=30)
combo_5 = simulate(slope=4.5, sd=3, number_of_cells=30)
cells_30 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_30$slope_over_var = "1.5"
cells_30$cells = "cells_30"

slope_var_15 = do.call("rbind", list(cells_7, cells_16, cells_20 ,cells_30))

write.csv(slope_var_15, "ten_reps_fig_2_second_run/slope_var_15")

####slope/var=2####
# 1/.5    1.5/.75   2/1   4/2   6/3
combo_1 = simulate(slope=1, sd=.5, number_of_cells=7)
combo_2 = simulate(slope=1.5, sd=.75, number_of_cells=7)
combo_3 = simulate(slope=2, sd=1, number_of_cells=7)
combo_4 = simulate(slope=4, sd=2, number_of_cells=7)
combo_5 = simulate(slope=6, sd=3, number_of_cells=7)
cells_7 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_7$slope_over_var = "2"
cells_7$cells = "cells_7"

combo_1 = simulate(slope=1, sd=.5, number_of_cells=16)
combo_2 = simulate(slope=1.5, sd=.75, number_of_cells=16)
combo_3 = simulate(slope=2, sd=1, number_of_cells=16)
combo_4 = simulate(slope=4, sd=2, number_of_cells=16)
combo_5 = simulate(slope=6, sd=3, number_of_cells=16)
cells_16 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_16$slope_over_var = "2"
cells_16$cells = "cells_16"

combo_1 = simulate(slope=1, sd=.5, number_of_cells=20)
combo_2 = simulate(slope=1.5, sd=.75, number_of_cells=20)
combo_3 = simulate(slope=2, sd=1, number_of_cells=20)
combo_4 = simulate(slope=4, sd=2, number_of_cells=20)
combo_5 = simulate(slope=6, sd=3, number_of_cells=20)
cells_20 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_20$slope_over_var = "2"
cells_20$cells = "cells_20"

combo_1 = simulate(slope=1, sd=.5, number_of_cells=30)
combo_2 = simulate(slope=1.5, sd=.75, number_of_cells=30)
combo_3 = simulate(slope=2, sd=1, number_of_cells=30)
combo_4 = simulate(slope=4, sd=2, number_of_cells=30)
combo_5 = simulate(slope=6, sd=3, number_of_cells=30)
cells_30 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_30$slope_over_var = "2"
cells_30$cells = "cells_30"

slope_var_2 = do.call("rbind", list(cells_7, cells_16, cells_20 ,cells_30))

write.csv(slope_var_2, "ten_reps_fig_2_second_run/slope_var_2")

####slope/var=3####
# 1.5/.5    2.25/.75    3/1   6/2   9/3
combo_1 = simulate(slope=1.5, sd=.5, number_of_cells=7)
combo_2 = simulate(slope=2.25, sd=.75, number_of_cells=7)
combo_3 = simulate(slope=3, sd=1, number_of_cells=7)
combo_4 = simulate(slope=6, sd=2, number_of_cells=7)
combo_5 = simulate(slope=9, sd=3, number_of_cells=7)
cells_7 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_7$slope_over_var = "3"
cells_7$cells = "cells_7"

combo_1 = simulate(slope=1.5, sd=.5, number_of_cells=16)
combo_2 = simulate(slope=2.25, sd=.75, number_of_cells=16)
combo_3 = simulate(slope=3, sd=1, number_of_cells=16)
combo_4 = simulate(slope=6, sd=2, number_of_cells=16)
combo_5 = simulate(slope=9, sd=3, number_of_cells=16)
cells_16 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_16$slope_over_var = "3"
cells_16$cells = "cells_16"

combo_1 = simulate(slope=1.5, sd=.5, number_of_cells=20)
combo_2 = simulate(slope=2.25, sd=.75, number_of_cells=20)
combo_3 = simulate(slope=3, sd=1, number_of_cells=20)
combo_4 = simulate(slope=6, sd=2, number_of_cells=20)
combo_5 = simulate(slope=9, sd=3, number_of_cells=20)
cells_20 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_20$slope_over_var = "3"
cells_20$cells = "cells_20"

combo_1 = simulate(slope=1.5, sd=.5, number_of_cells=30)
combo_2 = simulate(slope=2.25, sd=.75, number_of_cells=30)
combo_3 = simulate(slope=3, sd=1, number_of_cells=30)
combo_4 = simulate(slope=6, sd=2, number_of_cells=30)
combo_5 = simulate(slope=9, sd=3, number_of_cells=30)
cells_30 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_30$slope_over_var = "3"
cells_30$cells = "cells_30"

slope_var_3 = do.call("rbind", list(cells_7, cells_16, cells_20 ,cells_30))

write.csv(slope_var_3, "ten_reps_fig_2_second_run/slope_var_3")

####slop/var=4####
# 2/.5    3/.75   4/1   8/2   12/3
combo_1 = simulate(slope=2, sd=.5, number_of_cells=7)
combo_2 = simulate(slope=3, sd=.75, number_of_cells=7)
combo_3 = simulate(slope=4, sd=1, number_of_cells=7)
combo_4 = simulate(slope=8, sd=2, number_of_cells=7)
combo_5 = simulate(slope=12, sd=3, number_of_cells=7)
cells_7 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_7$slope_over_var = "4"
cells_7$cells = "cells_7"

combo_1 = simulate(slope=2, sd=.5, number_of_cells=16)
combo_2 = simulate(slope=3, sd=.75, number_of_cells=16)
combo_3 = simulate(slope=4, sd=1, number_of_cells=16)
combo_4 = simulate(slope=8, sd=2, number_of_cells=16)
combo_5 = simulate(slope=12, sd=3, number_of_cells=16)
cells_16 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_16$slope_over_var = "4"
cells_16$cells = "cells_16"

combo_1 = simulate(slope=2, sd=.5, number_of_cells=20)
combo_2 = simulate(slope=3, sd=.75, number_of_cells=20)
combo_3 = simulate(slope=4, sd=1, number_of_cells=20)
combo_4 = simulate(slope=8, sd=2, number_of_cells=20)
combo_5 = simulate(slope=12, sd=3, number_of_cells=20)
cells_20 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_20$slope_over_var = "4"
cells_20$cells = "cells_20"

combo_1 = simulate(slope=2, sd=.5, number_of_cells=30)
combo_2 = simulate(slope=3, sd=.75, number_of_cells=30)
combo_3 = simulate(slope=4, sd=1, number_of_cells=30)
combo_4 = simulate(slope=8, sd=2, number_of_cells=30)
combo_5 = simulate(slope=12, sd=3, number_of_cells=30)
cells_30 = do.call("rbind", list(combo_1, combo_2, combo_3 ,combo_4, combo_5))
cells_30$slope_over_var = "4"
cells_30$cells = "cells_30"

slope_var_4 = do.call("rbind", list(cells_7, cells_16, cells_20 ,cells_30))

write.csv(slope_var_4, "ten_reps_fig_2_second_run/slope_var_4")



# If I give you slope it will tell my how much dev there needs to be in order to be accuracte
#Generating data: 


#this should return df with accuracy of all cell numbers in it
simulate(slope=.5, sd=1, number_of_subsamples)
simulate(slope=.75, sd=1.5, number_of_subsamples)



#slope 1  ::    1/1   2/2     3/3   4/4   5/5 #do we want var to be what we're controling for?
#slope 2  ::  2/1   
#slope 3
#slope 4