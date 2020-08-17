library(ggplot2)
library(dplyr)

#Figure 1

slope5 = read.csv("~/ALS/data_code_for_fig_1_and_2/fig1/slope0p5_mean_accuraacy")
slope1 = read.csv("~/ALS/data_code_for_fig_1_and_2/fig1/slope1_mean_accuraacy")
slope2 = read.csv("~/ALS/data_code_for_fig_1_and_2/fig1/slope2_mean_accuraacy")
slope4 = read.csv("~/ALS/data_code_for_fig_1_and_2/fig1/slope4_mean_accuraacy")

cells_7 = slope5[slope5[, "cells"] == 7,]
cells_16 = slope5[slope5[, "cells"] == 16,]
cells_20 = slope5[slope5[, "cells"] == 20,]
cells_30 = slope5[slope5[, "cells"] == 30,]
cells_100 = slope5[slope5[, "cells"] == 100,]

cells_7 = slope1[slope1[, "cells"] == 7,]
cells_16 = slope1[slope1[, "cells"] == 16,]
cells_20 = slope1[slope1[, "cells"] == 20,]
cells_30 = slope1[slope1[, "cells"] == 30,]
cells_100 = slope1[slope1[, "cells"] == 100,]
# # 
cells_7 = slope2[slope2[, "cells"] == 7,]
cells_16 = slope2[slope2[, "cells"] == 16,]
cells_20 = slope2[slope2[, "cells"] == 20,]
cells_30 = slope2[slope2[, "cells"] == 30,]
cells_100 = slope2[slope2[, "cells"] == 100,]

cells_7 = slope4[slope4[, "cells"] == 7,]
cells_16 = slope4[slope4[, "cells"] == 16,]
cells_20 = slope4[slope4[, "cells"] == 20,]
cells_30 = slope4[slope4[, "cells"] == 30,]
cells_100 = slope4[slope4[, "cells"] == 100,]


plot(cells_7$deviation, cells_7$mean_accuracy, main="Slope = 4",
     xlab="Variation ", ylab="% Accuracy ", pch=19, type="o",
     ylim=c(0.45,1.0), lty=1) 

arrows(cells_7$deviation, cells_7$mean_accuracy, cells_7$deviation,cells_7$mean_accuracy+cells_7$dev_accuracy, length=0.05, angle=90)
arrows(cells_7$deviation, cells_7$mean_accuracy, cells_7$deviation,cells_7$mean_accuracy-cells_7$dev_accuracy, length=0.05, angle=90)



par(new=T)
plot(cells_16$deviation, cells_16$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0.45,1.0), lty=2)
arrows(cells_16$deviation, cells_16$mean_accuracy, cells_16$deviation,cells_16$mean_accuracy+cells_16$dev_accuracy, length=0.05, angle=90)
arrows(cells_16$deviation, cells_16$mean_accuracy, cells_16$deviation,cells_16$mean_accuracy-cells_16$dev_accuracy, length=0.05, angle=90)

par(new=T)
plot(cells_20$deviation, cells_20$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0.45,1.0), lty=3)
arrows(cells_20$deviation, cells_20$mean_accuracy, cells_20$deviation,cells_20$mean_accuracy+cells_20$dev_accuracy, length=0.05, angle=90)
arrows(cells_20$deviation, cells_20$mean_accuracy, cells_20$deviation,cells_20$mean_accuracy-cells_20$dev_accuracy, length=0.05, angle=90)

par(new=T)
plot(cells_30$deviation, cells_30$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0.45,1.0), lty=4)
arrows(cells_30$deviation, cells_30$mean_accuracy, cells_30$deviation,cells_30$mean_accuracy+cells_30$dev_accuracy, length=0.05, angle=90)
arrows(cells_30$deviation, cells_30$mean_accuracy, cells_30$deviation,cells_30$mean_accuracy-cells_30$dev_accuracy, length=0.05, angle=90)

par(new=T)
plot(cells_100$deviation, cells_100$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0.45,1.0), lty=5)
arrows(cells_100$deviation, cells_100$mean_accuracy, cells_100$deviation,cells_100$mean_accuracy+cells_100$dev_accuracy, length=0.05, angle=90)
arrows(cells_100$deviation, cells_100$mean_accuracy, cells_100$deviation,cells_100$mean_accuracy-cells_100$dev_accuracy, length=0.05, angle=90)


legend("bottomleft", legend = c("7 cells", "16 cells", "20 cells", "30 cells", "100 cells"),
       lty = 1:5, cex = 0.8)



