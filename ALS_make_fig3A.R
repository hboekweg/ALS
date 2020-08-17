library(ggplot2)
library(dplyr)

#FIGURE 3A
cells_7= read.csv("~/ALS/data_code_for_fig_1_and_2/fig2/cells_7")
cells_16 = read.csv("~/ALS/data_code_for_fig_1_and_2/fig2/cells_16")
cells_20= read.csv("~/ALS/data_code_for_fig_1_and_2/fig2/cells_20")
cells_30= read.csv("~/ALS/data_code_for_fig_1_and_2/fig2/cells_30")
cells_100= read.csv("~/ALS/data_code_for_fig_1_and_2/fig2/cells_100")


plot(cells_7$slope_over_var, cells_7$mean_accuracy, main="Slope Over Variation",
     xlab="Slope Over Variation ", ylab="% Accuracy ", pch=19, type="o",
     ylim=c(.5,1.0), lty=1) 

arrows(cells_7$slope_over_var, cells_7$mean_accuracy, cells_7$slope_over_var,cells_7$mean_accuracy+cells_7$deviation, length=0.05, angle=90)
arrows(cells_7$slope_over_var, cells_7$mean_accuracy, cells_7$slope_over_var,cells_7$mean_accuracy-cells_7$deviation, length=0.05, angle=90)



par(new=T)
plot(cells_16$slope_over_var, cells_16$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(.5,1.0), lty=2)
arrows(cells_16$slope_over_var, cells_16$mean_accuracy, cells_16$slope_over_var,cells_16$mean_accuracy+cells_16$deviation, length=0.05, angle=90)
arrows(cells_16$slope_over_var, cells_16$mean_accuracy, cells_16$slope_over_var,cells_16$mean_accuracy-cells_16$deviation, length=0.05, angle=90)


par(new=T)
plot(cells_20$slope_over_var, cells_20$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(.5,1.0), lty=3)
arrows(cells_20$slope_over_var, cells_20$mean_accuracy, cells_20$slope_over_var,cells_20$mean_accuracy+cells_20$deviation, length=0.05, angle=90)
arrows(cells_20$slope_over_var, cells_20$mean_accuracy, cells_20$slope_over_var,cells_20$mean_accuracy-cells_20$deviation, length=0.05, angle=90)

par(new=T)
plot(cells_30$slope_over_var, cells_30$mean_accuracy, pch=19, type="o", ann=F, axes=F, 
     ylim=c(.5,1.0), lty=4)
arrows(cells_30$slope_over_var, cells_30$mean_accuracy, cells_30$slope_over_var,cells_30$mean_accuracy+cells_30$deviation, length=0.05, angle=90)
arrows(cells_30$slope_over_var, cells_30$mean_accuracy, cells_30$slope_over_var,cells_30$mean_accuracy-cells_30$deviation, length=0.05, angle=90)

par(new=T)
plot(cells_100$slope_over_var, cells_100$mean_accuracy, pch=19, type="o", ann=F, axes=F,
     ylim=c(.5,1.0), lty=5)
arrows(cells_100$slope_over_var, cells_100$mean_accuracy, cells_100$slope_over_var,cells_100$mean_accuracy+cells_100$deviation, length=0.05, angle=90)
arrows(cells_100$slope_over_var, cells_100$mean_accuracy, cells_100$slope_over_var,cells_100$mean_accuracy-cells_100$deviation, length=0.05, angle=90)


legend("bottomright", legend = c("7 cells", "16 cells", "20 cells", "30 cells", "100 cells"),
        lty = 1:5, cex = 0.8)


