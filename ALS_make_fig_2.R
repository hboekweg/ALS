library(ggplot2)
library(babynames) # provide the dataset: a dataframe called babynames
library(dplyr)

slope5 = read.csv("~/ALS/jupyter_notebooks/data/fig_1/slope5.csv")
slope1 = read.csv("~/ALS/jupyter_notebooks/data/fig_1/slope1.csv")
slope2 = read.csv("~/ALS/jupyter_notebooks/data/fig_1/slope2.csv")
slope4 = read.csv("~/ALS/jupyter_notebooks/data/fig_1/slope4.csv")


cells_7 = slope4[slope4[, "cell"] == 7,]
cells_16 = slope4[slope4[, "cell"] == 16,]
cells_20 = slope4[slope4[, "cell"] == 20,]
cells_30 = slope4[slope4[, "cell"] == 30,]
cells_100 = slope4[slope4[, "cell"] == 100,]


plot(cells_7$deviation, cells_7$mean, main="Slope = 4",
     xlab="Variation ", ylab="% Accuracy ", pch=19, type="o",
     ylim=c(0,1.0), lty=1)

arrows(cells_7$deviation, cells_7$mean, cells_16$deviation,cells_7$mean+cells_7$standard.deviation, length=0.05, angle=90)
arrows(cells_7$deviation, cells_7$mean, cells_16$deviation,cells_7$mean-cells_7$standard.deviation, length=0.05, angle=90)

par(new=T)
plot(cells_16$deviation, cells_16$mean, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0,1.0), lty=2)
arrows(cells_16$deviation, cells_16$mean, cells_16$deviation,cells_16$mean+cells_16$standard.deviation, length=0.05, angle=90)
arrows(cells_16$deviation, cells_16$mean, cells_16$deviation,cells_16$mean-cells_16$standard.deviation, length=0.05, angle=90)

par(new=T)
plot(cells_20$deviation, cells_20$mean, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0,1.0), lty=3)
arrows(cells_20$deviation, cells_20$mean, cells_20$deviation,cells_20$mean+cells_20$standard.deviation, length=0.05, angle=90)
arrows(cells_20$deviation, cells_20$mean, cells_20$deviation,cells_20$mean-cells_20$standard.deviation, length=0.05, angle=90)

par(new=T)
plot(cells_30$deviation, cells_30$mean, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0,1.0), lty=4)
arrows(cells_30$deviation, cells_30$mean, cells_30$deviation,cells_30$mean+cells_30$standard.deviation, length=0.05, angle=90)
arrows(cells_30$deviation, cells_30$mean, cells_30$deviation,cells_30$mean-cells_30$standard.deviation, length=0.05, angle=90)

par(new=T)
plot(cells_100$deviation, cells_100$mean, pch=19, type="o", ann=F, axes=F, 
     ylim=c(0,1.0), lty=5)
arrows(cells_100$deviation, cells_100$mean, cells_100$deviation,cells_100$mean+cells_100$standard.deviation, length=0.05, angle=90)
arrows(cells_100$deviation, cells_100$mean, cells_100$deviation,cells_100$mean-cells_100$standard.deviation, length=0.05, angle=90)


legend("bottomleft", legend = c("7 cells", "16 cells", "20 cells", "30 cells", "100 cells"),
        lty = 1:5, cex = 0.8)
