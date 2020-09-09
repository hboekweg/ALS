##Figure 4
library(ggplot2)
df = read.csv("~/ALS/jupyter_notebooks/data/figure4_data.csv")


ggplot(df, aes(x = C10_stdev, y = abs_C10.SVEC)) +
  geom_point() +
  xlim(0, 1.5)+
  ylim(0, 1.5)+
  geom_abline(aes(intercept = 0, slope = 1, color="1")) +
  geom_abline(aes(intercept = 0, slope = 2, color="2")) +
  geom_abline(aes(intercept = 0, slope = 4, color="4"))+
  labs(color = "Slope/Varitaion") 

plot.new()

#need to add lbels for each slope line
#set x and y to be the same so it's square


