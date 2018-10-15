#title: "plot_result"
#author: "Zhiyuan_Jim_Li"
#date: "September 17, 2018"

# 1. read in data and build data frame
rm(list = ls()) # clear workspace
cat("\014") #clear console
sims_node_power = read.csv(file="/home/zhiyuan/Desktop/go_sim/result/sims_node_specific_power.csv", header=FALSE, sep=",") # change it to your path!
hyper_node_power = read.csv(file="/home/zhiyuan/Desktop/go_sim/result/hypergeometric_node_specific_power.csv", header=FALSE, sep=",") # change it to your path!
sims_FDP_power = read.csv(file="/home/zhiyuan/Desktop/go_sim/result/sims_test_result.csv", header=FALSE, sep=",") # change it to your path!
hyper_FDP_power = read.csv(file="/home/zhiyuan/Desktop/go_sim/result/hypergeometric_test_result.csv", header=FALSE, sep=",") # change it to your path!

hyper_FDP_power = data.frame(hyper_FDP_power)
colnames(hyper_FDP_power) <- c("regime", "experiment","FDP","power")
sims_FDP_power = data.frame(sims_FDP_power)
colnames(sims_FDP_power) <- c("regime", "experiment","FDP","power")

# 2. Plot FDP and Power of two tests
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
}
library(ggplot2)

p1 <- ggplot(hyper_FDP_power, aes(x=regime, y=FDP)) +
  geom_boxplot(aes(group = cut_width(regime, 10)),outlier.colour = "blue", outlier.shape = 1)+
  labs(title="FDP of hypergeometric test",x="n_regimes", y = "False Discovery Proportion")+
  coord_cartesian(xlim = c(10, 120),ylim=c(0,1))+ 
  scale_x_continuous(breaks = seq(10, 120, by = 10)) + 
  geom_hline(yintercept = 0.05,col="red",linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5))
p1 # plot FDP of hypergeometric test

p2 <- ggplot(sims_FDP_power, aes(x=regime, y=FDP)) +
  geom_boxplot(aes(group = cut_width(regime, 10)),outlier.colour = "blue", outlier.shape = 1)+
  labs(title="FDP of sims test",x="n_regimes", y = "False Discovery Proportion")+
  coord_cartesian(xlim = c(10, 120),ylim=c(0,1))+ 
  scale_x_continuous(breaks = seq(10, 120, by = 10)) + 
  geom_hline(yintercept = 0.05,col="red",linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5))
p2 # plot FDP of sims test

p3 <- ggplot(hyper_FDP_power, aes(x=regime, y=power)) +
  geom_boxplot(aes(group = cut_width(regime, 10)),outlier.colour = "blue", outlier.shape = 1)+
  labs(title="Power of hypergeometric test",x="n_regimes", y = "False Discovery Proportion")+
  coord_cartesian(xlim = c(10, 120),ylim=c(0,1))+ 
  scale_x_continuous(breaks = seq(10, 120, by = 10)) +
  theme(plot.title = element_text(hjust = 0.5))
p3 # plot power of hypergeometric test

p4 <- ggplot(sims_FDP_power, aes(x=regime, y=power)) +
  geom_boxplot(aes(group = cut_width(regime, 10)),outlier.colour = "blue", outlier.shape = 1)+
  labs(title="Power of sims test",x="n_regimes", y = "False Discovery Proportion")+
  coord_cartesian(xlim = c(10, 120),ylim=c(0,1))+ 
  scale_x_continuous(breaks = seq(10, 120, by = 10)) +
  theme(plot.title = element_text(hjust = 0.5))
p4 # plot power of sims test

# 3. Plot node-specific power of two tests
competitive_non_null = c(3,7,8,19,20,24,28,40,43,48,69,96)+1 # 
self_contained_non_null = c(0,1,3,4,7,8,13,19,20,24,28,40,43,48,69,96)+1
competitive=hyper_node_power[competitive_non_null,]
contained=sims_node_power[self_contained_non_null,]
Col = sort(unique(hyper_FDP_power[,1]))  # get the unique regimes from smallest to largest
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
}
library("gplots")

p5 <- heatmap.2(as.matrix(competitive), labCol= Col, tracecol=NA ,Rowv= FALSE, Colv = FALSE, dendrogram = "none", col = cm.colors, margins=c(5,5),xlab = "control samples", ylab =  "non-null gene",main="Node-specific power for\n hypergeometric test", key.title ="colorbar")

p6 <-  heatmap.2(as.matrix(competitive), labCol= Col, tracecol=NA ,Rowv= FALSE, Colv = FALSE, dendrogram = "none", col = cm.colors, margins=c(5,5),xlab = "control samples", ylab =  "non-null gene",main="Node-specific power for\n sims test", key.title ="colorbar")

