library(ggplot2)
data <- read.csv('~/Desktop/data/Pi/plot/IGS_sort_as_cp_order.txt', 
                 sep = '\t', header = FALSE) 
data$V1 <- factor(data$V1, levels = data$V1)

p1 <- ggplot(data, aes(x = V1, y = V2, group=1)) + geom_line(size = 0.3) + 
  ylim(0, 0.4) + theme_bw() + 
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle = 75, vjust = 0.9, hjust=1, size = 7, 
                                   face="bold.italic", margin = margin(-0.2,0,0,0, 'cm'), 
                                   colour = 'black'),  
        plot.margin = margin(1,1,1,1,'cm'),
        text=element_text(face = "bold"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  ylab("Nucleotide diversity")  + 
  geom_hline(aes(yintercept = 0.25), colour = "#990000", linetype = "dashed")  +
  geom_point(size = 1, shape = 21, fill = 'blue')

data$V3 <- ifelse(data$V2 > 0.25, data$V2, NA)

p1 <- ggplot(data, aes(x = V1, y = V2, group=1)) + geom_line(size = 0.3) + 
  ylim(0, 0.4) + theme_bw() + 
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle = 75, vjust = 0.9, hjust=1, size = 7, 
                                   face="bold.italic", margin = margin(-0.2,0,0,0, 'cm'), 
                                   colour = 'black'),  
        plot.margin = margin(1,1,1,1,'cm'),
        text=element_text(face = "bold"), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  ylab("Nucleotide diversity")  + 
  geom_hline(aes(yintercept = 0.25), colour = "#990000", linetype = "dashed")  +
  geom_point(size = 1, shape = 21, fill = 'blue')+
  geom_text(size = 2.5,  aes(label = V3), vjust=-0.35)
p1

# 拼图
library(cowplot)
combined_plot <- plot_grid(p1, p2, nrow = 2, align = "v", labels = c("A", "B"))
combined_plot
