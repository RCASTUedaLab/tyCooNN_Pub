dat = read.table('history.csv',sep=',',header=T)
dat2 = read.table('history_arg.csv',sep=',',header=T)
dat$aug = F
dat2$aug = T
d = rbind(dat,dat2)
library(ggplot2)

loss = (d$loss - min(d$loss)) / (max(d$loss) - min(d$loss))
val_loss = (d$val_loss - min(d$val_loss)) / (max(d$val_loss) - min(d$val_loss))
d1 = data.frame(val=d$accuracy,label='Accuracy,Train',x=1:110)
d2 = data.frame(val=d$val_accuracy,label='Accuracy,Test',x=1:110)
d3 = data.frame(val=loss,label='Loss,Train',x=1:110)
d4 = data.frame(val=val_loss,label='Loss,Test',x=1:110)
D = rbind(d1,d2,d3,d4)
ggplot(data=D,aes(y=val,x=x,fill=label,color=label)) + geom_line(size=1.5) + 
  xlab('Epoch') + ylab('Measures') + 
  geom_vline(xintercept = 100, linetype="dotted", color='black',size=1.5) + 
  theme_classic(base_size=22,base_family='Times')
ggsave('rcc_train.png',dpi=300,width = 11.4, height = 6.69,units = "in", device='png')
