
MAIN = "Training of tyCooNN (IVT, fMet, v0)"

dat = read.table('history.csv',header=TRUE,sep=',')
dat_aug = read.table('history_arg.csv',header=TRUE,sep=',')
library(repr)

df = rbind(dat,dat_aug)
n = nrow(df)

png('training.png',width=10,height=7,unit='in',res=600)
par(mar=par()$mar + c(0,0,0,12),lwd=2,family='Times')
plot(df$val_accuracy,type='l',col='black',lwd=2,
     xlab='Epoch',ylab='Accuracy',ylim=c(0,1),
     bty='u',cex.lab=1.5,cex.axis=1.3,
     main=MAIN,
     cex.main=1.7,font.main=2,
     font.axis=2,font.lab=2)
lines(df$accuracy,col='blue',lty=2,lwd=2)
lines(0.3 * df$val_loss,type='l',col='red',lwd=2)
lines(0.3 * df$loss,col='orange',lty=2,lwd=2)
axis(side=4,font=2,cex.axis=1.3)
mtext(side=4,line=4,'0.3 * Loss',col='red',cex=1.5,font=2)
abline(v=100,col='magenta',lwd=4)
rect(0,0,100,1,col=rgb(230/256,230/256,230/256,alpha=0.5),border=NA)
rect(100,0,n,1,col=rgb(230/256,0,230/256,alpha=0.1),border=NA)
par(xpd=TRUE)
legend(x=c(20,90),y=c(1.11,1.11),
       legend=c('without augmentation','with augmentation'),horiz=TRUE,bty='n',
       fill=c(rgb(230/256,230/256,230/256,alpha=0.5),rgb(230/256,0,230/256,alpha=0.1)),
       border=NA)
legend(135,1,legend=c('val_accuracy','accuracy','val_loss','loss'),bty='n',
      col=c('black','blue','red','orange'),lty=c(1,2,1,2))
o = dev.off()
