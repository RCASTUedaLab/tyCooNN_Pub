args = commandArgs(trailingOnly=TRUE)
trna1 = args[1]
trna2 = args[2]
trna  = args[3]
setwd(trna)
trna1.1 = paste0(trna1,'_merge')
trna2.1 = paste0(trna2,'_merge')
d_nn = read.table(paste0(trna,'_distr.label'))
d_map= read.table(paste0(trna,'_map.lab'))
d_nn[d_nn$V2 == trna1.1,]$V2 = trna1
d_nn[d_nn$V2 == trna2.1,]$V2 = trna2
d = merge(d_nn,d_map,by='V1')
names(d) = c('rid','nn','x','map')
d = d[,c('rid','nn','map')]

d_concordant = d[d$nn == d$map,]
d_concordant = d_concordant[,c('rid','nn')]
write.table(file=paste0(trna,'_concordant.lab'),d_concordant,quote=FALSE,row.names=FALSE,col.names=FALSE,sep=' ')

d_discordant = d[d$nn != d$map,]
write.table(file=paste0(trna,'_discordant.lab'),d_discordant,quote=FALSE,row.names=FALSE,col.names=FALSE,sep=' ')
