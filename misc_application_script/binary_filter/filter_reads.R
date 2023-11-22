
set.seed(10)

args = commandArgs(trailingOnly=TRUE)
trna = args[1]
trna1 = args[2]

inp_file = paste0(trna,"/",trna,".label")

dat0 = read.table(inp_file)
dat = dat0[dat0$V2 != -1,]

search_species = paste0(trna1,"_ivt")
index_1 = which(dat$V5 != trna1)

dat_fixed  = data.frame('read_id'=dat[index_1,]$V1,'name'=trna,'type'='default')

out_file = paste0(trna,"/",trna,"_filter.label")
write.table(file=out_file,dat_fixed,quote=FALSE,col.names=FALSE,row.names=FALSE,sep=' ')
