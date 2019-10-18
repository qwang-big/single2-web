nuc=c('A','C','G','T')

bnorm <- function(d){
x=split(d, substr(d[,1],1,1))
m=t(do.call(cbind,lapply(x,function(v) v[order(v[,1]),2])))
colnames(m)=nuc
melt(m/rowSums(m))
}

calcPWM <- function(d){
x=split(d, nchar(d[,1]))
lapply(x,function(y){
t(do.call(rbind,lapply(1:nchar(y[1,1]),function(i){
m=data.frame(substr(y[,1],i,i),y[,2])
sapply(nuc,function(d) sum(m[m[,1]==d,2])/sum(m[,2]))
})))})
}

logV2 <- function(d){
d[,2]=log(d[,2]+1)
d
}

odd <- function(d) as.numeric(d[seq(1,5,2)])
even <- function(d) as.numeric(d[seq(2,6,2)])
swap <- function(d,i,j){
tmp=d[i]
d[i]=d[j]
d[j]=tmp
d
}

readModel <- function(f){
x=readLines(f)
i=grep('#',x)
df=data.frame(from=i+1, to=c(i[2:length(i)]-1,length(x)))
y=apply(df,1,function(d) x[seq(d[1],d[2])])
names(y)=substr(x[i],2,nchar(x[i]))
lapply(y, function(d){
df2 <- do.call(rbind,strsplit(d,'\t'))
setNames(
    as.data.frame(lapply(1:ncol(df2), function(i) {
      type.convert(df2[,i], as.is = TRUE)}), stringsAsFactors = FALSE), 
  paste0('V', 1:ncol(df2)))
})
}

plotSunburst <- function(f){
x=read.table(f,fill=T)
d=c(odd(x[x[,1]=="Processed",3]), 0, odd(x[x[,1]=="short",3]), as.numeric(x[x[,1]=="written",3]))[-c(5,8)]
d[1]=d[1]/2
d[4]=d[9]=d[1]-d[2]-d[3]
d=swap(d,6,7)
plot_ly(
  labels = c("Total reads", "Forward adapter", "Reverse adapter", "No adapter", "Forward barcode", "No barcode", "Reverse barcode", "No barcode ", ""),
  parents = c("", "Total reads", "Total reads", "Total reads", "Forward adapter", "Forward adapter", "Reverse adapter", "Reverse adapter", "No adapter"),
  values = d,
  type = 'sunburst',
  branchvalues = 'total'
)
}
