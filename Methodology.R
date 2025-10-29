library(ggplot2)
library(latex2exp)

err <- seq(0.00001,0.5,by=0.0001)


# Functions for 
precifun <- Vectorize(function(prev,err){
  out <- 0
  if( (err < prev) & (err + prev < 1) ) {
    out <- prev/(prev + err)
  }
  out
},"err")

specififun <- Vectorize(function(prev,err){
  out <- 1
  if( (err < prev) & (err + prev < 1) ) {
    out <- (1-prev-err)/(1 - prev)
  }
  out
},"err")

recallfun <- Vectorize(function(prev,err){
  out <- 0
  if( (err < prev) & (err + prev < 1) ) {
    out <- 1
  }
  out
},"err")

f1fun <- Vectorize(function(prev,err){
  out <- 0
  if( (err < prev) & (err + prev < 1) ) {
    out <- (2*prev)/(2*prev+err)
  }
  out
},"err")

accfun <- Vectorize(function(prev,err){
  if( (err < prev) & (err + prev < 1) ) {
    out <- 1 - err
  }
  if( (err > prev) & (err + prev < 1) ) {
    out <- (1-prev-err)/(1-err)
  }
  out
},"err")

preci1 <- precifun(0.025,err) 
preci2 <- precifun(0.05,err)
preci3 <- precifun(0.1,err)
preci4 <- precifun(0.2,err)

speci1 <- specififun(0.025,err) 
speci2 <- specififun(0.05,err)
speci3 <- specififun(0.1,err)
speci4 <- specififun(0.2,err)

sensi1 <- recallfun(0.025,err) 
sensi2 <- recallfun(0.05,err)
sensi3 <- recallfun(0.1,err)
sensi4 <- recallfun(0.2,err)

fscore1 <- f1fun(0.025,err) 
fscore2 <- f1fun(0.05,err)
fscore3 <- f1fun(0.1,err)
fscore4 <- f1fun(0.2,err)

acc1 <- sapply(err, function(x)accfun(0.025,x)) 
acc2 <- sapply(err, function(x)accfun(0.05,x)) 
acc3 <- sapply(err, function(x)accfun(0.1,x)) 
acc4 <- sapply(err, function(x)accfun(0.2,x)) 


precifun2 <- Vectorize(function(prev,err){
  out <- 0
  if( (err < 1-prev) & (err + prev < 1) ) {
    out <- prev/(prev + err)
  }
  out
},"err")

specififun2 <- Vectorize(function(prev,err){
  out <- 1
  if( (err < 1-prev) & (err + prev < 1) ) {
    out <- (1-prev-err)/(1 - prev)
  }
  out
},"err")

recallfun2 <- Vectorize(function(prev,err){
  out <- 0
  if( (err < 1-prev) & (err + prev < 1) ) {
    out <- 1
  }
  out
},"err")

f1fun2 <- Vectorize(function(prev,err){
  out <- 0
  if( (err < 1-prev) & (err + prev < 1) ) {
    out <- (2*prev)/(2*prev+err)
  }
  out
},"err")

accfun2 <- Vectorize(function(prev,err){
  if( (err < 1-prev) & (err + prev < 1) ) {
    out <- 1 - err
  }
  if( (err > 1-prev) & (err + prev < 1) ) {
    out <- (1-prev-err)/(1-err)
  }
  out
},"err")

preci1b <- precifun2(0.025,err) 
preci2b <- precifun2(0.05,err)
preci3b <- precifun2(0.1,err)
preci4b <- precifun2(0.2,err)
plot(err,preci1b)
points(err,preci2b,col=2)
points(err,preci3b,col=3)
points(err,preci4b,col=4)




speci1b <- specififun2(0.025,err) 
speci2b <- specififun2(0.05,err)
speci3b <- specififun2(0.1,err)
speci4b <- specififun2(0.2,err)
plot(err,speci1b)
points(err,speci2b,col=2)
points(err,speci3b,col=3)
points(err,speci4b,col=4)

sensi1b <- recallfun2(0.025,err) 
sensi2b <- recallfun2(0.05,err)
sensi3b <- recallfun2(0.1,err)
sensi4b <- recallfun2(0.2,err)
plot(err,sensi1b)
points(err,sensi2b,col=2)
points(err,sensi3b,col=3)
points(err,sensi4b,col=4)

fscore1b <- f1fun2(0.025,err) 
fscore2b <- f1fun2(0.05,err)
fscore3b <- f1fun2(0.1,err)
fscore4b <- f1fun2(0.2,err)
plot(err,fscore1b)
points(err,fscore2b,col=2)
points(err,fscore3b,col=3)
points(err,fscore4b,col=4)

acc1b <- sapply(err, function(x)accfun2(0.025,x)) 
acc2b <- sapply(err, function(x)accfun2(0.05,x)) 
acc3b <- sapply(err, function(x)accfun2(0.1,x)) 
acc4b <- sapply(err, function(x)accfun2(0.2,x)) 
plot(err,acc1b,ylim=c(0.0,1))
points(err,acc2b,col=2)
points(err,acc3b,col=3)
points(err,acc4b,col=4)


plot(err,acc1,ylim=c(0.0,1),type="l")
lines(err,acc2,col=2)
lines(err,acc3,col=3)
lines(err,acc4,col=4)

dataplot <- as.data.frame(cbind(err,preci1,preci1b))
ggplot(dataplot, aes(err)) + 
  geom_line(aes(y = preci1)) + 
  geom_line(aes(y = preci1b),linetype="dashed")

ERR <- (rep(err,8))
PREV <- factor((rep(rep(c("2.5%","5%","10%","20%"),each=length(err)),2)),
                   levels=c("2.5%","5%","10%","20%"))
#PREV <- ((rep(rep(c(0.025,0.05,0.1,0.2),each=length(err)),2)))
PRECI <- c(preci1,preci2,preci3,preci4,
           preci1b,preci2b,preci3b,preci4b)
METHOD <- rep(c(1,2),each=length(err)*4)
dataplot <- as.data.frame(cbind(ERR,PRECI,PREV,METHOD))
ggplot(subset(dataplot), aes(x=ERR,y = PRECI,group = PREV,colour=PREV)) + 
  geom_line()+
  facet_wrap(~METHOD,  ncol=2)


ERR <- rep(err,16)
PREV <- rep(rep(c("2.5%","5%","10%","20%"),each=length(err)),4)
yvalues<- c(preci1,preci2,preci3,preci4,
            preci1b,preci2b,preci3b,preci4b,
            fscore1,fscore2,fscore3,fscore4,
            fscore1b,fscore2b,fscore3b,fscore4b,
            sensi1, sensi2, sensi3, sensi4,
            sensi1b, sensi2b, sensi3b, sensi4b,
            speci1, speci2, speci3, speci4,
            speci1b, speci2b, speci3b, speci4b)
METHOD <- rep(rep(c('Standard 0.5-threshold','Prevalence-threshold'),each=length(err)*4),2)
PERF <- rep(c("Precision","F1-score", "Recall", "Specificity"),each=length(err)*8)

dataplot <- as.data.frame(cbind(ERR,yvalues,PREV,METHOD,PERF))

# Factor the necessary columns.
dataplot$PERF <- as.factor(dataplot$PERF)
dataplot$METHOD <- as.factor(dataplot$METHOD)
dataplot$PREV <- as.factor(dataplot$PREV)
dataplot$ERR <- as.numeric(dataplot$ERR)
dataplot$yvalues <- as.numeric(dataplot$yvalues)
dataplot$Prevalence <- factor(dataplot$PREV, levels = c("2.5%","5%","10%","20%"))

# Now to make the plots better as a line plot, we make the point
# of discontinuation NA
dataplot$yvalues[which(diff(yvalues) < -0.1) ] <- NA
dataplot$yvalues[which(diff(yvalues) > 0.001)] <- NA

# Now add some nice PRIME-ROSE colours
PRIME_ROSE_colors <- c(
  "#001F5F", # navy blue
  "#3A0D73", # dark purple
  '#702572', # PRIME-ROSE blue
  '#BD1A8D', # PRIME-ROSE purple
  "#D1554F", # warm red
  "#E17A39", # burnt orange
  "#F4991E", # orange
  '#F37021'  # PRIME-ROSE orange
)
PRIME_ROSE_gradient <- colorRampPalette(PRIME_ROSE_colors)

ggplot(dataplot, aes(x= ERR, y = yvalues, colour = Prevalence)) + 
  geom_line(size = 1.2)+
  facet_grid(PERF ~ METHOD) +
  scale_color_manual(values = PRIME_ROSE_gradient(4)) +
  theme(legend.position = "top") +
  ylab('Performance measure') + xlab(TeX("$\\frac{n_{ 0 1}}{n}$"))

ggsave(filename = 'Figure1.pdf', device = 'pdf', width = 8, height = 10)
