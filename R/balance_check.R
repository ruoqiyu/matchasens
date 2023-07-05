balance_check<-function(fdata,mdata,fz,mz,measure='smd'){
  stopifnot(dim(fdata)[2]==dim(mdata)[2])
  stopifnot(colnames(fdata)==colnames(mdata))
  stopifnot(all(fz%in%c(0,1)))
  stopifnot(all(mz%in%c(0,1)))

  if (is.vector(fdata)) fdata<-matrix(fdata,ncol=1)
  if (is.vector(mdata)) mdata<-matrix(mdata,ncol=1)
  fd0<-fdata[which(fz==0),]
  fd1<-fdata[which(fz==1),]
  md0<-mdata[which(mz==0),]
  md1<-mdata[which(mz==1),]
  if (is.vector(fd0)) fd0<-matrix(fd0,length(fd0),1)
  if (is.vector(md0)) md0<-matrix(md0,length(md0),1)
  if (is.vector(fd1)) fd1<-matrix(fd1,length(fd1),1)
  if (is.vector(md1)) md1<-matrix(md1,length(md1),1)

  if (measure=='smd'){
    varf0<-apply(fd0,2,var,na.rm=TRUE)
    varf1<-apply(fd1,2,var,na.rm=TRUE)
    meanf0<-apply(fd0,2,mean,na.rm=TRUE)
    meanf1<-apply(fd1,2,mean,na.rm=TRUE)
    meanm0<-apply(md0,2,mean,na.rm=TRUE)
    meanm1<-apply(md1,2,mean,na.rm=TRUE)
    smdf<-(meanf1-meanf0)/sqrt((varf0+varf1)/2)
    smdm<-(meanm1-meanm0)/sqrt((varf0+varf1)/2)
    
    r<-cbind(meanf1,meanm0,meanf0,smdm,smdf)
    rownames(r)<-colnames(fdata)
    colnames(r)<-c('Treated Mean','Control Match Mean','Control All Mean',
                   'Control Match SMD','Control All SMD')
  }else if (measure=='t.test'){
    nvar=ncol(fdata)
    treatmat.before=fdata[fz==1,]
    controlmat.before=fdata[fz==0,]
    treatmat.after=mdata[mz==1,]
    controlmat.after=mdata[mz==0,]
    t.test.pval.before=rep(0,nvar)
    t.test.pval.after=rep(0,nvar)
    for (i in 1:nvar){
      t.test.pval.before[i]=
        stats::t.test(treatmat.before[,i],controlmat.before[,i])$p.value
      t.test.pval.after[i]=
        stats::t.test(treatmat.after[,i],controlmat.after[,i])$p.value
    }
    r=cbind(t.test.pval.before,t.test.pval.after)
    rownames(r)<-colnames(fdata)
  }else{
    warning('Consider choosing a different covariate balance measure.')
  }
  r
}
