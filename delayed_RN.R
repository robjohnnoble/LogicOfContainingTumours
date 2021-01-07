#delayed treatment tumour size until Nref at Cmax=C
delayed_treat<-function(C,Nref,parms=pars_fig,times=x,logged=FALSE,func=MonroGaffney, max_time = 2e3){
  yini<-init_yini(pars_fig)
  untreated<-ode(y=yini,parms = parms, times = times, func = func)
  #time JUST BEFORE we cross Nref (to simulate reaching exactly Nref before MTD and not changing the time to prog or else)
  cutoff <- dic_find(Nref,untreated[,2]+untreated[,3],0,length(untreated[,2]))
  t <- times[cutoff]
  delay<-NULL
  #Now its just maximum dose from t to max_time
  yini["C"]<-C
  xb <- times[(cutoff + 1):length(times)]
  #Stitching the untreated until N0 to the MTD part
  if(t==0){
    #if t = 0, then it's mtd from start to finish
    mtd<-ode(y=yini,parms=parms,times = xb,func=func)
    delay <- mtd
  }else if(t<max_time){
    #tumour starts where the untreated stopped
    yini["y1"]<-untreated[cutoff,2]
    yini["y2"]<-untreated[cutoff,3]
    yini["C"]<-C
    mtd<-ode(y=yini,parms=parms,times = xb,func=func)
    #then stitching with the untreated part
    delay <- rbind(untreated[1:cutoff,],mtd)
  }else{
    #if t == max_time, then there's no mtd
    delay<-untreated
  }
  return(delay)
}
#plotting delayed treatment. Cseq can be a sequence up to 5 c values
delayed_graph<-function(Cseq,Nref,parms=pars_fig,logged=FALSE,func=MonroGaffney, max_time = 2e3){
  #Cseq = sequence of doses
  #create a vector of color
  col<-c("red","blue","green","magenta","yellow")
  x<-seq(from=0,to=max_time,by=1)
  title<-"not logged"
  log<-""
  if(logged==TRUE){
    title<-"logged"
    log<-"y"
  }
  coul<-NULL
  legend<-NULL
  lty<-NULL
  for(i in 1:length(Cseq)){
    u <-delayed_treat(Cseq[i],Nref,times=x,parms=parms,logged=logged,func=func, max_time = max_time)
    if(i==1){
      plot(x=x,(u[,2]+u[,3]),log=log,type="l",col=col[i],main=paste(title,sprintf("delayed treatment at %s",Nref)),
           xlim=c(0,2000),ylim=c(parms["N0"]/2, 2e12))
    }else{
      lines(x=x,(u[,2]+u[,3]),type="l",col=col[i])
    }
    #color
    coul<-c(coul,col[i])
    legend<-c(legend,sprintf("MTD;C=%s",Cseq[i]))
    #type of line
    lty<-c(lty,"solid")
  }
  abline(h=(Nref),lty="dashed")
  legend("bottomright",lty=lty, col=coul,legend = legend,lwd=2);
  
}
#Time for the delayed treatment at Ndelay to reach Nref
time_delayed<-function(Ndelay,Nref,parms=pars_fig){
  #untreated until Ndelay
  t1<-time_from_to(parms["N0"],Ndelay)
  #calculate r1 = resistant population at Ndelay
  r1<-parms["r0"]*Ndelay/parms["N0"]
  #growth from r1 to Nref
  t2<-t1+time_from_to(r1,Nref)
  #if we delay too late, it's the same as an untreated tumor
  t2[which(Ndelay>=Nref)]<-time_from_to(parms["N0"],Nref)
  t2
}
#time to prog/fail/surv for various ndelay
time_delayed_plot<-function(parms=pars_fig,logged=FALSE){
  ndelay <-seq(from=parms["N0"],to=parms["Ncrit"]*1.2,length=10000)
  prog<-time_delayed(ndelay,parms["N0"],parms=parms)
  fail<-time_delayed(ndelay,parms["Nacc"],parms=parms)
  surv<-time_delayed(ndelay,parms["Ncrit"],parms=parms)
  title<-"not logged"
  log<-""
  if(logged==TRUE){
    title<-"logged"
    log<-"x"
  }
  plot(x=(ndelay),fail,log=log,type="l",col="red",ylim=c(0,250),main=paste(title,"delayed treatment"))
  lines(x=(ndelay),surv,type="l",col="blue")
  lines(x=(ndelay),prog,type="l",col ="green")
  legend("topleft",lty=c("solid","solid","solid"), col=c("blue","red","green"),legend = c("Tsurv","Tfail","Tprog"),lwd=2);
  
  abline(v=((parms["N0"])),lty=2)
  abline(v=((parms["Ncrit"])),lty=2)
  abline(v=((parms["Nacc"])),lty=2)
}