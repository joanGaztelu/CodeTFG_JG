##Code by Joan Gaztelu Camps
###########################################################

#load necessary libraries
library(ggsignif);library(ggplot2);library(plyr);

#Print subject residual plot
#sound options: "Flat" or "Loom"
analysis_plots<-function(sound){
  setwd("directory to store the plots")
  
  #unique times in all dataset
  times<-c(unique((ST$Time)))
  
  
  #adding factors information to pxqnts
  #fit$xqnts[9:11]<-fit$vars[9:11]
  fit$xqnts[10:12]<-fit$vars[9:11]
  fit$pxqnts[10:12]<-fit$xqnts[10:12]
  
  
  #changing column names of pqxnt data frame
  names(fit$pxqnts)[1]<-"10."
  names(fit$pxqnts)[2]<-"20."
  names(fit$pxqnts)[3]<-"30."
  names(fit$pxqnts)[4]<-"40."
  names(fit$pxqnts)[5]<-"50."
  names(fit$pxqnts)[6]<-"60."
  names(fit$pxqnts)[7]<-"70."
  names(fit$pxqnts)[8]<-"80."
  names(fit$pxqnts)[9]<-"90."
  
  
  i<-1
  n<-1
  r<-8
  #CO
  #for(sound in sounds){
  
  #selecting specific set of data for variable sound
  
  var4<-fit$vars[fit$vars$ST.Sound== sound,]
  xqnts<-fit$xqnts[fit$xqnts$ST.Sound==sound,][1:9]
  pxqnts<-fit$pxqnts[fit$pxqnts$ST.Sound==sound,][1:9]
  
  #joining both observed and predicted quantile dataframes  
  quantiles<- cbind(xqnts, pxqnts)
  #checking dimensions   
  ndat <- dim(var4)[1]
  #computing standard deviations for postierior residual standarization
  xse <- pxse <- sapply(1:ndat, function(ind)  sqrt(var4$alpha[ind]/(var4$gamma[ind]^3)))
  #computing residuals
  resid <- abs((xqnts/xse)-(pxqnts/pxse))
  #changing to long format
  residuals_long <- gather(resid[1:9], Quantile, Residual, `0.1`:`0.9`, factor_key=TRUE)
    
  cresid <- apply(abs(residuals_long[2]),1,sum) #vectorizing
  cresid <- data.frame(cresid, residuals_long[,1])
  cresid2 <- apply(abs(resid),1,sum)
  
  da0<-c(1:length(cresid[,1]))#creating range
    
  lencr<- length(cresid[,1])
  new_f2<- data.frame(cresid, da0 )
  names(new_f2)[names(new_f2) == "residuals_long...1."] <- "Quantiles"
    
  interval<- unique(new_f2$Quantiles)
  #dataframe to store deltas
  panel_fig2<-data.frame(delta = double(), quantile = character())
  for(k in interval){
    x<-new_f2[new_f2$Quantiles == k,][1]
    panel_fig2<- panel_fig2 %>% add_row(delta = round(mean(x$cresid),2), quantile = k)
  }
  
  plot_deltas<-ggplot()+
    geom_col(data=panel_fig2, aes(x=(quantile), y=(delta), group=quantile, fill=quantile))+theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey50"), axis.text.x = element_blank(), axis.ticks.x = element_blank())+ scale_fill_viridis_d()+labs( y="∆", x="Quantiles")+
    ggtitle(paste("∆ for sound ", sound)) + expand_limits(y=c(0, 0.3))
  #ggsave(filename=paste0("deltas ", sound, ".png"), plot = plot_deltas)
  print(plot_deltas)
}


#Plot mean plots
mean_plots("Flat")
mean_plots<-function(sound){
  #setwd("set working directory")
  setwd("/Volumes/PANDORA/Courses/PR/PLOTS/Paleta1/means")
  times<-c(unique((ST$Time)))
  fit$pxqnts[10:12]<-fit$xqnts[10:12]
  sequence<-c(seq(0.1, 0.9, 0.1))
  new_means_f<-data.frame(sequence)
  new_sd<- data.frame(sequence)
  for(time in times){
    sd <- c()
    vec<-c()
    xqnts<-fit$xqnts[fit$xqnts$ST.Time==time & fit$xqnts$ST.Sound==sound,][1:9]
    pxqnts<-fit$pxqnts[fit$pxqnts$ST.Time==time & fit$pxqnts$ST.Sound==sound,][1:9]
    res <- abs(xqnts-pxqnts)
    SD <- c(apply(res,2,sd))
    for(k in 1:9){
      vec[k]<-round(mean(xqnts[,k]))
      
    }
    new_means_f<- cbind(new_means_f, vec)
    names(new_means_f)[names(new_means_f) == "vec"]<-time
    new_sd[,paste0(time, "_SD")]<- SD
  }
  
  lev<-c(-700,300,800, 1500, 2200, 2700)
  means_long<-gather(new_means_f, Time, Mean, `2200`:`1500`)
  sd_long<-gather(new_sd, Time, SD, `2200_SD`:`1500_SD`)
  table_M_S<-cbind(means_long, sd_long[,3])
  names(table_M_S)[names(table_M_S) == "sd_long[, 3]"] <- "SD"
  
  if(sound=="Loom"){
    plot_loom<-ggplot(table_M_S) + geom_line(aes(x = factor(Time, level=lev), y=Mean, group= sequence))+
      geom_errorbar(aes(x=factor(Time, level=lev), ymin=Mean-SD, ymax=Mean+SD, group=(sequence), colour=factor(sequence)), width=0.2) +
      geom_point(aes(x = factor(Time, level=lev), y=Mean, group= sequence, color=factor(sequence))) + 
      theme_bw() + theme(axis.text.x = element_text(angle = 90))+labs(color="Quantiles" ,x="Temporal delay [ms]", y = "Means  [ms]", title ="Mean Quantiles", subtitle = "Sound Loom" ) + scale_x_discrete(labels = c("No sound", "300", "800", "1500", "2200", "1700"))+ scale_color_viridis_d()
    print(plot_loom) 
    ggsave(filename = paste0("Mean", sound, ".png",sep=""), plot = plot_loom, units = "cm")
  }
  else if(sound=="Flat"){
    plot_flat<-ggplot(table_M_S) + geom_line(aes(x = factor(Time, level=lev), y=Mean, group= sequence))+geom_errorbar(aes(x=factor(Time, level=lev), ymin=Mean-SD, ymax=Mean+SD, group=(sequence), colour=factor(sequence)), width=0.2)+ geom_point(aes(x = factor(Time, level=lev), y=Mean, group= sequence, color=factor(sequence))) + 
      theme_bw()+theme(axis.text.x = element_text(angle = 90)) + labs(color="Quantiles" ,x="Temporal Delay [ms]", y = "Means [ms]", title ="Mean Quantiles", subtitle = "Sound Flat")+scale_x_discrete(labels = c("No sound", "300", "800", "1500", "2200", "1700")) + scale_color_viridis_d() 
    print(plot_flat)
    ggsave(filename = paste0("Mean", sound, ".png",sep=""), plot = plot_flat, units = "cm")
  }
}

#Plots for the parameters
par_plot<- function(sound){
  setwd("/Volumes/PANDORA/Courses/PR/PLOTS/Paleta1/parameters")
  
  times<-c("300","800","1500","2200","2700")
  soundf<-c("Flat","Flat","Flat","Flat","Flat")
  soundl<-c("Loom","Loom","Loom","Loom","Loom")
  
  Mgamma_f<-c();Malpha_f<-c();Mtheta_f<-c()
  
  Mgamma_l<-c();Malpha_l<-c();Mtheta_l<-c()
  
  MgammaF_var<-c();MalphaF_var<-c();MthetaF_var<-c()
  
  MgammaL_var<-c();MalphaL_var<-c();MthetaL_var<-c()
  
  MgammaF_sd<-c();MalphaF_sd<-c();MthetaF_sd<-c()
  
  MgammaL_sd<-c();MalphaL_sd<-c();MthetaL_sd<-c()
  
  
  i<-1
  for(time in times){
    n <- length(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound==sound,][,1])
    Mgamma_f[i]<-mean(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,1])
    Malpha_f[i]<-mean(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,2])
    Mtheta_f[i]<-mean(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,3])
    
    Mgamma_l[i]<-mean(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,1])
    Malpha_l[i]<-mean(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,2])
    Mtheta_l[i]<-mean(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,3])
    
    MgammaF_var<-var(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,1])
    MalphaF_var<-var(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,2])
    MthetaF_var<-var(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,3])
    
    MgammaL_var<-var(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,1])
    MalphaL_var<-var(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,2])
    MthetaL_var<-var(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,3])
    
    
    MgammaF_sd[i]<-sd(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,1])/sqrt(n)
    MalphaF_sd[i]<-sd(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,2])/sqrt(n)
    MthetaF_sd[i]<-sd(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Flat",][,3])/sqrt(n)
    
    MgammaL_sd[i]<-sd(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,1])/sqrt(n)
    MalphaL_sd[i]<-sd(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,2])/sqrt(n)
    MthetaL_sd[i]<-sd(fit$vars[fit$vars$ST.Time==time & fit$vars$ST.Sound=="Loom",][,3])/sqrt(n)
    
    i<-i+1
  }
  
  mpF_df<-data.frame(times,soundf,Mgamma_f, Malpha_f, Mtheta_f,MgammaF_sd,MalphaF_sd,MthetaF_sd, MgammaF_var, MalphaF_var, MthetaF_var)
  mpL_df<-data.frame(times,soundl,Mgamma_l, Malpha_l, Mtheta_l,MgammaL_sd,MalphaL_sd,MthetaL_sd, MgammaL_var, MalphaL_var, MthetaL_var)
  joined<-cbind(mpF_df, mpL_df[2:11])
  lev<-c(300,800, 1500, 2200, 2700)
  
  gmplot<-ggplot(data = joined) +
    geom_errorbar(aes(x=factor(times, level = lev), ymin=Mgamma_f-MgammaF_sd,ymax=Mgamma_f+MgammaF_sd, colour=factor(soundf))) +
    geom_point(aes(x=factor(times, level = lev),Mgamma_f, group =MgammaF_var, colour=factor(soundf)))+
    geom_line(aes(x=factor(times, level = lev),Mgamma_f, group =MgammaF_var, colour=factor(soundf)))+ 
    geom_errorbar(aes(x=factor(times, level = lev), ymin=Mgamma_l-MgammaL_sd,ymax=Mgamma_l+MgammaL_sd, colour=factor(soundl))) +
    geom_point(aes(x=factor(times, level = lev),Mgamma_l, group = MgammaL_var,colour=factor(soundl)))+
    geom_line(aes(x=factor(times, level = lev),Mgamma_l, group =MgammaL_var, colour=factor(soundl)))+
    scale_color_manual(values = c("blue","darkgreen"))+
    labs(x="Onset Times [ms]", y="Means", title = "Gamma", colour = "Sound") + theme_bw()
  ggsave(filename = paste0("_Mean_gamma.png"), plot = gmplot, units = "cm")
  
  
  alplot<-ggplot(data = joined) +
    geom_signif(data =fit$vars ,aes(x=factor(fit$vars$ST.Time, levels = lev),y=fit$vars$alpha, y_position = 26, tip_length= 0.005, annotations = "*"),comparisons = list(c("800", "1500")), map_signif_level=TRUE, test = "t.test", manual = TRUE) +
    geom_signif(data =fit$vars ,aes(x=factor(fit$vars$ST.Time, levels = lev),y=fit$vars$alpha, y_position = 25, tip_length= 0.005, annotations = "*"),comparisons = list(c("300", "800")), map_signif_level=TRUE, test = "t.test", manual = TRUE) +
    geom_errorbar(aes(x=factor(times, level = lev), ymin=Malpha_f-MalphaF_sd,ymax=Malpha_f+MalphaF_sd, colour=factor(soundf))) +
    geom_point(aes(x=factor(times, level = lev),Malpha_f, group =MalphaF_var, colour=factor(soundf),show.legend = FALSE))+
    geom_line(aes(x=factor(times, level = lev),Malpha_f, group =MalphaF_var,colour=factor(soundf) ))+ 
    geom_errorbar(aes(x=factor(times, level = lev), ymin=Malpha_l-MalphaL_sd,ymax=Malpha_l+MalphaL_sd, colour=factor(soundl))) +
    geom_point(aes(x=factor(times, level = lev),Malpha_l, group = MalphaL_var, colour=factor(soundl)))+
    geom_line(aes(x=factor(times, level = lev),Malpha_l, group =MalphaL_var, colour=factor(soundl)))+
    scale_color_manual(values = c("blue","darkgreen", "blue", "darkgreen"))+
    labs(x="Onset Times [ms]", y="Means", title = "Alpha", colour = "Sound") + theme_bw()
  
  ggsave(filename = paste0("_Mean_Alpha.png"), plot = alplot, units = "cm")
  
  thplot<- ggplot() +  geom_signif(data =fit$vars ,aes(x=factor(fit$vars$ST.Time, levels = lev),y=fit$vars$theta, y_position = 270, tip_length= 0.005, annotations = "*"),comparisons = list(c("300", "2200")), map_signif_level=TRUE, test = "t.test", manual = TRUE) +
    geom_signif(data =fit$vars ,aes(x=factor(fit$vars$ST.Time, levels = lev),y=fit$vars$theta, y_position = 265, tip_length= 0.005, annotations = "*"),comparisons = list(c("300", "1500")), map_signif_level=TRUE, test = "t.test", manual = TRUE) +
    geom_signif(data =fit$vars ,aes(x=factor(fit$vars$ST.Time, levels = lev),y=fit$vars$theta, y_position = 275, tip_length= 0.005, annotations = "*"),comparisons = list(c("300", "2700")), map_signif_level=TRUE, test = "t.test", manual = TRUE) +
    geom_signif(data =fit$vars ,aes(x=factor(fit$vars$ST.Time, levels = lev),y=fit$vars$theta, y_position = 260, tip_length= 0.005, annotations = "*"),comparisons = list(c("800", "2700")), map_signif_level=TRUE, test = "t.test", manual = TRUE) +
    geom_errorbar(data = joined, aes(x=factor(times, level = lev), ymin=Mtheta_f-MthetaF_sd,ymax=Mtheta_f+MthetaF_sd, colour=factor(soundf)))+
    geom_point(data = joined,aes(x=factor(times, level = lev),Mtheta_f, group =MthetaF_var, colour=factor(soundf)))+
    geom_line(data = joined,aes(x=factor(times, level = lev),Mtheta_f, group =MthetaF_var, colour=factor(soundf)))+
    geom_errorbar(aes(x=factor(times, level = lev), ymin=Mtheta_l-MthetaL_sd,ymax=Mtheta_l+MthetaL_sd, colour=factor(soundl))) +
    geom_point(data = joined,aes(x=factor(times, level = lev),Mtheta_l, group = MthetaL_var, colour=factor(soundl)))+
    geom_line(data = joined,aes(x=factor(times, level = lev),Mtheta_l, group =MthetaL_var, colour=factor(soundl)))+
    scale_color_manual(values = c("blue","darkgreen"))+
    labs(x="Onset Times", y="Means", title = "Theta", colour = "Sound") + theme_bw()
  
  ggsave(filename = paste0("_Mean_Theta.png"), plot = thplot, units = "cm")
}