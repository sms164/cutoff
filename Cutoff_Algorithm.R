#Master Cutoff Algorithm Functions 

#The cutoff algorithm contained in this R script allows the user to descriminate between a positive and negative sample on a continuous diagnostic without a good gold standard via mixture modeling.
#The algorithm runs 10 different potential types of mixture models which are based on the normal and skew-normal distributions, looking at 1-5 components each. The optimal distribution and component combination is chosen by BIC (which should be minimized)
#Note: the fitloops function takes some time (typically 2-3 minutes) to run and can fail, if it does, read the error messages and decide whether to try and run it again, or restrict the maximum number of components tested.
#The user then checks which distribution has been chosen and can manually chose a different distribution if desired (there should be good scientific or clinical reasoning to choose a different combination)
#If a model with more than two components is optimal by BIC, then the user must allot components to the positive and negative side using their knowlege of the epidemiological conditions where the data comes from.
#After the decision of where to cut is made, the user can then output graphs and tables that detail what proportion of the population falls into what category, with various degrees of certainty.


#To complete this process, there are 5 necessary steps
#I)   Install/load packages and load in functions (run lines 67 - 1270)
#II)  Load in data (run lines: 1288)
#III) Run each line of code separately, from 1314 to the end of the program, following the directions in the comments
#     The code will: 1)determine which type of distribution and number of sub-populations is optimal; 2) let the user  
#     determine which cutpoint is best (if a model with more than two components is optimal) and 3) display the cutoff results
 
 
#-----------------------------------------------
#Table of Contents
#-----------------------------------------------
#[I] Install/load packages and load in functions
#[IA] Install Packages
#[IB] Load Packages
#[IC] Load Functions
#[IC1]  fitloops
#[IC2]  bicgraph
#[IC3]  modelpick
#[IC4]  rawuncertgraph
#[IC5]  rawdistgraph
#[IC6]  rawhistcuts
#[IC7]  cutoff
#[IC8] cutuncertgraph
#[IC9] cutdistgraph
#[IC10] summaryout
#-----------------------------------------------
#[II] Load and clean data
#[IIA] Import data 
#[IIB] Create dataframe with two variables, id (must be unique identifier) and data
#-----------------------------------------------
#[III] Determine which type and number of distributions is optimal
#[IIIA] Check data for distributional assumptions and extreme outliers
#[IIIB] Run fitting function and determine which distribution and number of components is best by BIC
#[IIIC] Pick the best combination of distribution and number of components (usually the one that is optimal by BIC unless there is scientific rationale)
#-----------------------------------------------
#[IV]  Determine which cutpoint is best if more than two components are optimal
#[IVA] Investigate proportions of positivity and compare to scientific context
#[IVB] Utilize visualizations of uncertainty, distributions and cutpoint to aid in decision-making
#-----------------------------------------------
#[V]   Create the cutpoint and investigate what resulting distributions look like
#[VA]  Implement Cutpoint with Standard and/or Non-standard Certainty Levels for Indeterminate Ranges
#[VB]  Investigate resulting distributions and classifications
#-----------------------------------------------





#-----------------------------------------------
#[I] Install/load packages and load in functions
#-----------------------------------------------

#-----------------------------------------------
#[IA] Install Packages
#-----------------------------------------------

install.packages("mixsmsn")
install.packages("sn")

#-----------------------------------------------
#[IB] Load Packages
#-----------------------------------------------

library("mixsmsn")
library("sn")


#-----------------------------------------------
#[IC] Load Functions
#-----------------------------------------------

#-----------------------------------------------
#[IC1]  fitloops
#-----------------------------------------------

#Runs loops of mixture model over different numbers of components and the normal and skew-normal distribution. This function accepts a data frame of the data and identifiers (variables must be entitled data and id respectively). This is the longest function to run and may fail, but can be re-run with adjustments using various options based on the error message output when the function fails. The default values are to do both skew-normal and normal approximation, with up to 5 components and 10 loops per combination of distribution and number of components. It is not reccommended to change these unless the default function fails.
#accepts a dataframe with columns id and data
fitloops<-function(datawithids,loops=10,maxcp=5){
  dists="Both"
  datawithids2<-unique(datawithids)
  if (nrow(unique(datawithids))!=nrow(datawithids)){
    stop("The input data (datawithids) must not include duplicates. Please remove duplicates.")
  }
  data<-datawithids$data
  datapts<-length(data)
  norms<-vector("list")
  errors<-matrix(0, nrow=maxcp, ncol=2)
  loopErr = function(data,i,j,family) {
    tryCatch(suppressWarnings(smsn.mix(data, nu =3, g = i, get.init = TRUE, criteria = TRUE, group = TRUE, family = family, calc.im=FALSE, obs.prob=TRUE, iter.max=5000, kmeans.param=list(iter.max = 50, n.start = 5))),
             error = function(e) {message(paste("Error in loop", j, "of", family, "with", i, "components", sep=" ")); 
               NA}) 
  }
  
  for (i in 1:maxcp){
    norms$mc[[i]]<-vector("list", maxcp)
    norms$ms[[i]]<-vector("list", maxcp)
    j=1
    while (max(errors)<=2 & j<=loops) {
      
      norms$mc[[i]][[j]]<-vector("list", 15)
      norms$ms[[i]][[j]]<-vector("list", 15)
      if (dists=="Both") {
        norms$mc[[i]][[j]]<-loopErr(data,i,j,"Normal")  
        norms$ms[[i]][[j]]<-loopErr(data,i,j,"Skew.normal")
      } else if (dists=="Normal"){ 
        norms$mc[[i]][[j]]<-loopErr(data,i,j,"Normal")
      } else if (dists=="Skew.normal"){
        norms$ms[[i]][[j]]<-loopErr(data,i,j,"Skew.normal")
      } else {
        stop(print("dists must be specified as 'Both' or 'Normal' or 'Skew-normal'"))
      }
      
      message(paste(i, "Component(s), Loop", j))
      if (is.na(norms$mc[[i]][[j]][1])==T) {
        norms$mc[[i]][[j]]<-vector("list")
        norms$mc[[i]][[j]]$bic<-NA
        errors[i,1]<-errors[i,1]+1
      } 
      if (is.na(norms$ms[[i]][[j]][1])==T) {
        norms$ms[[i]][[j]]<-vector("list")
        norms$ms[[i]][[j]]$bic<-NA
        errors[i,2]<-errors[i,2]+1
      } 
      if (dists=="Normal"){ 
        norms$ms[[i]][[j]]<-vector("list")
        norms$ms[[i]][[j]]$bic<-NA
      } else if (dists=="Skew.normal"){
        norms$mc[[i]][[j]]<-vector("list")
        norms$mc[[i]][[j]]$bic<-NA
      }
      j=j+1
      ic=i
    }
  }

  if (i==maxcp & j==loops+1 & max(errors[maxcp, 2])<=2 & max(errors[maxcp, 1])<=2 ) {
    nobj<-vector("list")
    nobj$mc<-matrix(nrow=maxcp, ncol=loops)
    nobj$ms<-matrix(nrow=maxcp, ncol=loops)
    
    for (i in 1:maxcp){
      for (j in 1:loops){
        nobj$mc[i,j]<-norms$mc[[i]][[j]]$bic
        nobj$ms[i,j]<-norms$ms[[i]][[j]]$bic
      }
    }
    
    for (i in 1:maxcp){
      if (dists=="Both" | dists=="Normal"){
        nobj$lmin$mc[i]<-which.min(nobj$mc[i,])
      } else {
        nobj$lmin$mc[i]<-1
      }
      if (dists=="Both" | dists=="Skew.normal"){
        nobj$lmin$ms[i]<-which.min(nobj$ms[i,])
      } else {
        nobj$lmin$ms[i]<-1
      }
    }
    
    fitres<-vector("list")
    for (i in 1:maxcp){
      fitres$no[[i]]<-vector("list", 15)
      fitres$sn[[i]]<-vector("list", 15)
      fitres$no[[i]]<-norms$mc[[i]][[nobj$lmin$mc[i]]]
      fitres$sn[[i]]<-norms$ms[[i]][[nobj$lmin$ms[i]]]
    }
  
    nobj<-vector("list")
    
    for (i in 1:maxcp){
      for (j in 1:loops){
        nobj$no[i]<-fitres$no[[i]]$bic
        nobj$sn[i]<-fitres$sn[[i]]$bic
      }
    }
    
    nobj$minind<-data.frame(no=as.numeric(c(NA)), sn=as.numeric(c(NA)))
    
    
    
    if (is.na(nobj$no[1])==T){
      nobj$minind$no<-1
    } else {
      nobj$minind$no <- (which.min(nobj$no))
    }
    
    if (is.na(nobj$sn[1])==T){
      nobj$minind$sn<-1
    } else {
      nobj$minind$sn <- (which.min(nobj$sn))
    }
    
    nobj$min<-data.frame(type=c("no", "sn"), bic=as.numeric(c(NA, NA)))
    
    
    nobj$min$type<-as.character(c("no", "sn"))
    type<-c("no", "sn")
    nobj$min$type<-type
    
    nobj$min[1,2]<-fitres$no[[nobj$minind$no[1]]]$bic
    nobj$min[2,2]<-fitres$sn[[nobj$minind$sn[1]]]$bic
    
    nobj$distnum<-which.min(nobj$min$bic)
    nobj$dist<-nobj$min$type[nobj$distnum]
    
    nobj$best<-fitres[[nobj$dist]][[nobj$minind[nobj$dist][1,1]]]
    
    if (nobj$dist=="no"){
      nobj$distp="Normal"
      nobj$distf="Normal"
    } else if (nobj$dist=="sn") {
      nobj$distp="Skew-normal"
      nobj$distf="Skew-normal"
    } else {
      nobj$distp=""
    }
    
    model<-vector("list")
    
    model$loops<-loops
    model$maxcp<-maxcp
    model$datawithids<-datawithids
    model$errors<-errors
    model$loopfits<-fitres
    
    model$sn2<-fitres$sn[[2]]
    model$snbest<-fitres$sn[[nobj$minind$sn[1]]]
    model$norm2<-fitres$no[[2]]
    model$normbest<-fitres$no[[nobj$minind$no[1]]]
    model$best<-nobj$best
    model$bictab<-data.frame(normal=c(NA)[rep(c(1), times=maxcp)], skew.normal=c(NA)[rep(c(1), times=maxcp)])
    model$bictab$normal<-nobj$no
    model$bictab$skew.normal<-nobj$sn
    componentnames<-vector()
    for (i in 1:maxcp){
      componentnames[i]<-paste(i, " Component(s):  ", sep="")  
    }    
    rownames(model$bictab)<-componentnames
    colnames(model$bictab)<-c("Normal", "Skew-normal")
    model$summary<-data.frame(desc=c(paste(nobj$distp, "with", nobj$minind[nobj$dist][1,1], "components ", sep=" "), paste("Skew-normal with", nobj$minind$sn, "components ", sep=" "), paste("Normal with", nobj$minind$no, "components ", sep=" "), "Skew-normal with 2 components ", "Normal with 2 components "), bic=c(nobj$best$bic, nobj$sn[nobj$minind$sn], nobj$no[nobj$minind$no], nobj$sn[2], nobj$no[2]))
    rownames(model$summary)<-c("Best Overall", "Best Skew-normal  ", "Best Normal", "Two Skew-normal", "Two Normal")
    colnames(model$summary)<-c("Description         ", "BIC  ")
    
    model$bestdesc<-list(ncomp=nobj$minind[nobj$dist][1,1], dist=nobj$distf)
    
    print("")
    print(paste ("The best model by BIC is", nobj$distp, "with", nobj$minind[nobj$dist][1,1], "components", sep = " "))
    print("")
    print("Bayesian Information Criterion (BIC) Matrix")
    print("BIC should be minimized and a difference of 10 BIC indicates strong evidence that the model with lower BIC is superior")
    print(model$bictab)
    print("Below is a table of the BIC's of the most common distribution and number of component combinations to base a cutpoint on")
    print(model$summary)
  
  return(model)
    
  } else {
    if (which.max(errors)<=maxcp){
      dist<-"Normal"
    } else {
      dist<-"Skew.normal"
    }
    if (ic<=2){
      stop(paste("More than two loops of the mixture modeling algorithm failed to converge for the", dist, "distribution with", ic, "component(s). First ensure that the necessary packages have been loaded correctly. If the packages are loaded correctly, check for extreme outliers (which may be dropped with extreme caution) or violations of skew-normal or normal distribution assumptions."))
    } else {
      stop(paste("More than two loops of the mixture modeling algorithm failed to converge for the", dist, "distribution with", ic, "component(s). Consider choosing a lower number of components (the default is 5) via the maxcp option of the fitloops function."))
    }
  }

}
  

#-----------------------------------------------
#[IC2]  bicgraph
#-----------------------------------------------
#Plots BIC to show which combination of distributions and number of components is optimal by BIC. BIC is a criterion which should be minimized (ie we are looking for the most negative value), and a decrease in BIC of 10 or more is "strong evidence" that the model with lower BIC is superior. However in some cases you may want to choose a different distribution and this graph allows for an easy compairison of relative BIC.

bicgraph<-function(fitobj,title="BIC by type and number of components",setcolor=c("magenta", "blue", "black")){
  model<-fitobj
  color<-setcolor
  buffer<-(max(model$bictab)-min(model$bictab))/20
  ylim<-c(min(model$bictab)-buffer, max(model$bictab)+buffer)
  if (is.na(model$bictab$Normal[1])==F) {
    plot(model$bictab$Normal, pch=1, col=color[1], main=title, xlab="Number of components", ylab="BIC", type="o",  xaxt="n", ylim=ylim)
    lines(model$bictab$`Skew-normal`, type="o", pch=2, lty=2, col=color[2])
  } else {
    plot(model$bictab$Normal, type="o", pch=2, lty=2, col=color[2], main=title, xlab="Number of components", ylab="BIC",  xaxt="n", ylim=ylim)
    lines(model$bictab$`Skew-normal`, pch=1, col=color[1], type="o")    
  }  
  points(model$bestdesc$ncomp,model$best$bic,type="o", pch="O", col=color[3],cex=2)
  axis(1, at = seq(1, length(model$bictab$Normal), by = 1))
  suppressWarnings(legend("topright", col=c(color[1], color[2], color[3]), lty=c(1,2,0), pch = c(1, 2, 1), legend=c("Normal", "Skew-normal", "Best by BIC"), cex=c(1,1,1), xjust=1, seg.len=3, title="Distribution"))
  
}

#-----------------------------------------------
#[IC3]  modelpick
#-----------------------------------------------
#Selects a specific distribution and format so that it can be read into the uncertainty functions. Default picks the optimized BIC function from bestfits, but this can be adjusted if the user chooses otherwise. The function also performs uncertainty calculations and calculates raw percentages and counts positive and negative for all possible cut points between distributions

modelpick<-function(fitobj, dist="", ncomp=NA){
  #fitpick<-function(fitres,bestfits,dist="",ncomp=NA){
  fitres<-fitobj
  if (dist=="" & is.na(ncomp)==T) {
    dist=fitres$bestdesc$dist
    ncomp=fitres$bestdesc$ncomp
  }
  if (ncomp==1) {
    stop("Single component specified so there is no cut-point")
  }
  
  data=fitres$datawithids$data
  n=ncomp
  if (dist!="Normal" & dist!="Skew-normal") {
    stop("Distribution must be specified as either 'Normal' or 'Skew-normal'")
  } else if ((dist=="Normal" & is.na(fitres$loopfits$no[[1]]$bic)==T)| (dist=="Skew-normal" & is.na(fitres$loopfits$sn[[1]]$bic)==T) | ncomp>length(fitres$loopfits$no) ){
    stop("Must choose a distribution and number of components that were fitted in the fitloops function")
  } else {
    if (dist=="Normal"){
      dist="no"
    } else {
      dist="sn"
    }
    
    if (dist=="no"){
      distp="Normal"
    } else if (dist=="sn") {
      distp="Skew-Normal"
    } else {
      distp=""
    }
    
    singlefit<-vector("list")
    singlefit<-fitres$loopfits[[dist]][[ncomp]]  
    
    comp<-data.frame(ogroup=singlefit$group, data=data)
    comp$index <- as.numeric(row.names(comp))
    comp <- comp[order(comp$data),] 
    comp$diff<-c(0, diff(comp$ogroup))
    nl=as.numeric(dim(comp)[1])
    
    #Will recode the distributions so that they are in order with 1 being the lowest
    comp$c<-vector(mode="numeric", length=nl)
    comp$c[1]<-1
    for (i in 2:nl){
      if (comp$diff[i]==0) {
        comp$c[i]=comp$c[i-1]
      } else {
        comp$c[i]=comp$c[i-1]+1
      }
    }
    comp <- comp[order(comp$index),]
    singlefit$group2<-comp$c
    comp <- comp[order(comp$data),] 
    
    intDistErr = function() {
      tryCatch(data.frame(c=unique(comp$c), o=unique(comp$ogroup)),
               error = function(e) {data.frame(c=c(NA)[rep(c(1), times=length(comp$c))], o=c(NA)[rep(c(1), times=length(comp$c))])}) 
    }
    
    cross<-intDistErr()
    if (suppressWarnings(is.na(cross[1,1])==T)){
      stop("This combination of distribution and number of components yeilds interrupting distributions, which make determining a cut-point impossible. Please choose a different combination of distribution and number of components")
    }
    
    for (i in 1:n){
      singlefit$mu[i]<-fitres$loopfits[[dist]][[ncomp]]$mu[cross$o[i]]
      singlefit$sigma2[i]<-fitres$loopfits[[dist]][[ncomp]]$sigma2[cross$o[i]]
      singlefit$shape[i]<-fitres$loopfits[[dist]][[ncomp]]$shape[cross$o[i]]
      singlefit$pii[i]<-fitres$loopfits[[dist]][[ncomp]]$pii[cross$o[i]]
      singlefit$mu[i]<-fitres$loopfits[[dist]][[ncomp]]$mu[cross$o[i]]
    }
    
    singlefit$obs.prob<-NULL
    
    #Detect if there are any interupting distributions (can also be recoded to kick someone out and tell them to pick something else:
    if (max(comp$c)>max(comp$ogroup)) {
      print("Warning: there is at least one interupting distribution, which may make the determination of a cutoff and the calculation of an indeterminate range difficult in later functions. If later functions fail, try a different combination of distribution and number of components (ideally the next most optimal by BIC).")
    }
    
    par<-vector("list")
    par<-data.frame(mean=c(NA)[rep(c(1), times=n)], sd=c(NA)[rep(c(1), times=n)], shape=c(NA)[rep(c(1), times=n)], pii=c(NA)[rep(c(1), times=n)])
    for (i in 1:n){
      par$mean[i]=singlefit$mu[i]
      par$sd[i]=sqrt(singlefit$sigma2[i])
      par$sigma2[i]=singlefit$sigma2[i]
      par$shape[i]=singlefit$shape[i]
      par$pii[i]=singlefit$pii[i]
    }
  }

  #Performs uncertainty calculations.
  #uncert<-function(fitpickobj){
    
  fit<-singlefit
  #Create Uncertainties
  v<-c(seq(min(data),max(data),0.0001)) 
  
  nl<-length(v)
  PD<-vector(mode="numeric", length=n)
  PD[1]<-1
  g1<-0
  g2<-0
  Pv.D<-matrix(0, nl, n)
  Pv.DPD<-matrix(0, nl, n)
  PDT<-matrix(0, nl, 1)
  PD.v<-matrix(0, nl, n)
  for (i in 2:n) {
    PD[i]<-par$pii[i]
    PD[1]<-PD[1]-PD[i]
  }
  for (i in 1:n){
    if (dist=="Normal"){
      Pv.D[,i]<-dnorm(v,mean=par$mean[i], sd=par$sd[i])
    } else {
      Pv.D[,i]<-dsn(v,xi=par$mean[i],omega=par$sd[i],alpha=par$shape[i])
    }
  }
  for (i in 1:n){
    for (j in 1:nl){
      Pv.DPD[j,i]<-Pv.D[j,i]*PD[i]
    }
  }
  PDT<-rowSums (Pv.DPD)
  for (i in 1:n){
    for (j in 1:nl){
      PD.v[j,i] <- Pv.DPD[j,i] / (PDT[j])
    }
  }
  
  mcluncert<-1-apply(PD.v, 1, max) 
  uncertmat<-PD.v
  restrmat<-cbind(PD.v, v)
  uncertout<-vector("list")

  #multcut<-function(uncertobj){
  
  closest<-function(uncertdf,min,max,level){
    uncertdf2 <- uncertdf[ which(uncertdf$v>=min & uncertdf$v<max), ]
    diff.cut<-abs(uncertdf2$uncertainty.v - level)
    min.diffc<-min(diff.cut)
    y<-which(min.diffc==diff.cut)
    #take first one
    z<-min(y)
    cut<-uncertdf2$v[z]  
    return(cut) 
  }
  g1<-vector("list")
  g2<-vector("list")
  uncertdf<-vector("list")
  uncertainty.v<-vector("list")
  cutpoint<-c(0)[rep(c(1), times=ncomp-1)]
  for (j in 1:(ncomp-1)){
    g1[[j]]<-0
    g2[[j]]<-0
    for (i in 1:j){
      g1[[j]]<-g1[[j]]+PD.v[,i]
    }
    for (i in (j+1):n){
      g2[[j]]<-g2[[j]]+PD.v[,i]
    }
    uncertainty.v[[j]]<-1-pmax(g1[[j]], g2[[j]])
  }
  
  for (j in 1:(ncomp-1)){
    uncertdf[[j]]<-data.frame(v, uncertainty.v=uncertainty.v[[j]])
    cutpoint[j]<-closest(uncertdf[[j]],min(data),max(data),1.0)
  }
  class2<-vector("list")
  for (j in 1:(ncomp-1)){
    class2[[j]]<-ifelse(data>=cutpoint[j], "positive", "negative")
  }
  class<- data.frame(matrix(unlist(class2), nrow=singlefit$n, byrow=F),stringsAsFactors=TRUE)
  class$data<-data
  
  table<-vector("list", ncomp-1)
  for (i in 1:(ncomp-1)){
    table[[i]]<-table(class2[[i]])
    
  }
  npos<-vector("numeric")
  nneg<-vector("numeric")
  dpos<-vector("numeric")
  dneg<-vector("numeric")
  ppos<-vector("numeric")
  pneg<-vector("numeric")
  grpname<-vector("numeric")
  for (i in 1:(ncomp-1)){
    npos[i]<-table[[i]][2]
    nneg[i]<-table[[i]][1]
    dpos[i]<-npos[i]/(npos[i]+nneg[i])
    dneg[i]<-1-dpos[i]
    ppos[i]<-paste(round(100*dpos[i], 2), "%", sep="")
    pneg[i]<-paste(round(100*dneg[i], 2), "%", sep="")
    grpname[i]<-paste("Cut between components ", i, " and ", i+1, sep="")
  }
  classtab<-data.frame( nneg, npos, dneg, dpos)
  outtab<-data.frame(nneg, npos, pneg, ppos)
  rownames(outtab) <- grpname
  colnames(outtab) <-c("Number Negative", "Number Positive", "Percent Negative", "Percent Positive")
  
  cutobj<-vector("list")
  
  cutobj$singlefit<-singlefit
  cutobj$par<-par
  cutobj$datawithids<-fitres$datawithids
  cutobj$desc<-vector("list")
  cutobj$desc$ncomp<-ncomp
  cutobj$desc$dist<-distp
  
  cutobj$uncertmat<-uncertmat
  cutobj$mcluncert<-mcluncert
  cutobj$restrmat<-restrmat
  cutobj$v<-v
  
  cutobj$cutpoint<-cutpoint
  cutobj$class<-class
  cutobj$uncertdf<-uncertdf
  cutobj$type<-"Standard"
  cutobj$g1<-g1
  cutobj$g2<-g2
  cutobj$classtab<-classtab
  cutobj$outtab<-outtab
  cutobj$cutnames<-grpname

  print(outtab)  
  
  return(cutobj)
}
  
#-----------------------------------------------
#[IC4]  rawuncertgraph
#-----------------------------------------------
#Displays the overall uncertainty function (which is 1-max(probability of membership to component i)) with the possible cut-points overlaid. Note: Sometimes the cutpoint does not line up exactly with a peak in the uncertainty function, which looks like a flaw, but is not. In a two component setting, the cutpoint will always occur at the peak of the uncertainty but this is not always the case when there are more than two components.
rawuncertgraph<-function(modelpickobj,xlab="X Value",xlim=c(NA,NA),suppresslegend=F,title="default",setcolor=c(
  "green4",          #1
  "turquoise3",      #2
  "royalblue2",      #3
  "darkorchid2",     #4
  "firebrick3",      #5
  "darkorange2",     #6
  "white",           #7
  "grey",            #8
  "black"            #9
)){
  
  if (title!="default"){
    title=title
  } else {
    title=paste("Uncertainty Plot with Possible Cutpoints for", modelpickobj$desc$dist, "with", modelpickobj$desc$ncomp, "components", sep = " ")
  }
  grobj<-list()
  if (is.na(xlim[1])==T | is.na(xlim[2])==T) {   
    grobj$minval<-(1/10)*trunc(10*min(modelpickobj$datawithids$data)) 
    grobj$maxval<-max(modelpickobj$datawithids$data)
  } else {
    grobj$minval<-xlim[1]
    grobj$maxval<-xlim[2]
  }

  plot(modelpickobj$v, modelpickobj$mcluncert, main=title, type="l", ylim=c(0,1), xlim=c(grobj$minval,grobj$maxval), xlab=xlab, ylab="Uncertainty")
    
  col<-setcolor
  ncomp<-modelpickobj$desc$ncomp
  for (i in 1:(modelpickobj$desc$ncomp-1)){
    abline(v =modelpickobj$cutpoint[i], untf = T, col=col[i], lwd=2)
  }
  
  if (suppresslegend==F) {
    legtxt<-vector("character", ncomp)
    colors<-vector("character", ncomp)
    for (i in 2:ncomp){
      legtxt[i]<-paste(modelpickobj$cutnames[i-1], " , x=", round(modelpickobj$cutpoint[i-1], 2), sep="")
      colors[i]<-col[i-1]
    }  
    legtxt[1]<-"Uncertainty"
    colors[1]<-col[length(col)]
    linewidth<-c(2)[rep(c(1), times=ncomp)]
    legend("topright", col=colors, lwd=c(1,linewidth), legend=legtxt, xjust=1, seg.len=3)
  }
}


#-----------------------------------------------
#[IC5]  rawdistgraph
#-----------------------------------------------
#Displays the distributions of all the components for the chosen combination of distribution and number of components. 
rawdistgraph<-function(modelpickobj,xlim=c(NA,NA),xlab="Optical Density",setbreaks=100,suppresslegend=F,title="default",setcolor=c(
  "green4",          #1
  "turquoise3",      #2
  "royalblue2",      #3
  "darkorchid2",     #4
  "firebrick3",      #5
  "darkorange2",     #6
  "white",           #7
  "grey",            #8
  "black"            #9
)){
  #Parameter specification
  ncomp=modelpickobj$desc$ncomp
  grobj<-vector("list")
  grobj$gpar<-vector("list")
  data<-modelpickobj$datawithids$data
  color<-setcolor
  
  if (is.na(xlim[1])==T | is.na(xlim[2])==T) {   
    grobj$minval<-(1/10)*trunc(10*min(modelpickobj$datawithids$data)) 
    grobj$maxval<-max(data)
  } else {
    grobj$minval<-xlim[1]
    grobj$maxval<-xlim[2]
  }
  
  if (title!="default"){
    title=title
  } else {
    title=paste("Distribution of", ncomp, modelpickobj$desc$dist,"Components", sep = " ")
  }
  
  grobj$t<-seq(grobj$minval,grobj$maxval,length=1000) 
  
  grobj$gpar<-data.frame(mean=c(NA)[rep(c(1), times=ncomp)], sd=c(NA)[rep(c(1), times=ncomp)], shape=c(NA)[rep(c(1), times=ncomp)], pii=c(NA)[rep(c(1), times=ncomp)])
  grobj$gt<-vector("list", ncomp+1)
  grobj$maxdens<-vector("numeric", ncomp)
  
  grobj$coloropt<-vector("list", 5)
  grobj$coloropt[[1]]<-c(color[1], color[9])
  grobj$coloropt[[2]]<-c(color[1], color[5], color[9])
  grobj$coloropt[[3]]<-c(color[1], color[3], color[5], color[9])
  grobj$coloropt[[4]]<-c(color[1], color[3], color[4], color[5], color[9])
  grobj$coloropt[[5]]<-c(color[1], color[2], color[3], color[4], color[5], color[9])
  grobj$colors<-grobj$coloropt[[ncomp]]
  
  for (i in 1:ncomp){
    grobj$gpar$mean[i]=modelpickobj$par$mean[i]
    grobj$gpar$sd[i]=sqrt(modelpickobj$par$sigma2[i])
    grobj$gpar$shape[i]=modelpickobj$par$shape[i]
    grobj$gpar$pii[i]=modelpickobj$par$pii[i]
  }
  #Creation of points for curves of Graph at full scale (components scaled to add to overall density)
  for (i in 1:ncomp) {
    grobj$gt[[i]] <- grobj$gpar$pii[i]*dsn(grobj$t,xi=grobj$gpar$mean[i],omega=grobj$gpar$sd[i],alpha=grobj$gpar$shape[i])
  }
  grobj$gt[[ncomp+1]]<-vector("numeric", 1000)
  for (i in 1:ncomp) {
    grobj$gt[[ncomp+1]]<-grobj$gt[[ncomp+1]]+grobj$gt[[i]]
  }
  grobj$ylim<-c(0, max(grobj$gt[[ncomp+1]], max(hist(data,setbreaks,plot=F)$density,na.rm=TRUE) ))
  grobj$xlim<-c(grobj$minval,grobj$maxval)
  
  dist<-c(0)[rep(c(0), times=ncomp)]
  for (i in 1:ncomp){
    if (modelpickobj$singlefit$shape[i]>0) {
      dist[i]<-1
    }
  }
  
  #Plot Graph
  
  hist(data, freq=F, breaks=setbreaks, xlab=xlab, ylab="Density", main=title, border="grey", ylim=grobj$ylim, xlim=grobj$xlim)
  for (i in 1:ncomp) {
    lines(grobj$t, grobj$gt[[i]], lwd=2, col=grobj$colors[i])
  }
  lines(grobj$t, grobj$gt[[length(grobj$gt)]], lwd=2, col=grobj$colors[length(grobj$colors)], lty=2)
  legtxt<-vector("character", ncomp+1)
  colors<-vector("character", ncomp+1)
  for (i in 1:ncomp){
    legtxt[i]<-paste("Distribution of Component ", i, sep="")
    colors[i]<-grobj$colors[i]
  }  
  legtxt[ncomp+1]<-"Overall Distribution"
  colors[ncomp+1]<-grobj$colors[length(grobj$colors)]
  linetype<-c(1)[rep(c(1), times=ncomp)]
  if (suppresslegend==F) {
    legend("topright", col=colors, lwd=2, legend=legtxt, xjust=1, seg.len=3, title="Distributions", lty=c(linetype,2))
  }
    
  #return(grobj)
}

#-----------------------------------------------
#[IC6]  rawhistcuts
#-----------------------------------------------
#Displays which component each bar of the histogram would be assigned to. 

rawhistcuts<-function(modelpickobj,xlab="Optical Density",xlim=c(NA,NA),setbreaks=250,suppresslegend=F,title="default",setcolor=c(
  "green4",          #1
  "turquoise3",      #2
  "royalblue2",      #3
  "darkorchid2",     #4
  "firebrick3",      #5
  "darkorange2",     #6
  "white",           #7
  "grey",            #8
  "black"            #9
)){
  data<-modelpickobj$datawithids$data
  grobj<-list()
  
  if (is.na(xlim[1])==T | is.na(xlim[2])==T) {   
    grobj$minval<-(1/10)*trunc(10*min(data)) 
    grobj$maxval<-max(data)
    hist<-hist(data, breaks=setbreaks, plot=F)
  } else {
    grobj$minval<-xlim[1]
    grobj$maxval<-xlim[2]
    datalim<-data[ which(data>=grobj$minval & data<=grobj$maxval)]
    hist<-hist(datalim, breaks=setbreaks, plot=F)
  }
  
  ncomp<-modelpickobj$desc$ncomp
  
  if (title!="default"){
    title=title
  } else {
    title=paste("Classification by", modelpickobj$desc$dist, "with", ncomp, "components", sep = " ")
  }
  
  cuts <- cut(hist$breaks, c(-Inf,modelpickobj$cutpoint,Inf))
  
  color<-setcolor
  coloropt<-vector("list", 5)
  coloropt[[1]]<-c(color[1], color[9])
  coloropt[[2]]<-c(color[1], color[5], color[9])
  coloropt[[3]]<-c(color[1], color[3], color[5], color[9])
  coloropt[[4]]<-c(color[1], color[3], color[4], color[5], color[9])
  coloropt[[5]]<-c(color[1], color[2], color[3], color[4], color[5], color[9])
  col<-coloropt[[ncomp]]
  
  cuts<-as.numeric(cuts)
  colors<-vector("character", length(cuts))
  for (i in 1:max(cuts)){
    for (j in 1:length(cuts)){
      if (cuts[j]==i) {
        colors[j]<-col[i]
      }
    }
  }
  plot(hist, col=colors, xlab=xlab, ylab="Density", main=title, freq=F, lty="blank")
  
  if (suppresslegend==F){
    legtxt<-vector("character", ncomp)
    colors<-vector("character", ncomp)
    for (i in 1:ncomp){
      legtxt[i]<-paste("Component ", i, sep="")
      colors[i]<-col[i]
    }  
    linetype<-c(1)[rep(c(1), times=ncomp)]
    legend("topright", col=colors, lwd=3, legend=legtxt, xjust=1, seg.len=3, title="Components", lty=linetype)
  }
}

#-----------------------------------------------
#[IC7]  cutoff
#-----------------------------------------------
#After deciding which distribution to cut after, this function provides cutpoint and standard indeterminate ranges (80 and 90), as well as classification into negative, indeterminate and positive components. It should be noted that the function is written so that it will still run if a lower or upper bound is not found, in which case all subjects below the upper bound (if lower bound not found) or above the lower bound (if upper bound not found) will be classified as indeterminate. In this case it is reccommeded that the non-standard cuts function be used to find a more descriptive indeterminate range.
cutoff<-function(modelpickobj, cutcomp=0, standardcert=T, newcertlevel=0){
  data<-modelpickobj$datawithids$data
  n=modelpickobj$desc$ncomp
  certlevel=newcertlevel
  if (n==2){
    cutcomp=1
  } else if (cutcomp==0) {
    stop("Component to cut after must be specified. For example if you wanted to cut between components 2 and 3, set cutcomp=2")
  }
  if ((certlevel>=1 & certlevel!=0) | (certlevel<0.5 & certlevel!=0)){
    stop("Certainty level must be between 0.5 and 1, for example if you wanted a certainty level of 80% set certlevel=0.80")
  }
  if (standardcert==F & newcertlevel==0){
    stop("Neither standard nor non-standard certainty levels are requested. For standard certainty levels set standardcert=T. For non-standard certainty levels set newcertlevel to desired certainty. Both standard and nonstandard certainties can be generated in the same function call.")
  }
  uncertdf<-modelpickobj$uncertdf[[cutcomp]]
  closest<-function(uncertdf,min,max,level){
    uncertdf2 <- uncertdf[ which(uncertdf$v>=min & uncertdf$v<max), ]
    diff.cut<-abs(uncertdf2$uncertainty.v - level)
    min.diffc<-min(diff.cut)
    y<-which(min.diffc==diff.cut)
    #take first one
    z<-min(y)
    cut<-uncertdf2$v[z]  
    return(cut) 
  }
  crosses<-function(uncertdf,min,max,level){
    uncertdf2 <- uncertdf[ which(uncertdf$v>=min & uncertdf$v<max), ]
    diff.cut<-abs(uncertdf2$uncertainty.v - level)
    min.diffc<-min(diff.cut)
    if (min.diffc<=0.001){
      y<-which(min.diffc==diff.cut)
      #take first one
      z<-min(y)
      cut<-uncertdf2$v[z]  
    } else {
      cut<-NA
    }
    return(cut)
  }
  
  cutpoint<-modelpickobj$cutpoint[cutcomp]  
  minunclb<-closest(uncertdf,min(data)+0.0001,cutpoint,0.0001)
  minuncub<-closest(uncertdf,cutpoint,max(data),0)
  class<-data.frame(id=modelpickobj$datawithids$id, data=modelpickobj$datawithids$data)
  class$groupdi<-ifelse(class$data>=cutpoint, "positive", "negative")
  class$groupdi<-as.factor(class$groupdi)

  if (standardcert==T){
    lb80<-crosses(uncertdf,minunclb,cutpoint,0.2)  
    ub80<-crosses(uncertdf,cutpoint,minuncub,0.2)
    lb90<-crosses(uncertdf,minunclb,cutpoint,0.1)  
    ub90<-crosses(uncertdf,cutpoint,minuncub,0.1)
    bound80<-c(lb80, ub80)
    bound90<-c(lb90, ub90)
    if (is.na(lb80)==T){
      message("Lower bound of the 80% indeterminate range not found because uncertainty is never below 20% to the left of the cutpoint. Consider adjusting the non-standard certainty option (newcertlevel) to find a more descriptive indeterminate range.")
    }
    if (is.na(ub80)==T){
      message("Upper bound of the 80% indeterminate range not found because uncertainty is never below 20% to the right of the cutpoint. Consider adjusting the non-standard certainty option (newcertlevel) to find a more descriptive indeterminate range.")
    }
    if (is.na(lb90)==T){
      message("Lower bound of the 90% indeterminate range not found because uncertainty is never below 10% to the left of the cutpoint. Consider adjusting the non-standard certainty option (newcertlevel) to find a more descriptive indeterminate range.")
    }
    if (is.na(ub90)==T){
      message("Upper bound of the 90% indeterminate range not found because uncertainty is never below 10% to the right of the cutpoint. Consider adjusting the non-standard certainty option (newcertlevel) to find a more descriptive indeterminate range.")
    }
    class$group90<-"indeterminate"
    class$group90[class$data<=bound90[1]]<-"negative"
    class$group90[class$data>=bound90[2]]<-"positive"
    class$group80<-"indeterminate"
    class$group80[class$data<=bound80[1]]<-"negative"
    class$group80[class$data>=bound80[2]]<-"positive"
    class$group90<-as.factor(class$group90)
    class$group80<-as.factor(class$group80)
    class$group80 = factor(class$group80,levels(class$group80)[c(2,3,1)])
    class$group90 = factor(class$group90,levels(class$group90)[c(2,3,1)])
    class$group<-NA
  }
  
  if (newcertlevel>0){
    uncertlevel<-1-certlevel
    uncertlevelpct<-certlevel
    lbuns<-crosses(uncertdf,minunclb,cutpoint,uncertlevel)  
    ubuns<-crosses(uncertdf,cutpoint,minuncub,uncertlevel)
    bound<-c(lbuns, ubuns)
    if (is.na(lbuns)==T){
      message(paste("Lower bound of the ", 100*certlevel, "% indeterminate range not found because uncertainty is never below ", 100*(1-certlevel), "% to the left of the cutpoint. Consider adjusting the non-standard certainty option (newcertlevel) to find a more descriptive indeterminate range.", sep=""))
    }
    if (is.na(ubuns)==T){
      message(paste("Upper bound of the ", 100*certlevel, "% indeterminate range not found because uncertainty is never below ", 100*(1-certlevel), "% to the right of the cutpoint. Consider adjusting the non-standard certainty option (newcertlevel) to find a more descriptive indeterminate range.", sep=""))
    }
    if (standardcert==F){
      class$group90<-NA
      class$group80<-NA
    }
    class$group<-"indeterminate"
    class$group[class$data<=bound[1]]<-"negative"
    class$group[class$data>=bound[2]]<-"positive"
    class$group<-as.factor(class$group)
    class$group = factor(class$group,levels(class$group)[c(2,3,1)])
  }
  
  cutobj<-vector("list")
  cutobj$cutpoint<-cutpoint
  cutobj$bound80<-c(NA,NA)
  cutobj$bound90<-c(NA,NA)
  cutobj$bound<-c(NA,NA)
  cutobj$class<-class
  cutobj$uncertdf<-uncertdf
  
  if (standardcert==T & newcertlevel==0){
    cutobj$type<-"Standard"
    cutobj$bound80<-bound80
    cutobj$bound90<-bound90
  } else if (standardcert==F & newcertlevel>0) {
    cutobj$type<-"Non.Standard"
    cutobj$bound<-bound
  } else {
    cutobj$type<-"Standard.and.Non"
    cutobj$bound80<-bound80
    cutobj$bound90<-bound90
    cutobj$bound<-bound
  }
  
  cutobj$par<-modelpickobj$par
  cutobj$uncertobj<-modelpickobj$uncertmat
  cutobj$density<-data.frame(g1=modelpickobj$g1[[cutcomp]], g2=modelpickobj$g2[[cutcomp]])
  cutobj$desc<-modelpickobj$desc
  cutobj$desc$cutcomp<-cutcomp
  cutobj$desc$certlevel<-newcertlevel
  cutobj$datawithids<-modelpickobj$datawithids
  
  return(cutobj) 
}

#-----------------------------------------------
#[IC8] cutuncertgraph
#-----------------------------------------------
#Displays the uncertainty function after components have been combined to create positive and negative components with cut-points and bounds of indeterminate range(s) overlaid. Accepts results from both standard and nonstandard cutpoints (from the functions standindet and specindet respectively).
cutuncertgraph<-function(cutobj, xlab="Optical Density", xlim=c(NA,NA), suppresslegend=F, title="default", setcolor=c(
  "green4",          #1
  "turquoise3",      #2
  "royalblue2",      #3
  "darkorchid2",     #4
  "firebrick3",      #5
  "darkorange2",     #6
  "white",           #7
  "grey",            #8
  "black"            #9
)){
  color=setcolor
  grobj<-vector("list")
  
  if (title!="default"){
    title=title
  } else {
    title=paste("Uncertainty Plot of", cutobj$desc$dist, "with", cutobj$desc$ncomp, "components, cut between distributions", cutobj$desc$cutcomp, "and", (cutobj$desc$cutcomp+1), sep = " ")
  }
  
  if (is.na(xlim[1])==T | is.na(xlim[2])==T) {   
    plot(cutobj$uncertdf$v, cutobj$uncertdf$uncertainty.v, main=title, type="l", ylim=c(0,1), xlab=xlab, ylab="Uncertainty")
  } else {
    grobj$minval<-xlim[1]
    grobj$maxval<-xlim[2]
    plot(cutobj$uncertdf$v, cutobj$uncertdf$uncertainty.v, main=title, type="l", ylim=c(0,1), xlab=xlab, ylab="Uncertainty", xlim=c(grobj$minval,grobj$maxval))
  }
  
  abline(v =cutobj$cutpoint, untf = T, col=color[6], lwd=1)
  if (cutobj$type=="Standard") {
    abline(h =0.1, untf = T, col=color[4], lwd=1, lty=2)
    abline(h =0.2, untf = T, col=color[2], lwd=1, lty=2)
    abline(v =cutobj$bound80, untf = T, col=color[2], lwd=1)
    abline(v =cutobj$bound90, untf = T, col=color[4], lwd=1)
    legtxt<-c("Uncertainty",
              paste("Cutpoint=", round(cutobj$cutpoint, 2), sep =""), 
              "Uncertainty=0.20",
              "Indeterminate Range with", 
              paste("80% Certainty: (", round(cutobj$bound80[1], 2), " , ", round(cutobj$bound80[2], 2), ")", sep =""),
              "Uncertainty=0.10",
              "Indeterminate Range with", 
              paste("90% Certainty: (", round(cutobj$bound90[1], 2), " , ", round(cutobj$bound90[2], 2), ")", sep ="")
    )
    if (suppresslegend==F){
      legend("topright", col=c(color[9], color[6], color[2], color[2], "white", color[4], color[4], "white"), lwd=1, lty=c(1,1,2,1,1,2,1,1),legend=legtxt, xjust=1, seg.len=1, title="Uncertainties, Cutpoint, and Boundaries")
      }
  } else if (cutobj$type=="Non.Standard") {
    abline(h =(1-cutobj$desc$certlevel), untf = T, col=color[3], lwd=1, lty=2)
    abline(v =cutobj$bound, untf = T, col=color[3], lwd=1)
    legtxt<-c("Uncertainty",
              paste("Cutpoint=", round(cutobj$cutpoint, 2), sep =""), 
              paste("Uncertainty=", (1-cutobj$desc$certlevel), sep =""),
              "Indeterminate Range with", 
              paste((100*cutobj$desc$certlevel), "% Certainty: (", round(cutobj$bound[1], 2), " , ", round(cutobj$bound[2], 2), ")", sep ="")
    )
    if (suppresslegend==F){
      legend("topright", col=c(color[9],color[6], color[3], color[3], "white"), lwd=1, lty=c(1,1,2,1), legend=legtxt, xjust=1, seg.len=1, title="Cutpoint, Uncertainties, and Boundaries")
    }
  } else if (cutobj$type=="Standard.and.Non") {
    abline(h =0.1, untf = T, col=color[4], lwd=1, lty=2)
    abline(h =0.2, untf = T, col=color[2], lwd=1, lty=2)
    abline(h =(1-cutobj$desc$certlevel), untf = T, col=color[3], lwd=1, lty=2)
    abline(v =cutobj$bound80, untf = T, col=color[2], lwd=1)
    abline(v =cutobj$bound90, untf = T, col=color[4], lwd=1)
    abline(v =cutobj$bound, untf = T, col=color[3], lwd=1)
    legtxt<-c("Uncertainty",
              paste("Cutpoint=", round(cutobj$cutpoint, 2), sep =""), 
              "Uncertainty=0.20",
              "Indeterminate Range with", 
              paste("80% Certainty: (", round(cutobj$bound80[1], 2), " , ", round(cutobj$bound80[2], 2), ")", sep =""),
              "Uncertainty=0.10",
              "Indeterminate Range with", 
              paste("90% Certainty: (", round(cutobj$bound90[1], 2), " , ", round(cutobj$bound90[2], 2), ")", sep =""),
              paste("Uncertainty=", (1-cutobj$desc$uncertlevel), sep =""),
              "Indeterminate Range with", 
              paste((100*cutobj$desc$certlevel), "% Certainty: (", round(cutobj$bound[1], 2), " , ", round(cutobj$bound[2], 2), ")", sep ="")
    )
    if (suppresslegend==F){
      legend("topright", col=c(color[9], color[6], color[2], color[2], "white", color[4], color[4], "white", color[3], color[3], "white"), lwd=1, lty=c(1,1,2,1,1,2,1,1,2,1,1),legend=legtxt, xjust=1, seg.len=1, title="Uncertainties, Cutpoint, and Boundaries")
    }
  }
}

#-----------------------------------------------
#[IC9]  cutdistgraph
#-----------------------------------------------
#Displays the distributions of the positive and negative components with cut-points and bounds of indeterminate range(s) overlaid. Accepts results from both standard and nonstandard cutpoints (from the functions standindet and specindet respectively). 
cutdistgraph<-function(cutobj,xlim=c(NA,NA),xlab="Optical Density",setbreaks=100,suppresslegend=F,title="default",setcolor=c(
  "green4",          #1
  "turquoise3",      #2
  "royalblue2",      #3
  "darkorchid2",     #4
  "firebrick3",      #5
  "darkorange2",     #6
  "white",           #7
  "grey",            #8
  "black"            #9
)){
  #Parameter specification
  ncomp<-cutobj$desc$ncomp
  cutcomp<-cutobj$desc$cutcomp
  cutobj<-cutobj
  data<-cutobj$datawithids$data
  color<-setcolor
  
  grobj<-vector("list")
  grobj$gpar<-vector("list")
  
  if (is.na(xlim[1])==T | is.na(xlim[2])==T) { 
    grobj$minval<-(1/10)*trunc(10*min(data))     
    grobj$maxval<-max(data)
  } else {
    grobj$minval<-xlim[1]
    grobj$maxval<-xlim[2]
  }
  
  if (title!="default"){
    title=title
  } else {
    title=paste(cutobj$desc$dist, "with", cutobj$desc$ncomp, "components, cut between distributions", cutobj$desc$cutcomp, "and", (cutobj$desc$cutcomp+1), sep = " ")
  }

  grobj$t<-seq(grobj$minval,grobj$maxval,length=1000) 
  
  grobj$gpar<-data.frame(mean=c(NA)[rep(c(1), times=ncomp)], sd=c(NA)[rep(c(1), times=ncomp)], shape=c(NA)[rep(c(1), times=ncomp)], pii=c(NA)[rep(c(1), times=ncomp)])
  grobj$gs<-vector("list", ncomp+1)
  grobj$gt<-vector("list", 2+1)
  grobj$maxdens<-vector("numeric", ncomp)
  grobj$colors<-c(setcolor[1], setcolor[5], setcolor[9])
  grobj$t<-cutobj$uncertdf$v
  
  for (i in 1:ncomp){
    grobj$gpar$mean[i]=cutobj$par$mean[i]
    grobj$gpar$sd[i]=sqrt(cutobj$par$sigma2[i])
    grobj$gpar$shape[i]=cutobj$par$shape[i]
    grobj$gpar$pii[i]=cutobj$par$pii[i]
  }
  #Creation of points for curves of Graph at full scale (components scaled to add to overall density)
  for (i in 1:ncomp) {
    grobj$gs[[i]] <- grobj$gpar$pii[i]*dsn(grobj$t,xi=grobj$gpar$mean[i],omega=grobj$gpar$sd[i],alpha=grobj$gpar$shape[i])
  }
  grobj$gt[[1]]<-0
  grobj$gt[[2]]<-0
  grobj$gt[[3]]<-0
  for (i in 1:cutcomp){
    grobj$gt[[1]]<-grobj$gt[[1]]+grobj$gs[[i]]
  }
  for (i in (cutcomp+1):ncomp){
    grobj$gt[[2]]<-grobj$gt[[2]]+grobj$gs[[i]]
  }
  for (i in 1:ncomp) {
    grobj$gt[[3]]<-grobj$gt[[3]]+grobj$gs[[i]]
  }
  
  grobj$ylim<-c(0, max(grobj$gt[[2+1]], max(hist(data,setbreaks,plot=F)$density,na.rm=TRUE) ))
  grobj$xlim<-c(round(min(data),digits=1),grobj$maxval)
  
  #Plot Graph
  hist(data, freq=F, breaks=setbreaks, xlab=xlab, ylab="Density", main=title, border=color[8], ylim=grobj$ylim, xlim=grobj$xlim)
  
  for (i in 1:2) {
    lines(cutobj$uncertdf$v, grobj$gt[[i]], lwd=2, col=grobj$colors[i])
  }
  lines(grobj$t, grobj$gt[[length(grobj$gt)]], lwd=2, col=grobj$colors[length(grobj$colors)], lty=2)
  abline(v =cutobj$cutpoint, untf = T, col=color[6], lwd=1)
  if (cutobj$type=="Standard") {
    abline(v =cutobj$bound80, untf = T, col=color[2], lwd=1)
    abline(v =cutobj$bound90, untf = T, col=color[4], lwd=1)
    if (suppresslegend==F) {
      legtxt<-c("Negative Distribution",
                "Positive Distribution",
                "Overall Distribution",
                paste("Cutpoint=", round(cutobj$cutpoint, 2), sep =""), 
                "Indeterminate Range with", 
                paste("80% Certainty: (", round(cutobj$bound80[1], 2), " , ", round(cutobj$bound80[2], 2), ")", sep =""),
                "Indeterminate Range with", 
                paste("90% Certainty: (", round(cutobj$bound90[1], 2), " , ", round(cutobj$bound90[2], 2), ")", sep ="")
      )
      legend("topright", col=c(color[1], color[5], color[9], color[6], color[2], "white", color[4], "white"), lwd=c(2,2,2,1,1,1,1,1), legend=legtxt, xjust=1, seg.len=3, title="Distributions, Cutpoint, and Boundaries", lty=c(1,1,2,1,1,1,1,1))
    }
  } else if (cutobj$type=="Non.Standard") {
    abline(v =cutobj$bound, untf = T, col=color[3], lwd=1)
    if (suppresslegend==F){
      legtxt<-c("Negative Distribution",
                "Positive Distribution",
                "Overall Distribution",
                paste("Cutpoint=", round(cutobj$cutpoint, 2), sep =""), 
                "Indeterminate Range with", 
                paste((cutobj$desc$certlevel*100), "% Certainty: (", round(cutobj$bound[1], 2), " , ", round(cutobj$bound[2], 2), ")", sep ="")
      )
      legend("topright", col=c(color[1], color[5], color[9], color[6], color[3], "white"), lwd=c(2,2,2,1,1,1), legend=legtxt, xjust=1, seg.len=3, title="Distributions, Cutpoint, and Boundaries", lty=c(1,1,2,1,1,1))
    }
  } else {
    abline(v =cutobj$bound80, untf = T, col=color[2], lwd=1)
    abline(v =cutobj$bound90, untf = T, col=color[4], lwd=1)
    abline(v =cutobj$bound, untf = T, col=color[3], lwd=1)
    if (suppresslegend==F){
    legtxt<-c("Negative Distribution",
              "Positive Distribution",
              "Overall Distribution",
              paste("Cutpoint=", round(cutobj$cutpoint, 2), sep =""), 
              "Indeterminate Range with", 
              paste("80% Certainty: (", round(cutobj$bound80[1], 2), " , ", round(cutobj$bound80[2], 2), ")", sep =""),
              "Indeterminate Range with", 
              paste("90% Certainty: (", round(cutobj$bound90[1], 2), " , ", round(cutobj$bound90[2], 2), ")", sep =""),
              "Indeterminate Range with", 
              paste((100*cutobj$desc$certlevel), "% Certainty: (", round(cutobj$bound[1], 2), " , ", round(cutobj$bound[2], 2), ")", sep ="")
      )
    legend("topright", col=c(color[1], color[5], color[9], color[6], color[2], "white", color[4], "white", color[3], "white"), lwd=c(2,2,2,1,1,1,1,1,1,1), legend=legtxt, xjust=1, seg.len=3, title="Distributions, Cutpoint, and Boundaries", lty=c(1,1,2,1,1,1,1,1,1,1))
    }
  }

  #return(grobj)
  
}

#-----------------------------------------------
#[IC10]  summaryout
#-----------------------------------------------

#Creates table which summarizes the results of the cutting functions, yeilding the cutpoint, indeterminate ranges, counts and percentages for the dichotomous cut, 80 and 90 indeterminate ranges and (if desired) one non-standard uncertainty level.
summaryout<-function(cutobj,fileandpathname=NULL){
  
  if (is.null(cutobj)==T) {
    print("Cutoff object must be specified")
  } 
  if (cutobj$type=="Standard"){
    dataout<-cutobj$class[c("id", "data", "group80", "group90","groupdi")]
    names(dataout)<-c(names(cutobj$class[c(1,2,5,4,3)]))
  } else if (cutobj$type=="Non.Standard"){
    dataout<-cutobj$class[c("id", "data", "group", "groupdi")]
    names(dataout)<-c(names(cutobj$class[c(1,2)]),paste('group', (100*cutobj$desc$certlevel), sep=''), names(cutobj$class[3]))
  } else {
    dataout<-cutobj$class[c("id", "data","group80", "group90", "group", "groupdi")]
    names(dataout)<-c(names(cutobj$class[c(1,2,5,4)]),paste('group', (100*cutobj$desc$certlevel), sep=''), names(cutobj$class[3]))
  }

  if (is.null(fileandpathname)==F){
    write.csv(dataout, file =paste(fileandpathname, ".csv", sep=""), row.names=F)
  }

  outdataobj<-dataout

  
  comptable <- function(table){
    
    dftable<-as.data.frame(table)
    len<-nrow(dftable)
    dftable$orders<-c(NA)[rep(c(1), times=len)]
    for (i in 1:len){
      if (as.character(dftable$Var1[i])=="negative"){
        dftable$orders[i]<-0
      } else if (as.character(dftable$Var1[i])=="indeterminate"){
        dftable$orders[i]<-1
      } else if (as.character(dftable$Var1[i])=="positive"){
        dftable$orders[i]<-2
      }
    }
    row.names(dftable)<-dftable$Var1
    dftable$Var1<-NULL
    if (is.na(dftable["indeterminate", 1])==T) {
      Freq<-c(dftable$Freq,0)
      grps2<-as.data.frame(Freq)
      x<-c(row.names(dftable),"indeterminate")
      row.names(grps2)<-x
      grps2$orders<-c(dftable$orders,1)
      dftable<-grps2
      Freq<-NULL
      x<-NULL
    }
    if (is.na(dftable["negative", 1])==T) {
      Freq<-c(dftable$Freq,0)
      grps2<-as.data.frame(Freq)
      x<-c(row.names(dftable),"negative")
      row.names(grps2)<-x
      grps2$orders<-c(dftable$orders,0)
      dftable<-grps2
      Freq<-NULL
      x<-NULL
    }
    if (is.na(dftable["positive", 1])==T) {
      Freq<-c(dftable$Freq,0)
      grps2<-as.data.frame(Freq)
      x<-c(row.names(dftable),"positive")
      row.names(grps2)<-x
      grps2$orders<-c(dftable$orders,2)
      dftable<-grps2    
      Freq<-NULL
      x<-NULL
    }
    dftable<-dftable[order(dftable$orders),]
    return(dftable)
  }
  
  grps<-table(outdataobj$groupdi)
  grps<-comptable(grps)
  grps$group<-row.names(grps)
  grps$groupdi<-grps$Freq
  grps$Freq<-NULL
  
  if (ncol(outdataobj)!=5 ){
    grpsns<-table(outdataobj[,4])
    grpsns<-comptable(grpsns)
    grps$gns<-grpsns$Freq
    names(grps)<-c(names(grps[1:3]),paste(names(outdataobj)[4]))
  } 
  
  if (ncol(outdataobj)!=4 ){ 
    grps80<-table(outdataobj$group80)
    grps80<-comptable(grps80)
    grps$group80<-grps80$Freq
    grps90<-table(outdataobj$group90)
    grps90<-comptable(grps90)
    grps$group90<-grps90$Freq
  } 
  
  names(grps)
  
  for (i in names(grps)){
    nneg<-vector("numeric")
    nind<-vector("numeric")
    npos<-vector("numeric")
    
    dneg<-vector("numeric")  
    dind<-vector("numeric")
    dpos<-vector("numeric")
    
    pneg<-vector("numeric")
    pind<-vector("numeric")
    ppos<-vector("numeric")
  }
  
  vars<-c("groupdi", "group80", "group90", "group80")
  
  grpname<-vector("numeric")
  for (i in 1:(ncol(grps)-2)){
    nneg[i]<-grps[1,i+2]
    nind[i]<-grps[2,i+2] 
    npos[i]<-grps[3,i+2]
  }
  for (i in 1:(ncol(grps)-2)){
    dneg[i]<-nneg[i]/(npos[i]+nind[i]+nneg[i])
    dind[i]<-nind[i]/(npos[i]+nind[i]+nneg[i])
    dpos[i]<-npos[i]/(npos[i]+nind[i]+nneg[i])
    ppos[i]<-paste(round(100*dpos[i], 2), "%", sep="")
    pind[i]<-paste(round(100*dind[i], 2), "%", sep="")
    pneg[i]<-paste(round(100*dneg[i], 2), "%", sep="")
  }
  if (ncol(grps)==4){
    grpname<-c("Raw Cutoff", paste("Cutoff with ", cutobj$desc$certlevel*100, "% Certainty", sep=""))
  } else if (ncol(grps)==5) {
    grpname<-c("Raw Cutoff", "Cutoff with 80% Certainty", "Cutoff with 90% Certainty")
  } else if (ncol(grps)==6) {
    grpname<-c("Raw Cutoff", paste("Cutoff with ", cutobj$desc$certlevel*100, "% Certainty", sep=""), "Cutoff with 80% Certainty", "Cutoff with 90% Certainty")
  }
  
  classtab<-data.frame(nneg, nind, npos, dneg, dind, dpos)
  outtab<-data.frame(nneg, nind, npos, pneg, pind, ppos)
  rownames(outtab) <- grpname
  colnames(outtab) <-c("No. Neg", "No. Indet", "No. Pos", "% Neg", "% Indet", "% Pos")
  
  model<-vector("list")
  if (is.null(cutobj)==F){
    model$cutpoint<-cutobj$cutpoint
    model$desc<-cutobj$desc
    model$bounds<-vector("list")
    model$bounds$boundns<-cutobj$bound
  } else {
    model$cutpoint<-cutobj$cutpoint     
    model$desc<-cutobj$desc
    model$bounds<-vector("list")
    model$bounds$bound80<-cutobj$bound80
    model$bounds$bound90<-cutobj$bound90
  }
  if (is.null(cutobj)==F & is.null(cutobj)==F){
    model$bounds$bound80<-cutobj$bound80
    model$bounds$bound90<-cutobj$bound90
  }
  model$classtab<-classtab
  model$outtab<-outtab
  model$classification<-outdataobj
  
  message(paste("The absolute cutoff is", round(cutobj$cutpoint,3), sep=" "))
  message("Boundries of Indeterminates:")
  if (cutobj$type!="Standard"){
    message(paste(cutobj$desc$certlevel*100, "% Certainty Range: (", round(cutobj$bound[1], 3), ", ", round(cutobj$bound[2], 3), ")", sep=""))
  }
  if (cutobj$type!="Non.Standard"){  
    message(paste("80% Certainty Range: (", round(cutobj$bound80[1], 3), ", ", round(cutobj$bound80[2], 3), ")", sep=""))
    message(paste("90% Certainty Range: (", round(cutobj$bound90[1], 3), ", ", round(cutobj$bound90[2], 3), ")", sep=""))
  }
  
  return(model)
}


###########################################################################################
###########################################################################################
################################  START HERE   ############################################
###########################################################################################

#-----------------------------------------------
#[II] Load and clean data
#-----------------------------------------------
#-----------------------------------------------
#[IIA] Import data 
#
#!!!Note: the full path pointing to you data must replace the text 
#         "C:/myfolders/mydata/ExampleData.csv" below with the exact file location.
#         R does not recognize backslashes ("\") so be sure to use "/" instead.
#         For this code to work your data must be a csv file
#
#-----------------------------------------------
input<-read.table("/Users/sarahsullivan/Google Drive/Sarah Sullivan/Data/Working/ExampleData.csv", sep = ",", quote = "", header=TRUE)


#-----------------------------------------------
#[IIB] Create dataframe with two variables, id (must be unique identifier) and data
#
#!!!Note: replace the names "identified" and "value" with the variable names
#         used in your data to indicate the unique ID for each observation, and
#         the data (e.g. OD value) on which you wish to base your cutoff determination.
#
#-----------------------------------------------

examp<-data.frame(id=input$identifier, data=input$value)



#-----------------------------------------------
#[III] Determine which type and number of distributions is optimal
#-----------------------------------------------
#-----------------------------------------------
#[IIIA] Check data for distributional assumptions and extreme outliers
#
#!!!Note: the user may modify the titles by replacing the text following "main = " 
#         to make them more descriptive of the data
#
#-----------------------------------------------
hist(examp$data, breaks=100, main="Histogram of Data")
dens<-density(examp$data)
plot(dens, main="Density of Data")

#-----------------------------------------------
#[IIIB] Run fitting function and determine which distribution and number of components is best by BIC
#
#!!!#Note: the fitloops function may take a decent amount of time. 
#     Pay attention to the error messages to decide what to do if it fails 
#     (With the example data it should run eventually, though sometimes bad seed values 
#     mean you may have to try it a couple times)
#-----------------------------------------------

fitloopsobject<-fitloops(datawithids = examp)

bicgraph(fitobj = fitloopsobject)

#-----------------------------------------------
#[IIIC] Pick the best combination of distribution and components (usually the one that is optimal by BIC unless there is scientific rationale)
#
#     To default to the optimal fit (recommended!) run the code as it is written below.
#     If you would like to specify a distribution or number of sub-populations the function call line is:
#
#              modelpick(fitobj = , dist = "", ncomp = NA){
#
#                 dist: either "Normal" or "Skew-normal"
#                 ncomp: number from 1 to 5
#-----------------------------------------------

modelpickobject<-modelpick(fitobj = fitloopsobject)


#-----------------------------------------------
#[IV]  Determine which cutpoint is best if more than two components are optimal
#-----------------------------------------------
#-----------------------------------------------
#[IVA] Investigate proportions of positivity and compare to scientific context
#-----------------------------------------------
#-----------------------------------------------
#[IVB] Utilize visualizations of uncertainty, distributions and cutpoint to aid in decision-making
#-----------------------------------------------
rawuncertgraph(modelpickobj = modelpickobject)
rawdistgraph(modelpickobj = modelpickobject)
rawhistcuts(modelpickobj = modelpickobject)

#-----------------------------------------------
#[V]   Create the cutpoint and investigate what resulting distributions look like
#-----------------------------------------------
#-----------------------------------------------
#[VA]  Implement Cutpoint with Standard and/or Non-standard Certainty Levels for Indeterminate Ranges
#
# !!! If there are more than 2 sub-populations identified in the chosen model then the user will need to 
#     tell the program between which two sub-populations to set the cutoff.  
#
#     The function call line is:
#
#           cutoff<-function(modelpickobj = , cutcomp = 0, standardcert = T, newcertlevel = NA){
#
#     cutcomp = the lower of the two distributions between which you want to set the cutoff. 
#               So, for example, set the cutoff between 1 and 2, let cutcomp = 1. To set the 
#               To set the cutoff between distributions 4 and 5, let cutcomp = 4.
#     standardcert = T to generate the standard indeterminate ranges (80 % and 90% certainty).
#                  = F to repress generating the standard indeterminate ranges (80 % and 90% certainty).
#     newcertlevel = the certainty level of a new indeterminate range which you want to create. It should lie on the range 0.5<x<1
#
#-----------------------------------------------
cutoffobject<-cutoff(modelpickobj = modelpickobject, cutcomp = 1, standardcert = T, newcertlevel = 0.95)

#-----------------------------------------------
#[VB]  Investigate resulting distributions and classifications
#-----------------------------------------------
cutuncertgraph(cutobj = cutoffobject)
cutdistgraph(cutobj = cutoffobject)

summaryobject<-summaryout(cutobj = cutoffobject)
outtable<-summaryobject$outtab
outtable




