

#################################################################################################
######### Finlay's ongoing contribution to Outbreaker2, working on incoporating contact ######### 
######### tracing data into the inference of posterior transmission trees.              ######### 
#################################################################################################

#for(i in list.files("R",".R")) source(paste("R\\",i,sep=""))

##################################################
######### LOADING LIBRARIES AND DATASETS ######### 
##################################################

#library(visNetwork)
#library(reshape2)

######################################
######### DEFINING FUNCTIONS #########
######################################

## A function to analyse the accuracy and precision of outbreaker inference 
result.analysis <- function(result,true.outbreak,plot=TRUE,print=FALSE){
  id <- seq_len(true.outbreak$n)
  adder <- which(names(result)=="alpha.1")-1
  samples <- length(result$step)
  
  #Determine the modal transmission network
  network <- data.frame(from=do.call(rbind,lapply(id, function(i) ({
    modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
    if(length(modal.ances)==0) return(NA) else return(modal.ances)
  }))), to=id) 
  
  
  #Plot proposed transmission network
  if(plot){
    
    import <- which(is.na(network$from))
    
    #Define the indices to call the times of infection
    t.inf.index <- which(names(result)=="t.inf.1")-1+id
    
    #Determine the median posterior time of onset for plotting
    onset <- unlist(lapply(result[t.inf.index],median))
    
    #Scale onset to begin at 0
    onset <- onset - min(onset)
    
    nodes <- data.frame(id=id,label=id)
    nodes$color <- gray(1/3+2*onset/(3*tail(onset,1)))
    nodes$color[import] <- "#e60000"
    
    plot.network <- visNetwork(nodes,network) %>%
      visNodes(shape="ellipse",font=list("size"=25),borderWidth=2) %>%
      visEdges(arrows="to")
    print(plot.network)
  }
  
  #Determine confidence in our results

  transmission.id <- id[-unlist(lapply(result[id+adder],function(i) any(is.na(i))))]
  
  square.conf <- round(mean(unlist(lapply(result[transmission.id+adder],function(i) sum(table(i)^2)/samples^2))),2)
  
  mean.conf <- round(mean(unlist(lapply(transmission.id,function(i) mean(result[[i+adder]]==true.outbreak$ances[i])))),2)
  
  entropy <- round(mean(unlist(lapply(result[transmission.id+adder],function(i) {fk <- table(i)/sum(table(i))
                                                                                 -sum(log(fk)*fk)}))),2)
   
  #Determine the proportion of correctly inferred ancestries
  num.correct <-  sum(true.outbreak$ances==network$from,na.rm=TRUE)
  num.correct <- num.correct + sum(is.na(true.outbreak$ances[is.na(network$from)]))
  prop.correct <- round(num.correct/nrow(network),2)
  
  if(print) print(paste(100*prop.correct,"% of transmission tree correctly inferred | Average confidence: ",mean.conf,sep=""))
  
  out <- list(transmission=network,square.conf=square.conf,mean.conf=mean.conf,
              prop.correct=prop.correct,entropy=entropy)
  return(out)
}

## A function to simulate contact tracing data (CTD)
simCTD <- function(temp.outbreak,eps=1,ksi=0,plot=FALSE,print.ratio=FALSE){
  
  if(temp.outbreak$n==1) return("No transmission observed")
  
  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
  
  accept.reject <- function(pair){
    if(pair[3]) return(runif(1,0,1) < eps)
    else return(runif(1,0,1) < ksi*eps)
  }
  
  import <- which(is.na(temp.outbreak$ances))
  infec.contact <- cbind(temp.outbreak$ances[-import],temp.outbreak$id[-import])
  
  potent.CTD <- as.data.frame(t(combn(temp.outbreak$id,2)))
  colnames(potent.CTD) = c("i","j")
  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))
  
  CTD <- potent.CTD[potent.CTD$accept,1:2]
  rownames(CTD) <- NULL

  if(plot){
    plot.CTD <- potent.CTD[potent.CTD$accept | potent.CTD$contact,]
    
    nodes.color <- gray(1/3+2*temp.outbreak$onset/(3*tail(temp.outbreak$onset,1)))
    nodes.color[import] <- "#e60000"
    
    nodes.label = temp.outbreak$id
    
    edges.color <- rep("#e60000",nrow(plot.CTD))
    edges.color[!plot.CTD$contact] <- "blue"
    edges.arrows <- rep("to",nrow(plot.CTD))
    edges.dashes <- plot.CTD$contact & !plot.CTD$accept
    edges.width <- rep(3,nrow(plot.CTD))
    edges.width[!plot.CTD$accept]  <- 0.5
    
    nodes <- data.frame(id=temp.outbreak$id,label=nodes.label,color=nodes.color)
    edges <- data.frame(from=plot.CTD$i,to=plot.CTD$j,dashes=edges.dashes,arrows=edges.arrows,
                        color=edges.color,width=edges.width)
    network <- visNetwork(nodes,edges,main="Simulated Contact Network") %>%
      visNodes(shape="ellipse",font=list("size"=25),borderWidth=2)
    print(network)
  }
  
  return(CTD)
}

## A function to compare accuracy and precision of outbreaker and CTD.outbreaker
compare.outbreakers <- function(runs=10,min.hosts=10,max.hosts=15){
  library(ggplot2)
  
  temp.w.dens <- dgamma(1:20,2,0.05)

  out <- matrix(NA,nrow=runs,ncol=8)
  colnames(out) <- c("accuracy","CTD.accuracy","entropy","CTD.entropy","mean.conf","CTD.mean.conf","time","CTD.time")
  counter <- 1
  
  while(counter<=runs){
    true.outbreak <- outbreaker::simOutbreak(2,temp.w.dens,n.hosts=15,mu.transi=1e-05,rate.import.case=0)
    if(true.outbreak$n<min.hosts) next
    
    print(counter)
    
    sim.contact <- simCTD(true.outbreak,plot=TRUE,ksi=0.01,eps=0.8)
    
    temp.time <- system.time(their.result <- outbreaker(data=list(dates=true.outbreak$onset,w.dens=temp.w.dens,
                                             dna=true.outbreak$dna),config=list(n.iter=2e4, sample.every=200)))
    analysis <- result.analysis(their.result,true.outbreak,print=TRUE,plot=FALSE)
    out[counter,c(1,3,5,7)] <- c(analysis$prop.correct,analysis$entropy,analysis$mean.conf,temp.time[1])
    
    temp.time <- system.time(our.result <- CTD.outbreaker(data=list(dates=true.outbreak$onset,
                                           dna=true.outbreak$dna,w.dens=temp.w.dens,CTD=sim.contact),
                                           config=list(n.iter=2e4, sample.every=200)))
    analysis <- result.analysis(our.result,true.outbreak,print=TRUE,plot=FALSE)
    out[counter,c(2,4,6,8)] <- c(analysis$prop.correct,analysis$entropy,analysis$mean.conf,temp.time[1])
    
    counter <- counter + 1
    
    print(paste("~",round(mean(c(out[,7],out[,8]),na.rm=TRUE)*(runs-counter)*2/60,2)," minutes left",sep=""))
  }
  
  plot.out <- out
  plot.out[,7:8] <- out[,7:8]/max(out[,7:8])
  plot.out <- as.data.frame(melt(plot.out))
  
  plot.out$variable <- c(rep("Accuracy",2*runs),rep("Entropy",2*runs),rep("Confidence",2*runs),rep("Time",2*runs))
  plot.out$Model <- factor(rep(c(rep("Original",runs),rep("CTD",runs)),4),levels=c("Original","CTD"))
    
  p <- ggplot(plot.out) + geom_violin(aes(x=variable,y=value,fill=Model)) + ylim(0.,1) +
       ggtitle("80% coverage | 0.01 false positive rate") + ylab("Value") +
       theme_set(theme_gray(base_size = 18)) + theme(axis.title.x=element_blank())
  
  print(p)
  
  ggsave("imperfect.CTD.png",p)
  
  return(list("analysis"=as.data.frame(out),"plot"=p))
}