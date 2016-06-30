

#################################################################################################
######### Finlay's ongoing contribution to Outbreaker2, working on incoporating contact ######### 
######### tracing data into the inference of posterior transmission trees.              ######### 
#################################################################################################


##################################################
######### LOADING LIBRARIES AND DATASETS ######### 
##################################################

library(visNetwork)
library(reshape2)


######################################
######### DEFINING FUNCTIONS #########
######################################

#This function analyses the accuracy and precision of outbreaker inference 
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
    nodes <- data.frame(id=id,label=id)
    plot.network <- visNetwork(nodes,edges=network) %>%
      visNodes(shape="ellipse",font=list("size"=25),borderWidth=2)
    print(plot.network)
  }
  
  #Determine confidence in our results
  id <- id[-unlist(lapply(id,function(i) any(is.na(result[[i+adder]]))))]
  confidence <- round(mean(unlist(lapply(id,function(i) sum(table(result[[i+adder]])^2)/samples^2))),2)
  
  #Determine the proportion of correctly inferred ancestries
  num.correct <-  sum(true.outbreak$ances==network$from,na.rm=TRUE)
  num.correct <- num.correct + sum(is.na(true.outbreak$ances[is.na(network$from)]))
  prop.correct <- round(num.correct/nrow(network),2)
  
  if(print) print(paste(100*prop.correct,"% of transmission tree correctly inferred | Average confidence: ",confidence,sep=""))
  
  out <- list(transmission=network,confidence=confidence,prop.correct=prop.correct)
  return(out)
}

#A function to simulate contact tracing data (CTD)
simCTD <- function(temp.outbreak,eps=1,chi=1,ksi=0,plot=FALSE,print.ratio=FALSE){
  
  if(temp.outbreak$n==1) return("No transmission observed")
  
  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
  
  accept.reject <- function(pair){
    if(pair[3]) return(runif(1,0,1) < chi*eps)
    else return(runif(1,0,1) < ksi*chi*eps)
  }
  
  import <- which(is.na(temp.outbreak$ances))
  infec.contact <- cbind(temp.outbreak$ances[-import],temp.outbreak$id[-import])
  
  potent.CTD <- as.data.frame(t(combn(temp.outbreak$id,2)))
  colnames(potent.CTD) = c("i","j")
  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))
  
  CTD <- potent.CTD[potent.CTD$accept,1:2]
  rownames(CTD) <- NULL
  
  if(print.ratio) print(paste("un/in:",round(sum(!CTD$contact)/sum(CTD$contact),2)))
  
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

###############################################
######### SIMULATING OUTBREAK AND CTD #########
###############################################

temp.w.dens <- dgamma(1:20,2,0.05)
true.outbreak <- outbreaker::simOutbreak(2,temp.w.dens,n.hosts=40,mu.transi=1e-05,rate.import.case=0)
sim.contact<- simCTD(true.outbreak,plot=TRUE,ksi=0.01,eps=0.8)

result <- outbreaker(data=list(dates=true.outbreak$onset,dna=true.outbreak$dna,w.dens=temp.w.dens,CTD=sim.contact))
analysis <- result.analysis(result,true.outbreak,print=TRUE)

runs <- 50
out2 <- matrix(NA,nrow=runs,ncol=4)
colnames(out2) <- c("accuracy","confidence","CTD.accuracy","CTD.confidence")
counter = 1

while(counter<=runs){
  true.outbreak <- outbreaker::simOutbreak(2,temp.w.dens,n.hosts=40,mu.transi=1e-05,rate.import.case=0)
  if(true.outbreak$n<15) next
  
  print(counter)
  
  sim.contact<- simCTD(true.outbreak,plot=TRUE,ksi=0.01,eps=0.8)
  
  their.result <- outbreaker(data=list(dates=true.outbreak$onset,w.dens=temp.w.dens,dna=true.outbreak$dna))
  analysis <- result.analysis(their.result,true.outbreak,print=TRUE)
  out2[counter,1:2] <- c(analysis$prop.correct,analysis$confidence)
  
  our.result <- CTD.outbreaker(data=list(dates=true.outbreak$onset,dna=true.outbreak$dna,w.dens=temp.w.dens,CTD=sim.contact))
  analysis <- result.analysis(our.result,true.outbreak,print=TRUE)
  out2[counter,3:4] <- c(analysis$prop.correct,analysis$confidence)
  
  counter <- counter + 1
}

out2 <- out2[,c(1,3,2,4)]
out2 <- as.data.frame(melt(out2))
out2$model <- c(rep("original",50),rep("CTD",50),rep("original",50),rep("CTD",50))

p <- ggplot(out2) + geom_violin(aes(x=factor(Var2),y=value,fill=factor(model))) + ylim(0.35,1) +
  ggtitle("80% coverage | 0.01 false positive rate") + xlab("") + ylab("Value") + guides(fill=FALSE) +
  theme_set(theme_gray(base_size = 18))
p

ggsave("imperfect.CTD.png",p)