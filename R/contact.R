#' Simulate contact tracing data from a simOutbreak object
#'
#' This function takes a simOutbreak object and returns a dataframe of contact pairs. The probabilites of reporting infectious and non-infectious contacts are provided as parameters.
#' 
#' @importFrom magrittr "%>%"
#' 
#' @param outbreak a simOutbreak object
#'
#' @param eps the reporting probability of infectious contacts
#' 
#' @param xi the scaling factor relative to eps, describing the probability of reporting contact between non-infectiously related pairs
#'
#' @param plot a logical indicating if the contact tracing data should be plotted
#'
#' @author Finlay Campbell (\email{f.campbell15@@imperial.ac.uk})
#' 
#' @export

simCTD <- function(outbreak,eps=1,xi=0,plot=FALSE){
  
  if(outbreak$n==1) return("No transmission observed")
  
  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
  
  accept.reject <- function(pair){
    if(pair[3]) return(stats::runif(1,0,1) < eps)
    else return(stats::runif(1,0,1) < xi*eps)
  }
  
  import <- which(is.na(outbreak$ances))
  
  infec.contact <- cbind(outbreak$ances[-import],outbreak$id[-import])
  
  potent.CTD <- as.data.frame(t(utils::combn(outbreak$id,2)))
  colnames(potent.CTD) = c("i","j")
  
  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))
  
  CTD <- potent.CTD[potent.CTD$accept,1:2]
  rownames(CTD) <- NULL
  
  if(plot){
    plot.CTD <- potent.CTD[potent.CTD$accept | potent.CTD$contact,]
    
    nodes.color <- grDevices::gray(1/3+2*outbreak$onset/(3*tail(outbreak$onset,1)))
    nodes.color[import] <- "#e60000"
    
    nodes.label = outbreak$id
    
    edges.color <- rep("#e60000",nrow(plot.CTD))
    edges.color[!plot.CTD$contact] <- "blue"
    
    edges.arrows <- rep("to",nrow(plot.CTD))
    edges.arrows[!plot.CTD$contact] <- FALSE
    
    edges.dashes <- plot.CTD$contact & !plot.CTD$accept
    
    edges.width <- rep(3,nrow(plot.CTD))
    edges.width[!plot.CTD$accept]  <- 0.5
    
    nodes <- data.frame(id=outbreak$id,label=nodes.label,color=nodes.color)
    edges <- data.frame(from=plot.CTD$i,to=plot.CTD$j,dashes=edges.dashes,arrows=edges.arrows,
                        color=edges.color,width=edges.width)
    
    network <- magrittr::`%>%`(visNetwork::visNetwork(nodes,edges,main="Simulated Contact Network"),
                               visNetwork::visNodes(shape="ellipse",font=list("size"=25),borderWidth=2))
    
    print(network)
  }
  
  return(CTD)
}
