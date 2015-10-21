\name{outbreaker graphics}
\alias{plotChains}
\alias{transGraph}
\alias{plotOutbreak}
\title{Plot outbreaker's results}
\description{
  These are the main functions used for generating graphics from the raw
  output of \code{outbreaker} and \code{outbreaker.parallel}.

  \itemize{
    \item \code{plotChains} is used for plotting MCMCs

    \item \code{transGraph} plots a graph of inferred ancestries

    \item \code{plotOutbreak} attempts to synthetize the reconstruction
    of small outbreaks
  }
}
\usage{
plotChains(x, what="post", type=c("series","density"), burnin=0,
           dens.all=TRUE, col=funky(x$n.runs), lty=1, lwd=1,
           main=what, legend=TRUE, posi="bottomleft", \dots)

transGraph(x, labels=NULL, burnin=x$burnin, threshold=0.2, col.pal=NULL,
           curved.edges=TRUE, annot=c("dist","support"), sep="/", \dots)

plotOutbreak(x, burnin=x$burnin, thres.hide=0.2, col=NULL,
             col.pal=colorRampPalette(c("blue","lightgrey")),
             edge.col.pal=NULL, col.edge.by="prob",
             annot=c("dist","prob"), sep="/", cex.bubble=1,
             edge.max.dist = 10, lwd.arrow=2, xlim=NULL, \dots)
}
\arguments{
  \item{x}{the output of \code{outbreaker} or
    \code{outbreaker.parallel}.}
  
  \item{what}{a character chains giving the name of the item to be
    plotted. See \code{names(x$chains)} for possible values. By default,
    log-posterior values are plotted}
  
  \item{type}{a character indicating if the chains should be plotted as
    time series ("series"), or as density ("density").}

  \item{burnin}{an integer indicating the number of MCMC steps to
    discard before plotting chains.}

  \item{dens.all}{a logical indicating if, in the case of multiple runs,
    the overall density of the different chains should be plotted in
    addition to individual densities.}

  \item{col}{a vector of colors to be used to plot different chains.}

  \item{lty}{a vector of integers specifying line types for the
    different chains.}

  \item{lwd}{same as \code{lty}, but for line width.}

  \item{main}{the title to be added to the plot.}

  \item{labels}{the labels to be used to name the nodes of the graph (cases).}

  \item{threshold}{the minimum support for ancestries to be plotted;
    'support' is defined as the frequency of a given ancestor in the
    posterior distribution; defaults to 0.2.}

  \item{thres.hide}{a threshold of posterior support for displaying
    ancestries; ancestries with less than this frequency in the
    posterior are hidden.}

  \item{col.pal,edge.col.pal}{the color palette to be used for the edges
    (ancestries).}

  \item{curved.edges}{a logical indicating whether edges should be
    curved.}

  \item{col.edge.by}{a character string indicating which information
    should be used to color the edges ('dist': genetic distance; 'prob':
    support for the ancestry)}

  \item{annot}{a character indicating which information should be used
    to annotate the edges; this can be the distances between ancestors
    and descendents ("dist") and the posterior support for ancestries
    ("support"); if both are requested, fields will be concatenated.}

  \item{sep}{a character indicating the separator to be used when
    concatenating several types of annotation.}

  % \item{arr.col.pal}{the color palette to be used for the arrows,
  %   representing the support for each ancestry; if NULL, a grey scale
  %   is used, with larger support in darker shades.}

  \item{cex.bubble}{a numeric value indicating the size factor for the
    bubbles representing the generation time distribution.}

  \item{edge.max.dist}{a number indicating the threshold distance
    bounding the color palette used for the edges; useful to avoid
    showing edges corresponding to distances larger than a given
    number. }
  
  \item{lwd.arrow}{a numeric value indicating the size factor for the
    arrows.}

  \item{xlim}{the limits of the X axis; if NULL, determined from the
    data.}

  \item{legend}{a logical indicating if a legend should be plotted for
    the different runs.}

  \item{posi}{a character string indicating the position of the legend
    (see \code{?legend}).}

  \item{\dots}{further arguments to be passed to other functions.}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\examples{
data(fakeOutbreak)
attach(fakeOutbreak)

## examine MCMC
plotChains(res)
plotChains(res,type="dens")
plotChains(res,type="dens", what="mu1", burnin=2e4)

## represent posterior ancestries
transGraph(res, annot="", main="Posterior ancestries")
transGraph(res, annot="", main="Posterior ancestries - support > 0.5",
   threshold=0.5)
if(require(adegenet)){
transGraph(res, annot="", main="Posterior ancestries - support > 0.01",
   threshold=0.01, col.pal=spectral)
}
## summary plot
plotOutbreak(res,cex.bubble=0.5, thres.hide=0.5,
   main="Outbreak reconstruction")


detach(fakeOutbreak)

}
