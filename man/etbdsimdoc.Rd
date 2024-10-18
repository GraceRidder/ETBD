\name{simETBD}
\alias{simETBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generates ETBD simulation 
%%  ~~function to do ... ~~
}
\description{ Simulates population size dependnet speciation or extinction dynamics in species richness and phylogenetic trees.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simETBD( t = NULL, DIST = "", JmaxV = NULL, NegExpEx = F, exparm = NULL,  psymp = NULL, ExpSp = T, ExpSpParm = NULL,  constantEX = NULL, splitparm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{Number of timesteps.
%%     ~~Describe \code{x} here~~
}
  \item{DIST}{A string either "SRS", "NORM", "GEO", "NO", describing the SAD.
%%     ~~Describe \code{x} here~~
}
  \item{JmaxV}{The level of resources.
%%     ~~Describe \code{x} here~~
}
  \item{NegExpEx}{If True turns on population size dependent extinction.
%%     ~~Describe \code{x} here~~
}
  \item{exparm}{Paramater of population size dependent extinction.
%%     ~~Describe \code{x} here~~
}
  \item{psymp}{Paramater of constant speciation.
%%     ~~Describe \code{x} here~~
}
  \item{ExpSp}{If True turns on population size dependent speciaiton.
%%     ~~Describe \code{x} here~~
}
  \item{constantEX}{Paramater of constant extinction.
%%     ~~Describe \code{x} here~~
}
  \item{splitparm}{Paramater of speciation symmetry.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

### runs a simulation for population size dependent speciation for 200 time steps.
res1 = simETBD(
  t =200,
  DIST = "SRS",     
  JmaxV = 2000,
  NegExpEx = F,  
  exparm = -.7,  
  psymp = .35,    
  ExpSp = T,     
  ExpSpParm = 1.6,  
  constantEX = .005, 
  splitparm = .5, 
)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
