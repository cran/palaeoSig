\name{NEWS}
\title{palaeoSig News}
\encoding{UTF-8}

\section{Changes in version 1.1}{

  \subsection{GENERAL}{
    \itemize{
      
      \item The package \pkg{autocorTF} has been merged with \pkg{palaeoSig} for ease of maintenance.
      
    }
  } % end GENERAL
  
  \subsection{NEW FUNCTIONS}{
    \itemize{      
      \item \code{identify.obscor} allows species names to be added to plots interactively.
      
      
    }
  } % end NEW FUNCTIONS

  \subsection{MODIFIED FUNCTIONS}{
    \itemize{      
      \item \code{obs.cor} Now includes correlations with several different species weights
      \item \code{plot.obscor} upgraded to allow the a choice of which abundance weighting is used, and the code for scaling points has been improved.
      
      
    }
  } % end NEW FUNCTIONS

  \subsection{BUG FIXES}{
    \itemize{      
      \item \code{randomTF} bug when partialling out reconstruction other than MAT now fixed.
    }
  } % end NEW FUNCTIONS
  
  

}% end VERSTION 1.1-1




\section{Changes in version 1.1}{

  \subsection{GENERAL}{
    \itemize{
      
      \item The package \pkg{autocorTF} has been merged with \pkg{palaeoSig} for ease of maintenance.
      
    }
  } % end GENERAL
  
  \subsection{NEW and REWRITTEN FUNCTIONS}{
    \itemize{
      
      \item \code{RNE} replaces function \code{mat.rne}, \code{mat.he} and \code{mat.rd}. \code{RNE} allows any of the transfer function methods in \pkg{rioja} to be used to find the dependence of transfer function performance on spatially close observations.
      \item \code{plot.RNE} replaces \code{plot.rne} to work with the output of \code{RNE}.
      \item \code{obs.cor} replaces the previous version of \code{obs.cor}, \code{sim.cor} and \code{obscor.sig} to make the function easier to use and more similar to \code{randomTF}
      \item \code{plot.obscor} replaces the previous version of \code{plot.obscor} and \code{plot.simcor} to use the new output of \code{obs.cor}. The two plot types can be selected with the argument \code{which}.
      \item \code{jointsig}: tests if two environmental variables have joint control on fossil assemblage composition. 
      \item \code{plot.js}: a \code{plot} function for \code{jointsig}
      \item \code{Hill.N2.core} now calculates the Minimum, first quartile and median effective number of species for all fossil observations. This makes much more sense than the previous version, which conflated diversity in individual levels and turnover between levels.
      \item \code{plot.palaeoSig} has been improved to give neater figures.
      
    }
  } % end NEW FUNCTIONS
  
  \subsection{REMOVED FUNCTIONS}{
    
    \itemize{
      
      \item \code{wajack} was removed as \code{WA} is now implemented in \pkg{rioja}. 
      
    }
  } % end REMOVED FUNCTIONS
  
  \subsection{FORTHCOMMING CHANGES}{
    \itemize{
      
      \item \code{mat.h} will be removed in the next release, when code for h-block resampling should be implimented in \pkg{rioja}.
      
    }
  } % end FORTHCOMMING CHANGES
}% end VERSTION 1.1
   	  