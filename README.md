Notice
=========

HitWalker2 now supersedes this package and is freely available at https://github.com/biodev/HitWalker2.  Setup is simplified for common datatypes though the implementation has switched to python/javascript as opposed to R.  For those interested, a streamlined R version of the prioritization code (no visualizations though) is available as part of the `hwhelper` package bundled with HitWalker2.

Installation
==========

Multiple packages are either required or suggested for HitWalker to be installed and the main ones can be retrieved from the following two sources:

From inside R use the following commands:

The Comprehensive R Archive Network  using:  
install.packages(c("igraph", "RSQLite", "cgdsr", "reshape2"))

Bioconductor using:  
source("http://bioconductor.org/biocLite.R")  
biocLite(c("graph", "Rgraphviz", "Streamer", "biomaRt"))

Note that of these packages Rgraphviz is the most difficult to install due to its dependence on the graphviz external library (though I believe this is 
may no longer be an issue in current and future Bioconductor releases). 

Please see the documentation at: http://www.bioconductor.org/packages/release/bioc/readmes/Rgraphviz/README for installation instructions for your platform.

When these dependencies are met, the release version of HitWalker can be downloaded and installed from within R using the following commands:

Linux, Mac OSX or Windows (assuming you are in the appropriate directory):
install.packages("HitWalker_0.99.1.tar.gz", repos=NULL, type="source")

Also, the accompanying HitWalkerData package should be installed as well.

install.packages("HitWalkerData_0.99.0.tar.gz", repos=NULL, type="source")

Note that HitWalker has currently only been tested on R versions R-2.15.1 R-3.0.0 using a 64-bit Ubuntu Linux machine.

Once installed, the package can be loaded and the user can familiarize themselves with its use by reading any of the four supplied vignettes.  The vignettes build on
each other and so the suggested reading order is:

Hitwalker_Create_DB  
Hitwalker_Add_metadata  
Hitwalker_glio  
HitWalker  

First load the package:

library(HitWalker)  
library(HitWalkerData)

Vignettes can be examined through a command similar to:

vignette("Hitwalker_Create_DB")

Similarly, all the vignettes can be listed by typing:

vignette(package="HitWalker")

To generate an executable file called 'HitWalker.R' containing all the code in the vignette carry out a command similar to:

Stangle(system.file("doc/HitWalker.Rnw", package="HitWalker"))

Where 'HitWalker.Rnw' can be any of the supplied vignettes. 

From the Stangle command a file called 'HitWalker.R' will be produced and can be executed section by section or all at once by:

source("HitWalker.R")

See the R documentation at http://cran.r-project.org/ for general use of the R programming language.

Deployment
===========

HitWalker requires a database backend and therefore a fair amount of work is required before use of the R package.  

Contributors
============

Several extensions/modifications of HitWalker are in the works (and will be posted on Github soon).  Though this code
base will continue to be maintained, anyone interested in new features is encouraged to contribute them and to keep me 
in the loop.

Contact
=============

Dan Bottomly
bottomly@ohsu.edu
