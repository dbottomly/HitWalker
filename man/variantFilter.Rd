\name{variantFilter-class}
\Rdversion{1.1}
\docType{class}
\alias{variantFilter-class}
\alias{variantFilter}
\alias{filterBestCategory,variantFilter-method}
\alias{filterBestCategory}
\alias{filterCategory,variantFilter-method}
\alias{filterCategory}
\alias{filterCollapseCol,variantFilter-method}
\alias{filterCollapseCol}
\alias{filterDta,variantFilter-method}
\alias{filterDta}
\alias{filterIgnoreCol,variantFilter-method}
\alias{filterIgnoreCol}
\alias{filterLabels,variantFilter-method}
\alias{filterLabels}
\alias{filterSep,variantFilter-method}
\alias{filterSep}
\alias{show,variantFilter-method}
\alias{makeVariantFilter}

\title{Class \code{"variantFilter"}}
\description{
This class provides a way of working with text-based variant filters such as those provided in the Genome Analysis Toolkit.  Most importantly, users can specify labels for the filters and the type of variant to which they should be applied.
}
\section{Objects from the Class}{
Objects should be created from a call to the function \code{makeVariantFilter}.  See the example section.

}
\section{Slots}{
  \describe{
    \item{\code{filter.dta}:}{Object of class \code{"data.frame"} containing three columns: 'category', 'label' and 'type', all of \code{character} mode.
        'Category' refers to how the end user should view the filter value (e.g. PASS or FAIL).  The \code{data.frame} should be sorted based on decreasing category value.
        'Label' refers to the actual filters, note multiple filters applied to a variant will be parsed using \code{filter.sep}.  All combinations of the filters will be considered internally so only one combination needs to be specified.
        'Type' refers to the type of variant that the filter(s) supplied in label should be applied to (e.g. INDEL or SNV).  Again can refer to multiple variants using \code{filter.sep}.
}
    \item{\code{filter.sep}:}{Object of class \code{"character"} delimiting multiple labels or types.}
    \item{\code{filter.collapse.col}:}{Object of class \code{"character"} containing the name of the column from the retrieved variant data that the values in label would apply to.}
    \item{\code{filter.ignore.col}:}{Object of class \code{"character"} indicating redudant columns with \code{filter.collapse.col}.  This could be an auto-incremeted integer value produced by a DBMS.}
    \item{\code{filter.best.category}:}{Object of class \code{"character"} containing the 'category' that should be considered the best quality.}
  }
}
\section{Methods}{
  \describe{
    \item{filterBestCategory}{\code{signature(obj = "variantFilter")}: A method for retrieving the best category} 
    \item{filterCategory}{\code{signature(obj = "variantFilter")}: A method for retrieving a \code{character} vector of categories named by label}
    \item{filterCollapseCol}{\code{signature(obj = "variantFilter")}: A method for retrieving the column of the variant \code{data.frame} that contains the variant filters}
    \item{filterDta}{\code{signature(obj = "variantFilter")}: A method to retrieve the filter \code{data.frame}}
    \item{filterIgnoreCol}{\code{signature(obj = "variantFilter")}: A method to retrieve the columns of the variant \code{data.frame} that should be ignored}
    \item{filterLabels}{\code{signature(obj = "variantFilter")}: A method to retrieve a \code{character} vector containing filter types named by filter labels}
    \item{filterSep}{\code{signature(obj = "variantFilter")}: A method for retrieving the character used to seperate the label and type values }
	 }
}
\references{
Genome Analysis Toolkit Website: http://www.broadinstitute.org/gatk/
}
\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\link{get.sample.variants}}, \code{\link{filter.variant.annots}}
}
\examples{
category <- c("OK", "OK", "WARN", "WARN", "WARN", "WARN")
label <- c("PASS", "StrandBiasFilter", "GATKStandardFilter", "HARD_TO_VALIDATE", "HARD_TO_VALIDATE;GATKStandardFilter", "GATKStandardFilter;StrandBiasFilter")
type <- c("SNV;INDEL", "SNV;INDEL", "SNV;INDEL", "SNV;INDEL", "SNV;INDEL", "INDEL")
    
filter.dta <- data.frame(category, label, type, stringsAsFactors=FALSE)
var.filter <- makeVariantFilter(filter.dta, filter.sep=";", filter.collapse.col="filter", filter.ignore.col="variant_qual_id", filter.best.category="OK")
}
\keyword{classes}
