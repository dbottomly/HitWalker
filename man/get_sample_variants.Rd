\name{get.sample.variants}
\alias{get.sample.variants}
\alias{filter.variant.annots}
\title{
Retrieving and filtering variants
}
\description{
These functions provide the mechanism for retrieval and processing of variants. The assumption is that the data retrieved has essentially been limited to 
variants of interest.  Further filtering of text labels from the genotyping procedure is provided by information specified in the \code{priorDbParams} class.
}
\usage{
get.sample.variants(sample.id.list, db.con, param.obj)
filter.variant.annots(var.dta, filter.list, var.id.col, filter.sep, collapse.col="filter", ignore.cols="variant_qual_id")
}
\arguments{
  \item{sample.id.list}{
	A list with three elements: 'recon' is the human readable name, 'ind' is the database index value and 'num.pats' is the number of samples in the database.
}
  \item{db.con}{
	An object inherited from \code{DBIConnection} which represents the datasource.
}
  \item{param.obj}{
	A \code{priorDbParams} object.
}
  \item{var.dta}{
	A \code{data.frame} containing the retrieved variant data.

}
  \item{filter.list}{
	A list containing the names of the filters to be applied.  Each list element should contain a vector with filter labels.
}
  \item{var.id.col}{
	The column of var.dta uniquely identifying variants (and potentially alleles) for a given sample.
}

  \item{filter.sep}{
	The character used to seperate filter values as specified in \code{variantFilter}.
}

  \item{collapse.col}{
	The column to be collapsed if necessary, (e.g. the filter column).
}
  \item{ignore.cols}{
	The columns redundant with the collapse.col that should be removed as well (e.g. an autoincremented DB index).  NULL if none are present.
}
}

\value{
\describe{
\item{get.sample.variants}{returns a list containing two elements:}
\item{var.dta}{The filtered \code{data.frame} containing the variant results.}
\item{summary.id}{The core or summary ID used for this function}

\item{filter.variant.annots}{returns a \code{data.frame} containing the remaining elements of \code{var.dta} after filtering}
}
}

\author{
Daniel Bottomly
}
\note{
\code{filter.variant.annots} will likely never be called directly by the user, instead it will be called according to the \code{variantFilter} class by \code{get.sample.variants}.
}

\seealso{
\code{\linkS4class{priorDbParams}}, \code{\linkS4class{variantFilter}}
}
\examples{
if (require(RSQLite) && require(HitWalkerData))
{
        set.seed(123)
	data(params)
        db.con <- dbConnect("SQLite", hitwalker.db.path())
        graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
        sample.id.list <- reconcile.sample.name("08-00102", db.con, examp.prior.param)
	raw.vars <- dbGetQuery(db.con, HitWalker:::testVariantQuery(sample.id.list, examp.prior.param))
	filter.values <- getFilterValues(examp.prior.param)		
	filt.list <- strsplit(names(filter.values)[filter.values == "SNV;INDEL"], getFilterSep(examp.prior.param))
	pat.var.id <- patVarID(examp.prior.param)
	filt.vars <- filter.variant.annots(var.dta=raw.vars, filter.list=filt.list, var.id.col=pat.var.id, filter.sep=";", collapse.col="filter", ignore.cols="variant_qual_id")
	nrow(raw.vars)##286
	nrow(filt.vars)##206
	
}
}

\keyword{methods}
