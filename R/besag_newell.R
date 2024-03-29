#' Besag-Newell Cluster Detection Method
#' 
#' @description Besag-Newell cluster detection method.  There are differences with the original paper and our implementation:
#' \itemize{
#' \item{we base our analysis on \eqn{k} cases, rather than \eqn{k} *other* cases as prescribed in the paper.} 
#' \item we do not subtract 1 from the *accumulated numbers of other cases* and *accumulated numbers of others at risk*, as was prescribed in the paper to discount selection bias
#' \item M is the total number of areas included, not the number of additional areas included.  i.e. \eqn{M} starts at 1, not 0.
#' \item p-values are not based on the original value of \eqn{k}, rather the actual number of cases observed until we view \eqn{k} or more cases.  Ex:  if \eqn{k = 10}, but as we consider neighbors we encounter 1, 2, 9 then 12 cases, we base our \eqn{p}-values on \eqn{k=12}
#' \item we do not provide a Monte-Carlo simulated \eqn{R}:  the number of tests that attain significance at a fixed level \eqn{\alpha}
#' }
#' The first two and last differences are because we view the testing on an area-by-area level, rather than a case-by-case level.
#' 
#' @param geo an `n x 2` table of the (x,y)-coordinates of the area centroids
#' @param population aggregated population counts for all `n` areas
#' @param cases aggregated case counts for all `n` areas
#' @param expected.cases expected numbers of disease for all `n` areas
#' @param k number of cases to consider
#' @param alpha.level alpha-level threshold used to declare significance
#'
#'
#'
#'
#'
#'
#' @details 
#' For the `population` and `cases` tables, the rows are bunched by areas first, and then for each area, the counts for each strata are listed.  It is important that the tables are balanced:  the strata information are in the same order for each area, and counts for each area/strata combination appear exactly once (even if zero). 
#'
#' @references Besag J. and Newell J. (1991) The Detection of Clusters in Rare Diseases *Journal of the Royal Statistical Society. Series A (Statistics in Society)*, **154**, 143--155
#' @author Albert Y. Kim
#'
#' @return
#' List containing
#' \item{clusters}{information on all clusters that are \eqn{\alpha}-level significant, in decreasing order of the \eqn{p}-value}
#' \item{p.values}{for each of the \eqn{n} areas, \eqn{p}-values of each cluster of size at least \eqn{k}}
#' \item{m.values}{for each of the \eqn{n} areas, the number of areas need to observe at least \eqn{k} cases}
#' \item{observed.k.values}{based on `m.values`, the actual number of cases used to compute the \eqn{p}-values}
#' 
#' @note The `clusters` list elements are themselves lists reporting:\cr\cr
#'\tabular{ll}{
#'  `location.IDs.included` \tab ID's of areas in cluster, in order of distance\cr
#'  `population` \tab population of cluster\cr
#'  `number.of.cases` \tab number of cases in cluster\cr
#'  `expected.cases` \tab expected number of cases in cluster\cr
#'  `SMR` \tab estimated SMR of cluster\cr
#'  `p.value` \tab \eqn{p}-value\cr
#' }
#' 
#' 
#' 
#' @export
#'
#' @examples
#' ## Load Pennsylvania Lung Cancer Data
#'data(pennLC)
#'data <- pennLC$data
#'
#' ## Process geographical information and convert to grid
#'geo <- pennLC$geo[,2:3]
#'geo <- latlong2grid(geo)
#'
#' ## Get aggregated counts of population and cases for each county
#' population <- tapply(data$population,data$county,sum)
#' cases <- tapply(data$cases,data$county,sum)
#'
#' ## Based on the 16 strata levels, computed expected numbers of disease
#' n.strata <- 16
#' expected.cases <- expected(data$population, data$cases, n.strata)
#'
#' ## Set Parameters
#' k <- 1250
#' alpha.level <- 0.05
#'
#' # not controlling for stratas
#' results <- besag_newell(geo, population, cases, expected.cases=NULL, k, 
#'                        alpha.level)
#'
#' # controlling for stratas
#' results <- besag_newell(geo, population, cases, expected.cases, k, alpha.level)
#' 
besag_newell <-
function(geo, population, cases, expected.cases=NULL, k, alpha.level){

#-------------------------------------------------------------------------------
# Initialization 
#-------------------------------------------------------------------------------
# If no expected.cases provided, set them if there are no expected counts
if(is.null(expected.cases)){
	p <- sum(cases)/sum(population)
	expected.cases <- population*p
}

# geographical information computation
geo.results <- zones(geo, population, 1)
nearest.neighbors <- geo.results$nearest.neighbors
distance <- geo.results$dist
n.zones <- length(unlist(nearest.neighbors))


#-------------------------------------------------------------------------------
# Observed statistic computation
#-------------------------------------------------------------------------------
results <- besag_newell_internal(cases, expected.cases, nearest.neighbors, 
                                 n.zones, k)

# observed p.values for each areas
p.values <- results$observed.p.values	
# observed number of neighbors needed to observe k cases
m.values <- results$observed.m.values	
# actual observed number of cases
k.values <- results$observed.k.values	

# pick out areas that were significant and order them by p-value
signif.indices <- order(p.values)[1:sum(p.values <= alpha.level)]

# order remaining values
signif.p.values <- p.values[signif.indices]
signif.m.values <- m.values[signif.indices]
signif.k.values <- k.values[signif.indices]


# Create object to output
# If none are significant, return NULL
if(length(signif.indices) == 0){
	clusters <- NULL
} else {
	clusters <- vector("list", length=length(signif.indices))

	for( i in 1:length(clusters) ){	
		# find areas included in cluster
		cluster <- order(distance[signif.indices[i],])[1:signif.m.values[i]]
		
		new.cluster <- list(
			location.IDs.included = cluster,
			population = sum(population[cluster]),
			number.of.cases = sum(cases[cluster]),
			expected.cases = sum(expected.cases[cluster]),
			SMR = sum(cases[cluster])/sum(expected.cases[cluster]),
			p.value = signif.p.values[i]	
		)
		clusters[[i]] <- new.cluster
	}
}


#-------------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------------
results <- list(
	clusters=clusters,
	p.values=p.values,
	m.values=m.values,
	observed.k.values=k.values
)	
return(results)
}
