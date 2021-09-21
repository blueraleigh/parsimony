.fitch.data = function(data, levs, ambig)
{
    if (missing(levs)) {
        if (is.factor(data)) {
            stopifnot(nlevels(data) < 32)
            f = bitwShiftL(1L, match(data, levels(data))-1L)
        } else {
            stopifnot(is.integer(data))
            stopifnot(all(data > 0L))
            stopifnot(max(data) < 32)

            f = bitwShiftL(1L, data-1L)
        }
    } else {
        f = match(data, levs)
        missing = is.na(f)
        has_missing = any(missing)
        if (has_missing) {
            stopifnot(!missing(ambig))
            stopifnot(is.list(ambig))
            stopifnot(!is.null(names(ambig)))
            stopifnot(length(intersect(levs, names(ambig))) == 0L)
            for (i in which(missing)) {
                p = match(ambig[[ as.character(data[i]) ]], levs)
                if (anyNA(p))
                    stop("invalid ambiguous character state encoding")
                f[i] = bitwShiftL(1L, p[1L]-1L)
                for (j in p[-1L])
                    f[i] = bitwOr(f[i], bitwShiftL(1L, j-1L))
            }
        }

        f[!missing] = bitwShiftL(1L, f-1L)

    }

    return (f)
}


#' Ancestral character state reconstruction using maximum parsimony
#'
#' Perform a maximum parsimony reconstruction of ancestral character
#' states using Fitch's algorithm
#'
#' @param phy An object of class \code{tree}.
#' @param data A named vector of character state data. If \code{levels} is
#' this should be a \code{factor} or \code{integer} vector.
#' @param levels A vector containing the different states the character can
#' take. Optional.
#' @param ambig A named list containing ambiguous character state encodings.
#' The names of the list correspond to ambiguous character states. The values
#' for each list item refer to states in \code{levels}. Thus, if one of the
#' names in \code{ambig} appears in \code{data}, that data entry is assumed to
#' potentially map to any of the states in \code{levels} that appear in the
#' corresponding \code{ambig} entry. Optional.
#' @note A maximum of 31 character states is allowed.
#' @return A list with five components:
#' \describe{
#' \item{score}{The minimum number of character state changes needed to explain
#' \code{data}.}
#' \item{scores}{The minimum number of character state changes needed to explain
#' the data in each clade of \code{phy}. The i-th entry in this array corresponds
#' to the subtree rooted at the node with index i.}
#' \item{mpr}{A matrix containing optimal state sets for each internal node.
#' Each column represents a state and is either 0 (the optimal state set does
#' not include the state) ) or 1 (the optimal state set includes the state).}
#' \item{mpr_count}{The number of maximum parsimony reconstructions.}
#' \item{simulate}{A function to randomly sample maximum parsimony 
#' reconstructions.}
#' }
mpr.fitch = function(phy, data, levels, ambig) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(!is.null(names(data)))

    data = data[tiplabels(phy)]
    downpass = integer(Nnode(phy))
    uppass = integer(Nnode(phy))
    pscore = integer(Nnode(phy))

    downpass[1:Ntip(phy)] = if (missing(levels))
         .fitch.data(data)
    else
        .fitch.data(data, levels, ambig)

    r = if (missing(levels)) length(unique(data)) else length(levels)

    mpr = .Call(parsimony_fitch_mpr, phy, r, uppass, downpass, pscore)
    mpr_count = .Call(parsimony_fitch_mpr_count, phy, r, uppass, downpass)

    list(
        score = pscore[root(phy)]

        , scores = pscore

        , mpr = mpr

        , mpr_count = mpr_count

        , simulate = function(n) {
            .Call(parsimony_fitch_mpr_sample, phy, r, as.integer(n),
                uppass, downpass)
        }
    )
}
