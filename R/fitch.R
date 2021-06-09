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


# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Fitch's algorithm
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
