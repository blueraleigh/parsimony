.sankoff.data = function(data, levs, ambig)
{
    n = length(data)
    if (missing(levs)) {
        if (is.factor(data)) {
            r = nlevels(data)
            g = matrix(Inf, n, r)
            data = match(data, levels(data))
            idx = (seq_len(n)-1L) + (data-1L)*n
            g[idx+1] = 0
        } else {
            stopifnot(is.integer(data))
            stopifnot(all(data > 0L))
            stopifnot(all(tabulate(data, max(data)) > 0L))
            r = length(unique(data))
            g = matrix(Inf, n, r)
            idx = (seq_len(n)-1L) + (data-1L)*n
            g[idx+1] = 0
        }
    } else {
        r = length(levs)
        g = matrix(Inf, n, r)
        f = match(data, levs)
        missing = is.na(f)

        idx = (seq_len(n)[!missing]-1) + (f[!missing]-1L)*n
        g[idx+1] = 0

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
                g[i, p] = 0
            }
        }
    }

    return (g)
}


# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Sankoff's algorithm
mpr.sankoff = function(phy, data, cost, levels, ambig) {
    stopifnot(is.tree(phy))
    stopifnot(!is.null(names(data)))

    data = data[tiplabels(phy)]
    g = .sankoff.data(data, levels, ambig)
    r = ncol(g)
    g = rbind(g, matrix(0, Nnode(phy)-Ntip(phy), r))
    h = matrix(0, Nnode(phy), r)
    f = matrix(0, Nnode(phy), r)

    if (missing(cost))
        cost = matrix(1, r, r)
    diag(cost) = 0
    stopifnot(is.matrix(cost) && is.numeric(cost))
    stopifnot(nrow(cost) == ncol(cost))
    stopifnot(nrow(cost) == r)

    .Call(parsimony_sankoff_mpr_downpass, phy, r, g, h, cost)
    .Call(parsimony_sankoff_mpr_uppass, phy, r, g, h, f, cost)

    mpr_count = .Call(parsimony_sankoff_mpr_count, phy, r, g, h, f, cost)

    list(

        mpr = f

        , mpr_count = mpr_count

        , simulate = function(n) {
            .Call(parsimony_sankoff_mpr_sample, phy, r, g, h, f, cost,
                as.integer(n), TRUE)
        }
    )
}
