#include <R.h>
#include <Rinternals.h>
#include <phy.h>

/*
** Fitch parsimony routines. We represent state sets
** as bit masks using 32 bit signed integers. For example
** the state set {1, 3, 7} is represented as
**
** 00000000000000000000000001000101
**
** corresponding to the integer value 138.
**
** Consequently, these routines can only handle a maximum
** of 31 states.
**
** States are numbered starting from 1 on the R side but
** starting from 0 on the C side.
*/


static int stateset_contains(int state, int stateset) {
    return (stateset & (1<<state)) != 0 ? 1 : 0;
}


static int stateset_add(int state, int stateset) {
    return stateset | 1<<state;
}


static int stateset_remove(int state, int stateset) {
    return stateset & ~(1<<state);
}


static void fitch_downpass(struct phy *phy, int *g, int *pscores)
{
    int parent;
    int lfchild;
    int rtchild;
    int pscore;
    struct phy_node *node;
    struct phy_cursor *cursor;

    memset(pscores, 0, phy_nnode(phy) * sizeof(int));

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        pscore = 0;
        parent = phy_node_index(node);
        lfchild = phy_node_index(phy_node_lfdesc(node));
        rtchild = phy_node_index(phy_node_rtdesc(node));

        g[parent] = g[lfchild] & g[rtchild];

        if (!g[parent]) {
            g[parent] = g[lfchild] | g[rtchild];
            ++pscore;
        }

        pscores[parent] += pscores[lfchild];
        pscores[parent] += pscores[rtchild];
        pscores[parent] += pscore;
    }
}


static void fitch_uppass(struct phy *phy, int *f, int *g)
{
    int focal;
    int parent;
    int lfchild;
    int rtchild;
    struct phy_node *node;
    struct phy_cursor *cursor;

    memcpy(f, g, phy_ntip(phy) * sizeof(int));

    f[phy_node_index(phy_root(phy))] = g[phy_node_index(phy_root(phy))];

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        PREORDER);
    phy_cursor_step(cursor);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        focal = phy_node_index(node);
        parent = phy_node_index(phy_node_anc(node));
        lfchild = phy_node_index(phy_node_lfdesc(node));
        rtchild = phy_node_index(phy_node_rtdesc(node));

        if (f[parent] == (g[focal] & f[parent])) {
            f[focal] = f[parent];
        } else if ((g[lfchild] & g[rtchild]) == 0) {
            f[focal] = g[focal] | f[parent];
        } else {
            f[focal] = g[lfchild] | g[rtchild];
            f[focal] &= f[parent];
            f[focal] |= g[focal];
        }
    }
}


/* Compute MPR state sets */
SEXP parsimony_fitch_mpr(SEXP rtree, SEXP r, SEXP up, SEXP down, SEXP pscore)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    fitch_downpass(phy, ddata, INTEGER(pscore));
    fitch_uppass(phy, udata, ddata);

    SEXP ans = PROTECT(allocMatrix(INTSXP, phy_nnode(phy), *INTEGER(r)));

    for (int i = 0; i < phy_nnode(phy); ++i) {
        for (int j = 0; j < *INTEGER(r); ++j) {
            if (stateset_contains(j, udata[i]))
                INTEGER(ans)[i + j*phy_nnode(phy)] = 1;
            else
                INTEGER(ans)[i + j*phy_nnode(phy)] = 0;
        }
    }

    UNPROTECT(1);
    return ans;
}


/* Uniformly sample a 0-valued index from an r-length
** binary vector (reservoir sampling) */
static int choose_state(int r, int stateset)
{
    int i;
    int j;
    double k = 1;
    for (i = 0; i < r; ++i) {
        if (stateset_contains(i, stateset)) {
            if (unif_rand() < 1/k)
                j = i;
            k += 1;
        }
    }
    return j;
}


/* Sample a MPR history */
static void fitch_history(struct phy *phy, int r, int nsample, int *nodestate,
    int *f, int *g)
{
    int k;
    int nnode = phy_nnode(phy);
    int root = phy_node_index(phy_root(phy));
    int focal;
    int parent;
    int cstate;
    int pstate;

    struct phy_node *node;
    struct phy_cursor *cursor;

    for (k = 0; k < nsample; ++k)
    {
        nodestate[root + k * nnode] = choose_state(r, f[root]) + 1;

        cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, PREORDER);
        phy_cursor_step(cursor);

        while ((node = phy_cursor_step(cursor)) != NULL)
        {
            focal = phy_node_index(node);
            parent = phy_node_index(phy_node_anc(node));
            pstate = nodestate[parent + k * nnode] - 1;

            if (stateset_contains(pstate, g[focal]))
            {
                cstate = pstate;
            }
            else
            {
                if (stateset_contains(pstate, f[focal]))
                {
                    // temporarily modify the downpass state set
                    // to include the parent state, which is a
                    // valid MPR state-to-state transition in this case.
                    g[focal] = stateset_add(pstate, g[focal]);
                    cstate = choose_state(r, g[focal]);
                    g[focal] = stateset_remove(pstate, g[focal]);
                }
                else
                {
                    cstate = choose_state(r, g[focal]);
                }
            }
            nodestate[focal + k*nnode] = cstate + 1;
        }
    }
}


SEXP parsimony_fitch_mpr_sample(SEXP rtree, SEXP r, SEXP nsample,
    SEXP up, SEXP down)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    SEXP nodestate = PROTECT(allocMatrix(
        INTSXP, phy_nnode(phy), *INTEGER(nsample)));

    GetRNGstate();

    fitch_history(phy, *INTEGER(r), *INTEGER(nsample), INTEGER(nodestate),
        INTEGER(up), INTEGER(down));

    PutRNGstate();

    UNPROTECT(1);
    return nodestate;
}


/* Count the number of MPR reconstructions
**
** Algorithm modified from Mesquite project
** https://github.com/MesquiteProject/MesquiteCore/blob/master/Source/mesquite/parsimony/lib/MPRProcessor.java
*/
static double fitch_count(struct phy *phy, int r, int *f, int *g)
{
    int i;
    int j;
    int k;
    int parent;
    int child;
    int nnode = phy_nnode(phy);
    double np;
    double nbr = 0;
    double nh[nnode * r];
    struct phy_node *d;
    struct phy_node *node;
    struct phy_cursor *cursor;

    memset(nh, 0, (nnode * r) * sizeof(double));

    for (i = 0; i < phy_ntip(phy); ++i) {
        for (j = 0; j < r; ++j) {
            if (stateset_contains(j, g[i]))
                nh[i + j * nnode] = 1;
        }
    }

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0) {
        parent = phy_node_index(node);
        for (j = 0; j < r; ++j) {
            if (stateset_contains(j, f[parent]))
                nh[parent + j * nnode] = 1;
        }
        for (d = phy_node_lfdesc(node); d; d = phy_node_next(d)) {
            child = phy_node_index(d);
            for (j = 0; j < r; ++j) {
                if (stateset_contains(j, f[parent])) {
                    np = 0;
                    if (stateset_contains(j, g[child])) {
                        np += nh[child + j * nnode];
                    } else {
                        for (k = 0; k < r; ++k) {
                            if (stateset_contains(k, g[child]) ||
                                stateset_contains(k, f[child]))
                            {
                                np += nh[child + k * nnode];
                            }
                        }
                    }
                    nh[parent + j * nnode] *= np;
                }
            }
        }
    }

    for (j = 0; j < r; ++j)
        nbr += nh[phy_node_index(phy_root(phy)) + j * nnode];

    return nbr;
}


SEXP parsimony_fitch_mpr_count(SEXP rtree, SEXP r, SEXP up, SEXP down)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    return ScalarReal(fitch_count(phy, *INTEGER(r), udata, ddata));
}
