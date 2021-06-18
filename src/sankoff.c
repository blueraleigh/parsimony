#include <R.h>
#include <Rinternals.h>
#include <phy.h>

/* Calculate node and stem state costs. */
static void sankoff_downpass(struct phy *phy, int r, double *g, double *h,
    double *cost)
{
    int i;
    int j;
    int k;
    int u;
    int v;
    int ntip = phy_ntip(phy);
    int nnode = phy_nnode(phy);
    double t;
    double min_t;
    struct phy_node *node;
    struct phy_node *desc;
    struct phy_cursor *cursor;

#define C(i, j) cost[(i) + r * (j)]
#define G(i, j) g[(i) + nnode * (j)]
#define H(i, j) h[(i) + nnode * (j)]

    for (u = 0; u < ntip; ++u)
    {
        for (i = 0; i < r; ++i)
        {
            min_t = R_PosInf;
            for (j = 0; j < r; ++j)
            {
                t = C(i, j) + G(u, j);
                if (t < min_t)
                    min_t = t;
            }
            H(u, i) = min_t;
        }
    }

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        u = phy_node_index(node);
        for (j = 0; j < r; ++j)
        {
            desc = phy_node_lfdesc(node);
            v = phy_node_index(desc);
            G(u, j) = H(v, j);
            for (desc = phy_node_next(desc); desc; desc = phy_node_next(desc))
            {
                v = phy_node_index(desc);
                G(u, j) += H(v, j);
            }
        }

        for (i = 0; i < r; ++i)
        {
            min_t = R_PosInf;
            for (j = 0; j < r; ++j)
            {
                t = C(i, j) + G(u, j);
                if (t < min_t)
                    min_t = t;
            }
            H(u, i) = min_t;
        }
    }
}


/* Calculate final state costs */
static void sankoff_uppass(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost)
{
    int i;
    int j;
    int focal;
    int parent = phy_node_index(phy_root(phy));
    int nnode = phy_nnode(phy);
    double t;
    double min_t;
    struct phy_node *node;
    struct phy_cursor *cursor;

#define F(i,j) f[(i) + nnode * (j)]

    for (j = 0; j < r; ++j)
        F(parent, j) = G(parent, j);

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, PREORDER);
    phy_cursor_step(cursor);
    while ((node = phy_cursor_step(cursor)) != 0)
    {
        focal = phy_node_index(node);
        parent = phy_node_index(phy_node_anc(node));

        for (i = 0; i < r; ++i)
        {
            min_t = R_PosInf;
            for (j = 0; j < r; ++j)
            {
                t = (F(parent, j) - H(focal, j)) + C(j, i) + G(focal, i);
                if (t < min_t)
                    min_t = t;
            }
            F(focal, i) = min_t;
        }
    }
}


/* Count the number of MPR histories */
static double sankoff_count(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost)
{
    int j;
    int k;
    int u;
    int v;
    int ntip = phy_ntip(phy);
    int nnode = phy_nnode(phy);
    double np;
    double L;
    double nbr = 0;
    double nh[nnode * r];
    struct phy_node *node;
    struct phy_node *child;
    struct phy_cursor *cursor;

    memset(nh, 0, nnode * r * sizeof(double));

#define NH(i, j) nh[(i) + (j) * nnode]

    for (u = 0; u < ntip; ++u)
    {
        for (j = 0; j < r; ++j)
        {
            if (G(u, j) < R_PosInf)
                NH(u, j) = 1;
            else
                NH(u, j) = 0;
        }
    }

    u = phy_node_index(phy_root(phy));
    L = R_PosInf;
    for (j = 0; j < r; ++j)
    {
        if (F(u, j) < L)
            L = F(u, j);
    }

    cursor = phy_cursor_prepare(phy, phy_root(phy), INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        u = phy_node_index(node);
        for (j = 0; j < r; ++j)
        {
            if (F(u, j) > L)
                continue;
            NH(u, j) = 1;
        }
        for (child = phy_node_lfdesc(node); child; child = phy_node_next(child))
        {
            v = phy_node_index(child);
            for (j = 0; j < r; ++j)
            {
                if (F(u, j) > L)
                    continue;
                np = 0;
                for (k = 0; k < r; ++k)
                {
                    if ((F(u, j) - H(v, j) + C(j, k) + G(v, k)) > L)
                        continue;
                    np += NH(v, k);
                }
                NH(u, j) *= np;
            }
        }
    }

    u = phy_node_index(phy_root(phy));
    for (j = 0; j < r; ++j)
        nbr += NH(u, j);

    return nbr;
}


static int sankoff_choose(int r, int stride, double L, double *f, int only_mpr)
{
    int i;
    int j = 0;
    double g;
    double w;
    double max_w = R_NegInf;
    for (i = 0; i < r; ++i, f += stride)
    {
        if (*f > L && only_mpr)
            continue;
        w = L - *f;
        g = -log(-log(unif_rand()));
        if ((g + w) > max_w)
        {
            max_w = g + w;
            j = i;
        }
    }
    return j;
}


/* Sample histories of character evolution */
static void sankoff_sample(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost, int nsample, int *nodestate, int only_mpr)
{
    int i;
    int j;
    int k;
    int u;
    int v;
    int nnode = phy_nnode(phy);
    double L;
    double w[r];
    struct phy_node *node;
    struct phy_cursor *cursor;

    u = phy_node_index(phy_root(phy));
    L = R_PosInf;
    for (j = 0; j < r; ++j)
    {
        if (F(u, j) < L)
            L = F(u, j);
    }

#define NS(i,j) nodestate[(i)+(j)*nnode]

    for (k = 0; k < nsample; ++k)
        NS(u, k) = sankoff_choose(r, nnode, L, f+u, only_mpr) + 1;

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, PREORDER);
    phy_cursor_step(cursor);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        u = phy_node_index(phy_node_anc(node));
        v = phy_node_index(node);
        for (k = 0; k < nsample; ++k)
        {
            i = NS(u, k) - 1;
            for (j = 0; j < r; ++j)
                w[j] = F(u, i) - H(v, i) + C(i, j) + G(v, j);
            NS(v, k) = sankoff_choose(r, 1, F(u, i), w, only_mpr) + 1;
        }
    }
}


SEXP parsimony_sankoff_mpr_downpass(SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP cost)
{
    sankoff_downpass((struct phy *)R_ExternalPtrAddr(rtree), INTEGER(r)[0],
        REAL(g), REAL(h), REAL(cost));
    return R_NilValue;
}


SEXP parsimony_sankoff_mpr_uppass(
    SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP f, SEXP cost)
{
    sankoff_uppass((struct phy *)R_ExternalPtrAddr(rtree), INTEGER(r)[0],
        REAL(g), REAL(h), REAL(f), REAL(cost));
    return R_NilValue;
}


SEXP parsimony_sankoff_mpr_count(
    SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP f, SEXP cost)
{
    return ScalarReal(sankoff_count((struct phy *)R_ExternalPtrAddr(rtree),
        INTEGER(r)[0], REAL(g), REAL(h), REAL(f), REAL(cost)));
}


SEXP parsimony_sankoff_mpr_sample(
    SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP f, SEXP cost, SEXP nsample,
    SEXP only_mpr)
{
    int n = INTEGER(nsample)[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    SEXP nodestate = PROTECT(allocMatrix(INTSXP, phy_nnode(phy), n));
    GetRNGstate();
    sankoff_sample(phy, INTEGER(r)[0], REAL(g), REAL(h), REAL(f), REAL(cost),
        n, INTEGER(nodestate), INTEGER(only_mpr)[0]);
    PutRNGstate();
    UNPROTECT(1);
    return nodestate;
}
