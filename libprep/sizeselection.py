import numpy as nu
import pandas as pd

def size_selection(frags,std):
    a_frags =[]
    acpt=[]
    post_frg = [x[1][2] for x in frags]
    post_frg = nu.array(post_frg)
    count = len(post_frg)
    rng = nu.ptp(post_frg)
    intervals = nu.round(nu.sqrt(count), decimals=0)
    spacing = rng/(intervals-1)
    bin_list = nu.arange(post_frg.min(),post_frg.max()+spacing,spacing)
    pfs = pd.Series(post_frg).value_counts(bins = bin_list, normalize = True).sort_index().reset_index().values #as_matrix()
    cum = nu.cumsum([x[1] for x in pfs])
    bns = [x[0].left for x in pfs] #[x[0] for x in pfs] <- pandas 0.19 in 0.23 use left to access left interval
    acpt = [x[0] for x in zip(bns,cum) if 0.06 <= x[1]]

    for i in frags:
        if acpt[0] <= i[1][2] <=  acpt[-1]:
            a_frags.append(i)
    return a_frags