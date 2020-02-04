import numpy as nu
from libprep.fragmentation import fragment_dna
from libprep.pcr import PCR
from libprep.ligation import ligate_dna
from libprep.sizeselection import size_selection
import multiprocessing as mp
import psutil


def mp_count(pfrags, i):
    pfrags = [[x[2][0], x[2][1]] for x in pfrags]
    cov_count = sum(list(
        1 if j[0] <= i <= j[0]+80 else 1 if j[1]-80 <= i <= j[1] else 0 for j in pfrags))

    return cov_count


def mp_count_wrapper(args):

    return mp_count(*args)


def raw_cov(dna, pfrags):

    # single process coverage counter
    # raw_cv = [mp_count(pfrags, i) for i in xrange(len(dna))]

    # multi process coverage counter
    pool = mp.Pool(processes=psutil.cpu_count(logical=False))
    args_gen = [[pfrags, x] for x in range(len(dna))]
    raw_cv = pool.map(mp_count_wrapper, args_gen)
    pool.close()

    print("Coverage analysyis 1 complete!")

    return raw_cv


def mv_avg_cov(cov, n):

    mov_cv = []
    for i in range(len(cov)-n):
        start1 = min(i, len(cov)-n)
        stop2 = min(len(cov), i + n)
        start2 = int(stop2-40)
        stop1 = int(start1+40)
        window1 = nu.mean(cov[start1:stop1])
        window2 = nu.mean(cov[start2:stop2])
        window = [window1, window2]
        avg = nu.mean(window)
        mov_cv.append(avg)
    print("Coverage analysis 2 complete!")

    return mov_cv


def gc_mv_avg(var_dict):
    result = []
    for i in range(len(var_dict['dna'])-var_dict['window']):
        start1 = min(i, len(var_dict['dna'])-var_dict['window'])
        stop2 = min(len(var_dict['dna']), i + var_dict['window'])
        start2 = int(stop2-40)
        stop1 = int(start1+40)
        window1 = var_dict['dna'][start1:stop1]
        GC1 = (window1.count('G') + window1.count('C')) / \
            float(len(window1))*100
        window2 = var_dict['dna'][start2:stop2]
        GC2 = (window2.count('G') + window2.count('C')) / \
            float(len(window2))*100
        window = [GC1, GC2]
        GC = nu.mean(window)
        result.append(GC)

    print("GC analysis complete!")

    return result


def seq(var_dict):
    nu.random.seed()
    frags, frg = fragment_dna(var_dict['mu_frags'], var_dict['sd_frags'], var_dict[
        'no_frags'], var_dict['psplit'], var_dict['dna'])
    frags = size_selection(frags, var_dict['sd_frags'])
    lfrags = ligate_dna(var_dict['pligate'], frags)
    pfrags, pfrags_no_dup = PCR(var_dict['d_temp'], var_dict[
        'el_temp'], var_dict['cycles'], lfrags, var_dict['sd_pcr'])
    raw = raw_cov(var_dict['dna'], pfrags_no_dup)

    return raw


def mv_coverage(var_dict, raw):
    mov = mv_avg_cov(raw, var_dict['window'])

    return mov


def evenness(D):
    C = nu.around(nu.mean(D))
    D2 = [x for x in D if x <= C]
    E = 1-(len(D2)-sum(D2)/C)/len(D)

    return E
