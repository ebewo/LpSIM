from __future__ import division
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import copy
from scipy.stats import norm
import numpy as nu

def PCR(tmp, el_tmp, cycles, ligated_frags, sd_pcr):
    n_pcrd = 0
    n_elnd = 0
    PCR.pcrd_fragments_cov = []
    PCR.pcrd_frags = []
    ligated_frags = copy.deepcopy(ligated_frags)
    PCR.dd = []
    PCR.pcrd = 0
    PCR.elnd = 0
    PCR.n_adp = 0
    PCR.h_adp = 0
    PCR.f_adp = 0
    PCR.sd = sd_pcr

    def PCR_process(i):
        f_lig = ["TAD", "TND", "NAD", "NND"]
        h_lig = ["TAD5", 'TND5', "NAD5", "NND5",
                 "TAD3", "NAD3", "TND3", "NND3"]
        n_lig = ["TA", "TN", "NA", "NN"]
        
        i.insert(1, i[0].complement())
        bb = [copy.deepcopy(i)]
        c_count = 0
        d_count = 0

        for j in xrange(cycles):
            cc = []
            c_count += 1

            for k in bb:
                a = k[0].complement()
                b = k[1].complement()
                d_count += 1

                if k[2][7] in f_lig:
                    k[2][8] = "PCRD"
                    cc.extend([[k[0], a, k[2]], [b, k[1], k[2]]])
                    PCR.f_adp += 1
                    if j == 0:
                        PCR.pcrd_fragments_cov.append(k)
                elif k[2][7] in h_lig:
                    k[2][8] = "PCR_FAIL_2"
                    cc.extend([[k[0], a, k[2]]])
                    PCR.h_adp += 1
                elif k[2][7] in n_lig:
                    k[2][8] = "PCR_FAIL_3"
                    PCR.n_adp += 1

            bb = []
            bb.extend(cc)
        PCR.dd.extend(bb)
        PCR.pcrd += 1
        PCR.pcrd_frags = PCR.dd
        return PCR.pcrd_frags, PCR.pcrd_fragments_cov

    for i in ligated_frags:
        i[0] = Seq(i[0])
        dt_diff = mt.Tm_NN(i[0]) - tmp
        el_diff = el_tmp - mt.Tm_NN(i[0])

        if dt_diff < 0:
            if el_diff < 0:
                PCR.pcrd_frags = PCR_process(i)
            else:
                pd = norm.pdf(el_diff, 0, PCR.sd)/norm.pdf(0, 0, PCR.sd)
                if nu.random.uniform(0,1) <= pd:
                    PCR.pcrd_frags = PCR_process(i)
                else:
                    n_elnd += 1
        else: 
            pd = norm.pdf(dt_diff, 0, PCR.sd)/norm.pdf(0, 0, PCR.sd)
            if nu.random.uniform(0,1) <= pd:
                if el_diff < 0:
                    PCR.pcrd_frags = PCR_process(i)    
                else:
                    pd = norm.pdf(el_diff, 0, PCR.sd)/norm.pdf(0, 0, PCR.sd)
                    if nu.random.uniform(0,1) <= pd:
                        PCR.pcrd_frags = PCR_process(i)
                    else:
                        n_elnd += 1
            else:
                n_pcrd += 1

    return PCR.pcrd_frags, PCR.pcrd_fragments_cov