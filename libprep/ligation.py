from __future__ import division
import numpy as nu
import copy

def ligate_dna(p,frags):
    frags=copy.deepcopy(frags)
    ligated_frags = []
    adpC = 0
    bdpC = 0
    cdpC1 = 0
    cdpC2 = 0

    adp = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
    bdp = "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"

    for i in frags:
        i = i[:]
        i[0] = i[0][:]
        i[1] = i[1][:]
        i[0] = "".join([i[0], "A"])

        if i[0][0] == "T" and nu.random.random() >= p:

            if i[0][len(i[0]) - 2] == "A":
                i[1][7] = "TAD5"
            elif i[0][len(i[0]) - 2] != "A":
                i[1][7] = "TND5"
            i[0] = "".join([adp, i[0]])
            adpC += 1

        elif i[0][0] != "T":

            if i[0][len(i[0]) - 2] == "A":
                i[1][7] = "NAD5"
            elif i[0][len(i[0]) - 2] != "A":
                i[1][7] = "NND5"
            i[0] = "".join([adp, i[0]])
            cdpC1 += 1

        if i[0][len(i[0]) - 2] == "A" and nu.random.random() >= p:

            if i[1][7] == "!ADP":
                if i[0][0] == "T":
                    i[1][7] = "TAD3"
                elif i[0][0] != "T":
                    i[1][7] = "NAD3"

            elif i[1][7] == "TAD5" or i[1][7] == "TND5":
                i[1][7] = "TAD"
            elif i[1][7] == "NAD5" or i[1][7] == "NND5":
                i[1][7] = "NAD"
            i[0] = "".join([i[0], bdp])
            bdpC += 1

        elif i[0][len(i[0]) - 2] != "A":

            if i[1][7] == "!ADP":
                if i[0][0] == "T":
                    i[1][7] = "TND3"
                elif i[0][0] != "T":
                    i[1][7] = "NND3"

            elif i[1][7] == "TAD5" or i[1][7] == "TND5":
                i[1][7] = "TND"
            elif i[1][7] == "NAD5" or i[1][7] == "NND5":
                i[1][7] = "NND"
            i[0] = "".join([i[0], bdp])
            cdpC2 += 1

        if i[1][7] == "!ADP":
            if i[0][0] == "T" and i[0][len(i[0]) - 2] == "A":
                i[1][7] = "TA"
            elif i[0][0] == "T" and i[0][len(i[0]) - 2] != "A":
                i[1][7] = "TN"
            elif i[0][0] != "T" and i[0][len(i[0]) - 2] == "A":
                i[1][7] = "NA"
            elif i[0][0] != "T" and i[0][len(i[0]) - 2] != "A":
                i[1][7] = "NN"

        ligated_frags.append(i)

    f_lig = ["TAD", "TND", "NAD", "NND"]
    ligated_frags = [x for x in ligated_frags if x[1][7] in f_lig]
    print "Ligation complete!"
    return ligated_frags
