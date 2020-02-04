
import numpy as nu

def bias_counter(x):
    # check gc and at content of a fragment
    gc = ((x.count("G") + x.count("C")) / len(x))
    at = ((x.count("A") + x.count("T")) / len(x))
    return gc, at

def fragment_dna(m, s, Nfrags, PsplitCG,dna):
    nu.random.seed()
    # fragment DNA string with fragment sizes from a known distribution of mean(m), standard deviation(s) and
    # size(Nfrags). Split DNA prefrentially at CpG sites using a probability PsplitCG

    # calculate parameters mu and std of a lognormal distribution given the mean(m) and standard deviation of
    # the lognormal distribution
    v = s**2
    mu = nu.log(((m**2) / (nu.sqrt(v + m**2))))
    std = nu.sqrt(nu.log((v / (m**2) + 1)))

    # generate fragment sizes from a lognormal distribution of size Nfrags
    frg = nu.random.lognormal(mu, std, Nfrags)
    frg = nu.round(frg)
    frg = frg.tolist()
    flist_size= len(frg)

    # z DNA size
    # CpG counter for splits in CpG sites
    # nCpG counter for splits not in CpG sites
    # frags list of fragments
    # a controls while loop
    z = len(dna)
    CpG = 0
    nCpG = 0
    frags = []
    for h,i in  enumerate(frg): 
        a = 0
        # i fragment sizes
        i = int(i)
        # calculate latest fragment start
        #b = z - i
        # while condition ensures every fragment is used
        while a == 0:
            if h < flist_size/2:
                d = nu.random.random_integers(0, z-i)
                c = i + d
                e = d - 1
                if c < len(dna): #10000 b4
                    c=c
                else:
                    c=len(dna)-1 #9999
            else:
                c = nu.random.random_integers(i, z)
                d = c - i
                e = d - 1
                if d > 0:
                    d=d
                else:
                    d=0
                    e=0
            sz=c-d
            # check nucleotide identity of breakage point and the preceding nucleotide and allow breakage
            # if breakage point is G and preceeding nucleotide is C and probabilty PsplitCG is TRUE
            if dna[e] == "C" and dna[d] == "G" and nu.random.random() <= PsplitCG:
                # select fragment from dna string
                x = dna[d:c]
                # GC/AT counter
                gc, at = bias_counter(x)
                # append fragment and other neccasary information to list of fragments
                frags.append(
                    [x, [d, c, sz, gc, at, "R1", "a", "!ADP", "!PCRD"]])
                # satisfy while condition if DNA is split
                a = 1
                # CpG split counter
                CpG += 1
            # check identity of breakage point and preceding nucleotide and allow breakage if breakage point is not G
            # or preceding breakage point is not C and reduced probability of PsplitCG is TRUE
            elif (dna[e] != "C" or dna[d] != "G") and nu.random.random() >= PsplitCG:
                x = dna[d:c]
                gc, at = bias_counter(x)
                frags.append(
                    [x, [d, c, sz, gc, at, "R1", "b", "!ADP", "!PCRD"]])
                a = 1
                nCpG += 1
    print("Fragmentation complete!")
    return frags, frg
