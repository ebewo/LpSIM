from __future__ import division
import numpy as nu

def shuffle_dna(Dsize, i):
    if i < (0.3 * Dsize):
        Pa, Pc, Pg, Pt = 40, 50, 60, 100
    elif i < (0.6 * Dsize):
        Pa, Pc, Pg, Pt = 25, 50, 75, 100
    elif i < Dsize:
        Pa, Pc, Pg, Pt = 10, 50, 90, 100

    return Pa, Pc, Pg, Pt

def shuffle_dna2(Dsize, i):
    if i < (0.2 * Dsize):
        Pa, Pc, Pg, Pt = 25, 50, 75, 100
    elif i < (0.5 * Dsize):
        Pa, Pc, Pg, Pt = 40, 50, 60, 100
    elif i < (0.8 * Dsize):
        Pa, Pc, Pg, Pt = 10, 50, 90, 100
    elif i < Dsize:
        Pa, Pc, Pg, Pt = 25, 50, 75, 100

    return Pa, Pc, Pg, Pt

def generate_dna(dna_parameters):
    # generate random DNA string with different probabilities(Pa,Pc,Pg,Pt)
    # for each nucleotide
    dna = []
    Pa = (100 - dna_parameters['GC']) / 2
    Pc = (dna_parameters['GC'] / 2) + Pa
    Pg = (dna_parameters['GC'] / 2) + Pc
    Pt = Pa + Pg

    for i in xrange(dna_parameters['seq_size']):
        if dna_parameters['shuffle_opt'] == 1:
            Pa, Pc, Pg, Pt = shuffle_dna(dna_parameters['seq_size'], i)
        elif dna_parameters['shuffle_opt'] == 2:
            Pa, Pc, Pg, Pt = shuffle_dna2(dna_parameters['seq_size'], i)
        r = nu.random.random()
        # create clumps of CG to create CpG islands, degree of clumpness is
        # based on pclump
        if nu.random.random() <= dna_parameters['pclump'] and len(dna) != 0 and dna[i - 1] == "C":
            dna.append("G")
        elif r <= Pa / 100:
            dna.append("A")
        elif r <= Pc / 100:
            dna.append("C")
        elif r <= Pg / 100:
            dna.append("G")
        elif r <= Pt / 100:
            dna.append("T")

    dna = "".join(dna)
    # save dna to file
    fname = "sequences/simulated_dna.txt"
    file = open(fname, "w")
    file.write(dna)
    file.close()
    return dna
