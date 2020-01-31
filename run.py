import numpy as nu
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from libprep.coverage import seq, mv_coverage, gc_mv_avg, evenness
from Bio import SeqIO
import time

# parameters and files
var_dict = {'mu_frags': 300, 'sd_frags': 30, 'no_frags': 5000, 'psplit': 0, 'pligate': 0,
                'd_temp': 98, 'el_temp': 50, 'cycles': 1, 'window': 80, 'fasta': "sequences/phix.fasta", 'sd_pcr': 1}

# utilise fasta file
# record = SeqIO.read(var_dict['fasta'], "fasta")
# dna = str(record.seq)

# utilise simulated DNA
dna_filename = 'sequences/simulated_dna.txt'
file = open(dna_filename, "r")
dna = file.read()


# Add DNA to parameter dictionary
var_dict['dna'] = dna

# run simulator
coverage = seq(var_dict)
moving_coverage = mv_coverage(var_dict, coverage)
GC = gc_mv_avg(var_dict)
E = evenness(coverage)
df_coverage = pd.DataFrame(nu.array([list(dna),coverage]).transpose(), columns = ["base position", "coverage"])
df_moving_coverage = pd.DataFrame(nu.array([list(dna[:-80]),GC,moving_coverage]).transpose(), columns = ["base position","GC","moving coverage"])

# plot coverage
gs = gridspec.GridSpec(2, 2, height_ratios=[6, 1])
fig1 = []
fig1 = plt.figure(1)
ax = fig1.add_subplot(gs[0,:])
lnse1 = ax.plot(moving_coverage, 'y', label = '100*94')
ax.set_ylabel("Coverage")
ax.set_xlabel("GC level")
ax.xaxis.labelpad =10
plt.title("Evenness of coverage = "+str(E))
ax5 = fig1.add_subplot(gs[1,:])
ax.tick_params(labelbottom=False)  
x = range(len(moving_coverage))
extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
ax5.imshow(nu.asarray(GC)[nu.newaxis,:], cmap="jet", aspect="auto", extent=extent)
ax5.set_yticks([])
ax.set_xlim(extent[0], extent[1])
ax5.set_xlabel("nucleotide position")
fig1.set_size_inches(10, 5)

# save files
timestr = time.strftime("%Y%m%d-%H%M%S")
df_coverage.to_csv('results/coverage_'+timestr+'.csv', index = False)
df_moving_coverage.to_csv('results/moving_coverage_'+timestr+'.csv', index = False)
fig1.savefig('results/coverage_'+timestr+'.png')
print "Coverage files generated"

# print evenness of coverage
print "Evenness of coverage = "+str(E)
