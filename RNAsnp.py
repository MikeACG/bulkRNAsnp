import sys
import ctypes
import numpy as np
import pandas as pd
from ViennaRNA import RNA
#import time
inputPath = sys.argv[1]
outputPath = sys.argv[2]

# # # debug
# inputPath = "~/projects/translateSelection/MC3/producedData/fold/ACC/jobs/h1-1/h2-1/h3-1/batch1/input.tsv"
# outputPath = "~/projects/translateSelection/MC3/producedData/fold/ACC/jobs/h1-1/h2-1/h3-1/batch1/ouput.tsv"

print("INIT", flush = True)

# setup rnasnp with ctypes
clibrary = ctypes.CDLL("./rnasnp.so")
rnasnp = clibrary.compute_bpd
rnasnp.restype = None
rnasnp.argtypes = [
        np.ctypeslib.ndpointer(dtype = np.uintp, ndim = 1, flags = 'C'),
        np.ctypeslib.ndpointer(dtype = np.uintp, ndim = 1, flags = 'C'),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype = np.int32, ndim = 1, flags = 'C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype = np.float64, ndim = 1, flags = 'C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype = np.int32, ndim = 1, flags = 'C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype = np.float64, ndim = 1, flags = 'C_CONTIGUOUS'),]

def getbppm(rna, n, W = 200, L = 150, c = 0.0):

    bppm = np.zeros((n + 1, n + 1), 'float64')

    bppl = RNA.pfl_fold(rna, W, L, c)
    if len(bppl) == 0: return bppm
    
    i = np.array([ep.i for ep in bppl])
    j = np.array([ep.j for ep in bppl])
    p = np.array([ep.p for ep in bppl])

    bppm[i, j] = p
    bppm[j, i] = p

    return bppm

def dirRNAsnp(mutbppm, wtbppm, n, regionX = 20, regionY = 150):

    m = (mutbppm.__array_interface__['data'][0] + np.arange(mutbppm.shape[0]) * mutbppm.strides[0]).astype(np.uintp)
    w = (wtbppm.__array_interface__['data'][0] + np.arange(wtbppm.shape[0]) * wtbppm.strides[0]).astype(np.uintp)
    start = np.zeros(1, 'int32')
    dmax = np.zeros(1, 'float64')
    end = np.zeros(1, 'int32')
    d = np.zeros(1, 'float64')
    rnasnp(w, m, 1, n, regionX, regionY, start, dmax, end, d)

    start = start[0]
    dmax = dmax[0]
    end = end[0]
    d = d[0]

    # matrix bounds are violated when summary stats are not computable (some edge cases like matrices being all zeros or having very few small differences and the seqs being only a few nucleotides long)
    # represent this case with missing values in start/end as -1
    if (start < 1 or start > n) or (end < 1 or end > n) or (start > end):
        # if both matrices are equal we can say the difference between them is 0, else we can't say anything
        if ((mutbppm == wtbppm).all()):
            return -1, -1, 0.0, 0.0, 0.0
        return -1, -1, float('nan'), float('nan'), float('nan')

    dmedp = deltaMedianProb(mutbppm, wtbppm, start, end)

    return start, end, dmax, d, dmedp

def deltaMedianProb(mutbppm, wtbppm, start, end1):

    end0 = end1 + 1 # 0-based end to include end1 index itself

    return np.median(
        mutbppm[start:end0, start:end0].sum(1) \
        - wtbppm[start:end0, start:end0].sum(1))

def compare(mutrna, n, wtbppm):

    mutbppm = getbppm(mutrna, n)
    stats = dirRNAsnp(mutbppm, wtbppm, n)

    return stats

def processRow(wtrna, mutrnaList, n, i, l):

    # verbosity
    if i % 1000 == 0: 
        print(str(i) + "/" + str(l) + "...", flush = True)

    # get wild-type base pair probability matrix
    wtbppm = getbppm(wtrna, n)

    # compare each mutant to the wild type
    r = [compare(m, n, wtbppm) for m in mutrnaList]

    return r

# load the inpput
print("READING INPUT...", flush = True)
indf = pd.read_table(inputPath, header = 0)
indf['n'] = indf['seq'].str.len() # length of each sequence

# Process each row
print("MAIN...", flush = True)
nrows = indf.shape[0]
R = [processRow(w, [m1, m2, m3], n, i, nrows) for w, m1, m2, m3, n, i
     in zip(indf['seq'], indf['mutant1'], indf['mutant2'], indf['mutant3'], indf['n'], range(1, nrows + 1))]

# write result to disk
print("WRITING OUTPUT...", flush = True)
snps = [[m1, m2, m3] for m1, m2, m3 in zip(indf['snp1'], indf['snp2'], indf['snp3'])] # get mutated bases in the right order
outdf = pd.concat(
    [
        pd.DataFrame([tpl for nl in R for tpl in nl], columns = ['start', 'end', 'dmax', 'd', 'dmedp']), # from flattened results list
        pd.Series([snp for nl in snps for snp in nl], name = 'snp'), # flattened snps list
        indf[['transcript_id', 'position', 'wt', 'position.abs', 'chr', 'n']] \
            .loc[indf.index.repeat(3)].reset_index(drop = True)
    ],
    axis = 1
)
outdf.to_csv(outputPath, sep = "\t", index = False, na_rep = 'NA')

print("DONE", flush = True)

#### some usage considerations
# it takes ~ 5.33 hours to run this script for 14,590 records of the input data frame, the median seq length of the records was 401
# the max memory used for this setup was about 453.969M
# these timings are for a computing cluster which unfortunately has a lot of variance between nodes
# certain nodes won't be able to complete 14,000 in 6 hours, probably its safer to settle for 12,000 max
