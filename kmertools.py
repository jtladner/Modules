from collections import defaultdict
import fastatools as ft
import itertools as it
try:
    from Bio.Seq import Seq
except:
    print("You must have Biopython installed to use the Revcomp functions")

# Pattern should be a string of 0s and 1s with 0s indicating positions to discard and 1s indicating positions to include
def patternedKmers(seq, pattern, outType="set",filter=[]):
    out=[]
    k=len(pattern)
    for i in range(len(seq)-k+1):
        this = seq[i:i+k]
        flag = sum([1 for p in this if p in filter])
        if not flag:
            this = "".join([a for i,a in enumerate(this) if int(pattern[i])])
            out.append(this)
    
    if outType == "set":
        return set(out)
    elif outType == "dict":
        outD = defaultdict(int)
        for k in out:
            outD[k]+=1
        return outD
    else:
        print(f"{outType} is not a supported output type. Options are 'set' and 'dict'.")

# Function to generate relevant patterns for use with patternedKmers function
#def generatePatterns(k, span):
    

# Returns a dictionary with one key for each seq in a fasta (sequence name)
# Values will be sets of all kmers contained
def kmerDictSetFasta(fasta,k,filter=[]):
    fD = ft.read_fasta_dict_upper(fasta)
    outD = kmerDictSet(fD,k,filter)
    return outD

# Returns a dictionary with one key for each seq in a fasta dict (sequence name)
# Values will be sets of all kmers contained
# The same as kmerDictSetFasta, but starting with a fasta file that has already been loaded into a dictionary
def kmerDictSet(fD,k,filter=[]):
    outD = {}
    for n,s in fD.items():
        outD[n] = kmerSet(s,k, filter)
    return outD


#Returns set containing all unique kmers.
def kmerSetFasta(fasta,k,filter=[]):
    names, seqs = ft.read_fasta_lists(fasta)
    total = set()
    for s in seqs:
        total.update(kmerSet(s,k, filter))
    return total

#Returns set containing all unique kmers, including those for the reverse complement
def kmerSetFastaRevcomp(fasta,k,filter=[], upper=True):
    names, seqs = ft.read_fasta_lists(fasta)
    if upper:
        seqs = [s.upper() for s in seqs]
    total = set()
    for s in seqs:
        total.update(kmerSet(s,k, filter))
        dna = Seq(s)
        rcs = str(dna.reverse_complement())
        total.update(kmerSet(rcs,k, filter))

    return total



#Returns dictionary containing counts for unique kmers.
def kmerDictCountFasta(fasta,k,filter=[]):
    names, seqs = ft.read_fasta_lists(fasta)
    cD = defaultdict(int)
    for s in seqs:
        ks = kmerSet(s,k, filter)
        for kmer in ks:
            cD[kmer]+=1
    return cD

#Returns dictionary containing counts for unique kmers.
def kmerDictCount(seqs,k,filter=[]):
    cD = defaultdict(int)
    for s in seqs:
        ks = kmerSet(s,k, filter)
        for kmer in ks:
            cD[kmer]+=1
    return cD


#Returns set containing all unique kmers.
def kmerSet(seq,k, filter=[]):
    out=[]
    for i in range(len(seq)-k+1):
        this = seq[i:i+k]
        flag = sum([1 for p in this if p in filter])
        if not flag:
            out.append(this)
    return set(out)

#Returns list containing all unique kmers.
def kmerList(seq,k):
    out=[]
    for i in range(len(seq)-k+1):
        out.append(seq[i:i+k])
    return list(set(out))

#Returns dict with keys for each unique kmer. All values will be empty strings.
def kmerEmptyDict(seq,k, circular=False):
    out=[]
    for i in range(len(seq)-k+1):
        out.append(seq[i:i+k])
    
    if circular:
        bridge = seq[-k:] + seq[:k]
        for i in range(len(bridge)-k+1):
            out.append(bridge[i:i+k])
    
    return {x:"" for x in set(out)}

#Returns proportion of identical kmers
def compSeqs(s1, s2, k, filter=[]):
    s1k = kmerSet(s1,k, filter)
    s2k = kmerSet(s2,k, filter)
    if len(s1k) == 0 or len(s2k)==0:
        print("At least one of the seqs have 0 %dmers. Returning 0 overlap." % (k))
        return 0
    elif len(s1k)<=len(s2k):
        return(len(s1k.intersection(s2k))/len(s1k))
    else:
        return(len(s1k.intersection(s2k))/len(s2k))

#####----Below here, haven't checked over, just copied from another script

    
def pairwiseComps(seqs, k):
    comps=[]
    for a,b in it.combinations(seqs, 2):
        comps.append(compSeqs(a,b,k))
    return comps

def btwnComps(g1, g2, k):
    comps=[]
    for a in g1:
        for b in g2:
            comps.append(compSeqs(a,b,k))
    return comps

def cluster(seqs, k, thresh, names):
    longest=[]
    clusts=[]
    clustNames=[]
    for i, each in enumerate(seqs):
        match=0
        if i == 0:
            longest=[each]
            clusts=[[each]]
            clustNames.append([names[i]])
        else:
            for j, rep in enumerate(longest):
                if compSeqs(each, rep,k)>=thresh:
                    clusts[j].append(each)
                    clustNames[j].append(names[i])
                    match=1
                    break
            if not match:
                longest.append(each)
                clusts.append([each])
                clustNames.append([names[i]])
    return clusts, clustNames
    
#Look at epitope merging script for algorithm to combine these clusters
def combPairs(pairs):
    d={}
    groups = pairs[::]
    for each in pairs:
        for x in each:
            d[x] = d.get(x, 0) + 1
    for each, count in d.items():
        if count>1:
            groups = combine(each, groups)
    return groups

def combine(focal, gList):
    newL = []
    combo = []
    for each in gList:
        if focal in each:
            combo+=each
        else:
            newL.append(each)
    newL.append(list(set(combo)))
    return newL


def merge(clusts, names, k, maxThresh, meanThresh):
    pairs2merge = []
    for a,b in it.combinations(range(len(clusts)), 2):
        cs = btwnComps(clusts[a],clusts[b],k)
        if max(cs)>=maxThresh or max(np.mean(cs),np.median(cs))>=meanThresh:
            pairs2merge.append([a,b])
        
    groups2merge = combPairs(pairs2merge)
    
    newClusts = []
    newNames = []
    merged = []
    for each in groups2merge:
        newClusts.append([])
        newNames.append([])
        for a in each:
            merged.append(a)
            newClusts[-1] += clusts[a]
            newNames[-1] += names[a]
    for i in range(len(clusts)):
        if i not in merged:
            newClusts.append(clusts[i])
            newNames.append(names[i])
            
    return newClusts, newNames

