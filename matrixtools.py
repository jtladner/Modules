import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt

def writeCounts(cd, outname):
    probeNames = sorted(cd[list(cd.keys())[0]].keys())
    sampNames =  sorted(list(cd.keys()))
    with open(outname, "w") as fout:
        fout.write("Probe\t%s\n" % ("\t".join(sampNames)))
        for each in probeNames:
            fout.write("%s\t%s\n" % (each, "\t".join([str(cd[x][each]) for x in sampNames])))

def parseCounts(countFile, delim="\t", valType="float"):
    counts={}
    with open(countFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols=line.rstrip("\n").split(delim)
            if lc == 1:
                names = cols[1:]
                for n in names:
                    counts[n]={}
            else:
                for i, count in enumerate(cols[1:]):
                    if valType=="float":
                        counts[names[i]][cols[0]] = float(count)
                    elif valType=="str":
                        counts[names[i]][cols[0]] = str(count)
                    elif valType=="int":
                        counts[names[i]][cols[0]] = int(count)
    return counts

def readMatrix(countFile, delim="\t"):
    cmat={}
    with open(countFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols=line.rstrip("\n").split(delim)
            if lc == 1:
                names = cols[1:]
                for n in names:
                    cmat[n]=[]
            else:
                for i, count in enumerate(cols[1:]):
                    cmat[names[i]].append(float(count))
    return cmat

def readMatrixFlip(countFile, delim="\t"):
    cmat={}
    with open(countFile, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols=line.rstrip("\n").split(delim)
            if lc > 1:
                n = cols[0]
                cmat[n]=[float(c) for c in cols[1:]]
    return cmat


def corrMat(matFile, outfile, cMap="viridis", setdpi=300, meth='pearson'):
    pyMat = readMatrix(matFile)
    df = pd.DataFrame(pyMat,columns=sorted(list(pyMat.keys())))
    corrMatrix = df.corr(method=meth)
    
    height = 0.3*len(pyMat)
    width = 0.4*len(pyMat)
    fig, ax = plt.subplots(figsize=(width,height))
    sn.heatmap(corrMatrix, annot=False, ax=ax, cmap=cMap)
    fig.savefig(outfile,dpi=setdpi,bbox_inches='tight')
    