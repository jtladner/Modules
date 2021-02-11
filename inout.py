from collections import defaultdict

def fileDictHeader(file, key, val, delim="\t", splitVal=False, splitKey=False, valType="str"):
    fileD={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols = line.rstrip("\n").split(delim)
            if lc==1:
                headD={}
                for i,x in enumerate(cols):
                    headD[x] = i
            elif cols != ['']:
                if splitKey:
                    for each in cols[headD[key]].split(splitKey):
                        if splitVal:
                            if valType=="float":
                                fileD[each] = float(cols[headD[val].split(splitDelim)])
                            elif valType=="str":
                                fileD[each] = str(cols[headD[val].split(splitDelim)])
                            elif valType=="int":
                                fileD[each] = int(cols[headD[val].split(splitDelim)])

                        else:
                            if valType=="float":
                                fileD[each] = float(cols[headD[val]])
                            elif valType=="str":
                                fileD[each] = str(cols[headD[val]])
                            elif valType=="int":
                                fileD[each] = int(cols[headD[val]])
                
                else: 
                    if splitVal:
                        if valType=="float":
                            fileD[cols[headD[key]]] = float(cols[headD[val].split(splitDelim)])
                        elif valType=="str":
                            fileD[cols[headD[key]]] = str(cols[headD[val].split(splitDelim)])
                        elif valType=="int":
                            fileD[cols[headD[key]]] = int(cols[headD[val].split(splitDelim)])
                    else:
                        if valType=="float":
                            fileD[cols[headD[key]]] = float(cols[headD[val]])
                        elif valType=="str":
                            fileD[cols[headD[key]]] = str(cols[headD[val]])
                        elif valType=="int":
                            fileD[cols[headD[key]]] = int(cols[headD[val]])

    return fileD

def fileDictHeaderLists(file, key, val, delim="\t"):
    fileD=defaultdict(list)
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            cols = line.rstrip("\n").split(delim)
            if lc==1:
                headD={}
                for i,x in enumerate(cols):
                    headD[x] = i
            elif cols != ['']:
                fileD[cols[headD[key]]].append(cols[headD[val]])
    return fileD


def fileList(file, col=0, delim="\t", header=True):
    l=[]
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                l.append(thisCols[col])
    return l

def fileListHeader(file, name, delim="\t", valType="str"):
    l=[]
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            thisCols = line.rstrip("\n").split(delim)
            if lc == 1:
                col = thisCols.index(name)
            else:
                if valType=="float":
                    l.append(float(thisCols[col]))
                elif valType=="str":
                    l.append(str(thisCols[col]))
                elif valType=="int":
                    l.append(int(thisCols[col]))

                
    return l


def fileEmptyDict(file, col=0, delim="\t", header=True):
    l={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                l[thisCols[col]]=""
    return l
    
def fileDict(file, key=0, val=1, delim="\t", header=True):
    l={}
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            if lc>1 or header==False:
                thisCols = line.rstrip("\n").split(delim)
                try: l[thisCols[key]]=thisCols[val]
                except: print(thisCols)
    return l

def fileDictFull(file, delim="\t"):
    l=defaultdict(list)
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            thisCols = line.rstrip("\n").split(delim)
            if lc==1:
                headMap = {i:x for i,x in enumerate(thisCols)}
            else:
                for i,x in enumerate(thisCols):
                    l[headMap[i]].append(x)
    return l

def fileDictFullRowNames(file, delim="\t"):
    l=defaultdict(dict)
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            thisCols = line.rstrip("\n").split(delim)
            if lc==1:
                headMap = {i:x for i,x in enumerate(thisCols)}
            else:
                for i,x in enumerate(thisCols):
                    if i>0:
                        l[headMap[i]][thisCols[0]] = x
    return l


def writeList(lst, outname, delim="\n"):
    with open(outname, "w") as fout:
        fout.write(delim.join(lst))
