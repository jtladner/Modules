from collections import defaultdict
import os

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

def fileDictHeaderLists(file, key, val, delim="\t", valType="str"):
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
				if valType=="float":
					fileD[cols[headD[key]]].append(float(cols[headD[val]]))
				elif valType=="str":
					fileD[cols[headD[key]]].append(str(cols[headD[val]]))
				elif valType=="int":
					fileD[cols[headD[key]]].append(int(cols[headD[val]]))
				
	return fileD

def fileHeaderDictDictLists(file, k1, k2, val, delim="\t", valType="str"):
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
				if cols[headD[k1]] not in fileD:
					fileD[cols[headD[k1]]] = {}
					fileD[cols[headD[k1]]][cols[headD[k2]]] = []
				elif cols[headD[k2]] not in fileD[cols[headD[k1]]]:
					fileD[cols[headD[k1]]][cols[headD[k2]]] = []
				
				if valType=="float":
					fileD[cols[headD[k1]]][cols[headD[k2]]].append(float(cols[headD[val]]))
				elif valType=="str":
					fileD[cols[headD[k1]]][cols[headD[k2]]].append(str(cols[headD[val]]))
				elif valType=="int":
					fileD[cols[headD[k1]]][cols[headD[k2]]].append(int(cols[headD[val]]))
				
	return fileD

def fileHeaderDictDict(file, k1, k2, val, delim="\t", valType="str"):
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
				if cols[headD[k1]] not in fileD:
					fileD[cols[headD[k1]]] = {}
				
				if valType=="float":
					fileD[cols[headD[k1]]][cols[headD[k2]]] = float(cols[headD[val]])
				elif valType=="str":
					fileD[cols[headD[k1]]][cols[headD[k2]]] = str(cols[headD[val]])
				elif valType=="int":
					fileD[cols[headD[k1]]][cols[headD[k2]]] = int(cols[headD[val]])
				
	return fileD


def fileDictLists(file, key=0, val=1, delim="\t", header=True):
	fileD=defaultdict(list)
	with open(file, "r") as fin:
		lc=0
		for line in fin:
			lc+=1
			rows = line.rstrip("\n").split(delim)
			if lc>1 or header==False:
				if rows != ['']:
					try:
						fileD[rows[key]].append(rows[val])
					except:
						print(rows)
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

def fileDictFull(file, delim="\t", valType="str", rowNames=False):
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
					if i>0 or not rowNames:
						if valType=="str":
							l[headMap[i]].append(x)
						elif valType=="int":
							l[headMap[i]].append(int(x))
						elif valType=="float":
							l[headMap[i]].append(float(x))
						else:
							print("%s is an invalid valType" % (valType))
	return l

def fileDictFullRowNames(file, delim="\t", valType="str"):
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
						if valType=="str":
							l[headMap[i]][thisCols[0]] = x
						elif valType=="int":
							l[headMap[i]][thisCols[0]] = int(x)
						elif valType=="float":
							l[headMap[i]][thisCols[0]] = float(x)

	return l


def writeList(lst, outname, delim="\n"):
	with open(outname, "w") as fout:
		fout.write(delim.join(lst))

def writeDict(dct, keyHeader, valHeader, outname, delim="\t"):
	with open(outname, "w") as fout:
		fout.write(f"{keyHeader}{delim}{valHeader}\n")
		for k, v in dct.items():
			fout.write(f"{k}{delim}{v}\n")


def makeDir(dirName):
	if os.path.isdir(dirName):
		print("Warning: the directory %s already exists!" % (dirName))
	else:
		os.mkdir(dirName)