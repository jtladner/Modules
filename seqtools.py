from collections import defaultdict

import kmertools as kt

std_nt = {"A":2, "T":2, "C":3, "G":3, "U":3}
std_dna = {"A":2, "T":2, "C":3, "G":3}
comp_dna = {"A":"T", "T":"A", "C":"G", "G":"C", "Y":"R", "R":"Y", "W":"W", "S":"S", "M":"K", "K":"M", "B":"V", "V":"B", "D":"H", "H":"D", "N":"N"}
std_rna = {"A":2, "C":3, "G":3, "U":3}
ambig_nt = {"Y":{"T":2, "C":3}, 
			"R":{"A":2, "G":3},
			"W":{"A":2, "T":2},
			"S":{"C":3, "G":3},
			"M":{"A":2, "C":3},
			"K":{"T":2, "G":3},
			"B":{"C":3, "G":3, "T":2},
			"D":{"A":2, "G":3, "T":2},
			"H":{"A":2, "C":3, "T":2},
			"V":{"A":2, "C":3, "G":3},
			"N":{"A":2, "C":3, "G":3, "T":2},
			}

pair2ambig = {("C", "T"):"Y", 
			  ("A", "G"):"R",
			  ("A", "T"):"W",
			  ("C", "G"):"S",
			  ("A", "C"):"M",
			  ("G", "T"):"K",
			}

trip2ambig = {("C", "G", "T"):"B", 
			  ("A", "G", "T"):"D",
			  ("A", "C", "G"):"V",
			  ("A", "C", "T"):"H",
			}


# Dictionary linking a DNA codon to the 1-letter amino acid to which it corresponds
codon_table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
		'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
	}

# Dictionary linking a PHRED character to the integer quality score, +33
phred={chr(x+33):x for x in range(0,94)}

# Calculate GC content
def calcGC(seq):
	seq = seq.upper()
	countD = {nt:seq.count(nt) for nt in std_dna}
	return ((countD["G"]+countD["C"])/sum(countD.values()))*100

# Reverse complement a DNA sequence
def revCompDNA(seq):
	revSeq = seq[::-1].upper()
	revComp = [comp_dna[b] for b in revSeq]
	return "".join(revComp)

def mutType(codon1, codon2):
	
	aa1 = codon_table[codon1]
	aa2 = codon_table[codon2]
	
	if aa1 == aa2:
		return "Synonymous"
	elif aa2 == "*":
		return "Nonsense"
	elif aa1 == "*":
		return "Extend ORF"
	else:
		return "Non-synonymous"

def propMutTypes(codingSeq):
	mutCounts = defaultdict(int)

	for i in range(0,len(codingSeq)-2, 3):
		this = codingSeq[i:i+3]
		poss = possCodonMuts(this)
		for p in poss:
			mutCounts[mutType(this,p)] += 1
	
	total = sum(mutCounts.values())
	return {k:v/total for k,v in mutCounts.items()}
	
def possCodonMuts(codon):
	poss = []
	for i,b in enumerate(codon):
		other = [x for x in std_dna if x != b]
		for each in other:
			poss.append(codon[:i] + each + codon[i+1:])
	return poss

# Calculate % identity for two aligned nucleotide sequences
def percID_nt(s1, s2, indels=False, skipN=False):
	if len(s1) != len(s2):
		print("Sequences must be the same length to run percID_nt()")
		return False
	
	else:
		total=0
		same=0
	
		#Convert both seqs to upper case
		s1 = s1.upper()
		s2 = s2.upper()
	
		gap1 = False
		gap2 = False
	
		for i, b1 in enumerate(s1):
			b2 = s2[i]
			# If seq 1 has a gap
			if b1=="-":
				# If both sequences have gaps
				if b2 == "-":
					continue
				else:
					gap2=False
			
				if indels and not gap1:
					total+=1
					gap1=True
					continue
		
			else:
				gap1=False
			
				if b2 == "-":
					if indels and not gap2:
						total+=1
						gap2=True
						continue
				else:
					gap2=False
				
					total+=1

					if not skipN or (b1 != "N" and b2 != "N"):

						#If the positions are the same
						if b1 == b2:
							same+=1
							continue
						#If both positions are non-ambiguous
						elif b1 in std_nt and b2 in std_nt:
							continue
						else:				
							if b1 in ambig_nt:
								b1 = set(ambig_nt[b1].keys())
							else:
								b1 = set([b1])

							if b2 in ambig_nt:
								b2 = set(ambig_nt[b2].keys())
							else:
								b2 = set([b2])

							if len(b1.intersection(b2)) > 0:
								same+=1
		if total:
			return same/total*100
		else:
			return False

# Calculate % identity for two aligned amino acid sequences
def percID_aa(s1, s2, indels=False, skipX=False):
	if len(s1) != len(s2):
		print("Sequences must be the same length to run percID_aa()")
		return False
	
	else:
		total=0
		same=0
	
		#Convert both seqs to upper case
		s1 = s1.upper()
		s2 = s2.upper()
	
		gap1 = False
		gap2 = False
	
		for i, b1 in enumerate(s1):
			b2 = s2[i]
			# If seq 1 has a gap
			if b1=="-":
				# If both sequences have gaps
				if b2 == "-":
					continue
				else:
					gap2=False
			
				if indels and not gap1:
					total+=1
					gap1=True
					continue
		
			else:
				gap1=False
			
				if b2 == "-":
					if indels and not gap2:
						total+=1
						gap2=True
						continue
				else:
					gap2=False
				
					total+=1
					
					if not skipX or (b1 != "X" and b2 != "X"):
						#If the positions are the same
						if b1 == b2:
							same+=1
							continue

		if total:
			return same/total*100
		else:
			return False

		
def rmv_common_gaps(*argv):
	toFill = [""]*len(argv)
	for i in range(len(argv[0])):
		basesPresent = [s[i] for s in argv]
		uniqBases = set(basesPresent)
		if uniqBases != set(["-"]):
			for j,b in enumerate(basesPresent):
				toFill[j]+=b
	return toFill

# Output coordinates are 0-indexed and relative to seq1
def findIndels(seq1, seq2):
	insD = {}
	delD = {}
	
	currDel = 0
	currIns = 0
	
	refPos = -1
	
	for i, b1 in enumerate(seq1):
		b2 = seq2[i]
		
		# Deletion
		if b1 != "-" and b2 == "-":
			refPos+=1
			if currIns:
				insD[refPos-currIns] = currIns
				currIns = 0
			currDel+=1

		# Insertion
		elif b1 == "-" and b2 != "-":
			if currDel:
				delD[refPos-currDel] = currDel
				currDel = 0
			currIns+=1

		else:
			refPos+=1
			if currIns:
				insD[refPos-currIns] = currIns
				currIns = 0

			if currDel:
				delD[refPos-currDel] = currDel
				currDel = 0

	# In case there's an indel at the very end of the alignment
	if currIns:
		insD[i-currIns+1] = currIns

	if currDel:
		delD[i-currDel+1] = currDel
		
	return insD, delD

def findSNPs(seq1, seq2):
	snpD = {}
	
	for i, b1 in enumerate(seq1):
		b2 = seq2[i]
		
		# If not an indel
		if b1 != "-" and b2 != "-":
			# If different
			if b1 != b2:
				snpD[i] = b2

	return snpD

# Remove insertions in seq2 relative to seq1
def rmv_ins(seq1, seq2):
	out1 = ""
	out2 = ""
	for i, b1 in enumerate(seq1):
		if b1 != "-":
			out1+=b1
			out2+=seq2[i]
	return out1, out2
	
# Use a sliding window approach to trim poor qulaity sequence from the beginning and/or end of sequences
def sw_qual_trim(seqD, qualD, avgQ = 25, winS=5, beg=True, end=True):
	
	for name in list(qualD.keys()):
		#If trimming poor quality from the beginning of the seqs
		if beg:
			begTrim = 0
			for i in range(0,len(qualD[name])-winS+1, 1):
				these = [phred[q] for q in qualD[name][i:i+winS]]
#				print(these)
				if sum(these)/len(these) < avgQ:
					begTrim = i+winS
				else:
					break
			if begTrim:
#				print(name, "beg", begTrim)
				seqD[name] = seqD[name][begTrim:]
				qualD[name] = qualD[name][begTrim:]
			
		#If trimming poor quality from the end of the seqs
		if end:
			revQuals = qualD[name][::-1]
			endTrim = 0
			for i in range(0,len(revQuals)-winS+1, 1):
				these = [phred[q] for q in revQuals[i:i+winS]]
				if sum(these)/len(these) < avgQ:
					endTrim = i+winS
				else:
					break
			if endTrim:
#				print(name, "end", endTrim)
				seqD[name] = seqD[name][:-endTrim]
				qualD[name] = qualD[name][:-endTrim]
					
	return seqD, qualD

# Merge forward and reverse Sanger sequences, if they overlap by a kmer of at least a certain size
def merge_FandR(forS, revS, minK=15):
	k = min([len(forS), len(revS)])
	while k >= minK:
		kf = kt.kmerSet(forS, k, filter=["N"])
		kr = kt.kmerSet(revS, k, filter=["N"])
		ovlp = kf.intersection(kr)
		if len(ovlp)>0:
#			print(k)
			firstKmer = list(ovlp)[0]
			indF = forS.index(firstKmer)
			indR = revS.index(firstKmer)
			break
		else:
			k-=1
	
	if k<minK:
		print("No match found.")
	else:
		aligned= []
		if indF>indR:
			aligned.append(forS)
			aligned.append(" "*(indF-indR) + revS)
		elif indF<indR:
			aligned.append(" "*(indR-indF) + forS)
			aligned.append(revS)
		else:
			aligned.append(forS)
			aligned.append(revS)
		
		# If the sequences aren't the same length
		if len(aligned[0]) < len(aligned[1]):
			aligned[0] = aligned[0] + " "*(len(aligned[1]) - len(aligned[0]))
		elif len(aligned[1]) < len(aligned[0]):
			aligned[1] = aligned[1] + " "*(len(aligned[0]) - len(aligned[1]))
	
	return consensus_2seqs(aligned), aligned

def consensus_2seqs(seqL):
	cons = ""
	for i in range(len(seqL[0])):
		bases = [s[i] for s in seqL if s[i] != " "]
		if len(set(bases)) == 1:
			cons+=bases[0]
		else:
			cons+=pair2ambig[tuple(sorted(bases))]
	return cons

def consensus(seqL):
	cons = ""
	for i in range(len(seqL[0])):
		bases = [s[i] for s in seqL if s[i] != " " and s[i] != "-"]
		uniqBases = set(bases)
		if len(uniqBases) == 1:
			cons+=bases[0]
		elif len(uniqBases) > 1:
			countD = {b:bases.count(b) for b in uniqBases}
			countDrev = defaultdict(list)
			for k,v in countD.items():
				countDrev[v].append(k)
			
			mostCommon = max(countDrev.keys())
			#If there is one base with >50% frequency
			if mostCommon/len(bases) >0.5:
				cons+=countDrev[mostCommon][0]
			elif len(countDrev[mostCommon]) == 2:
				cons+=pair2ambig[tuple(sorted(countDrev[mostCommon]))]
			elif len(countDrev[mostCommon]) == 3:
				cons+=trip2ambig[tuple(sorted(countDrev[mostCommon]))]
			else:
				cons+="N"
	return cons


