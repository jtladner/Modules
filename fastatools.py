#!/usr/bin/env python

def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

def read_fasta_dict_upper(file):
    names, seqs = read_fasta_lists(file)
    seqs = [x.upper() for x in seqs]
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

#writes a new fasta file
def write_fasta_dict(fD, new_filename):
    fout=open(new_filename, 'w')
    for k,v in fD.items():
        fout.write(">%s\n%s\n" % (k, v))
    fout.close()


#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    if len(names) != len(seqs):
        print(f"Warning: the number of names ({len(names)}) and the number of sequences ({len(seqs)}) provided to write_fasta are NOT the same!")
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


def seqLenList(seqs):
    lens = [len(s) for s in seqs]
    return sorted(list(set(lens)))

def seqLenDict(fD):
    seqs = list(fD.values())
    lens = [len(s) for s in seqs]
    return sorted(list(set(lens)))

def combine_fastafiles(fileList, outname):
    allNames=[]
    allSeqs=[]
    
    for f in fileList:
        n,s = read_fasta_lists(f)
        # Makes sure that some sequences were found
        if n:
            allNames+=n
            allSeqs+=s
    
    write_fasta(allNames, allSeqs, outname)