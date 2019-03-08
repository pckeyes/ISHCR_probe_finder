#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:51:59 2019

@author: piperkeyes
"""
import copy

#Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)

def main():
    gene = input("Enter your gene of interest (enclosed in double quotes): ").upper()
    gene_name = input ("Enter the gene's name (enclosed in double quotes): ")
    gene_len = len(gene)
    original_pairs = list()
    rev_comp_pairs = list()
    accepted_pairs = list()
    not_accepted_pairs = list()
    while gene_len >= 52:
        strand1 = gene[0:25]
        spacer = gene[25:27]
        strand2 = gene[27:52]
        (tm1, tm2) = getMeltingTemp(strand1, strand2)
        print("passed melting temp")
        if tm1 < 65 and tm1 > 55 and tm2 < 65 and tm2 > 55:
            rev_comp_pair = getReverseComplement(strand1, strand2, original_pairs, rev_comp_pairs)
            print("passed getting reverse complement")
            accepted = acceptOrRejectPair(rev_comp_pair, accepted_pairs, not_accepted_pairs)
            print("passed accepting or rejecting")
            if (accepted): gene = gene[79:]
            else: gene = gene[1:]
        else: gene = gene[1:]
        gene_len = len(gene)
        print(gene_len)
    writeToFile(gene_name, accepted_pairs)

def getMeltingTemp(strand1, strand2):
    tm1 = 0
    tm2 = 0
    A1 = 0
    T1 = 0
    G1 = 0
    C1 = 0
    A2 = 0
    T2 = 0
    G2 = 0
    C2 = 0
    #calculate tm of strand1
    for i in strand1:
        if (i == "A"): A1 = A1 + 1
        elif (i == "T"): T1 = T1 + 1
        elif (i == "G"): G1 = G1 + 1
        else: C1 = C1 + 1 
    tm1 = 64.9 + 41 * (G1 + C1 - 16.4)/(A1+T1+G1+C1)
    #calculate tm of strand2
    for i in strand2:
        if (i == "A"): A2 = A2 + 1
        elif (i == "T"): T2 = T2 + 1
        elif (i == "G"): G2 = G2 + 1
        else: C2 = C2 + 1 
    tm2 = 64.9 + 41 * (G2 + C2 - 16.4)/(A2+T2+G2+C2)
    return(tm1, tm2)

def getReverseComplement(strand1, strand2, original_pairs, rev_comp_pairs):
    original_pair = (strand1, strand2)
    original_pairs.append(original_pair)
    rev_comp1 = ""
    rev_comp2 = ""
    #generare first reverse complement
    for i in strand1:
        if (i == "A"): rev_comp1 = "T" + rev_comp1
        elif (i == "T"): rev_comp1 = "A" + rev_comp1
        elif (i == "G"): rev_comp1 = "C" + rev_comp1
        else: rev_comp1 = "G" + rev_comp1
    #generate second reverse complement
    for i in strand2:
        if (i == "A"): rev_comp2 = "T" + rev_comp2
        elif (i == "T"): rev_comp2 = "A" + rev_comp2
        elif (i == "G"): rev_comp2 = "C" + rev_comp2
        else: rev_comp2 = "G" + rev_comp2
    rev_comp_pairs.append((rev_comp1, rev_comp2))
    return (rev_comp1, rev_comp2)
    
def acceptOrRejectPair(rev_comp_pair, accepted_pairs, not_accepted_pairs):
    for strand in rev_comp_pair:
        #loops over all possible secondary folding structures with no base in the bulge
        for x in range(6,21):
            short_side = strand[:x]
            long_side = strand[x:]
            #make sure shorter side of the fold is assigned correctly
            if len(short_side) > len(long_side):
                temp = copy.deepcopy(short_side)
                short_side = copy.deepcopy(long_side)
                long_side = temp
            (consecutive_pairs, GC_pairs) = getSecondaryStructure(short_side, long_side, x)
            if GC_pairs >= 3 or consecutive_pairs >=4:
                not_accepted_pairs.append(rev_comp_pair)
                return False
            
        #loops over all possible secondary folding structures with one base in the bulge
        for x in range(6,21):
            short_side = strand[:x]
            long_side = strand[x + 1:]
            #make sure shorter side of the fold is assigned correctly
            if len(short_side) > len(long_side):
                temp = copy.deepcopy(short_side)
                short_side = copy.deepcopy(long_side)
                long_side = temp
            (consecutive_pairs, GC_pairs) = getSecondaryStructure(short_side, long_side, x)
            if GC_pairs >= 3 or consecutive_pairs >=4:
                not_accepted_pairs.append(rev_comp_pair)
                return False
    accepted_pairs.append(rev_comp_pair)
    return True

def getSecondaryStructure(short_side, long_side, x):
    consecutive_pairs = 0
    GC_pairs = 0
    secondary_pairs = list()
    #accumulate totals of secondary structure for this folding
    for c in range(len(short_side)):
        if short_side[c] == "A" and long_side[x-2-c] == "T" or short_side[c] == "T" and long_side[x-2-c] == "A":
            secondary_pairs.append("AT")
        elif short_side[c] == "G" and long_side[x-2-c] == "C" or short_side[c] == "C" and long_side[x-2-c] == "G":
            secondary_pairs.append("GC")
        else: 
            secondary_pairs.append("N/A")
    for y in range(1,len(secondary_pairs)):
        #find consecutive strong pairs in secondary structure
        if (secondary_pairs[y-1] != "N/A" and secondary_pairs [y] != "N/A"):
            consecutive_pairs = consecutive_pairs + 1
        else: consecutive_pairs = 0
        #find consecutive GC pairs
        if (secondary_pairs[y-1] == "GC" and secondary_pairs [y] == "GC"):
            GC_pairs = GC_pairs + 1
        else: GC_pairs = 0
    return (consecutive_pairs, GC_pairs)

def writeToFile(gene_name, accepted_pairs):
    with open(gene_name + "_pairs","w") as f:
        for line in accepted_pairs:
            strs=" ".join(str(x) for x in line)
            f.write(strs+"\n")
  
if __name__== "__main__":
  main()
        
            
#write to .txt file

