#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:51:59 2019

@author: piperkeyes
"""

#Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)

gene = input("Enter your gene of interest (enclosed in double quotes): ")
gene_len = len(gene)
temps = list()
original_pairs = list()
rev_comp_pairs = list()



#generate pairs of primer dna
while gene_len >= 52:
    strand1 = gene[0:25]
    spacer = gene[26:28]
    strand2 = gene[29:54]
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
    #if we found a good pair
    if tm1 < 65 and tm1 > 55 and tm2 < 65 and tm2 > 55:
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
        rev_comp_pair = (rev_comp1, rev_comp2)
        rev_comp_pairs.append(rev_comp_pair)
        gene = gene[79:]
        gene_len = len(gene)
    #if we didn't find a good pair
    else:
        gene = gene [1:]
        gene_len = len(gene)

#write to .txt file
with open("pairs","w") as f:
    for line in rev_comp_pairs:
        strs=" ".join(str(x) for x in line)
        f.write(strs+"\n")
