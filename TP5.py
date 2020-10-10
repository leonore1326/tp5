#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 14:21:46 2020

@author: leonore
"""
import numpy as np
import Levenshtein as lv

## QUESTION 1

DNA_samples = ['ACCATACCTTCGATTGTCGTGGCCACCCTCGGATTACACGGCAGAGGTGC',
               'GTTGTGTTCCGATAGGCCAGCATATTATCCTAAGGCGTTACCCCAATCGA',
               'TTTTCCGTCGGATTTGCTATAGCCCCTGAACGCTACATGCACGAAACCAC',
               'AGTTATGTATGCACGTCATCAATAGGACATAGCCTTGTAGTTAACAG',
               'TGTAGCCCGGCCGTACAGTAGAGCCTTCACCGGCATTCTGTTTG',
               'ATTAAGTTATTTCTATTACAGCAAAACGATCATATGCAGATCCGCAGTGCGCT',
               'GGTAGAGACACGTCCACCTAAAAAAGTGA',
               'ATGATTATCATGAGTGCCCGGCTGCTCTGTAATAGGGACCCGTTATGGTCGTGTTCGATCAGAGCGCTCTA',
               'TACGAGCAGTCGTATGCTTTCTCGAATTCCGTGCGGTTAAGCGTGACAGA',
               'TCCCAGTGCACAAAACGTGATGGCAGTCCATGCGATCATACGCAAT',
               'GGTCTCCAGACACCGGCGCACCAGTTTTCACGCCGAAAGCATC',
               'AGAAGGATAACGAGGAGCACAAATGAGAGTGTTTGAACTGGACCTGTAGTTTCTCTG',
               'ACGAAGAAACCCACCTTGAGCTGTTGCGTTGTTGCGCTGCCTAGATGCAGTGG',
               'TAACTGCGCCAAAACGTCTTCCAATCCCCTTATCCAATTTAACTCACCGC',
               'AATTCTTACAATTTAGACCCTAATATCACATCATTAGACACTAATTGCCT',
               'TCTGCCAAAATTCTGTCCACAAGCGTTTTAGTTCGCCCCAGTAAAGTTGT',
               'TCAATAACGACCACCAAATCCGCATGTTACGGGACTTCTTATTAATTCTA',
               'TTTTTCGTGGGGAGCAGCGGATCTTAATGGATGGCGCCAGGTGGTATGGA']

def hamming_distance(S1,S2):
    d=0
    for k in range (len(S1)):
        if S1[k]!=S2[k]:
            d+=1
    return(d)

distance1=hamming_distance(DNA_samples[0], DNA_samples[1])
print(distance1)

##QUESTION 2

def levenshtein_distance(S1,S2):
    M=np.zeros((len(S1)+1,len(S2)+1))
    for i in range(1,len(S1)+1):
        M[i,0]=i
    for j in range(1,len(S2)+1):
        M[0,j]=j
    for j in range (1,len(S2)+1):
        for i in range(1,len(S1)+1):
            insertCost=M[i-1][j]+1
            deleteCost=M[i][j-1]+1
            if S1[i-1]==S2[j-1]:
                subCost=M[i-1][j-1]
            else:
                subCost=M[i-1][j-1]+1
            M[i][j]=min(insertCost,deleteCost,subCost)
    return(M,M[-1][-1])


print(levenshtein_distance('kryptonite','python'))

# tableau des distances entre séquences d'ADN

M_leo=[]
M_str=[]
M_module=[]
for k in range (len(DNA_samples)):
    M_leo.append([])
    M_module.append([])
    M_str.append([])
    for m in range (len(DNA_samples)):
        trace,d=levenshtein_distance(DNA_samples[k], DNA_samples[m])
        M_leo[k].append(d)
        M_str[k].append(str(int(d)))
        M_module[k].append(lv.distance(DNA_samples[k], DNA_samples[m]))

np.array(M_leo)
M_module=np.array(M_module)

print(M_leo, M_module, sep='\n')

# création du fichier

fichier=open("/Users/leonore/Desktop/Programmation/matrix.txt","w")
for i in range (len(M_str)):
    ligne=' '.join(M_str[i])+'\n'
    fichier.write(ligne)
fichier.close()

### QUESTION 3 et 4

def smith_waterman(a, b, alignment_score=2, gap_cost=1):
  # """
  # Compute the Smith-Waterman alignment score for two strings.
  # See https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Algorithm
  # This implementation has a fixed gap cost (i.e. extending a gap is considered
  # free). In the terminology of the Wikipedia description, W_k = {c, c, c, ...}.
  # This implementation also has a fixed alignment score, awarded if the relevant
  # characters are equal.
  # Kinda slow, especially for large (50+ char) inputs.
  # """
  # H holds the alignment score at each point, computed incrementally
  H = np.zeros((len(a) + 1, len(b) + 1))
  for i in range(1, len(a) + 1):
    for j in range(1, len(b) + 1):
      # The score for substituting the letter a[i-1] for b[j-1]. Generally low
      # for mismatch, high for match.
      if a[i-1] == b[j-1]:       
          match = H[i-1,j-1] + alignment_score 
      else:
          match = H[i-1,j-1] - alignment_score
    
      # The scores for for introducing extra letters in one of the strings (or
      # by symmetry, deleting them from the other).
      if i>1:         
          delete = H[i-1,j] - gap_cost 
      else:
          delete = 0
          
      if j>1:          
          insert = H[i,j-1] - gap_cost 
      else :
          insert= 0
      H[i,j] = max(match, delete, insert, 0)
      
  # The highest score is the best local alignment.
  # For our purposes, we don't actually care _what_ the alignment was, just how
  # aligned the two strings were.
  return (H, H.max())

print(smith_waterman('actagag','gacatat'))           

###   QUESTION 5

def smith_waterman_visu(a, b, alignment_score=2, gap_cost=1):
  # """
  # Compute the Smith-Waterman alignment score for two strings.
  # See https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Algorithm
  # This implementation has a fixed gap cost (i.e. extending a gap is considered
  # free). In the terminology of the Wikipedia description, W_k = {c, c, c, ...}.
  # This implementation also has a fixed alignment score, awarded if the relevant
  # characters are equal.
  # Kinda slow, especially for large (50+ char) inputs.
  # """
  # H holds the alignment score at each point, computed incrementally
  H = np.zeros((len(a) + 1, len(b) + 1))
  V=[[None for i in range(len(b)+1)]for i in range(len(a)+1)]  
  for i in range(1, len(a) + 1):
    for j in range(1, len(b) + 1):
      # The score for substituting the letter a[i-1] for b[j-1]. Generally low
      # for mismatch, high for match.
      if a[i-1] == b[j-1]:       
          match = H[i-1,j-1] + alignment_score 
      else:
          match = H[i-1,j-1] - alignment_score
    
      # The scores for for introducing extra letters in one of the strings (or
      # by symmetry, deleting them from the other).
      if i>1:         
          delete = H[i-1,j] - gap_cost 
      else:
          delete = 0
          
      if j>1:          
          insert = H[i,j-1] - gap_cost 
      else :
          insert= 0
      H[i,j] = max(match, delete, insert, 0)
#creation de la matrice trace
      if max(match, delete, insert, 0)==match:
          V[i][j]="&"
      elif max(match, delete, insert, 0)==delete:
          V[i][j]="|"
      elif max(match, delete, insert, 0)==insert:
          V[i][j]="-"
          
  indmax=np.where(H==H.max())
  Ord=indmax[0][0]          
  Abs=indmax[1][0]
#remonter la trace : 
  Aff1=[]
  Aff2=[]
  while V[Ord][Abs]!=None:
      #print(Abs,Ord)
      obis=Ord
      abis=Abs
      if V[Ord][Abs]=="&":
          Aff1.append(b[Abs-1])
          Aff2.append(a[Ord-1])
          abis-=1
          obis-=1
      if V[Ord][Abs]=="|":
          obis-=1
      if V[Ord][Abs]=="-":
          Aff1.append(b[Abs-1])
          Aff2.append("-")
          abis-=1
      Ord=obis
      Abs=abis
  Aff1.reverse()
  Aff2.reverse()
  middle=[]
  for i in range(len(Aff2)):
      if Aff2[i]!='-':
          middle.append('| ')
      else:
          middle.append('  ')
  
  str1=' '.join(Aff2)
  str2=''.join(middle)
  str3=' '.join(Aff1)   
  print('\n',str1,str2,str3,sep='\n')
      
              
  
smith_waterman_visu('actagag','gacatat')

