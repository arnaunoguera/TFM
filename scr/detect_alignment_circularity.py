#!/usr/bin/python

import sys
start = 0
data = { }
circular = 0

while True:
    l = sys.stdin.readline() # read each line from stdin
    if len(l)==0:
        break
    query,target,fident,qstart,qend,tstart,tend,lenQ,lenT,cigar,evalue,covQ,covT = l.split() # get the columns
    data.setdefault((target,query), [ ]).append([qstart, qend, tstart, tend, fident, lenQ, covQ]) #use as dict key the combination of ref and query

for (ref, query), values in sorted(data.items(), key=lambda item: item[0][0]):

    if ref == query: # if the reference and query are the same, next
        continue
    elif len(values) == 1: # if there is only one local alignment, next
        continue
    else:
        lenQ = int(values[0][5])
        # I will make an object like: [(qmin, qmax, tmin, tmax, identity), (...), ...]
        # (making sure to order qstart and qend)
        alignments = [(min(int(v[0]), int(v[1])), max(int(v[0]), int(v[1])), min(int(v[2]), int(v[3])), max(int(v[2]), int(v[3])), v[4]) for v in values]

        # Keep only the alignemtns where qmin is 1 or qmax is lenQ
        alignments = [a for a in alignments if a[0] == 1 or a[1] == lenQ]
        # if only one alignment is left (or there are more than 2, which shouldn't be the case), next
        if len(alignments) != 2:
            continue

        # Ensure that the identity is at least 0.9
        if float(alignments[0][4]) < 0.9 or float(alignments[1][4]) < 0.9:
            continue

        # Now sort the two alignments by qmin
        alignments.sort(key=lambda x: x[0])


        # Check if the two alignments are adjacent: tmin of the first alignment = tmax of the second + 1 (or vice versa)
        if alignments[0][2] == alignments[1][3] + 1 or alignments[0][3] + 1 == alignments[1][2]:
            circular = 1
            break # if the two alignments are adjacent, the plasmid is circular, we can stop iterating 

print(circular) # print the result