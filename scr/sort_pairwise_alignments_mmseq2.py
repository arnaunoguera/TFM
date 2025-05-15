#!/usr/bin/python

import sys
start = 0
data = { }

while True:
    l = sys.stdin.readline() # read each line from stdin
    if start == 1:
        if len(l)==0:
            break
        query,target,fident,qstart,qend,tstart,tend,lenQ,lenT,cigar,evalue,covQ,covT = l.split() # get the columns
        data.setdefault((target,query), [ ]).append([qstart, qend, fident, lenQ, covQ]) #use as dict key the combination of ref and query
    else:
        if l.startswith("query"):
            start = 1
            continue

# Print the header
print("""# C : Proportion of query positions covered
# I_local : Local sequence identity
# I_global : Global sequence identity
Reference_Sequence\tQuery_Sequence\tC\tI_local\tI_global""")

for (ref, query), values in sorted(data.items(), key=lambda item: item[0][0]):
    if ref == query: # if the reference and query are the same, print 1.0 for all values
        print(f"{ref}\t{query}\t1.0\t1.0\t1.0")
    elif len(values) == 1:
        s2, e2, ident, lenQ, covQ = values[0]
        I_local = round(float(ident), 4)
        c = round(float(covQ), 4)
        I_global = round(c * I_local, 4)
        print(f"{ref}\t{query}\t{c}\t{I_local}\t{I_global}")
    else:
        # Get the summary of the alignment if there are multiple local alignments
        lenQ = int(values[0][3])

        #first, check if start (s2) and end (e2) coordinates overlap
        # get a list of all start and end coordinates
        intervals = [(min(int(v[0]), int(v[1])), max(int(v[0]), int(v[1]))) for v in values]
        # Sort the intervals by the start coordinate
        intervals.sort()

        # Merge overlapping intervals and calculate the total length
        merged_intervals = []
        current_start, current_end = intervals[0]
        total_length = 0
        sum_block_lengths = current_end - current_start + 1

        for start, end in intervals[1:]:
            sum_block_lengths += end - start + 1
            if start <= current_end:  # Overlapping intervals
                current_end = max(current_end, end)
            else:  # Non-overlapping interval
                merged_intervals.append((current_start, current_end))
                total_length += current_end - current_start + 1
                current_start, current_end = start, end

        # Add the last interval
        merged_intervals.append((current_start, current_end))
        total_length += current_end - current_start + 1

        # Calculate the total coverage
        total_c = total_length/lenQ

        # Calculate total local identity
        # For this, calculate ifrst the number of mismatches / gaps in the query alignment: (1 - identity) * (qstart - qend)
        delta = round(sum([(1 - float(v[2])) * (max(int(v[0]), int(v[1])) - min(int(v[0]), int(v[1]))) for v in values]))
        
        I_local = (sum_block_lengths - delta) / sum_block_lengths

        I_global = total_c * I_local

        # Round the results
        total_c = round(total_c, 4)
        I_local = round(I_local, 4)
        I_global = round(I_global, 4)
        print(f"{ref}\t{query}\t{total_c}\t{I_local}\t{I_global}")