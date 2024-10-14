import sys
import argparse
from pafpy import PafFile
import os 


parser = argparse.ArgumentParser()
parser.add_argument("--paf", required=True ,help="alignment paf of reference unique sequence to query haplotype.")
parser.add_argument("--bed", required=True ,help="output bed file.")
parser.add_argument("--sample", required=True ,help="sample name passed as wildcard.")
parser.add_argument("--hap", required=True ,help="hap name passed as wildcard.")
parser.add_argument("--rgn", required=True ,help="rgn name passed as wildcard.")

# Process PAF file with pysam
args = parser.parse_args()
path =  args.paf
out_bed = args.bed
print(path)
samp = args.sample
hap = args.hap
rgn = args.rgn



records = list()
with PafFile(path) as paf:
    for record in paf:
        records.append(record)

def write_output(out_line):
    '''write output to bed file - whether actual bed or failure'''
    with open(out_bed, 'w') as file:  # Use 'a' to append to the file instead of overwriting it
        file.write(f'{out_line}\n')  # Write a line to the file

if(len(records) != 2):
    write_output(f"Fail: {rgn}:{samp}_{hap}: more than 2 alignments")
    exit()

if( len( set( rec.tname for rec in records ) )!= 1):
    write_output(f"Fail: {rgn}:{samp}_{hap}: unique regions aligned to different contigs")
    exit()

if( len( set( str(rec.strand) for rec in records ) )!= 1):
    write_output(f"Fail: {rgn}:{samp}_{hap} : unique regions aligned on different strands") 
    exit()

rgn_start = min([rec.tstart for rec in records])
rgn_end = max([rec.tend for rec in records])
rgn_strand = str(records[0].strand)
rgn_tig = records[0].tname

out_string = f"{rgn_tig}\t{rgn_start}\t{rgn_end}\t{rgn}_{samp}_{hap}\t.\t{rgn_strand}"
write_output( out_string )
