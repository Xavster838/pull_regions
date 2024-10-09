from pafpy import PafFile
import os 

# Process PAF file with pysam
for record in pysam.AlignmentFile("./tmp/exp1_GM19129_1_mapping.paf", 'r', format = 'paf'):
    print(f"Query: {record.query_name}, Target: {record.target_name}")

path = snakemake.input.paf

records = list()
with PafFile(path) as paf:
    for record in paf:
        records.append(record)

def write_output(out_line):
    '''write output to bed file - whether actual bed or failure'''
    with open(snakemake.output.bed, 'w') as file:  # Use 'a' to append to the file instead of overwriting it
        file.write(f'{out_line}\n')  # Write a line to the file

if(len(records) != 2):
    write_output(f"Fail: {snakemake.wildcards.rgn}_{snakemake.wildcards.samp}_{snakemake.wildcards.hap}:more than 2 alignments")
    exit()

if( len( set( rec.tname for rec in records ) )!= 1):
    write_output(f"Fail: {snakemake.wildcards.rgn}_{snakemake.wildcards.samp}_{snakemake.wildcards.hap} : unique regions aligned to different contigs")
    exit()

if( len( set( str(rec.strand) for rec in records ) )!= 1):
    write_output(f"Fail: {snakemake.wildcards.rgn}_{snakemake.wildcards.samp}_{snakemake.wildcards.hap} : unique regions aligned on different strands")
    exit()

rgn_start = min([rec.tstart for rec in records])
rgn_end = max([rec.tend for rec in records])
rgn_strand = str(records[0].strand)
rgn_tig = records[0].tname

out_string = f"{rgn_tig}\t{rgn_start}\t{rgn_end}\t{snakemake.wildcards.rgn}_{snakemake.wildcards.samp}_{snakemake.wildcards.hap}\t.\t{rgn_strand}\n"
write_output( out_string )