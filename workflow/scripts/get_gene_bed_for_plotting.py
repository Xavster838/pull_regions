import sys
import argparse
import os 
import pandas as pd

is_rev = lambda bed_line : bed_line[5] == "-" #check if loc_bed flipped

def get_rev_gene_bed(g_bed, l_bed):
    '''get flipped gene_bed to annotate flipped loc_bed for plotting'''
    l_end = l_bed[2][0]
    g_bed[1], g_bed[2] = l_end - g_bed[2] , l_end - g_bed[1]  #NOTE: flipping start and end of gene beds to accomodate orientation change
    g_bed[5] = g_bed[5].apply(lambda x: '-' if x == '+' else '+')
    assert (g_bed[1] < g_bed[2]).all(), "Not all values in column 'x' are greater than those in column 'y'."
    return(g_bed)

def get_gene_bed(g_bed, l_bed):
    '''return gene_bed with loc_bed start subtracted form each start and end to get coordinates for plotting'''
    g_bed[1] = [ int( x - l_bed.iloc[0,1]) for x in g_bed[1] ]  #subtract start of loc_bed
    g_bed[2] = [ int( x - l_bed.iloc[0,1]) for x in g_bed[2] ] 
    return(g_bed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--loc_bed", required=True ,help="bed file of specific region. output of rules.get_loc_bed")
    parser.add_argument("--loc_gene_bed", required=True ,help="bed file of gene_oi intersect with loc_bed. output of rules.get_loc_bed")
    parser.add_argument("--out_bed", required=True ,help="output bed file name")
    parser.add_argument("--sample", required=True ,help="sample name passed as wildcard.")
    parser.add_argument("--hap", required=True ,help="hap name passed as wildcard.")
    parser.add_argument("--rgn", required=True ,help="rgn name passed as wildcard.")

    #load parser vars
    args = parser.parse_args()
    loc_bed =  args.loc_bed
    gene_bed = args.loc_gene_bed
    out_bed_path = args.out_bed
    samp = args.sample
    hap = args.hap
    rgn = args.rgn

    #turn beds into dataframes
    loc_bed = pd.read_csv(loc_bed, sep = "\t", header = None )
    gene_bed = pd.read_csv(gene_bed, sep = "\t", header = None)

    #process gene_bed
    if is_rev(loc_bed.loc[0,:]):
        out_bed = get_rev_gene_bed(g_bed = gene_bed, l_bed = loc_bed)
    else:
        out_bed = get_gene_bed(g_bed = gene_bed, l_bed = loc_bed)

    out_bed.iloc[:,0] = [f"{rgn}_{samp}_{hap}__{x}" for x in out_bed.iloc[:,0] ] #add region, samp, hap info to name.

    out_bed.to_csv(out_bed_path, sep = "\t", header = False, index = False)

if __name__ == "__main__":
    main()
