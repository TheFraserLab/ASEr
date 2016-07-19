from __future__ import print_function
from pysam import Samfile
from argparse import ArgumentParser, FileType
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
from sys import stdout
from os import path
import ASEr.logme as lm
from math import log2
import pickle

try:
    from progressbar import ProgressBar as pbar
except ImportError:
    print("Not loading progress bar")
    pbar = lambda: lambda x: iter(x)

def get_phase(read, snps):
    phase = None
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos + 1 in snps:
            if phase == None:
                try:
                    # 1 if alternate, -1 if reference
                    phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0 # This SNP isn't in the dataset
            else:
                try:
                    new_phase = -1 + 2*snps[ref_pos + 1].index(read.seq[read_pos])
                except ValueError:
                    return 0
                if new_phase != phase:
                    return 0 # read seems misphased
    return phase



def get_snps(snpfile):
    snps = defaultdict(dict)
    if path.exists(path.join(path.dirname(snpfile), 'true_hets.tsv')):
        print("using true hets")
        true_hets = {tuple(line.strip().split()):True
                     for line in open(path.join(path.dirname(snpfile), 'true_hets.tsv'))
                    }
    else:
        true_hets = defaultdict(lambda x: True)
    if snpfile.endswith('.bed'):
        for line in open(snpfile):
            chrom, _, start, refalt = line.strip().split()
            if true_hets.get((chrom, start), True):
                snps[chrom][int(start)] = refalt.split('|')

    return snps

def get_gene_coords(gff_file, id_name, feature_type='exon'):
    if path.exists(gff_file + '.pkl'):
        return pickle.load(open(gff_file + '.pkl', 'rb'))
    gene_coords = defaultdict(lambda : [None, set()])
    for line in open(gff_file):
        chrom, _, feature, left, right, _, _, _, annot = (
                line.split('\t'))
        if feature != feature_type: continue
        annot = annot.strip().strip(';').split('; ')
        feature_id = 'MISSING'
        for a in annot:
            if a.startswith(id_name):
                feature_id = a.split()[1].strip('"').strip("'")
                break
        else:
            feature_id = 'MISSING'
            lm.log("Can't find {} in line: '{}'".format(
                id_name, line.strip()),
                level='warn'
                )
            continue
        gene_coords[feature_id][0] = chrom
        gene_coords[feature_id][1].add((int(left), int(right)))
    gene_coords_out = {}
    for entry in gene_coords:
        gene_coords_out[entry] = gene_coords[entry]
    pickle.dump(gene_coords_out, open(gff_file + '.pkl', 'wb'))
    return gene_coords

def get_ase_by_coords(chrom, coords, samfile, snp_dict):
    left_most = 1e99
    right_most = 0

    for left, right in coords:
        left_most = min(left, left_most)
        right_most = max(right, right_most)
    assert left_most < right_most

    read_results = Counter()
    phases = defaultdict(set)
    snps_on_chrom = snp_dict[chrom]

    for read in samfile.fetch(chrom, left_most, right_most, multiple_iterators=True):
        read_left = read.reference_start
        read_right = read.reference_end

        # Make sure both ends are in an exon for this gene
        left_hit = False
        right_hit = False
        for left, right in coords:
            if left <= read_left <= right:
                left_hit = True
            if left <= read_right <= right:
                right_hit = True
            if left_hit and right_hit:
                break
        else:
            # One of the ends falls outside of the gene, probably indicating an
            # exon of a different gene, but also possibly an expanded 5' or 3'
            # UTR
            read_results['missed exon boundaries'] += 1
            continue
        phase = get_phase(read, snps_on_chrom)
        phases[read.qname].add(phase)

    read_counts = Counter()

    for phase_set in phases.values():
        phase_set.discard(None)
        if len(phase_set) == 1:
            # Unambiguously phased
            read_counts[phase_set.pop()] += 1
        elif len(phase_set) == 0:
            read_results['no phasing information'] += 1
            read_counts[None] += 1
        else:
            read_results['discordant phase information'] += 1
            read_counts[0] += 1
    return read_counts



def pref_index(ref, alt):
    return (alt-ref)/(alt+ref)

def log2ase(ref, alt):
    return log2(alt/ref)

def ratio(ref, alt):
    return alt/ref

ase_fcns = {
        'pref_index': pref_index,
        'log2': log2ase,
        'ratio': ratio,
        }

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('snp_file')
    parser.add_argument('gff_file')
    parser.add_argument('reads', type=Samfile)
    parser.add_argument('--max-jobs', '-p', default=0, type=int,
            help=''
            )
    parser.add_argument('--id-name', '-n', default='gene_id', type=str)
    parser.add_argument('--outfile', '-o', default=stdout,
            type=FileType('w'),
            )
    parser.add_argument('--min-reads-per-gene', '-m', default=20, type=int)
    parser.add_argument('--ase-function', '-f', default='log2', type=str)
    parser.add_argument('--min-reads-per-allele', '-M', default=0, type=int)

    args = parser.parse_args()
    if args.ase_function not in ase_fcns:
        print("Unrecognized function: {}".format(args.ase_fcn))
        raise ValueError
    print(args)
    return args

if __name__ == "__main__":
    args = parse_args()
    snp_dict = get_snps(args.snp_file)
    gene_coords = get_gene_coords(args.gff_file, args.id_name)

    ase_vals = {}
    if False and args.max_jobs != 1:
        with Pool(args.max_jobs or cpu_count()) as pool:
            for gene in gene_coords:
                ase_vals[gene] = pool.apply_async(
                        get_ase_by_coords,
                        (
                            gene_coords[gene][0],
                            gene_coords[gene][1],
                            args.reads,
                            snp_dict,
                            ))

            prog = pbar()
            for gene in prog(ase_vals):
                ase_vals[gene] = ase_vals[gene].get()
            if 'finish' in dir(prog):
                prog.finish()
    else:
        prog = pbar()
        for gene in prog(gene_coords):
            ase_vals[gene] = get_ase_by_coords(
                    gene_coords[gene][0],
                    gene_coords[gene][1],
                    args.reads,
                    snp_dict
                    )
        if 'finish' in dir(prog):
            prog.finish()
    print('gene', 'chrom', 'ref_counts', 'alt_counts', 'no_ase_counts',
            'ambig_ase_counts', 'ase_value',
            file=args.outfile, sep='\t', end='\n'
            )
    for gene in sorted(ase_vals):
        avg = ase_vals[gene]
        # Not "average", "ase vals for gene"
        if (min(avg[1], avg[-1]) > args.min_reads_per_allele) and (avg[1] + avg[-1] > args.min_reads_per_gene):
            ase_fcn = ase_fcns[args.ase_function]
            ase_val = ase_fcn(avg[-1], avg[1])
        else:
            ase_val = 'NA'
        print(
                gene, gene_coords[gene][0],
                ase_vals[gene][-1],
                ase_vals[gene][1],
                ase_vals[gene][None],
                ase_vals[gene][0],
                ase_val,
                file=args.outfile, sep='\t', end='\n'
                )

