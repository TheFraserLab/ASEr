####
ASEr
####

Get ASE counts from BAMs or raw fastq data -- repackage of pipeline by Carlo Artieri

Scripts: MaskReferencefromBED.pl, create_individual_snp_files, CountSNPASE.py, create_phased_bed, GetGeneASE.py

Required software installed in PATH:
  - samtools
  - STAR (for mapping, but can use others)

Required python libraries:
  - pysam
  - pybedtools

.. contents:: **Contents**

Installation
------------

Simple install procedure::

    git clone https://github.com/TheFraserLab/ASEr.git
    cd ASEr
    python ./setup.py install

Alternatively, to install without root privileges::

    python ./setup.py install --user

Pipeline Flow
-------------

- First generate a FASTA formatted file containing the genome where each SNP position has
  been masked by 'N'. An existing genome file can be masked using the 
  MaskReferencefromBED.pl script::
       
    usage: MaskReferencefromBED.pl <SNP BED FILE> <GENOME FASTA FILE> <MASKED OUTPUT FASTA>
    
    A list of SNPs in BED format must be supplied as follows:
    
    CHR \t 0-POSITION \t 1-POSITION \t REF|ALT
    
    e.g.
    
    chr02  1242  1243  A|G

- The pipeline requires that reads mapped to the masked genome be supplied in SAM or BAM
  format. Assuming that reads will be mapped with STAR 
  (http://bioinformatics.oxfordjournals.org/content/29/1/15): The masked reference must 
  be used to create a STAR index. STAR's efficiency at mapping spliced transcripts is 
  strongly aided by supplying an annotation file in GTF format. The command to generate a
  STAR index is:
     
  ``STAR --runThreadN <NUMBER OF CORES> --runMode genomeGenerate --genomeDir <LOCATION FOR INDEX OUTPUT> --genomeFastaFiles <FIXED MASKED GENOME>.fa --sjdbGTFfile  <ANNOTATION>.gtf --sjdbOverhang <READ LENGTH - 1>``
     

- Now map the FASTQ files to the genome using STAR with the following flags. It is 
  critical that SAM/BAM files contain the MD flag for the pipeline to identify SNPs. 
  Also, because of the increased incidence of errors in the first 6 bp of reads, we trim 
  them off.
     
  ``STAR --runThreadN <NUMBER OF CORES> --genomeDir <LOCATION OF INDEX> --outFilterMultimapNmax 1 --outFileNamePrefix <PREFIX FOR OUTPUT> --outSAMtype BAM SortedByCoordinate --outSAMattributes MD NH --clip5pNbases 6``
     

- Next, we must remove duplicate reads from the mapped output. If we use the the samtools
  or the PICARD tools, we'll create a slight reference allele bias, therefore we should 
  use the rmdup.py program in the WASP package 
  (see: http://biorxiv.org/content/early/2014/11/07/011221)


- The duplicate-removed BAM file then needs to be sorted by mate-pair name rather than 
  coordinates::

  samtools sort -n [DUPLICATES REMOVED].bam [SORTED PREFIX]

***************************
create_individual_snp_files
***************************  
  
- Next we need to generate custom bed files for every individual. These files must contain
  only exonic SNPs that are heterozygous in that individual. This is necessary because
  CountSNPASE will randomly choose a SNP to count if more than one SNP is present in a read.
  As some RNA-seq reads span non-coding sequence, this can lead to SNP counts that cannot
  contribute to gene-level counts. In addition, no homozygous SNP can possibly contribute
  to ASE counts, and so randomly picking a homozygous SNP over a heterozygous SNP will lead
  to lower counts.
  
  These bed files can be generated easily by the ``create_individual_snp_files`` you should
  provide the same SNP bed file that you used in the mapping step script::
  
    usage:  create_individual_snp_files [-g] [-o outdir] [-j threads] plink snps
            create_individual_snp_files [-g] [-o outdir] [-j threads] --snpfile snp_bed --exonfile exonfile plink
            create_individual_snp_files --help

    Reference Files (positional):
      plink                The plink file, can be a prefix or full path
      snps                 A bed file of exonic SNPs. If not available, use the
                          --snpfile and --exonfile flags to generate this list on
                          the fly

    Generate exonic SNP list:
      [Use if exon_snp file doesn't already exist]

      --snpfile snp_bed    A bed file of snps (gzipped OK)
      --exonfile exon_bed  A bed file of exons (gzipped OK).

    Filter Individuals:
      -f , --filter        List of individuals to keep. Either comma-separated or
                          as a file (newline separated
      --split_name         Split the individual name in plink on this character
                          (must be a single character)
      --split_index        Index of split name to compare individuals to

    Multi(plex) mode arguments:
      -j , --jobs          Divide into # of jobs
      -w , --walltime      Walltime for each job
      -k , --mem           Memory for each job in MB (int)
      --queue              Queue to submit jobs to
      --cluster            Which cluster to use

    Optional Arguments:
      -o , --outdir        The output directory to write files to
      -g, --gzip           gzip compress the output files
      -q, --quiet          Quiet output
      -v, --verbose        Verbose output
      -h, --help           Show this help and exit.

  
***********
CountSNPASE
***********

- Now we can count reads overlapping each SNP. The CountSNPASE.py script does this::
  
    usage: CountSNPASE.py -m mode -s <BED> -r <[S/B]AM> [-p] [-b] [-n] [-q] [-v]

    Count number of reads overlapping each SNP in a sam/bam file.

    Required arguments:
      -m mode, --mode mode  Operation mode (default: None)
      -s <BED>, --snps <BED>
                            SNP BED file (default: None)
      -r <[S/B]AM>, --reads <[S/B]AM>
                            Mapped reads file [sam or bam] (default: None)

    Universal optional arguments:
      -p , --prefix         Prefix for temp files and output (default: TEST)
      -b, --bam             Mapped read file type is bam (auto-detected if *.bam)
                            (default: False)
      -n, --noclean         Do not delete intermediate files (for debuging)
                            (default: False)
      -h, --help            show this help message and exit

    Multi(plex) mode arguments:
      -j , --jobs           Divide into # of jobs (default: 100)
      -w , --walltime       Walltime for each job (default: 3:00:00)
      -k , --mem            Memory for each job (default: 5000MB)
      --queue               Queue to submit jobs to (default: batch)
      --cluster {torque,slurm,normal}
                            Which cluster to use, normal uses threads on this
                            machine (default: torque)
      --threads             Max number of threads to run at a time (normal mode
                            only). (default: 24)

    Single mode arguments:
      -f , --suffix         Suffix for multiplexing [set automatically] (default:
                            )

    Logging options:
      -q, --quiet           Quiet mode, only prints warnings. (default: False)
      -v, --verbose         Verbose mode, prints debug info too. (default: False)
      --logfile LOGFILE     Logfile to write messages too, default is STDERR
                            (default: None)

    Detailed description of inputs/outputs follows:

    -s/--snps
        A tab-delimited BED file with positions of masked SNPs of interest as follows:

        [CHR]    [0 POSITION]    [1 POSITION]

        Additional columns are ignored.

    -r/--reads
        A SAM or BAM file containing all of the reads masked to the masked genome. The file
        shound have all duplicates removed and MUST be sorted by read name
        (i.e. samtools sort -n ).

    -m/--mode
        The script can be run in two modes. In 'single' mode, the entire SNP counting is performed
        locally. In 'multi' mode, the read file will be split up by the number of specified jobs on
        the cluster. This is much faster for large SAM/BAM files.

    OUTPUT:

    The output of the script is a tab-delimited text file, [PREFIX]_SNP_COUNTS.txt, which contains the
    following columns:

    CHR		        Chromosome where SNP is found
    POSITION	    1-based position of SNP
    POS_A|C|G|T	  Count of reads containing A|C|G|T bases at the SNP position on the POSITIVE strand
    NEG_A|C|G|T	  Count of reads containing A|C|G|T bases at the SNP position on the NEGATIVE strand
    SUM_POS_READS	Sum of all reads assigned to the SNP on POSITIVE strand
    SUM_NEG_READS	Sum of all reads assigned to the SNP on NEGATIVE strand
    SUM_READS	    Sum of all reads assigned to the SNP

Note: this script has a multiplexing mode that can dramatically accelerate its performance by splitting
sam/bam files and running in parallel on all the fragments. This mode will can be enabled with the
``-m multi`` argument. On a simple system it will just use threads up to a maximum of ``--threads``. On
a system with torque or slurm, it will submit its jobs to those systems. The cluster system is
auto-detected, but you will need to provide the queue/partition to run in and other submission variables.

*****************
create_phased_bed
*****************

- To run GetGeneASE, a bed file of phased SNPs is required. This can be created by running shapeit on the
  plink data of your individuals. Note: you may wish to run impute2 on your data also ot increase your
  power to detect SNPs.

- To create the bed file from the .haps file output by shapeit, run the ``create_phase_bed`` script::

    usage:  create_phased_bed -i hap_file [hap_file...] -o bed_file
            cat hap_file | create_phased_bed > bed_file
            create_phased_bed --help

    Create a phased SNP bed file from haplotype data.

    Files:
      -i infile [infile ...], --hap_files infile [infile ...]
                            List of haps files, default STDIN
      -o , --bed_file       Output bed file, default STDOUt, gzipped OK

    Optional Arguments:
      --chr_format {num,chr}
                            Convert chromsome to number only (num) or to chr#
                            (chr)
      -h, --help            Show this help and exit.

**********  
GetGeneASE
**********

- Once we've determined the counts at individual SNPs, we can then obtain the gene/
  transcript-level counts with GetGeneASE.py::
     
    usage: GetGeneASE.py -c  -p  -g  -o  [-w] [-i] [-t] [-m MIN] [-s] [-h]

    This script takes the output of CountSNPASE.py and generates gene level ASE counts.

    Required arguments::
      -c , --snpcounts      SNP-level ASE counts from CountSNPASE.py (default:
                None)
      -p , --phasedsnps     BED file of phased SNPs (default: None)
      -g , --gff            GFF/GTF formatted annotation file (default: None)
      -o , --outfile        Gene-level ASE counts output (default: None)

    Optional arguments::
      -w, --writephasedsnps
                Write a phased SNP-level ASE output file
                [OUTFILE].snps.txt (default: False)
      -i , --identifier     ID attribute in information column (default: gene_id)
      -t , --type           Annotation feature type (default: exon)
      -m MIN, --min MIN     Min reads to calculate proportion ref/alt biased
                (default: 10)
      -s, --stranded        Data are stranded? [Default: False] (default: False)
      -h, --help            Show this help message and exit

    NOTE:  SNPs that overlap multiple features on the same strand (or counting from 
        unstranded libraries) will be counted in EVERY feature that they overlap. It is
        important to filter the annotation to count features of interest!  

    Detailed description of inputs/outputs follows:

    -p/--phasedsnps 
      A tab-delimited BED file with positions of masked SNPs of interest as follows:

      [CHR]  [0 POSITION]  [1 POSITION]  [REF|ALT]

      The fourth column MUST contain the phased SNPs alleles. 

    -g/--gff
      The script accepts both GTF and GFF annotation files. This should be combined with
      the -i/--identifier option specifying the identifier in the info column (column 9) 
      that will be used for grouping counts. For example, in a GTF 'gene_id' will group
      counts by gene with 'transcript_id' with group counts by transcript. In addition,
      the -t/--type option sets the feature type (column 3) from which to pull features
      typically you'd want to count from 'exon', but many annotations may use non-
      standard terms.

    -m/--min
      This sets the minimum # of reads required to include a SNP in the calculation of 
      the fraction of SNPs agreeing in allelic direction.

    -w/--writephasedsnps
      If this is specified, then the program will output an additional output file named
      [OUTFILE].snp.txt with phased SNP-level ASE calls. This can be useful for checking
      SNP consistency across samples. See below for a description of the output.

    -s/--stranded
      If the data come from a stranded library prep, then this option will only count 
      reads mapped to the corresponding strand.
     
    OUTPUT:

    The output of the script is a tab-delimited text file set by -o/--outfile, which 
    contains the following columns:

    FEATURE            Name of the counted feature  
    CHROMOSOME         Chromosome where feature is found
    ORIENTATION        Orientation of feature (+/-)
    START-STOP         Ultimate 5' and 3' 1-based start and stop positions
    REFERENCE_COUNTS   Total reference allele counts across SNPS (or first allele in the REF|ALT phasing)
    ALT_COUNTS         Total alternate allele counts across SNPs (or second allele in the REF|ALT phasing)
    TOTAL_SNPS         The total number of SNPs overlapped by the feature 
    REF_BIASED         Number of REF biased SNPs passing the -m/--min threshold
    ALT_BIASED         Number of ALT biased SNPs passing the -m/--min threshold
    REF-ALT_RATIO      The proportion of SNPs agreeing in direction (0.5 - 1)
    SNPS               A list of all SNPs overlapped by the feature separated by ';' and of the format:

      [1-based position],[REF_ALLELE]|[ALT_ALLELE],[REF_COUNTS]|[ALT_COUNTS];

    If the -w/--writephasedsnps option has been set, it will produce a tab-delimited table 
    with the following columns:

    CHROMOSOME         Chromosome where SNP is found
    POSITION           1-based position
    FEATURE            Feature in which SNP is found
    ORIENTATION        Orientation of feature (if stranded only reads on this strand are counted)
    REFERENCE_ALLELE   Reference base
    ALTERNATE_ALLELE   Alternate base
    REF_COUNTS         Reference base counts
    ALT_COUNTS         Alternate base counts

