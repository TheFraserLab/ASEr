
test: testdir/pass_countsnpase testdir/pass_getgenease |testdir

testdir:
	mkdir $@

testdir/pass_countsnpase: bin/CountSNPASE.py testdir/ase.bam testdir/variants.bed testdir/reference_SNP_COUNTS.txt
	mkdir -p testdir/countsnpase_tmp
	python bin/CountSNPASE.py --mode single \
		--random-seed 0 \
		--reads testdir/ase.bam \
		--snps testdir/variants.bed \
		--prefix testdir/countsnpase_tmp/
	mv testdir/countsnpase_tmp/_SNP_COUNTS.txt testdir/SNP_COUNTS.txt
	diff -q testdir/SNP_COUNTS.txt testdir/reference_SNP_COUNTS.txt
	touch $@

testdir/pass_getgenease: testdir/pass_countsnpase testdir/ref.gtf testdir/variants.bed testdir/reference_gene_ase.tsv
	python bin/GetGeneASE.py \
		--snpcounts testdir/SNP_COUNTS.txt \
		--phasedsnps testdir/variants.bed \
		--gff testdir/ref.gtf \
		-o testdir/gene_ase.tsv \
		--writephasedsnps
	diff -q testdir/gene_ase.tsv testdir/reference_gene_ase.tsv
	touch $@

testdir/ref.gtf : | testdir
	wget -O $@ http://web.stanford.edu/~pcombs/asetest/ref.gtf

	
testdir/ase.bam : | testdir
	wget -O $@ http://web.stanford.edu/~pcombs/asetest/ase.bam

testdir/variants.bed : | testdir
	wget -O $@ http://web.stanford.edu/~pcombs/asetest/variants.bed

testdir/reference_SNP_COUNTS.txt : | testdir
	wget -O $@ http://web.stanford.edu/~pcombs/asetest/reference_SNP_COUNTS.txt

testdir/reference_gene_ase.tsv :  | testdir
	wget -O $@ http://web.stanford.edu/~pcombs/asetest/reference_gene_ase.tsv
