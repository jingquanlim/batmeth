Description:
DNA methylation plays a crucial role in higher organisms. Coupling bisulfite treatment with next-generation sequencing enables the interrogation of 5-methylcytosine sites in the genome. However, bisulfite conversion introduces mismatches between the reads and the reference genome, which makes mapping of Illumina and SOLiD reads slow and inaccurate. BatMeth is an algorithm that integrates novel mismatch counting, list filtering, mismatch stage filtering and fast mapping onto two indexes to improve unique mapping rate, speed and precision. Experimental results show that BatMeth is faster and more accurate than existing tools.

In short: BatMeth maps BS-reads (Illumina base-reads and SOLiD color-reads) onto a reference genome efficiently and accurately.

Current version: Google-SVN (Illumina reads), 1.04b in download section (ABi color reads)

Pre-Built Index: hg19-base-space, hg19-color-space http://compbio.ddns.comp.nus.edu.sg/~limjingq/BATMETH_INDEX/

"BatMeth: Improved mapper for bisulfite sequencing reads on DNA methylation" was published in Genome Biology's Special [Issue 2012](https://code.google.com/p/batmeth/issues/detail?id=2012) on Epigenomics.