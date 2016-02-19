# GenomeResearch_2016_scripts
Bash and Perl scripts associated with research paper "Sister chromatid, but not NHEJ-mediated inter-chromosomal telomere fusions can occur independently of DNA ligases 3 and 4" (Genome Research, 2016).

Scripts rely on the following software.

* BWA (https://sourceforge.net/projects/bio-bwa/files/).
* samtools (https://sourceforge.net/projects/samtools/files/)
* circos (http://circos.ca/software/download/circos/)
* WGP-Toolkit (https://github.com/nestornotabilis/WGP-Toolkit)

To run pipeline:

`$ ./run.sh <ID> <PATH TO R1 FASTQ> <PATH TO R2 FASTQ> <PATH TO TELOMERE BWA INDEX> <PATH TO CUSTOM GENOME BWA INDEX> <PATH TO CUSTOM GENOME FASTA>`
