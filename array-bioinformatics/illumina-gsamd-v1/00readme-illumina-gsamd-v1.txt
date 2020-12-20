readme

PF Sullivan & Jie Song
12/2020

Output file:	gsamd-v1-snp-informatics-2020-12.tsv.gz
md5:			903b1dcfb114a4c47946c93ab6a3b86c

processing Illumina GSAMD v1 manifest files
map probe sequences using bowtie (required exactly one perfect match)
add allele frequencies
add + strand annotations (REF/ALT from 1000Genomes and HRC)
note SNV under probe (AF > 0.0005 from TOPMed freeze 8)

Illumina sources:
From the genotyping lab, required by Robert Karlsson <robert.karlsson@ki.se>, used for the new Swedish schizophrenia/control data genotyping.
The M(ulti)D(isease) content (illumina-gsamd-v1/2/3, hg19 versions only) makes these into custom products, so the manifests cannot be found on Illumina's public web

Directory contents:
- this readme
- R code
- array manifest files from genotyping lab
- bowtie (v1) mapping results
- final array infomation file


#=== column definitions gsamd-v1-snp-informatics-2020-12.tsv.gz
IlmnID		illumina variant ID
Name		illumina variant name
SNP			illumina SNP alleles
hg38chr		hg38 chromosome (bowtie probe sequence mapping)
hg38bp		hg38 bp
hg38ProbeStart	hg38 start position of illumina 50bp probe (0-based)
hg38ProbeEnd	hg38 end
hg38str			hg38 strand
hg38PosDup		flag for positional duplicate
ilmnChrHg19		illumina hg19 chromosome
ilmnBpHg19		illumina hg19 bp
ilmnStrHg19		illumina hg19 strand
hg19chr		hg19 chromosome (bowtie probe sequence mapping)
hg19bp		hg19 bp
hg19str		hg19 strand
kg1type		1000 Genomes variant type
kg1REF		1000 Genomes REF base (+ strand) 
kg1ALT		1000 Genomes ALT base (+ strand)
hrctype		HRC variant type
hrcREF		HRC REF base (+ strand) 
hrcALT		HRC ALT base (+ strand)
hrcSNPid	HRC rs ID (if available)
hrcAF		HRC allele frequency
UnderN		Number of SNVs under hg38 probe position (TOPMed freeze 8, AF > 0.0005)
UnderBases	Number of bases affected
UnderMaxAF	Maximum allele frequency of SNVs under probe



#=== first 10 lines of output file
IlmnID	Name	SNP	hg38chr	hg38bp	hg38ProbeStart	hg38ProbeEnd	hg38str	hg38PosDup	ilmnChrHg19	ilmnBpHg19	ilmnStrHg19	hg19chr	hg19bp	hg19str	kg1type	kg1REF	kg1ALT	hrctype	hrcREF	hrcALT	hrcSNPid	hrcAF	UnderN	UnderBases	UnderMaxAF
1:100292476-0_T_F_2346041289	1:100292476	[A/G]	chr1	99826920	99826869	99826919	+	FALSE	chr1	100292476	+	chr1	100292476	+	biSNP	A	G
1:101064936-0_T_F_2346041295	1:101064936	[A/G]							chr1	101064936	+
1:103380393-0_B_R_2346041316	1:103380393	[T/C]	chr1	102914837	102914837	102914887	-	FALSE	chr1	103380393	-	chr1	103380393	-	biSNP	G	A
1:104303716-0_B_R_2346041324	1:104303716	[T/C]	chr1	103761094	103761094	103761144	-	FALSE	chr1	104303716	-	chr1	104303716	-	biSNP	G	A
1:104864464-0_T_R_2346041326	1:104864464	[A/G]	chr1	104321842	104321842	104321892	-	FALSE	chr1	104864464	-	chr1	104864464	-	biSNP	C	T
1:106737318-0_T_R_2346041332	1:106737318	[A/C]							chr1	106737318	-
1:109439680-0_T_F_2348625138	1:109439680	[A/G]	chr1	108897058	108897007	108897057	+	FALSE	chr1	109439680	+	chr1	109439680	+
1:111119214-0_B_F_2346041355	1:111119214	[T/C]	chr1	110576592	110576541	110576591	+	FALSE	chr1	111119214	+	chr1	111119214	+	biSNP	C	T
1:114483147-0_T_F_2346041367	1:114483147	[A/C]	chr1	113940525	113940474	113940524	+	FALSE	chr1	114483147	+	chr1	114483147	+	biSNP	A	C