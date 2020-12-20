readme

PF Sullivan & Jie Song
12/2020

Output file:	gsa-v1-snp-informatics-2020-12.tsv.gz
md5:			b02f4302845cb8f5a47478da5b685fd4

processing Illumina GSA v1 manifest files
map probe sequences using bowtie (required exactly one perfect match)
add allele frequencies
add + strand annotations (REF/ALT from 1000Genomes and HRC)
note SNV under probe (AF > 0.0005 from TOPMed freeze 8)

Illumina sources:
https://emea.illumina.com/products/by-type/microarray-kits/infinium-global-screening.html (for gsa-v3)
https://grcf.jhmi.edu/wp-content/uploads/2017/12/infinium-commercial-gsa-data-sheet-370-2016-016.pdf

Directory contents:
- this readme
- R code
- array manifest files from Illumina
- Illumina product sheet pdf
- bowtie (v1) mapping results
- final array infomation file


#=== column definitions gsa-v1-snp-informatics-2020-12.tsv.gz
IlmnID		illumina variant ID
Name		illumina variant name
SNP			illumina SNP alleles
hg38chr		hg38 chromosome (bowtie probe sequence mapping)
hg38bp		hg38 bp
hg38ProbeStart	hg38 start position of illumina 50bp probe (0-based)
hg38ProbeEnd	hg38 end
hg38str			hg38 strand
hg38PosDup		flag for positional duplicate
ilmnChrHg38		illumina hg38 chromosome
ilmnBpHg38		illumina hg38 bp
ilmnStrHg38		illumina hg38 strand
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
IlmnID	Name	SNP	hg38chr	hg38bp	hg38ProbeStart	hg38ProbeEnd	hg38str	hg38PosDup	ilmnChrHg38	ilmnBpHg38	ilmnStrHg38	ilmnChrHg19	ilmnBpHg19	ilmnStrHg19	hg19chr	hg19bp	hg19str	kg1type	kg1REF	kg1ALT	hrctype	hrcREF	hrcALT	hrcSNPid	hrcAF	UnderN	UnderBases	UnderMaxAF
1:103380393-0_B_R_2346041316	1:103380393	[T/C]	chr1	102914837	102914837	102914887	-	FALSE	chr1	102914837	-	chr1	103380393	-	chr1	103380393	-	biSNP	G	A
1:106737318-0_T_R_2346041332	1:106737318	[A/C]							chr1	106194696	-	chr1	106737318	-
1:109439680-0_T_F_2348625138	1:109439680	[A/G]	chr1	108897058	108897007	108897057	+	FALSE	chr1	108897058	+	chr1	109439680	+	chr1	109439680	+
1:118227370-0_T_R_2346041384	1:118227370	[A/G]	chr1	117684748	117684748	117684798	-	FALSE	chr1	117684748	-	chr1	118227370	-	chr1	118227370	-	biSNP	T	C
1:1183442-0_T_F_2346041385	1:1183442	[A/G]	chr1	1248062	1248011	1248061	+	FALSE	chr1	1248062	+	chr1	1183442	+	chr1	1183442	+	biSNP	A	G						1	1	0.0555904257270404
1:118933200-0_B_F_2346041391	1:118933200	[T/G]	chr1	118390577	118390526	118390576	+	FALSE	chr1	118390577	+	chr1	118933200	+	chr1	118933200	+	multiallelic	G,GT	T,G
1:11907740-0_B_R_2348625055	1:11907740	[T/C]	chr1	11847683	11847683	11847733	-	FALSE	chr1	11847683	-	chr1	11907740	-	chr1	11907740	-
1:119872141-0_B_F_2346041397	1:119872141	[T/C]	chr1	119329518	119329467	119329517	+	FALSE	chr1	119329518	+	chr1	119872141	+	chr1	119872141	+	biSNP	T	C
1:120608075-0_B_R_2346041408	1:120608075	[T/C]	chr1	120065461	120065461	120065511	-	FALSE	chr1	120065461	-	chr1	120608075	-	chr1	120608075	-	biSNP	G	A	biSNP	G	A	rs61200250	0.0217893	1	1	0.000907820530178186