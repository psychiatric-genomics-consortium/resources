readme

PF Sullivan & Jie Song
12/2020

Output file:	gsa-v3-snp-informatics-2020-12.tsv.gz
md5:			86c5b059279d5a7d6f6d4e791e50b057

processing Illumina GSA v3 manifest files
map probe sequences using bowtie (required exactly one perfect match)
add allele frequencies
add + strand annotations (REF/ALT from 1000Genomes and HRC)
note SNV under probe (AF > 0.0005 from TOPMed freeze 8)

Illumina sources:
https://emea.illumina.com/products/by-type/microarray-kits/infinium-global-screening.html
https://emea.support.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/infinium-global-screening-array-data-sheet-370-2016-016.pdf

Directory contents:
- this readme
- R code
- array manifest files from Illumina
- Illumina product sheet pdf
- bowtie (v1) mapping results
- final array infomation file


#=== column definitions gsa-v3-snp-informatics-2020-12.tsv.gz
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
Intensity_Only	illumina flag
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
IlmnID	Name	SNP	hg38chr	hg38bp	hg38ProbeStart	hg38ProbeEnd	hg38str	hg38PosDup	ilmnChrHg38	ilmnBpHg38	Intensity_Only	ilmnStrHg38	ilmnChrHg19	ilmnBpHg19	ilmnStrHg19	hg19chr	hg19bp	hg19str	kg1typekg1REF	kg1ALT	hrctype	hrcREF	hrcALT	hrcSNPid	hrcAF	UnderN	UnderBases	UnderMaxAF
1:103380393-0_B_R_2346041316	1:103380393	[T/C]	chr1	102914837	102914837	102914887	-	FALSE	chr1	102914837	0	-	chr1	103380393	-	chr1	103380393	-	biSNP	""	""	""	""	""			
1:109439680-0_T_F_2348625138	1:109439680	[A/G]	chr1	108897058	108897007	108897057	+	FALSE	chr1	108897058	0	+	chr1	109439680	+	chr1	109439680	+	""	""	""	""	""	""	""	""			
ilmnseq_1:110198788-3_B_F_2599938692	1:110198788	[T/C]	chr1	109656166	109656115	109656165	FALSE	chr1	109656166	0	+	chr1	110198788	+	chr1	110198788	+	biSNP	""	""	""	""	""			
ilmnseq_1:110201112-3_B_R_2599950325	1:110201112	[G/C]	""				""		chr1	109658490	0	-	chr1	110201112	-	""		""	""	""	""	""	""	""	""	""			
ilmnseq_1:110201667-3_B_F_2599938696	1:110201667	[T/C]	chr1	109659045	109658994	109659044	FALSE	chr1	109659045	0	+	chr1	110201667	+	chr1	110201667	+	biSNP	""	""	""	""	""	1	1	0.000887156012037961
ilmnseq_1:110202904-3_T_R_2599938699	1:110202904	[A/G]	chr1	109660282	109660282	109660332	FALSE	chr1	109660282	0	-	chr1	110202904	-	chr1	110202904	-	biSNP	""	""	""	""	""			
ilmnseq_1:110203240_ilmndup1-3_T_R_2610792938	1:110203240	[A/C]	chr1	109660618	109660618	109660668	-	FALSE	chr1	109660618	0	-	chr1	110203240	-	chr1	110203240	biSNP	G	T	""	""	""	""	""	2	2	0.0163305194789433
ilmnseq_1:110203911-3_B_F_2599938716	1:110203911	[T/C]	chr1	109661289	109661238	109661288	FALSE	chr1	109661289	0	+	chr1	110203911	+	chr1	110203911	+	biSNP	""	""	""	""	""			
ilmnseq_1:110206675-3_B_R_2599938729	1:110206675	[T/C]	chr1	109664053	109664053	109664103	FALSE	chr1	109664053	0	-	chr1	110206675	-	chr1	110206675	-	biSNP	biSNP	A	G	rs78856997	0.00331075			
