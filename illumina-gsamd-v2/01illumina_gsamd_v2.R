#=== illumina_gsamd_v2.R by PF Sullivan & Jie Song 12/2020
# GSA array processing for Illumina GSAMD v2.0
# All the new Swedish schizophrenia/control data comes from illumina-gsamd-v1/2/3 (hg19) from the genotyping lab
# The M(ulti)D(isease) content makes these into custom products, so the manifests cannot be found on Illumina's public web

# input file: GSAMD-24v2-0_20024620_A1.csv
# output file: gsamd-v2-snp-informatics-2020-12.tsv.gz
# longleaf: /nas/depts/007/sullilab/shared/bioinf/arrays/illumina-gsamd-v2

#=== output file column definitions 
#IlmnID		illumina variant ID
#Name		illumina variant name
#SNP			illumina SNP alleles
#hg38chr		hg38 chromosome (bowtie probe sequence mapping)
#hg38bp		hg38 bp
#hg38ProbeStart	hg38 start position of illumina 50bp probe (0-based)
#hg38ProbeEnd	hg38 end
#hg38str			hg38 strand
#hg38PosDup		flag for positional duplicate
#ilmnChrHg19		illumina hg19 chromosome
#ilmnBpHg19		illumina hg19 bp
#ilmnStrHg19		illumina hg19 strand
#hg19chr		hg19 chromosome (bowtie probe sequence mapping)
#hg19bp		hg19 bp
#hg19str		hg19 strand
#kg1type		1000 Genomes variant type
#kg1REF		1000 Genomes REF base (+ strand) 
#kg1ALT		1000 Genomes ALT base (+ strand)
#hrctype		HRC variant type
#hrcREF		HRC REF base (+ strand) 
#hrcALT		HRC ALT base (+ strand)
#hrcSNPid	HRC rs ID (if available)
#hrcAF		HRC allele frequency
#UnderN		Number of SNVs under hg38 probe position (TOPMed freeze 8, AF > 0.0005)
#UnderBases	Number of bases affected
#UnderMaxAF	Maximum allele frequency of SNVs under probe

#=== initialize
library(tidyverse)
library(data.table)
library(skimr)
setwd("/nas/depts/007/sullilab/shared/bioinf/arrays/illumina-gsamd-v2")

#=== unzip array manifests 
# GSAMD-24v2-0_20024620_A1.zip
system2("unzip", args = c("GSAMD-24v2-0_20024620_A1.zip"))

#=== rename files (manual step)
system("mv GSAMD-24v2-0_20024620_A1.csv GSAMD-24v2-0_20024620_A1.hg19.csv")

#=== read 
# get read start and end lines
c <- "GSAMD-24v2-0_20024620_A1.hg19.csv"
d <- paste("head -100", c, "| grep -n 'IlmnStrand' | cut -f 1 -d ':'")
e <- as.integer(system(d, intern = TRUE))-1
f <- paste("grep -n 'Controls'", c, "| cut -f 1 -d ':'")
g <- as.integer(system(f, intern = TRUE))
h <- g-e-2
i19 <- fread(file = c, header = TRUE, skip = e, nrows = h) %>%
  mutate(ilmnChrHg19 = paste0("chr", Chr)) %>%
  select(IlmnID, Name, SNP, AlleleA_ProbeSeq, AlleleB_ProbeSeq,
         GenomeBuild, ilmnChrHg19, ilmnBpHg19=MapInfo, ilmnStrHg19=RefStrand)
# mutate(i19, buildCheck = GenomeBuild==37) %>% skim(buildCheck)
rm(list=ls(pattern="^[a-h]"))


#=== use bowtie (v1) to map probes onto hg38 and hg19
# run on longleaf
# map only probeA (probeB is for A/T or C/G SNPs due to 2 dye system)
j <- "gsamd-v2.fasta"
i19 %>% 
  mutate(ll = paste0(">", Name, "\n", AlleleA_ProbeSeq)) %>%
  select(ll) %>%
  fwrite(file = j, col.names = FALSE, quote = FALSE)
system2("gzip", args = c("-f", j))

#=== upload the fasta file to longleaf, run these commands:
# TO=/nas/depts/007/sullilab/shared/bioinf/arrays/illumina-gsamd-v2
# cd ${TO}
# srun -t 8:00:00 -p interact -N 1 -n 4 --mem=16g --pty /bin/bash
# module load bowtie/1.2.3
# gunzip gsamd-v2.fasta.gz
# bowtie -m 10000 -v 0 -a -f --suppress 5,6,7,8 /nas/depts/007/sullilab/shared/bin/bowtie/grch38_1kgmaj gsamd-v2.fasta gsamd-v2.hg38.txt
# bowtie -m 10000 -v 0 -a -f --suppress 5,6,7,8 /nas/depts/007/sullilab/shared/bin/bowtie/hg19_1kgmaj gsamd-v2.fasta gsamd-v2.hg19.txt
# gzip gsamd*
# sh /nas/depts/007/sullilab/shared/bin/PERMISSON.sh ${TO}

#=== bowtie logic
# report all valid alignments (-a) up to 10000 (-m 10000)
# must be a perfect match (-v 0)
# specify fasta input (-f)
# omit redundant columns (--suppress 5,6,7,8)

#=== bowtie output hg38
# reads processed: 759993
# reads with at least one reported alignment: 727488 (95.72%)
# reads that failed to align: 32495 (4.28%)
# Reported 3039531 alignments
#=== bowtie output hg19
# reads processed: 759993
# reads with at least one reported alignment: 727623 (95.74%)
# reads that failed to align: 32360 (4.26%)
# Reported 3071869 alignments


#=== continue processing (download gsamd-v2.hg38.txt.gz and gsamd-v2.hg19.txt.gz)
# rules: + strand bp=bowtieBP+51, - strand bp=bowtieBP. probe: bowtieBP to bowtieBP+50 (0-based)
# keep those present once
k <- fread("gsamd-v2.hg38.txt.gz") %>%
  rename(Name=V1, hg38str=V2, hg38chr=V3, hg38bowtie=V4) %>%
  filter(str_length(hg38chr) <= 5) %>%
  mutate(hg38bp = as.integer(ifelse(hg38str=="+", hg38bowtie+51, hg38bowtie)),
         hg38ProbeStart = as.integer(hg38bowtie),
         hg38ProbeEnd = as.integer(hg38ProbeStart+50)) %>%
  group_by(Name) %>%
  filter(n()==1) %>%
  ungroup() %>%
  select(Name, hg38chr, hg38bp, hg38ProbeStart, hg38ProbeEnd, hg38str) #717596

l <- fread("gsamd-v2.hg19.txt.gz") %>%
  mutate(V3 = paste0("chr", V3)) %>%    # chr as digit
  rename(Name=V1, hg19str=V2, hg19chr=V3, hg19bowtie=V4) %>%
  filter(str_length(hg19chr) <= 5) %>%
  mutate(hg19bp = as.integer(ifelse(hg19str=="+", hg19bowtie+51, hg19bowtie)),
         hg19ProbeStart = as.integer(hg19bowtie),
         hg19ProbeEnd = as.integer(hg19ProbeStart+50)) %>%
  group_by(Name) %>%
  filter(n()==1) %>%
  ungroup() %>%
  select(Name, hg19chr, hg19bp, hg19str) #717128

m <- i19 %>% select(-AlleleA_ProbeSeq, -AlleleB_ProbeSeq, -GenomeBuild) %>%
  left_join(k, by = "Name") %>%
  left_join(l, by = "Name") 

# recodes after exploration: PAR1 and PAR2 mishandled
# illumina uses chrMT, bowtie chrM
# PAR hg38 chrX ≤ 2781479 or ≥ 155701382
# PAR hg19 chrX ≤ 2699520 or ≥ 154931044 
n <- m %>%
  mutate(hg38chr = ifelse(hg38chr=="chrX" & (hg38bp<=2781479 | hg38bp>=155701382), "chrXY", hg38chr),
         hg19chr = ifelse(hg19chr=="chrX" & (hg19bp<=2699520 | hg19bp>=154931044), "chrXY", hg19chr)) %>%
  group_by(hg38chr, hg38bp) %>%
  mutate(hg38PosDup = n()>1) %>%
  ungroup() %>%
  mutate(hg38PosDup = ifelse(is.na(hg38chr), as.logical(NA), hg38PosDup))


#=== get 1kg REF ALT
o <- fread("/nas/depts/007/sullilab/shared/bioinf/resources/1000GenomesV3.sites.type.freq.tsv.gz") %>% 
  select(hg38chr=chr, hg38bp, kg1type=type, kg1REF=REF, kg1ALT=ALT)
p <- n %>%
  left_join(o, by = c("hg38chr", "hg38bp")) %>%
  distinct(IlmnID, .keep_all = TRUE)     # kludge, adds 3 rows, prob dups in o
rm(list=ls(pattern="^[i-o]"))


#=== get HRC REF ALT
q <- fread("/nas/depts/007/sullilab/shared/bioinf/resources/HRC.sites.type.tsv.gz") %>% 
  select(hg38chr=chr, hg38bp, hrctype=type, hrcREF=REF, hrcALT=ALT, hrcSNPid=ID, hrcAF=AF)
r <- p %>%
  left_join(q, by = c("hg38chr", "hg38bp")) 
rm(q)


#=== variants under probe, TOPMed freeze 8, AF > 0.0005, chr1-X 
# NB! some probe positions are identical
system("mkdir tmp/")
for(ii in c(paste0("chr", seq(1:22)), "chrX")) {
  print(ii)
   r %>%
     filter(ii==hg38chr) %>%
     select(hg38chr, hg38ProbeStart, hg38ProbeEnd) %>%
     arrange(hg38chr, hg38ProbeStart, hg38ProbeEnd) %>%
     fwrite(paste0("tmp/gsa.", ii, ".bed"), col.names = FALSE, sep="\t")
  system(paste("bedtools intersect -wao",
               "-a", paste0("tmp/gsa.", ii, ".bed"),
               "-b", paste0("/nas/depts/007/sullilab/shared/topmed/08maf/snv.loc.", ii, ".hg38.bed.gz"),
               ">", paste0("tmp/out", ii, ".txt")))
}

# read files (^ means begins with)
s <- list.files(path="tmp", pattern="^out", full.names=TRUE) %>% map_df(function(x) fread(x))

# head(s)
#      V1     V2     V3 V4 V5 V6 V7 V8
#1: chr1  58763  58813  . -1 -1  .  0
#2: chr1 629241 629291  . -1 -1  .  0
#3: chr1 629861 629911  . -1 -1  .  0
#4: chr1 630053 630103  . -1 -1  .  0
#5: chr1 630077 630127  . -1 -1  .  0
#6: chr1 631712 631762  . -1 -1  .  0
# nrows are close but not exact (gsa has chrXY chrY chrM, not in topmed)
# there are also positional dups
t <- s %>%
  filter(V4 != ".") %>% 
  rename(hg38chr=V1, hg38ProbeStart=V2, hg38ProbeEnd=V3) %>%
  group_by(hg38chr, hg38ProbeStart, hg38ProbeEnd) %>%
  summarise(UnderN = n(),
            UnderBases = sum(V8),
            UnderMaxAF = max(10^-as.numeric(V7))) %>%
  ungroup()

#=== combine, save
u <- left_join(r, t, by = c("hg38chr", "hg38ProbeStart", "hg38ProbeEnd")) %>%
  select(IlmnID, Name, SNP, starts_with("hg38"), everything())
t <- "gsamd-v2-snp-informatics-2020-12.tsv"
fwrite(u, file=t, sep="\t")
system2("gzip", args = c("-f", t))

# cleanup
rm(list=ls(pattern="^[a-z]"))
system("rm GSAMD-24v2-0_20024620_A1.hg19.csv gsamd-v2.fasta.gz")
system("rm -r tmp/")

#to verify the integrity
system("md5sum gsamd-v2-snp-informatics-2020-12.tsv.gz")


