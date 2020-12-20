#=== illumina_gsa_v3.R by PF Sullivan & Jie Song 12/2020
# GSA array processing for Illumina GSA v3.0
# output files with the following info:
#   Illumina name (ILMN)
#   rsid : getting from illumina 
#   illmuchr illmubp (both for 38/37): chr.illumina.hg19/hg38, bp.illumina.hg19/hg38
#   bowtie chr /bp unique map, chr.bowtie.hg19/hg38, bp.bowtie.hg19/38
#   flag1: unique mapping in both hg19/38 (not equal to 1) could be 0 or more than 1
#   flag2: all.pos.agree (TRUE/FALSE), if bp.illumina.hg19/hg38= chr.bowtie.hg19/hg38, code to present
#   flag3: tbd
#   TOPMed info (filter for biallelic SNPs), REF and ALT at that position (is + strand), allele freqs in TOPMed pops
# longleaf: /nas/depts/007/sullilab/shared/bioinf/arrays/illumina-gsa-v3


#=== initialize
library(tidyverse)
library(data.table)
library(skimr)
# setwd("/nas/depts/007/sullilab/shared/bioinf/arrays/IlluminaGSA3.0")
setwd("/Volumes/dataMidden/genomicReferenceData/snpArrays/illumina-gsa-v3")


#=== download array manifests 
# grch38: GSA-24v3-0_A2.csv
#a <- "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A2-manifest-file-csv.zip"
#b <- word(string = a, start = -1, sep = "/")
#download.file(url = a, destfile = b)
#system2("unzip", args = c(b))

# grch37: GSA-24v3-0_A1.csv
#a <- "https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/v3-0/GSA-24v3-0-A1-manifest-file-csv.zip"
#b <- word(string = a, start = -1, sep = "/")
#download.file(url = a, destfile = b)
#system2("unzip", args = c(b))

#=== rename files (manual step)
# system("mv GSA-24v3-0_A2.csv GSA-24v3-0_A2.hg38.csv")
# system("mv GSA-24v3-0_A1.csv GSA-24v3-0_A1.hg19.csv")


#=== read hg38 version
# get read start and end lines
c <- "GSA-24v3-0_A2.hg38.csv"
d <- paste("head -100", c, "| grep -n 'IlmnStrand' | cut -f 1 -d ':'")
e <- as.integer(system(d, intern = TRUE))-1
f <- paste("grep -n 'Controls'", c, "| cut -f 1 -d ':'")
g <- as.integer(system(f, intern = TRUE))
h <- g-e-2
i38 <- fread(file = c, header = TRUE, skip = e, nrows = h) %>%
  mutate(ilmnChrHg38 = paste0("chr", Chr)) %>%
  select(IlmnID, Name, SNP, AlleleA_ProbeSeq, AlleleB_ProbeSeq, 
         GenomeBuild, ilmnChrHg38, ilmnBpHg38=MapInfo, Intensity_Only, ilmnStrHg38=RefStrand)
# mutate(i38, buildCheck = GenomeBuild==38) %>% skim(buildCheck)
rm(list=ls(pattern="^[a-h]"))


#=== read hg19 version - keep only position
# get read start and end lines
c <- "GSA-24v3-0_A1.hg19.csv"
d <- paste("head -100", c, "| grep -n 'IlmnStrand' | cut -f 1 -d ':'")
e <- as.integer(system(d, intern = TRUE))-1
f <- paste("grep -n 'Controls'", c, "| cut -f 1 -d ':'")
g <- as.integer(system(f, intern = TRUE))
h <- g-e-2
i19 <- fread(file = c, header = TRUE, skip = e, nrows = h) %>%
  mutate(ilmnChrHg19 = paste0("chr", Chr)) %>%
  select(IlmnID, SNP, GenomeBuild, ilmnChrHg19, ilmnBpHg19=MapInfo, ilmnStrHg19=RefStrand)
# mutate(I19, buildCheck = GenomeBuild==37) %>% skim(buildCheck)
rm(list=ls(pattern="^[a-h]"))


#=== combine
# both IlmnID and Name are unique
i <- inner_join(i38, i19, by = "IlmnID") %>%
  rename(SNP=SNP.x) %>% 
  select(-AlleleB_ProbeSeq, -GenomeBuild.x, -SNP.y, -GenomeBuild.y)
# nrow(i)==nrow(I38) & nrow(i)==nrow(I19)
rm(i38,i19)



#=== use bowtie (v1) to map probes onto hg38 and hg19
# run on longleaf
# map only probeA (probeB is for A/T or C/G SNPs due to 2 dye system)
j <- "gsa-v3.fasta"
i %>% 
  mutate(ll = paste0(">", Name, "\n", AlleleA_ProbeSeq)) %>%
  select(ll) %>%
  fwrite(file = j, col.names = FALSE, quote = FALSE)
system2("gzip", args = c("-f", j))

#=== upload the fasta file to longleaf, run these commands:
# TO=/nas/depts/007/sullilab/shared/bioinf/arrays/illumina-gsa-v3
# cd ${TO}
# srun -t 8:00:00 -p interact -N 1 -n 4 --mem=16g --pty /bin/bash
# module load bowtie/1.2.3
# gunzip gsa-v3.fasta.gz
# bowtie -m 10000 -v 0 -a -f --suppress 5,6,7,8 /nas/depts/007/sullilab/shared/bin/bowtie/grch38_1kgmaj gsa-v3.fasta gsa-v3.hg38.txt
# bowtie -m 10000 -v 0 -a -f --suppress 5,6,7,8 /nas/depts/007/sullilab/shared/bin/bowtie/hg19_1kgmaj   gsa-v3.fasta gsa-v3.hg19.txt
# gzip gsa*
# sh /nas/depts/007/sullilab/shared/bin/PERMISSON.sh ${TO}

#=== bowtie logic
# report all valid alignments (-a) up to 10000 (-m 10000)
# must be a perfect match (-v 0)
# specify fasta input (-f)
# omit redundant columns (--suppress 5,6,7,8)

#=== bowtie output hg38
# reads processed: 654027
# reads with at least one reported alignment: 628272 (96.06%)
# reads that failed to align: 25755 (3.94%)
# Reported 632709 alignments
#=== bowtie output hg19
# reads processed: 654027
# reads with at least one reported alignment: 628689 (96.13%)
# reads that failed to align: 25338 (3.87%)
# Reported 672448 alignments


#=== continue processing (download gsa-v3.hg38.txt.gz and gsa-v3.hg19.txt.gz)
# rules: + strand bp=bowtieBP+51, - strand bp=bowtieBP. probe: bowtieBP to bowtieBP+50 (0-based)
# keep those present once
k <- fread("gsa-v3.hg38.txt.gz") %>%
  rename(Name=V1, hg38str=V2, hg38chr=V3, hg38bowtie=V4) %>%
  filter(str_length(hg38chr) <= 5) %>%
  mutate(hg38bp = as.integer(ifelse(hg38str=="+", hg38bowtie+51, hg38bowtie)),
         hg38ProbeStart = as.integer(hg38bowtie),
         hg38ProbeEnd = as.integer(hg38ProbeStart+50)) %>%
  group_by(Name) %>%
  filter(n()==1) %>%
  ungroup() %>%
  select(Name, hg38chr, hg38bp, hg38ProbeStart, hg38ProbeEnd, hg38str)

l <- fread("gsa-v3.hg19.txt.gz") %>%
  mutate(V3 = paste0("chr", V3)) %>%    # chr as digit
  rename(Name=V1, hg19str=V2, hg19chr=V3, hg19bowtie=V4) %>%
  filter(str_length(hg19chr) <= 5) %>%
  mutate(hg19bp = as.integer(ifelse(hg19str=="+", hg19bowtie+51, hg19bowtie))) %>%
  group_by(Name) %>%
  filter(n()==1) %>%
  ungroup() %>%
  select(Name, hg19chr, hg19bp, hg19str)

m <- i %>%
  left_join(k, by = "Name") %>%
  left_join(l, by = "Name") %>%
  select(-AlleleA_ProbeSeq)

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
# o <- fread("/nas/depts/007/sullilab/shared/bioinf/resource/1000GenomesV3.sites.type.freq.tsv.gz") %>% 
o <- fread("/Volumes/dataMidden/genomicReferenceData/1000Genomes/1000GenomesV3.sites.type.freq.tsv.gz") %>%
  select(hg38chr=chr, hg38bp, kg1type=type, kg1REF=REF, kg1ALT=ALT)
p <- n %>%
  left_join(o, by = c("hg38chr", "hg38bp")) %>%
  distinct(IlmnID, .keep_all = TRUE)     # kludge, adds 3 rows, prob dups in o
rm(list=ls(pattern="^[i-o]"))


#=== get HRC REF ALT
# q <- fread("/nas/depts/007/sullilab/shared/bioinf/resource/HRC.sites.type.tsv.gz") %>% 
q <- fread("/Volumes/dataMidden/genomicReferenceData/HRC/HRC.sites.type.tsv.gz") %>%
  select(hg38chr=chr, hg38bp, hrctype=type, hrcREF=REF, hrcALT=ALT, hrcSNPid=ID, hrcAF=AF)
r <- p %>%
  left_join(q, by = c("hg38chr", "hg38bp")) %>%
rm(q)


#=== variants under probe, TOPMed freeze 8, AF > 0.0005, chr1-X
# NB! some probe positions are identical
for(ii in c(paste0("chr", seq(1:22)), "chrX")) {
  print(ii)
  # r %>%
  #   filter(ii==hg38chr) %>%
  #   select(hg38chr, hg38ProbeStart, hg38ProbeEnd) %>%
  #   arrange(hg38chr, hg38ProbeStart, hg38ProbeEnd) %>%
  #   fwrite(paste0("tmp/gsa3.", ii, ".bed"), col.names = FALSE, sep="\t")
  system(paste("bedtools intersect -wao",
               "-a", paste0("tmp/gsa3.", ii, ".bed"),
               # "-b", paste0("/nas/depts/007/sullilab/shared/TOPMed/08maf/snv.loc.", ii, ".hg38.bed.gz"),
               "-b", paste0("/Volumes/dataMidden/genomicReferenceData/TOPMed/08maf/snv.loc.", ii, ".hg38.bed.gz"),
               ">", paste0("tmp/out", ii, ".txt")))
}

# read files (^ means begins with)
s <- list.files(path="tmp", pattern="^out", full.names=TRUE) %>% map_df(function(x) fread(x))

# head(s)
#   V1     V2     V3   V4     V5     V6    V7 V8
# chr1 632828 632878    .     -1     -1     .  0
# chr1 633147 633197    .     -1     -1     .  0
# chr1 792461 792511    .     -1     -1     .  0
# chr1 817290 817340    .     -1     -1     .  0
# chr1 823656 823706 chr1 823669 823670 3.023  1
# chr1 858952 859002 chr1 858961 858962 3.143  1
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
t <- "gsa-v3-snp-informatics-2020-12.tsv"
fwrite(u, file=t, sep="\t")
system2("gzip", args = c("-f", t))

# cleanup
rm(list=ls(pattern="^[a-z]"))
system("rm GSA-24v3-0_A2.hg38.csv GSA-24v3-0_A1.hg19.csv")
system("rm -r tmp/")



