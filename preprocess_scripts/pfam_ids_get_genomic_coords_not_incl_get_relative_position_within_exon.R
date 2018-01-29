#Load libraries.

library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)

#Load Rdata.

installation_dir_FusionVis <- getwd()
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_exons_as_data_frame.Rdata"))

#Get a version of exons_as_data_frame that includes only exons with non-NA CDS.

exons_with_CDS_dat <- exons_as_data_frame[which(is.na(exons_as_data_frame$CDSLENGTH) == FALSE),]

rm(exons_as_data_frame)

#For each transcript, get the amino acids of the transcript contained within the exon.
#For each transcript, we need to get the coordinates within the whole CDS for each exon.
#For example, if we have three exons, CDS length 67, 81, and 110, we want the cumulative CDS intervals to be 1-67, 68-148, and 149-258.

#To do this, we'll want to add to each exon a cumulative sum of the CDS lengths of the exons before the current exon.
#So we make a dummy exon for each transcript (we'll call it EXONRANK=0), where the CDSLENGTH will be 0.

transcripts_in_exons_with_CDS_dat <- unique(exons_with_CDS_dat$TXNAME)
num_transcripts_in_exons_with_CDS_dat <- length(transcripts_in_exons_with_CDS_dat)

dummy_exons <- data.frame(TXNAME = unique(exons_with_CDS_dat$TXNAME),
CDSSTART = rep(NA,times=num_transcripts_in_exons_with_CDS_dat),
CDSEND = rep(NA,times=num_transcripts_in_exons_with_CDS_dat),
EXONCHROM = rep(NA,times=num_transcripts_in_exons_with_CDS_dat),
EXONSTRAND = rep("+",times=num_transcripts_in_exons_with_CDS_dat),
EXONSTART = rep(NA,times=num_transcripts_in_exons_with_CDS_dat),
EXONEND = rep(NA,times=num_transcripts_in_exons_with_CDS_dat),
EXONRANK = rep(0,times=num_transcripts_in_exons_with_CDS_dat),
CDSLENGTH = rep(0,times=num_transcripts_in_exons_with_CDS_dat),
EXONLENGTH = rep(0,times=num_transcripts_in_exons_with_CDS_dat),
stringsAsFactors=FALSE)

#Add the dummy exon data frame to exons_with_CDS_dat, then order exons_with_CDS_dat by transcript then exon rank so the dummy exon always appears first for each transcript.

exons_with_CDS_dat <- rbind(exons_with_CDS_dat,dummy_exons)

exons_with_CDS_dat <- exons_with_CDS_dat[order(exons_with_CDS_dat$TXNAME,exons_with_CDS_dat$EXONRANK),]

rownames(exons_with_CDS_dat) <- paste0(exons_with_CDS_dat$TXNAME,"/",exons_with_CDS_dat$EXONRANK)

#Calculate the cumulative sum of the CDS lengths as you move along each transcript.

cumsum_CDS_length_per_transcript <- aggregate(exons_with_CDS_dat$CDSLENGTH ~ exons_with_CDS_dat$TXNAME,FUN=cumsum)
colnames(cumsum_CDS_length_per_transcript) <- c("TXNAME","CUMCDS")

#First few rows look like this:

#           TXNAME                                                 CUMCDS
#1 ENST00000000233                         0, 67, 148, 258, 330, 456, 543
#2 ENST00000000412                        0, 176, 343, 453, 584, 711, 834
#3 ENST00000000442                      0, 325, 442, 571, 742, 1012, 1272
#4 ENST00000001008 0, 105, 250, 393, 514, 671, 762, 846, 1032, 1272, 1380
#5 ENST00000001146                      0, 204, 429, 705, 861, 1146, 1539
#6 ENST00000002125   0, 55, 216, 297, 408, 622, 681, 792, 936, 1110, 1326

#For each vector in CUMCDS, we want all but the last item.

#Convert so it looks like this.

#TXNAME	CUMCDS
#ENST00000000233	0
#ENST00000000233	67
#..
#ENST00000000233	456
#ENST00000000233	543
#ENST00000000412	0
#ENST00000000412	176
#..
#ENST00000000412	711
#ENST00000000412	834

#Then, remove the row with the highest CUMCDS for each transcript from this table to remove the last item.

cumsum_CDS_length_per_transcript_vector <- as.vector(unlist(cumsum_CDS_length_per_transcript[,2]))
exons_per_transcript_plus_one <- as.vector(sapply(cumsum_CDS_length_per_transcript[,2],FUN=function(x)length(x)))

cumsum_CDS_length_per_transcript <- data.frame(TXNAME = rep(cumsum_CDS_length_per_transcript$TXNAME,times=exons_per_transcript_plus_one),
CUMCDS = cumsum_CDS_length_per_transcript_vector,
stringsAsFactors=FALSE)

#To remove row with highest CUMCDS, reverse sort by CUMCDS, take only rows that are the second or later occurence of TXNAME.
#Then, convert back to regular (not reverse) sort.

cumsum_CDS_length_per_transcript <- cumsum_CDS_length_per_transcript[order(cumsum_CDS_length_per_transcript$TXNAME,cumsum_CDS_length_per_transcript$CUMCDS,decreasing=TRUE),]
cumsum_CDS_length_per_transcript <- cumsum_CDS_length_per_transcript[which(duplicated(cumsum_CDS_length_per_transcript$TXNAME) == TRUE),]
cumsum_CDS_length_per_transcript <- cumsum_CDS_length_per_transcript[order(cumsum_CDS_length_per_transcript$TXNAME,cumsum_CDS_length_per_transcript$CUMCDS),]

#Remove the dummy exons (EXONRANK = 0) from exons_with_CDS_dat, then merge with the cumulative CDS length info.

exons_with_CDS_dat <- exons_with_CDS_dat[which(exons_with_CDS_dat$EXONRANK > 0),]

exons_with_CDS_dat <- data.frame(exons_with_CDS_dat,
CUMCDS = cumsum_CDS_length_per_transcript$CUMCDS,
stringsAsFactors=FALSE)

#Add columns CDS.base.within.transcript.start and CDS.base.within.transcript.end.
#CDS.base.within.transcript.start will be one base beyond the CDS bases we have gone through so far as we move along the transcript.
#CDS.base.within.transcript.end will add the length of the CDS of the current exon to the cumulative CDS length from all previous exons.

exons_with_CDS_dat <- data.frame(exons_with_CDS_dat,
CDS.base.within.transcript.start = exons_with_CDS_dat$CUMCDS  + 1,
CDS.base.within.transcript.end = exons_with_CDS_dat$CUMCDS + exons_with_CDS_dat$CDSLENGTH,
stringsAsFactors=FALSE)

#Now, let's get amino acid coordinates per protein domain per transcript from Ensembl Biomart.
#Then, to convert relative amino acids to relative bases, multiply by three.
#Then for the start, subtract 2, so we can start with the first base of the amino acid instead of the third.

library(biomaRt)

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",dataset="hsapiens_gene_ensembl")
pfam_start_and_end_per_transcript <- getBM(attributes = c("ensembl_transcript_id","pfam","pfam_start","pfam_end"),mart=ensembl)
pfam_start_and_end_per_transcript <- data.frame(pfam_start_and_end_per_transcript,
CDS.base.within.transcript.start = (pfam_start_and_end_per_transcript$pfam_start*3) - 2,
CDS.base.within.transcript.end = pfam_start_and_end_per_transcript$pfam_end*3,
stringsAsFactors=FALSE)

#Order just for ease of viewing.

pfam_start_and_end_per_transcript <- pfam_start_and_end_per_transcript[order(pfam_start_and_end_per_transcript$ensembl_transcript_id,pfam_start_and_end_per_transcript$pfam_start),]

#Remove transcripts that do not have any Pfam ID.

pfam_start_and_end_per_transcript <- pfam_start_and_end_per_transcript[which(is.na(pfam_start_and_end_per_transcript$pfam_start) == FALSE),]

#Now, let's check the max CDS base per transcript in both data sets.
#Do the domains ever go beyond what the CDS length is supposed to be?

#pfam_max_CDS_base_within_transcript_end = aggregate(CDS.base.within.transcript.end ~ ensembl_transcript_id,data=pfam_start_and_end_per_transcript,FUN=max)

#overall_max_CDS_base_within_transcript_end = aggregate(CDS.base.within.transcript.end ~ TXNAME,data=exons_with_CDS_dat,FUN=max)

#overall_max_CDS_base_within_transcript_end <- overall_max_CDS_base_within_transcript_end[match(pfam_max_CDS_base_within_transcript_end[,1],overall_max_CDS_base_within_transcript_end[,1]),]

#pfam_vs_overall_max_CDS_base_within_transcript_end <- data.frame(Pfam = pfam_max_CDS_base_within_transcript_end[,2],Overall = overall_max_CDS_base_within_transcript_end[,2],
#row.names=pfam_max_CDS_base_within_transcript_end[,1],
#stringsAsFactors=FALSE)

#range(pfam_vs_overall_max_CDS_base_within_transcript_end[,2] - pfam_vs_overall_max_CDS_base_within_transcript_end[,1])

#Looks like in some cases, there is a Pfam domain that goes one beyond the amino acids that are supposed to be in the transcript.
#I checked and these are cases like this:

#                           TXNAME EXONRANK CDSLENGTH EXONLENGTH CUMCDS CDS.base.within.transcript.start CDS.base.within.transcript.end
#ENST00000341400/1 ENST00000341400        1       120        399      0                                1                            120
#ENST00000341400/2 ENST00000341400        2       129        129    120                              121                            249
#ENST00000341400/3 ENST00000341400        3       100        100    249                              250                            349
#ENST00000341400/4 ENST00000341400        4       139        139    349                              350                            488

#       ensembl_transcript_id    pfam pfam_start pfam_end CDS.base.within.transcript.start CDS.base.within.transcript.end
#165646       ENST00000341400 PF06409         41      163                              121                            489

#Shouldn't be a big deal. In these cases, just trim the domain a bit so it does not go beyond the CDS when we intersect with exons in a bit.

#Turn both exons_with_CDS_dat and pfam_start_and_end_per_transcript into special GRanges (with each transcript ID as a sequence) so we can intersect.

exons_with_CDS_GRanges_of_CDS_coord_within_transcript <- GRanges(seqnames = exons_with_CDS_dat$TXNAME,
ranges = IRanges(exons_with_CDS_dat$CDS.base.within.transcript.start,exons_with_CDS_dat$CDS.base.within.transcript.end),
name = rownames(exons_with_CDS_dat),
mcols = exons_with_CDS_dat[,c("TXNAME","CDSSTART","CDSEND","EXONCHROM","EXONSTRAND","EXONSTART","EXONEND","EXONRANK","CDSLENGTH","EXONLENGTH")])

pfam_domains_GRanges_of_CDS_coord_within_transcript <- GRanges(seqnames = pfam_start_and_end_per_transcript$ensembl_transcript_id,
ranges = IRanges(pfam_start_and_end_per_transcript$CDS.base.within.transcript.start,pfam_start_and_end_per_transcript$CDS.base.within.transcript.end),
name = paste0(pfam_start_and_end_per_transcript$ensembl_transcript_id,"/",pfam_start_and_end_per_transcript$pfam,"/",pfam_start_and_end_per_transcript$pfam_start),
mcols = pfam_start_and_end_per_transcript[,c("pfam","pfam_start","pfam_end")])

overlaps_exons_and_pfam_domains_from_CDS_coord_within_transcript <- mergeByOverlaps(exons_with_CDS_GRanges_of_CDS_coord_within_transcript,
pfam_domains_GRanges_of_CDS_coord_within_transcript,maxgap=-1L)

exons_overlapping_pfam_domains_from_CDS_coord_within_transcript <- BiocGenerics::as.data.frame(overlaps_exons_and_pfam_domains_from_CDS_coord_within_transcript[,"exons_with_CDS_GRanges_of_CDS_coord_within_transcript"])

pfam_domains_overlapping_exons_from_CDS_coord_within_transcript <- BiocGenerics::as.data.frame(overlaps_exons_and_pfam_domains_from_CDS_coord_within_transcript[,"pfam_domains_GRanges_of_CDS_coord_within_transcript"])

#head(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript)
#head(pfam_domains_overlapping_exons_from_CDS_coord_within_transcript)

#nrow(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript)
#nrow(pfam_domains_overlapping_exons_from_CDS_coord_within_transcript)

#Save progress.

save.image("pfam_ids_get_genomic_coords.Rdata")

#Now, if we want to convert these to genomic coordinates, we can easily do so.

#         seqnames start end width strand              name    mcols.TXNAME
#1 ENST00000000233     1  67    67      * ENST00000000233/1 ENST00000000233
#2 ENST00000000233    68 148    81      * ENST00000000233/2 ENST00000000233
#3 ENST00000000233   149 258   110      * ENST00000000233/3 ENST00000000233
#4 ENST00000000233   259 330    72      * ENST00000000233/4 ENST00000000233
#5 ENST00000000233   331 456   126      * ENST00000000233/5 ENST00000000233
#6 ENST00000000233   457 543    87      * ENST00000000233/6 ENST00000000233
#  mcols.CDSSTART mcols.CDSEND mcols.EXONCHROM mcols.EXONSTRAND mcols.EXONSTART
#1      127228553    127228619               7                +       127228399
#2      127229137    127229217               7                +       127229137
#3      127229539    127229648               7                +       127229539
#4      127230120    127230191               7                +       127230120
#5      127231017    127231142               7                +       127231017
#6      127231267    127231353               7                +       127231267
#  mcols.EXONEND mcols.EXONRANK mcols.CDSLENGTH mcols.EXONLENGTH
#1     127228619              1              67              221
#2     127229217              2              81               81
#3     127229648              3             110              110
#4     127230191              4              72               72
#5     127231142              5             126              126
#6     127231759              6              87              493
#         seqnames start end width strand                      name mcols.pfam
#1 ENST00000000233    19 528   510      * ENST00000000233/PF00025/7    PF00025
#2 ENST00000000233    19 528   510      * ENST00000000233/PF00025/7    PF00025
#3 ENST00000000233    19 528   510      * ENST00000000233/PF00025/7    PF00025
#4 ENST00000000233    19 528   510      * ENST00000000233/PF00025/7    PF00025
#5 ENST00000000233    19 528   510      * ENST00000000233/PF00025/7    PF00025
#6 ENST00000000233    19 528   510      * ENST00000000233/PF00025/7    PF00025
#  mcols.pfam_start mcols.pfam_end
#1                7            176
#2                7            176
#3                7            176
#4                7            176
#5                7            176
#6                7            176
#[1] 466090
#[1] 466090

overlap_exons_with_pfam_to_get_pfam_genomic_coord <- data.frame(Exon.CDS.base.within.transcript.start = exons_overlapping_pfam_domains_from_CDS_coord_within_transcript$start,
Exon.CDS.base.within.transcript.end = exons_overlapping_pfam_domains_from_CDS_coord_within_transcript$end,
Pfam.CDS.base.within.transcript.start = pfam_domains_overlapping_exons_from_CDS_coord_within_transcript$start,
Pfam.CDS.base.within.transcript.end = pfam_domains_overlapping_exons_from_CDS_coord_within_transcript$end,
EXONSTRAND = exons_overlapping_pfam_domains_from_CDS_coord_within_transcript$mcols.EXONSTRAND,
CDSSTART = exons_overlapping_pfam_domains_from_CDS_coord_within_transcript$mcols.CDSSTART,
CDSEND = exons_overlapping_pfam_domains_from_CDS_coord_within_transcript$mcols.CDSEND,
Pfam.genomic.start = rep(NA,times=nrow(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript)),
Pfam.genomic.end = rep(NA,times=nrow(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript)),
stringsAsFactors=FALSE)

#Where Pfam.CDS.base.within.transcript.end is larger than Exon.CDS.base.within.transcript.end, set it equal to Exon.CDS.base.within.transcript.end.
#Where Pfam.CDS.base.within.transcript.start is less than Exon.CDS.base.within.transcript.start, set it equal to Exon.CDS.base.within.transcript.start.

overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.base.within.transcript.end[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.base.within.transcript.end > overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.base.within.transcript.end)] <- overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.base.within.transcript.end[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.base.within.transcript.end > overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.base.within.transcript.end)]

overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.base.within.transcript.start[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.base.within.transcript.start < overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.base.within.transcript.start)] <- overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.base.within.transcript.start[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.base.within.transcript.start < overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.base.within.transcript.start)]

#Shorten some column names so easier to look at.

colnames(overlap_exons_with_pfam_to_get_pfam_genomic_coord)[1:4] <- c("Exon.CDS.start","Exon.CDS.end","Pfam.CDS.start","Pfam.CDS.end")

#Now, only problem is that Pfam.CDS.start and Pfam.CDS.end are in overall CDS coordinates, not CDS coordinates within this exon.
#Need to change this by switching to be in comparison to Exon.CDS.start.

overlap_exons_with_pfam_to_get_pfam_genomic_coord <- data.frame(overlap_exons_with_pfam_to_get_pfam_genomic_coord,
Pfam.this.exon.CDS.start = overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.start - overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.start + 1,
Pfam.this.exon.CDS.end = overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.CDS.end - overlap_exons_with_pfam_to_get_pfam_genomic_coord$Exon.CDS.start + 1,
stringsAsFactors=FALSE)

#Now, if on positive strand, get genomic start and end using CDSSTART + Pfam start/end - 1.
#If on negative strand, instead do CDSEND - Pfam start/end + 1.
#Also if on negative strand, genomic start will use Pfam end, while genomic end will use Pfam start.

overlap_exons_with_pfam_to_get_pfam_genomic_coord_pos_strand <- overlap_exons_with_pfam_to_get_pfam_genomic_coord[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$EXONSTRAND == "+"),]

overlap_exons_with_pfam_to_get_pfam_genomic_coord_neg_strand <- overlap_exons_with_pfam_to_get_pfam_genomic_coord[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$EXONSTRAND == "-"),]

overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.genomic.start[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$EXONSTRAND == "+")] <- overlap_exons_with_pfam_to_get_pfam_genomic_coord_pos_strand$CDSSTART + overlap_exons_with_pfam_to_get_pfam_genomic_coord_pos_strand$Pfam.this.exon.CDS.start - 1

overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.genomic.end[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$EXONSTRAND == "+")] <- overlap_exons_with_pfam_to_get_pfam_genomic_coord_pos_strand$CDSSTART + overlap_exons_with_pfam_to_get_pfam_genomic_coord_pos_strand$Pfam.this.exon.CDS.end - 1

overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.genomic.start[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$EXONSTRAND == "-")] <- overlap_exons_with_pfam_to_get_pfam_genomic_coord_neg_strand$CDSEND - overlap_exons_with_pfam_to_get_pfam_genomic_coord_neg_strand$Pfam.this.exon.CDS.end + 1

overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.genomic.end[which(overlap_exons_with_pfam_to_get_pfam_genomic_coord$EXONSTRAND == "-")] <- overlap_exons_with_pfam_to_get_pfam_genomic_coord_neg_strand$CDSEND - overlap_exons_with_pfam_to_get_pfam_genomic_coord_neg_strand$Pfam.this.exon.CDS.start + 1

#Think this is finally correct. Let's save progress again.

save.image("pfam_ids_get_genomic_coords.Rdata")
