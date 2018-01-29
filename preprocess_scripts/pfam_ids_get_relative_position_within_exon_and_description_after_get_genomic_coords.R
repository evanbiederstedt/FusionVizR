#Load libraries and Rdata.

library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)

load("pfam_ids_get_genomic_coords.Rdata")

#This script is a followup script to pfam_ids_get_genomic_coords_not_incl_get_relative_position_within_exon.R.
#Here, we combine Pfam.genomic.start and Pfam.genomic.end information with whatever info we want from exons_overlapping_pfam_domains_from_CDS_coord_within_transcript and pfam_domains_overlapping_exons_from_CDS_coord_within_transcript.
#We also get descriptions for the Pfam IDs.

#Let's start by removing some columns from the intersection of the exons with Pfam domains.

exons_overlapping_pfam_domains_from_CDS_coord_within_transcript <- exons_overlapping_pfam_domains_from_CDS_coord_within_transcript[,grep('mcols',colnames(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript))]

colnames(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript) <- sapply(strsplit(colnames(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript),"mcols."),"[[",2)

pfam_domains_overlapping_exons_from_CDS_coord_within_transcript <- pfam_domains_overlapping_exons_from_CDS_coord_within_transcript[,grep('mcols',colnames(pfam_domains_overlapping_exons_from_CDS_coord_within_transcript))]

colnames(pfam_domains_overlapping_exons_from_CDS_coord_within_transcript) <- sapply(strsplit(colnames(pfam_domains_overlapping_exons_from_CDS_coord_within_transcript),"mcols."),"[[",2)

#Now, combine the information in these two tables with the Pfam domain genomic coordinates.

pfam_match_up_to_exons_plus_genomic_coords <- data.frame(exons_overlapping_pfam_domains_from_CDS_coord_within_transcript,
pfam_domains_overlapping_exons_from_CDS_coord_within_transcript,
Pfam.genomic.start = overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.genomic.start,
Pfam.genomic.end = overlap_exons_with_pfam_to_get_pfam_genomic_coord$Pfam.genomic.end,
stringsAsFactors=FALSE)

#Finally, let's get descriptions per Pfam ID using the PFAM.db package. Then match these up with the table.

library(PFAM.db)
pfam_id_vs_description <- data.frame(PFAMID = keys(PFAMDE),Description = as.vector(unlist(mget(keys(PFAMDE),PFAMDE))),stringsAsFactors=FALSE)

#pfam_id_vs_description includes all Pfam IDs (including some that may be non-human).
#So will have IDs not found in pfam_match_up_to_exons_plus_genomic_coords.
#Question is, will all IDs in pfam_match_up_to_exons_plus_genomic_coords be found in pfam_id_vs_description? 
#Let's check.

#length(unique(pfam_match_up_to_exons_plus_genomic_coords$pfam))
#length(intersect(unique(pfam_match_up_to_exons_plus_genomic_coords$pfam),pfam_id_vs_description$PFAMID))

#Looks like there are a few Pfam IDs not in pfam_id_vs_description, but not too many. Probably a version issue.
#Let's add these to the table with description "No description in Jan. 2018 version of Pfam DB".

pfam_id_vs_description <- rbind(pfam_id_vs_description,data.frame(PFAMID = setdiff(unique(pfam_match_up_to_exons_plus_genomic_coords$pfam),pfam_id_vs_description$PFAMID),Description = "No description in Jan. 2018 version of Pfam DB",stringsAsFactors=FALSE))

colnames(pfam_id_vs_description)[1] <- "pfam"
colnames(pfam_id_vs_description)[2] <- "Pfam.description"

#Add domain descriptions to pfam_match_up_to_exons_plus_genomic_coords.

pfam_match_up_to_exons_plus_genomic_coords <- merge(pfam_match_up_to_exons_plus_genomic_coords,pfam_id_vs_description,by="pfam")
pfam_match_up_to_exons_plus_genomic_coords <- pfam_match_up_to_exons_plus_genomic_coords[order(pfam_match_up_to_exons_plus_genomic_coords$TXNAME,pfam_match_up_to_exons_plus_genomic_coords$EXONRANK),]

#Now, we only want to keep rows in pfam_match_up_to_exons_plus_genomic_coords corresponding to domain boundaries.
#For example, if a domain spans exons 2-8, keep only exons 2 and 8 and discard 3-7.

#Now, we may sometimes have multiple occurences of the same domain for the same transcript.
#For example, ENST00000001008, domain PF00254 occurs from amino acids 45-134 and 162-249. We do NOT want to merge these two intervals.

#Create a dummary variable column with TXNAME, pfam, and pfam_start pasted together.
#Then, regular and reverse sort by EXONRANK, and take only the first occurence of the dummy variable.

pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- pfam_match_up_to_exons_plus_genomic_coords
pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- pfam_match_up_to_exons_plus_genomic_coords

pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- pfam_match_up_to_exons_plus_genomic_coords_domain_starts[order(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONRANK),]
pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- pfam_match_up_to_exons_plus_genomic_coords_domain_starts[which(duplicated(paste0(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$TXNAME,"/",pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam,"/",pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam_start)) == FALSE),]

pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- pfam_match_up_to_exons_plus_genomic_coords_domain_ends[order(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONRANK,decreasing=TRUE),]
pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- pfam_match_up_to_exons_plus_genomic_coords_domain_ends[which(duplicated(paste0(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$TXNAME,"/",pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam,"/",pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam_start)) == FALSE),]

#pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- pfam_match_up_to_exons_plus_genomic_coords_domain_starts[order(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$TXNAME,pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONRANK),]
#pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- pfam_match_up_to_exons_plus_genomic_coords_domain_ends[order(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$TXNAME,pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONRANK),]

pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- pfam_match_up_to_exons_plus_genomic_coords_domain_starts[order(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$TXNAME,pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam,pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam_start),]
pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- pfam_match_up_to_exons_plus_genomic_coords_domain_ends[order(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$TXNAME,pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam,pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam_start),]

#Now we are going to create a new variable Pfam.genomic.single.coord in pfam_match_up_to_exons_plus_genomic_coords_domain_starts and pfam_match_up_to_exons_plus_genomic_coords_domain_ends.
#Set it based on either Pfam.genomic.start or Pfam.genomic.end based on strand information.

pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- data.frame(pfam_match_up_to_exons_plus_genomic_coords_domain_starts,
Pfam.genomic.single.coord = NA,
stringsAsFactors=FALSE)

pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- data.frame(pfam_match_up_to_exons_plus_genomic_coords_domain_ends,
Pfam.genomic.single.coord = NA,
stringsAsFactors=FALSE)

pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Pfam.genomic.single.coord <- ifelse(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONSTRAND == "+",
pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Pfam.genomic.start,
pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Pfam.genomic.end)

pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Pfam.genomic.single.coord <- ifelse(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONSTRAND == "+",
pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Pfam.genomic.end,
pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Pfam.genomic.start)

#Add another column "Strand.corrected.exon.start" which will contain the EXONSTART adjusted for strand.

pfam_match_up_to_exons_plus_genomic_coords_domain_starts <- data.frame(pfam_match_up_to_exons_plus_genomic_coords_domain_starts,
Strand.corrected.exon.start = ifelse(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONSTRAND == "+",
pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONSTART,
pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONEND),
stringsAsFactors=FALSE)

pfam_match_up_to_exons_plus_genomic_coords_domain_ends <- data.frame(pfam_match_up_to_exons_plus_genomic_coords_domain_ends,
Strand.corrected.exon.start = ifelse(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONSTRAND == "+",
pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONSTART,
pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONEND),
stringsAsFactors=FALSE)

#Final step after this will be to check distance of Pfam.genomic.single.coord from Strand.corrected.exon.start (abs(Pfam.genomic.single.coord - Strand.corrected.exon.start).
#If 0, then start.
#If abs(Pfam.genomic.single.coord - Strand.corrected.exon.start) + 1 = exon length, then end.
#Otherwise length into exon is abs(Pfam.genomic.single.coord - Strand.corrected.exon.start) + 1, then divide by exon length to get percent into exon.

#Then set up a table similar to this one:

#exon_num_per_transcript_per_breakpoint_geneA <- data.frame(
#Breakpoint = exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.coord,
#Transcript = exon_num_per_transcript_per_breakpoint_geneA$TXNAME,
#ExonicOrNot = "yes",
#ExonInOrBefore = exon_num_per_transcript_per_breakpoint_geneA$EXONRANK,
#Breakpoint.name = exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.name,
#StartMiddleEnd = geneA_matching_exon_start_middle_end,LengthIntoExon = geneA_length_into_exon,PercentIntoExon = geneA_percent_into_exon,
#Transcript.other.info = paste0("Transcript_gene_name=",exon_num_per_transcript_per_breakpoint_geneA$external_gene_name,";Transcript_biotype=",exon_num_per_transcript_per_breakpoint_geneA$transcript_biotype,";Transcript_CDS_length=",exon_num_per_transcript_per_breakpoint_geneA$cds_length,";Transcript_total_length=",exon_num_per_transcript_per_breakpoint_geneA$transcript_length),
#stringsAsFactors=FALSE)

pfam_domain_starts_start_middle_end <- rep("Middle",times=nrow(pfam_match_up_to_exons_plus_genomic_coords_domain_starts))
pfam_domain_ends_start_middle_end <- rep("Middle",times=nrow(pfam_match_up_to_exons_plus_genomic_coords_domain_ends))

pfam_domain_starts_lengths_into_exon <- abs(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Pfam.genomic.single.coord - pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Strand.corrected.exon.start) + 1
pfam_domain_starts_start_middle_end <- ifelse(pfam_domain_starts_lengths_into_exon == 1,"Start","Middle")
pfam_domain_starts_start_middle_end <- ifelse(pfam_domain_starts_lengths_into_exon == pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONLENGTH,"End",pfam_domain_starts_start_middle_end)

pfam_domain_ends_lengths_into_exon <- abs(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Pfam.genomic.single.coord - pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Strand.corrected.exon.start) + 1
pfam_domain_ends_start_middle_end <- ifelse(pfam_domain_ends_lengths_into_exon == 1,"Start","Middle")
pfam_domain_ends_start_middle_end <- ifelse(pfam_domain_ends_lengths_into_exon == pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONLENGTH,"End",pfam_domain_ends_start_middle_end)

exon_num_per_pfam_domain_start <- data.frame(
Breakpoint = pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Pfam.genomic.single.coord,
Transcript = pfam_match_up_to_exons_plus_genomic_coords_domain_starts$TXNAME,
ExonicOrNot = "yes",
ExonInOrBefore = pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONRANK,
Breakpoint.name = paste0(pfam_match_up_to_exons_plus_genomic_coords_domain_starts$TXNAME,"/",pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam,"(aa",pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam_start,"-aa",pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam_end,")"),
StartMiddleEnd = pfam_domain_starts_start_middle_end,
LengthIntoExon = pfam_domain_starts_lengths_into_exon,
PercentIntoExon = ifelse(pfam_domain_starts_start_middle_end == "Middle",
pfam_domain_starts_lengths_into_exon*100/pfam_match_up_to_exons_plus_genomic_coords_domain_starts$EXONLENGTH,
NA),
stringsAsFactors=FALSE)

exon_num_per_pfam_domain_end <- data.frame(
Breakpoint = pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Pfam.genomic.single.coord,
Transcript = pfam_match_up_to_exons_plus_genomic_coords_domain_ends$TXNAME,
ExonicOrNot = "yes",
ExonInOrBefore = pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONRANK,
Breakpoint.name = paste0(pfam_match_up_to_exons_plus_genomic_coords_domain_ends$TXNAME,"/",pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam,"(aa",pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam_start,"-aa",pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam_end,")"),
StartMiddleEnd = pfam_domain_ends_start_middle_end,
LengthIntoExon = pfam_domain_ends_lengths_into_exon,
PercentIntoExon = ifelse(pfam_domain_ends_start_middle_end == "Middle",
pfam_domain_ends_lengths_into_exon*100/pfam_match_up_to_exons_plus_genomic_coords_domain_ends$EXONLENGTH,
NA),
stringsAsFactors=FALSE)

#Now, we just want to add the Pfam ID and description to exon_num_per_pfam_domain_start and exon_num_per_pfam_domain_end.

exon_num_per_pfam_domain_start <- data.frame(exon_num_per_pfam_domain_start,
Pfam.ID = pfam_match_up_to_exons_plus_genomic_coords_domain_starts$pfam,
Pfam.description = pfam_match_up_to_exons_plus_genomic_coords_domain_starts$Pfam.description,
stringsAsFactors=FALSE)

exon_num_per_pfam_domain_end <- data.frame(exon_num_per_pfam_domain_end,
Pfam.ID = pfam_match_up_to_exons_plus_genomic_coords_domain_ends$pfam,
Pfam.description = pfam_match_up_to_exons_plus_genomic_coords_domain_ends$Pfam.description,
stringsAsFactors=FALSE)

#Now, we have a table each for domain start and end that looks like this:

#  Breakpoint      Transcript ExonicOrNot ExonInOrBefore
#1  127228571 ENST00000000233         yes              1
#2    9099000 ENST00000000412         yes              2
#3   64082314 ENST00000000442         yes              5
#                       Breakpoint.name StartMiddleEnd LengthIntoExon
#1   ENST00000000233/PF00025(aa7-aa176)         Middle            173
#2   ENST00000000412/PF02157(aa1-aa277)         Middle              2
#3 ENST00000000442/PF00104(aa225-aa398)         Middle            102
#  PercentIntoExon Pfam.ID                                  Pfam.description
#1       78.280543 PF00025                    ADP-ribosylation factor family
#2        1.129944 PF02157                      Mannose-6-phosphate receptor
#3       59.649123 PF00104 Ligand-binding domain of nuclear hormone receptor

#Save objects.

save.image("pfam_ids_get_relative_position_within_exon_and_description.Rdata")

#Also, index the unique domains within each transcript by amino acid start.
#For example if a transcript has domain A from amino acids 1-5 and 12-15 and B from 6-10, A would be index 1.
#Index 1 = A, index 2 = B, and so on.

#First, for each transcript/domain combo in exon_num_per_pfam_domain_start, get only the rows where ExonInOrBefore is lowest.

temp_exon_num_per_pfam_domain_start <- data.frame(exon_num_per_pfam_domain_start,
transcript.plus.domain = paste0(exon_num_per_pfam_domain_start$Transcript,"/",exon_num_per_pfam_domain_start$Pfam.ID),
stringsAsFactors=FALSE)

temp_exon_num_per_pfam_domain_start <- temp_exon_num_per_pfam_domain_start[order(temp_exon_num_per_pfam_domain_start$Transcript,temp_exon_num_per_pfam_domain_start$ExonInOrBefore,temp_exon_num_per_pfam_domain_start$LengthIntoExon),]

#Now, take only unique occurences of transcript/domain to get the first instance of each domain.

temp_exon_num_per_pfam_domain_start <- temp_exon_num_per_pfam_domain_start[which(duplicated(temp_exon_num_per_pfam_domain_start$transcript.plus.domain) == FALSE),]

#Create a running count of each transcript ID in a new column.

temp_running_count_transcript <- aggregate(rep(1,times=nrow(temp_exon_num_per_pfam_domain_start)) ~ temp_exon_num_per_pfam_domain_start$Transcript,FUN=cumsum)

temp_running_count_transcript <- as.vector(unlist(temp_running_count_transcript[,2]))

temp_exon_num_per_pfam_domain_start$Pfam.short.label <- LETTERS[temp_running_count_transcript]

#Add these short labels to the original table.

rm(temp_running_count_transcript)

exon_num_per_pfam_domain_start <- merge(exon_num_per_pfam_domain_start,temp_exon_num_per_pfam_domain_start[,c("Transcript","Pfam.ID","Pfam.short.label")],by=c("Transcript","Pfam.ID"))
exon_num_per_pfam_domain_start <- exon_num_per_pfam_domain_start[order(exon_num_per_pfam_domain_start$Transcript,exon_num_per_pfam_domain_start$ExonInOrBefore,exon_num_per_pfam_domain_start$LengthIntoExon),]

exon_num_per_pfam_domain_end <- merge(exon_num_per_pfam_domain_end,temp_exon_num_per_pfam_domain_start[,c("Transcript","Pfam.ID","Pfam.short.label")],by=c("Transcript","Pfam.ID"))
exon_num_per_pfam_domain_end <- exon_num_per_pfam_domain_end[order(exon_num_per_pfam_domain_end$Transcript,exon_num_per_pfam_domain_end$ExonInOrBefore,exon_num_per_pfam_domain_end$LengthIntoExon),]

#Save again.

save.image("pfam_ids_get_relative_position_within_exon_and_description.Rdata")

#Also save just the last objects we made, as these are the only ones we need to plot.

save(list=c("exon_num_per_pfam_domain_start","exon_num_per_pfam_domain_end"),file="annotation/relative_exon_coords_for_all_annotated_pfam_ids.Rdata")
