#Load libraries.

library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)

#Placeholder - Evan to set this dynamically in the future.
installation_dir_FusionVis <- "/nethome/hmgeiger/FusionViz"

#This script will accept as arguments the full path to a fusion detection program output, and arguments for genome (GRCh37/GRCh38) and program (STAR-Fusion/FusionCatcher).
#Also an output directory.

args <- commandArgs(trailingOnly=T)

fusion_file <- args[1]
fusion_program <- args[2]
genome <- args[3]
output_dir <- args[4]

#Load different annotation Rdata files depending on genome version specified.
#Still need to make the ones for GRCh38.

if(genome == "GRCh37")
{
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_exons_as_data_frame.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_exons_as_GRanges.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_introns_as_GRanges.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_total_exons_per_transcript.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_rank_per_transcript.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_utrs_for_plotting.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_relative_exon_coords_for_all_annotated_pfam_ids.Rdata")) #This line added 02/02/18 to add loading protein domain information.
}

if(genome == "GRCh38")
{
print("GRCh38 annotation not available yet")
quit()
}

#Read in both STAR-Fusion and FusionCatcher fusion files the same way.
fusion_file <- read.table(fusion_file,header=TRUE,sep="\t",comment.char="",quote="",check.names=FALSE,stringsAsFactors=FALSE)

#If FusionCatcher, breakpoints are called "Fusion_point_for_gene_1(5end_fusion_partner)" and "Fusion_point_for_gene_2(3end_fusion_partner)".
#If STAR-Fusion, breakpoints are called "LeftBreakpoint" and "RightBreakpoint".

#Both are in format chr:breakpoint:strand.

#Other information of interest:
#Fusioncatcher geneA and geneB are called "Gene_1_symbol(5end_fusion_partner)" and "Gene_2_symbol(3end_fusion_partner)".
#Other potentially interesting information to display as text:
#"Fusion_description" - gene combination specific
#"Spanning_pairs" and "Spanning_unique_reads" - breakpoint specific
#"Predicted_effect" - breakpoint specific

#STAR-Fusion:
#"#FusionName" = GeneA--GeneB
#JunctionReadCount
#SpanningFragCount
#SpliceType

if(fusion_program == "FusionCatcher")
{
breakpoints_geneA <- fusion_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"]
breakpoints_geneB <- fusion_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"]
geneA_plus_geneB_names <- paste0(fusion_file[,"Gene_1_symbol(5end_fusion_partner)"],"--",fusion_file[,"Gene_2_symbol(3end_fusion_partner)"])
in_frame_vs_not <- rep("Not in-frame",times=nrow(fusion_file))
in_frame_vs_not[which(fusion_file[,"Predicted_effect"] == "in-frame")] <- "in-frame"
descriptions_per_fusion <- fusion_file[,"Fusion_description"]
#Here put information we do not technically need to make the plot, but that may be useful to know (like number of supporting reads).
#This contains breakpoint-specific information.
other_info_to_possibly_display <- paste0("Spanning.pairs=",fusion_file[,"Spanning_pairs"],";Spanning.unique.reads=",fusion_file[,"Spanning_unique_reads"],";Predicted.effect=",fusion_file[,"Predicted_effect"])
}

if(fusion_program == "STAR-Fusion")
{
breakpoints_geneA <- fusion_file[,"LeftBreakpoint"]
breakpoints_geneB <- fusion_file[,"RightBreakpoint"]
geneA_plus_geneB_names <- fusion_file[,"#FusionName"]
#STAR-Fusion does not output anything equivalent to in-frame from FusionCatcher. So, will not be able to color by in-frame vs. not for STAR-Fusion, so just make NA.
in_frame_vs_not <- rep(NA,times=nrow(fusion_file))
#STAR-Fusion also does not output any description, so just leave a blank for description for those outputs.
descriptions_per_fusion <- rep("No description",times=nrow(fusion_file))
#STAR-Fusion uses different definitions for read support values.
#Also has whether or not fusion has non-ref splice sites instead of predicted effect.
other_info_to_possibly_display <- paste0("JunctionReadCount=",fusion_file[,"JunctionReadCount"],";SpanningFragCount=",fusion_file[,"SpanningFragCount"],";SpliceType=",fusion_file[,"SpliceType"])
}

fusion_file_main_info <- data.frame(GeneA.GeneB = geneA_plus_geneB_names,FusionDescription = descriptions_per_fusion,GeneA.Breakpoint = breakpoints_geneA,GeneB.Breakpoint = breakpoints_geneB,In.frame.status = in_frame_vs_not,Other.breakpoint.specific.info = other_info_to_possibly_display,stringsAsFactors=FALSE)

#Get chromosome, coordinate, and strand from breakpoints_geneA and breakpoints_geneB using the strsplit function to separate the fields separated by ":".

breakpoints_geneA_data_frame <- data.frame(chr = sapply(strsplit(breakpoints_geneA,":"),"[[",1),coordinate = sapply(strsplit(breakpoints_geneA,":"),"[[",2),strand = sapply(strsplit(breakpoints_geneA,":"),"[[",3),Breakpoint.string = breakpoints_geneA,stringsAsFactors=FALSE)
breakpoints_geneA_data_frame$coordinate <- as.numeric(as.vector(breakpoints_geneA_data_frame$coordinate))

breakpoints_geneB_data_frame <- data.frame(chr = sapply(strsplit(breakpoints_geneB,":"),"[[",1),coordinate = sapply(strsplit(breakpoints_geneB,":"),"[[",2),strand = sapply(strsplit(breakpoints_geneB,":"),"[[",3),Breakpoint.string = breakpoints_geneB,stringsAsFactors=FALSE)
breakpoints_geneB_data_frame$coordinate <- as.numeric(as.vector(breakpoints_geneB_data_frame$coordinate))

#Previously, annotation had "chr" as prefix to the chromosome names.
#In this version of the script, I am now using Ensembl annotation, which does not have this.
#So we instead do the opposite, and remove "chr" prefix if it exists.

if(length(grep('chr',breakpoints_geneA_data_frame$chr[1])) == 1)
{
breakpoints_geneA_data_frame$chr <- sapply(strsplit(breakpoints_geneA_data_frame$chr,"chr"),"[[",2)
breakpoints_geneB_data_frame$chr <- sapply(strsplit(breakpoints_geneB_data_frame$chr,"chr"),"[[",2)
}

#Now, make a GRanges object out of each of the two breakpoint data frames.

breakpoints_geneA_GRanges <- GRanges(seqnames = breakpoints_geneA_data_frame$chr,
ranges = IRanges(breakpoints_geneA_data_frame$coordinate,breakpoints_geneA_data_frame$coordinate),
name = breakpoints_geneA_data_frame$Breakpoint.string,
strand = breakpoints_geneA_data_frame$strand)

breakpoints_geneB_GRanges <- GRanges(seqnames = breakpoints_geneB_data_frame$chr,
ranges = IRanges(breakpoints_geneB_data_frame$coordinate,breakpoints_geneB_data_frame$coordinate),
name = breakpoints_geneB_data_frame$Breakpoint.string,
strand = breakpoints_geneB_data_frame$strand)

#Overlap these with exon and intron annotation.

overlaps_transcript_exon_coords_and_breakpoints_geneA <- mergeByOverlaps(exons_as_GRanges,breakpoints_geneA_GRanges,maxgap=-1L)
transcript_exon_coords_overlapping_breakpoints_geneA <- BiocGenerics::as.data.frame(overlaps_transcript_exon_coords_and_breakpoints_geneA[,"exons_as_GRanges"],row.names=NULL)
breakpoint_geneA_coords_overlapping_transcript_exons <- BiocGenerics::as.data.frame(overlaps_transcript_exon_coords_and_breakpoints_geneA[,"breakpoints_geneA_GRanges"],row.names=NULL)

exon_num_per_transcript_per_breakpoint_geneA <- data.frame(Breakpoint.name = overlaps_transcript_exon_coords_and_breakpoints_geneA[,"name"], #Formerly grabbed column "name.1" instead.
TXNAME = names(overlaps_transcript_exon_coords_and_breakpoints_geneA$exons_as_GRanges), #This is now obtained a completely different way. Also now name TXNAME.
EXONRANK = as.numeric(as.vector(transcript_exon_coords_overlapping_breakpoints_geneA$exon_rank)), #Column name used to be different. Also now name EXONRANK instead of Exon.num.
EXONSTART = as.numeric(as.vector(transcript_exon_coords_overlapping_breakpoints_geneA$start)), #Same as before.
EXONEND = as.numeric(as.vector(transcript_exon_coords_overlapping_breakpoints_geneA$end)), #Same as before.
Breakpoint.coord = as.numeric(as.vector(breakpoint_geneA_coords_overlapping_transcript_exons$start)), #Same as before.
stringsAsFactors=FALSE)

overlaps_transcript_exon_coords_and_breakpoints_geneB <- mergeByOverlaps(exons_as_GRanges,breakpoints_geneB_GRanges,maxgap=-1L)
transcript_exon_coords_overlapping_breakpoints_geneB <- BiocGenerics::as.data.frame(overlaps_transcript_exon_coords_and_breakpoints_geneB[,"exons_as_GRanges"],row.names=NULL)
breakpoint_geneB_coords_overlapping_transcript_exons <- BiocGenerics::as.data.frame(overlaps_transcript_exon_coords_and_breakpoints_geneB[,"breakpoints_geneB_GRanges"],row.names=NULL)

exon_num_per_transcript_per_breakpoint_geneB <- data.frame(Breakpoint.name = overlaps_transcript_exon_coords_and_breakpoints_geneB[,"name"],
 TXNAME = names(overlaps_transcript_exon_coords_and_breakpoints_geneB$exons_as_GRanges),
 EXONRANK = as.numeric(as.vector(transcript_exon_coords_overlapping_breakpoints_geneB$exon_rank)),
 EXONSTART = as.numeric(as.vector(transcript_exon_coords_overlapping_breakpoints_geneB$start)),
 EXONEND = as.numeric(as.vector(transcript_exon_coords_overlapping_breakpoints_geneB$end)),
 Breakpoint.coord = as.numeric(as.vector(breakpoint_geneB_coords_overlapping_transcript_exons$start)),
 stringsAsFactors=FALSE)

overlaps_transcript_intron_coords_and_breakpoints_geneA <- mergeByOverlaps(introns_as_GRanges,breakpoints_geneA_GRanges,maxgap=-1L)
transcript_intron_coords_overlapping_breakpoints_geneA <- BiocGenerics::as.data.frame(overlaps_transcript_intron_coords_and_breakpoints_geneA[,"introns_as_GRanges"],row.names=NULL)
breakpoint_geneA_coords_overlapping_transcript_introns <- BiocGenerics::as.data.frame(overlaps_transcript_intron_coords_and_breakpoints_geneA[,"breakpoints_geneA_GRanges"],row.names=NULL)

intron_per_transcript_per_breakpoint_geneA <- data.frame(Breakpoint.name = overlaps_transcript_intron_coords_and_breakpoints_geneA[,"name"],
 TXNAME = names(overlaps_transcript_intron_coords_and_breakpoints_geneA$introns_as_GRanges),
 INTRONSTART = as.numeric(as.vector(transcript_intron_coords_overlapping_breakpoints_geneA$start)),
 INTRONEND = as.numeric(as.vector(transcript_intron_coords_overlapping_breakpoints_geneA$end)),
 Breakpoint.coord = as.numeric(as.vector(breakpoint_geneA_coords_overlapping_transcript_introns$start)),
 stringsAsFactors=FALSE)

overlaps_transcript_intron_coords_and_breakpoints_geneB <- mergeByOverlaps(introns_as_GRanges,breakpoints_geneB_GRanges,maxgap=-1L)
transcript_intron_coords_overlapping_breakpoints_geneB <- BiocGenerics::as.data.frame(overlaps_transcript_intron_coords_and_breakpoints_geneB[,"introns_as_GRanges"],row.names=NULL)
breakpoint_geneB_coords_overlapping_transcript_introns <- BiocGenerics::as.data.frame(overlaps_transcript_intron_coords_and_breakpoints_geneB[,"breakpoints_geneB_GRanges"],row.names=NULL)

intron_per_transcript_per_breakpoint_geneB <- data.frame(Breakpoint.name = overlaps_transcript_intron_coords_and_breakpoints_geneB[,"name"],
 TXNAME = names(overlaps_transcript_intron_coords_and_breakpoints_geneB$introns_as_GRanges),
 INTRONSTART = as.numeric(as.vector(transcript_intron_coords_overlapping_breakpoints_geneB$start)),
 INTRONEND = as.numeric(as.vector(transcript_intron_coords_overlapping_breakpoints_geneB$end)),
 Breakpoint.coord = as.numeric(as.vector(breakpoint_geneB_coords_overlapping_transcript_introns$start)),
 stringsAsFactors=FALSE)

#Now, choose the best transcript to represent each breakpoint.
#First, we need to separate breakpoints that can be represented within exons, from those that are only in introns, from those that do not overlap either (and thus are "incompatible annotation").

incompatible_annotation_rows <- c()

geneA_breakpoints_overlapping_exons <- unique(exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.name)
geneA_breakpoints_overlapping_introns_only <- setdiff(unique(intron_per_transcript_per_breakpoint_geneA$Breakpoint.name),geneA_breakpoints_overlapping_exons)
geneA_breakpoints_incompatible_annotation <- setdiff(unique(breakpoints_geneA),c(geneA_breakpoints_overlapping_exons,geneA_breakpoints_overlapping_introns_only))
if(length(geneA_breakpoints_incompatible_annotation) > 0){incompatible_annotation_rows <- c(incompatible_annotation_rows,grep(paste0(geneA_breakpoints_incompatible_annotation,collapse="|"),breakpoints_geneA))}

geneB_breakpoints_overlapping_exons <- unique(exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.name)
geneB_breakpoints_overlapping_introns_only <- setdiff(unique(intron_per_transcript_per_breakpoint_geneB$Breakpoint.name),geneB_breakpoints_overlapping_exons)
geneB_breakpoints_incompatible_annotation <- setdiff(unique(breakpoints_geneB),c(geneB_breakpoints_overlapping_exons,geneB_breakpoints_overlapping_introns_only))
if(length(geneB_breakpoints_incompatible_annotation) > 0){incompatible_annotation_rows <- c(incompatible_annotation_rows,grep(paste0(geneB_breakpoints_incompatible_annotation,collapse="|"),breakpoints_geneB))}

incompatible_annotation_rows <- unique(incompatible_annotation_rows)

#So we will not visualize for incompatible_annotation_rows. For the remainder, we will, using exons when possible.
#Choose transcript prioritizing coding transcripts, then longer CDS, then longer transcript.

exon_num_per_transcript_per_breakpoint_geneA <- merge(exon_num_per_transcript_per_breakpoint_geneA,transcript_info,by.x="TXNAME",by.y = "ensembl_transcript_id")
exon_num_per_transcript_per_breakpoint_geneB <- merge(exon_num_per_transcript_per_breakpoint_geneB,transcript_info,by.x="TXNAME",by.y = "ensembl_transcript_id")

exon_num_per_transcript_per_breakpoint_geneA <- exon_num_per_transcript_per_breakpoint_geneA[order(exon_num_per_transcript_per_breakpoint_geneA$coding,exon_num_per_transcript_per_breakpoint_geneA$cds_length,exon_num_per_transcript_per_breakpoint_geneA$transcript_length,decreasing=TRUE),]
exon_num_per_transcript_per_breakpoint_geneB <- exon_num_per_transcript_per_breakpoint_geneB[order(exon_num_per_transcript_per_breakpoint_geneB$coding,exon_num_per_transcript_per_breakpoint_geneB$cds_length,exon_num_per_transcript_per_breakpoint_geneB$transcript_length,decreasing=TRUE),]

intron_per_transcript_per_breakpoint_geneA <- merge(intron_per_transcript_per_breakpoint_geneA,data.frame(Breakpoint.name = geneA_breakpoints_overlapping_introns_only,stringsAsFactors=FALSE),by = "Breakpoint.name")
intron_per_transcript_per_breakpoint_geneB <- merge(intron_per_transcript_per_breakpoint_geneB,data.frame(Breakpoint.name = geneB_breakpoints_overlapping_introns_only,stringsAsFactors=FALSE),by = "Breakpoint.name")

intron_per_transcript_per_breakpoint_geneA <- merge(intron_per_transcript_per_breakpoint_geneA,transcript_info,by.x="TXNAME",by.y = "ensembl_transcript_id")
intron_per_transcript_per_breakpoint_geneB <- merge(intron_per_transcript_per_breakpoint_geneB,transcript_info,by.x="TXNAME",by.y = "ensembl_transcript_id")

intron_per_transcript_per_breakpoint_geneA <- intron_per_transcript_per_breakpoint_geneA[order(intron_per_transcript_per_breakpoint_geneA$coding,intron_per_transcript_per_breakpoint_geneA$cds_length,intron_per_transcript_per_breakpoint_geneA$transcript_length,decreasing=TRUE),]
intron_per_transcript_per_breakpoint_geneB <- intron_per_transcript_per_breakpoint_geneB[order(intron_per_transcript_per_breakpoint_geneB$coding,intron_per_transcript_per_breakpoint_geneB$cds_length,intron_per_transcript_per_breakpoint_geneB$transcript_length,decreasing=TRUE),]

exon_num_per_transcript_per_breakpoint_geneA <- exon_num_per_transcript_per_breakpoint_geneA[which(duplicated(exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.name) == FALSE),]
exon_num_per_transcript_per_breakpoint_geneB <- exon_num_per_transcript_per_breakpoint_geneB[which(duplicated(exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.name) == FALSE),]
intron_per_transcript_per_breakpoint_geneA <- intron_per_transcript_per_breakpoint_geneA[which(duplicated(intron_per_transcript_per_breakpoint_geneA$Breakpoint.name) == FALSE),]
intron_per_transcript_per_breakpoint_geneB <- intron_per_transcript_per_breakpoint_geneB[which(duplicated(intron_per_transcript_per_breakpoint_geneB$Breakpoint.name) == FALSE),]

#For introns, can get one coordinate of the nearest exon.
#If on positive strand, subtract one from the intron start. This will give EXONEND of the exon before.
#If on negative strand, add one to the intron end. This will give EXONSTART of the exon before.

introns_geneA_find_exon_before <- merge(intron_per_transcript_per_breakpoint_geneA,exons_as_data_frame[,c("EXONSTRAND","EXONSTART","EXONEND","EXONRANK","TXNAME")],by="TXNAME")
introns_geneB_find_exon_before <- merge(intron_per_transcript_per_breakpoint_geneB,exons_as_data_frame[,c("EXONSTRAND","EXONSTART","EXONEND","EXONRANK","TXNAME")],by="TXNAME")

introns_geneA_find_exon_before <- introns_geneA_find_exon_before[which((introns_geneA_find_exon_before$EXONSTRAND == "+" & (introns_geneA_find_exon_before$INTRONSTART - 1) == introns_geneA_find_exon_before$EXONEND) | (introns_geneA_find_exon_before$EXONSTRAND == "-" & (introns_geneA_find_exon_before$INTRONEND + 1) == introns_geneA_find_exon_before$EXONSTART)),]
introns_geneB_find_exon_before <- introns_geneB_find_exon_before[which((introns_geneB_find_exon_before$EXONSTRAND == "+" & (introns_geneB_find_exon_before$INTRONSTART - 1) == introns_geneB_find_exon_before$EXONEND) | (introns_geneB_find_exon_before$EXONSTRAND == "-" & (introns_geneB_find_exon_before$INTRONEND + 1) == introns_geneB_find_exon_before$EXONSTART)),]

#Add additional info to exon table as well.

exon_num_per_transcript_per_breakpoint_geneA <- merge(exon_num_per_transcript_per_breakpoint_geneA,exons_as_data_frame[,c("EXONSTRAND","TXNAME","EXONLENGTH","EXONRANK")],by=c("TXNAME","EXONRANK"))
exon_num_per_transcript_per_breakpoint_geneB <- merge(exon_num_per_transcript_per_breakpoint_geneB,exons_as_data_frame[,c("EXONSTRAND","TXNAME","EXONLENGTH","EXONRANK")],by=c("TXNAME","EXONRANK"))

#Switch exon start and end if on negative strand.

geneA_new_exon_start_and_end <- data.frame(EXONSTART = exon_num_per_transcript_per_breakpoint_geneA$EXONSTART,EXONEND = exon_num_per_transcript_per_breakpoint_geneA$EXONEND)
geneA_new_exon_start_and_end$EXONSTART <- ifelse(exon_num_per_transcript_per_breakpoint_geneA$EXONSTRAND == "+",exon_num_per_transcript_per_breakpoint_geneA$EXONSTART,exon_num_per_transcript_per_breakpoint_geneA$EXONEND)
geneA_new_exon_start_and_end$EXONEND <- ifelse(exon_num_per_transcript_per_breakpoint_geneA$EXONSTRAND == "+",exon_num_per_transcript_per_breakpoint_geneA$EXONEND,exon_num_per_transcript_per_breakpoint_geneA$EXONSTART)

geneB_new_exon_start_and_end <- data.frame(EXONSTART = exon_num_per_transcript_per_breakpoint_geneB$EXONSTART,EXONEND = exon_num_per_transcript_per_breakpoint_geneB$EXONEND)
geneB_new_exon_start_and_end$EXONSTART <- ifelse(exon_num_per_transcript_per_breakpoint_geneB$EXONSTRAND == "+",exon_num_per_transcript_per_breakpoint_geneB$EXONSTART,exon_num_per_transcript_per_breakpoint_geneB$EXONEND)
#exon_num_per_transcript_per_breakpoint_geneB$EXONEND <- ifelse(exon_num_per_transcript_per_breakpoint_geneB$EXONSTRAND == "+",exon_num_per_transcript_per_breakpoint_geneB$EXONEND,exon_num_per_transcript_per_breakpoint_geneB$EXONSTART)
geneB_new_exon_start_and_end$EXONEND <- ifelse(exon_num_per_transcript_per_breakpoint_geneB$EXONSTRAND == "+",exon_num_per_transcript_per_breakpoint_geneB$EXONEND,exon_num_per_transcript_per_breakpoint_geneB$EXONSTART) #This line added 02/02/18 to replace line just before it, which was a bug.

#Now, determine start, middle, or end and percent into exon based on distance from start or end.

geneA_matching_exon_start_middle_end <- rep("Middle",times=nrow(exon_num_per_transcript_per_breakpoint_geneA))
geneB_matching_exon_start_middle_end <- rep("Middle",times=nrow(exon_num_per_transcript_per_breakpoint_geneB))

geneA_matching_exon_start_middle_end[which(exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.coord == geneA_new_exon_start_and_end$EXONSTART)] <- "Start"
geneA_matching_exon_start_middle_end[which(exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.coord == geneA_new_exon_start_and_end$EXONEND)] <- "End"

geneB_matching_exon_start_middle_end[which(exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.coord == geneB_new_exon_start_and_end$EXONSTART)] <- "Start"
geneB_matching_exon_start_middle_end[which(exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.coord == geneB_new_exon_start_and_end$EXONEND)] <- "End"

geneA_length_into_exon <- abs(exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.coord - geneA_new_exon_start_and_end$EXONSTART) + 1
geneB_length_into_exon <- abs(exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.coord - geneB_new_exon_start_and_end$EXONSTART) + 1

geneA_length_into_exon[which(geneA_matching_exon_start_middle_end != "Middle")] <- NA
geneB_length_into_exon[which(geneB_matching_exon_start_middle_end != "Middle")] <- NA

geneA_percent_into_exon <- geneA_length_into_exon/exon_num_per_transcript_per_breakpoint_geneA$EXONLENGTH
geneB_percent_into_exon <- geneB_length_into_exon/exon_num_per_transcript_per_breakpoint_geneB$EXONLENGTH

exon_num_per_transcript_per_breakpoint_geneA <- data.frame(Breakpoint = exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.coord,
Transcript = exon_num_per_transcript_per_breakpoint_geneA$TXNAME,
ExonicOrNot = "yes",
ExonInOrBefore = exon_num_per_transcript_per_breakpoint_geneA$EXONRANK,
Breakpoint.name = exon_num_per_transcript_per_breakpoint_geneA$Breakpoint.name,
Transcript.other.info = paste0("Transcript_gene_name=",exon_num_per_transcript_per_breakpoint_geneA$external_gene_name,";Transcript_biotype=",exon_num_per_transcript_per_breakpoint_geneA$transcript_biotype,";Transcript_CDS_length=",exon_num_per_transcript_per_breakpoint_geneA$cds_length,";Transcript_total_length=",exon_num_per_transcript_per_breakpoint_geneA$transcript_length),
StartMiddleEnd = geneA_matching_exon_start_middle_end,LengthIntoExon = geneA_length_into_exon,PercentIntoExon = geneA_percent_into_exon,
stringsAsFactors=FALSE)

exon_num_per_transcript_per_breakpoint_geneA <- merge(exon_num_per_transcript_per_breakpoint_geneA,total_exons_per_transcript,by.x="Transcript",by.y="TXNAME")

exon_num_per_transcript_per_breakpoint_geneB <- data.frame(Breakpoint = exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.coord,
Transcript = exon_num_per_transcript_per_breakpoint_geneB$TXNAME,
ExonicOrNot = "yes",
ExonInOrBefore = exon_num_per_transcript_per_breakpoint_geneB$EXONRANK,
Breakpoint.name = exon_num_per_transcript_per_breakpoint_geneB$Breakpoint.name,
Transcript.other.info = paste0("Transcript_gene_name=",exon_num_per_transcript_per_breakpoint_geneB$external_gene_name,";Transcript_biotype=",exon_num_per_transcript_per_breakpoint_geneB$transcript_biotype,";Transcript_CDS_length=",exon_num_per_transcript_per_breakpoint_geneB$cds_length,";Transcript_total_length=",exon_num_per_transcript_per_breakpoint_geneB$transcript_length),
StartMiddleEnd = geneB_matching_exon_start_middle_end,LengthIntoExon = geneB_length_into_exon,PercentIntoExon = geneB_percent_into_exon,
stringsAsFactors=FALSE)

exon_num_per_transcript_per_breakpoint_geneB <- merge(exon_num_per_transcript_per_breakpoint_geneB,total_exons_per_transcript,by.x="Transcript",by.y="TXNAME")

#Set table up similarly for introns.

introns_geneA_find_exon_before <- data.frame(Breakpoint = introns_geneA_find_exon_before$Breakpoint.coord,Transcript = introns_geneA_find_exon_before$TXNAME,ExonicOrNot = "no",ExonInOrBefore = introns_geneA_find_exon_before$EXONRANK,Breakpoint.name = introns_geneA_find_exon_before$Breakpoint.name,Transcript.other.info = paste0("Transcript_gene_name=",introns_geneA_find_exon_before$external_gene_name,";Transcript_biotype=",introns_geneA_find_exon_before$transcript_biotype,";Transcript_CDS_length=",introns_geneA_find_exon_before$cds_length,";Transcript_total_length=",introns_geneA_find_exon_before$transcript_length),StartMiddleEnd = "Intronic",LengthIntoExon = NA,PercentIntoExon = NA,stringsAsFactors=FALSE)

introns_geneB_find_exon_before <- data.frame(Breakpoint = introns_geneB_find_exon_before$Breakpoint.coord,Transcript = introns_geneB_find_exon_before$TXNAME,ExonicOrNot = "no",ExonInOrBefore = introns_geneB_find_exon_before$EXONRANK,Breakpoint.name = introns_geneB_find_exon_before$Breakpoint.name,Transcript.other.info = paste0("Transcript_gene_name=",introns_geneB_find_exon_before$external_gene_name,";Transcript_biotype=",introns_geneB_find_exon_before$transcript_biotype,";Transcript_CDS_length=",introns_geneB_find_exon_before$cds_length,";Transcript_total_length=",introns_geneB_find_exon_before$transcript_length),StartMiddleEnd = "Intronic",LengthIntoExon = NA,PercentIntoExon = NA,stringsAsFactors=FALSE)

introns_geneA_find_exon_before <- merge(introns_geneA_find_exon_before,total_exons_per_transcript,by.x="Transcript",by.y="TXNAME")
introns_geneB_find_exon_before <- merge(introns_geneB_find_exon_before,total_exons_per_transcript,by.x="Transcript",by.y="TXNAME")

#Rbind the exon and intron tables.

geneA_representative_transcripts_per_breakpoint <- rbind(exon_num_per_transcript_per_breakpoint_geneA,introns_geneA_find_exon_before)
geneB_representative_transcripts_per_breakpoint <- rbind(exon_num_per_transcript_per_breakpoint_geneB,introns_geneB_find_exon_before)

#Bind the match-up to transcripts to the other info.

fusion_file_main_info_incompatible_annotation <- fusion_file_main_info[incompatible_annotation_rows,]

fusion_file_main_info <- fusion_file_main_info[setdiff(1:nrow(fusion_file_main_info),incompatible_annotation_rows),]

fusion_file_main_info_match_up_geneA <- merge(fusion_file_main_info,geneA_representative_transcripts_per_breakpoint,by.x = "GeneA.Breakpoint",by.y = "Breakpoint.name")
fusion_file_main_info_match_up_geneB <- merge(fusion_file_main_info,geneB_representative_transcripts_per_breakpoint,by.x = "GeneB.Breakpoint",by.y = "Breakpoint.name")

#Last minor detail: multiply PercentIntoExon by 100, as it is actually currently a proportion in range 0-1.

fusion_file_main_info_match_up_geneA$PercentIntoExon <- fusion_file_main_info_match_up_geneA$PercentIntoExon*100
fusion_file_main_info_match_up_geneB$PercentIntoExon <- fusion_file_main_info_match_up_geneB$PercentIntoExon*100

#Also sort geneA and geneB info so they match up.

fusion_file_main_info_match_up_geneA <- fusion_file_main_info_match_up_geneA[match(paste0(fusion_file_main_info$GeneA.Breakpoint,"/",fusion_file_main_info$GeneB.Breakpoint),paste0(fusion_file_main_info_match_up_geneA$GeneA.Breakpoint,"/",fusion_file_main_info_match_up_geneA$GeneB.Breakpoint)),]
fusion_file_main_info_match_up_geneB <- fusion_file_main_info_match_up_geneB[match(paste0(fusion_file_main_info$GeneA.Breakpoint,"/",fusion_file_main_info$GeneB.Breakpoint),paste0(fusion_file_main_info_match_up_geneB$GeneA.Breakpoint,"/",fusion_file_main_info_match_up_geneB$GeneB.Breakpoint)),]

#save(list=c("fusion_file","fusion_program","genome","output_dir","fusion_file_main_info_incompatible_annotation","fusion_file_main_info","fusion_file_main_info_match_up_geneA","fusion_file_main_info_match_up_geneB"),file=paste0(output_dir,"/overlaps_breakpoints_and_transcripts.Rdata"))
#save(list=c("fusion_file","fusion_program","genome","output_dir","fusion_file_main_info_incompatible_annotation","fusion_file_main_info","fusion_file_main_info_match_up_geneA","fusion_file_main_info_match_up_geneB","exon_num_per_pfam_domain_start","exon_num_per_pfam_domain_end","utrs_data_frame_for_plotting"),
#file=paste0(output_dir,"/overlaps_breakpoints_and_transcripts.Rdata"))

#Final step in annotation is to match up UTR table.
#utrs_data_frame_for_plotting contains 5' and 3' UTRs for all annotated transcripts.
#Subset it to include only transcripts being used to represent our fusions.

fiveprime_utrs <- utrs_data_frame_for_plotting[order(utrs_data_frame_for_plotting$ExonInOrBefore),]
fiveprime_utrs <- fiveprime_utrs[which(duplicated(fiveprime_utrs$Transcript) == FALSE),]
rownames(fiveprime_utrs) <- fiveprime_utrs$Transcript

threeprime_utrs <- utrs_data_frame_for_plotting[order(utrs_data_frame_for_plotting$ExonInOrBefore,decreasing=TRUE),]
threeprime_utrs <- threeprime_utrs[which(duplicated(threeprime_utrs$Transcript) == FALSE),]
rownames(threeprime_utrs) <- threeprime_utrs$Transcript

fiveprime_utrs <- fiveprime_utrs[unique(c(fusion_file_main_info_match_up_geneA$Transcript,fusion_file_main_info_match_up_geneB$Transcript)),]
threeprime_utrs <- threeprime_utrs[unique(c(fusion_file_main_info_match_up_geneA$Transcript,fusion_file_main_info_match_up_geneB$Transcript)),]

#Now ready to save all objects.
#This way if want to run just the plotting step after this, can load overlaps_breakpoints_and_transcripts.Rdata and have all the objects we need.

save(list=c("fusion_file","fusion_program","genome","output_dir","fusion_file_main_info_incompatible_annotation","fusion_file_main_info","fusion_file_main_info_match_up_geneA","fusion_file_main_info_match_up_geneB",
"exon_num_per_pfam_domain_start","exon_num_per_pfam_domain_end","fiveprime_utrs","threeprime_utrs"),
file=paste0(output_dir,"/overlaps_breakpoints_and_transcripts.Rdata"))

#Commenting this out for now, but here are some lines to better understand the objects used for plotting.

#dim(fusion_file_main_info)
#dim(fusion_file_main_info_match_up_geneA)
#dim(fusion_file_main_info_match_up_geneB)

#head(fusion_file_main_info,n=3)
#head(fusion_file_main_info_match_up_geneA,n=3)
#head(fusion_file_main_info_match_up_geneB,n=3)

#length(unique(fusion_file_main_info_match_up_geneA$Transcript))
#length(unique(fusion_file_main_info_match_up_geneB$Transcript))
#length(unique(c(fusion_file_main_info_match_up_geneA$Transcript,fusion_file_main_info_match_up_geneB$Transcript)))

#dim(fiveprime_utrs)
#dim(threeprime_utrs)

#head(fiveprime_utrs,n=3)
#head(threeprime_utrs,n=3)

#dim(exon_num_per_pfam_domain_start)
#dim(exon_num_per_pfam_domain_end)

#head(exon_num_per_pfam_domain_start,n=3)
#head(exon_num_per_pfam_domain_end,n=3)

#Now, we source the functions from the script that contains the functions specifically for plotting.
#See that script for more detailed comments.

source('scripts/plot_fusions_incl_domains_after_process_and_match_up_to_annotation.R')

#Run through each row that does not have descriptions we want to remove, and put the row number into the function final_plot_for_given_row_num.
#As we do this, output to PDF.

descriptions_to_remove <- c("readthrough","healthy","banned")

row_numbers_to_plot <- grep(paste(descriptions_to_remove,collapse="|"),fusion_file_main_info$FusionDescription,perl=TRUE,invert=TRUE)

pdf(paste0(output_dir,"/fusion_visualizations.pdf"),width=12,height=8)

for(rownum in row_numbers_to_plot)
{
final_plot_for_given_row_num(rownum)
}

dev.off()

