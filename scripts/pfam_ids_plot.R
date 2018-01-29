#Evan to set this dynamically in the future.
#For now, set to current working directory, and have to run script from where it is installed.

installation_dir_FusionVis <- getwd()

#First, load libraries and Rdata files.

library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)

#Next, load annotation files.
#Here, we load one file we already had before (biomart_db_GRCh37_total_exons_per_transcript.Rdata),
#along with a new one (biomart_db_GRCh37_relative_exon_coords_for_all_annotated_pfam_ids.Rdata).

#The new annotation file gives coordinates within exons of transcripts for Pfam domains.
#This file was created with these scripts:
#preprocess_scripts/pfam_ids_get_genomic_coords_not_incl_get_relative_position_within_exon.R
#preprocess_scripts/pfam_ids_get_relative_position_within_exon_and_description_after_get_genomic_coords.R

#Will document those scripts in more detail later on to make sure their output is correct.
#For now, we trust the objects in biomart_db_GRCh37_relative_exon_coords_for_all_annotated_pfam_ids.Rdata.

load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_relative_exon_coords_for_all_annotated_pfam_ids.Rdata"))
load(paste0(installation_dir_FusionVis,"/annotation/biomart_db_GRCh37_total_exons_per_transcript.Rdata"))

#Let's look at what objects we have now.

#ls()
#[1] "exon_num_per_pfam_domain_end"   "exon_num_per_pfam_domain_start"
#[3] "installation_dir_FusionVis"     "total_exons_per_transcript"  

#Look at each object.

#head(exon_num_per_pfam_domain_start[,c(1,2,5:7,9:11)],n=4)

#       Transcript Pfam.ID ExonInOrBefore                      Breakpoint.name
#1 ENST00000000233 PF00025              1   ENST00000000233/PF00025(aa7-aa176)
#2 ENST00000000412 PF02157              2   ENST00000000412/PF02157(aa1-aa277)
#4 ENST00000000442 PF00105              2  ENST00000000442/PF00105(aa78-aa145)
#3 ENST00000000442 PF00104              5 ENST00000000442/PF00104(aa225-aa398)
#  StartMiddleEnd PercentIntoExon
#1         Middle       78.280543
#2         Middle        1.129944
#4         Middle       72.403561
#3         Middle       59.649123
#                                   Pfam.description Pfam.short.label
#1                    ADP-ribosylation factor family                A
#2                      Mannose-6-phosphate receptor                A
#4                Zinc finger, C4 type (two domains)                A
#3 Ligand-binding domain of nuclear hormone receptor                B

# head(exon_num_per_pfam_domain_end[,c(1,2,5:7,9:11)],n=4)

#       Transcript Pfam.ID ExonInOrBefore                      Breakpoint.name
#1 ENST00000000233 PF00025              6   ENST00000000233/PF00025(aa7-aa176)
#2 ENST00000000412 PF02157              7   ENST00000000412/PF02157(aa1-aa277)
#4 ENST00000000442 PF00105              3  ENST00000000442/PF00105(aa78-aa145)
#3 ENST00000000442 PF00104              7 ENST00000000442/PF00104(aa225-aa398)
#  StartMiddleEnd PercentIntoExon
#1         Middle       14.604462
#2         Middle        7.614213
#4         Middle       94.017094
#3         Middle       17.635659
 #                                  Pfam.description Pfam.short.label
#1                    ADP-ribosylation factor family                A
#2                      Mannose-6-phosphate receptor                A
#4                Zinc finger, C4 type (two domains)                A
#3 Ligand-binding domain of nuclear hormone receptor                B

#To plot this, we will go row-by-row for a given transcript and draw a rectangle from the exon coordinates specified by exon_num_per_pfam_domain_start to those specified in exon_num_per_pfam_domain_end.

#E.g. transcript ENST00000000442 will have rectangles going from 72.4% into exon 2 to 94% into exon 3, and another rectangle (with a different color since different domain) going from 59.65% into exon 5 to 17.64% into exon 7.

#Now, let's write functions to set up the fusion plot (draw boxes for number of exons).
#Then, another function to get the x coordinate on this plot given the "breakpoint info" table.
#These functions are similar if not identical to those in visualize_incl_GenomicRanges.R.

setup_fusion_plot <- function(exon_num){
plot_xmax <- 120*(exon_num - 1)+200

plot(0,pch='',ylab='',xlab='',xlim=c(0,plot_xmax),ylim=c(0,10),xaxt = 'n',yaxt = 'n',bty='n')

for(exon in 0:(exon_num - 1))
{
rect(120*exon,8.5,(120*exon)+100,9.25)
end_this_exon <- (120*exon)+100
if(exon < (exon_num - 1)){lines(c(end_this_exon,end_this_exon+20),c(8.875,8.875),type="l")}
}

}

x_coord_on_fusion_plot <- function(breakpoint_info_data_frame){

exon_index <- breakpoint_info_data_frame$ExonInOrBefore - 1
breakpoint_type <- breakpoint_info_data_frame$StartMiddleEnd

exon_of_interest_start_coordinate <- 120*exon_index

if(breakpoint_type == "Start"){breakpoint_x_coordinate_on_plot <- exon_of_interest_start_coordinate}
if(breakpoint_type == "End"){breakpoint_x_coordinate_on_plot <- exon_of_interest_start_coordinate + 100}
if(breakpoint_type == "Intronic"){breakpoint_x_coordinate_on_plot <- exon_of_interest_start_coordinate + 110}
if(breakpoint_type == "Middle")
{
breakpoint_x_coordinate_on_plot <- exon_of_interest_start_coordinate + (breakpoint_info_data_frame$PercentIntoExon/100)*101
}

return(breakpoint_x_coordinate_on_plot)

}

#Let's plot for an example transcript.

pdf("test_data_and_output/Pfam_IDs_visualization_using_new_Pfam_annotation_01.28.18.pdf")

#First, set up a colorblind-friendly vector to differentiate the different protein domains.
#In this annotation, transcripts have at most 9 unique domains, so this should be enough colors.

mycol <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")

#We chose this transcript because it has multiple unique domains, one of which occurs in more than one place in the transcript.

transcript <- "ENST00000011653"

#Draw boxes for the exons, getting the total number of exons from total_exons_per_transcript.

setup_fusion_plot(total_exons_per_transcript[which(total_exons_per_transcript$TXNAME == transcript),"TotalExonNum"])

#Get domain info for just this transcript in a separate set of data frames.

this_transcript_domain_starts <- exon_num_per_pfam_domain_start[which(exon_num_per_pfam_domain_start$Transcript == transcript),]
this_transcript_domain_ends <- exon_num_per_pfam_domain_end[which(exon_num_per_pfam_domain_end$Transcript == transcript),]

#Look at these in more detail.

#data.frame(this_transcript_domain_starts[,c("Pfam.ID","Pfam.short.label")],
#ExonInOrBefore.start = this_transcript_domain_starts$ExonInOrBefore,
#Percent.into.exon.start = round(this_transcript_domain_starts$PercentIntoExon,digits=2),
#ExonInOrBefore.end = this_transcript_domain_ends$ExonInOrBefore,
#PercentIntoExon.end = round(this_transcript_domain_ends$PercentIntoExon,digits=2))

#    Pfam.ID Pfam.short.label ExonInOrBefore.start Percent.into.exon.start ExonInOrBefore.end PercentIntoExon.end
#128 PF07686                A                    3                   12.73                  4               89.94
#126 PF05790                B                    5                    1.28                  5               99.57
#129 PF09191                C                    6                    1.72                  6               97.99
#127 PF05790                B                    6                   98.28                  8                6.56
#130 PF12104                D                    8                   98.36                 10                0.86

#Currently, domains are labeled with letters (eg A-D if there are 4 unique domains).
#Map these to which color to make them.

#color_per_rectangle <- plyr::mapvalues(x = this_transcript_domain_starts$Pfam.short.label,from = unique(this_transcript_domain_starts$Pfam.short.label)[order(unique(this_transcript_domain_starts$Pfam.short.label))],to = mycol[1:(length(unique(this_transcript_domain_starts$Pfam.short.label)))])

color_per_rectangle <- plyr::mapvalues(x = this_transcript_domain_starts$Pfam.short.label,
from = LETTERS[1:9],
to = mycol)

#We draw the boxes with y-coordinates from 7.75 to 8.
#x-coordinates based on the exon of the domain start and end, converted to plot coordinates using x_coord_on_fusion_plot.

for(i in 1:nrow(this_transcript_domain_starts))
{
rect(x_coord_on_fusion_plot(this_transcript_domain_starts[i,]),7.75,x_coord_on_fusion_plot(this_transcript_domain_ends[i,]),8,col=color_per_rectangle[i])
}

#Now we just need to add a legend saying what each color represents.
#Make a table with description per Pfam ID by taking only the first occurence per Pfam ID.
#Also make sure we order by A/B/C label.

this_transcript_description_per_domain <- data.frame(Pfam.short.label = this_transcript_domain_starts$Pfam.short.label,
Pfam.ID = this_transcript_domain_starts$Pfam.ID,
Pfam.description = this_transcript_domain_starts$Pfam.description,
stringsAsFactors=FALSE)

this_transcript_description_per_domain <- this_transcript_description_per_domain[which(duplicated(this_transcript_description_per_domain$Pfam.ID) == FALSE),]
this_transcript_description_per_domain <- this_transcript_description_per_domain[order(this_transcript_description_per_domain$Pfam.short.label),]

legend("center",legend=transcript,col="white",bty="n")

legend("bottom",
paste0(this_transcript_description_per_domain$Pfam.ID,"/",this_transcript_description_per_domain$Pfam.description),
col=mycol[1:(length(unique(this_transcript_domain_starts$Pfam.short.label)))],
lwd=3,bty="n")

#I also tried testing with the MGA-NUTM1 transcripts (ENST00000219905 and ENST00000537011) with similarly good results.



