#This file contains functions to actually make the plots after creating all the necessary intermediate objects like fusion_file_main_info_match_up_geneA and fusion_file_main_info_match_up_geneB.

#First, here is the function to make the initial plot, incl. draw the boxes and make the plot the appropriate width.
#Only thing is I no longer leave room for the coverage bars, so the y coordinates are different than they were in earlier versions of the script.

setup_fusion_plot <- function(exon_num){

	#Plot will go from y=1 to y=10, and x=0 to x=120*(exon_num - 1)+200.
	#The plot x-limit is based on the fact that the exons will be 101 bp, and the introns will be 20bp, plus a little extra room to put text, legends, etc.
	#Then we can use the standard coordinate system to place things. For example, if we want something near the top but not quite, and in the middle, setting the x coordinate as half the plot max,
	#and the y coordinate at 9.5 out of 10, should work well.

	plot_xmax <- 120*(exon_num - 1)+200

	plot(0,pch='',ylab='',xlab='',xlim=c(0,plot_xmax),ylim=c(0,10),xaxt = 'n',yaxt = 'n',bty='n')

	#Go through each exon and draw a rectangle at the appropriate x coordinate, as well as an intron line at the end of the rectangle if we are not up to the last exon yet.

	for(exon in 0:(exon_num - 1))
	{
		rect(120*exon,8.5,(120*exon)+100,9.25)
		end_this_exon <- (120*exon)+100
		if(exon < (exon_num - 1)){lines(c(end_this_exon,end_this_exon+20),c(8.875,8.875),type="l")}
	}
}

#Added 02/02/18: Instead of having a function called draw_breakpoint, we are going to have a function called x_coord_on_fusion_plot.
#This is almost the same, but stops at getting breakpoint_x_coordinate_on_plot.
#Instead of plotting this as part of the function, we return it, then can do the actual drawing of the line outside of the function.

#Outside of the function, we want to note that we will add 0 to the x coordinate of the exon if breakpoint is at exact start, 100 if at exact end, 110 is breakpoint is in the intron after the exon.

amount_to_add_to_exon_coordinate <- data.frame(StartEndIntronic = c("Start","End","Intronic"),Amount.to.add = c(0,100,110),stringsAsFactors=FALSE,row.names=c("Start","End","Intronic"))

x_coord_on_fusion_plot <- function(breakpoint_info_data_frame){
	
	exon_index <- breakpoint_info_data_frame$ExonInOrBefore - 1
	breakpoint_type <- breakpoint_info_data_frame$StartMiddleEnd

	exon_of_interest_start_coordinate <- 120*exon_index
	
	#If breakpoint is at the exact start or end of an exon, or in an intron, we can get the x coordinate without using PercentIntoExon from breakpoint_info_data_frame.
	#Instead, we add 0 if start coordinate, 100 if end, or 110 if in intron after, as specified in amount_to_add_to_exon_coordinate data frame.

	#If breakpoint is in the middle of an exon, we add the appropriate amount based on PercentIntoExon.

	breakpoint_x_coordinate_on_plot <- ifelse(breakpoint_type == "Middle",
		exon_of_interest_start_coordinate + (breakpoint_info_data_frame$PercentIntoExon/100)*101,
		exon_of_interest_start_coordinate + amount_to_add_to_exon_coordinate[breakpoint_type,"Amount.to.add"])

	return(breakpoint_x_coordinate_on_plot)
}

#Also added 02/02/18: function to plot the protein domains for a given transcript.
#After drawing the actual rectangles with color, function also returns domain descriptions.
#We do not add these descriptions to the plot as part of the function.
#Instead after run this function, save the return value so can add the domain descriptions legend later on.

domaincolors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000") #Before run function, set up a colorblind-friendly vector called domaincolors that we can use to give each unique domain a color.

plot_domains_of_function <- function(transcript){
	
	this_transcript_domain_starts <- exon_num_per_pfam_domain_start[which(exon_num_per_pfam_domain_start$Transcript == transcript),]
	this_transcript_domain_ends <- exon_num_per_pfam_domain_end[which(exon_num_per_pfam_domain_end$Transcript == transcript),]

	color_per_rectangle <- plyr::mapvalues(x = this_transcript_domain_starts$Pfam.short.label,
		from = LETTERS[1:9],
		to = domaincolors,
		warn_missing=FALSE)

	#We draw the boxes with y-coordinates from 7.75 to 8.
	#x-coordinates based on the exon of the domain start and end, converted to plot coordinates using x_coord_on_fusion_plot.

	for(i in 1:nrow(this_transcript_domain_starts))
	{
		rect(x_coord_on_fusion_plot(this_transcript_domain_starts[i,]),7.75,x_coord_on_fusion_plot(this_transcript_domain_ends[i,]),8,col=color_per_rectangle[i])
	}

	#Now for our return value, the information to put in the domain legend.

	this_transcript_description_per_domain <- data.frame(Pfam.short.label = this_transcript_domain_starts$Pfam.short.label,
		Pfam.ID = this_transcript_domain_starts$Pfam.ID,
		Pfam.description = this_transcript_domain_starts$Pfam.description,
		stringsAsFactors=FALSE)

	this_transcript_description_per_domain <- this_transcript_description_per_domain[which(duplicated(this_transcript_description_per_domain$Pfam.ID) == FALSE),] #Take only the first occurence of each unique domain.
	this_transcript_description_per_domain <- this_transcript_description_per_domain[order(this_transcript_description_per_domain$Pfam.short.label),] #Order domain descriptions in same order as they appear in plot.

	#For the legend, will want the Pfam ID and description, along with the appropriate colors.
	#Return this as a data frame.

	info_for_legend <- data.frame(Text = paste0(this_transcript_description_per_domain$Pfam.ID,"/",
		this_transcript_description_per_domain$Pfam.description),
		Color = domaincolors[1:(length(unique(this_transcript_domain_starts$Pfam.short.label)))],
		stringsAsFactors=FALSE)

	return(info_for_legend)
}

#Finally, let's write a function to generate the plot including boxes and breakpoints for each gene plus appropriate text for a given row number in our "main info" objects.

final_plot_for_given_row_num <- function(i){
	
	#Make a two-panel plot, with appropriate margins.

	par(mfrow=c(1,2))
	par(mai = c(1, 0.5, 0.1, 0.1))

	#Get transcript ID for geneA and geneB.

	geneA_transcript  <- fusion_file_main_info_match_up_geneA$Transcript[i]
	geneB_transcript <- fusion_file_main_info_match_up_geneB$Transcript[i]

	#Set the color of the breakpoint and UTR lines.
	#Eventually, will set the color of the breakpoint line dynamically depending on in-frame vs. not in-frame.
	#And will not have UTR lines (will plot boxes smaller for the UTRs instead).

	breakpoint_color <- "orange"
	UTR_color <- "blue"

	#Draw the boxes for geneA.

	setup_fusion_plot(fusion_file_main_info_match_up_geneA$TotalExonNum[i])

	#Draw the breakpoint for geneA.

	geneA_x_coordinate_on_fusion_plot <- x_coord_on_fusion_plot(fusion_file_main_info_match_up_geneA[i,c("ExonInOrBefore","StartMiddleEnd","PercentIntoExon")])
	lines(c(geneA_x_coordinate_on_fusion_plot,geneA_x_coordinate_on_fusion_plot),c(8.35,9.4),col=breakpoint_color)

	#If geneA has a 5' and/or 3' UTR, draw these as well.

	if(geneA_transcript %in% rownames(fiveprime_utrs)){
		fiveprime_utr_coordinate_on_fusion_plot <- x_coord_on_fusion_plot(fiveprime_utrs[geneA_transcript,c("ExonInOrBefore","StartMiddleEnd","PercentIntoExon")])
		lines(c(fiveprime_utr_coordinate_on_fusion_plot,fiveprime_utr_coordinate_on_fusion_plot),c(8.35,9.4),col=UTR_color)
	}

	if(geneA_transcript %in% rownames(threeprime_utrs)){
		threeprime_utr_coordinate_on_fusion_plot <- x_coord_on_fusion_plot(threeprime_utrs[geneA_transcript,c("ExonInOrBefore","StartMiddleEnd","PercentIntoExon")])
		lines(c(threeprime_utr_coordinate_on_fusion_plot,threeprime_utr_coordinate_on_fusion_plot),c(8.35,9.4),col=UTR_color)
	}

	#Draw the rectangles for the protein domains, and save descriptions of the domains in an object.

	if(geneA_transcript %in% exon_num_per_pfam_domain_start$Transcript){
		geneA_domain_descriptions <- plot_domains_of_function(geneA_transcript)
	}

	#Now, we are still on the left side of the panel plot, which is where we want to put all our explanatory text like:
	#Gene names
	#Both geneA and geneB transcripts
	#etc.

	#We can get the x coordinate of the middle of the plot by taking the median of par('usr')[1] and par('usr')[2].
	#Middle of the plot on y-coordinate will always be y=5.
	
	gene_names <- fusion_file_main_info$GeneA.GeneB[i]
	fusion_description <- fusion_file_main_info$FusionDescription[i]
	other_breakpoint_specific_info <- gsub('\\;','\n',fusion_file_main_info_match_up_geneA$Other.breakpoint.specific.info[i])
	geneA_transcript_ID_plus_breakpoint_text <- paste0("GeneA.transcript=",geneA_transcript," (GeneA.Breakpoint=",fusion_file_main_info$GeneA.Breakpoint[i],")")
	geneB_transcript_ID_plus_breakpoint_text <- paste0("GeneB.transcript=",geneB_transcript," (GeneB.Breakpoint=",fusion_file_main_info$GeneB.Breakpoint[i],")")

	all_descriptive_text <- paste(c(gene_names,fusion_description,other_breakpoint_specific_info,"\n\n",geneA_transcript_ID_plus_breakpoint_text,geneB_transcript_ID_plus_breakpoint_text),collapse="\n")

	text(median(c(par('usr')[1],par('usr')[2])),5,all_descriptive_text,cex=0.75)

	#Now we just need to plot for geneB, then add descriptions of protein domains to the geneB plot so they will be on the right side.

	setup_fusion_plot(fusion_file_main_info_match_up_geneB$TotalExonNum[i])
	
	geneB_x_coordinate_on_fusion_plot <- x_coord_on_fusion_plot(fusion_file_main_info_match_up_geneB[i,c("ExonInOrBefore","StartMiddleEnd","PercentIntoExon")])
	lines(c(geneB_x_coordinate_on_fusion_plot,geneB_x_coordinate_on_fusion_plot),c(8.35,9.4),col=breakpoint_color)

	if(geneB_transcript %in% rownames(fiveprime_utrs)){
		fiveprime_utr_coordinate_on_fusion_plot <- x_coord_on_fusion_plot(fiveprime_utrs[geneB_transcript,c("ExonInOrBefore","StartMiddleEnd","PercentIntoExon")])
		lines(c(fiveprime_utr_coordinate_on_fusion_plot,fiveprime_utr_coordinate_on_fusion_plot),c(8.35,9.4),col=UTR_color)
	}

	if(geneB_transcript %in% rownames(threeprime_utrs)){
		threeprime_utr_coordinate_on_fusion_plot <- x_coord_on_fusion_plot(threeprime_utrs[geneB_transcript,c("ExonInOrBefore","StartMiddleEnd","PercentIntoExon")])
		lines(c(threeprime_utr_coordinate_on_fusion_plot,threeprime_utr_coordinate_on_fusion_plot),c(8.35,9.4),col=UTR_color)
	}

	if(geneB_transcript %in% exon_num_per_pfam_domain_start$Transcript){
		geneB_domain_descriptions <- plot_domains_of_function(geneB_transcript)
	}

	#Add descriptions of protein domains as legends.

	if(geneA_transcript %in% exon_num_per_pfam_domain_start$Transcript){
		legend(par('usr')[1],7,legend=geneA_domain_descriptions$Text,col=geneA_domain_descriptions$Color,lwd=3,title="GeneA domains")
	}

	if(geneB_transcript %in% exon_num_per_pfam_domain_start$Transcript){
		legend(par('usr')[1],3.5,legend=geneB_domain_descriptions$Text,col=geneB_domain_descriptions$Color,lwd=3,title="GeneB domains")
	}

}


