#!/bin/bash

ID=$1
TELOMERE_REFERENCE=$2
CUSTOM_GENOME_REFERENCE=$3

# Map to telomere sequences
bwa mem \
	-t24 \
	-M \
	-R "@RG\tID:$ID\tPL:illumina\tPU:0\tLB:$ID\tSM:$ID" \
	$TELOMERE_REFERENCE \
	trimmed.r1.fq.gz \
	trimmed.r2.fq.gz \
	> TMP_telomere.sam


# Create exclude list
samtools view -S -f67 TMP_telomere.sam | cut -f1 > exclude.list

# Now map to full reference
bwa mem \
	-t24 \
	-M \
	-R "@RG\tID:$ID\tPL:illumina\tPU:0\tLB:$ID\tSM:$ID" \
	$CUSTOM_GENOME_REFERENCE \
	trimmed.r1.fq.gz \
	trimmed.r2.fq.gz \
	> TMP.sam

# Convert to BAM.
samtools view -bt $CUSTOM_GENOME_REFERENCE TMP.sam | samtools sort - $ID

# Identify bad-pairs that straddle the two references
java -jar WGP-Toolkit/BadPairsOnly.jar \
	-b $ID.bam \
	-o $ID\_badPairsOnly.bam

# Remove ambiguous secondary alignments.
samtools view -bF \
	256 $ID\_badPairsOnly.bam \
	> $ID\_badPairsOnly_unique.bam

samtools index $ID\_badPairsOnly_unique.bam

# Filter full genome mapping with telomere exclude list
java -jar WGP-Toolkit/SubsampleBAMWithReadNames.jar \
        -b $ID\_badPairsOnly_unique.bam \
        -e true \
        -l exclude.list \
        -o $ID\_badPairsOnly_unique_filtered.bam

samtools index $ID\_badPairsOnly_unique_filtered.bam

# Create Circos plot input for all bad pairs.
create_linkage_file_v3.pl \
	-i DB018_badPairsOnly_unique_filtered.bam \
	> all.links

# Convert to Circos one-line link format for downstream compatibility.
to_one_line_format.pl -i all.links > all.oneline.links

# Filter for each telomere reference.
filter_oneline_links_by_ref.pl -i all.oneline.links -r 17p > 17p.links
filter_oneline_links_by_ref.pl -i all.oneline.links -r 21q > 21q.links
filter_oneline_links_by_ref.pl -i all.oneline.links -r XpYp > XpYp.links

# Apply the Circos bundling script, bundlelinks. 
bundlelinks \	
	-links 17p.links \
	-max_gap 100 \
	-min_bundle_membership 5 \
	> 17p.bundled.100.5.links

bundlelinks \
	-links 21q.links \
	-max_gap 100 \
	-min_bundle_membership 5 \
	> 21q.bundled.100.5.links

bundlelinks \
	-links XpYp.links \
	-max_gap 100 \
	-min_bundle_membership 5 \
	> XpYp.bundled.100.5.links

At this point links files are ready to view in Circos.

# Subsampling BAM for viewing
# e.g. 17p
# (same procedure applied to reads mapping to XpYp and 21q).
samtools view \
	$ID\_badPairsOnly_unique_filtered.bam | \
	perl -ne '@x=split(/\t/,$_); print if $x[2] eq "17p" || $x[6] eq "17p"' | \
	samtools view -bt $CUSTOM_GENOME_REFERENCE - | \
	samtools sort - 17p_only
 
samtools index 17p_only.bam

# Creating BED file to annotate mapping hot-spots
samtools mpileup 17p_only.bam -ABd100000 | \
	pileupToBed.pl -i - -m 3 -g 10 -r true > 17p_only.m3_g10.bed

# Creating spreadsheet of mapping hot-spots.
# e.g. 17p
# (same procedure applied to reads mapping to XpYp and 21q).
annotate_hotspots_basic.pl \
	-f 17p_only.m3_g10.bed \
	-b 17p_only.bam \
	-r  $CUSTOM_GENOME_REFERENCE |\
 	annotate_hotspots__add_feature_distance.pl \ 
	-i - \
	-f 21q_only.m3_g10.bed \
	-l 'closest 21q hotspot' |\
	annotate_hotspots__add_feature_distance.pl \
	-i - \
	-f XpYp_only.m3_g10.bed \
	-l 'closest XpYp hotspot' |\
     annotate_hotspots__add_feature_distance.pl \
	-i - \
	-f merged_telomere.bed \
	-l 'closest mask' |\
      annotate_hotspots__add_telomere_distance.pl \
	-i - \
	-f centromere_telomere_locations.bed \
	-l "closest telomere" |\
	annotate_hotspots__add_consensus.pl \
	-i - \
	-b 17p_only.bam |\
	textToExcel.pl \
	-i - \
	-o 17p_hotspot_table.xlsx

