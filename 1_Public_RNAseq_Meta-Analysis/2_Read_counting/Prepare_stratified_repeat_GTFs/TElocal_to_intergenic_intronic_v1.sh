####################################################### INPUTS
# Define the output directory
output_dir='/Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE_subsets_vs_GENCODE44/'

# Define gene and TE GTF files
# NOTE: THE TE BED FILE WAS MADE WITH A CUSTOM SCRIPT 'TElocal_locations_to_BED.R'
my_gene_gtf='/Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/homo_sapiens_gencode_v44/gencode.v44.primary_assembly.annotation.gtf'
my_TE_BED='/Users/juanb/Library/CloudStorage/Dropbox/Research/Shared_Bioinformatic_Resources/STAR_genome_indices/TElocal_rmsk_files/GRCh38_GENCODE_rmsk_TE.GTF'

####################################################### PROCESS INPUTS

# Go to the output_dir
cd $output_dir

# Split gene GTF
grep '	transcript	' $my_gene_gtf > GRCh38.genes.gtf
grep 'exon' $my_gene_gtf > GRCh38.exons.gtf
#grep '	CDS	' $my_gene_gtf > GRCh38.CDS.gtf

# Extract intergenic TEs
# The subtract -A flag removes a feature entirely if there is any overlap
bedtools subtract -nonamecheck -A -a $my_TE_BED -b GRCh38.genes.gtf > Repeats_intergenic.gtf

    # Add 5 kb to both sides of each intergenic TE and output repeats that overlap a gene (-u flag)
    bedtools window -w 5000 -u -a Repeats_intergenic.gtf -b GRCh38.genes.gtf > Repeats_intergenic_within_5kb.gtf

    # Add 5 kb to both sides of each intergenic TE and output repeats that DO NOT overlap a gene (-v flag)
    bedtools window -w 5000 -v -a Repeats_intergenic.gtf -b GRCh38.genes.gtf > Repeats_intergenic_farther_than_5kb.gtf

# Extract exonic TEs
# The -u flag reports an original entry in A once if any overlaps are found in B.
bedtools intersect -nonamecheck -u -a $my_TE_BED -b GRCh38.exons.gtf > Repeats_exonic.gtf

# Remove intergenic and exonic repeats from the repeat list
# The subtract -A flag removes a feature entirely if there is any overlap
bedtools subtract -nonamecheck -A -a $my_TE_BED -b Repeats_intergenic.gtf > Repeats_minus_intergenic.gtf
bedtools subtract -nonamecheck -A -a Repeats_minus_intergenic.gtf -b Repeats_exonic.gtf > Repeats_minus_intergenic_exonic.gtf

# Concatenate the exon gtf with exonic repeats and sort.
#The 'f' in the second command ignores upper/lowercases in chromosome names.
cat GRCh38.exons.gtf Repeats_exonic.gtf > GRCh38.exons_with_exonic_repeats.gtf
sort -k1,1f -k4,4n GRCh38.exons_with_exonic_repeats.gtf > GRCh38.exons_with_exonic_repeats.sorted.gtf

# Generate merged exon+repeat BED file
bedtools merge -i GRCh38.exons_with_exonic_repeats.sorted.gtf -c 2,9 -o distinct > GRCh38.exons_with_exonic_repeats.merged.bed

# Subtract merged exons+exonic_repeat features from the gene GTF, to get introns
bedtools subtract -nonamecheck -a GRCh38.genes.gtf -b GRCh38.exons_with_exonic_repeats.merged.bed > GRCh38.introns.gtf

# Extract intronic TEs
# The -u flag reports an original entry in A once if any overlaps are found in B.
bedtools intersect -nonamecheck -u -a Repeats_minus_intergenic_exonic.gtf -b GRCh38.introns.gtf > Repeats_intronic.gtf

# Remove intronic repeats from the repeat list
bedtools subtract -nonamecheck -A -a Repeats_minus_intergenic_exonic.gtf -b Repeats_intronic.gtf > Repeats_minus_intergenic_exonic_intronic.gtf











# Delete unneeded files
rm GRCh38.exons.gtf
rm GRCh38.genes.gtf
rm Repeats_minus_intergenic.gtf
rm Repeats_minus_intergenic_exonic.gtf
rm GRCh38.exons_with_exonic_repeats.gtf
rm GRCh38.exons_with_exonic_repeats.sorted.gtf
rm GRCh38.exons_with_exonic_repeats.merged.bed
rm GRCh38.introns.gtf
