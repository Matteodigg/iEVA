import re
import pysam

def Extract_Variant_Type(variant, H_CHR):

	CHR = variant[H_CHR.index('#CHROM')]
	POS = variant[H_CHR.index('POS')]
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]

	Variant_Class = '.'

	if re.match("[A|T|C|G]+$",REF) and re.match("[A|T|C|G]+$",ALT):
		if len(REF) > len(ALT) and len(ALT) == 1:
			Variant_Class = 'Del'

		elif len(ALT) > len(REF)  and len(REF) == 1:
			Variant_Class = 'Ins'

		elif len(ALT) == len(REF) and len(REF) == 1 and len(ALT) == 1:
			Variant_Class = 'snv'

		elif len(ALT) == len(REF) or len(REF) != len(ALT) and len(REF) > 1 and len(ALT) > 1:
			Variant_Class = 'Alt'

		else:
			Variant_Class = '.'

	return Variant_Class
	

def Extract_Reds_Info(read, Reads_Info):

	Reads_Info['Total_Reads_Unfilter'] += 1

	Reads_Info['Total_Mapping_Quality'] += read.alignment.mapping_quality

	if read.alignment.is_duplicate:
		Reads_Info['Duplicate_reads'] += 1

	if read.alignment.is_duplicate == False:
		Reads_Info['Total_Reads_No_Dup'] += 1
		Reads_Info['Alignment_Score'] += read.alignment.get_tag('AS')
		Reads_Info['Suboptimal_Alignment_Score'] += read.alignment.get_tag('XS')

		if read.alignment.mapping_quality == 0:
			Reads_Info['Mapping_Quality_Zero'] += 1

		if read.alignment.is_unmapped:
			Reads_Info['Unmapped_reads'] += 1

		if read.alignment.is_paired == False:
			Reads_Info['Not_Paired_Reads'] += 1

		if read.alignment.is_proper_pair == False:
			Reads_Info['Not_Proper_Paired_Reads'] += 1

		if read.alignment.is_supplementary:
			Reads_Info['Supplementary_Align'] += 1

		if read.alignment.is_secondary:
			Reads_Info['Not_Primary_Align'] += 1

		if read.alignment.is_unmapped == False:

			if read.alignment.is_read1:
				Reads_Info['is_read1'] += 1
				if read.alignment.mate_is_reverse:
					Reads_Info['is_read1_forward'] += 1
				elif read.alignment.is_reverse:
					Reads_Info['is_read1_reverse'] += 1
				
			elif read.alignment.is_read2:
				Reads_Info['is_read2'] += 1
				if read.alignment.is_reverse:
					Reads_Info['is_read2_reverse'] += 1
				if read.alignment.mate_is_reverse:
					Reads_Info['is_read2_forward'] += 1

	return Reads_Info


def Genotype_Reads_Info(read, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match):

	if Allele == 'ALT':

		if Variant_Class == 'snv':
			Reads_Info['Alt_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['Base_Alt_Qual'] += read.alignment.query_qualities[read.query_position]
		
		elif Variant_Class == 'Del':
			Reads_Info['Alt_Mapping_Quality'] += read.alignment.mapping_quality
			Alt_Qual = round(float(sum([read.alignment.query_qualities[read.query_position]] + [read.alignment.query_qualities[read.query_position+1]]))/float(2),2)
			Reads_Info['Base_Alt_Qual'] += Alt_Qual

		elif Variant_Class == 'Ins':
			Reads_Info['Alt_Mapping_Quality'] += read.alignment.mapping_quality
			Alt_Qual = round(float(sum(read.alignment.query_qualities[Exact_Match[1][0]:Exact_Match[-1][0]]))/float(variant_lenght-1),2)
			Reads_Info['Base_Alt_Qual'] += Alt_Qual

		Reads_Info['Read_Alt'] += 1
		Reads_Info['AS_Alt'] += read.alignment.get_tag('AS')
		Reads_Info['XS_Alt'] += read.alignment.get_tag('XS')

		if read.alignment.get_tag('XS') == 0:
			Reads_Info['XS0_Alt'] += 1
		
		if read.alignment.is_unmapped:
			Reads_Info['UnMap_Alt'] += 1
		
		if read.alignment.mapping_quality == 0:
			Reads_Info['MapQ0_Alt'] += 1
		
		if read.alignment.is_secondary:
			Reads_Info['NPA_Alt'] += 1
		
		if read.alignment.is_supplementary:
			Reads_Info['Suppl_Alt'] += 1
		
		if read.alignment.is_paired == False:
			Reads_Info['NP_Alt'] += 1
		
		if read.alignment.is_proper_pair == False:
			Reads_Info['NPP_Alt'] += 1
		
		if 'S' in read.alignment.cigarstring or 'H' in read.alignment.cigarstring:
			Reads_Info['Clipped_Reads_Alt'] += 1

	elif Allele == 'REF':

		if Variant_Class == 'snv':
			Reads_Info['Ref_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['Base_Ref_Qual'] += read.alignment.query_qualities[read.query_position]

		elif Variant_Class == 'Del':
			Reads_Info['Ref_Mapping_Quality'] += read.alignment.mapping_quality
			Base_Qual = round(float(sum(read.alignment.query_qualities[read.query_position:read.query_position+variant_lenght]))/float(variant_lenght),2)
			Reads_Info['Base_Ref_Qual'] += Base_Qual

		elif Variant_Class == 'Ins':
			Reads_Info['Ref_Mapping_Quality'] += read.alignment.mapping_quality
			Base_Qual = round(float(sum(read.alignment.query_qualities[Exact_Match[0][0]:(Exact_Match[-1][0])+1]))/float(2),2)
			Reads_Info['Base_Ref_Qual'] += Base_Qual
		
		Reads_Info['Read_Ref'] += 1
		Reads_Info['AS_Ref'] += read.alignment.get_tag('AS')
		Reads_Info['XS_Ref'] += read.alignment.get_tag('XS')

		if read.alignment.get_tag('XS') == 0:
			Reads_Info['XS0_Ref'] += 1

		if read.alignment.is_unmapped:
			Reads_Info['UnMap_Ref'] += 1
		
		if read.alignment.mapping_quality == 0:
			Reads_Info['MapQ0_Ref'] += 1
		
		if read.alignment.is_secondary:
			Reads_Info['NPA_Ref'] += 1
		
		if read.alignment.is_supplementary:
			Reads_Info['Suppl_Ref'] += 1
		
		if read.alignment.is_paired == False:
			Reads_Info['NP_Ref'] += 1
		
		if read.alignment.is_proper_pair == False:
			Reads_Info['NPP_Ref'] += 1
		
		if 'S' in read.alignment.cigarstring or 'H' in read.alignment.cigarstring:
			Reads_Info['Clipped_Reads_Ref'] += 1

	return Reads_Info