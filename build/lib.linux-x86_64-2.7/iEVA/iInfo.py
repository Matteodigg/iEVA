import re
import pysam
import math

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

	#Reads_Info['Total_Mapping_Quality'] += read.alignment.mapping_quality

	if read.alignment.is_duplicate:
		Reads_Info['Duplicate_reads'] += 1

	else:
		Reads_Info['Total_Mapping_Quality'] += read.alignment.mapping_quality
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

		if read.alignment.mate_is_unmapped:
			Reads_Info['Mate_Unmapped'] += 1

	return Reads_Info


def Genotype_Reads_Info(read, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match, STOP):

	if Allele == 'ALT':

		if Variant_Class == 'snv':
			Reads_Info['Alt_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['MappingRankAlt'] += [read.alignment.mapping_quality]
			Reads_Info['Base_Alt_Qual'] += read.alignment.query_qualities[read.query_position]
			Reads_Info['QualRankAlt'] += [read.alignment.query_qualities[read.query_position]]
			Reads_Info['PosRankAlt'] += [(read.query_position)+1]
			Around_Qual = read.alignment.query_qualities[read.query_position-3:read.query_position] + read.alignment.query_qualities[read.query_position+1:read.query_position+4]
			Reads_Info['BQAltAround'] += sum(Around_Qual)/len(Around_Qual)
			Reads_Info['AltNormPos'] += math.fabs(round(float(((read.query_position)+1)-(read.alignment.query_alignment_length/2))/float(read.alignment.query_alignment_length/2),3))

		elif Variant_Class == 'Del':
			Reads_Info['Alt_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['MappingRankAlt'] += [read.alignment.mapping_quality]
			Alt_Qual = round(float(sum([read.alignment.query_qualities[read.query_position]] + [read.alignment.query_qualities[read.query_position+1]]))/float(2),2)
			Reads_Info['Base_Alt_Qual'] += Alt_Qual
			Reads_Info['QualRankAlt'] += [Alt_Qual]
			Reads_Info['PosRankAlt'] += [(read.query_position)+1]
			Around_Qual = read.alignment.query_qualities[read.query_position-2:read.query_position+1] + read.alignment.query_qualities[read.query_position+1:read.query_position+4]
			Reads_Info['BQAltAround'] += sum(Around_Qual)/len(Around_Qual)
			Reads_Info['AltNormPos'] += math.fabs(round(float(((read.query_position)+1)-(read.alignment.query_alignment_length/2))/float(read.alignment.query_alignment_length/2),3))

		elif Variant_Class == 'Ins':
			Reads_Info['Alt_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['MappingRankAlt'] += [read.alignment.mapping_quality]
			Alt_Qual = round(float(sum(read.alignment.query_qualities[Exact_Match[1][0]:Exact_Match[-1][0]]))/float(variant_lenght-1),2)
			Reads_Info['Base_Alt_Qual'] += Alt_Qual
			Reads_Info['QualRankAlt'] += [Alt_Qual]
			Reads_Info['PosRankAlt'] += [(read.query_position)+1]
			Around_Qual = read.alignment.query_qualities[(Exact_Match[1][0])-3:(Exact_Match[1][0])] + read.alignment.query_qualities[(Exact_Match[1][0])+variant_lenght-1:(Exact_Match[1][0])+variant_lenght-1+3]
			Reads_Info['BQAltAround'] += sum(Around_Qual)/len(Around_Qual)
			Reads_Info['AltNormPos'] += math.fabs(round(float(((read.query_position)+1)-(read.alignment.query_alignment_length/2))/float(read.alignment.query_alignment_length/2),3))

		Reads_Info['Read_Alt'] += 1
		Reads_Info['AS_Alt'] += read.alignment.get_tag('AS')
		Reads_Info['XS_Alt'] += read.alignment.get_tag('XS')

		if read.alignment.get_tag('XS') == 0:
			Reads_Info['XS0_Alt'] += 1
		
		if read.alignment.mapping_quality == 0:
			Reads_Info['MapQ0_Alt'] += 1

		if 'S' in read.alignment.cigarstring or 'H' in read.alignment.cigarstring:
			Reads_Info['Clipped_Reads_Alt'] += 1
			text = re.findall(r'[0-9]+[S|H]', read.alignment.cigarstring)
			Reads_Info['ClipRankAlt'] += [sum(int(s.replace('H','').replace('S','')) for s in text)]

		if read.alignment.mate_is_reverse:
			Reads_Info['ALT+'] += 1

		if read.alignment.is_reverse:
			Reads_Info['ALT-'] += 1

	elif Allele == 'REF':

		if Variant_Class == 'snv':
			Reads_Info['Ref_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['MappingRankRef'] += [read.alignment.mapping_quality]
			Reads_Info['Base_Ref_Qual'] += read.alignment.query_qualities[read.query_position]
			Reads_Info['QualRankRef'] += [read.alignment.query_qualities[read.query_position]]
			Reads_Info['PosRankRef'] += [(read.query_position)+1]

		elif Variant_Class == 'Del':
			Reads_Info['Ref_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['MappingRankRef'] = [read.alignment.mapping_quality]
			Base_Qual = round(float(sum(read.alignment.query_qualities[read.query_position:read.query_position+variant_lenght]))/float(variant_lenght),2)
			Reads_Info['Base_Ref_Qual'] += Base_Qual
			Reads_Info['QualRankRef'] = [Base_Qual]
			Reads_Info['PosRankRef'] += [(read.query_position)+1]

		elif Variant_Class == 'Ins':
			Reads_Info['Ref_Mapping_Quality'] += read.alignment.mapping_quality
			Reads_Info['MappingRankRef'] = [read.alignment.mapping_quality]
			Base_Qual = round(float(sum(read.alignment.query_qualities[Exact_Match[0][0]:(Exact_Match[-1][0])+1]))/float(2),2)
			Reads_Info['Base_Ref_Qual'] += Base_Qual
			Reads_Info['QualRankRef'] = [Base_Qual]
			Reads_Info['PosRankRef'] += [(read.query_position)+1]

		Reads_Info['Read_Ref'] += 1
		Reads_Info['AS_Ref'] += read.alignment.get_tag('AS')
		Reads_Info['XS_Ref'] += read.alignment.get_tag('XS')

		if read.alignment.get_tag('XS') == 0:
			Reads_Info['XS0_Ref'] += 1
		
		if read.alignment.mapping_quality == 0:
			Reads_Info['MapQ0_Ref'] += 1
				
		if 'S' in read.alignment.cigarstring or 'H' in read.alignment.cigarstring:
			Reads_Info['Clipped_Reads_Ref'] += 1
			text = re.findall(r'[0-9]+[S|H]', read.alignment.cigarstring)
			Reads_Info['ClipRankRef'] += [sum(int(s.replace('H','').replace('S','')) for s in text)]

		if read.alignment.mate_is_reverse:
			Reads_Info['REF+'] += 1

		if read.alignment.is_reverse:
			Reads_Info['REF-'] += 1

	return Reads_Info