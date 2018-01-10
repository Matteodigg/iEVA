import iInfo
import iSequence
import iCheck
import iBam

def iHeader(Header_File, H_CHR, opts):

	with open(opts.input) as vcf:
		for line in vcf:
			line = line.rstrip()
			if line.startswith('##'):
				Header_File += [line]
			elif line.startswith('#CHROM'):
				H_CHR += [line]
				break

	if opts.SimpleRepeat:
		Header_File += ['##INFO=<ID=iSR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">']
	if opts.SimpleRepeatLength:
		Header_File += ['##INFO=<ID=iSRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">']
	if opts.SimpleRepeatUnit:
		Header_File += ['##INFO=<ID=iSRU,Number=1,Type=String,Description="Simple repeated sequence unit composing repeated sequence (SR)">']
	if opts.PseudoNucleotidesComposition:
		Header_File += ['##INFO=<ID=iPNC,Number=16,Type=Float,Description="Pseudo Nucleotide sequence Composition using Kmer size of 2. Reported as: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT">']
	if opts.RepeatMasker:
		Header_File += ['##INFO=<ID=iRM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">']
	if opts.gcContent:
		Header_File += ['##INFO=<ID=iGC,Number=1,Type=Float,Description="Percentage of GC content in sequence">']
	if opts.VariantClass:
		Header_File += ['##INFO=<ID=iVC,Number=1,Type=String,Description="Annotated variant class: SNV=snv, Insertion=Ins, Deletion=Del, SequenceAlteration=Alt">']
	if opts.UnMappedReads:
		Header_File += ['##FORMAT=<ID=iUnMap,Number=1,Type=Float,Description="Fraction of unmapped reads">']
	if opts.MeanMappingQuality:
		Header_File += ['##FORMAT=<ID=iMQ,Number=1,Type=Float,Description="Mean mapping quality for reads mapping variant position">']
	if opts.MappingQualityZero:
		Header_File += ['##FORMAT=<ID=iMQ0,Number=1,Type=Float,Description="Fraction of reads mapping position with Mapping Quaility=0">']
	if opts.NotPrimaryAlignment:
		Header_File += ['##FORMAT=<ID=iNPA,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not primary alignment">']
	if opts.SupplementaryAlignment:
		Header_File += ['##FORMAT=<ID=iSA,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as supplementary alignment">']
	if opts.NotPairedReads:
		Header_File += ['##FORMAT=<ID=iNP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not paired">']
	if opts.NotProperPairedReads:
		Header_File += ['##FORMAT=<ID=iNPP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not proper paired">']
	if opts.AlignmentScore:
		Header_File += ['##FORMAT=<ID=iAS,Number=1,Type=Float,Description="Mean alignment score of reads mapping variant position">']
	if opts.SuboptimalAlignmentScore:
		Header_File += ['##FORMAT=<ID=iXS,Number=1,Type=Float,Description="Mean suboptimal alignment score of reads mapping variant position">']
	if opts.TotalDupReads:
		Header_File += ['##FORMAT=<ID=iDUP,Number=1,Type=Integer,Description="Fraction of total reads mapping position marked as duplicate">']
	if opts.TotalDPUnfilter:
		Header_File += ['##FORMAT=<ID=iTDP,Number=1,Type=Integer,Description="Total read depth. No filter applied">']
	if opts.iEvaDepth:
		Header_File += ['##FORMAT=<ID=iDP,Number=1,Type=Integer,Description="iEVA read depth. Duplicate reads are excluded.">']
	if opts.AlleleDepth:
		Header_File += ['##FORMAT=<ID=iAD,Number=R,Type=Integer,Description="Allelic depth reported by iEVA as Ref,Alt">']
	if opts.StrandBiasDepth:
		Header_File += ['##FORMAT=<ID=iSBD,Number=4,Type=Integer,Description="Depth of bases supporting REF and ALT allele on forward and reverse strand for strand bias detection (R+, R-, A+, A-)">']
	if opts.StrandBias:
		Header_File += ['##FORMAT=<ID=iSB,Number=1,Type=Float,Description="Fisher exact test to detect strand bias (R+,R-,A+,A-)">']
	if opts.AlleleDuplicateDepth:
		Header_File += ['##FORMAT=<ID=iADup,Number=R,Type=Integer,Description="iEVA Allele Depth of read marked as duplicate. Reported as NumberDupRef, NumberDupAlt">']
	if opts.DeltaDuplicate:
		Header_File += ['##FORMAT=<ID=iDDup,Number=1,Type=Float,Description="Difference between fraction of duplicate reads for REF and ALT alleles (DupREF-DupALT) for Het variants">']
	if opts.AlleleQscore:
		Header_File += ['##FORMAT=<ID=iQual,Number=R,Type=Float,Description="Mean Q-score for REF and ALT allele reported as MeanRefQscore, MeanAltQscore">']
	if opts.AlleleMeanMappingQuality:
		Header_File += ['##FORMAT=<ID=iAMQ,Number=R,Type=Float,Description="Mean mapping quality for reads supporting REF and ALT allele reported as MeanMappingQualRef,MeanMappingQualAlt">']
	if opts.AlleleMeanAlignmentScore:
		Header_File += ['##FORMAT=<ID=iAAS,Number=R,Type=Float,Description="Mean alignment score for reads supporting REF and ALT allele. Reported as MeanAlignmentScoreRef,MeanAlignmentScoreAlt">']
	if opts.AlleleSuboptimalAlignmentScore:
		Header_File += ['##FORMAT=<ID=iAXS,Number=R,Type=Float,Description="Mean suboptimal alignment score for reads supporting REF and ALT allele. Reported as MeanSuboptimalAlignmentScoreRef,MeanSuboptimalAlignmentScoreAlt">']
	if opts.AlleleSuboptimalAlignmentScoreZero:
		Header_File += ['##FORMAT=<ID=iAXS0,Number=R,Type=Integer,Description="Number of reads supporting REF and ALT allele with Suboptimal Alignment Score = 0. Reported as NumberReadsXSscore0Ref,NumberReadsXSscore0Alt">']
	if opts.AlleleUnMappedReads:
		Header_File += ['##FORMAT=<ID=iAUnMap,Number=R,Type=Integer,Description="Number reads supporting REF and ALT allele flagged as unmapped. Reported as NumberUnmappedRef,NumberUnmappedAlt">']
	if opts.AlleleMappingQualityZero:
		Header_File += ['##FORMAT=<ID=iAMQ0,Number=R,Type=Integer,Description="Number of reads supporting REF and ALT allele with mapping quality = 0. Reported as NumberMappingQuality0Ref,NumberMappingQuality0Alt">']
	if opts.AlleleNotPrimaryAlignment:
		Header_File += ['##FORMAT=<ID=iANPA,Number=R,Type=Integer,Description="Number of reads supporting REF and ALT allele flagged as secondary alignment. Reported as NumberSecondaryAlignmentRef,NumberSecondaryAlignmentAlt">']
	if opts.AlleleSupplementaryAlignment:
		Header_File += ['##FORMAT=<ID=iASA,Number=R,Type=Integer,Description="Number reads supporting REF and ALT allele flagged as supplementary alignment. Reported as NumberUnmappedRef,NumberUnmappedAlt">']
	if opts.AlleleNotPaired:
		Header_File += ['##FORMAT=<ID=iANP,Number=R,Type=Integer,Description="Number of reads supporting REF and ALT allele flagged as not paired. Reported as NumberNotPairedRef,NumberNotPairedAlt">']
	if opts.AlleleNotProperPaired:
		Header_File += ['##FORMAT=<ID=iANPP,Number=R,Type=Integer,Description="Number of reads supporting REF and ALT allele flagged as not proper paired. Reported as NumberNotProperPairedRef,NumberNotProperPairedAlt">']
	if opts.AlleleClippedReads:
		Header_File += ['##FORMAT=<ID=iACR,Number=R,Type=Integer,Description="Number of clipped reads mapping REF and ALT allele reported as NumberClippedRef, NumberClippedAlt">']


	return Header_File, H_CHR


def iAnnotation(Sample_dict, out, opts, Sequence_opts, Bam_opts, Genotype_opts, H_CHR, Reference, Sample_tot, Repeat_list):

	CHR_Counter = ''
	
	with open(opts.input) as vcf:

		for variant in vcf:
			if variant.startswith('#'):
				continue
			else:
				variant = variant.rstrip()
				variant = variant.split('\t')
				Variant_Class = iInfo.Extract_Variant_Type(variant, H_CHR)

				if any(Sequence_opts):
					RepeatSeq, RepeatSeq_Lenght, RepeatSeq_Unit, Psuedo_Nucleotide, RM, GC = iSequence.SimpleRepeats_Finder(Reference, variant, H_CHR, Variant_Class, opts, Repeat_list)
				if opts.SimpleRepeat:
					variant[H_CHR.index('INFO')] += ';' + 'SR=' + str(RepeatSeq)
				if opts.SimpleRepeatLength:
					variant[H_CHR.index('INFO')] += ';' + 'SRL=' + str(RepeatSeq_Lenght)
				if opts.SimpleRepeatUnit:
					variant[H_CHR.index('INFO')] += ';' + 'SRU=' + str(RepeatSeq_Unit)
				if opts.PseudoNucleotidesComposition:
					variant[H_CHR.index('INFO')] += ';' + 'PNC=' + ','.join(str(Din) for Din in Psuedo_Nucleotide)
				if opts.RepeatMasker:
					variant[H_CHR.index('INFO')] += ';' + 'RM=' + RM
				if opts.gcContent:
					variant[H_CHR.index('INFO')] += ';' + 'GC=' + GC
				if opts.VariantClass:
					variant[H_CHR.index('INFO')] += ';' + 'VC=' + str(Variant_Class)

				if any(Genotype_opts) or any(Bam_opts):
					Sample_Stat = iBam.Sequence_Annotator(variant, Sample_dict, H_CHR, Reference, Variant_Class, opts.SNVMinBaseQuality, opts.SNVMinMappingQuality, opts.IndelMinBaseQuality, opts.IndelMinMappingQuality, Bam_opts, Genotype_opts)

				if opts.UnMappedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'iUnMap'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Percentage_Unmapped_Reads'))) if sample in Sample_dict.keys() else ':.'
				if opts.MeanMappingQuality:
					variant[H_CHR.index('FORMAT')] += ':' + 'iMQ'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Total_Mean_Mapping_Quality'))) if sample in Sample_dict.keys() else ':.'
				if opts.MappingQualityZero:
					variant[H_CHR.index('FORMAT')] += ':' + 'iMQ0'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Mapping_Quality_Zero'))) if sample in Sample_dict.keys() else ':.'
				if opts.NotPrimaryAlignment:
					variant[H_CHR.index('FORMAT')] += ':' + 'iNPA'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Percentage_Not_Primary_Alignment_Reads'))) if sample in Sample_dict.keys() else ':.'
				if opts.SupplementaryAlignment:
					variant[H_CHR.index('FORMAT')] += ':' + 'iSA'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Percentage_Supplementary_Align'))) if sample in Sample_dict.keys() else ':.'
				if opts.NotPairedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'iNP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Percentage_Not_Paired_Reads'))) if sample in Sample_dict.keys() else ':.'
				if opts.NotProperPairedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'iNPP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Percentage_Not_Proper_Paired_Reads'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlignmentScore:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAS'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Mean_Alignment_Score'))) if sample in Sample_dict.keys() else ':.'
				if opts.SuboptimalAlignmentScore:
					variant[H_CHR.index('FORMAT')] += ':' + 'iXS'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Mean_Suboptimal_Alignment_Score'))) if sample in Sample_dict.keys() else ':.'
				if opts.TotalDupReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'iDUP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Total_Duplicate_Reads'))) if sample in Sample_dict.keys() else ':.'
				if opts.TotalDPUnfilter:
					variant[H_CHR.index('FORMAT')] += ':' + 'iTDP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Total_Reads_Unfilter')) if sample in Sample_dict.keys() else ':.'
				if opts.iEvaDepth:
					variant[H_CHR.index('FORMAT')] += ':' + 'iDP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Coverage')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleDepth:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAD'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Read_Ref')) + ',' + str(Sample_Stat.get(sample).get('Read_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.StrandBiasDepth:
					variant[H_CHR.index('FORMAT')] += ':' + 'iSBD'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('REF+')) + ',' + str(Sample_Stat.get(sample).get('REF-')) + ',' + str(Sample_Stat.get(sample).get('ALT+')) + ',' + str(Sample_Stat.get(sample).get('ALT-')) if sample in Sample_dict.keys() else ':.'
				if opts.StrandBias:
					variant[H_CHR.index('FORMAT')] += ':' + 'iSB'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Strand_Bias_Reads'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleDuplicateDepth:
					variant[H_CHR.index('FORMAT')] += ':' + 'iADup'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Number_Read_Dup_Ref')) + ',' + str(Sample_Stat.get(sample).get('Number_Read_Dup_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.DeltaDuplicate:
					variant[H_CHR.index('FORMAT')] += ':' + 'iDDup'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Delta_Duplicate'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleQscore:
					variant[H_CHR.index('FORMAT')] += ':' + 'iQual'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Qual_Ref')) + ',' + str(Sample_Stat.get(sample).get('Qual_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleMeanMappingQuality:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAMMQ'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Ref_Mean_Mapping_Quality'))) + ',' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('Alt_Mean_Mapping_Quality'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleMeanAlignmentScore:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAAS'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('AS_Ref'))) + ',' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('AS_Alt'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleMeanAlignmentScore:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAXS'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('XS_Ref'))) + ',' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('XS_Alt'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleMeanAlignmentScore:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAXS0'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('XS0_Ref'))) + ',' + str(iCheck.Check_Zero(Sample_Stat.get(sample).get('XS0_Alt'))) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleUnMappedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAUnMap'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('UnMap_Ref')) + ',' + str(Sample_Stat.get(sample).get('UnMap_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleMappingQualityZero:
					variant[H_CHR.index('FORMAT')] += ':' + 'iAMQ0'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('MapQ0_Ref')) + ',' + str(Sample_Stat.get(sample).get('MapQ0_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleNotPrimaryAlignment:
					variant[H_CHR.index('FORMAT')] += ':' + 'iANPA'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Not_Primary_Align_Ref')) + ',' + str(Sample_Stat.get(sample).get('Not_Primary_Align_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleSupplementaryAlignment:
					variant[H_CHR.index('FORMAT')] += ':' + 'iASA'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Suppl_Alignment_Ref')) + ',' + str(Sample_Stat.get(sample).get('Suppl_Alignment_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleNotPaired:
					variant[H_CHR.index('FORMAT')] += ':' + 'iANP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Not_Paired_Ref')) + ',' + str(Sample_Stat.get(sample).get('Not_Paired_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleNotProperPaired:
					variant[H_CHR.index('FORMAT')] += ':' + 'iANPP'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Not_Proper_Paired_Ref')) + ',' + str(Sample_Stat.get(sample).get('Not_Proper_Paired_Alt')) if sample in Sample_dict.keys() else ':.'
				if opts.AlleleClippedReads:
					variant[H_CHR.index('FORMAT')] += ':' + 'iACR'
					for sample in Sample_tot:
						variant[H_CHR.index(sample)] += ':' + str(Sample_Stat.get(sample).get('Clipped_Reads_Ref')) + ',' + str(Sample_Stat.get(sample).get('Clipped_Reads_Alt')) if sample in Sample_dict.keys() else ':.'

				if opts.verbose:
					if variant[H_CHR.index('#CHROM')] != CHR_Counter and CHR_Counter == '':
						CHR_Counter = variant[H_CHR.index('#CHROM')]
						print '\n' + 'Extracting attributes on: ' + str(CHR_Counter)
					if variant[H_CHR.index('#CHROM')] != CHR_Counter and CHR_Counter != '':
						CHR_Counter = variant[H_CHR.index('#CHROM')]
						print '\n' + 'Extracting attributes on: ' + str(CHR_Counter)

				out.write('\t'.join(variant) + '\n')