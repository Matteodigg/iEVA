import sys
import argparse
import iCheck
import vcfwriter

def main():

	parser = argparse.ArgumentParser(description='iEVA: identification and Extraction of Variant Attributes.')
	parser.add_argument('-v','--verbose', action="store_true", help="Increase output verbosity. Enable iEVA extraction progress information")
	
	Required = parser.add_argument_group(description='Required arguments')
	Required.add_argument('-I','--input',required=True,metavar='path/to/input',help="Input file in vcf format.")
	Required.add_argument('-Ref','--reference',required=True,metavar='path/to/fasta',help="Reference file in fasta format.")
	Required.add_argument('-O','--outfile',required=True,metavar='path/to/outfile',help="Output vcf file.")
	Required.add_argument('-L','--list',metavar='path/to/bamlist',help="REQUIRED ONLY TO EXTRACT FEATURES FROM BAM FILE. Bam list file, including path, one for each raw. Sample name in vcf genotype field has to be the same of Read Group Sample Name (SM) tag in bam file header. Bam Index file required in same path of bam file.")
	
	Sequence = parser.add_argument_group(description='Reference.fasta extraction arguments')
	Sequence.add_argument('-WS','--WindowSize',default=40, type=int, metavar='[20-600]',help="Set window size for -SR -SRL -GC and -PNC command. Default=40, Min=20, Max=600")
	Sequence.add_argument('-SR','--SimpleRepeat',action="store_true",help="Enable Simple Repeat Finder. Check if a variant fall in a Simple Repeated Sequence. 0 for None, 1 for Simple Repeat Sequence, 2 for Homopolymer.")
	Sequence.add_argument('-SRList','--SimpleRepeatList', metavar='path/to/SampleRepeatList.txt',help="path to list containing simple repeated sequence (one for each raw) to extract with iEVA. If enabled, only reported sequences in file will be extracted. -SR option is required.")
	Sequence.add_argument('-SRL','--SimpleRepeatLength',action="store_true",help="Length of repeated sequence (expressed as number of nucleotides) in -SR option")
	Sequence.add_argument('-SRU','--SimpleRepeatUnit',action="store_true",help="Repeated sequence unit composing simple repeated sequence")
	Sequence.add_argument('-PNC','--PseudoNucleotidesComposition',action="store_true",help="Describe nucleotide sequence using Pseudo Nucleotide Composition with Kmer size of 2. Values reported as: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT")
	Sequence.add_argument('-RM','--RepeatMasker',action="store_true",help="Check from fasta reference if a variant falls in a sequence flagged as repeated sequence by RepeatMasker tool.")
	Sequence.add_argument('-GC','--gcContent',action="store_true",help="Percentage of sequence GC content.")
	Sequence.add_argument('-VC','--VariantClass',action="store_true",help="Annotated variant class: SNV (snv), Insertion (Ins), Deletion (Del), Sequence Alteration (Alt)")	
	
	Bam = parser.add_argument_group(description='Bam extraction arguments')
	Bam.add_argument('-SBR','--StrandBiasReads',action="store_true",help="Fisher's exact test based on read orientation (R1+,R1-,R2+,R2-) to detect strand bias. Only for paired-end sequencing experiments.")
	Bam.add_argument('-UnMap','--UnMappedReads', action="store_true",help="Fraction of unmapped reads for variant position.")
	Bam.add_argument('-MMQ','--MeanMappingQuality',action="store_true",help="Mean mapping quality for reads mapping variant position")
	Bam.add_argument('-MQ0','--MappingQualityZero',action="store_true",help="Fraction of reads with Mapping Quaility=0.")
	Bam.add_argument('-NPA','--NotPrimaryAlignment',action="store_true",help="Fraction of reads mapping position flagged as not primary alignment.")
	Bam.add_argument('-SA','--SupplementaryAlignment',action="store_true",help="Fraction of reads mapping position flagged as supplementary alignment.")
	Bam.add_argument('-NP','--NotPairedReads',action="store_true",help="Fraction of reads mapping position flagged as not paired.")
	Bam.add_argument('-NPP','--NotProperPairedReads',action="store_true",help="Fraction of reads mapping position flagged as not proper paired.")
	Bam.add_argument('-AS','--AlignmentScore',action="store_true",help="Mean Alignment Score (AS tag in bam file) of reads mapping position.")
	Bam.add_argument('-XS','--SuboptimalAlignmentScore',action="store_true",help="Mean Suboptimal Alignment Score (XS tag in bam file) of reads mapping position")
	Bam.add_argument('-DUP','--TotalDupReads',action="store_true",help="Fraction of total reads mapping position marked as duplicate.")
	
	Genotype = parser.add_argument_group(description='Genotype extraction arguments')
	Genotype.add_argument('-SNVmbq','--SNVMinBaseQuality',default=12,type=int,metavar='[0-66]',help="Minimum Base Quality threshold for base supporting SNV position. Used on Genotype exatrction arguments. Default=12")
	Genotype.add_argument('-SNVmmq','--SNVMinMappingQuality',default=30,type=int,metavar='[0-60]',help="Minimum Mapping Quality threshold for reads supporting SNV position. Used on Genotype exatrction arguments. Default=30")
	Genotype.add_argument('-INDELmbq','--IndelMinBaseQuality',default=10,type=int,metavar='[0-66]',help="Minimum Base Quality threshold for base supporting InDel position. Used on Genotype exatrction arguments. Default=10")
	Genotype.add_argument('-INDELmmq','--IndelMinMappingQuality',default=20,type=int,metavar='[0-60]',help="Minimum Mapping Quality threshold for reads supporting InDel position. Used on Genotype exatrction arguments. Default=20")
	Genotype.add_argument('-TDP','--TotalDPUnfilter',action="store_true",help="Total read depth. No filter applied, include duplicate reads. Look at -iDP option for filtered iDP.")
	Genotype.add_argument('-iDP','--iEvaDepth',action="store_true",help="iEVA read depth for variant position. Only proper paired, proper mapped and not duplicate reads are included.")
	Genotype.add_argument('-iAD','--iAlleleDepth',action="store_true",help="iEVA Allele Depth reported as Ref,Alt")
	Genotype.add_argument('-iADup','--AlleleDuplicateDepth',action="store_true",help="iEVA Allele Depth of read marked as duplicate. Reported as NumberDupRef, NumberDupAlt")
	Genotype.add_argument('-iDDup','--DeltaDuplicate',action="store_true",help="Difference between fraction duplicate reads in REF and ALT alleles (REF-ALT) for Het variants. Negative values show preference in duplicate for REF allele, positive for ALT.")
	Genotype.add_argument('-iQual','--AlleleQscore',action="store_true",help="Mean Q-score for REF and ALT allele. Reported as MeanQscoreRef, MeanQscoreAlt")
	Genotype.add_argument('-iMMQ','--AlleleMeanMappingQuality',action="store_true",help="Mean mapping quality for reads supporting REF and ALT allele. Reported as MeanMappingQualRef,MeanMappingQualAlt")
	Genotype.add_argument('-iAS','--AlleleMeanAlignmentScore',action="store_true",help="Mean alignment score for reads supporting REF and ALT allele. Reported as MeanAlignmentScoreRef,MeanAlignmentScoreAlt")
	Genotype.add_argument('-iXS','--AlleleSuboptimalAlignmentScore',action="store_true",help="Mean Suboptimal Alignment Score for reads supporting REF and ALT allele. Reported as MeanXSscoreRef,MeanXSscoreAlt")
	Genotype.add_argument('-iXS0','--AlleleSuboptimalAlignmentScoreZero',action="store_true",help="Number of reads supporting REF and ALT allele with Suboptimal Alignment Score = 0. Reported as NumberReadsXSscore0Ref,NumberReadsXSscore0Alt")
	Genotype.add_argument('-iUnMap','--AlleleUnMappedReads',action="store_true",help="Number reads supporting REF and ALT allele flagged as unmapped. Reported as NumberUnmappedRef,NumberUnmappedAlt")
	Genotype.add_argument('-iMQ0','--AlleleMappingQualityZero',action="store_true",help="Number of reads supporting REF and ALT allele with mapping quality = 0. Reported as NumberMappingQuality0Ref,NumberMappingQuality0Alt")
	Genotype.add_argument('-iNPA','--AlleleNotPrimaryAlignment',action="store_true",help="Number of reads supporting REF and ALT allele flagged as secondary alignment. Reported as NumberSecondaryAlignmentRef,NumberSecondaryAlignmentAlt")
	Genotype.add_argument('-iSA','--AlleleSupplementaryAlignment',action="store_true",help="Number reads supporting REF and ALT allele flagged as supplementary alignment. Reported as NumberUnmappedRef,NumberUnmappedAlt")
	Genotype.add_argument('-iNP','--AlleleNotPaired',action="store_true",help="Number of reads supporting REF and ALT allele flagged as not paired. Reported as NumberNotPairedRef,NumberNotPairedAlt")
	Genotype.add_argument('-iNPP','--AlleleNotProperPaired',action="store_true",help="Number of reads supporting REF and ALT allele flagged as not proper paired. Reported as NumberNotProperPairedRef,NumberNotProperPairedAlt")
	Genotype.add_argument('-iCR','--AlleleClippedReads',action="store_true",help="Number of clipped reads mapping REF and ALT allele. Reported as NumberClippedReadsRef,NumberClippedReadsAlt")


	global opts
	opts = parser.parse_args()

	Genotype_opts = [opts.iAlleleDepth, opts.TotalDPUnfilter, opts.AlleleQscore, opts.AlleleMeanMappingQuality, opts.AlleleSuboptimalAlignmentScore,
					opts.AlleleUnMappedReads, opts.AlleleMappingQualityZero, opts.AlleleNotPrimaryAlignment, opts.AlleleSupplementaryAlignment,
					opts.AlleleNotPaired, opts.AlleleNotProperPaired, opts.DeltaDuplicate, opts.iEvaDepth, opts.AlleleClippedReads,
					opts.AlleleDuplicateDepth, opts.AlleleMeanAlignmentScore, opts.AlleleSuboptimalAlignmentScoreZero]

	Bam_opts = [opts.StrandBiasReads, opts.NotPrimaryAlignment, opts.MeanMappingQuality, opts.SuboptimalAlignmentScore,
				opts.SupplementaryAlignment, opts.AlignmentScore, opts.MappingQualityZero,
				opts.NotProperPairedReads, opts.NotPairedReads, opts.UnMappedReads, opts.TotalDupReads]

	Sequence_opts = [opts.SimpleRepeat, opts.SimpleRepeatUnit, opts.SimpleRepeatLength,
					opts.RepeatMasker, opts.PseudoNucleotidesComposition, opts.gcContent, opts.SimpleRepeatList]

	iCheck.Check_Parameters('-WS, --WindowSize',opts.WindowSize,20,600)
	iCheck.Check_SimpleRepeatList('-SRList, --SimpleRepeatList',opts)
	iCheck.Check_list(opts, Bam_opts, Genotype_opts)

	out = open(opts.outfile,'w')
	Simple_Repeat_list = []
	Reference = opts.reference
	bam_list = opts.list
	Header = []
	H_CHR = []

	iEVA_Header, Head_CHR = vcfwriter.iHeader(Header, H_CHR, opts)

	if opts.SimpleRepeatList:
		with open(opts.SimpleRepeatList) as repeat_list:
			for rep in repeat_list:
				Simple_Repeat_list += [rep.rstrip()]

	out.write('\n'.join(iEVA_Header) + '\n' + '\n'.join(Head_CHR) + '\n')
	Head_CHR = Head_CHR[0].split('\t')

	Sample_tot = Head_CHR[Head_CHR.index('FORMAT')+1:]

	if opts.list:
		Sample_list = Head_CHR[Head_CHR.index('FORMAT')+1:]
		Sample_dict, Sample_excluded = iCheck.Check_Bam(Sample_list, bam_list)
	if not opts.list:
		Sample_dict = {}

	variant = vcfwriter.iAnnotation(Sample_dict, out, opts, Sequence_opts, Bam_opts, Genotype_opts, Head_CHR, Reference, Sample_tot, Simple_Repeat_list)

	out.close()

	print '\nExatraction: Done\n'