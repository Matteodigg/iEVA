import pysam
import sys

def Check_Parameters(arg,opts,min,max):

	if opts > max or opts < min:
		sys.exit('\nArgument %s : expected value in range %s - %s\n' % (arg,str(min),str(max)))


def Check_SimpleRepeatList(arg,opts):

	if opts.SimpleRepeatList and not opts.SimpleRepeat:
		sys.exit('\nTo use %s please enable -SR, --SimpleRepeats option' % (arg))


def Check_Zero(number):

	if number == 0.0:
		number = int(number)
	else:
		pass
	return number


def Check_Bam(Sample_list, bam_list):

	Sample_dict = {}
	File_list = []
	File_match = []
	path_list = open(bam_list, 'r')

	for path in path_list:
		path = path.rstrip()
		File_list += [path]
		bamfile = pysam.AlignmentFile(path, "rb")
		bam_header = bamfile.header
		Header_Sample_Name = bam_header.get('RG')[0].get('SM')

		if Header_Sample_Name in Sample_list:
			Sample_dict[Header_Sample_Name] = path
			Sample_list.remove(Header_Sample_Name)
			File_match += [path]

		else:
			continue

	path_list.close()

	unmatched_file = set(File_list) - set(File_match)

	if len(unmatched_file) != 0:
		for elem in unmatched_file:
			print '\nBam file %s does not match any sample in genotype vcf field.\n' % (elem)
	if len(Sample_list) != 0:
		for elem in Sample_list:
			print '\nSample %s will not be annotated. Missing bam file.\n' % (elem)

	return Sample_dict,Sample_list


def Check_list(opts,Bam_opts,Genotype_opts):
	
	if any(Bam_opts) or any(Genotype_opts):
		if opts.list:
			pass
		else:
			sys.exit('\nTo Enable bam features insert path to bam list with command -L (--list)\n')

def Check_All(opts):	

	if opts.AllSequence:
		opts.SimpleRepeat = True
		opts.SimpleRepeatUnit = True
		opts.SimpleRepeatLength = True
		opts.RepeatMasker = True
		opts.PseudoNucleotidesComposition = True
		opts.gcContent = True
		opts.VariantClass = True

	if opts.AllBam:
		opts.NotPrimaryAlignment = True
		opts.MeanMappingQuality = True
		opts.SuboptimalAlignmentScore = True
		opts.SupplementaryAlignment = True
		opts.AlignmentScore = True
		opts.MappingQualityZero = True
		opts.MateIsUnmapped = True
		opts.NotProperPairedReads = True
		opts.NotPairedReads = True
		opts.UnMappedReads = True
		opts.TotalDupReads = True

	if opts.AllGenotype:	
		opts.AlleleDepth = True
		opts.AlleleQscore = True
		opts.AlleleMeanMappingQuality = True
		opts.AlleleSuboptimalAlignmentScore = True
		opts.AlleleMappingQualityZero = True
		opts.iEvaDepth = True
		opts.AlleleClippedReads = True
		opts.AlleleMeanAlignmentScore = True
		opts.AlleleSuboptimalAlignmentScoreZero = True
		opts.StrandBias = True
		opts.StrandBiasDepth = True
		opts.iBaseQualRankSumTest = True
		opts.iMapQualRankSumTest = True
		opts.iReadPosRankSumTest = True
		opts.iAltNormReadPos = True
		opts.iClipRankSumTest = True
		opts.iBaseQualValAround = True
		opts.AlleleFrequency = True
		

	