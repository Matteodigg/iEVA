import pysam
import sys

def Check_Parameters(arg,opts,min,max):

	if opts > max or opts < min:
		sys.exit('\nArgument ' + arg + ' :' + ' expecated value in range ' + str(min) + '-' + str(max) + '\n')


def Check_SimpleRepeatList(arg,opts):

	if opts.SimpleRepeatList and not opts.SimpleRepeat:
		sys.exit('\nTo use ' + arg + ' please enable -SR, --SimpleRepeats option' + '\n')


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
			print '\n' + "Bam file " + str(elem) + " does not match any sample in genptype vcf field."
	if len(Sample_list) != 0:
		for elem in Sample_list:
			print '\n' + "Sample " + str(elem) + " will not be annotated. Missing bam file." + '\n'

	return Sample_dict,Sample_list

def Check_list(opts,Bam_opts,Genotype_opts):
	
	if any(Bam_opts) or any(Genotype_opts):
		if opts.list:
			pass
		else:
			sys.exit('\n' + 'To Enable bam features insert path to bam list with command -L (--list)' + '\n')