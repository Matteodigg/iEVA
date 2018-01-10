import pysam
import scipy.stats as stats
import iInfo
import multiprocessing
import warnings

def Sequence_Annotator(variant, Sample_dict, H_CHR, Reference, Variant_Class, SNVMinBaseQuality, SNVMinMappingQuality, IndelMinBaseQuality, IndelMinMappingQuality, Bam_opts, Genotype_opts):

	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]

	variant_lenght = max(len(REF), len(ALT))

	STOP = int(variant[H_CHR.index('POS')]) + variant_lenght - 1

	Sample_Stat = {}

	# pool = multiprocessing.Pool(processes=4)
	# results = [pool.apply(Prova_multirpocess, args=(CHR, POS, REF, ALT, variant_lenght, STOP, Sample_Stat, sample, variant, Sample_dict, H_CHR, Reference, Variant_Class, SNVMinBaseQuality, SNVMinMappingQuality, IndelMinBaseQuality, IndelMinMappingQuality, Bam_opts, Genotype_opts,)) for sample in Sample_dict.keys()]
	# print results
	# return Sample_Stat

	#DA PROVARE UN MULTIPROCESS SULLE KEYS
	#from multiprocessing import Pool
	#p = Pool()                                   #number of processes = number of CPUs
	#keys, values= zip(*d.iteritems())            #ordered keys and values
	#processed_values= p.map( f, values )         #apply the function f to each set and wait for result
	#then proceed to join the three sets

	for sample in Sample_dict.keys():

		Allele = ''
		Exact_Match = []
		Variant_Stat = {}
		Reads_Info = {}
		Genotype_Info = {}
		Reads_Info['Total_Reads_No_Dup'] = 0
		Reads_Info['Total_Reads_Unfilter'] = 0
		Reads_Info['Coverage'] = 0
		Reads_Info['Duplicate_reads'] = 0
		Reads_Info['Not_Paired_Reads'] = 0
		Reads_Info['Not_Proper_Paired_Reads'] = 0
		Reads_Info['Dup_Read_Alt'] = 0
		Reads_Info['Dup_Read_Ref'] = 0
		Reads_Info['Alignment_Score'] = 0
		Reads_Info['Suboptimal_Alignment_Score'] = 0
		Reads_Info['Unmapped_reads'] = 0
		Reads_Info['Total_Mapping_Quality'] = 0
		Reads_Info['Mapping_Quality_Zero'] = 0
		Reads_Info['Ref_Mapping_Quality'] = 0
		Reads_Info['Alt_Mapping_Quality'] = 0
		Reads_Info['Read_Ref'] = 0
		Reads_Info['Read_Alt'] = 0
		Reads_Info['Base_Alt_Qual'] = 0
		Reads_Info['Base_Ref_Qual'] = 0
		Reads_Info['Clipped_Reads_Ref'] = 0
		Reads_Info['Clipped_Reads_Alt'] = 0
		Reads_Info['Supplementary_Align'] = 0
		Reads_Info['Not_Primary_Align'] = 0
		Reads_Info['MeanRefPos'] = 0
		Reads_Info['MeanAltPos']=0
		Reads_Info['REF+'] = 0
		Reads_Info['REF-'] = 0
		Reads_Info['ALT+'] = 0
		Reads_Info['ALT-'] = 0
		Reads_Info['QCfail'] = 0

		Reads_Info['UnMap_Ref'] = 0
		Reads_Info['UnMap_Alt'] = 0
		Reads_Info['MapQ0_Ref'] = 0
		Reads_Info['MapQ0_Alt'] = 0
		Reads_Info['NPA_Ref'] = 0
		Reads_Info['NPA_Alt'] = 0
		Reads_Info['Suppl_Ref'] = 0
		Reads_Info['Suppl_Alt'] = 0
		Reads_Info['NP_Ref'] = 0
		Reads_Info['NP_Alt'] = 0
		Reads_Info['NPP_Ref'] = 0
		Reads_Info['NPP_Alt'] = 0
		Reads_Info['AS_Ref'] = 0
		Reads_Info['AS_Alt'] = 0
		Reads_Info['XS_Ref'] = 0
		Reads_Info['XS_Alt'] = 0
		Reads_Info['XS0_Ref'] = 0
		Reads_Info['XS0_Alt'] = 0

		Reads_Info['QualRankRef'] = []
		Reads_Info['QualRankAlt'] = []
		Reads_Info['MappingRankRef'] = []
		Reads_Info['MappingRankAlt'] = []
		Reads_Info['PosRankRef'] = []
		Reads_Info['PosRankAlt'] = []



		if Variant_Class == 'Alt' or Variant_Class == '.':

			Reads_Info['Coverage'] = '.'
			Reads_Info['Dup_Read_Alt'] = '.'
			Reads_Info['Dup_Read_Ref'] = '.'
			Reads_Info['Read_Ref'] = '.'
			Reads_Info['Read_Alt'] = '.'
			Reads_Info['Base_Alt_Qual'] = '.'
			Reads_Info['Base_Ref_Qual'] = '.'
			Reads_Info['Clipped_Reads_Ref'] = '.'
			Reads_Info['Clipped_Reads_Alt'] = '.'
			Reads_Info['UnMap_Ref'] = '.'
			Reads_Info['UnMap_Alt'] = '.'
			Reads_Info['MapQ0_Ref'] = '.'
			Reads_Info['MapQ0_Alt'] = '.'
			Reads_Info['NPA_Ref'] = '.'
			Reads_Info['NPA_Alt'] = '.'
			Reads_Info['Suppl_Ref'] = '.'
			Reads_Info['Suppl_Alt'] = '.'
			Reads_Info['NP_Ref'] = '.'
			Reads_Info['NP_Alt'] = '.'
			Reads_Info['NPP_Ref'] = '.'
			Reads_Info['NPP_Alt'] = '.'
			Reads_Info['AS_Ref'] = '.'
			Reads_Info['AS_Alt'] = '.'
			Reads_Info['XS_Ref'] = '.'
			Reads_Info['XS_Alt'] = '.'
			Reads_Info['XS0_Ref'] = '.'
			Reads_Info['XS0_Alt'] = '.'

		bampile = pysam.AlignmentFile(Sample_dict.get(sample), "rb")

		if any(Bam_opts) or any(Genotype_opts):
			for pileupcolumn in bampile.pileup(CHR, POS, POS+1, stepper='nofilter'):
				if pileupcolumn.pos == POS:
					for pileupread in pileupcolumn.pileups:
						Reads_Info = iInfo.Extract_Reds_Info(pileupread, Reads_Info)
						if any(Genotype_opts):

							if Variant_Class == 'snv':
								if not pileupread.is_del and not pileupread.is_refskip:
									if pileupread.alignment.query_sequence[pileupread.query_position] == ALT:

										Allele = 'ALT'
										if pileupread.alignment.query_qualities[pileupread.query_position] < SNVMinBaseQuality \
										or pileupread.alignment.mapping_quality < SNVMinMappingQuality:
											continue
										elif pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Alt'] += 1
										else:
											Reads_Info = iInfo.Genotype_Reads_Info(pileupread, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match)

									elif pileupread.alignment.query_sequence[pileupread.query_position] == REF:
										
										Allele = 'REF'
										if pileupread.alignment.query_qualities[pileupread.query_position] < SNVMinBaseQuality \
										or pileupread.alignment.mapping_quality < SNVMinMappingQuality:
											continue
										elif pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Ref'] += 1
										else:
											Reads_Info = iInfo.Genotype_Reads_Info(pileupread, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match)

							elif Variant_Class == 'Del':

								Read_Tuple = pileupread.alignment.get_aligned_pairs(with_seq=True)
								Start_Index = 0
								Stop_Index = 0

								for index, val in enumerate(Read_Tuple):
									if val[1] == POS and val[0] is not None and val[2].isupper() and val[2] == REF[0]:
										Start_Index = index
									if val[1] == STOP+1 and val[0] is not None and val[2].isupper():
										Stop_Index = index
										break

								Exact_Match = Read_Tuple[Start_Index:Stop_Index]

								try:
									Indel_Sequence = pileupread.alignment.query_sequence[Exact_Match[0][0]:(Exact_Match[-1][0])]
								except:
									Indel_Sequence = ''

								if Exact_Match == [] or len(Exact_Match) > variant_lenght+1:
									continue

								else:
									if pileupread.indel == -(variant_lenght-1) and not pileupread.is_refskip \
									and not pileupread.is_tail and len(Indel_Sequence) == 1:

										Allele = 'ALT'
										if pileupread.alignment.query_qualities[pileupread.query_position] < IndelMinBaseQuality \
										or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
											continue

										elif pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Alt'] += 1
										else:
											Reads_Info = iInfo.Genotype_Reads_Info(pileupread, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match)

									elif not pileupread.is_refskip and not pileupread.query_position == None \
									and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:

										Allele = 'REF'
										for values in Exact_Match[1:-1]:
											if values[0] is not None and values[1] is not None and values[2].isupper():
												Not_Ref = False
											else:
												Not_Ref = True
												break

										if pileupread.alignment.query_qualities[pileupread.query_position] < IndelMinBaseQuality \
										or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
											Not_Ref = True

										if Not_Ref is False:
											if pileupread.alignment.is_duplicate:
												Reads_Info['Dup_Read_Ref'] += 1

											else:
												Reads_Info = iInfo.Genotype_Reads_Info(pileupread, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match)


							elif Variant_Class == 'Ins':


								Read_Tuple = pileupread.alignment.get_aligned_pairs(with_seq=True)
								Start_Index = 0
								Stop_Index = 0
								for index, val in enumerate(Read_Tuple):
									if val[1] == POS and val[0] is not None and val[2].isupper():
										Start_Index = index

									if val[1] == (POS+2) and val[0] is not None and val[2].isupper():
										Stop_Index = index
										break

								Exact_Match = Read_Tuple[Start_Index:Stop_Index]

								try:
									Indel_Sequence = pileupread.alignment.query_sequence[Exact_Match[0][0]:(Exact_Match[-1][0])]
								except:
									Indel_Sequence = ''

								if Exact_Match == []:
									continue

								elif len(Exact_Match) == variant_lenght+1:
									if pileupread.indel == (variant_lenght-1) and not pileupread.is_refskip \
									and not pileupread.is_tail and Indel_Sequence == ALT \
									and Exact_Match[-1][0] is not None and Exact_Match[-1][2].isupper() \
									and Exact_Match[0][0] is not None and Exact_Match[0][2].isupper():

										Allele = 'ALT'
										if float(sum(pileupread.alignment.query_qualities[Exact_Match[1][0]:Exact_Match[-1][0]]))/float(variant_lenght-1) < IndelMinBaseQuality \
										or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
											continue

										elif pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Alt'] += 1

										else:
											Reads_Info = iInfo.Genotype_Reads_Info(pileupread, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match)

								elif len(Exact_Match) == 2:
									if not pileupread.is_refskip and Exact_Match[0][0] is not None and Exact_Match[0][2].isupper() \
									and Exact_Match[-1][0] is not None and Exact_Match[-1][2].isupper():

										Allele = 'REF'
										if float(sum(pileupread.alignment.query_qualities[Exact_Match[0][0]:(Exact_Match[-1][0])+1]))/float(2) < IndelMinBaseQuality \
										or pileupread.alignment.mapping_quality < IndelMinMappingQuality:
											continue

										elif pileupread.alignment.is_duplicate:
											Reads_Info['Dup_Read_Ref'] += 1

										else:
											Reads_Info = iInfo.Genotype_Reads_Info(pileupread, Reads_Info, Allele, Variant_Class, variant_lenght, Exact_Match)

								else:
									continue

		try:
			Variant_Stat['Mean_Alignment_Score'] = round(float(Reads_Info['Alignment_Score'])/float(Reads_Info['Total_Reads_No_Dup']),3)
		except:
			Variant_Stat['Mean_Alignment_Score'] = '.'
		try:
			Variant_Stat['Mean_Suboptimal_Alignment_Score'] = round(float(Reads_Info['Suboptimal_Alignment_Score'])/float(Reads_Info['Total_Reads_No_Dup']),3)
		except:
			Variant_Stat['Mean_Suboptimal_Alignment_Score'] = '.'
		try:
			Variant_Stat['Percentage_Unmapped_Reads'] = round(float(Reads_Info['Unmapped_reads'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Unmapped_Reads'] = '.'
		try:
			Variant_Stat['Percentage_Dup_Ref'] = round(float(Reads_Info['Dup_Read_Ref']) / float((Reads_Info['Dup_Read_Ref']+Reads_Info['Read_Ref'])),4)
		except:
			Variant_Stat['Percentage_Dup_Ref'] = '.'
		try:
			Variant_Stat['Percentage_Dup_Alt'] = round(float(Reads_Info['Dup_Read_Alt']) / float((Reads_Info['Dup_Read_Alt']+Reads_Info['Read_Alt'])),4)
		except:
			Variant_Stat['Percentage_Dup_Alt'] = '.'
		try:
			Variant_Stat['Delta_Duplicate'] = round(Variant_Stat['Percentage_Dup_Ref']-Variant_Stat['Percentage_Dup_Alt'],4)
		except:
			Variant_Stat['Delta_Duplicate'] = '.'
		try:
			Variant_Stat['Total_Duplicate_Reads'] = round(float(Reads_Info['Duplicate_reads'])/float(Reads_Info['Total_Reads_Unfilter']),4)
		except:
			Variant_Stat['Total_Duplicate_Reads'] = '0'
		try:
			Variant_Stat['Total_Reads_Unfilter'] = Reads_Info['Total_Reads_Unfilter']
		except:
			Variant_Stat['Total_Reads_Unfilter'] = '.'
		try:
			Variant_Stat['Percentage_Supplementary_Align'] = round(float(Reads_Info['Supplementary_Align'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Supplementary_Align'] = '0'
		try:
			Variant_Stat['Percentage_Not_Primary_Alignment_Reads'] = round(float(Reads_Info['Not_Primary_Align'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Not_Primary_Alignment_Reads'] = '0'
		try:
			Variant_Stat['Percentage_Not_Paired_Reads'] = round(float(Reads_Info['Not_Paired_Reads'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Not_Paired_Reads'] = '0'
		try:
			Variant_Stat['Percentage_Not_Proper_Paired_Reads'] = round(float(Reads_Info['Not_Proper_Paired_Reads'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Percentage_Not_Proper_Paired_Reads'] = '0'
		try:
			Variant_Stat['Mapping_Quality_Zero'] = round(float(Reads_Info['Mapping_Quality_Zero'])/float(Reads_Info['Total_Reads_No_Dup']),4)
		except:
			Variant_Stat['Mapping_Quality_Zero'] = '0'
		try:
			Variant_Stat['Total_Mean_Mapping_Quality'] = round(float(Reads_Info['Total_Mapping_Quality'])/float(Reads_Info['Total_Reads_Unfilter']),4)
		except:
			Variant_Stat['Total_Mean_Mapping_Quality'] = '.'
		try:
			Variant_Stat['Coverage'] = int(Reads_Info['Read_Ref'])+int(Reads_Info['Read_Alt'])
		except:
			Variant_Stat['Coverage'] = '.'
		try:
			Variant_Stat['Read_Ref'] = Reads_Info['Read_Ref']
		except:
			Variant_Stat['Read_Ref'] = '.'
		try:
			Variant_Stat['Read_Alt'] = Reads_Info['Read_Alt']
		except:
			Variant_Stat['Read_Alt'] = '.'
		try:
			oddsratio, SBR = stats.fisher_exact([[Reads_Info['REF+'], Reads_Info['REF-']], [Reads_Info['ALT+'], Reads_Info['ALT-']]])
		except:
			SBR = '.'
		try:
			Variant_Stat['REF+'] = Reads_Info['REF+']
		except:
			Variant_Stat['REF+'] = '.'
		try:
			Variant_Stat['REF-'] = Reads_Info['REF-']
		except:
			Variant_Stat['REF-'] = '.'
		try:
			Variant_Stat['ALT+'] = Reads_Info['ALT+']
		except:
			Variant_Stat['ALT+'] = '.'
		try:
			Variant_Stat['ALT-'] = Reads_Info['ALT-']
		except:
			Variant_Stat['ALT-'] = '.'
		try:
			Variant_Stat['Strand_Bias_Reads'] = round(SBR,4)
		except:
			Variant_Stat['Strand_Bias_Reads'] = '.'
		try:
			Variant_Stat['Qual_Alt'] = round(float(Reads_Info['Base_Alt_Qual'])/float(Reads_Info['Read_Alt']),2)
		except:
			Variant_Stat['Qual_Alt'] = '.'
		try:
			Variant_Stat['Qual_Ref'] = round(float(Reads_Info['Base_Ref_Qual'])/float(Reads_Info['Read_Ref']),2)
		except:
			Variant_Stat['Qual_Ref'] = '.'
		try:
			Variant_Stat['Number_Read_Dup_Alt'] = Reads_Info['Dup_Read_Alt']
		except:
			Variant_Stat['Number_Read_Dup_Alt'] = '.'
		try:
			Variant_Stat['Number_Read_Dup_Ref'] = Reads_Info['Dup_Read_Ref']
		except:
			Variant_Stat['Number_Read_Dup_Ref'] = '.'
		try:
			Variant_Stat['Clipped_Reads_Ref'] = Reads_Info['Clipped_Reads_Ref']
		except:
			Variant_Stat['Clipped_Reads_Ref'] = '.'
		try:
			Variant_Stat['Clipped_Reads_Alt'] = Reads_Info['Clipped_Reads_Alt']
		except:
			Variant_Stat['Clipped_Reads_Alt'] = '.'
		try:
			Variant_Stat['Ref_Mean_Mapping_Quality'] = round(float(Reads_Info['Ref_Mapping_Quality'])/float(Reads_Info['Read_Ref']),4)
		except:
			Variant_Stat['Ref_Mean_Mapping_Quality'] = '.'
		try:
			Variant_Stat['Alt_Mean_Mapping_Quality'] = round(float(Reads_Info['Alt_Mapping_Quality'])/float(Reads_Info['Read_Alt']),4)
		except:
			Variant_Stat['Alt_Mean_Mapping_Quality'] = '.'
		try:
			Variant_Stat['UnMap_Ref'] = Reads_Info['UnMap_Ref']
		except:
			Variant_Stat['UnMap_Ref'] = '.'
		try:
			Variant_Stat['UnMap_Alt'] = Reads_Info['UnMap_Alt']
		except:
			Variant_Stat['UnMap_Alt'] = '.'
		try:
			Variant_Stat['MapQ0_Ref'] = Reads_Info['MapQ0_Ref']
		except:
			Variant_Stat['MapQ0_Ref'] = '.'
		try:
			Variant_Stat['MapQ0_Alt'] = Reads_Info['MapQ0_Alt']
		except:
			Variant_Stat['MapQ0_Alt'] = '.'
		try:
			Variant_Stat['Not_Primary_Align_Ref'] = Reads_Info['NPA_Ref']
		except:
			Variant_Stat['Not_Primary_Align_Ref'] = '.'
		try:
			Variant_Stat['Not_Primary_Align_Alt'] = Reads_Info['NPA_Alt']
		except:
			Variant_Stat['Not_Primary_Align_Alt'] = '.'
		try:
			Variant_Stat['Not_Paired_Ref'] = Reads_Info['NP_Ref']
		except:
			Variant_Stat['Not_Paired_Ref'] = '.'
		try:
			Variant_Stat['Not_Paired_Alt'] = Reads_Info['NP_Alt']
		except:
			Variant_Stat['Not_Paired_Alt'] = '.'
		try:
			Variant_Stat['Not_Proper_Paired_Ref'] = Reads_Info['NPP_Ref']
		except:
			Variant_Stat['Not_Proper_Paired_Ref'] = '.'
		try:
			Variant_Stat['Not_Proper_Paired_Alt'] = Reads_Info['NPP_Alt']
		except:
			Variant_Stat['Not_Proper_Paired_Alt'] = '.'
		try:
			Variant_Stat['Suppl_Alignment_Ref'] = Reads_Info['Suppl_Ref']
		except:
			Variant_Stat['Suppl_Alignment_Ref'] = '.'
		try:
			Variant_Stat['Suppl_Alignment_Alt'] = Reads_Info['Suppl_Alt']
		except:
			Variant_Stat['Suppl_Alignment_Alt'] = '.'
		try:
			Variant_Stat['AS_Ref'] = round(float(Reads_Info['AS_Ref'])/float(Reads_Info['Read_Ref']),3)
		except:
			Variant_Stat['AS_Ref'] = '.'
		try:
			Variant_Stat['AS_Alt'] = round(float(Reads_Info['AS_Alt'])/float(Reads_Info['Read_Alt']),3)
		except:
			Variant_Stat['AS_Alt'] = '.'
		try:
			Variant_Stat['XS_Ref'] = round(float(Reads_Info['XS_Ref'])/float(Reads_Info['Read_Ref']),3)
		except:
			Variant_Stat['XS_Ref'] = '.'
		try:
			Variant_Stat['XS_Alt'] = round(float(Reads_Info['XS_Alt'])/float(Reads_Info['Read_Alt']),3)
		except:
			Variant_Stat['XS_Alt'] = '.'
		try:
			Variant_Stat['XS0_Ref'] = Reads_Info['XS0_Ref']
		except:
			Variant_Stat['XS0_Ref'] = '.'
		try:
			Variant_Stat['XS0_Alt'] = Reads_Info['XS0_Alt']
		except:
			Variant_Stat['XS0_Alt'] = '.'

		#if opts.QualRankTest:
		try:
			QualRank = stats.ranksums(Reads_Info['QualRankAlt'], Reads_Info['QualRankRef'])
			Variant_Stat['QualRankTest'] = QualRank.statistic
		except:
			Variant_Stat['QualRankTest'] = '.'

		try:
			MappingRank = stats.ranksums(Reads_Info['MappingRankAlt'], Reads_Info['MappingRankRef'])
			Variant_Stat['MappingRankTest'] = MappingRank.statistic
		except:
			Variant_Stat['MappingRankTest'] = '.'

		try:
			PosRank = stats.ranksums(Reads_Info['PosRankAlt'], Reads_Info['PosRankRef'])
			Variant_Stat['PosaRankTest'] = PosRank.statistic
		except:
			Variant_Stat['PosaRankTest'] = '.'

		#if Variant_Stat['Total_Duplicate_Reads'] != 0 and Variant_Stat['Total_Duplicate_Reads'] != 0.0:
		#	print CHR
		#	print POS
		#	print sample
		#	#print Reads_Info['PosRankAlt']
		#	print Reads_Info['Dup_Read_Alt']
		#	print Reads_Info['Dup_Read_Ref']
		#	print Variant_Stat['Percentage_Dup_Ref']
		#	print Variant_Stat['Percentage_Dup_Alt']
		#	print Variant_Stat['Total_Duplicate_Reads']
		#	#print Reads_Info['QCfail']
		#	#print Variant_Stat['QualRankTest']
		#	#print Variant_Stat['MappingRankTest']
		#	print '\n'

		#print Variant_Stat['Total_Duplicate_Reads']

		#RISOLVI DUP

		#INERISCI is_qcfail per failing vendor

		Sample_Stat[sample] = Variant_Stat
		Reads_Info = {}

	return Sample_Stat