from pyfasta import Fasta
from repDNA.nac import Kmer
from Bio import SeqUtils
import regex

def SimpleRepeats_Finder(Reference, variant, H_CHR, Variant_Class, opts, Repeat_list):

	WindowSize = opts.WindowSize / 2

	Ref = Fasta(Reference)
	sorted(Ref.keys())

	CHR = variant[H_CHR.index('#CHROM')]
	POS = int(variant[H_CHR.index('POS')])-1
	REF = variant[H_CHR.index('REF')]
	ALT = variant[H_CHR.index('ALT')]
	GC = '.'
	RM = '.'
	new_str = []
	Dinucloetides = []
	
	if WindowSize >= POS:
		WindowSize = POS +1

	sequence = Ref[CHR][POS-WindowSize:POS+WindowSize+1]

	seq = sequence.upper()

	if opts.PseudoNucleotidesComposition:

		try:
			kmer = Kmer(k=2,normalize=True)
			Dinuc = kmer.make_kmer_vec([seq])
			Dinucloetides = Dinuc[0]
		except:
			Dinucloetides = '.'

		Combo_Seq = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']

	if opts.SimpleRepeat:

		Homopolymer_pattern = regex.compile(r'(A){4,}|(G){4,}|(C){4,}|(T){4,}')

		Homopolymer_list = [(m.group(), m.groups(), m.start(), m.end(), len(m.group())) for m in regex.finditer(Homopolymer_pattern, seq)]

		RepeatSeq = '0'
		RepeatSeq_Lenght = '.'
		RepeatSeq_Unit = '.'

		for poly in Homopolymer_list:

			End_Hom = poly[3]-1

			if poly[2] <= WindowSize+1 and End_Hom >= WindowSize:
				RepeatSeq = '2'
				RepeatSeq_Lenght = poly[-1]
				RepeatSeq_Unit = [i for i in poly[1] if i is not None]
				break

		if RepeatSeq == '0':
			if opts.SimpleRepeatList:
				for string in Repeat_list:
					new_str += ['('+string+')'+'{3,}']
				Simple_Repeat_pattern = regex.compile('|'.join(new_str))

			else:
				Simple_Repeat_pattern = regex.compile('(AC){3,}|(AG){3,}|(AT){3,}|(CA){3,}|(GA){3,}|(TA){3,}|(CG){3,}|(GC){3,}|(CT){3,}|(TC){3,}|'
					'(GT){3,}|(TG){3,}|(AAC){3,}|(ACA){3,}|(ACC){3,}|(CAA){3,}|(CAC){3,}|(CCA){3,}|(AAG){3,}|(AGA){3,}|(AGG){3,}|(GAA){3,}|(GAG){3,}|'
					'(GGA){3,}|(AAT){3,}|(ATA){3,}|(ATT){3,}|(TAA){3,}|(TAT){3,}|(TTA){3,}|(CCG){3,}|(CGC){3,}|(CGG){3,}|(GCC){3,}|(GCG){3,}|(GGC){3,}|'
					'(CCT){3,}|(CTC){3,}|(CTT){3,}|(TCC){3,}|(TCT){3,}|(TTC){3,}|(GGT){3,}|(GTG){3,}|(GTT){3,}|(TGG){3,}|(TGT){3,}|(TTG){3,}|'
					'(AAAC){3,}|(AACA){3,}|(AACC){3,}|(ACAA){3,}|(ACAC){3,}|(ACCA){3,}|(ACCC){3,}|(CAAA){3,}|(CAAC){3,}|(CACA){3,}|(CACC){3,}|'
					'(CCAA){3,}|(CCAC){3,}|(CCCA){3,}|(AAAG){3,}|(AAGA){3,}|(AAGG){3,}|(AGAA){3,}|(AGAG){3,}|(AGGA){3,}|(AGGG){3,}|(GAAA){3,}|'
					'(GAAG){3,}|(GAGA){3,}|(GAGG){3,}|(GGAA){3,}|(GGAG){3,}|(GGGA){3,}|(AAAT){3,}|(AATA){3,}|(AATT){3,}|(ATAA){3,}|(ATAT){3,}|'
					'(ATTA){3,}|(ATTT){3,}|(TAAA){3,}|(TAAT){3,}|(TATA){3,}|(TATT){3,}|(TTAA){3,}|(TTAT){3,}|(TTTA){3,}|(CCCG){3,}|(CCGC){3,}|'
					'(CCGG){3,}|(CGCC){3,}|(CGCG){3,}|(CGGC){3,}|(CGGG){3,}|(GCCC){3,}|(GCCG){3,}|(GCGC){3,}|(GCGG){3,}|(GGCC){3,}|(GGCG){3,}|'
					'(GGGC){3,}|(CCCT){3,}|(CCTC){3,}|(CCTT){3,}|(CTCC){3,}|(CTCT){3,}|(CTTC){3,}|(CTTT){3,}|(TCCC){3,}|(TCCT){3,}|(TCTC){3,}|'
					'(TCTT){3,}|(TTCC){3,}|(TTCT){3,}|(TTTC){3,}|(GGGT){3,}|(GGTG){3,}|(GGTT){3,}|(GTGG){3,}|(GTGT){3,}|(GTTG){3,}|(GTTT){3,}|'
					'(TGGG){3,}|(TGGT){3,}|(TGTG){3,}|(TGTT){3,}|(TTGG){3,}|(TTGT){3,}|(TTTG){3,}|(AAAAC){3,}|(AAACA){3,}|(AAACC){3,}|(AACAA){3,}|'
					'(AACAC){3,}|(AACCA){3,}|(AACCC){3,}|(ACAAA){3,}|(ACAAC){3,}|(ACACA){3,}|(ACACC){3,}|(ACCAA){3,}|(ACCAC){3,}|(ACCCA){3,}|'
					'(ACCCC){3,}|(CAAAA){3,}|(CAAAC){3,}|(CAACA){3,}|(CAACC){3,}|(CACAA){3,}|(CACAC){3,}|(CACCA){3,}|(CACCC){3,}|(CCAAA){3,}|'
					'(CCAAC){3,}|(CCACA){3,}|(CCACC){3,}|(CCCAA){3,}|(CCCAC){3,}|(CCCCA){3,}|(AAAAG){3,}|(AAAGA){3,}|(AAAGG){3,}|(AAGAA){3,}|'
					'(AAGAG){3,}|(AAGGA){3,}|(AAGGG){3,}|(AGAAA){3,}|(AGAAG){3,}|(AGAGA){3,}|(AGAGG){3,}|(AGGAA){3,}|(AGGAG){3,}|(AGGGA){3,}|'
					'(AGGGG){3,}|(GAAAA){3,}|(GAAAG){3,}|(GAAGA){3,}|(GAAGG){3,}|(GAGAA){3,}|(GAGAG){3,}|(GAGGA){3,}|(GAGGG){3,}|(GGAAA){3,}|'
					'(GGAAG){3,}|(GGAGA){3,}|(GGAGG){3,}|(GGGAA){3,}|(GGGAG){3,}|(GGGGA){3,}|(AAAAT){3,}|(AAATA){3,}|(AAATT){3,}|(AATAA){3,}|'
					'(AATAT){3,}|(AATTA){3,}|(AATTT){3,}|(ATAAA){3,}|(ATAAT){3,}|(ATATA){3,}|(ATATT){3,}|(ATTAA){3,}|(ATTAT){3,}|(ATTTA){3,}|'
					'(ATTTT){3,}|(TAAAA){3,}|(TAAAT){3,}|(TAATA){3,}|(TAATT){3,}|(TATAA){3,}|(TATAT){3,}|(TATTA){3,}|(TATTT){3,}|(TTAAA){3,}|'
					'(TTAAT){3,}|(TTATA){3,}|(TTATT){3,}|(TTTAA){3,}|(TTTAT){3,}|(TTTTA){3,}|(CCCCG){3,}|(CCCGC){3,}|(CCCGG){3,}|(CCGCC){3,}|'
					'(CCGCG){3,}|(CCGGC){3,}|(CCGGG){3,}|(CGCCC){3,}|(CGCCG){3,}|(CGCGC){3,}|(CGCGG){3,}|(CGGCC){3,}|(CGGCG){3,}|(CGGGC){3,}|'
					'(CGGGG){3,}|(GCCCC){3,}|(GCCCG){3,}|(GCCGC){3,}|(GCCGG){3,}|(GCGCC){3,}|(GCGCG){3,}|(GCGGC){3,}|(GCGGG){3,}|(GGCCC){3,}|'
					'(GGCCG){3,}|(GGCGC){3,}|(GGCGG){3,}|(GGGCC){3,}|(GGGCG){3,}|(GGGGC){3,}|(CCCCT){3,}|(CCCTC){3,}|(CCCTT){3,}|(CCTCC){3,}|'
					'(CCTCT){3,}|(CCTTC){3,}|(CCTTT){3,}|(CTCCC){3,}|(CTCCT){3,}|(CTCTC){3,}|(CTCTT){3,}|(CTTCC){3,}|(CTTCT){3,}|(CTTTC){3,}|'
					'(CTTTT){3,}|(TCCCC){3,}|(TCCCT){3,}|(TCCTC){3,}|(TCCTT){3,}|(TCTCC){3,}|(TCTCT){3,}|(TCTTC){3,}|(TCTTT){3,}|(TTCCC){3,}|'
					'(TTCCT){3,}|(TTCTC){3,}|(TTCTT){3,}|(TTTCC){3,}|(TTTCT){3,}|(TTTTC){3,}|(GGGGT){3,}|(GGGTG){3,}|(GGGTT){3,}|(GGTGG){3,}|'
					'(GGTGT){3,}|(GGTTG){3,}|(GGTTT){3,}|(GTGGG){3,}|(GTGGT){3,}|(GTGTG){3,}|(GTGTT){3,}|(GTTGG){3,}|(GTTGT){3,}|(GTTTG){3,}|'
					'(GTTTT){3,}|(TGGGG){3,}|(TGGGT){3,}|(TGGTG){3,}|(TGGTT){3,}|(TGTGG){3,}|(TGTGT){3,}|(TGTTG){3,}|(TGTTT){3,}|(TTGGG){3,}|'
					'(TTGGT){3,}|(TTGTG){3,}|(TTGTT){3,}|(TTTGG){3,}|(TTTGT){3,}|(TTTTG){3,}|(AAAACA){3,}|(AAAACC){3,}|(AAACAA){3,}|(AAACAC){3,}|(AAACCA){3,}|(AAACCC){3,}|(AACAAA){3,}|(AACAAC){3,}|(AACACA){3,}|(AACACC){3,}|(AACCAA){3,}|(AACCAC){3,}|(AACCCA){3,}|(AACCCC){3,}|'
					'(ACAAAA){3,}|(ACAAAC){3,}|(ACAACA){3,}|(ACAACC){3,}|(ACACAA){3,}|(ACACAC){3,}|(ACACCA){3,}|(ACACCC){3,}|(ACCAAA){3,}|(ACCAAC){3,}|'
					'(ACCACA){3,}|(ACCACC){3,}|(ACCCAA){3,}|(ACCCAC){3,}|(ACCCCA){3,}|(ACCCCC){3,}|(CAAAAA){3,}|(CAAAAC){3,}|(CAAACA){3,}|(CAAACC){3,}|'
					'(CAACAA){3,}|(CAACAC){3,}|(CAACCA){3,}|(CAACCC){3,}|(CACAAA){3,}|(CACAAC){3,}|(CACACA){3,}|(CACACC){3,}|(CACCAA){3,}|(CACCAC){3,}|'
					'(CACCCA){3,}|(CACCCC){3,}|(CCAAAA){3,}|(CCAAAC){3,}|(CCAACA){3,}|(CCAACC){3,}|(CCACAA){3,}|(CCACAC){3,}|(CCACCA){3,}|(CCACCC){3,}|'
					'(CCCAAA){3,}|(CCCAAC){3,}|(CCCACA){3,}|(CCCACC){3,}|(CCCCAA){3,}|(CCCCAC){3,}|(AAAAGA){3,}|(AAAAGG){3,}|(AAAGAA){3,}|(AAAGAG){3,}|'
					'(AAAGGA){3,}|(AAAGGG){3,}|(AAGAAA){3,}|(AAGAAG){3,}|(AAGAGA){3,}|(AAGAGG){3,}|(AAGGAA){3,}|(AAGGAG){3,}|(AAGGGA){3,}|(AAGGGG){3,}|'
					'(AGAAAA){3,}|(AGAAAG){3,}|(AGAAGA){3,}|(AGAAGG){3,}|(AGAGAA){3,}|(AGAGAG){3,}|(AGAGGA){3,}|(AGAGGG){3,}|(AGGAAA){3,}|(AGGAAG){3,}|'
					'(AGGAGA){3,}|(AGGAGG){3,}|(AGGGAA){3,}|(AGGGAG){3,}|(AGGGGA){3,}|(AGGGGG){3,}|(GAAAAA){3,}|(GAAAAG){3,}|(GAAAGA){3,}|(GAAAGG){3,}|'
					'(GAAGAA){3,}|(GAAGAG){3,}|(GAAGGA){3,}|(GAAGGG){3,}|(GAGAAA){3,}|(GAGAAG){3,}|(GAGAGA){3,}|(GAGAGG){3,}|(GAGGAA){3,}|(GAGGAG){3,}|'
					'(GAGGGA){3,}|(GAGGGG){3,}|(GGAAAA){3,}|(GGAAAG){3,}|(GGAAGA){3,}|(GGAAGG){3,}|(GGAGAA){3,}|(GGAGAG){3,}|(GGAGGA){3,}|(GGAGGG){3,}|'
					'(GGGAAA){3,}|(GGGAAG){3,}|(GGGAGA){3,}|(GGGAGG){3,}|(GGGGAA){3,}|(GGGGAG){3,}|(AAAATA){3,}|(AAAATT){3,}|(AAATAA){3,}|(AAATAT){3,}|'
					'(AAATTA){3,}|(AAATTT){3,}|(AATAAA){3,}|(AATAAT){3,}|(AATATA){3,}|(AATATT){3,}|(AATTAA){3,}|(AATTAT){3,}|(AATTTA){3,}|(AATTTT){3,}|'
					'(ATAAAA){3,}|(ATAAAT){3,}|(ATAATA){3,}|(ATAATT){3,}|(ATATAA){3,}|(ATATAT){3,}|(ATATTA){3,}|(ATATTT){3,}|(ATTAAA){3,}|(ATTAAT){3,}|'
					'(ATTATA){3,}|(ATTATT){3,}|(ATTTAA){3,}|(ATTTAT){3,}|(ATTTTA){3,}|(ATTTTT){3,}|(TAAAAA){3,}|(TAAAAT){3,}|(TAAATA){3,}|(TAAATT){3,}|'
					'(TAATAA){3,}|(TAATAT){3,}|(TAATTA){3,}|(TAATTT){3,}|(TATAAA){3,}|(TATAAT){3,}|(TATATA){3,}|(TATATT){3,}|(TATTAA){3,}|(TATTAT){3,}|'
					'(TATTTA){3,}|(TATTTT){3,}|(TTAAAA){3,}|(TTAAAT){3,}|(TTAATA){3,}|(TTAATT){3,}|(TTATAA){3,}|(TTATAT){3,}|(TTATTA){3,}|(TTATTT){3,}|'
					'(TTTAAA){3,}|(TTTAAT){3,}|(TTTATA){3,}|(TTTATT){3,}|(TTTTAA){3,}|(TTTTAT){3,}|(CCCCGC){3,}|(CCCCGG){3,}|(CCCGCC){3,}|(CCCGCG){3,}|'
					'(CCCGGC){3,}|(CCCGGG){3,}|(CCGCCC){3,}|(CCGCCG){3,}|(CCGCGC){3,}|(CCGCGG){3,}|(CCGGCC){3,}|(CCGGCG){3,}|(CCGGGC){3,}|(CCGGGG){3,}|'
					'(CGCCCC){3,}|(CGCCCG){3,}|(CGCCGC){3,}|(CGCCGG){3,}|(CGCGCC){3,}|(CGCGCG){3,}|(CGCGGC){3,}|(CGCGGG){3,}|(CGGCCC){3,}|(CGGCCG){3,}|'
					'(CGGCGC){3,}|(CGGCGG){3,}|(CGGGCC){3,}|(CGGGCG){3,}|(CGGGGC){3,}|(CGGGGG){3,}|(GCCCCC){3,}|(GCCCCG){3,}|(GCCCGC){3,}|(GCCCGG){3,}|'
					'(GCCGCC){3,}|(GCCGCG){3,}|(GCCGGC){3,}|(GCCGGG){3,}|(GCGCCC){3,}|(GCGCCG){3,}|(GCGCGC){3,}|(GCGCGG){3,}|(GCGGCC){3,}|(GCGGCG){3,}|'
					'(GCGGGC){3,}|(GCGGGG){3,}|(GGCCCC){3,}|(GGCCCG){3,}|(GGCCGC){3,}|(GGCCGG){3,}|(GGCGCC){3,}|(GGCGCG){3,}|(GGCGGC){3,}|(GGCGGG){3,}|'
					'(GGGCCC){3,}|(GGGCCG){3,}|(GGGCGC){3,}|(GGGCGG){3,}|(GGGGCC){3,}|(GGGGCG){3,}|(CCCCTC){3,}|(CCCCTT){3,}|(CCCTCC){3,}|(CCCTCT){3,}|'
					'(CCCTTC){3,}|(CCCTTT){3,}|(CCTCCC){3,}|(CCTCCT){3,}|(CCTCTC){3,}|(CCTCTT){3,}|(CCTTCC){3,}|(CCTTCT){3,}|(CCTTTC){3,}|(CCTTTT){3,}|'
					'(CTCCCC){3,}|(CTCCCT){3,}|(CTCCTC){3,}|(CTCCTT){3,}|(CTCTCC){3,}|(CTCTCT){3,}|(CTCTTC){3,}|(CTCTTT){3,}|(CTTCCC){3,}|(CTTCCT){3,}|'
					'(CTTCTC){3,}|(CTTCTT){3,}|(CTTTCC){3,}|(CTTTCT){3,}|(CTTTTC){3,}|(CTTTTT){3,}|(TCCCCC){3,}|(TCCCCT){3,}|(TCCCTC){3,}|(TCCCTT){3,}|'
					'(TCCTCC){3,}|(TCCTCT){3,}|(TCCTTC){3,}|(TCCTTT){3,}|(TCTCCC){3,}|(TCTCCT){3,}|(TCTCTC){3,}|(TCTCTT){3,}|(TCTTCC){3,}|(TCTTCT){3,}|'
					'(TCTTTC){3,}|(TCTTTT){3,}|(TTCCCC){3,}|(TTCCCT){3,}|(TTCCTC){3,}|(TTCCTT){3,}|(TTCTCC){3,}|(TTCTCT){3,}|(TTCTTC){3,}|(TTCTTT){3,}|'
					'(TTTCCC){3,}|(TTTCCT){3,}|(TTTCTC){3,}|(TTTCTT){3,}|(TTTTCC){3,}|(TTTTCT){3,}|(GGGGTG){3,}|(GGGGTT){3,}|(GGGTGG){3,}|(GGGTGT){3,}|'
					'(GGGTTG){3,}|(GGGTTT){3,}|(GGTGGG){3,}|(GGTGGT){3,}|(GGTGTG){3,}|(GGTGTT){3,}|(GGTTGG){3,}|(GGTTGT){3,}|(GGTTTG){3,}|(GGTTTT){3,}|'
					'(GTGGGG){3,}|(GTGGGT){3,}|(GTGGTG){3,}|(GTGGTT){3,}|(GTGTGG){3,}|(GTGTGT){3,}|(GTGTTG){3,}|(GTGTTT){3,}|(GTTGGG){3,}|(GTTGGT){3,}|'
					'(GTTGTG){3,}|(GTTGTT){3,}|(GTTTGG){3,}|(GTTTGT){3,}|(GTTTTG){3,}|(GTTTTT){3,}|(TGGGGG){3,}|(TGGGGT){3,}|(TGGGTG){3,}|(TGGGTT){3,}|'
					'(TGGTGG){3,}|(TGGTGT){3,}|(TGGTTG){3,}|(TGGTTT){3,}|(TGTGGG){3,}|(TGTGGT){3,}|(TGTGTG){3,}|(TGTGTT){3,}|(TGTTGG){3,}|(TGTTGT){3,}|'
					'(TGTTTG){3,}|(TGTTTT){3,}|(TTGGGG){3,}|(TTGGGT){3,}|(TTGGTG){3,}|(TTGGTT){3,}|(TTGTGG){3,}|(TTGTGT){3,}|(TTGTTG){3,}|(TTGTTT){3,}|'
					'(TTTGGG){3,}|(TTTGGT){3,}|(TTTGTG){3,}|(TTTGTT){3,}|(TTTTGG){3,}|(TTTTGT){3,}')

			Simple_Repeat_list = [(m.group(), m.groups(), m.start(), m.end(), len(m.group())) for m in regex.finditer(Simple_Repeat_pattern, seq)]

			for pattern in Simple_Repeat_list:

				End_Rep = pattern[3]-1

				if pattern[2] <= WindowSize+1 and End_Rep >= WindowSize:
					RepeatSeq = '1'
					RepeatSeq_Lenght = pattern[-1]
					RepeatSeq_Unit = [i for i in pattern[1] if i is not None]
					break

	if opts.RepeatMasker:
		if sequence[WindowSize].isupper():
			RM = '0'
		else:
			RM = '1'

	if opts.gcContent:
		try:
			GC = str(round(SeqUtils.GC(seq), 3))
		except:
			GC = '.'

	return RepeatSeq, RepeatSeq_Lenght, RepeatSeq_Unit[0], Dinucloetides, RM, GC