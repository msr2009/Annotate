"""
annotate.py

takes input of BED formatted file with mutation locations (from mutation_caller.py) and annotates SNPs

python annotate.py -f <BED file containing mutations> -s <ORF sequences> -n <non-coding GFF file for annotations> -g <Genome sequence as FASTA>

"""

def main(BED, orfs, noncoding_file, genome_file):

	"""GET DATA TOGETHER"""

	#gather data about genes:
	#orf id, start, stop (both 5'-3'), exons, introns, protein length, chr
	genes = {}
	for record in SeqIO.parse( open(orfs, 'r'), 'fasta' ):
		start = 0
		stop = 0
		exons = map(lambda x: x.split('-'), record.description.split(', ')[1].split(' ')[-1].split(','))
		ch = chromosome_conversion(record.description.split('Chr ')[1].split(' ')[0])
		
		#find introns in genes that have them
		introns = []
		if record.id.split('-')[0][-1] == 'C':
			exons = exons[::-1]
			start = int(exons[0][0])
			stop = int(exons[-1][1])
			if len(exons) > 1:
				for e in range(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		else:
			start = int(exons[0][0])
			stop = int(exons[-1][1])
			if len(exons) > 1:
				for e in range(len(exons)-1):
					introns.append([int(exons[e][1])+1, int(exons[e+1][0])-1])
		
		genes[record.id] = [record.id] + [start, stop] + [exons] + [introns] + [len(record.seq)/3.0] + [ch]
		
		
	#populate list of chromosomes in genome	
	genome = {}
	for record in SeqIO.parse( open(genome_file, 'r'), 'fasta' ):
		genome[chromosome_conversion(record.description.split(' [')[4].split('=')[1][0:-1])] = record.seq.tostring()
	
	#populate second dictionary of non-coding annotations
	#noncoding[ID] = [ ID, chrom, regiontype, start, stop ]
	noncoding = {}
	for line in open(noncoding_file, 'r'):
		if line.startswith('#') != True:
			l = line.strip().split('\t')
			noncoding[randint(1, 9999999)] = [ l[8].split(';')[0].split('=')[1], chromosome_conversion(l[0][3:]), l[2], int(l[3]), int(l[4]) ]

	"""ANNOTATE SNPS"""
	
	#open output file
	f_out = open(BED+'.annotated', 'w')
	
	#start reading BED file
	#chr, start, stop, ref, obs
	for line in open(BED, 'r'):
		#if it's a header, then print out more header to the next line
		if 'start' in line:
			print >> f_out, line.strip() + '\t' + '\t'.join(['annotation', 'region', 'protein'])
			continue
		l = line.strip().split('\t')
		
		snp_pos = int(l[1])
		
		#make copy of chromosome with mutation
		mut_chr = list(copy(genome[chromosome_conversion(l[0])]))
		try:
			mut_chr[snp_pos-1] = l[4]
		except IndexError:
			print >> f_out, "\t".join(l + ["POSITION", "OUTOF", "RANGE"])
			continue

		#is mutation an indel
		indel = False
		if len(l[5]) != len(l[6]):
			indel = True

		annotation = False
		
		#loop through genes, trying to find one containing snp
		for g in genes:
			#if gene on correct chromosome
			if genes[g][6] == chromosome_conversion(l[0]):
				#CRICK STRAND
				if genes[g][0].split('-')[0][-1] == 'C':
					#found gene containing snp
					if snp_pos <= genes[g][1] and snp_pos >= genes[g][2]:
						#check if snp is in an intron:
						for intron in genes[g][4]:
							if snp_pos <= intron[0] and snp_pos >= intron[1]:
								#found an intronic snp - check if it's near a splice site
								if snp_pos in range(intron[0]-2, intron[0]+1) or snp_pos in range(intron[1],intron[1]+3):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
								#otherwise it's just intronic (boring!)
								else:
									print >> f_out, '\t'.join(l + ['intron', genes[g][0], 'NA'])
						
						#if the snp is within gene start-stop, but not in an intron, then it's in an exon
						#remove intronic sequences from gene
						mut_gene = []
						wt_gene = []
						
						#or, if it's an indel within an exon, mark it as such
						#but don't do anything else, because this code doesn't
						#do indels

						if indel:
							print >> f_out, '\t'.join(l + ['INDEL-exon', genes[g][0], 'NA'])
							annotation = True
							break

						for exon in genes[g][3]:
							#check first if the snp is near a splice site
							if snp_pos in range(int(exon[0])-2,int(exon[0])+1) or snp_pos in range(int(exon[1]),int(exon[1])+3):
								#if it's not in the start/stop regions (these aren't splice-sites)
								if snp_pos not in range(genes[g][1]-2,genes[g][1]+1) and snp_pos not in range(genes[g][2],genes[g][2]+3):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
								
							mut_gene += list(complement(mut_chr[ int(exon[1])-1:int(exon[0]) ])[::-1])
							wt_gene += list(complement(genome[genes[g][6]][ int(exon[1])-1:int(exon[0]) ])[::-1])
													
						#loop through codons, find mismatch
						for codon in range(0, len(mut_gene), 3):
							if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
								mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
								wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
								if mut_aa != wt_aa:
									print >> f_out, '\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
								else:
									print >> f_out, '\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
					
						annotation = True
						
					#if snp_pos isn't in a gene, check if it's upstream of a gene
					elif snp_pos in range(genes[g][1]+1,genes[g][1]+201):
						print >> f_out, '\t'.join(l + ["5'-upstream", genes[g][0], 'NA'])
						annotation = True	

				#WATSON STRAND
				else:
					#found gene containing snp
					if snp_pos >= genes[g][1] and snp_pos <= genes[g][2]:

						#check if snp is in an intron:
						for intron in genes[g][4]:
							if snp_pos >= intron[0] and snp_pos <= intron[1]:
								#found an intronic snp - check if it's splice site
								if snp_pos in range(intron[0],intron[0]+3) or snp_pos in range(intron[1]-2,intron[1]+1):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
								#otherwise, it's just intronic
								else:
									print >> f_out, '\t'.join(l + ['intron', genes[g][0], 'NA'])
								annotation = True	
						
						#if the snp is within the gene start-stop, but not in an intron, then it's in an exon
						#remove intronic sequences
						mut_gene = []
						wt_gene = []
						
						#or, if it's an indel within an exon, mark it as such
						#but don't do anything else, because this code doesn't
						#do indels

						if indel:
							print >> f_out, '\t'.join(l + ['INDEL-exon', genes[g][0], 'NA'])
							annotation = True
							break

						for exon in genes[g][3]:
							#first, check if in a splice-site
							if snp_pos in range(int(exon[0]),int(exon[0])+3) or snp_pos in range(int(exon[1])-2,int(exon[1])+1):
								#if it's not the start/stop
								if snp_pos not in range(genes[g][1],genes[g][1]+3) and snp_pos not in range(genes[g][2]-2,genes[g][2]+1):
									print >> f_out, '\t'.join(l + ['splice-site', genes[g][0], 'NA'])
								
							mut_gene += list(mut_chr[ int(exon[0])-1:int(exon[1]) ])
							wt_gene += list(genome[genes[g][6]][int(exon[0])-1:int(exon[1])])

						#loop through codons, find mismatch
						for codon in range(0, len(mut_gene), 3):
							if mut_gene[codon:codon+3] != wt_gene[codon:codon+3]:
								mut_aa = lookup_codon(''.join(mut_gene[codon:codon+3]))
								wt_aa = lookup_codon(''.join(wt_gene[codon:codon+3]))
								#check if synonymous, non-synonymous
								if mut_aa != wt_aa:
									print >> f_out, '\t'.join(l + ['coding-nonsynonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
								else:
									print >> f_out, '\t'.join(l + ['coding-synonymous', genes[g][0], wt_aa + str(codon/3+1) + mut_aa])
								annotation = True

						annotation = True
											
					#if snp isn't in a gene, check if it's in the upstream region of a gene
					elif snp_pos in range(genes[g][1]-200,genes[g][1]):
						print >> f_out, '\t'.join(l + ["5'-upstream", genes[g][0], 'NA'])
						annotation = True						
	
		#after checking all genes for genic or non-coding snps, check non-coding elements
		if annotation == False:
			#if snp isn't upstream of a gene, check if it's in a non-coding region of the genome
			for n in noncoding:
				#if correct chromosome
				if noncoding[n][1] == chromosome_conversion(l[0]):
					if snp_pos >= noncoding[n][3] and snp_pos <= noncoding[n][4]:
						print >> f_out, '\t'.join(l + [noncoding[n][2], noncoding[n][0], 'NA'])
						annotation = True		
			#if annotation is still false after check non-coding elements, then the snp is just intergenic
			if annotation == False:
				print >> f_out, '\t'.join(l + ['intergenic', 'NA', 'NA'])
									
	"""OTHER USEFUL FUNCTIONS"""

def chromosome_conversion(chrom_number):
	chrom_conv = {'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 
				'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16,
				'Mito':17, 'mitochondrion':17, 'sp1':1, 'sp2':2, 'sp3':3,
				'sp4':4, 'sp5':5, 'sp6':6, 'sp7':7, 'sp8':8, 'sp9':9, 'sp10':10,
				'sp11':11, 'sp12':12, 'sp13':13, 'sp14':14, 'sp15':15, 'sp16':16}
	if chrom_number.startswith('chr'):
		chrom_number = chrom_number[3:]
	
	try:
		if int(chrom_number) in chrom_conv.values():
			return int(chrom_number)
	except ValueError:			
		return chrom_conv[chrom_number]

def complement(base):
	compbase = []
	comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	for i in range(len(base)):
		compbase.append(comp[base[i].upper()])
	return ''.join(compbase)

def lookup_codon(codon):
	lookup = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
             'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
             'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
             'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
             'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
             'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
             'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
             'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
             'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
             'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
             'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
             'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
             'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
             'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
             'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
             'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G' }
	return lookup[codon.lower()]

# translate DNA -> amino acid
def translate_sequence(seq):
	translated_seq = ''
	i = 0
	while i <= len(seq)-3:
		translated_seq += lookup_codon(seq[i:i+3])
		i += 3
	return translated_seq

if __name__ == "__main__":
	from optparse import OptionParser
	import sys
	from Bio import SeqIO
	from copy import copy
	from random import randint
	
	parser = OptionParser()
	parser.add_option('-f', '--input', action = 'store', type = 'string', dest = 'inputfile', help = 'file with mutations')
	parser.add_option('-s', '--sequences', action = 'store', type = 'string', dest = 'sequences', help = 'fasta file of coding sequences')
	parser.add_option('-n', '--non-coding', action = 'store', type = 'string', dest = 'noncoding', help = 'gff file containing non-coding regions')
	parser.add_option('-g', '--genome', action = 'store', type = 'string', dest = 'genome', help = 'fasta file containing genome sequence')
	(option, args) = parser.parse_args()

	main(option.inputfile, option.sequences, option.noncoding, option.genome)
