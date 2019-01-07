
# NOTE: some of this code has been cleaned for readability purposes. 
# original code was run on 30 Mb windows to speed up the computation and paralellize
# contact maxime.rotival@pasteur.fr if experiencing trouble running the code.

# define the scope of junctions that will be considered in the analysis
TotalCoverage=0 # this feature is legacy and was used to remove low coverage junction. 
L=100000 # max intron length to consider a junction

import gzip
import re
import sys                        
import time
import numpy as np
import pandas as pd

from StringIO import StringIO
from Bio import SeqIO, Phylo, AlignIO
from treetime import TreeAnc
from dendropy import Tree, TaxonNamespace


evo_immuno_pop='/Volumes/evo_immuno_pop'
home='/Volumes/@home'

# define the list of Samples:
SampleFile=evo_immuno_pop+'/Maxime/Splicing/HISAT2/SampleFile.txt'
AllSamples=pd.read_csv(SampleFile,sep='\t',names=['samples','batch','merge'])['samples']

# extract sequence for chromosome chr 
def read_seq_chrom(chr):
	print 'reading sequence chr'+str(chr)
	chrSeq_file=open(evo_immuno_pop+'/Maxime/humanGenome_b37/chr'+str(chr)+'.fa')
	chr_seq=SeqIO.parse(chrSeq_file, 'fasta').next().seq
	return chr_seq

# extract GerpRS for chromosome chr from analysis_start to analysis_end (need total chr length as input). the rest is set to 0.

def read_Gerp_chr (chr,chr_len,analysis_start,analysis_end):
	print 'reading Gerp++ chr'+str(chr)
	Gerp_path=home+'/Annotation/Conservation/hg19.GERP_scores/chr'+str(chr)+'.maf.rates'
	GerpFile=open(Gerp_path,'r') # open new GerpFile
	GerpScore=[0]*chr_len # 0 used for missing data
	Gerp_pos=-1
	for Gerp_line in GerpFile :
		Gerp_pos=Gerp_pos+1
		if Gerp_pos<analysis_start-1:
					continue
		if Gerp_pos%1000000 == 0:
			print str(float(Gerp_pos)/1000000.)+'Mb'
		if Gerp_pos>analysis_end:
			break 
		try:
			score=Gerp_line.strip().split('\t')[1];
			GerpScore[Gerp_pos]=float(score)
		except:
			print 'error during Gerp reading'
			print Gerp_line
			break
	GerpFile.close()
	return GerpScore

chr_seq=read_seq_chrom(cur_chr)
chr_seq_rev=chr_seq.reverse_complement()
chr_len=len(chr_seq_rev)
analysis_start=0
analysis_end=chr_len
GerpRS_chr=read_Gerp_chr(cur_chr,chr_len,0,chr_len)

## read all junctions from a file
juncFile=evo_immuno_pop+"/Maxime/Splicing/HISAT2/Results/All_junctions_aggregated_bycond.txt"
#names=['junc_id','chrom','start','end','strand','start_id','end_id','Total_count_NS','Total_count_LPS','Total_count_PAM3CSK4','Total_count_R848','Total_count_IAV','Total_count'])
junctions=pd.read_csv(juncFile,sep='\t')
junctions=junctions.loc[junctions['chrom'] == cur_chr]

######################################################
##	annotate conservation at donor/acceptor sites   ##
######################################################

# the next bit of code is legacy and was based on start values at on the last nucleotide of the exon.
# we thus shift the start for compatibility reasons
# TODO: update the code 

junctions['start']=junctions['start']-1
junctions['GerpRS_start']=[GerpRS_chr[int(i)] for i in junctions.start.values]
junctions['GerpRS_start2']=[GerpRS_chr[int(i)+1] for i in junctions.start.values]
junctions['GerpRS_end']=[GerpRS_chr[int(i)-1] for i in junctions.end.values]
junctions['GerpRS_end2']=[GerpRS_chr[int(i)-2] for i in junctions.end.values]

######################################################
##	  annotate sequence at donor/acceptor sites     ##
######################################################

junctions['acceptor_plus']=[''.join([str(chr_seq)[i] for i in range(int(end)-2,int(end))]).upper() for end in junctions.end.values]
junctions['acceptor_minus']=[''.join([str(chr_seq_rev[i]) for i in range(chr_len-int(start)-2,chr_len-int(start))]).upper() for start in junctions.start.values]
junctions['acceptor_site']=np.where(junctions.strand.values=='-',junctions['acceptor_minus'],junctions['acceptor_plus'])
del junctions['acceptor_plus']
del junctions['acceptor_minus']
junctions['donor_plus']=[''.join([str(chr_seq)[i] for i in range(int(start),int(start)+2)]).upper() for start in junctions.start.values]
junctions['donor_minus']=[''.join([str(chr_seq_rev[i]) for i in range(chr_len-int(end),chr_len-int(end)+2)]).upper() for end in junctions.end.values]
junctions['donor_site']=np.where(junctions.strand.values=='-',junctions['donor_minus'],junctions['donor_plus'])
del junctions['donor_plus']
del junctions['donor_minus']

###########################
##	add intropolis info  ##
###########################
intropolis_file=home+'/Annotation/intropolis/ByChrSummary/intropolis.v1.hg19_summary_chr'+str(cur_chr)+'.tsv.gz'
intropolis_chr=pd.read_csv(intropolis_file,sep='\t',names=['chrom','start','end','strand','donor','acceptor','Nbread','NbSample'])
intropolis_chr['chrom']=[x[3:] for x in intropolis_chr['chrom']]
intropolis_chr['start']=intropolis_chr['start']-1

junctions=pd.merge(junctions, intropolis_chr, how='left', on=['chrom','start','end'])

## write junctions to a file
csv_file=evo_immuno_pop+'/Maxime/Splicing/Splicing/HISAT2/Results/All_junctions_aggregated_bycond_annotated_chr'+cur_chr+'.txt'
junctions.to_csv(csv_file,sep='\t')

######################################################
##	 read sequence from Multiple alignment file     ##
######################################################

with gzip.open(home+"/Annotation/Conservation/SEqAlign/chr"+str(cur_chr)+".maf.gz",'rb') as maf_file:
	line=maf_file.readline() # discard header line 
	#	print line
	line=maf_file.readline()
	species={}
	hg19chrSize=0
	posMin=analysis_start
	skip=0
	pos=0
	posMax=analysis_end
	i=0
	has_hg19=False
	while 1:
		i=i+1
		line = maf_file.readline()
		if not line:
			print 'stop: empty line' 
			break
		if pos > (posMax+1e6) :
			print 'posmax reached'
			break
	# 	if i > 5000:
	# 		print 'stop: all lines done' 
	# 		break
		if line[0]=='a':
			score=float(line.strip().split('=')[1])
			skip=0
			print line,
		if line[0]=='s':
			align=re.split(' +', line.strip())
			species_chr=align[1].split('.')
			if hg19chrSize==0 :
				if species_chr[0]=='hg19':
					hg19chrSize=int(align[5])
			if not species_chr[0] in species:
				species[species_chr[0]]=['-']*(analysis_end-analysis_start+1)
			if species_chr[0]=='hg19':
				print align[2],
				has_hg19=True
				pos=int(align[2])
				seq_hg19=align[6]
				if pos < posMin-2e6 :
					print '.'
					skip=1
					continue
				seq_length=int(align[3])
				seqhg19_noVoid=[seq_hg19[i] for i in range(0,len(seq_hg19)) if seq_hg19[i] != '-']
				species['hg19'][(pos-analysis_start):(pos+seq_length-analysis_start)]=seqhg19_noVoid
			else :
				if not has_hg19 or skip == 1:
					continue
				seq_species=align[6]
				if len(seq_species)!=len(seq_hg19):
					print "oups"
					break
				seq_noVoid=[seq_species[i] for i in range(0,len(seq_hg19)) if seq_hg19[i] != '-']
				species[species_chr[0]][(pos-analysis_start):(pos+seq_length-analysis_start)]=seq_noVoid


SpeciesCommon='Human:Chimp:Gorilla:Orangutan:Rhesus:Baboon:Marmoset:Tarsier:lemur:Bushbaby:TreeShrew:Mouse:Rat:Kangaroo_rat:Guinea_Pig:Squirrel:Rabbit:Pika:Alpaca:Dolphin:Cow:Horse:Cat:Dog:Microbat:Megabat:Hedgehog:Shrew:Elephant:Rock_hyrax:Tenrec:Armadillo:Sloth:Opossum:Wallaby:Platypus:Chicken:Zebra_finch:Lizard:X_tropicalis:Tetraodon:Fugu:Stickleback:Medaka:Zebrafish:Lamprey'.split(':')
SpeciesScientific='hg19:panTro2:gorGor1:ponAbe2:rheMac2:papHam1:calJac1:tarSyr1:micMur1:otoGar1:tupBel1:mm9:rn4:dipOrd1:cavPor3:speTri1:oryCun2:ochPri2:vicPac1:turTru1:bosTau4:equCab2:felCat3:canFam2:myoLuc1:pteVam1:eriEur1:sorAra1:loxAfr3:proCap1:echTel1:dasNov2:choHof1:monDom5:macEug1:ornAna1:galGal3:taeGut1:anoCar1:xenTro2:tetNig2:fr2:gasAcu1:oryLat2:danRer6:petMar1'.split(':')

def getAge(start,end,species,mytree,margin=50):
	from StringIO import StringIO
	from Bio import Phylo,AlignIO
	from treetime import TreeAnc
	myfasta=''
	for spec in species.keys():
		myfasta=myfasta+'>'+spec+'\nA'+str(''.join(species[spec][(start-margin-analysis_start):(end+margin-analysis_start)]))+'A\n'
	my_aln = AlignIO.read(StringIO(myfasta), 'fasta')
	ta = TreeAnc(tree=mytree, aln=my_aln, gtr='JC69')
	ta.infer_ancestral_sequences(method = 'ml', infer_gtr=True, marginal=False)
	seq=ta.get_reconstructed_alignment()
	print 'ancestral seq inferred'
	# 	for i in range(len(seq)):
	# 		print seq[i].id+' '+seq[i].seq[(margin+1):(margin+1+end-start)]
	Phylo.write(ta.tree, home+'/Annotation/Conservation/SEqAlign/test_'+id+'.nwk', 'newick')
	mytreeAnnot=Phylo.read(home+'/Annotation/Conservation/SEqAlign/test_'+id+'.nwk', "newick")
	AncestralSeq={}
	for y in seq:
		AncestralSeq[y.id]=str(y.seq[(margin+1):(margin+1+end-start)])
	    AppearedInCommonAncestor=['human','chimp','gorilla','orangutan','macaque','marmoset','tarsier','lemur','treeShrew','mouse','cow','elephant','opossum','platypus','chicken','frog','zebrafish','lamprey']
    	AppearedInCommonAncestorScientific=['hg19','panTro2','gorGor1','ponAbe2','rheMac2','calJac1','tarSyr1','micMur1','tupBel1','mm9','bosTau4','loxAfr3','monDom5','ornAna1','galGal3','xenTro2','danRer6','petMar1']
    	CommonAncestor=['hg19','NODE_0000044','NODE_0000043','NODE_0000042','NODE_0000040','NODE_0000039','NODE_0000038','NODE_0000036','NODE_0000035','NODE_0000028','NODE_0000018','NODE_0000013','NODE_0000011','NODE_0000010','NODE_0000007','NODE_0000006','NODE_0000001','NODE_0000000']
	    CommonAncestorName=['Humans','Chimpanzees','Gorillas','Orangutans','Gibbons','Monkeys','Tarsiers','Lemurs','Shrews','Rodents','Boroeutheria','Afroeutheria','Marsupials','Monotremes','Sauropsids','Amphibians','Vertebrates','Chordates']
	    PhyloGeneticDistance=[ mytreeAnnot.distance('hg19', CommonAncestor[i]) for i in range(0,len(CommonAncestor))]
    	SeqTime=[ AncestralSeq[i] for i in CommonAncestor]
    	i=0
#	print SeqTime
	while i< len(SeqTime):
		print i, SeqTime[i]
		if SeqTime[i]==SeqTime[0] :
			i=i+1
		else:
			break
	ageSpliceSite1=i-1
	print "res1:", ageSpliceSite1
	while i< len(SeqTime):
		print i, SeqTime[i]
		if SeqTime[i]==SeqTime[0] or SeqTime[i]=='--' :
			i=i+1
		else:
			break
	ageSpliceSite2=i-1
	print "res2:", ageSpliceSite2
	ageSpliceSite3=max([i for i in range(0,len(SeqTime)) if SeqTime[i]==SeqTime[0]])
	print "res3:",ageSpliceSite3
	
	# age 1: first occurence of the current splice site in the ancestral sequences (not allowing missing data or loss and recovery of function)
	# age 2: first occurence of the current splice site in the ancestral sequences (allowing missing data)
	# age 3: first occurence of a splice site in the ancestral sequences (allowing loss and recovery of function, and missing data) # Final definition used in the manuscript)
	 
	all_ancetral_seq=[AppearedInCommonAncestor[i]+'_'+AncestralSeq[CommonAncestor[i]] for i in range(0,len(CommonAncestor))]
	all_species_seq=[SpeciesCommon[i]+'_'+AncestralSeq[SpeciesScientific[i]] for i in range(0,len(SpeciesCommon))]
	print all_species_seq,len(CommonAncestorName),len(PhyloGeneticDistance)
	return [str(ageSpliceSite1),str(ageSpliceSite2),str(ageSpliceSite3),CommonAncestorName[ageSpliceSite1],str(PhyloGeneticDistance[ageSpliceSite1]),':'.join(all_ancetral_seq),':'.join(all_species_seq)]

mytree = Phylo.read(home+'/Annotation/Conservation/SEqAlign/tree46vertebrates.nwk', "newick")

####################
JunctionFileName=evo_immuno_pop+'/Maxime/Splicing/Splicing/HISAT2/Results/All_junctions_aggregated_bycond_annotated_chr'+cur_chr+'.txt'
JunctionFile=pd.read_csv(JunctionFileName,sep='\t')

SpliceSiteStart=JunctionFile.pivot_table(index=['chrom','start','strand'],values=['Total_count'], aggfunc=sum).fillna(0)
print len(SpliceSiteStart)
mychr=[SpliceSiteStart.index.values[i][0] for i in range(0,len(SpliceSiteStart))]
mystart=[int(SpliceSiteStart.index.values[i][1]) for i in range(0,len(SpliceSiteStart))]
mystrand=[SpliceSiteStart.index.values[i][2] for i in range(0,len(SpliceSiteStart))]

tic=time.time()
with open(evo_immuno_pop+'/Maxime/Splicing/Conservation/SpliceSiteStart_conserv_'+str(cur_chr)+'.txt','w') as spliceSiteFile :
	for i in range(0,len(SpliceSiteStart)):
		print i, str(mystart[i]),str(mystrand[i]),str(myGerpRS[i]),
		try:
			ageData=getAge(mystart[i],mystart[i]+2,species,mytree)
			spliceSiteFile.write('\t'.join([str(cur_chr),str(mystart[i]),str(mystrand[i]),'.','.']),'\t'+'\t'.join(ageData)+'\n')
		except:
			continue
		if i>100 and test:
			break

toc=time.time()
print toc-tic,'s'

tic=time.time()
SpliceSiteEnd=JunctionFile.pivot_table(index=['chrom','end','strand'],values=['Total_count'], aggfunc=sum).fillna(0)
print len(SpliceSiteEnd)
mychr=[SpliceSiteEnd.index.values[i][0] for i in range(0,len(SpliceSiteEnd))]
myend=[int(SpliceSiteEnd.index.values[i][1]) for i in range(0,len(SpliceSiteEnd))]
mystrand=[SpliceSiteEnd.index.values[i][2] for i in range(0,len(SpliceSiteEnd))]
 
with open(evo_immuno_pop+'/Maxime/Splicing/Conservation/SpliceSiteEnd_conserv_'+str(cur_chr)+'.txt','w') as spliceSiteFile :
	for i in range(0,len(SpliceSiteEnd)):
		print i, str(myend[i]),str(mystrand[i]),
		try:
			ageData=getAge(myend[i]-2,myend[i],species,mytree)
			spliceSiteFile.write('\t'.join([str(cur_chr),str(myend[i]),str(mystrand[i]),'.','.'])+'\t'+'\t'.join(ageData)+'\n')
		except:
			continue
		if i>100 and test:
			break
toc=time.time()
print toc-tic,'s'

