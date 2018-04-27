import sys
import getopt

def getSettings(argv):
#	----	This function is to obtain the settings for CTCF-binding processing 	---- #

#	---- default values ---- #
	Celltype = 'Gm12878';bind_flxb=100;cap=50;chrom_lst=[1];
#	------------------------ #

	try:
		opts,args = getopt.getopt(argv,'hC:b:a:c:',\
								['Cell=','bindflex=','cap=','chrom='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: genCTCFbinding.py -C <Celltype> -b <binding flexbility> -a <CTCF-cohesin nearest dist> -c <chromosome id> 
          or: genCTCFbinding.py --Cell <Celltype> --bindflex <binding flexbility> --cap <CTCF-cohesin nearest dist> --chrom <chromosome id> 
''')
				sys.exit()
			elif opt in ('-C','--Cell'):
				Celltype = arg
			elif opt in ('-b','--bindflex'):
				bind_flxb = int(arg)
			elif opt in ('-a','--cap'):
				cap = int(arg)
			elif opt in ('-c','--chrom'):
				chrom1st	= [int(arg)]
				chrom2te	= map(eval, args)
				chrom_lst	= chrom1st+chrom2te
		print('''
>>>> Calculating CTCF-binding sites of %s ......'''%Celltype)

		return Celltype,bind_flxb,cap,chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options for calculating the CTCF-binding sites.
   > To see the manual and change the settings:
     python genCTCFbinding.py -h
''')
		sys.exit()


def getMotifSettings(argv):
#	----	This function is to obtain the settings for motif processing 	---- #

#	---- default values ---- #
	motif_fi = './motifs/hg19.motifs.txt';option='lbm';chrom_lst=[1];
#	------------------------ #
	
	try:
		opts,args = getopt.getopt(argv,'hm:p:c:',\
								['motif=','option=','chrom='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: prepareMotif.py -m <motif_name> -p <motif_folder_name_option> -c <chromosome id>
          or: prepareMotif.py -motif <motif_name> -option <motif_folder_name_option> -chrom <chromosome id>
''')
				sys.exit()
			elif opt in ('-m','--motif'):
				motif_fi = arg
			elif opt in ('-p','--option'):
				option = arg
			elif opt in ('-c','--chrom'):
				chrom1st	= [int(arg)]
				chrom2te	= map(eval, args)
				chrom_lst	= chrom1st+chrom2te
		print('''
>>>> Preprocessing the motif file: %s ......'''%motif_fi)

		return motif_fi,option,chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options for processing the motif files.
   > To see the manual and change the settings:
     python prepareMotif.py -h
''')
		sys.exit()


def getNarrowPeakSettings(argv):
#	----	This function is to obtain the settings for narrow peak processing 	---- #

#	---- default values ---- #
	Celltype = 'Gm12878';tf='ctcf';chrom_lst=[1];
#	------------------------ #
	
	try:
		opts,args = getopt.getopt(argv,'hC:t:c:',\
								['Cell=','tf=','chrom='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: processingNarrowPeak.py -C <Celltype> -t <transcriptional factor> -c <chromosome id>
          or: processingNarrowPeak.py --Cell <Celltype> --tf <transcriptional factor> --chrom <chromosome id>
''')
				sys.exit()
			elif opt in ('-C','--Cell'):
				Celltype = arg
			elif opt in ('-t','--tf'):
				tf = arg
			elif opt in ('-c','--chrom'):
				chrom1st	= [int(arg)]
				chrom2te	= map(eval, args)
				chrom_lst	= chrom1st+chrom2te
		print('''
>>>> Processing narrow binding peak of %s, %s ......'''%(tf,Celltype))
		return Celltype,tf,chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options for processing the NarrowPeak files.
   > To see the manual and change the settings:
     python processingNarrowPeak.py -h
''')
		sys.exit()

