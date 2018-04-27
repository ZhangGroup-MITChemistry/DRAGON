import sys
import getopt

def getSettings(argv):

# ---- default values ---- #
	Celltype = 'Gm12878';chrom_lst=[1];
# ------------------------ #
	
	try:
		opts,args = getopt.getopt(argv,'hC:c:',\
								['Cell=','chrom='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: genChromState.py -C <Celltype> -c <chromosome id>
          or: genChromState.py --Cell <Celltype> --chrom <chromosome id>
''')
				sys.exit()
			elif opt in ('-C','--Cell'):
				Celltype = arg
			elif opt in ('-c','--chrom'):
				chrom1st	= [int(arg)]
				chrom2te	= map(eval, args)
				chrom_lst	= chrom1st+chrom2te
		print('''
>>>> Calculating chromatin states of %s ......'''%Celltype)
		return Celltype,chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options.
   > To see the manual and change the settings:
     python genChromState.py -h
''')
		sys.exit()

