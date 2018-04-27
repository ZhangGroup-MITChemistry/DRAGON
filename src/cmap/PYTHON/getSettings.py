import sys
import getopt

def getSettings(argv):

# ---- default values ---- #
	Celltype = 'Gm12878';runnum=8;
	jobname='cmap';usrname='user';ptn='mit';
	chrom_lst=[1];
# ------------------------ #

	try:
		opts,args = getopt.getopt(argv,'hC:n:j:u:i:c:',\
								['Cell=','runnum=','job','user','ptn','chrom='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: calContactMap.py -C <Celltype> -n <run number> -j <jobname> -u <username> -i <partition> -c <chromosome id>
          or: calContactMap.py --Cell <Celltype> --runnum <run number> --job <jobname> --user <username> --ptn <partition> --chrom <chromosome id>
''')
				sys.exit()
			elif opt in ('-C','--Cell'):
				Celltype = arg
			elif opt in ('-n','--runnum'):
				runnum = int(arg)
			elif opt in ('-j','--job'):
				jobname = arg
			elif opt in ('-u','--user'):
				usrname = arg
			elif opt in ('-i','--ptn'):
				ptn = arg
			elif opt in ('-c','--chrom'):
				chrom1st	= [int(arg)]
				chrom2te	= map(eval, args)
				chrom_lst	= chrom1st+chrom2te
		return Celltype,runnum,jobname,usrname,ptn,chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options.
   > To see the manual and change the settings:
     python calContactMap.py -h
''')
		sys.exit()
