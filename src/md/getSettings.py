import sys
import getopt

def getSettings(argv):
#
# ---- default values ---- #
	Celltype = 'Gm12878';runnum=8;
	nNode=1;ncpu=14;ptn='mit';simtime=48;lmpsdir='../../../../../lammps/src/';
	ctcfthres=4;step=40E6;
	bind_flxb=100;cap=50;
	chrom_lst=[1];
# ------------------------ #

	try:
		opts,args = getopt.getopt(argv,'hC:n:N:p:i:t:l:d:r:b:a:c:',\
						['Cell=','runnum=','nNode=','ncpu','ptn','time','lmpsdir','ctcfthres','step','bindflex=','cap=','chrom='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: parallel_cmap.py -C <Celltype> -n <run number> -N <number of Node> -p <number of cpu> -i <partition> 
                               -t <simulation time> -l <Lammps dir> -d <near CTCF threshold> -r <simulation steps> -b <binding flexbility> -a <CTCF-cohesin nearest dist> -c <chromosome id> 
          or: parallel_cmap.py --Cell <Celltype> --runnum <run number> --nNode <number of Node> --ncpu <number of cpu> --ptn <partition>
                               --time <simulation time> --lmpsdir <Lammps dir> --ctcfthres <near CTCF threshold> --step <simulation steps> --bindflex <binding flexbility> --cap <CTCF-cohesin nearest dist> --chrom <chromosome id>
''')
				sys.exit()
			elif opt in ('-C','--Cell'):
				Celltype = arg
			elif opt in ('-n','--runnum'):
				runnum = int(arg)
			elif opt in ('-N','--nNode'):
				nNode = int(arg)
			elif opt in ('-p','--ncpu'):
				ncpu = int(arg)
			elif opt in ('-i','--ptn'):
				ptn = arg
			elif opt in ('-t','--time'):
				simtime = int(arg)
			elif opt in ('-l','--lmpsdir'):
				lmpsdir = arg
			elif opt in ('-d','--ctcfthres'):
				ctcfthres = int(arg)
			elif opt in ('-r','--step'):
				step = int(arg)
			elif opt in ('-b','--bindflex'):
				bind_flxb = arg
			elif opt in ('-a','--cap'):
				cap = arg
			elif opt in ('-c','--chrom'):
				chrom1st	= [int(arg)]
				chrom2te	= map(eval, args)
				chrom_lst	= chrom1st+chrom2te
		return Celltype,runnum,nNode,ncpu,ptn,simtime,lmpsdir,ctcfthres,step,bind_flxb,cap,chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options amd intializing the simulation!
   > To see the manual and change the settings:
     python parallel_cmap.py -h
''')
		sys.exit()
