import sys
import re
import getopt
import json
from Params import *

def getSettings(argv):
#
# ---- default values ---- #
	Celltype = 'Gm12878';njob=1;
	nNode=1;ncpu=14;ptn='mit';simtime=48;lmpsdir='%s/../lammps/src/'%glb_path;
	ctcfthres=4;step=40E6;
	bind_flxb=100;cap=50;
	chrom_lst=[1];chr_reg=[[20,45]];
# ------------------------ #

	try:
		opts,args = getopt.getopt(argv,'hC:n:N:p:i:t:l:r:c:g:',\
								['Cell=','njob=','nNode=','ncpu=',\
								'ptn=','time=','lmpsdir=','step=',\
								'chrom=','region='])
		for opt,arg in opts:
			if opt=='-h':
				print('''
>>>> Options: parallel_cmap.py -C <Celltype>
                               -c <chromosome id>
                               -g <chromosome region>
                               -l <Lammps dir>
                               -r <simulation steps>
                               -n <number of job>
                               -N <number of Node>
                               -p <number of cpu>
                               -i <partition>
                               -t <simulation time> 

          or: parallel_cmap.py --Cell <Celltype>
                               --chrom <chromosome id>
                               --region <chromosome region>
                               --lmpsdir <Lammps dir>
                               --step <simulation steps>
                               --njob <number of job>
                               --nNode <number of Node>
                               --ncpu <number of cpu>
                               --ptn <partition>
                               --time <simulation time>
''')
				sys.exit()
			elif opt in ('-C','--Cell'):
				Celltype = arg
			elif opt in ('-c','--chrom'):
				chrom_lst = init(arg,',')
			elif opt in ('-g','--region'):
				chr_reg = initReg(arg,len(chrom_lst))
			elif opt in ('-n','--njob'):
				njob = int(arg)
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
			elif opt in ('-r','--step'):
				step = int(arg)

		writeReg(chrom_lst,chr_reg)
		return 	Celltype,njob,\
				nNode,ncpu,ptn,simtime,lmpsdir,\
				ctcfthres,step,\
				bind_flxb,cap,\
				chrom_lst

	except getopt.GetoptError:
		print('''
>>>> [Warning] Error in setting options amd intializing the simulation!
   > To see the manual and change the settings:
     python parallel_cmap.py -h
''')
		sys.exit()


def init(inpt,cond):
	tmp = re.split(r'%s'%cond,inpt)
	filtLst = filt(tmp,'')
	filtLst = map(int, filtLst)
	return filtLst


def initReg(inpt,lenc):
# 
# ---- This function is to 
	filtLst = init(inpt,'\[|,|\]')
	chr_reg = np.column_stack((filtLst[::2],filtLst[1::2])).tolist()
	leng = len(chr_reg)
	if lenc != leng:
		if lenc > 1 and leng == 1:
			chr_reg *= lenc
		else:
			print('''
>>>> Error in mismatching of chromosome id with genomic regions!
''')
			sys.exit()
	return chr_reg


def filt(inpt,cond):
# 
# ---- This function is to filter

	try:
		condition = lambda t: t != cond
		filtLst = map(int,filter(condition,inpt))
		return filtLst

	except ValueError:
		print('''
>>>> Error in input format!
   > To see the manual and change the settings:
     python parallel_cmap.py -h
''')
		sys.exit()


def writeReg(chrom_lst,chr_reg):
# 
# ---- Write the chrom region to the src folder

	lenc = len(chrom_lst)
	leng = len(chr_reg)
	if lenc != leng:
	 	if lenc > 1 and leng == 1:
	 		chr_reg *= lenc
	 	else:
	 		print('''
>>>> Error in mismatching of chromosome id with genomic regions!
''')
	 		sys.exit()
	chrDict = dict(zip(chrom_lst,chr_reg))
	chrFile = '%s/../src/chr_region.txt'%glb_path
	with open(chrFile,'w') as f:
	    json.dump(chrDict,f)

