from Ipt_module import *

def writein_2d(to_path,to_path_file,writeinfile):

#	----	templating the writing of a 2d matrix	---- #
	if not os.path.exists(to_path):
		os.makedirs(to_path)
	fw = open(to_path+to_path_file,'w')
	for ii in xrange(len(writeinfile)):
		fw.writelines(str(int(writeinfile[ii][0]))+'\t'\
						+str(int(writeinfile[ii][1]))+'\n')
	fw.close()


def writein_ctcf(celltype,chrId,cap,writeinfile,chr_region):
	
#	----	write in the matrix in the original form	---- #
	gSta = chr_region[chrId-1,1]
	gEnd = gSta+25
	to_path = './%s/rawCTCF.%dbp/'%(celltype,cap)
	to_path_file = '%s_chr%d_ctcf_%dMbTo%dMb.txt'%(celltype,chrId,gSta,gEnd)
	writein_2d(to_path,to_path_file,writeinfile)


def writein_motif(chrId,opt,writeinfile):

#	----	write in the motif orientation file		---- #
	to_path = './motif_%s/'%(glb_path,opt)
	to_path_file = 'motif_chr%d.txt'%chrId
	writein_2d(to_path,to_path_file,writeinfile)
	

def writein_narrowpeak(celltype,tf,chrId,writeinfile):

#	----	write in the motif orientation file		---- #
	to_path = './%s/%s/'%(celltype,tf)
	to_path_file = 'chip-seq_peak_%d.txt'%chrId
	writein_2d(to_path,to_path_file,writeinfile)