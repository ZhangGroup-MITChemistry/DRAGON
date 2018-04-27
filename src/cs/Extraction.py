# convert the ChromHMM output
from Ipt_module import *
from Params import *
Params()

class Extraction():

	# ----    Global path 	---- #
	_cs_folder = 'OUTPUTSAMPLE_5kb_6celltype_15states'
	_csdir = '%s/../../../../processEpigenomicsData/chromatinStates/%s/'\
														%(glb_path,_cs_folder)
	_csfile = ''
	_todir = ''


	def convert2raw(self,celltype,chrId):
	# ----    raw output from ChromHMM    ---- #
		infile = self._csdir+celltype+'_15_segments.bed'

		#   ----    raw State processed folder  ---- #
		global _todir
		_todir = '%s/%s/'%(glb_path,celltype)
		todirraw = _todir+'rawStates/'
		if os.path.exists(todirraw) is not True:
			os.makedirs(todirraw)

		#   ----    raw State processed file    ---- #
		global _csfile
		_csfile = '%s/%s_chr%d_chromatin_states_raw.txt'%(todirraw,celltype,chrId)

		#   ----    generate raw State file     ---- #     
		fo = open(_csfile, 'w')
		for line in fi.input(infile):
			items = line.split()
			if items[0] == 'chr%d'%chrId:
				gsta = int(items[1])
				gend = int(items[2])
				for gp in range(gsta, gend, resolution):
					fo.write('%d %s\n'%(gp, items[3]))
		fo.close()


	def raw2state(self,celltype,chrId,realPos=False):
	#   ----    Chromatin state file name     ---- #
		staid = chr_region[chrId-1,1]
		endid = staid+25
		csta = staid * Mb +1
		cend = endid * Mb
		outfile = '%s/%s_chr%d_chromatin_states_From%dMbTo%dMb.txt'\
								%(_todir,celltype,chrId,staid,endid)

		#   ----    generate chrom state file     ---- #
		fo = open(outfile, 'w')
		for line in fi.input(_csfile):
			items = line.split()
			gpos = int(items[0])
			#   ----    output the absolute position  ---- #
			if gpos >=csta and gpos <= cend:
				if realPos:
					fo.write('%8d %4d\n'%(gpos,int(items[1][1::])))
				else:
					fo.write('%8d %4d\n'%((gpos-csta)/resolution+1,\
											int(items[1][1::])))
		fo.close()
		print('''   > Chromatin state for %s, chromosome %d is successfully generated.'''\
																		%(celltype,chrId))
