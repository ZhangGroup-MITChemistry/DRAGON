from Ipt_module import *
from Params import *
Params()

from writeinfunc import *

class GenCTCFinput():

	def convert2sq(self,ctcf_path):
	#	----	processing into complete sequence	---- #
		tmp = np.loadtxt(ctcf_path)
		mat = np.zeros((nbead,2),dtype='int')
		for i in range(len(mat)):
			mat[i,0] = i+1
			mat[i,1] = 3
		
		for j in range(len(tmp)):
			mat[int(tmp[j][0])-1][1] = int(tmp[j][1])

		mat = mat.tolist()
		for i in range(len(mat)):
			mat[i][0] = int(mat[i][0])
			mat[i][1] = int(mat[i][1])
			if mat[i][1] == 0:
				mat[i][1] = 4
		return mat

	def extractCtcfConv(self,ctcfSeq):
	#	----	processing into CTCF index sequence	---- #
		ctcfL = []
		ctcfR = []
		ctcfInd = np.zeros((nbead,2),dtype='int')

		for ictcf in ctcfSeq:
			aid = int(ictcf[0])
			cs = int(ictcf[1])
			if (cs == 1):
				ctcfL.append(aid)
			elif (cs == 2):
				ctcfR.append(aid)
			elif (cs == 4):
				ctcfL.append(aid)
				ctcfR.append(aid)
			else:
				pass

		for aid in range(1, nbead+1):
			tmp = -1
			for cid in ctcfL:
				if cid <= aid:
					if cid > tmp:
						tmp = cid
				else:
					break
			if tmp == -1:
				ctcfInd[aid-1,0] = tmp
			else:
				ctcfInd[aid-1,0] = ctcfL.index(tmp)
			tmp = nbead+1
			for cid in ctcfR:
				if cid >= aid:
					if cid < tmp:
						tmp = cid
			if tmp <= nbead:
				ctcfInd[aid-1,1] = ctcfR.index(tmp)
			else:
				ctcfInd[aid-1,1] = -1

		return ctcfInd


	def generate(self,celltype,chrId,cap):
	#	----	generate the input data format for the model 	---- #
	#	----	common bead:3	CTCF+:1	CTCF-:2	CTCF+-:4 		---- #
		gSta = chr_region[chrId-1,1]
		gEnd = gSta+25

		#	----	path to raw CTCF list
		to_path = '%s/%s/rawCTCF.%dbp/'%(glb_path,celltype,cap)
		ctcf_path = '%s/%s_chr%d_ctcf_%dMbTo%dMb.txt'\
											%(to_path,celltype,chrId,gSta,gEnd)

		#	----	path to input CTCF list (index)
		gen_path = '%s/%s/'%(glb_path,celltype)
		gen_path_pos = '%s_chr%d_ctcf_position_From%dMbTo%dMb.txt'\
												%(celltype,chrId,gSta,gEnd)
		gen_path_idx = '%s_chr%d_ctcf_index_From%dMbTo%dMb.txt'\
												%(celltype,chrId,gSta,gEnd)

		ctcfSeq = self.convert2sq(ctcf_path)
		ctcfInd = self.extractCtcfConv(ctcfSeq)

		#	----	write to path for output	---- #
		writein_2d(gen_path,gen_path_pos,ctcfSeq)
		writein_2d(gen_path,gen_path_idx,ctcfInd)
