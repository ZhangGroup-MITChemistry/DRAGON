import sys
sys.path.append('../../../../src/ctcf/')
from Ipt_module import *
from Params import *
Params()

from getSettings import getSettings
from ProcessCTCFSites import *
from GenCTCFinput import *
from writeinfunc import *

if __name__ == '__main__':

	# ---- input settings ---- #
	celltype,\
	bind_flxb,cap,\
		chrom_lst=getSettings(sys.argv[1:])
	
	Pcs = ProcessCTCFSites()
	for chrId in chrom_lst:
		# ---- prepare ctcf motif ---- #
		orilst_lbm = Pcs.extractMotif('lbm',chrId)
		orilst_known = Pcs.extractMotif('known',chrId)
		orilst_disc = Pcs.extractMotif('disc',chrId)

		try:
			# ---- process CTCF site decide with near cohesin ---- #
			# ---- and orientation decided according to motif ---- #
			final_ctcf_states = Pcs.processingCTCFori(celltype,chrId,orilst_lbm,\
												orilst_known,orilst_disc,bind_flxb,cap)

			# ---- output the list of ctcf sites ---- #
			writein_ctcf(celltype,chrId,cap,final_ctcf_states,chr_region)

			# ---- generate the complete list of ctcf index as input to the model ---- #
			GCInput = GenCTCFinput()
			GCInput.generate(celltype,chrId,cap)

			print('''   > CTCF-binding sites for %s, chromosome %d is successfully constructed.'''\
							%(celltype,chrId))
		except IOError:
			print('''
>>>> [WARNING]: Error in calculating CTCF-binding sites of %s, chromosome %d!
     Please recheck the README file for detail.'''%celltype,chrId)
	print