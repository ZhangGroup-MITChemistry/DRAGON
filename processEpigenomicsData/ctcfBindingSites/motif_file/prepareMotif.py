import sys
sys.path.append('../../../src/ctcf/')
from Ipt_module import *

from getSettings import getMotifSettings
from preprocessing import prepMotif
from writeinfunc import *

if __name__ == '__main__':
#
#	----	extract the motif information with the start position and orientation	---- #
#
	# ---- input settings ---- #
	motif_fi,option,chrom_lst=getMotifSettings(sys.argv[1:])

	for chrId in chrom_lst:
		try:
			oriList = prepMotif(motif_fi,option,chrId)
			writein_motif(chrId,option,oriList)

		except IOError:
			print('''
>>>> [WARNING]: Error in preprocessing the motif files!
                Please provide a valid raw motif file, see README for detail, or to use the existing motif files to process.''')
	print
