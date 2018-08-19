# generate chromatin states from Epigenetic data
import sys
sys.path.append('../../../../src/cs/')
from Ipt_module import *
from Params import *
Params()

from getSettings import getSettings
from Extraction import *


if __name__=='__main__':

    # ---- input settings ---- #
    celltype,chrom_lst=getSettings(sys.argv[1:])
    Extn = Extraction()

    # ---- two-step extraction ---- #
    try:
        for chrId in chrom_lst:
            Extn.convert2raw(celltype, chrId)
            Extn.raw2state(celltype, chrId)
    except IOError:
		print('''
>>>> [Warning] Error in calculating chromatin states of %s!
     Please check the README file for detail.'''%celltype)
    print
    