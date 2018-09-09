# combine to generate ensemble average contact map
from Ipt_module import *
from Params import *

def combineMaps(celltype,njob,chrId):
# ---- combine the contact map of individual parallel running ---- #
	nf_tot = 0
	gSta = chr_region[str(chrId)][0]
	gEnd = chr_region[str(chrId)][1]
	sepDist = gEnd-gSta
	nbead = int(sepDist*Mb/resolution)
	nbead_cg = int(nbead/cg_fac*0.8);
	comb_map = np.zeros([nbead_cg,nbead_cg])

	cmap_path = '%s/%s/chr%d/'%(glb_path,celltype,chrId)

	if os.path.exists('%s/run%02d/contact_map_CG.txt'%(cmap_path,njob-1)):
		for runid in range(0,njob):
			in_temp = np.loadtxt('%s/run%02d/contact_map_CG.txt'%(cmap_path,runid))
			nf = np.loadtxt('%s/run%02d/nframes.txt'%(cmap_path,runid))

			comb_map += in_temp*nf
			nf_tot += nf

		comb_map /= nf_tot
		np.savetxt('./cmap/contact_map_CG_comb_%s_chr%d.txt'\
					%(celltype,chrId),comb_map,fmt = '%.8f')
		print('   > Contact_map_CG_%s_chrom%d is processed to be combined.'%(celltype,chrId))
	
	else:
		print('''
>>>> [Warning] Error in combining the contact maps of %s, chromosome %d!
               Individual contact maps are not correctly calculated.'''%(celltype,chrId))

