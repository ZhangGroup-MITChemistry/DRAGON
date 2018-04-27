from Ipt_module import *

def checkcell(celltype0):

	cellList = [\
'Gm12878',\
'H1hesc',\
'Helas3',\
'Hepg2',\
'Huvec',\
'K562']
	flag = 1
	celltype = celltype0
	if celltype0 not in cellList:
		flag = 0
		print('''>>>> [WARNING]: Please recheck the README about the available cell types.''')
	elif celltype0 == 'Helas3':
		celltype = 'Hela'

	return flag,celltype

