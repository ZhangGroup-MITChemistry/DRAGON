from Ipt_module import *

def prepMotif(file_name,option,chrId):
#	----	This function is to pre-process the motif file 	---- #

	oriList = []

	for lines in fi.input(file_name):
		every_line = lines.split()
		temp = []
		if option == 'lbm' and every_line[1] == '%d'%chrId:
			temp.append(every_line[2])		# starting position
			temp.append(every_line[4])		# orientation
			oriList.append(temp)			# with orientations as '+' or '-'
		elif option in ['known','disc'] and every_line[1] == 'chr%d'%chrId:
			temp.append(every_line[2])		# starting position
			temp.append(every_line[4])		# orientation
			oriList.append(temp)			# with orientations as '+' or '-'

	return oriList


def prepNarrowPeak(raw_mat,chrId):
#	----	This function is to pre-process the narrow peak file 	---- #
	
	return region = [[peak[1],peak[2]] for peak in raw_mat \
							if peak[0] == 'chr%d'%(chrId)]
