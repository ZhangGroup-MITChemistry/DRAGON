# check running status
from Ipt_module import *
from Params import *
Params()

def checkStatus(usr_name,job_name):
# ---- check the job status ---- #
	njob = -1
	while True:
		cmd = 'squeue |grep %s |grep %s |wc -l'%(usr_name,job_name)
		q = Popen(cmd, shell=True, stdout=PIPE)
		njob = int(q.communicate()[0])
		print("   > Number of jobs calculating cmap on the cluster: %d"%njob)
		if njob == 0:
			print('''
>>>> Start to combine the contact maps ......''')
			break
		else:
			time.sleep(4)

