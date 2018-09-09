
package require topotools

topo readlammpsdata /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/lmps/rdm//data.chromosome_40
animate read dcd /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/lmps/condense//DUMP_FILE.dcd waitfor all

set nf [molinfo top get numframes]
puts $nf
animate write pdb /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/pdb//initial_config_40.pdb beg [expr $nf-1]
animate write psf /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/psf//initial_config_40.psf
topo guessangles
