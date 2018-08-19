
package require topotools

mol load hoomd   /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/xml//rdm_polymer_150.xml
animate read dcd /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/xml//rdm_polymer_150.dcd waitfor all

set nf [molinfo top get numframes]
puts $nf
animate write pdb /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/pdb//initial_config_150.pdb beg [expr $nf-1]
animate write psf /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/psf//initial_config_150.psf
topo guessangles
topo writelammpsdata /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/lmps/rdm//data.chromosome_150
