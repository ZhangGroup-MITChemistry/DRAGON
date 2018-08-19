
package require topotools

mol load psf /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/psf//initial_config_40.psf pdb /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/pdb//initial_config_40.pdb
animate read dcd /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/lmps/condense//DUMP_FILE.dcd waitfor all
[atomselect top all] moveby [vecscale -1.0 [measure center [atomselect top all]]]

set sel [atomselect top all]
foreach sl [$sel get serial] {
	set sel [atomselect top "serial $sl"]
	$sel set name CA
	$sel set resname ALA
	$sel set resid $sl
	$sel delete 
}

set nf [molinfo top get numframes]
animate write pdb /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/pdb//initial_config_40_centered.pdb beg [expr $nf-1]

#   ---
mol delete all
mol load psf /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/psf//initial_config_40.psf pdb /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//initial_config/pdb//initial_config_40_centered.pdb
topo guessangles
topo writelammpsdata /nobackup1b/users/qiyf/project_3d_genome/swf-Dragon-extension/runMolecularDynamics/../src/md/lmps_input//data.chromosome.init40
