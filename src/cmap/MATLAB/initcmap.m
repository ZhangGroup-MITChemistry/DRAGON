
function [simmap,lbl,lbl_tick,celltype,chrId,resolution,nbead,hic_path,hic_norm_path,gpSta,gpEnd] = initcmap(myparam)
    
    % initialize the parameters 
    myparam = initparam(myparam);
    celltype = myparam.celltype;
    chrId = myparam.chrId;
    resolution = myparam.resolution;
    nbead = myparam.nbead;
    simmap_path = myparam.simmap_path;
    hic_path = myparam.hic_path;
    hic_norm_path = myparam.hic_norm_path;
    
    % Hi-C map file name index
    chr_seg = load('../../src/chr_region.txt');
    sepDist = chr_seg(chrId,3)-chr_seg(chrId,2);
    gpSta = (chr_seg(chrId,2)+sepDist*0.08)*1E6/resolution+1;
    gpEnd = (chr_seg(chrId,3)-sepDist*0.12)*1E6/resolution;
    
    % axis indicator
    lbl = [1];
    lbl_tick = cell(1,5);
    lbl_tick{1} = num2str(myparam.gpSta);
    step = myparam.sepDist*0.8/4;
    for x_lbl=2:1:5
        lbl = [lbl;(x_lbl-1)*myparam.nbead/4];
        lbl_tick{x_lbl} = num2str(myparam.gpSta+step*(x_lbl-1));
    end
    
    % load the simulated and hic map
    simmap = load(simmap_path);

    % exclude unreasonable values
    simmap(simmap > 1.0) = 1.0;
    
end