
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
    startnum = myparam.startnum;
    endnum = myparam.endnum;
    
    % Hi-C map file name index
    sepDist = endnum-startnum;
    gpSta = (startnum+sepDist*0.08)*1E6/resolution+1;
    gpEnd = (endnum-sepDist*0.12)*1E6/resolution;
    
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