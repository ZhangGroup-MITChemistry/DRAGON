
function [cmap,lbl,lbl_tick,celltype,chrId,resolution] = initcmap(myparam)
    % initialize the parameters 
    myparam = initparam(myparam);
    celltype = myparam.celltype;
    chrId = myparam.chrId;
    resolution = myparam.resolution;
    
    % axis indicator
    lbl = [1];
    lbl_tick = cell(1,5);
    lbl_tick{1} = num2str(myparam.startnum+2);
    for x_lbl=2:1:5
        lbl = [lbl;(x_lbl-1)*1000/(myparam.resolution/5E3)];
        lbl_tick{x_lbl} = num2str((x_lbl-1)*5+myparam.startnum+2);
    end
    
    % load the simulated and hic map
    simmap = load(myparam.simmap_path);
%     hicmap = load(myparam.hic_path);                            % hic map in .txt format
    hicmap = cell2mat(struct2cell(load(myparam.hic_path)));   % hic map in .mat format
    norm_const = load(myparam.hic_norm_path);
    hicmap = hicmap/norm_const;
    
    % exclude unreasonable values
    simmap(simmap > 1.0) = 1.0;
    hicmap(isnan(hicmap))=0;
    
    %%% combine into one map
    % simulated map on the upper diag
    % hic map on the lower diag
    cmap = zeros(myparam.nbead);
    for ii=1:1:myparam.nbead
        for jj=ii+1:1:myparam.nbead
            cmap(ii,jj) = simmap(ii,jj);
        end
        for kk=ii:1:myparam.nbead
            cmap(kk,ii) = hicmap(kk,ii);
        end
    end
    
    % logarithmize the map
    cmap = cmap+1e-12;
    cmap = log(cmap);
    
end