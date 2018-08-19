
function cmap = combcmap(nbead,hic_path,hic_norm_path,simmap)

%     hicmap = load(myparam.hic_path);                            % hic map in .txt format
    hicmap = cell2mat(struct2cell(load(hic_path)));   % hic map in .mat format
    norm_const = load(hic_norm_path);
    hicmap = hicmap/norm_const;
    
    % exclude unreasonable values
    hicmap(isnan(hicmap))=0;
    
    %%% combine into one map
    % simulated map on the upper diag
    % hic map on the lower diag
    cmap = zeros(nbead);
    for ii=1:1:nbead
        for jj=ii+1:1:nbead
            cmap(ii,jj) = simmap(ii,jj);
        end
        for kk=ii:1:nbead
            cmap(kk,ii) = hicmap(kk,ii);
        end
    end
    
    % logarithmize the map
    cmap = cmap+1e-12;
    cmap = log(cmap);
    
end