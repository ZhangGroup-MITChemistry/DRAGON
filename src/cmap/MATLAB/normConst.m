
function normConst(celltype,chrId,resolution,gpSta,gpEnd)

    resolution = resolution/1E3;
    resolutionStr = [num2str(resolution,'%d') 'kb'];
    chrIdStr = num2str(chrId, '%d');    
    load(['./hic/hicMat/',celltype,'_chr',chrIdStr,'_',resolutionStr,'_',num2str(gpSta,'%d'),...
                                            '_',num2str(gpEnd,'%d'),'.mat']);
    hicmat = hicmat + 1e-15;
    
    pii = nonzeros(triu(hicmat,1) - triu(hicmat,2));
    ind = isnan(pii); pii = pii(~ind); 
    nc = mean(pii);
    save(['./hic/normConst/',celltype,'_chr',chrIdStr,'_',resolutionStr,'.txt'],'nc','-ascii');
    
end