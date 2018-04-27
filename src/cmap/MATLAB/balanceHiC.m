
function balanceHiC(celltype,chrId,resolution)

    resolution = resolution/1E3;
    resolutionStr = [num2str(resolution,'%d') 'kb'];
    chr_seg = load('../../src/chr_region.txt');

    gpSta = (chr_seg(chrId,2)+2) * 1000/resolution + 1;
    gpEnd = (chr_seg(chrId,3)-3) * 1000/resolution;

    extract_balanced_hicmat_loop(celltype,chrId,resolutionStr,gpSta,gpEnd)

end

function extract_balanced_hicmat_loop(celltype,chrId,resolutionStr,gpSta,gpEnd)
    
    chrIdStr = num2str(chrId, '%d');
    load(['./hic/rawMap/sparse_matrix/',celltype,'_chr',chrIdStr,'_',resolutionStr,'_sparse.mat']);
    load(['./hic/rawMap/sparse_matrix/',celltype,'_chr',chrIdStr,'_',resolutionStr,'_krnorm.mat']);
        
    rawmat_seg = rawmat(gpSta:gpEnd, gpSta:gpEnd);
    
    nl = size(rawmat_seg,1);
    [i,j,v] = find(rawmat_seg);
    hicmat = sparse(i,j, v ./ krnorm(gpSta+i-1) ./ krnorm(gpSta+j-1),nl,nl);
    hicmat = hicmat + hicmat';
    for i = 1:nl
        hicmat(i,i) = hicmat(i,i)/2;
    end
    
    save(['./hic/hicMat/',celltype,'_chr',chrIdStr,'_',resolutionStr,'_',num2str(gpSta,'%d'),...
                        '_',num2str(gpEnd, '%d'), '.mat'],'hicmat','-v7.3')
end
