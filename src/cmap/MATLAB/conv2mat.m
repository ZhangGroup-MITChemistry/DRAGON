
function conv2mat(celltype,chrId,resolution)

    fprintf(['Working on processing Hi-C for ',celltype,' chrom ',num2str(chrId),' ...\n\n'])
    
    % initialize the parameters
    chrIdStr = num2str(chrId,'%d');
    resolutionStr = [num2str(resolution/1E3,'%d') 'kb'];
    mainFolder = ['./hic/rawMap/',celltype,'_combined/'];
    dataFolder = [mainFolder resolutionStr ...
            '_resolution_intrachromosomal/chr',chrIdStr '/'];
    
    % normalize
    tmp = load([dataFolder,'/MAPQGE30/chr',chrIdStr,'_',...
                            resolutionStr,'.RAWobserved']);
                        
    % convert genomic position to index
    tmp(:,1:2) = tmp(:,1:2) / resolution + 1;
    rawmat = spconvert(tmp);
    save(['./hic/rawMap/sparse_matrix/',celltype,'_chr',chrIdStr,...
                    '_',resolutionStr,'_sparse.mat'],'rawmat');
                
    % krnorm
    krnorm = load([dataFolder, '/MAPQGE30/chr', chrIdStr, '_',...
                                resolutionStr, '.KRnorm']);
    save(['./hic/rawMap/sparse_matrix/',celltype,'_chr',chrIdStr,...
                    '_', resolutionStr,'_krnorm.mat'],'krnorm');

end