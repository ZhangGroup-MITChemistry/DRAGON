
function myparam = initparam(myparam)
    
    Mb = 1E6;
    
    myparam.celltype=input('Enter the cell type (default Gm12878):\n','s');
    if isempty(myparam.celltype)
        myparam.celltype = 'Gm12878';
    end
    
    myparam.chrId=input('Enter the chromosome id (default chromosome 1):\n');
    if isempty(myparam.chrId)
        myparam.chrId = 1;
    end
    
    myparam.resolution=input('Enter the resolution of contact map (default 50kb):\n');
    if isempty(myparam.resolution)
        myparam.resolution = 50000;
    end
    
    resolutionStr = [num2str(myparam.resolution/1E3,'%d') 'kb'];
    
    myparam.startnum = input('Enter the start genomic position in Mb (default 20Mb):\n');
    if isempty(myparam.startnum)
        myparam.startnum = 20;
    end
    
    myparam.endnum = input('Enter the end genomic position in Mb (default 45Mb):\n');
    if isempty(myparam.endnum)
        myparam.endnum = 45;
    end
    
    myparam.sepDist = myparam.endnum-myparam.startnum;
    myparam.gpSta = myparam.startnum+myparam.sepDist*0.08;
    myparam.gpEnd = myparam.endnum-myparam.sepDist*0.12;
    myparam.nbead = (myparam.gpEnd-myparam.gpSta)*Mb/myparam.resolution;
    
    myparam.simmap_path = ['./cmap/contact_map_CG_comb_' myparam.celltype ...
        '_chr' num2str(myparam.chrId) '.txt'];
    
    if strcmp(myparam.celltype,'Gm12878')
        myparam.celltype='GM12878';
    elseif strcmp(myparam.celltype,'Hela')
        myparam.celltype='HeLa';
    end
    
    myparam.hic_path = input('Enter the path where the Hi-C map is located:\n');
    if isempty(myparam.hic_path)
        myparam.hic_path = ['./hic/hicMat/' myparam.celltype '_chr'...
            num2str(myparam.chrId,'%d'),'_' resolutionStr '_',...
            num2str(myparam.gpSta*Mb/myparam.resolution+1), '_', ...
            num2str(myparam.gpEnd*Mb/myparam.resolution),'.mat'];
    end

    myparam.hic_norm_path = input('Enter the path where the normalization constant is located:\n');
    if isempty(myparam.hic_norm_path)
        myparam.hic_norm_path = ['./hic/normConst/' myparam.celltype '_chr' ...
            num2str(myparam.chrId,'%d'),'_' resolutionStr '.txt'];
    end
    
end