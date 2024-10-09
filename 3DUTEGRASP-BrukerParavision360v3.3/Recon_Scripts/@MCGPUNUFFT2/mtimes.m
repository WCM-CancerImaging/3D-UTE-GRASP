function res = mtimes(a,bb)
if (a.adjoint)
    [~,nc,nt]=size(bb);
    bb=bb.*repmat(a.w,[1,nc,1]);
    for tt=1:nt
        a.op{tt}.sensChn = size(a.b1,3);
        a.op{tt}.sens=a.b1;
        res(:,:,:,tt) = gpuNUFFT_adj(a.op{tt},bb(:,:,tt));
    end
    res=res/a.normsqr;
else
    [nx,ny,nz,nt]=size(bb);
    for tt=1:nt
        a.op{tt}.sensChn = size(a.b1,3);
        a.op{tt}.sens=a.b1;
        res(:,:,tt) = gpuNUFFT_forw(a.op{tt},bb(:,:,:,tt));
    end
    size(res);
    size(a.w);
    size(a.b1);
    if size(res,1)==1
	res=permute(res,[2,1,3]);   
    end
    res=res.*repmat(a.w,[1,size(a.b1,3),1]);
    res=res/a.normsqr;
end

