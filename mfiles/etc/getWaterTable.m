function H=getWaterTable(H)
    % H=getWaterTable(H)
    % obtain water table from H struct (generated from readDat)
    % this function is of important when the aquifer is not within one layer
    % H.nonNanlay  -- the layer  index where unconfined surface locates.
    % H.HUnConf    -- the hydraulic head of the unconfined aquifer
    % H.ind_nonNanLay -- the linear index of the unconfined aquifer stored in
    %                    in [ny*nx] format.
    % at the moment H(3).HUnConf == H(3).values( H(3).ind_nonNanLay) 
    % so now it is clear that H.ind_nonNanLay is easier to be used than 
    % H.nonNanLay.
    % 
  fprintf(1,'now get the water Table from H object\n');
  for n=1:size(H,1)
%    nx=size(H(n).values,1);
%    ny=size(H(n).values,2);
%    nz=size(H(n).values,3);
%    H(n).nonNanLay=zeros(nx,ny);
%    H(n).HUnConf=zeros(nx,ny);
%    tmp=isnan(H(n).values);
%
%
%    for i=1:nx  
%      for j=1:ny
%         aa=find( ~isnan(H(n).values(i,j,:))>0,1);
%%	 H(n).mask_UnConf
%	 if isempty(aa)
%	   H(n).nonNanLay(i,j)=nan;
%	   H(n).HUnConf(i,j)=nan; 
%	 else
%           H(n).nonNanLay(i,j)=aa;
%	   H(n).ind_nonNanLay(i,j)=sub2ind(size( H(n).values) ,i,j,aa);
%	   H(n).HUnConf(i,j)=H(n).values(i,j,H(n).nonNanLay(i,j) );
%	 end
%
%      end  % i loop
%    end   % j loop
    H(n).mask_UnConf= ~isnan(H(n).values) ;
    H(n).mask_UnConf=cumsum(H(n).mask_UnConf,3)==1 ;
    H(n).HUnConf=  nansum (bsxfun(@times , H(n).mask_UnConf, H(n).values),3);
    H(n).HUnConf=  nansum (bsxfun(@times , H(n).mask_UnConf, H(n).values),3);
    % one-liner to make the dry column from 0 to nan
    H(n).HUnConf= H(n).values(:,:,end)*0.*H(n).HUnConf;
%         fprintf(1,[num2str(i),num2str(j),'\n']);
%     fprintf(1,[n,'\n']);
  end  % n loop
  fprintf(1,'Water Table from H object finished\n');



   



