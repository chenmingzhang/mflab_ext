function H=getWaterTable(H)
    % obtain water table from H struct (generated from readDat)
    % this function is of important when the aquifer is not within one layer
  for n=1:size(H,1)
    nx=size(H(n).values,1);
    ny=size(H(n).values,2);
    nz=size(H(n).values,3);
    H(n).nonnanlay=zeros(nx,ny);
    H(n).HUnConf=zeros(nx,ny);
    for i=1:nx  
      for j=1:ny
           
         aa=find( ~isnan(H(n).values(i,j,:))>0,1);
         
         %if j==28
         %    kk=1
         %end
	 if isempty(aa)
	   H(n).nonnanlay(i,j)=nan;
	   H(n).HUnConf(i,j)=nan; 
	 else
           H(n).nonnanlay(i,j)=aa;
	   H(n).HUnConf(i,j)=H(n).values(i,j,H(n).nonnanlay(i,j) );
	 end

%         fprintf(1,[num2str(i),num2str(j),'\n'])
      end  % i loop
    end   % j loop
%     fprintf(1,[n,'\n'])
  end  % n loop
