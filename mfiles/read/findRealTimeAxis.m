function o=findRealTimeAxis(H)
   % o=findRealTimeAxis(H)
   % this function finds realtime axis associated with H
   % the reason of having this function is because H.totim some time
   % is working as expected.
   % the stessperiod
   % for  example:
   % if o(1)=15, then B(15) is the end of the first stress period
    tmp=[H.period];
    time_sp=[H.pertim];
    last_sp=H(end).period
    o=zeros(1,tmp(end));
    for i=1:last_sp
       sp=find(tmp==i);
       if ~isempty(sp)
         if i==1
           o(sp)=time_sp(sp);
	 else
           o(sp)=o(sp(1)-1)+time_sp(sp);
	 end
       end
    end
