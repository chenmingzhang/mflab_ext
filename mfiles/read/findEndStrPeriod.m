function o=findEndStrPeriod(B)
   % o=findEndStrPeriod(B)
   % this function finds the sequence of the output that refers to the end of 
   % the stessperiod
   % for  example:
   % if o(1)=15, then B(15) is the end of the first stress period
    tmp=[B.period];
    o=zeros(1,tmp(end));
    for i=1: length(o)
       sp=find([B.period]==i);
       if ~isempty(sp)
           o(i)=sp(end);
       end
    end
