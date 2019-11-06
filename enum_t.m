function [nt,t_level,eo_index]=enum_t(n)
tmp=1;
t_level=[];
eo_index=[];
nt=0;
while tmp
    eo_index=[eo_index rem(n,2)];
 n=n+eo_index(end);
 n=n/2;
 if n==1;
     tmp=0;
 end
 t_level=[t_level n];
end
nt=sum(t_level);
 
