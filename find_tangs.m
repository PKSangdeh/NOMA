function [z_t,r_t,v_t]=find_tangs(H,gamma,sig,eps)
[n,m]=size(H);
% H=zeros(n,m);
max_iter=300;
% for i=1:n
%     H(i,:)=10/10^(m-1)*(randn(1,m)+1j*rand(1,m));
% end
z_t=zeros(1,n-1);
r_t=zeros(1,n-1);
v_t=zeros(m,n);
it=0;
%% %%%%%%%%%%%%%%%%% find tangent point %%%%%%%%%%%%%%%%%%%%%%%
flag=1;
it=0;
while flag
      tmp=1;
      while tmp
          if n<=5
            v_t(:,n)=1/(10^(n/2))*(randn(m,1)+1j*randn(m,1));
          elseif n>5
            v_t(:,n)=1/(10^(n/4))*(randn(m,1)+1j*randn(m,1));
          end
            if norm(v_t(:,n)*H(n,:))>=sqrt(gamma*sig);
               tmp=0;
            end
      end
      fprintf('v_t(:,%d) is found\n',n)
      flag3=1;
      for i=n-1:-1:1
          flag2=1;
          iter=0;
          while (flag2 && flag3)
                iter=iter+1;
                v_t(:,i)=(1/(10^(i/2)))*(randn(m,1)+1j*randn(m,1));
                if norm(v_t(:))^2<=(1-(i-1)/n)
                   c1=1;
                   for k=i+1:n
                       buf1=0;
                       for j=i+1:n
                           buf1=buf1+norm(H(k,:)*v_t(:,j))^2;
                       end
                       buf1=buf1+sig;
                       if norm(H(k,:)*v_t(:,i))^2<(gamma-eps)*buf1
                           c1=0;
                       end
                   end
                   if c1==1
                       flag2=0;
                   end
                end
                if iter>=max_iter*(n-i)
                    flag3=0;
                end
          end
          if flag3==1
             fprintf('v_t(:,%d) is found\n',i)
          end
      end
      if flag3==1
         flag=0;
      end 
end
%% double check feasibility
% double check feasibility
status=1;
for i=1:n-1;
    for k=i+1:n;
        intf=sig;
        for j=i+1:n
            intf=intf+norm(H(k,:)*v_t(:,j))^2;
        end
        str=norm(H(k,:)*v_t(:,i))^2/intf;
        if str<(gamma-2*eps)
            fprintf('Violation: Interfernce from user %d at user %d,equals %f, exceeds gamma\n',i,k,str)
            status=0;
        end
    end
end
if status==1 && norm(v_t(:))^2>1
    status=0;
end
if status==1
    for i=1:n-1
        intf=sig;
        for j=i+1:n
            intf=intf+norm(H(i,:)*v_t(:,j))^2;
        end
        z_t(i)=intf-eps;
        r_t(i)=(1-eps)*norm(H(i,:)*v_t(:,i))^2/(z_t(i)+eps);
    end
end 

