function [v,opt_value,z_out,r_vec,r_tot,status,err]=CVX_PA(v_t,r_t,z_t,sig,gamma,H,iter,tn,t_lv)
% if display_iter==0
% cvx_begin quiet
% else
%   cvx_begin
% end
[m,n]=size(v_t);
cvx_begin quiet
variable t(1,tn)
variable r(1,n)
variable z(1,(n-1))
variable v(m,n) complex
maximize (t(tn))
subject to
r>=zeros(1,n);
t>=zeros(1,tn);
for i=1:length(t_lv)
    if i==1
        for j=1:t_lv(i)
            if (j==t_lv(i))&&(rem(n,2)==1)
                norm([2*t(t_lv(i)),(r(n))])<=(r(n)+2);
            else
               norm([2*t(j),(r(2*j-1)-r(2*j))])<=(r(2*j-1)+r(2*j)+2);
            end
        end
    else
        for j=1:t_lv(i)
            if (j==t_lv(i))&&(rem(t_lv(i-1),2)==1)
                norm([2*t(sum(t_lv(1:i))),(t(sum(t_lv(1:(i-1))))-1)])<=(t(sum(t_lv(1:(i-1))))+1);
            else
                norm([2*t(sum(t_lv(1:(i-1)))+j),(t(sum(t_lv(1:(i-2)))+2*j-1)-t(sum(t_lv(1:(i-2)))+2*j))])<=(t(sum(t_lv(1:(i-2)))+2*j-1)+t(sum(t_lv(1:(i-2)))+2*j));
            end
        end
    end
end
% constraints in form 2
for i=1:n-1
0.25*((r(i)+z(i))^2-(r_t(i)-z_t(i))^2-2*(r_t(i)-z_t(i))*(r(i)-r_t(i)+z_t(i)-z(i)))<=2*real(H(i,:)*v(:,i)*v_t(:,i)'*H(i,:)')-real(H(i,:)*v_t(:,i)*v_t(:,i)'*H(i,:)');
end
%% constraints form 1
for i=1:n-1
    buf=[];
    for j=i+1:n
        buf=[buf,2*H(i,:)*v(:,j)];
    end
    norm([buf,z(i)-sig-1])<=z(i)-sig+1;
end
% constraint on last desired signal (constraint in form 3)
sig*r(n)<=2*real(H(n,:)*v(:,n)*v_t(:,n)'*H(n,:)')-real(H(n,:)*v_t(:,n)*v_t(:,n)'*H(n,:)');
% constraints in form 4-interferences
for i=1:n-1;
    for k=i+1:n;
        buf=[];
        for j=i+1:n
            buf=[buf,2*H(k,:)*v(:,j)];
        end
        norm([buf,(1/gamma*(2*real(H(k,:)*v_t(:,i)*v(:,i)'*H(k,:)')-H(k,:)*v_t(:,i)*v_t(:,i)'*H(k,:)')-sig-1)])<1/gamma*(2*real(H(k,:)*v_t(:,i)*v(:,i)'*H(k,:)')-real(H(k,:)*v_t(:,i)*v_t(:,i)'*H(k,:)'))-sig+1;
    end
end
norm(v(:))<=1;
cvx_end
z_out=z;
r_vec=r(1:n-1);
r_tot=r;
opt_value=t(tn);
err=norm(v_t(:)-v(:));
fprintf('\n Iteration#%d: ',iter)
disp(cvx_status)
fprintf('- optimal value: %f- error=%f',t(tn),err)
status=cvx_status;
