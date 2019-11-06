%% Double check the feasiblity of soultion
function [status]=double_check_sol(r,z,v,H,gamma,sigma,eps)
disp('======================= Solution is being checked==================')
status=1;
[m,n]=size(v);
while status
    % constraint on rate
    for i=1:length(r)-1
        intf=sigma;
        for j=i+1:n
            intf=intf+norm(H(i,:)*v(:,j))^2;
        end
        snr=norm(H(i,:)*v(:,i))^2/intf;
        if r(i)<0 || r(i)>(1+eps)*snr
            status=0;
            fprintf('Violation: Rate of user %d is negative or does not match with SNR\n',i);
            break;
        end
    end
    %user n
    if r(n)<0 || r(n)>(1+eps)*norm(H(n,:)*v(:,n))^2/sigma
       status=0;
       fprintf('Violation: Rate of user %d is negative or does not match with SNR\n',n);
       break;
    end 
    if status==1
        fprintf('Check point 1: All rates are positive and matched with SINRs\n')
    end
vp=v(:);
    if norm(vp)<=1
        fprintf('Check point 2: Total power constraint is met - Pt=%f\n',norm(vp)^2)
    else
        disp('Violation:sum_power exceeds available power budget')
    end
for i=1:length(z)
    if z(i)<0
        fprintf('Violation:auxiliary variable z(%d) is negative\n',i)
        status=0;
        break
    end
end
if status==1
    disp('Check point 3: All auxiliary variablales, z, are positive')
end
% chech interference from user i at user k
for i=1:n-1;
    for k=i+1:n;
        intf=sigma;
        for j=i+1:n
            intf=intf+norm(H(k,:)*v(:,j))^2;
        end
        str=norm(H(k,:)*v(:,i))^2/intf;
        if str<(gamma-2*eps)
            fprintf('Violation: Interfernce from user %d at user %d,equals %f, exceeds gamma\n',i,k,str)
            status=0;
        end
    end
end
if status==1
   disp('Check point 4: Interference res. requirment is met at all users') 
end
    status=0; 
disp('===================================================================')
end
