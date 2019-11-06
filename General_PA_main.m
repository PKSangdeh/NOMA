%%
clc
clear all
close all
warning('off','all')
run('C:\Users\pedra\Desktop\cvx\cvx\cvx_startup.m')
%% Parameters
m=2;n=3;
convg_err1=10^-2; % for solved status
convg_err2=10^-2; % for inaccurate/solved status
eps=10^-4;        % safty margin for checking feasibility
sig=10^-5;
gamma=7;          % Least required resolution for strong interferences in dB
gamma=10^(gamma/10);
num_block=13;
ng=4;
v_stack=zeros(num_block,m,n);
[tn,t_level,eo_index]=enum_t(n); % enumerate required auxiliary variable for relaxing objective function
%% channel stack over different block
load channel_u1;
load channel_u2;
load channel_u3;
H=zeros(3,2);
st=zeros(num_block,1);
for i=1:num_block
    clc
    fprintf('\n=================== Block %d of %d =======================\n',i,num_block)
    buf1=ch_u1(:,(i-1)*ng+1:i*ng);
    mag=mean(abs(buf1).');
    ang=mean(angle(buf1).');
    buf2=(mag.*exp(1j*ang));
    H(1,:)=buf2;

    buf1=ch_u2(:,(i-1)*ng+1:i*ng);
    mag=mean(abs(buf1).');
    ang=mean(angle(buf1).');
    buf2=(mag.*exp(1j*ang));
    H(2,:)=buf2;
    
    buf1=ch_u3(:,(i-1)*ng+1:i*ng);
    mag=mean(abs(buf1).');
    ang=mean(angle(buf1).');
    buf2=(mag.*exp(1j*ang));
    H(3,:)=buf2;
    tmp=1;
    iter=0;
    disp('=== Step 1 of 3: finding good tangent points')
while tmp
      iter=iter+1;
      [z_t,r_t,v_t]=find_tangs(H,gamma,sig,eps);
      [v,opt_value,z,r,r_tot,status,error]=CVX_PA(v_t,r_t,z_t,sig,gamma,H,iter,tn,t_level);
      if length(status)==6
         if prod((status=='Solved')) % the staus of solution for initial tangent point must be solved (other status in cvx are: failed, infeasible, and inaccurate/solved)
         tmp=0;
         v_t=v;
         z_t=z;
         r_t=r;
         fprintf('\n')
         disp('>>>>>> reached to a good initial point***')
         end
      end     
end
opt_val=[];
err=[];
capc=[];
%% find the solution with iterative algorithm
itt=0;
fprintf('\n')
disp('====Step 2 of 3: Iterative algorithm begins')
tmp2=1;
while tmp2
    tmp=1;
    if itt>=1
        disp('== Retry===')
        itt=0;
    end
while tmp
    itt=itt+1;
    [v,opt_value,z,r,r_tot,status,error]=CVX_PA(v_t,r_t,z_t,sig,gamma,H,itt,tn,t_level);
    v_t=v;
    z_t=z;
    r_t=r;
    capc=[capc log2(prod(r_tot))];
    err=[err,error];
    opt_val=[opt_val,opt_value];
    if length(status)==6
       if prod((status=='Failed')) 
          tmp=0;
          clc
          disp('=== Reached to a Failure point ===')
          tmp=0;
       elseif prod((status=='Solved')) && error<convg_err1
              fprintf('\n')
              disp('=== Reached convergence point ===')
              v
              st(i)=1;
              fprintf('Optimal sum-capacity: %f\n',capc(end))
              tmp=0; tmp2=0;
       end
    end
    if (error<convg_err2) && (length(status)~=6)
        fprintf('\n')
        disp('=== Reached convergence point ===')
        v
        st(i)=1;
        fprintf('Optimal sum-capacity: %f \n',capc(end))
        tmp=0;tmp2=0;
    end
end
end
v_stack(i,:,:)=v;
end
clc
disp('================= Status of precoding stack ====================')
fprintf('\n')
for i=1:num_block
    if st(i)==1;
        fprintf('Block %2d: found\n',i)
    else
        fprintf('Block %2d: failed\n',i)
    end
end
disp('=================================================================')
save('C:\Users\pedra\Desktop\NOMA_2X3\v_stack.mat','v_stack');