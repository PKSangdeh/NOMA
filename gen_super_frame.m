% generate_super_frame
clear all
clc
close all
load v_stack
num_data_stream = 1;
num_ant_ap = 1;
num_valid_sc = 52;
const_num_sc = 64;
num_symbol_per_frame = 100;
const_num_frame = 10;
[num_block,M,N]=size(v_stack);
ng=num_valid_sc/num_block;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STF 
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(1, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
%%% LTF
LTF_DAT_GRP = zeros(3, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1];
LTF_DAT_GRP(3,:) = [1 -1 -1 -1 -1 1 -1 -1 -1 1 1 -1 1 1 1 1 1 -1 -1 1 -1 -1 -1 -1 -1 1 -1 1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 1 1 1 -1 1 -1];
%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP=ones(4,num_symbol_per_frame);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % generate freq data of all users
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_data_grp = zeros(N,num_valid_sc, num_symbol_per_frame);
for indx=1:N
    data_freq = zeros(num_valid_sc,num_symbol_per_frame);
    %stf_dat = STF_DAT_GRP(1,:);
    ltf_dat = LTF_DAT_GRP(indx,:);
    % data_freq(:,1) = stf_dat.';
    % data_freq(:,2) = stf_dat.';
    data_freq(:,2*(indx-1)+1) = ltf_dat.';
    data_freq(:,2*(indx-1)+2) = ltf_dat.';
    for jj = 7:num_symbol_per_frame
        data_freq(1:4,jj) = (1/sqrt(2))*(sign(randn(4,1))+1i*sign(randn(4,1)));
        % data_freq(POS_PILOT_SC,jj) = PILOT_DAT_GRP(:,jj);
    end
    freq_data_grp(indx,:,:)=data_freq;
end
 freq_dat1=squeeze(freq_data_grp(1,:,:));   
 freq_dat2=squeeze(freq_data_grp(2,:,:)); 
 freq_dat3=squeeze(freq_data_grp(3,:,:)); 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Superimposition at antenna 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_ant_1=zeros(num_valid_sc, num_symbol_per_frame);
for i=1:num_block
 data_ant_1((ng*(i-1)+1:ng*i),7:end)=squeeze(v_stack(i,1,1))*freq_dat1((ng*(i-1)+1:ng*i),7:end)+squeeze(v_stack(i,1,2))*freq_dat2((ng*(i-1)+1:ng*i),7:end)+squeeze(v_stack(i,1,3))*freq_dat3((ng*(i-1)+1:ng*i),7:end);
end
for i=1:num_block
data_ant_1((ng*(i-1)+1:ng*i),1:2)=freq_dat1((ng*(i-1)+1:ng*i),1:2)*squeeze(v_stack(i,1,1));
data_ant_1((ng*(i-1)+1:ng*i),3:4)=freq_dat2((ng*(i-1)+1:ng*i),3:4)*squeeze(v_stack(i,1,2));
data_ant_1((ng*(i-1)+1:ng*i),5:6)=freq_dat3((ng*(i-1)+1:ng*i),5:6)*squeeze(v_stack(i,1,3));
end
data_ant_1(POS_PILOT_SC,7:end)=PILOT_DAT_GRP(:,7:end);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Superimposition at antenna 2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_ant_2=zeros(num_valid_sc, num_symbol_per_frame);
for i=1:num_block
 data_ant_2((ng*(i-1)+1:ng*i),7:end)=squeeze(v_stack(i,2,1))*freq_dat1((ng*(i-1)+1:ng*i),7:end)+squeeze(v_stack(i,2,2))*freq_dat2((ng*(i-1)+1:ng*i),7:end)+squeeze(v_stack(i,2,3))*freq_dat3((ng*(i-1)+1:ng*i),7:end);
end
for i=1:num_block
data_ant_2((ng*(i-1)+1:ng*i),1:2)=freq_dat1((ng*(i-1)+1:ng*i),1:2)*squeeze(v_stack(i,2,1));
data_ant_2((ng*(i-1)+1:ng*i),3:4)=freq_dat2((ng*(i-1)+1:ng*i),3:4)*squeeze(v_stack(i,2,2));
data_ant_2((ng*(i-1)+1:ng*i),5:6)=freq_dat3((ng*(i-1)+1:ng*i),5:6)*squeeze(v_stack(i,2,3));
end
% %%
data_ant_2(POS_PILOT_SC,7:end)=PILOT_DAT_GRP(:,7:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Convert to time domain
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     time_signal_grp = zeros(2, 80*num_symbol_per_frame);
     % antenna 1
     data_f_frm = zeros(const_num_sc, num_symbol_per_frame);
     data_f_frm(POS_VALID_SC,:) = data_ant_1;
     data_t_frm = ifft(data_f_frm)*sqrt(const_num_sc);
     data_s_frm = [data_t_frm(end-15:end,:); data_t_frm];
     data_s_frm = data_s_frm(:).';
     time_signal_grp(1,:) = data_s_frm;
     data1=data_ant_1(:).';
     % antenna 2
     data_f_frm = zeros(const_num_sc, num_symbol_per_frame);
     data_f_frm(POS_VALID_SC,:) = data_ant_2;
     data2=data_ant_2(:).';
     data_t_frm = ifft(data_f_frm)*sqrt(const_num_sc);
     data_s_frm = [data_t_frm(end-15:end,:); data_t_frm];
     data_s_frm = data_s_frm(:).';
     time_signal_grp(2,:) = data_s_frm;
     time_signal_grp = repmat(time_signal_grp, 1, const_num_frame);
     figure();
     scatter(real(data1),imag(data1),'ko');
     figure();
     scatter(real(data2),imag(data2),'ko');     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for indx=1:M
    fileID= fopen(['C:\Users\pedra\Desktop\NOMA_2X3\signals\super_tx_' num2str(indx) '.dat'], 'wb');       
    if (fileID < 0)
        disp('Error: fail to open files!');
        %pause;
    end
    data_tmp_c = time_signal_grp(indx,:);
    data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
    data_tmp_f = data_tmp_f(:);
    fwrite(fileID,data_tmp_f,'float');
    fclose(fileID);
end