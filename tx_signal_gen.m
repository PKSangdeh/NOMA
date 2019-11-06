clear;
clc;
close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_data_stream = 1;
num_ant_ap = 1;
num_valid_sc = 52;
const_num_sc = 64;
num_symbol_per_frame = 20;
const_num_frame = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STF 
STF_POS = [3  7  11  15  19  23  30  34  38  42  46  50];
STF_DAT_GRP = zeros(10, 52);    
STF_DAT_GRP(1,STF_POS) = sqrt(2)*[-1-1i  -1+1i  -1-1i  -1-1i  -1-1i  1-1i  -1-1i  -1-1i  1+1i  1-1i  1+1i  -1+1i];
STF_DAT_GRP(2,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1-1i  1-1i  -1-1i  1-1i  -1+1i  1+1i  1-1i  -1+1i  1-1i  -1+1i];
STF_DAT_GRP(3,STF_POS) = sqrt(2)*[-1+1i  1+1i  -1+1i  1+1i  1-1i  1+1i  -1+1i  -1+1i  1-1i  1-1i  -1+1i  -1+1i];
STF_DAT_GRP(4,STF_POS) = sqrt(2)*[1+1i  -1+1i  -1+1i  -1+1i  -1-1i  1+1i  1+1i  1-1i  1-1i  -1+1i  -1+1i  -1-1i];
STF_DAT_GRP(5,STF_POS) = sqrt(2)*[-1+1i  -1+1i  1+1i  -1+1i  -1-1i  1+1i  1-1i  -1-1i  -1-1i  -1+1i  -1+1i  1-1i];
STF_DAT_GRP(6,STF_POS) = sqrt(2)*[1+1i  1-1i  -1-1i  1+1i  -1+1i  1+1i  1-1i  -1-1i  1+1i  -1-1i  -1+1i  -1+1i];
STF_DAT_GRP(7,STF_POS) = sqrt(2)*[1+1i  1+1i  1+1i  1-1i  -1-1i  1-1i  -1-1i  1-1i  1+1i  1+1i  -1-1i  1-1i];
STF_DAT_GRP(8,STF_POS) = sqrt(2)*[-1+1i  1-1i  1+1i  -1+1i  1+1i  1-1i  -1-1i  -1-1i  -1+1i  1-1i  -1+1i  -1-1i];
STF_DAT_GRP(9,STF_POS) = sqrt(2)*[1-1i  -1+1i  1-1i  -1+1i  1-1i  1+1i  1+1i  -1-1i  -1-1i  -1+1i  1+1i  -1-1i];
STF_DAT_GRP(10,STF_POS) = sqrt(2)*[1-1i  1-1i  -1-1i  1+1i  1+1i  1+1i  1+1i  1-1i  -1-1i  -1+1i  -1+1i  -1-1i];  
%%% LTF
LTF_DAT_GRP = zeros(10, 52);
LTF_DAT_GRP(1,:) = [-1 1 -1 -1 -1 1 1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1 1 1 -1 -1 -1];
LTF_DAT_GRP(2,:) = [1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 1 -1];
LTF_DAT_GRP(3,:) = [1 -1 -1 -1 -1 1 -1 -1 -1 1 1 -1 1 1 1 1 1 -1 -1 1 -1 -1 -1 -1 -1 1 -1 1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 -1 1 1 -1 -1 1 1 1 -1 1 -1];
LTF_DAT_GRP(4,:) = [1 1 1 -1 -1 1 -1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 1 1 -1 -1 1 -1 -1 1 -1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 1 1 -1 -1 1 -1 -1 1 -1 1 -1 1 1];
LTF_DAT_GRP(5,:) = [1 -1 -1 1 1 1 1 1 -1 -1 -1 1 1 1 1 1 -1 -1 1 1 1 1 -1 -1 -1 1 1 1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 -1 1 -1 -1];
LTF_DAT_GRP(6,:) = [-1 1 -1 1 -1 1 1 1 1 1 -1 1 1 1 -1 -1 -1 1 -1 1 1 1 1 1 1 1 -1 1 1 1 1 1 -1 1 -1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 1 -1 -1 1 1 1];
LTF_DAT_GRP(7,:) = [-1 1 1 -1 1 -1 -1 1 1 -1 -1 -1 -1 -1 -1 1 -1 1 1 -1 1 1 -1 -1 1 -1 1 1 -1 1 1 -1 -1 1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 -1 -1 -1 -1 1 1 1];
LTF_DAT_GRP(8,:) = [1 -1 1 1 1 -1 -1 -1 1 -1 1 1 -1 -1 -1 1 1 -1 1 1 1 -1 1 -1 -1 -1 1 -1 1 1 1 1 -1 -1 -1 1 1 1 -1 -1 1 1 -1 1 -1 1 -1 1 -1 1 1 1];
LTF_DAT_GRP(9,:) = [-1 1 1 -1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 1 1 -1 -1 -1 -1 1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 1 1 1 -1 -1 -1 1 -1];
LTF_DAT_GRP(10,:) = [1 -1 -1 -1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 -1 1 -1 -1 -1 -1];
%%%% pilot
POS_VALID_SC = [39:64 2:27];
POS_PAYLOAD_SC = [1:5 7:19 21:32 34:46 48:52];
POS_PILOT_SC = [6 20 33 47];
PILOT_DAT_GRP = zeros(4,num_symbol_per_frame);
PILOT_DAT_GRP(1,:) = [-1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 -1 -1 1 -1 1];
PILOT_DAT_GRP(2,:) = [1 -1 -1 1 -1 1 1 1 1 1 1 -1 1 1 1 -1 -1 1 1 1];
PILOT_DAT_GRP(3,:) = [-1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 1 1 1 1 1];
PILOT_DAT_GRP(4,:) = [1 1 1 -1 1 1 -1 -1 -1 -1 1 -1 -1 -1 -1 -1 -1 1 1 -1];

   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_data_grp = zeros(num_data_stream, num_valid_sc, num_symbol_per_frame);
for ii = 1:num_data_stream
    stf_dat = STF_DAT_GRP(ii,:);
    ltf_dat = LTF_DAT_GRP(ii,:);
    data_freq = zeros(num_valid_sc,num_symbol_per_frame);
    data_freq(:,1) = stf_dat.';
    data_freq(:,2) = stf_dat.';
    data_freq(:,3) = ltf_dat.';
    data_freq(:,4) = ltf_dat.';
    for jj = 5:num_symbol_per_frame
        data_freq(:,jj) = (1/sqrt(2))*(sign(randn(num_valid_sc,1))+1i*sign(randn(num_valid_sc,1)));
        data_freq(POS_PILOT_SC,jj) = PILOT_DAT_GRP(:,jj);
    end
    freq_data_grp(ii,:,:) = data_freq;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_grp = zeros(num_ant_ap, 80*num_symbol_per_frame);
for ii = 1:num_ant_ap
    data_ant = squeeze(freq_data_grp(ii,:,:));
    data_f_frm = zeros(const_num_sc, num_symbol_per_frame);
    data_f_frm(POS_VALID_SC,:) = data_ant;
    data_t_frm = ifft(data_f_frm)*sqrt(const_num_sc);
    data_s_frm = [data_t_frm(end-15:end,:); data_t_frm];
    data_s_frm = data_s_frm(:).';
    time_signal_grp(ii,:) = data_s_frm;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how much data to write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_signal_grp = repmat(time_signal_grp, 1, const_num_frame);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID= fopen('signals/signal_tx.dat', 'wb');
    if (fileID < 0)
        disp('Error: fail to open files!');
        %pause;
    end
    data_tmp_c = time_signal_grp(ii,:);
    data_tmp_f = [real(data_tmp_c); imag(data_tmp_c)];
    data_tmp_f = data_tmp_f(:);
    fwrite(fileID,data_tmp_f,'float');
    fclose(fileID);




 


