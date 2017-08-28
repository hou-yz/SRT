clear
clc
T_range=600;%an hour
delta_t=6;
J_range=3;
N=1;
P_i_max=10;P_j_max=1;
BW=20*10^6;%bandwidth=20MHz
%pos_i=zeros(2,1);pos_j=ones(2,t_range,J_range);
pos_j=srt_pos_j(J_range,T_range);%rand(J_range,T_range)+1i*rand(J_range,T_range);
C_qos=10^7*ones(1,J_range);
pos_i=[zeros(1,T_range);pos_j];
r=zeros(J_range+1,J_range,T_range);
beta=srt_beta(pos_i,pos_j);
gamma=0*ones(J_range,J_range);
[srt_E1,srt_E2]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
[ref_E1,ref_E2]=reference_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);