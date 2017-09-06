function [pos_j,C_qos] = srt_pos_and_qos( J_range,T_range,delta_t,interval,ratio )
%SRT_POS_J Summary of this function goes here
%   Detailed explanation goes here

%认为有两条航线，近似于彼此平行；第三条航线穿过前两条
%其余船只为渔船，随机分布

%航线1，航线2，航线3，渔船 速度
dir=[1,-1,-0.6+0.8i];
v=[36,36,36,0]/3.6;%单位：m/s


%航线1，航线2，航线3 发船间隔
interval=interval*ones(1,3);%单位：s

%航线1，航线2，航线3，渔船对应j
j_index_1=(1 : floor(J_range*ratio(1)));
j_index_2=(floor(J_range*ratio(1))+1 : floor(J_range*(ratio(1)+ratio(2))));
j_index_3=(floor(J_range*(ratio(1)+ratio(2)))+1 : floor(J_range*(ratio(1)+ratio(2)+ratio(3))));
j_index_rand=(floor(J_range*sum(ratio))+1 : J_range);
%航线1，航线2，航线3 起点 & 距离基站最近距离
init_dist=[-10+9i,15+12i,12]*10^3;


pos_j=zeros(J_range,T_range);

pos_j(j_index_1,:)=(1:length(j_index_1))'*ones(1,T_range)*(-interval(1))*v(1)*dir(1) + init_dist(1) + v(1)*dir(1)*delta_t*ones(length(j_index_1),1)*(1:T_range);
pos_j(j_index_2,:)=(1:length(j_index_2))'*ones(1,T_range)*(-interval(2))*v(2)*dir(2) + init_dist(2) + v(2)*dir(2)*delta_t*ones(length(j_index_2),1)*(1:T_range);
pos_j(j_index_3,:)=(1:length(j_index_3))'*ones(1,T_range)*(-interval(3))*v(3)*dir(3) + init_dist(3) + v(3)*dir(3)*delta_t*ones(length(j_index_3),1)*(1:T_range);
direction_rand=rand(length(j_index_rand),1);
direction_rand=direction_rand+(1-direction_rand.^2).^0.5*1i;
pos_j(j_index_rand,:)=(random('Uniform',-10*10^3,10*10^3,length(j_index_rand),1)+random('Uniform',2*10^3,10*10^3,length(j_index_rand),1)*1i)*ones(1,T_range) + v(3)*delta_t*direction_rand*(1:T_range);

%qos
C_qos=1*8*10^6*ones(1,J_range);%MB
C_qos(j_index_rand)=C_qos(j_index_rand)/10;%MB
end

