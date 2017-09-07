function [E_stage_1,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2)
%REFERENCE_ALGORITHM Summary of this function goes here
%   Detailed explanation goes here

r=srt_channel([P_i_max;ones(J_range,1)*P_j_max],BW,beta,sigma2);
S=[];%[i,j,t]->[0,j,t,0,0,0]
C_max=C_qos;

%����Ϊ�㷨��һ����
delta=zeros(J_range+1,J_range,T_range);
C_tmp=zeros(J_range,T_range+1);
%C_tmp(j,t) means have that much data at the time slot before t
unsatisfied_j=(1:J_range);
for t0=(1:T_range) %����ʱ϶
    if isempty(unsatisfied_j)
        break
    end
    [~,beta_sorted_j_pointer]=sort(beta(1,unsatisfied_j,t0),'descend');
    ijt_j_pointer=1;
    j=unsatisfied_j(beta_sorted_j_pointer(ijt_j_pointer));
    satisfied_cnt=0;
    
    while C_tmp(j,T_range+1)<C_qos(j) && ijt_j_pointer<=N %��������qos����������һ����  
        [~,beta_sorted_j_pointer]=sort(beta(1,unsatisfied_j,t0),'descend');
        if r(0+1,j,t0)==0
            %fprintf('srt:�û�%dδ������QoS����',j)
            break
        end
        if C_tmp(j,T_range+1)+r(0+1,j,t0)*delta_t<C_qos(j)%δ���
            delta(0+1,j,t0)=1;%index 0 for BS-i
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0)*delta_t*delta(0+1,j,t0);
            if C_tmp(j,T_range+1)==C_qos(j)
                [~,j_index]=ismember(j,unsatisfied_j);
                unsatisfied_j(j_index)=[];
                [~,beta_sorted_j_pointer]=sort(beta(1,unsatisfied_j,t0),'descend');
                satisfied_cnt=satisfied_cnt+1;
            end
        else%���
            delta(0+1,j,t0)=(C_qos(j)-C_tmp(j,T_range+1))/(r(0+1,j,t0)*delta_t);
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0)*delta_t*delta(0+1,j,t0);
             %��unsatisfied_j��ȥ��
            [~,j_index]=ismember(j,unsatisfied_j);
            unsatisfied_j(j_index)=[];
            [~,beta_sorted_j_pointer]=sort(beta(1,unsatisfied_j,t0),'descend');
            satisfied_cnt=satisfied_cnt+1;
        end
        if isempty(unsatisfied_j)
            break
        end
        S=[S;[0,j,t0,0,0,0]];
        ijt_j_pointer=ijt_j_pointer+1;
        if ijt_j_pointer-satisfied_cnt>length(unsatisfied_j)
            break
        else
            j=unsatisfied_j(beta_sorted_j_pointer(ijt_j_pointer-satisfied_cnt));
        end
    end
    
end
if ~isempty(unsatisfied_j)
    for j=unsatisfied_j
        %fprintf('ref:�û�%dδ������QoS����',j)
        C_max(j)=C_tmp(j,T_range+1);
    end
end
%��һ���ֽ���
E_stage_1=sum(sum(delta(1,:,:)))*P_i_max*delta_t;
%fprintf('��һ���ֽ�����ϵͳ�ܹ��ģ�%d\n',E_stage_1)
delta_1=reshape(delta(0+1,:,:),J_range,T_range);
C_unsatisfied=C_qos-C_max;

end

