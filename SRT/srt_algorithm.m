function [E_stage_2,E_stage_3]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW)
r=srt_channel([P_i_max;ones(J_range,1)*P_j_max],BW,beta);
S=[];%[i,j,t]->[0,j,t,0,0,0]
%����Ϊ�㷨��һ����
delta=zeros(J_range+1,J_range,T_range);
C_tmp=zeros(J_range,T_range+1);
%C_tmp(j,t) means have that much data at the time slot before t

for j=(1:J_range) %�����û�j
    [~,beta_j_sorted_t_index]=sort(beta(1,j,:),'descend');
    ijt_t_pointer=1;
    while C_tmp(j,T_range+1)<C_qos(j) %��������qos����������һ����
        t0=beta_j_sorted_t_index(ijt_t_pointer);
        if C_tmp(j,T_range+1)+r(0+1,j,t0)*delta_t<C_qos(j)%δ���
            delta(0+1,j,t0)=1;%index 0 for BS-i
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0)*delta_t*delta(0+1,j,t0);
        else%���
            delta(0+1,j,t0)=(C_qos(j)-C_tmp(j,T_range+1))/(r(0+1,j,t0)*delta_t);
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0)*delta_t*delta(0+1,j,t0);
        end
        
        S=[S;[0,j,t0,0,0,0]];
        ijt_t_pointer=ijt_t_pointer+1;
        if ijt_t_pointer>=T_range
            fprintf('�û�%dδ������QoS����',j)
            return
        end
    end
    
end
%��һ���ֽ���
E_stage_1=sum(sum(delta(1,:,:)))*P_i_max*delta_t;
%fprintf('��һ���ֽ�����ϵͳ�ܹ��ģ�%d\n',E_stage_1)
%reshape(delta(0+1,:,:),J_range,T_range)

%����Ϊ�㷨�ڶ�����
t_free=ones(1,T_range);%time slots that are availible ��ʼȫ��ʱ϶Ϊ����
t_OOC=zeros(1,T_range);%out of capacity ��ʼ��ʱ϶������������
for t=(1:T_range)
    t_free(t)=(sum(delta(0+1,:,t)>0)<N);
    t_OOC(t)=(sum(delta(0+1,:,t)>0)>N);
end
while sum(t_OOC)~=0 %������վ������Χ
    r_i_free=reshape(r(1,:,:).*(1-delta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free);% ����ʱ϶�� ���� ����
    r_i_OOC=reshape(r(1,:,:).*delta(0+1,:,:),J_range,T_range).*(ones(J_range,1)*t_OOC);r_i_OOC(r_i_OOC==0)=inf;% ��������ʱ϶��Ӧ����
    delta_C=zeros(1,J_range);% ��¼���û�j�������仯����ϵͳ����Ӱ�죩��С��
    t0_index=zeros(1,J_range);% ��¼���û�j�������仯����ϵͳ����Ӱ�죩��С���Ӧ ԭ ʱ϶
    for j=(1:J_range)
        [delta_C(j),t0_index(j)]=min(r_i_OOC(j,:)*delta_t);
        delta_C(j)=delta_C(j)-max(r_i_free(j,:))*delta_t;
    end
    [~,j]=min(delta_C); %�õ�ȥ���ߺ�Ч����õ��û�j
    t0=t0_index(j);
    %ȥ��ǰϵͳӰ����С�ıߣ������仯��С��.
    [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
    S(S_index,:)=[];
    C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r(0+1,j,t0)*delta(0+1,j,t0)*delta_t;
    delta(0+1,j,t0)=0;
    %fprintf('ȥ��һ���ߣ�[0 %d %d]\n',j,t0)
    
    while C_tmp(j,T_range+1)<C_qos(j)%ȥ���ñߺ�Ѱ��һ����Ž�����QoS����
        if sum(sum(r_i_free))==0
            fprintf('�û�%dδ������QoS����',j)
            return
        end
        [r_ijt0,t0]=max(r_i_free(j,:));
        if C_tmp(j,T_range+1)+r_ijt0*delta_t<C_qos(j)%δ���
            delta(0+1,j,t0)=1;% ��֮ǰdeltaȡֵ�޹أ�����ȫ������Ϊ1
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r_ijt0*delta_t*delta(0+1,j,t0);
        else%���
            delta(0+1,j,t0)=delta(0+1,j,t0)+(C_qos(j)-C_tmp(j,T_range+1))/(r_ijt0*delta_t);
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r_ijt0*delta_t*delta(0+1,j,t0);
        end
        
        S=[S;[0,j,t0,0,0,0]];
        %fprintf('����һ���ߣ�[0 %d %d]\n',j,t0)
        %����
        for t=(1:T_range)
            t_free(t)=(sum(delta(0+1,:,t)>0)<N);
            t_OOC(t)=(sum(delta(0+1,:,t)>0)>N);
        end
        r_i_free=reshape(r(1,:,:).*(1-delta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free);
    end    
end 
%�ڶ����ֽ���
E_stage_2=sum(sum(delta(1,:,:)))*P_i_max*delta_t;
fprintf('�ڶ����ֽ�����ϵͳ���ܺģ�%d\n',E_stage_2)
%reshape(delta(0+1,:,:),J_range,T_range)
delta2=delta;

%��������
r0_index=find(delta(0+1,:,:)>0);
[~,jj,tt]=ind2sub(size(delta(0+1,:,:)),r0_index);
for j_loop=(1:J_range)%for all user j_loop
    R=[];
    r0_j=r(0+1,j_loop,tt(jj==j_loop));
    [r0_min,~]=min(r0_j);
    for t2=(1:T_range)
        for j_prime=(1:J_range)
            r2=r(j_prime+1,j_loop,t2);%r2
            if j_prime~=j_loop && sum(delta(:,j_prime,t2)>0)+sum(delta(j_prime+1,:,t2)>0)-(delta(j_prime+1,j_loop,t2)>0)==0 && sum(delta(:,j_loop,t2)>0)-(delta(j_prime+1,j_loop,t2)>0)+sum(delta(j_loop+1,:,t2)>0)==0 && r2>r0_min % j' and j free
                for t1=(1:T_range)
                    if t1~=t2 && sum(delta(:,j_prime,t1)>0)-(delta(0+1,j_prime,t1)>0)+sum(delta(j_prime+1,:,t1)>0)==0 && sum(delta(0+1,:,t1)>0)-(delta(0+1,j_prime,t1)>0)<N % j' free and i availible
                        r1=r(1,j_prime,t1);%r1
                        if ((t1>t2 && C_tmp(j_prime,t2)>=r0_min*delta_t) || (t1<=t2 && C_tmp(j_prime,t2)+r1*delta_t>=r0_min*delta_t)) && r1>r0_min
                            R=[R;[0,j_prime,t1,j_prime,j_loop,t2]];
                        end
                    end
                end
            end
        end
    end
    while ~isempty(R)
        r0_index=find(delta(0+1,:,:));% ���¿�ѡr0��Ӧ����
        [~,jj,tt]=ind2sub(size(delta(0+1,:,:)),r0_index);
        r0_j=r(0+1,j_loop,tt(jj==j_loop));
        delta_P=zeros(length(tt(jj==j_loop))*size(R,1),3);edges=zeros(length(tt(jj==j_loop))*size(R,1),9);
        for i_0=(1:length(tt(jj==j_loop)))
            j=j_loop;ttt=tt(jj==j_loop);t0=ttt(i_0);
            edge0=[0,j,t0];
            r0_this=r0_j(i_0);
            for i_12=(1:size(R,1))
                edge1=R(i_12,1:3);edge2=R(i_12,4:6);
                j_prime=edge1(2);
                t1=edge1(3);t2=edge2(3);
                beta1=beta(edge1(1)+1,edge1(2),edge1(3));beta2=beta(edge2(1)+1,edge2(2),edge2(3));
                P1=srt_get_p(r0_this*(1-gamma(j_prime,j))*delta(0+1,j,t0),BW,beta1);
                P2=srt_get_p(r0_this*delta(0+1,j,t0),BW,beta2);
                ii=(i_0-1)*size(R,1)+i_12;
                if delta(0+1,j_prime,t1)+P1/P_i_max<1 && delta(j_prime+1,j,t2)+P2/P_j_max<1 && sum(delta(0+1,:,t1)>0)-(delta(0+1,j_prime,t1)>0)==0 && sum(delta(:,j_prime,t1)>0)-(delta(0+1,j_prime,t1)>0)+sum(delta(j_prime+1,:,t1)>0)==0 && sum(delta(:,j_prime,t2)>0)+sum(delta(j_prime+1,:,t2)>0)-(delta(j_prime+1,j,t2)>0)==0 && sum(delta(j+1,:,t2)>0)+sum(delta(:,j,t2)>0)-(delta(j_prime+1,j,t2)>0)==0
                    delta_P(ii,:)=[P1,P2,P_i_max-(P1+P2)];
                else
                    delta_P(ii,:)=[P1,P2,-inf];
                end
                edges(ii,:)=[edge0,edge1,edge2];
            end
        end
        [m,ii]=max(delta_P(:,3));%�õ����ŵ��滻��ʽ delta_P(ii,:)=[P1,P2,P_i_max-(P1+P2)];
        if m>0
            %�滻
            edge0=edges(ii,1:3);edge1=edges(ii,4:6);edge2=edges(ii,7:9);
            j_prime=edge1(2);
            t0=edge0(3);t1=edge1(3);t2=edge2(3);
            r0_this=r(0+1,j,t0);
            P1=delta_P(ii,1);P2=delta_P(ii,2);
            %ȥ��edge0��.
            delta(0+1,j,t0)=0;
            [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
            if S_index~=0
                S(S_index,:)=[];
            else
                break
            end
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r0_this*delta_t;
            %fprintf('ȥ��һ���ߣ�[0 %d %d]\n',j,t0)
            %����edge2��
            delta(j_prime+1,j,t2)=delta(j_prime+1,j,t2)+P2/P_j_max;
            C_tmp(j,t2:T_range+1)=C_tmp(j,t2:T_range+1)+r0_this*delta_t;
            C_tmp(j_prime,t2:T_range+1)=C_tmp(j_prime,t2:T_range+1)-r0_this*delta_t;%*(1-gamma(j_prime,j))
            %����edge1��
            delta(0+1,j_prime,t1)=delta(0+1,j_prime,t1)+P1/P_i_max;
            S=[S;[0,j_prime,t1,j_prime,j,t2]];
            C_tmp(j_prime,t1:T_range+1)=C_tmp(j_prime,t1:T_range+1)+r0_this*delta_t;%*(1-gamma(j_prime,j))
            %fprintf('���������ߣ�[0 %d %d] [%d %d %d]\n',j_prime,t1,j_prime,j,t2)
            
            %��R��ȥ��
            [~,R_index]=ismember([0,j_prime,t1,j_prime,j,t2],R,'rows');
            R(R_index,:)=[];
        else
            break
        end
    end
end
E_stage_3=sum(sum(sum(delta,3),2).*[P_i_max;ones(J_range,1)*P_j_max]*delta_t);
fprintf('�������ֽ�����ϵͳ���ܺģ�%d\n',E_stage_3)

if E_stage_3>E_stage_2
   delta2-delta 
end
end    