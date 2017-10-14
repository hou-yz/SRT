function [E_stage_2,E_stage_3,C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0)
P=reshape([P_i_max;ones(J_range,1)*P_j_max]*ones(1,J_range*T_range),J_range+1,J_range,T_range);
P1_tmp=linspace(0,P_i_max);P2_tmp=linspace(0,P_j_max);
%fast fading, size(h0,3)>T_range
T_range_h0=size(h0,3);
t_step=T_range_h0/T_range;
r=zeros(J_range+1,J_range,T_range);
for t=(1:T_range)
    t_h0=(t-1)*t_step;
    r(:,:,t)=sum(BW*log(1+P(:,:,t).*ones(1,1,t_step).*(beta(:,:,t).*ones(1,1,t_step)).*abs(h0(:,:,t_h0+1:t_h0+t_step)).^2/sigma2)/log(2),3);
end
r=r/t_step;
    
    
% r=BW*log(1+P.*beta.*h0.^2/sigma2)/log(2);

S=[];%[i,j,t]->[0,j,t,0,0,0]
%����Ϊ�㷨��һ����
% eta = ����ʱ�� / ʱ϶����
eta=zeros(J_range+1,J_range,T_range);
C_tmp=zeros(J_range,T_range+1);
C_max=C_qos;
%C_tmp(j,t) means have that much data at the time slot before t

for j=(1:J_range) %�����û�j
    [~,beta_j_sorted_t_index]=sort(beta(1,j,:),'descend');
    ijt_t_pointer=1;
    while C_tmp(j,T_range+1)<C_max(j) %��������qos����������һ����
        t0=beta_j_sorted_t_index(ijt_t_pointer);
        if r(0+1,j,t0)==0
            %fprintf('srt:�û�%dδ������QoS����',j)
            C_max(j)=C_tmp(j,T_range+1);
            break
        end
        if C_tmp(j,T_range+1)+r(0+1,j,t0)*delta_t<C_max(j)%δ���
            eta(0+1,j,t0)=1;%index 0 for BS-i
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0)*delta_t*eta(0+1,j,t0);
        else%���
            eta(0+1,j,t0)=(C_max(j)-C_tmp(j,T_range+1))/(r(0+1,j,t0)*delta_t);
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0)*delta_t*eta(0+1,j,t0);
        end
        
        S=[S;[0,j,t0,0,0,0]];
        ijt_t_pointer=ijt_t_pointer+1;
        if ijt_t_pointer>T_range
            %fprintf('srt:�û�%dδ������QoS����',j)
            C_max(j)=C_tmp(j,T_range+1);
        end
    end
    
end
%��һ���ֽ���
E_stage_1=sum(sum(eta(1,:,:)))*P_i_max*delta_t;
%fprintf('��һ���ֽ�����ϵͳ�ܹ��ģ�%d\n',E_stage_1)
delta_1=reshape(eta(0+1,:,:),J_range,T_range);

%����Ϊ�㷨�ڶ�����
t_free=ones(1,T_range);%time slots that are availible ��ʼȫ��ʱ϶Ϊ����
t_OOC=zeros(1,T_range);%out of capacity ��ʼ��ʱ϶������������
while sum(t_OOC)~=0 %������վ������Χ
    %����
    for t=(1:T_range)
        t_free(t)=(sum(eta(0+1,:,t)>0)<=N);
        t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
    end
    r_i_free=reshape(r(1,:,:).*(1-eta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free);% ����ʱ϶�� ���� ����
    r_i_OOC=reshape(r(1,:,:),J_range,T_range).*(ones(J_range,1)*t_OOC);r_i_OOC(eta(0+1,:,:)==0)=0;r_i_OOC(r_i_OOC==0)=NaN;% ��������ʱ϶��Ӧ����
    delta_C=zeros(1,J_range);% ��¼���û�j�������仯����ϵͳ����Ӱ�죩��С��
    t0_index=zeros(1,J_range);% ÿ���û�j�������仯����ϵͳ����Ӱ�죩��С���Ӧ ԭ ʱ϶t0
    for j=(1:J_range)
        [delta_C(j),t0_index(j)]=min(r_i_OOC(j,:)*delta_t);
        delta_C(j)=delta_C(j)-max(r_i_free(j,:))*delta_t;
    end
    [~,j]=min(delta_C); %ȥ��Ч�ʣ�beta���ϲ�ı�
    t0=t0_index(j);
    
    [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
    S(S_index,:)=[];
    C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r(0+1,j,t0)*eta(0+1,j,t0)*delta_t;
    eta(0+1,j,t0)=0;
    %fprintf('ȥ��һ���ߣ�[0 %d %d]\n',j,t0)
    %����
    for t=(1:T_range)
        t_free(t)=(sum(eta(0+1,:,t)>0)<=N);
        t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
    end
    r_i_free=reshape(r(1,:,:).*(1-eta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free);
    
    while C_tmp(j,T_range+1)<C_max(j)%ȥ���ñߺ�Ѱ��һ����Ž�����QoS����
        if sum(r_i_free(j,:))<=0
            %fprintf('srt:�û�%dδ������QoS����',j)
            C_max(j)=C_tmp(j,T_range+1);
            break;
        end
        [r_ijt0,t0]=max(r_i_free(j,:));
        if C_tmp(j,T_range+1)+r_ijt0*delta_t<C_max(j)%δ���
            eta(0+1,j,t0)=1;% ��֮ǰdeltaȡֵ�޹أ�����ȫ������Ϊ1
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r_ijt0*delta_t*eta(0+1,j,t0);
        else%���
            eta(0+1,j,t0)=eta(0+1,j,t0)+(C_max(j)-C_tmp(j,T_range+1))/(r_ijt0*delta_t);
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r_ijt0*delta_t*eta(0+1,j,t0);
        end
        
        S=[S;[0,j,t0,0,0,0]];
        %fprintf('����һ���ߣ�[0 %d %d]\n',j,t0)
        %����
        for t=(1:T_range)
            t_free(t)=(sum(eta(0+1,:,t)>0)<=N);
            t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
        end
        r_i_free=reshape(r(1,:,:).*(1-eta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free);
    end    
end 
%�ڶ����ֽ���
E_stage_2=sum(sum(eta(1,:,:)))*P_i_max*delta_t;
%fprintf('�ڶ����ֽ�����ϵͳ���ܺģ�%d\n',E_stage_2)
%reshape(delta(0+1,:,:),J_range,T_range)
delta_2=reshape(eta(0+1,:,:),J_range,T_range);
C_unsatisfied=C_qos-C_max;

%��������
r0_index=find(eta(0+1,:,:)>0);
[~,jj,tt]=ind2sub(size(eta(0+1,:,:)),r0_index);
for j=(1:J_range)%for all user j
    R=[];
    r0_j=r(0+1,j,tt(jj==j));
    if isempty(r0_j)
        continue
    end
    [r0_max,~]=max(r0_j);
    for j_prime_R=(1:J_range)
        
        t2s=(1:T_range);
        t2s=t2s(r(j_prime_R+1,j,:)>r0_max);
        for t2=t2s
            %disp([j_prime_R j sum(delta(:,j_prime_R,t2)>0)+sum(delta(j_prime_R+1,:,t2)>0)-(delta(j_prime_R+1,j,t2)>0) sum(delta(:,j,t2)>0)-(delta(j_prime_R+1,j,t2)>0)+sum(delta(j+1,:,t2)>0) t2])
            if j_prime_R~=j && sum(eta(:,j_prime_R,t2)>0)+sum(eta(j_prime_R+1,:,t2)>0)-2*(eta(j_prime_R+1,j,t2)>0)==0 && sum(eta(:,j,t2)>0)-2*(eta(j_prime_R+1,j,t2)>0)+sum(eta(j+1,:,t2)>0)==0 && sum(sum(eta(:,:,t2)>0))-(eta(j_prime_R+1,j,t2)>0)<N && eta(j_prime_R+1,j,t2)<1% j' and j and sys free
                
                t1s=(1:T_range);
                t1s=t1s(r(0+1,j_prime_R,:)>r0_max);
                for t1=t1s
                    if t1~=t2 && sum(eta(:,j_prime_R,t1)>0)-2*(eta(0+1,j_prime_R,t1)>0)+sum(eta(j_prime_R+1,:,t1)>0)==0 && sum(sum(eta(:,:,t1)>0))-(eta(0+1,j_prime_R,t1)>0)<N && eta(0+1,j_prime_R,t1)<1 % j' free and sys availible
                        r1=r(1,j_prime_R,t1);%r1
                        if ((t1>t2 && C_tmp(j_prime_R,t2)>=r0_max*delta_t) || (t1<=t2 && C_tmp(j_prime_R,t2)+r1*delta_t>=r0_max*delta_t))
                              R=[R;[0,j_prime_R,t1,j_prime_R,j,t2]];
                        end
                    end
                end
            end
        end
        
    end
    while ~isempty(R)
        r0_index=find(eta(0+1,:,:));% ���¿�ѡr0��Ӧ����
        [~,jj0,tt0]=ind2sub(size(eta(0+1,:,:)),r0_index);%��Ӧr0:(0+1,j,t0)
%         r0_j=r(0+1,j,tt0(jj0==j));
        delta_P=zeros(length(tt0(jj0==j))*size(R,1),3);edges=zeros(length(tt0(jj0==j))*size(R,1),9);
        for i_0=(1:length(tt0(jj0==j)))
            ttt=tt0(jj0==j);t0=ttt(i_0);
            edge0=[0,j,t0];
            
            r0=sum(BW*log(1+P(0+1,j,t0).*ones(1,1,t_step).*(eta(0+1,j,t0).*ones(1,1,t_step)).*(beta(0+1,j,t0).*ones(1,1,t_step)).*abs(h0(0+1,j,(t0-1)*t_step+1:t0*t_step)).^2/sigma2)/log(2),3)/t_step;
%             r0_prime=r(0+1,j,t0)*eta(0+1,j,t0);
            for i_12=(1:size(R,1))
                edge1=R(i_12,1:3);edge2=R(i_12,4:6);
                j_prime=edge1(2);
                t1=edge1(3);t2=edge2(3);
                beta1=beta(0+1,j_prime,t1);beta2=beta(j_prime+1,j,t2);
                if t_step==1
                    h0_1=1;h0_2=1;
                    P1=(2^(r0/BW)-1)*sigma2/(beta1*abs(h0_1)^2);%*eta(0+1,j,t0)
                    P2=(2^(r0/BW)-1)*sigma2/(beta2*abs(h0_2)^2);%*eta(0+1,j,t0)
                else
                    r1_tmp=sum(BW*log(1+P1_tmp.*ones(1,1,t_step).*(beta1.*ones(1,length(P1_tmp),t_step)).*(abs(h0(0+1,j,(t0-1)*t_step+1:t0*t_step).*ones(1,length(P1_tmp),t_step))).^2/sigma2)/log(2),3)/t_step;
                    [~,P1_index]=min(abs(r1_tmp-r0));
                    P1=P1_tmp(P1_index);
                    r2_tmp=sum(BW*log(1+P2_tmp.*ones(1,1,t_step).*(beta2.*ones(1,length(P2_tmp),t_step)).*(abs(h0(0+1,j,(t0-1)*t_step+1:t0*t_step).*ones(1,length(P2_tmp),t_step))).^2/sigma2)/log(2),3)/t_step;
                    [~,P2_index]=min(abs(r2_tmp-r0));
                    P2=P2_tmp(P2_index);
                end
                ii=(i_0-1)*size(R,1)+i_12;
                if eta(0+1,j_prime,t1)+P1/P_i_max<1 && eta(j_prime+1,j,t2)+P2/P_j_max<1 && sum(eta(:,j_prime,t1)>0)-2*(eta(0+1,j_prime,t1)>0)+sum(eta(j_prime+1,:,t1)>0)==0 && sum(eta(:,j_prime,t2)>0)+sum(eta(j_prime+1,:,t2)>0)-2*(eta(j_prime+1,j,t2)>0)==0 && sum(eta(j+1,:,t2)>0)+sum(eta(:,j,t2)>0)-2*(eta(j_prime+1,j,t2)>0)==0 && sum(sum(eta(:,:,t2)>0))-(eta(j_prime+1,j,t2)>0)<N && sum(sum(eta(:,:,t1)>0))-(eta(0+1,j_prime,t1)>0)<N
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
            r0=r(0+1,j,t0);
            P1=delta_P(ii,1);P2=delta_P(ii,2);
            %ȥ��edge0��.
            eta(0+1,j,t0)=0;
            [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
            if S_index~=0
                S(S_index,:)=[];
            else
                break
            end
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r0*delta_t;
            %fprintf('ȥ��һ���ߣ�[0 %d %d]\n',j,t0)
            %����edge2��
            eta(j_prime+1,j,t2)=eta(j_prime+1,j,t2)+P2/P_j_max;
            C_tmp(j,t2:T_range+1)=C_tmp(j,t2:T_range+1)+r0*delta_t;
            C_tmp(j_prime,t2:T_range+1)=C_tmp(j_prime,t2:T_range+1)-r0*delta_t;%*(1-gamma(j_prime,j))
            %����edge1��
            eta(0+1,j_prime,t1)=eta(0+1,j_prime,t1)+P1/P_i_max;
            S=[S;[0,j_prime,t1,j_prime,j,t2]];
            C_tmp(j_prime,t1:T_range+1)=C_tmp(j_prime,t1:T_range+1)+r0*delta_t;%*(1-gamma(j_prime,j))
            %fprintf('���������ߣ�[0 %d %d] [%d %d %d]\n',j_prime,t1,j_prime,j,t2)
            
            %��R��ȥ��
            %[~,R_index]=ismember([0,j_prime,t1,j_prime,j,t2],R,'rows');
            tmp_max=size(R,1);tmp=1;
            while tmp<tmp_max
                if R(tmp,6)==t2
                    R(tmp,:)=[];
                    tmp=tmp-1;
                    tmp_max=tmp_max-1;
                end
                tmp=tmp+1;
            end
        else
            break
        end
    end
end
E_stage_3=sum(sum(sum(eta,3),2).*[P_i_max;ones(J_range,1)*P_j_max]*delta_t);
%fprintf('�������ֽ�����ϵͳ���ܺģ�%d\n',E_stage_3)

end    