clear;yalmip('clear');
address=[];
for x=0:200:1800
    for y=0:200:1800
         address=[address;[x y]];
    end
end
% wind=[12,0,1];
wind=[];
for n=0:35
    wind=[wind;[12,2*pi/36*n,1/36]];
end
DD=Wake_effect_v1(address,wind);
I=size(address,1); %IΪ��Ҫѡȡ�ĵ�ַ����ɢ����,��ʽΪaddress(x1 y1;x2 y2...��
D=size(wind,1);%DΪ�����ɢ����,��ʽΪwind(v1 ang1 p1;v2 ang2 p2;.....�����Ƕ�Ϊ����������˳ʱ�����ĽǶȴ�С
vv=1-DD;  %vvΪ��Ӧ���ٶȱ���
y=binvar(size(address,1),1,'full');
p=sdpvar(size(address,1),size(wind,1),'full');
Constraints=[sum(y)==30,p>=0];  %Լ�����ܻ���
for dd=1:D
    Constraints=[Constraints,p(:,dd)<=0.3*wind(dd,3)*(wind(dd,1)^3)*y];
end
loss_num=1-vv.^3;
for dd=1:D
    for jj= 1:I
    Constraints=[Constraints,p(jj,dd)<=0.3*wind(dd,3)*(wind(dd,1)^3)*(1-y'*loss_num(:,jj,dd))];
    end
end
Objective=sum(p,'all');
ops = sdpsettings('solver', 'gurobi','verbose', 2 ,'gurobi.Timelimit',600,'gurobi.NonConvex',2,'usex0',1);%���������
% read= readmatrix('Donovan_p.txt');  %����ֵ
% assign(y,read);
solution=optimize(Constraints,-Objective,ops)
aaa=value(y);
for n=1:length(aaa)
    if(0.999<=aaa(n)&&aaa(n)<=1.001)
        plot(address(n,1),address(n,2),'pentagram');
        hold on;
    end
end
title('���΢��ѡַ���ַ���');xlabel('����λ��/m');ylabel('����λ��/m');
writematrix(aaa,'Donovan_p.txt');
loss=0;
turbine_number=0;
power=0;
for dd=1:D
    for jj= 1:I
        for ii=1:I
            loss=loss+(DD(ii,jj,dd)^2)*aaa(ii);
        end
        turbine_number=turbine_number+aaa(jj)*(1-sqrt(loss))^3;  %V_3���㹦��
        loss=0;
    end
    power=power+0.3*wind(dd,3)*(wind(dd,1)^3)*turbine_number;
    turbine_number=0;
end
% loss=0;
% turbine_number=0;
% power=0;
% for dd=1:D
%     for jj= 1:I
%         for ii=1:I
%             loss=loss+(1-vv(ii,jj,dd)^3)*aaa(ii)*aaa(jj);
%         end
%         turbine_number=turbine_number+aaa(jj)-loss;
%         loss=0;
%     end
%     power=power+0.3*wind(dd,3)*(wind(dd,1)^3)*turbine_number;
%     turbine_number=0;
% end
power
power_obj=value(Objective)