clear;yalmip('clear');
address=[];
for x=0:200:1800
    for y=0:200:1800
         address=[address;[x y]];
    end
end
% wind=[12,0,1];
wind=[];
for n=0:11
    wind=[wind;[12,2*pi/12*n,1/12]];
end
DD=Wake_effect_v1(address,wind);
vv=1-DD; 
y=binvar(size(address,1),1,'full');
% ww=binvar(size(address,1),size(address,1),'full');
Constraints=[sum(y)==19];  %约束，总机数
% tic
% for i=1:size(address,1)
%     for j=1:size(address,1)
%         if (i==j)
%         Constraints=[Constraints,ww(i,j)==y(i)];
%         else
%              Constraints=[Constraints,ww(i,j)<=y(i),ww(i,j)<=y(j),ww(i,j)>=y(i)+y(j)-1];
%         end
%     end
% end
% toc
Objective = 0;
I=size(address,1); %I为将要选取的地址的离散集合,格式为address(x1 y1;x2 y2...）
D=size(wind,1);%D为风的离散集合,格式为wind(v1 ang1 p1;v2 ang2 p2;.....），角度为与正北方向按顺时针计算的角度大小
KK=DD.^2;
tic
for dd=1:D
            Objective=Objective+wind(dd,3)*(wind(dd,1)^3)*y'*KK(:,:,dd)*y;
end
toc
ops = sdpsettings('solver', 'gurobi','verbose', 2,'gurobi.NonConvex',2,'gurobi.Timelimit',1200,'gurobi.Mipgap',0.05);%设置求解器
% read= readmatrix('A1.txt');  %赋初值
% assign(y,read);
solution=optimize(Constraints,Objective,ops)
aaa=value(y);
for n=1:length(aaa)
    if(aaa(n)==1)
        plot(address(n,1),address(n,2),'+');
        hold on;
    end
end
writematrix(aaa,'A1.txt')
read= readmatrix('A1.txt');
loss=0;
turbine_number=0;
power=0;
% for dd=1:D
%     for jj= 1:I
%         for ii=1:I
%             loss=loss+(DD(ii,jj,dd)^2)*aaa(ii);
%         end
%         turbine_number=turbine_number+aaa(jj)*(1-sqrt(loss))^3;
%         loss=0;
%     end
%     power=power+0.3*wind(dd,3)*(wind(dd,1)^3)*turbine_number;
%     turbine_number=0;
% end
for dd=1:D
    for jj= 1:I
        for ii=1:I
            loss=loss+(1-vv(ii,jj,dd)^3)*aaa(ii)*aaa(jj);
        end
        turbine_number=turbine_number+aaa(jj)-loss;
        loss=0;
    end
    power=power+0.3*wind(dd,3)*(wind(dd,1)^3)*turbine_number;
    turbine_number=0;
end
power