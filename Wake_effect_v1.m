   function DD = Wake_effect_v1(address,wind)
%   尾流损失生成函数 输入格式为(address,wind)
%   其中具体格式address(x1 y1;x2 y2...) wind(v1 ang1 p1;v2 ang2 p2;.....) ang为从正北方向顺时针旋转所需的角度 
%   在此模型中考虑重叠面积的影响，使用非高斯分布的jessen模型
I=size(address,1); %I为将要选取的地址的离散集合,格式为(x1 y1;x2 y2...）
D=size(wind,1);%D为风的离散集合,格式为(v1 ang1 p1 ;v2 ang2 p2;.....），角度为与正北方向按顺时针计算的角度大小
alpha=0.1;
rd=20;
DD=zeros(I,I,D);
Ct=0.88;
    for k=1:D  %每个风向、风速
        for j=1:I    %每个受影响的j
            for i=1:I   %对每个j风机来说不同的i风机对其的影响
                d=(address(i,2)-address(j,2))*cos(wind(k,2))+(address(i,1)-address(j,1))*sin(wind(k,2));%d为两个风机在风向上的距离
                d_v=abs((address(i,2)-address(j,2))*sin(wind(k,2))-(address(i,1)-address(j,1))*cos(wind(k,2)));%l为两个风机在与风向垂直方向上的距离
                if(d<=0)  %d小于等于0，尾流对其无影响
                    DD(i,j,k)=0;
                else     %d>0
                            if(d_v<=d*alpha)
                            DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2);
                            elseif(d*alpha<d_v&&d_v<=rd+d*alpha)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+(pi-acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v)))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v;
                               ita=A_shadow/(pi*rd^2);
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita; 
                            elseif(d_v>rd+d*alpha&&d_v<=rd+rd+d*alpha)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v;
                               ita=A_shadow/(pi*rd^2);
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita;   
                            elseif(d_v>rd+rd+d*alpha)  
                            DD(i,j,k)=0;
                            end
                end
            end
        end
    end
end

