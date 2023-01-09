# WFLO
==============================
### Wake_effect_v1.m
尾流效应因子矩阵计算模块，输入格式为(address,wind)，address为各个微观选址点的坐标，具体格式address=(x1 y1;x2 y2...)  
wind=(v1 ang1 p1;v2 ang2 p2;.....)  
其中v为风速，ang为从正北方向顺时针旋转所需的角度（例如北风的ang=0，东风的ang=π/2)
### Microscopic_siting.m 
复现*Turner, S. D. O., D. A. Romero, P. Y. Zhang, C. H. Amon, and T. C. Y. Chan. "A new mathematical programming approach to optimize wind farm layouts." Renewable Energy 63 (2014): 674-680.*
目标函数为 
$$∑_{𝑊𝑖𝑛𝑑}∑_{𝑗=1}^𝑁∑_{𝑖=1}^𝑁 p_{wind} v_{\infty}^3 d_{ij}^2 y_jy_i$$
### Microscopic_siting_Donovan_p.m
Donovan提出的优化模型，将目标函数做了近似处理，优化模型： 
![image](https://user-images.githubusercontent.com/57510093/211391564-37dbcd9b-5866-49c6-9dfe-36861598c6f7.png)
### Microscopic_siting_linecombine.m
认为多机尾流效应损失为单机尾流效应损失的线性和，优化模型：  
![image](https://user-images.githubusercontent.com/57510093/211392493-b6070387-bdb9-46dd-9d27-b9028d955cda.png) 
### Microscopic_siting_quarcombine.m
认为多机尾流效应损失为单机尾流效应损失的二次和，优化模型：  
![image](https://user-images.githubusercontent.com/57510093/211392341-d4891f8b-c5c0-4398-a9a0-eaab6defccd2.png)  
其他的txt文件为运算结果的储存文件
