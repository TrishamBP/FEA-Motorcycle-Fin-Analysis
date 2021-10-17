%% Applied Finite Element Method.
%% HW-3. Trisham Bharat Patil. 688519827.
%% 2) 2 gauss points
clc 
clear all
% Preliminary Given Data
D=0.00635; %Thickness of fin at x=0 in meters(m).
alpha=0.0174533; %Taper angle of fin in radians.
K=175;H=120;Length=0.075; X_start = 0; X_end =Length;

% Reading from the given text file
[fnam, path] = uigetfile('*.txt','InputFinite Element Mesh file');
filename = sprintf('%s%s', path, fnam);
fid = fopen(filename);
nn = fscanf(fid,'%d',1);
ne = fscanf(fid,'%d',1);
nBC = fscanf(fid,'%d',3);
BC1 = zeros(nBC(1),2);     
BC2 = zeros(nBC(2),2);
BC3 = zeros(nBC(3),3);

% Making some array space.
x = zeros(nn,1);  y = x;
El = zeros(ne,2);
ds = (X_end -X_start)/(nn-1);
Lhs = zeros(nn,nn);
Rhs= zeros(nn,1);

% Creating our finite element mesh.
for i = 1:nn
    x(i) = (i-1)*ds + X_start;
end
for L = 1:ne
    El(L,1) = L;
    El(L,2) = L+1;
end
ngp=input('How many gauss points do you want? ');
[Xi,W_GP]=GPoints(ngp);
% Calculate Lhs(i,j), Rhs(i)
for L=1:ne
    xL(1)=x(El(L,1)); xL(2)=x(El(L,2)); %a and b
    dx = x(El(L,2))-x(El(L,1));
    for g=1:ngp
    [N,DNDX,DJ]=gp1d(xL,Xi(g));
    W=N; DWDX=DNDX;
    for i = 1:2
        iG = El(L,i);% Need a Global Mapping
        for j = 1:2
            jG = El(L,j);% Same rationale as with “i” 
            Term_1=DWDX(i)*DNDX(j)*(D-(2*alpha*Xi(g)));
            Term_2=W(i)*N(j);
            Lhs(iG,jG) = Lhs(iG,jG)+(Term_1+Term_2*((2*H)/K))*DJ*W_GP(g);
        end
        Rhs(iG) = Rhs(iG) + 0;
    end
    end
end
% Printing out the output to a text file.
if nn == 5
    fid = fopen('HW3Output.txt','w');
    fprintf(fid,'The Lhs Matrix and Rhs Vector \n');
    fprintf(fid,'Lhs(i,1) Lhs(i,2) Lhs(i,3) Lhs(i,4) Lhs(i,5) Rhs(i)\n');
    for i = 1:nn
        fprintf(fid,'%8f %8f %8f %8f %8f %8f \n',Lhs(i,:),Rhs(i));
    end
end
%% 3).
% Applying the type 2 BC.
Bdy_Loc= nn;
Rhs(Bdy_Loc) = Rhs(Bdy_Loc) -( (1)*((D-(2*alpha*Length)))*5);
% Printing out the output to a text file.
if nn == 5
    fid = fopen('HW3Output_2.txt','w');
    fprintf(fid,'The Lhs Matrix and Rhs Vector after applying BC 1 \n');
    fprintf(fid,'Lhs(i,1) Lhs(i,2) Lhs(i,3) Lhs(i,4) Lhs(i,5) Rhs(i)\n');
    for i = 1:nn
        fprintf(fid,'%8f %8f %8f %8f %8f %8f \n',Lhs(i,:),Rhs(i));
    end
end
%% 4).
% Applying type 1 BC.
Bdy_Loc = 1; 
Lhs(Bdy_Loc,:) = 0;  
Lhs(Bdy_Loc,Bdy_Loc) = 1.0;
Rhs(Bdy_Loc) = 150;
if nn == 5
    fid = fopen('HW3Output_3.txt','w');
    fprintf(fid,'The Lhs Matrix and Rhs Vector after applying BC 1 and BC 2 \n');
    fprintf(fid,'Lhs(i,1) Lhs(i,2) Lhs(i,3) Lhs(i,4) Lhs(i,5) Rhs(i)\n');
    for i = 1:nn
        fprintf(fid,'%8f %8f %8f %8f %8f %8f \n',Lhs(i,:),Rhs(i));
    end
end

%% 5). Solving using “MCycle_Fin_20.txt” file.
% Calling the matrix solver.
temperature = Lhs\Rhs; %The solution.
% Display temperature solution
plot(x,temperature)
xlabel('Distance along fin')
ylabel('Temperature')
title('Thin fin heat conduction problem')

%% 6). Finding the temperature at the tip.
temp_at_tip=temperature(nn,1) %% T=89.5704 degC


function [N,DNDX,DJ] = gp1d(XLocal,XI)
% GAUSS POINT 1 DIMENSIONAL LINEAR Elements
% THIS PROGRAM WAS LAST UPDATED ON 10-27-83
% BY JOHN M. SULLIVAN, JR.   THE PROGRAM PERFORMS
% GAUSS QUADRATURE. LINEAR IN (XI) AND (ETA) DIRECTIONS
%
%      %1        2%
%      *----------*
%
% STANDARD SERINDEPITY FINITE FORMULATIONS
N = zeros(2,1);  DNDX = zeros(2,1);
%
% BASIS FUNCTIONS FOR UNKNOWNS (I.E. TEMPERATURE)
Range = (XLocal(2)-XLocal(1));
Slope = Range/2;  DJ = Slope;
Intercept = (XLocal(2)+XLocal(1))/2;
X = Slope*XI + Intercept;
N(1) = (XLocal(2)-X)/Range;
N(2) = (X - XLocal(1))/Range;
%
DNDX(1) = -1/Range;  DNDX(2) = 1/Range;

end
function [X, W]= GPoints(n)
%% This function receives the number of desired Gauss Points and
%   returns the specific locations "X" and weights "W"
%   John Sullivan  2014

X = zeros(n,1); W=zeros(n,1);
switch n
    case 1
        W(1) = 	2.0;
        X(1) = 0.0;
    case 2
        W(1) = 1; W(2) = 1;
        X(1) = -0.5773502691896257; X(2) = -X(1);
    case 3
        W(2) = 	0.8888888888888888;
        X(2) = 0.0000;
        W(1) = 0.5555555555555556; W(3)=W(1);
        X(1) = -0.7745966692414834; X(3) = -X(1);
    case 4
        W(1) =	0.6521451548625461; W(2) = W(1);
        X(1) = -0.3399810435848563; X(2) = -X(1);
        W(3) = 	0.3478548451374538; W(4) = W(3);
        X(3) = -0.8611363115940526; X(4) = -X(3);
    case 5
        W(3) =	0.5688888888888889;
        X(3) = 0.0000000000000000;
        W(2) =	0.4786286704993665; W(4) = W(2);
        X(2) = -0.5384693101056831; X(4) = -X(2);
        W(1) = 0.2369268850561891;  W(5) = W(1);
        X(1) = -0.9061798459386640; X(5) = -X(1);
    case 6
        W(1)= 0.3607615730481386; X(1)= 0.6612093864662645;
        W(2)= 0.3607615730481386; X(2)= -0.6612093864662645;
        W(3)= 0.4679139345726910; X(3)= -0.2386191860831969;
        W(4)= 0.4679139345726910; X(4)= 0.2386191860831969;
        W(5)= 0.1713244923791704; X(5)= -0.9324695142031521;
        W(6)= 0.1713244923791704; X(6)= 0.9324695142031521;
    case 7
        W(1)= 0.4179591836734694; X(1)=	0.0000000000000000;
        W(2)= 0.3818300505051189; X(2)=	0.4058451513773972;
        W(3)= 0.3818300505051189; X(3)=	-0.4058451513773972;
        W(4)= 0.2797053914892766; X(4)=	-0.7415311855993945;
        W(5)= 0.2797053914892766; X(5)=	0.7415311855993945;
        W(6)= 0.1294849661688697; X(6)=	-0.9491079123427585;
        W(7)= 0.1294849661688697; X(7)=	0.9491079123427585;
    case 8
        W(1) = 0.3626837833783620; X(1)= -0.1834346424956498;
        W(2) = 0.3626837833783620; X(2)= 0.1834346424956498;
        W(3) = 0.3137066458778873; X(3)= -0.5255324099163290;
        W(4) = 0.3137066458778873; X(4)= 0.5255324099163290;
        W(5) = 0.2223810344533745; X(5)= -0.7966664774136267;
        W(6) = 0.2223810344533745; X(6)= 0.7966664774136267;
        W(7) = 0.1012285362903763; X(7)= -0.9602898564975363;
        W(8) = 0.1012285362903763; X(8)= 0.9602898564975363;
    case 9
        W(1) = 0.3302393550012598; X(1) = 0.0000000000000000;
        W(2) = 0.1806481606948574; X(2) = -0.8360311073266358;
        W(3) = 0.1806481606948574; X(3) = 0.8360311073266358;
        W(4) = 0.0812743883615744; X(4) = -0.9681602395076261;
        W(5) = 0.0812743883615744; X(5) = 0.9681602395076261;
        W(6) = 0.3123470770400029; X(6) = -0.3242534234038089;
        W(7) = 0.3123470770400029; X(7) = 0.3242534234038089;
        W(8) = 0.2606106964029354; X(8) = -0.6133714327005904;
        W(9) = 0.2606106964029354; X(9) = 0.6133714327005904;
    case 10
        W(1) = 0.2955242247147529; X(1) = -0.1488743389816312;
        W(2) = 0.2955242247147529; X(2) = 0.1488743389816312;
        W(3) = 0.2692667193099963; X(3) = -0.4333953941292472;
        W(4) = 0.2692667193099963; X(4) = 0.4333953941292472;
        W(5) = 0.2190863625159820; X(5) = -0.6794095682990244;
        W(6) = 0.2190863625159820; X(6) = 0.6794095682990244;
        W(7) = 0.1494513491505806; X(7) = -0.8650633666889845;
        W(8) = 0.1494513491505806; X(8) = 0.8650633666889845;
        W(9) = 0.0666713443086881; X(9) = -0.9739065285171717;
        W(10) = 0.0666713443086881; X(10) = 0.9739065285171717;
    case 11
        W(1) = 0.2729250867779006; X(1) = 0.0000000000000000;
        W(2) = 0.2628045445102467; X(2) = -0.2695431559523450;
        W(3) = 0.2628045445102467; X(3) = 0.2695431559523450;
        W(4) = 0.2331937645919905; X(4) = -0.5190961292068118;
        W(5) = 0.2331937645919905; X(5) = 0.5190961292068118;
        W(6) = 0.1862902109277343; X(6) = -0.7301520055740494;
        W(7) = 0.1862902109277343; X(7) = 0.7301520055740494;
        W(8) = 0.1255803694649046; X(8) = -0.8870625997680953;
        W(9) = 0.1255803694649046; X(9) = 0.8870625997680953;
        W(10) = 0.0556685671161737; X(10) = -0.9782286581460570;
        W(11) = 0.0556685671161737; X(11) = 0.9782286581460570;
    otherwise
        warning('Invalid passing parameters');
end
end