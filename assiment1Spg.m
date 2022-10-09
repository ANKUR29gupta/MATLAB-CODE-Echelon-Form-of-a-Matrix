% assignment 1 of advanced signal processing  
% @@date@@ 13-10-2020
% submitted by-- ANkur gupta -- roll no- 220EE1282
%submitted to--  Dr. Supratim Gupta

% =============================================================================================
% assignment task--Write a function in MATLAB to get a general solution vector (x) a 
% set of linear equations given a matrix (A) of size MxN and of any rank and a vector
% b = Ax.  Write proper comments along your code. You should not use in-built functions
% of MATLAB at any stage of your code and no two 
% codes must be identical in terms of parameters' name or 
% style otherwise negative mark will be applied (50% of total points) to
% all answers those match. You need to submit the following m-files (*.m)
%  The function like x = my_func(A,b)
%  A script file demonstarting the use of the function  like A1 = []; b1 = []; x1 = my_func(A1,b1)
% =================================================================================================


clear;
clc;
% ==================================================================================================
% % matrix for cheking code uncomment them and use it
% % (1) row < colum
% A=[1 3 3 2  ; 2 6 9 7;-1 -3 3 4];
% B=[1;5;5];
% % (2)row <colum
% A=[1 2 3 5  ; 2 4 8 12;3 6 7 13]
% B=[0;6;-6]
% % (3) row =colum
% A=[1 2 -5  ; 3 -1 2;2 3 -1]
% B=[-9;5;3]
% % (4)row =colum
% A=[2 1 2 1  ; 6 -6 6 12;4 3 3 -3;2 2 -1 1]
% B=[6;36;-1;10]
% % (5)row>colum
% A=[1 2;2 4;3 6;4 8];
% B=[5;10;15;20];
% % (5)row>colum
% A=[3 -6;4 -8;0 1];
% B=[-1;7;2];
% A=[3 1 3;2 4 1];
% B=[3;4];
A=[1 3 3 2;2 6 9 7;-1 -3 3 4]
B=[1;5;5];
% =====================================================================================================
%this genrate my solution of Xp and Xn
[xpsol,xnsol]=my_funtion(A,B);
[p,t,rankofA,rankofAB]=myEchelon(A,B);
% =====================================================================================================
% this is funtion which solve X value
function [Xp,Xn]=my_funtion(k,t)
[p,t,rankofA,rankofAB]=myEchelon(k,t); %find echelon form
Xp=findXp(p,t,rankofA,rankofAB); %find xp
Xn=findXn(p,rankofA);  %find xn
XpValue=num2str(Xp);   %number to string convert for display 
XnValue=num2str(Xn);
disp('=====================================================================================')
disp('Xp =')
disp(XpValue)
disp('=====================================================================================')
disp('Xn =')
disp(XnValue)
disp('=====================================================================================')
end
% ===========================================================================================
% for finding xp
function xp=findXp(p,t,rankofA,rankofAB)
[Prow,Pcol]=size(p);   %here find size of p matrix
if Prow>Pcol           % here check row >colum and covert this matrix into colum*colum matrix 
    p=p(1:Pcol,:);     % for solving purpose
    t=t(1:Pcol,1);
    Prow=Pcol;
end
xp=ones(Pcol,1);      % here intial creat Xp values of ones matrix
Xnfree=zeros(Pcol,1); % this is use for checking which position have free variable
K=horzcat(p,t);        % this funtion use for joining two matrix in horizontal form


for q=1:Prow          %in these loop i am creat first non zero values of row is one 
    for w=1:Pcol
        if(K(q,w)~=0)
            K(q,:)= K(q,:)/K(q,w);
            break;
        end
    end
end
xp(end,1)=K(Prow,Pcol+1); % here i update end value of xp
freeVariable=Pcol-rankofA; % how many free variable is there
for h=1:freeVariable
    xp(Pcol-(2*h-2),1)=0;     % in these loop update the value of xp for free variable
    Xnfree(Pcol-(2*h-2),1)=1; % here fill the position where free variable is filled
end
% here apply back subsitution for calcultaing xp
valueofX=0;
No=Pcol;
for roW=Prow-1:-1:1
    I=0;
    S=0;
    for s=1:Pcol
        
        if K(roW,s)==1
            S=s;             % this is for how many time i run the loop
            break;
        end
    end
    
    for coL=Pcol:-1:S+1
        valueofX=valueofX+K(roW,coL)*xp(Pcol-I,1);  %this for calculating xp
        I=I+1;
    end
    if freeVariable~=0          % for skiping xp value which is free variable
        for Ho=No-1:-1:1
            if Xnfree(Ho,1)==0
                J=Ho;
                break;
            end
        end
        xp(Ho,1)=K(roW,end)-valueofX;    % these line for updating xp
        No=No-1;
    else
        xp(roW,1)=K(roW,end)-valueofX;
    end
    valueofX=0;
end
end
%=============================================================================================

% finding xn
function xn=findXn(p,rankofA)
[Prow,Pcol]=size(p); %here find size of p matrix
t=zeros(Prow,1);     % here intial creat t values of zeros matrix
if Prow>Pcol
    p=p(1:Pcol,:);  % here check row >colum and covert this matrix into colum*colum matrix 
    t=t(1:Pcol,1);  % for solving purpose
    Prow=Pcol;
end
xn=ones(Pcol,1);        % here intial creat Xn values of ones matrix
Xnfree=zeros(Pcol,1);  % this is use for checking which position have free variable
K=horzcat(p,t);        % this funtion use for joining two matrix in horizontal form
for q=1:Prow
    for w=1:Pcol        %in these loop i am creat first non zero values of row of K is one 
        if(K(q,w)~=0)
            K(q,:)= K(q,:)/K(q,w);
            break;
        end
    end
end


freeVariable=Pcol-rankofA;  % how many free variable is there
if freeVariable~=0
    for h=1:freeVariable   % in these loop update the value of xn for free variable
        xn(Pcol-(2*h-2),1)=input('please give a value for infinity many solution(for freevariable)');
        Xnfree(Pcol-(2*h-2),1)=1; % here fill the position where free variable is filled
    end
    % here apply back subsitution for calcultaing xn
    valueofX=0;
    No=Pcol;
    for roW=Prow-1:-1:1
        I=0;
        S=0;
        for s=1:Pcol
            
            if K(roW,s)==1
                S=s;            % this is for how many time i run the loop
                break;
            end
        end
        
        for coL=Pcol:-1:S+1
            valueofX=valueofX+K(roW,coL)*xn(Pcol-I,1);  %this for calculating xn
            I=I+1;
        end
        if freeVariable~=0
            for Ho=No-1:-1:1
                if Xnfree(Ho,1)==0      % for skiping xn value which is free variable
                    J=Ho;
                    break;
                end
            end
            xn(Ho,1)=K(roW,end)-valueofX;    % these line for updating xn
            No=No-1;
        else
            xn(roW+1,1)=K(roW,end)-valueofX;
        end
        valueofX=0;
    end
else
    xn=zeros(Pcol,1);
end
end
%=========================================================================================
% find echelon for
function [Unew,Bnew,RAnkforA,RAnkforAB]=myEchelon(U,B)
A=horzcat(U,B); 
[Arow,Acol]=size(A);
while(1)
    for n=1:Arow-1
        if A(n,n)~=0
            
            for c=n:Arow-1
                A(c+1,:)= A(c+1,:)-(A(c+1,n)/A(n,n))*A(n,:); % row operation conducted
                
            end
        elseif A(n+1,n)~=0
          
            A([n n+1],:)=A([n+1 n],:); % row swipping 
            
            break;
        else
            
            if n<min(Arow,Acol)
                if A(n,n+1)==0&&A(n+1,n+1)==0
                    continue;
                end
                A(n+1,:)=A(n+1,:)-(A(n+1,n+1)/A(n,n+1))*A(n,:); % for row operation when pivot is zero
                break;
            end
            
        end
        
    end
    if A(Arow,Acol-2)==0  % it is breaking for loop after completeing process
        break;
    end
end

Unew=A(:,1:Acol-1);  %for take out A after redution operation of A|B matrix (Ax=B)
Bnew=A(:,Acol);      %for take out B after redution operation of  A|B matrix
%for checking how many is not equal to zero B=A~=0
% for finding rank A|B
B=A~=0;
RAnkforAB=0;
%here find how many row is nonzero row
for j=1:Arow
    if (sum(B(j,:))~=0)
        RAnkforAB=RAnkforAB+1; 
    end
end
% for finding rank A
D=Unew~=0;
RAnkforA=0;
%here find how many row is nonzero row
for j=1:Arow
    if (sum(D(j,:))~=0)
        RAnkforA=RAnkforA+1;
    end
end

end
