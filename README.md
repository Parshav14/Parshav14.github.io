%% Algorithm of Simplex Method to solve LPP

% clc
% clear all
% format short
% %To solve the LPP by Simplex Method
% %Min z=x1-3x2+2x3
% %Subject to 3x1-x2+2x3<=7
% %-2x1+4x2<=12
% %-4x1+3x2+8x3<=10
% %x1,x2,x3>=0
% %First to change the objective function from minimization to maximization
% %Max z=-x1+3x2-2x3
% %To input parameters
% [info s b]

% To input Parameters
C = [-1 3 -2];
info = [3 -1 2;-2 4 0;-4 3 8];
b = [7; 12; 10];
NOVariables = size(info,2); % x1,x2,x3
s = eye(size(info,1));
A = [info s b];
Cost = zeros(1,size(A,2));
Cost(1:NOVariables) = C;
BV = NOVariables+1:size(A,2)-1;
%To calculate Z-Row(Zj - Cj)
ZRow = Cost(BV)*A-Cost;
% To print the table
ZjCj = [ZRow; A]; % First Column of Zj - Cj
SimpleTable = array2table(ZjCj);
SimpleTable.Properties.VariableNames(1:size(ZjCj,2))={'x_1','x_2','x_3','s_1','s_2','s_3','b'}
% Simplex Table Starts
Run = true;
while Run
if any(ZRow<0) % To check any negative value is there
fprintf('The current BFS is not optimal\n');
fprintf('\n=======The next iteration continues======\n');
disp('Old Basic Variable(BV)=')
disp(BV);
% To find the entering variable
ZC = ZRow(1:end-1);
[EnterCol, Pvt_Col] = min(ZC);
fprintf('The most negative element in Z-Row is %d Corresponding to Column %d \n',EnterCol,Pvt_Col);
fprintf('Entering Variable is %d \n', Pvt_Col);
% To find the leaving variable
sol = A(:,end);
Column = A(:,Pvt_Col);
if all(Column<=0)
error('LPP has unbounded solution. All entries<=0 in column %d',Pvt_Col)
else
% To check minimum ratio is with positive entering column
% entries
for i=1:size(Column,1)
if Column(i)>0
ratio(i)=sol(i)./Column(i)
else
ratio(i) = inf
end
end
% To find minimum ratio
[MinRatio, Pvt_Row] = min(ratio)
fprintf('Minimum ratio corresponding to pivot row is %d \n', Pvt_Row);
fprintf('Leaving Variable is %d \n', BV(Pvt_Row));
end
BV(Pvt_Row) = Pvt_Col;
disp('Newe Basic Variable (BV) =')
disp(BV)
% Pivot Key
Pvt_key = A(Pvt_Row,Pvt_Col);
% Update table for next iteration
A(Pvt_Row,:)=A(Pvt_Row,:)./Pvt_key
for i=1:size(A,1)
if i~=Pvt_Row
A(i,:)=A(i,:)-A(i,Pvt_Col).*A(Pvt_Row,:)
end
ZRow=ZRow-ZRow(Pvt_Col).*A(Pvt_Row,:)
%To print the updated table
ZjCj=[ZRow;A]
SimpTable=array2table(ZjCj)
SimpTable.Properties.VariableNames(1:size(ZjCj,2))={'x_1','x_2','x_3','s_1','s_2','s_3','Sol'}
BFS=zeros(1,size(A,2))
BFS(BV)=A(:,end)
BFS(end)=sum(BFS.*Cost)
CurrentBFS=array2table(BFS)
CurrentBFS.Properties.VariableNames(1:size(CurrentBFS,2))={'x_1','x_2','x_3','s_1','s_2','s_3','Sol'}
end
else
Run=false;
fprintf('The current BFS is optimal and Optimality is reached \n')
end
end

%% Algorithm of Big M to solve LPP

%% Phase I: Information Phase
Variables={'x_1','x_2','s_2','s_3','A_1','A_2','Sol'};
M=1000;
Cost=[-2,-1,0,0,-M,-M,0]; %cost of the LPP
A=[3, 1, 0, 0, 1, 0, 3;4, 3, -1, 0, 0, 1, 6;1, 2, 0, 1, 0, 0, 3]; %Constraints
s=eye(size(A,1));
%% To find the starting BFS
BV=[];
for j=1:size(s,2)
for i=1:size(A,2)
if A(:,i)==s(:,j)
BV=[BV i];
end
end
end
%% To compute Z-Row(zj-cj)
ZjCj=Cost(BV)*A-Cost;
%% To print the table
ZCj=[ZjCj;A];
SimpTable=array2table(ZCj);
SimpTable.Properties.VariableNames(1:size(ZCj,2))=Variables;
%% Simplex Table starts
Run=true;
while Run
ZC=ZjCj(:,1:end-1);
if any(ZC<0) %To check any negative value is there
fprintf('The current BFS is not Optimal \n');
fprintf('\n ============The Next Iteration Results========\n');
%To find the entering variable
[Entval,Pvt_Col]=min(ZC);
fprintf('Entering Column is %d \n', Pvt_Col);
%To find the leaving variable
sol=A(:,end);
Column=A(:,Pvt_Col);
if all(Column<=0)
fprintf('LPP has unbounded solution');
else
%To check minimum ration is with positive entering column entries
for i=1:size(Column,1)
if Column(i)>0
ratio(i)=sol(i)./Column(i);
else
ratio(i)=inf;
end
end
%To Finding the minimum Ratio
[MinRatio, Pvt_Row]=min(ratio);
fprintf('Leaving Row is %d \n', Pvt_Row);
end
BV(Pvt_Row)=Pvt_Col;
%Pivot Key
Pvt_Key=A(Pvt_Row,Pvt_Col);
%Update the Table for next Iteration
A(Pvt_Row,:)=A(Pvt_Row,:)./Pvt_Key;
for i=1:size(A,1)
if i~=Pvt_Row
A(i,:)=A(i,:)-A(i,Pvt_Col).*A(Pvt_Row,:);
end
ZjCj=ZjCj-ZjCj(Pvt_Col).*A(Pvt_Row,:);
%To print the table
ZCj=[ZjCj;A];
SimpTable=array2table(ZCj);
SimpTable.Properties.VariableNames(1:size(ZCj,2))=Variables;
end  
else
Run=false;
fprintf('==================\n');
fprintf('The current BFS is optimal and Optimality is reached \n');
fprintf('==================\n');
end
end
BFS=zeros(1,size(A,2));
BFS(BV)=A(:,end);
BFS(end)=sum(BFS.*Cost);
CurrentBFS=array2table(BFS);
CurrentBFS.Properties.VariableNames(1:size(CurrentBFS,2))=Variables

%% Algorithm of LCM to solve LPP

clc
clear all
c=[2 1 4; 3 1 2; 5 6 7]
a=[50 150 200]
b=[100 130 200]
z=0
if sum(a) == sum(b)
fprintf('Given Problem is balanced')
else
fprintf('Given problem is unbalanced')
if sum(a)>sum(b)
c(:,end+1)=zeros(length(a),1)
b(end+1)=sum(a)-sum(b)
else
c(end+1,:)=zeros(1,length(b))
a(end+1)=sum(b)-sum(a)
end
end
X=zeros(size(c,1),size(c,2))
initialc = c;
for i=1:size(c,1)
for j=1:size(c,2)
cpq=min(c(:))
if cpq==inf
break
end
[p1,q1]=find(cpq==c)

        possalloc = min(a(p1), b(q1))
        [maxalloc, indx] = max(possalloc)

        p = p1(indx)
        q = q1(indx)

        X(p,q) = min(a(p),b(q))
        if min(a(p),b(q))==a(p)
            b(q)=b(q)-a(p)
            a(p)=a(p)-X(p,q)
            c(p,:)=inf
        else
            a(p)=a(p)-b(q)
            b(q)=b(q)-X(p,q)
            c(:,q)=inf
        end
    end

end
for i=1:size(c,1)
for j=1:size(c,2)
z=z+ initialc(i,j)\*X(i,j);
end
end
array2table(X)
fprintf('transportation cost is %f',z)

%% Algorithm of Steepest Descent Method to solve LPP

syms x y
f1=x-y+2*x^2+2*x*y+y^2
f=inline(f1)
fobj=@(x)f(x(:,1),x(:,2))
grad=gradient(f1)
G=inline(grad)
grad_x=@(x)G(x(:,1),x(:,2))
H1=hessian(f1)
HX=inline(H1)
x_1=[1 1]
max_iteration=4
tol=0.003
iter=0
X=[]
while norm(grad_x(x_1))>tol && iter<max_iteration
X=[X,x_1]
S=-grad_x(x_1)
H=HX(x_1)
lam=(S'*S)./(S'*H*S)
x_new=x_1+lam.\*S'
x_1=x_new
iter=iter+1
end
fprintf('optimal sol X=[%f %f]',x_1(1),x_1(2))
fprintf('optimal value f(x)=%f',fobj(x_1))

<a href="https://colab.research.google.com/drive/1TYKESyV-NJmtT7DHVbzvryK0urJMUWd9?usp=drive_link">Colab Link</a><br><pre>Code here</pre>
