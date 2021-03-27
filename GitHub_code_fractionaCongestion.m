% Fractional congestion dynamics.
%  It assumes a binary/weighted and dir/un-directed matrix A
%  It requires of the ml_matrix.m code for computing the Caputo Derivative

clearvars; clc; close all;

A = rand(20,20);        
A = A - diag(diag(A));  % just as simple example of a directed weighted graph
%%     fractional SI model

n = length(A);
fexponent=0.75;
beta = 0.01;
c = 0.005;
gamma = 1-c/n;


tmax=100;
for y=1:tmax
    t=(y-1);
    
    G(:,:,y) = ml_matrix(t^(fexponent)*beta^(fexponent)*gamma*A,fexponent,1);
    
    v1(:,y) = (1/gamma-1)*G(:,:,y)*ones(n,1) - (1/gamma-1+log(gamma))*ones(n,1);    % Export analysis
    v2(:,y) = (1/gamma-1)*ones(1,n)*G(:,:,y) - (1/gamma-1+log(gamma))*ones(1,n);    % Import analysis
    
    vv1(y) = v1(:,y)'*ones(n,1);    % global congestion for Export analysis
    vv2(y) = v2(:,y)'*ones(n,1);
end
t=0:(tmax-1);

%% %Plotting ALL nodes 
x=(ones(1,tmax)-exp(-v1));
figure
plot(t,x)
axis square
xlabel('Time')
ylabel('fractional congestion dynamics for nodes')
