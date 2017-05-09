function [v lambda]=perronvector(P,method,tol,v0);
%
% [v lambda]=perronvector(P,method,tol,v0);
%
% gets the Perron vector (or the maximum real part one) of a positive matrix
% normalized s.t. norm(v,2)=1 and v>0
% v0 is a "guess" of the vector (that may or may not 
% be used by the algorithm)
% tol is the tolerance (that may or may not...)
%
% available methods:
% 'eig'
% 'eigs'
% 'power'
% 'squaring' -> power method + successive squaring acceleration
% 'squaring2'-> same but more aggressive acceleration
% % (c) f.poloni@sns.it 2009-2010

n=size(P,1);
if(not(exist('v0','var')))
    v0=rand(n,1);
end
if(not(all(P>=0)))
  warning('P is not a positive matrix');
end
if(strcmp(method,'eig'))
  [V Lambda]=eig(P);
  [useless,j]=max(real(diag(Lambda)));
  v=V(:,j);
  lambda=Lambda(j,j);
elseif(strcmp(method,'eigs'))
  opts.disp=0;
  opts.tol=tol;
  opts.v0=v0;
  [v lambda]=eigs(P,1,'LR',opts);
elseif(strcmp(method,'power'))
  res=inf;
  v=v0;
  while(res>tol)
    vold=v;
    v=P*v;
    v=v/norm(v,2);
    res=norm(vold-v);
    %disp(sprintf('   inner power it: %g ',res));
  end
  lambda=sum(P*v)/sum(v);
elseif(strcmp(method,'squaring'))
    %hybrid power method + successive quadratures
  res=inf;
  v=v0;
  P1=P;
  while(res>tol)
    vold=v;
    P1=P1*P1;
    v=P1*v;
    v=v/norm(v,2);
    res=norm(vold-v);
    %disp(sprintf('   inner power it: %g ',res));
  end
  lambda=sum(P*v)/sum(v);
elseif(strcmp(method,'squaring2'))
  %more aggressive version
  res=inf;
  v=v0;
  P1=P*P;
  P1=P1*P1;
  while(res>tol)
    vold=v;
    P1=P1*P1;
    v=P1*v;
    v=v/norm(v,2);
    res=norm(vold-v);
    %disp(sprintf('   inner power it: %g ',res));
  end
  lambda=sum(P*v)/sum(v);
else
  error('Unknown method')
end
v=v/v(1); %makes sure it's positive
v=v/norm(v,2); %normalizes
