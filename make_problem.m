function [alpha,v,R] = make_problem(kind, varargin)

switch kind
    case 'GleichR2'
        R = [
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1/2;
            0 0 0 0 0 1 0 1 0 1/2 0 0 0 1/2 0 0;
            0 0 0 0 0 0 1 0 0 1/2 1 1 0 0 0 0;
            1 1 1 1 1 0 0 0 1 0 0 0 1 1/2 1 1/2
            ];
        v = ones(4,1)/4;
        if length(varargin) < 1
            alpha = 0.9;
        else
            alpha = varargin{1};
        end
    otherwise 
        error('wrong kind');        
end
n = size(R,1);
e = ones(n,1);
assert(all(abs(e'*R-1)<sqrt(eps)));
assert(abs(e'*v-1) < sqrt(eps));
