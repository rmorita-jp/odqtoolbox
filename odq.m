function [Q E Hk gain] = odq(G,T,d,gamma,dim,solver)
%ODQ Optimal Dynamic Quantizer
%
%Q = odq(G,T,d) is ODQ for G with evaluation interval 'T'.
%Q = odq(G,T,d,gamma) is ODQ for G whose upper bound of gain is gamma, and
%quantizing interval is d.
%(gamma is structure which contains uv and wv.)
%Q = odq(G,T,d,gamma,dim) can specify the dimension of Q. dim is the
%dimension.
%[Q E] = odq(G,T,d,gamma,dim) shows Q's cost.
%[Q E Hk] = odq(G,T,d,gamma,dim) keeps Hankel matrix to workspace.
%
%Q = odq(G,T,d,gamma,dim,solver) specifies optimization solver.
%'linprog' uses linprog. (Optimization Toolbox by Mathworks)
%'cplex' uses CPLEX. (by ILOG http://www.ilog.com/)
%CPLEX needs CPLEX MEX INTERFACE (contains some bugs)
%
%See also compg, odqreal, odqgain, odqcost, odqstb.

% -------------------------------------------------------------------------
% Copyright is with the following author. 
% (C) 2008 Ryosuke Morita, 
%          Kyoto University;
%          Gokasho, Uji, Kyoto 611-0011, Japan
%          morita@robot.kuass.kyoto-u.ac.jp
% -------------------------------------------------------------------------
% Legal Note:
%           
%     (a)  This program is a free software. 
%          
%     (b)  This program is distributed according to GNU General Public
%          License, i.e., it is allowed to use WITHOUT ANY WARRANTY; 
%          without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%          See the GNU General Public License for more details.
% 
% -------------------------------------------------------------------------

if nargin==5
    if ischar(dim)==1
        solver = dim;
        dim=0;
    else
        solver='linprog';
    end
elseif nargin==4
    if ischar(gamma)==1
        solver   = gamma;
        clear gamma
        gamma.wv = inf;
    else
        solver='linprog';
    end
    dim=0;
elseif nargin==3
    gamma.wv = inf;
    dim=0;
    solver='linprog';
elseif nargin<3
    error('there are not enough input argument')
end

%%%%%set dimension%%%%%
m = size(G.c1*G.b2,2);
p = size(G.c1*G.b2,1);
%%%%%%%%%%%%%%%%%%%%%%%

G.aa = G.a+G.b2*G.c2;    %convert to closed loop

%%%%%%%%%%Compose Quantizer%%%%%%%%%%
endflag=0;
if T~=inf
    mstkn=1;
else
    mstkn=0;
end
while (endflag == 0)
    if (m >= p && mstkn == 0)
        Q = odqanly(G);
        %Hk = 'null';
        Hk=odqhnkl_anly(Q);
        if max(abs(eig(Q.a+Q.b2*Q.c)))<=1
            endflag=1;
        else
            mstkn = 0;
            endflag=1;
            disp('Try optimization...')
        end
    else
        if strcmp(solver,'linprog')==1
            [x exitflag,Tnew,fval] = odqlp(G,T,gamma);
        elseif strcmp(solver,'cplex')==1
            [x exitflag,Tnew,fval] = odqlp_cplex(G,T,gamma);
        elseif strcmp(solver,'sdpt3')==1
            [x exitflag,Tnew,fval] = odqlp_sdpt3(G,T,gamma,solver);
        elseif strcmp(solver,'sedumi')==1
            [x exitflag,Tnew,fval] = odqlp_sedumi(G,T,gamma,solver);
        else
            error('wrong solver')
        end
        if (exitflag==1)
            Hk = odqhnkl(x,Tnew,m,p);
            if dim~=0
                [Q Hk] = odqreal(G,Hk,dim);
            else
                [Q Hk] = odqreal(G,Hk);
            end
            endflag=2;
        else
            error('Optimization is failed. Please change T and gamma')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Check Quantizer%%%%%
if endflag==1
    E=odqcost(G,Q,d,T);
elseif endflag==2
    E=fval*d/2;
end
if (size(Q.a,1)>1000)
    disp('Q has too big dimension. Skip verifying')
    skipflag=1;
    gain.wv='skipped';
else
    gain=odqgain(Q,T);
    skipflag=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Display Result%%%%%
fprintf('\n%s\n','RESULT')
disp('Q')
disp(Q)
if (skipflag==1)
    fprintf('%s%d\n%s\n','T = ',T,'Computing gain and cost are skiped')
else
    fprintf('%s%d\n%s%d\n%s%d\n','T = ',T,'gain uv = ',gain.uv,'gain wv = ',gain.wv)
    if E~=0
        fprintf('%s%1.4e%s%d%s\n','E = ',E,' in ',T,' steps.') 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%