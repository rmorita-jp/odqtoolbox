function [x,exitflag,T,FOPT] = odqlp_sedumi(G,T,gamma,solver)
%ODQLP Optimal Dynamic Quantizer with LP
%
%This function is not to use alone.
%Please use 'ODQ'.
%
%See also odq.

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

%%%%%set parameters%%%%%
if (mod(T,1)~=0)
    error('steptime must be a natural number');
end
if (mod(T,2)==0)
    T=T+1;
end
if (gamma.wv<1)
    error('max gain w->v must be greater than 1');
end

m = size(G.c1*G.b2,2);
p = size(G.c1*G.b2,1);

G.aa = G.a+G.b2*G.c2;    %convert to closed loop
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%initialize%%%%%
if gamma.wv==inf
    vectorf = zeros( 1 + m*m*T + p*m*(T-1) , 1 );
else
    vectorf = sparse(zeros( 1 + m*m*T*2 + p*m*(T-1) + p + p*m*(T-1)*2 + m + m*m*T*2 , 1 ));
end
%vectorb = zeros( p + p*m*(T-1)*2 + m + m*m*T*2 , 1 );
%matrixA = zeros( p + p*m*(T-1)*2 + m + m*m*T*2 , 1 + m*m*T + p*m*(T-1) );
%%%%%%%%%%%%%%%%%%%%

vectorf(1)=1;

fprintf('Preprocessing for LP...\n')
pause(0.1)

%%%%%making PHI%%%%%
PHI     = zeros( m*p*(T-1) , m*m*T );   %initialize
% PHIdash = zeros( p , m );               %initialize
for i=1:T-1
    PHIdash = kron( G.c1*G.aa^(i-1)*G.b2 , eye(m) );
    for j=i:T-1
        PHI( (j-1)*m*p+1:j*m*p , (j-i)*m*m+1:(j-i+1)*m*m ) = PHIdash;
    end
end
%%%%%%%%%%%%%%%%%%%%

%%%%%making 'matrix -> vector' transposer%%%%%
eye_sumE = kron( ones(1,m*(T-1)) , eye(p) );
eye_sumH = kron( ones(1,m*T) , eye(m) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%finalize%%%%%%%%%%
if gamma.wv==inf
    matrixA = [ -ones( p , 1 )         zeros( p , m*m*T )  eye_sumE                 ; % 
                zeros( p*m*(T-1) , 1 )   PHI              -eye( m*p*(T-1) )         ;
                zeros( p*m*(T-1) , 1 )  -PHI              -eye( m*p*(T-1) )         ];
else
    matrixA = sparse([ -ones( p , 1 )         zeros( p , m*m*T ) zeros( p , m*m*T )           eye_sumE                             eye(p) zeros(p        ,p*m*(T-1)) zeros(p        ,p*m*(T-1)) zeros(p        ,m) zeros(p        ,m*m*T) zeros(p        ,m*m*T);  % sum(E)-G <-abs(CB) 
                       zeros( p*m*(T-1) , 1 )   PHI              zeros( p*m*(T-1) , m*m*T )  -eye( m*p*(T-1) )         zeros(p*m*(T-1),p)             eye(p*m*(T-1)) zeros(p*m*(T-1),p*m*(T-1)) zeros(p*m*(T-1),m) zeros(p*m*(T-1),m*m*T) zeros(p*m*(T-1),m*m*T);  %  PHI*H-E  <-C*A^k*B
                       zeros( p*m*(T-1) , 1 )  -PHI              zeros( p*m*(T-1) , m*m*T )  -eye( m*p*(T-1) )         zeros(p*m*(T-1),p) zeros(p*m*(T-1),p*m*(T-1))             eye(p*m*(T-1)) zeros(p*m*(T-1),m) zeros(p*m*(T-1),m*m*T) zeros(p*m*(T-1),m*m*T);  % -PHI*H-E <C*A^k*B
                       zeros( m , 1 )         zeros( m , m*m*T )   eye_sumH                 zeros( m , m*p*(T-1) )     zeros(m        ,p) zeros(m        ,p*m*(T-1)) zeros(m        ,p*m*(T-1))             eye(m) zeros(m        ,m*m*T) zeros(m        ,m*m*T);  % sum(Hbar)<gamma-1
                       zeros( m*m*T , 1 )       eye( m*m*T )      -eye( m*m*T )             zeros( m*m*T , m*p*(T-1) ) zeros(m*m*T    ,p) zeros(m*m*T    ,p*m*(T-1)) zeros(m*m*T    ,p*m*(T-1)) zeros(m*m*T    ,m)             eye(m*m*T) zeros(m*m*T    ,m*m*T);  % H-Hbar   <0
                       zeros( m*m*T , 1 )      -eye( m*m*T )      -eye( m*m*T )             zeros( m*m*T , m*p*(T-1) ) zeros(m*m*T    ,p) zeros(m*m*T    ,p*m*(T-1)) zeros(m*m*T    ,p*m*(T-1)) zeros(m*m*T    ,m) zeros(m*m*T    ,m*m*T)             eye(m*m*T)]);% -H-Hbar  <0
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%making C A^k B ,changing matrix to vector%%%%%
CAB = zeros( (T-1)*p , m );             %initialize
el_CAB = zeros( m*p*(T-1) , 1 );        %initialize

for i=1:T-1
    CAkB = G.c1*G.aa^i*G.b2;
    CAB( (i-1)*p+1:i*p , 1:m ) = CAkB;
end

for j=1:p*(T-1)
    for i=1:m
        el_CAB( i+(j-1)*m , 1 ) = CAB( j , i );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%finalize%%%%%%%%%%
if gamma.wv==inf
    vectorb = [ -abs( G.c1*G.b2 )*ones( m , 1 ) ;
                -el_CAB                         ;
                 el_CAB                         ];
else
    vectorb = sparse([ -abs( G.c1*G.b2 )*ones( m , 1 ) ;
                -el_CAB                         ;
                 el_CAB                         ;
                 (gamma.wv-1)*ones( m , 1 )     ;
                 zeros( m*m*T*2,1 )             ]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%liner programing problem%%%%%%%%%%
% switch solver
%     case 'sdpt3'
%         cvx_solver sdpt3
%         disp('Solving LP problem(SDPT3)...')
%     case 'sedumi'
%         cvx_solver sedumi
%         disp('Solving LP problem(SeDuMi)...')
%     otherwise
%         error('wrong solver')
% end
% %length(vectorf)
% pause(0.1)
% cvx_begin
%     variable x(length(vectorf));
%     minimize(vectorf'*x);
%     subject to
%         matrixA*x <= vectorb;
% cvx_end
% exitflag=1;
% FOPT=vectorf'*x;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%liner programing problem%%%%%%%%%%
disp('Solving LP problem(SeDuMi)...')
K.f=1 + m*m*T;
K.l= m*m*T + p*m*(T-1) + p + p*m*(T-1)*2 + m + m*m*T*2;
x=sedumi(matrixA,vectorb,vectorf,K);
length(vectorf)
exitflag=1;
FOPT=vectorf'*x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
