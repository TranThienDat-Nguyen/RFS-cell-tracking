function R = crr(s,varargin)
% CRR - Confidence Region Radius (R13)
%
%   R = CRR(S) computes the radius of the mean-centered interval, circle,
%   or sphere with 95% probability given S, which is either a vector of
%   standard deviations or a covariance matrix from a multivariate normal
%   distribution.  If S is a real, symmetric, positive semidefinite
%   matrix, CRR(S) is equivalent to CRR(SQRT(EIG(S))).  Scalar S is treated
%   as a standard deviation.
%
%   R = CRR(S,P) computes the confidence region radius with probability P
%   instead of the default, which is 0.95.
%
%   R = CRR(S,P,TOL) uses a quadrature tolerance of TOL instead of the
%   default, which is 1e-15.  Larger values of TOL may result in fewer
%   function evaluations and faster computation, but less accurate results.
%   Use [] as a placeholder to obtain the default value of P.
%
%   R = CRR(S,P,TOL,M) performs a bootstrap validation with M normally
%   distributed random samples of size 1e6.  Use [] as a placeholder to
%   obtain the default value of TOL.
%
%   R = CRR(S,P,TOL,[M N]) performs a bootstrap validation with M normally
%   distributed random samples of size N.
%
% EXAMPLES:
%
%   % Probable Error (PE)
%   r=crr(1,0.5) % = 0.6745
%
%   % 95% confidence interval
%   r=crr(1) % = 1.9600
%
%   % 90% confidence interval from data vector X
%   r=crr(std(X),0.9)
%
%   % Circular Error Probable (CEP)
%   r=crr([1 2],0.5) % = 1.7408
%
%   % 95% confidence circle with bootstrap validation
%   r=crr([1 2],[],[],[10 1e5]) % = 4.0717
%
%   % 90% confidence circle from covariance matrix V
%   V=[1 1; 1 4];
%   r=crr(V,0.9) % = 3.5263
%
%   % 95% confidence circle from data matrix [X Y]
%   r=crr(cov([X Y]))
%
%   % Spherical Error Probable (SEP)
%   r=crr([1 2 3],0.5) % = 3.1068
%
%   % 95% confidence sphere from covariance matrix V
%   V=[1 1 2; 1 4 -1; 2 -1 9];
%   r=crr(V) % = 6.5817
%
%   % 90% confidence sphere from data matrix [X Y Z]
%   r=crr(cov([X Y Z]),0.9)
%
%   % Inverse chi-square
%   % chi2inv(p,n) = crr(eye(n),p)^2 = crr(ones(1,n),p)^2;  n=1,2,3
%   r2=crr(1)^2       % = 3.8415 = chi2inv(0.95,1)
%   r2=crr([1 1])^2   % = 5.9915 = chi2inv(0.95,2)
%   r2=crr([1 1 1])^2 % = 7.8147 = chi2inv(0.95,3)
%
% REFERENCES:
%
%   Abramowitz, M. and Stegun, I. A., eds. (A&S 1972). Handbook of
%     mathematical functions with formulas, graphs, and mathematical
%     tables. National Bureau of Standards. 6.5.1, 7.1.23, 7.1.29, 9.6.12.
%     (http://www.math.sfu.ca/~cbm/aands/frameindex.htm)
%   Acklam, P. J. (2002). "Variance matrix checking: varchk."
%     (http://home.online.no/~pjacklam/matlab/software/util/statutil/
%     varchk.m)
%   Altham, P. M. E. (2006). "Applied multivariate analysis, notes."
%     (http://www.statslab.cam.ac.uk/~pat/AppMultNotes.pdf)
%   Davis, T. and Kleder, M. (2006). "Confidence region radius."
%     Mathworks Central File Exchange.
%     (http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?
%     objectId=10526&objectType=file)
%   Godfrey, P. and Acklam, P. J. (2003). "Error function for complex
%     inputs: erfz."  Mathworks Central File Exchange.
%     (http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?
%     objectId=3574&objectType=file)
%   Kleder, M. (2004). "An algorithm for converting covariance to
%     spherical error probable." Mathworks Central File Exchange.
%     (http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?
%     objectId=5688&objectType=file)
%   Koopman, R. (2004). "Joint normal distribution integral over the unit
%     sphere." The Math Forum.
%     (http://mathforum.org/kb/message.jspa?messageID=1552353&tstart=0)
%   National Institute of Standards and Technology (NIST 2006).
%     "The multivariate normal distribution."
%     (http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc542.htm)
%   Sloane, N. J. A. (2006). "The on-line encyclopedia of integer
%     sequences." AT&T Labs Research. A000984
%     (http://www.research.att.com/~njas/sequences/A000984)
%
% See also crr.pdf.
% Copyright(c)2006-2008
%   Tom Davis (tdavis@metzgerwillard.com)
%   Michael Kleder (mkleder@hotmail.com)
%
% Last revision: 03/16/2008
if nargin<1 || isempty(s)
  error('CRR requires at least one argument.')
end
if all(size(s)>1)                                    % covariance matrix
  if ~isa(s,'double')
    error('Covariance matrix must be of class ''double''.')
  elseif ndims(s)>2
    error('Covariance matrix cannot have more than two dimensions.')
  elseif size(s,1)~=size(s,2);
    error('Covariance matrix must be square.')
  elseif length(s)>3;
    error('Covariance matrix cannot have more than three rows.')
  elseif ~isreal(s)
    error('Covariance matrix must be real.')
  elseif any(any(s~=s.'))
    error('Covariance matrix must be symmetric.')
  end
  s=eig(s);
  if any(s<-eps)
    error('Covariance matrix must be positive semidefinite.')
  end
  s(abs(s)<=eps)=0; s=sqrt(s);
else                                                 % deviation vector
  if ~isa(s,'double')
    error('Deviation vector must be of class ''double''.')
  elseif ndims(s)>2
    error('Deviation vector cannot have more than two dimensions.')
  elseif length(s)>3;
    error('Deviation vector cannot have more than three elements.')
  elseif ~isreal(s)
    error('Deviation vector must be real.')
  elseif any(s<-eps)
    error('Deviation vector must be nonnegative.')
  end
  s(abs(s)<=eps)=0;
end
if all(s==0), R=0; return, end
s=s(s~=0); smax=max(s);
if smax>=1e7*min(s)
  error('Ratio of nonzero deviations must be less than 1e7.')
end
s=sort(s)/smax;
if nargin>1 && ~isempty(varargin{1}) && isa(varargin{1},'double') ...
  && numel(varargin{1})==1  
  p=real(varargin{1}); else p=0.95;                  % probability
end
if p<eps || p>1-eps, error('CRR requires 0 < P < 1.'), end
if nargin>2 && ~isempty(varargin{2}) && isa(varargin{2},'double') ...
  && numel(varargin{2})==1  && real(varargin{2})~=0
  tol=real(varargin{2}); else tol=1e-15;             % quadrature tolerance
end
if tol<0, tol=-tol; qflag=0; else qflag=1; end       % quadrature flag
options=optimset('display','off','tolx',eps);        % fzero tolerance
n=length(s);
if n==1                                              % crr.pdf (3)
  R=sqrt(2)*erfinv(p);
elseif n==2 && abs(diff(s))<eps                      % crr.pdf (10)
  R=sqrt(-2*log(1-p));
elseif n==3 && all(abs(diff(s))<eps)                 % crr.pdf (18)
  R=sqrt(-2*log(1-p)); R0=1e20;
  b=sqrt(pi)/2; count=0;
  while abs(R-R0)>tol && count<50                    % Newton's method
    a=R/sqrt(2); R0=R; count=count+1;
    R=R+(1+b*exp(a^2)*(p-erf(a))/a)/R;
  end
  if count==50, R=(R+R0)/2; end
else
  R0=sqrt(2)*erfinv(p);
  try                                                % crr.pdf (7,9,15,16)
    R=fzero(@difference,[R0 R0+1],options,s,p,tol,qflag);
  catch                                              % crr.pdf (7,15)
    s=flipud(s(:));
    R=abs(fzero(@difference,R0,options,s,p,tol,0));
  end
end
if nargin>3 && ~isempty(varargin{3}) && isa(varargin{3},'double')
  bootstrap(R,s,varargin{3})
end
R=R*smax;
%--------------------------------------------------------------------------
function f=difference(R,s,p,tol,qflag)
n=length(s);
if qflag && n==3 && any(abs(diff(s))<eps)             % crr.pdf (16)
  if abs(s(3)-s(2))<eps, s=[s(2) s(3) s(1)]; end
  a=R/sqrt(2); b=a/s(3);
  c=sqrt(1-(s(3)/s(1))^2);
  F=erf(b)-exp(-(a/s(1))^2)*erfir(c*b)/c;
elseif qflag && n==2 && max(s)/min(s)<=10 && p<=0.995 % crr.pdf (9)
  k=160; K=2*(0:k)+1;
  a=0.25*(s(1)^-2-s(2)^-2);
  b=0.25*(s(1)^-2+s(2)^-2);
  % central binomial coefficients (cbc)
  % C(2k,k) = (2k)!/(k!)^2 (Sloane A000984)
  cbc=ones(1,k+1);
  for n=1:k, cbc(n+1)=cbc(n)*(4-2/n); end
  F=sum(cbc.*(0.5*a/b).^K.*gammainc(b*R^2,K))/(a*s(1)*s(2));
else                                                  % crr.pdf (7,15)
  F=quadl(@integrand,0,R,tol,0,R,s)/(s(1)*s(2));
end
f=p-F;
%--------------------------------------------------------------------------
function f=integrand(r,R,s)
n=length(s);
a=abs(0.25*(s(1)^-2-s(2)^-2));
b=-0.25*(s(1)^-2+s(2)^-2);
r2=r.^2;
f=r.*exp((b+a)*r2);                                   % crr.pdf (7,15)
if a>eps
  % Io( z,1) = exp[-|Re(z)|] Io(z)
  % Io(-z,1) = Io(z,1)
  f=f.*besseli(0,a*r2,1);                             % crr.pdf (7,15)
end
if n==3
  f=f.*erf(sqrt(0.5*(R^2-r2))/s(3));                  % crr.pdf (15)
end
%--------------------------------------------------------------------------
function f=erfir(z)
% error function of pure imaginary or real argument;
% cf. erf(z) = gammainc(z^2,0.5)
y=imag(z);
if abs(y)<eps
  f=erf(real(z));
elseif abs(y)<=8                                      % A&S 7.1.29
  n=1:32;
  f=i/pi*(sum(exp(-0.25*n.*n)./n.*sinh(n*y))*2+y);
elseif abs(y)<27                                      % A&S 7.1.23
  m=193; s=1; y2=2*y^2;
  for n=m:-2:1, s=1+n*s/y2; end
  f=i*s*exp(y^2)/(y*sqrt(pi));
else
  f=sign(y)*i*inf;
end
%--------------------------------------------------------------------------
function bootstrap(R,s,mn)
disp('Running test of CRR.')
mn=fix(abs(real(mn)));
m=max(mn(1),1);
if length(mn)>1, n=max(mn(2),1); else n=1e6; end
j=length(s);
p=zeros(m,1);
for i=m:-1:1
  fprintf('%g ',i);
  X=randn(n,j)*diag(s);
  p(i)=sum(sum(X.^2,2)<R^2)/n;
end
disp(' ');
disp([num2str(mean(p)*100),'% of ',num2str(n*m),...
  ' random points were located within the computed radius.'])
