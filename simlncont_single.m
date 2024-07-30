function epsilon = simlncont_single(mu,R,p)
% epsilon = simlncont(mu,R,p)
% this function generates a stochastic simulation of 
% a lognormal multifractal continuous cascade, using the method of
%  exponential of a moving average over a cone
%   with parameter mu, R the total scale ratio and p number of iterations
%   meaning that we generate p*R data points in 1d.
%
% simlncont.m
%
% F. G. Schmitt, July 2016 - March 2017
% 
% input:  - scalar mu   (0<mu<1)
%         - scalar p (integer, size limited by the memory of the computer)
%         - scalar R (scale ratio, R>10, with new version, it is easy to do this up to R=4096)
%
% ouput:  - epsilon, 1-column matrix containing 
%                           the multifractal simulation, length p*R
%
% 
% Modified by Yongxiang Huang, 25/01/2019


% first, generate a gaussian noise on a rectange with size
%     (R,(p+2)R)
% and with mean -mu/(2r(2r+1)) and standard deviation sqrt(mu/(r(2r+1)))



%correct finite scale ratio effect based on the numerical experiment
R0=[128,256,512,1024,2048];
alpha0=[0.9703,0.9953,1.0049,1.0066,1.0069];
if R<2048
    alpha0=interp1(R0,alpha0,R,'spline');
else
    alpha0=mean(alpha0(end-2:end));
end
tt=(p+2)*R;
ir=1:R;
m=-1./(2*ir.*(2*ir+1));
sdev=sqrt(1./(ir.*(2*ir+1)));
Nu=length(mu);
if Nu~=1
    error('single mu')
end

epsilon=zeros(R*p,Nu);
if exist('convnfft','file') %using FFT-based convolution, it is faster than conv provided by matlab
    Flag=1;
else
    Flag=0;
end

for i=1:R
    II = zeros(2*R+1,1);%to save memory, we do not generate the full matrix
    a1=R-i+1;
    a2=R+i+1;
    II(a1:a2)=1;
    L=sum(II,1);
    x=randn(tt,1);
    if Flag==0
        tmp=conv(x,II,'valid');
    elseif Flag==1
        tmp=convnfft(x,II,'valid'); %fast version 
    end
    
    tu=mu/(0.0262*mu^2-0.1038*mu+alpha0);%finite scaling ratio effect
    % tu=mu;
    epsilon(:)=epsilon(:)+tmp*sdev(i)*sqrt(tu)+L*m(i)*tu;

end

epsilon=exp(epsilon);
% epsilon=bsxfun(@times,epsilon,mean(epsilon,1).^-1);
% for i=1:Nu
% epsilon=epsilon/mean(epsilon);
% end

end
            
