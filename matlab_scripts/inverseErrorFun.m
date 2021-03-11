%%%%%%%%%%%%%%%
%%%
%%%   routine for computing inverse error function
%%%


x = -3:0.01:3;
f = exp(-x.^2)/sqrt(pi);
F = cumtrapz(x,f);
F0 = 1/2+erf(x)/2;


close(figure(1));
f1=figure(1); set(f1,'position',[550 440 890 360]);
%
subplot(1,2,1);
plot(x,f); title('distibution');
axis('tight');

subplot(1,2,2);
plot(x,F); title('cumulative distribution');
hold on; plot(x,F0,'linestyle','--');
axis('tight');



%%%%%%%%%%%%%%%
%%%
%%%  compute inverse error function and sample to reproduce normal dfn
%%%
%%%%%%%%%%%%%%
%F0 = 1/2+erf(x)/2;


%%%   randomly sample
%
Nsamples = 1.0e5;
v = zeros(1,Nsamples);
for i=1:length(v)
    v(i) = erfinv(2*rand-1);
end
%v = randn(Nsamples,1)/sqrt(2);

%%%   build distribution
%
dfn = zeros(size(x));
dx = x(2)-x(1);
for i=1:length(v)
    thisv = v(i);
    [~,index]=min(abs(x-thisv));
    dfn(index) = dfn(index)+1;
end
normC = trapz(x,dfn);
dfn = dfn/normC;    

figure(1);
subplot(1,2,1); hold on;
plot(x,dfn);


%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   Now that we know matlabs routine works, test ours...
%%%   (https://stackoverflow.com/questions/5971830/need-code-for-inverse-error-function)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%

y = -1:0.01:1;
xinv = zeros(size(y));


a = [0.886226899,  -1.645349621,  0.914624893, -0.140543331];
b = [-2.118377725,  1.442710462, -0.329097515,  0.012229801];
c = [-1.970840454, -1.624906493,  3.429567803,  1.641345311];
d = [ 3.543889200,  1.637067800];
y0 = 0.7;
for i=1:length(y)
   thisy = y(i);
   if (abs(thisy)==1.0 )
      xinv(i) = -thisy*log(0.0);
   elseif (thisy < -y0)
      z = sqrt(-log((1.0+thisy)/2.0));
      xinv(i) = -(((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0);
   else
      if (thisy<y0)
          z = thisy^2;
          %xinv(i) = thisy*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(4))*z+b(2))*z+b(1))*z+1.0);
          xinv(i) = thisy*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(3))*z+b(2))*z+b(1))*z+1.0);
      else
          z = sqrt(-log((1.0-thisy)/2.0));
          xinv(i) = (((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0);
      end
%     //polish x to full accuracy
      xinv(i) = xinv(i) - (erf(xinv(i)) - thisy) / (2.0/sqrt(pi) * exp(-xinv(i).^2));
      xinv(i) = xinv(i) - (erf(xinv(i)) - thisy) / (2.0/sqrt(pi) * exp(-xinv(i).^2));
   end
   
end


close(figure(2));
f2=figure(2);
plot(y,erfinv(y),'displayName','matlab');
hold on; plot(y,xinv,'r--','displayName','JRA');
title('inverse error function');
hold on; legend('show');

error = abs(erfinv(y)-xinv);

figure(3); hold on;
plot(y,error);
set(gca,'yscale','log');
title('error in calculation of inverse error function');



