function results = exactPeriodicSoln(x,t,v,m)
numerator = (x - 1).*exp(-(x - pi).^2/(4*v*(t+1)));
denominator = exp(-(x - pi).^2/(4*v*(t+1)));

for r = 1:m
    numerator = numerator + (x - pi*(2*r + 1) ).*exp( -(x - pi*(2*r + 1)).^2 / (4*v*(t+1)) ) + ...
        (x - pi*(1 - 2*r)).*exp( -(x - pi*(1 - 2*r)).^2 / (4*v*(t+1)) );
    denominator = denominator + exp( -(x - pi*(2*r+1)).^2 / (4*v*(t+1)) ) + ...
        exp( -(x - pi*(1 - 2*r)).^2 / (4*v*(t+1)) );
end

results = numerator./denominator/(t+1);

% numerator = zeros(size(x));
% denominator = trapz( exp( - (1 - cos(pi*x))/(2*pi*v) ) );
% for k = m:n
%     a = 2*trapz( exp(-(1 - cos(pi*x))/(2*pi*v)).*cos(k*pi*x) );
%     numerator = numerator + a*exp(-k^2*pi^2*v*t)*k*sin(k*pi*x);
%     denominator = denominator + a*exp(-k^2*pi^2*v*t)*cos(k*pi*x);
% end
% results = 2*pi*v*numerator./denominator;