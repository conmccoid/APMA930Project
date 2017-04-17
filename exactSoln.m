function results = exactSoln(x,t,e)
e1 = exp((-x + 0.5 - 4.95*t)/(20*e));
e2 = exp((-x + 0.5 - 0.75*t)/(4*e));
e3 = exp((-x + 0.375)/(2*e));

results = (0.1*e1 + 0.5*e2 + e3)./(e1 + e2 + e3);