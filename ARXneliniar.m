load("datemotor.mat");
N = 200;
u = [zeros(50,1) ; idinput(N, 'prbs', [], [-0.7 0.7]) ; zeros(100,1) ; 0.4*ones(200,1)];
%[vel, alpha, t] = run(u, '2');
%plot(t, vel);
u_id = u(1:300);
u_val = u(300:550);
y_id = vel(1:300);
y_val = vel(300:550);

na = 7;
nb = 7;
phia = zeros(length(y_id),na);
phib = zeros(length(u_id),nb);

for i=1:length(y_id)
    for j=1:na
    if(i < j)
    phia(i,j) = 0;
    elseif(i == j)
    phia(i,j) = -y_id(1);
    else
    phia(i,j) = -y_id(i-j);
    end
    end
end

for i=1:length(u_id)
    for j=1:nb
    if(i < j)
    phib(i,j) = 0;
    elseif(i == j)
    phib(i,j) = u_id(1);
    else
    phib(i,j) = u_id(i-j);
    end
    end
end
phi = zeros(length(u_id),na+nb);

phi = [phia,phib];
teta = phi\y_id';

yaprox = zeros(1,length(y_val));
for i=1:length(y_val)
    for j =1:na
        if i>j
        yaprox(i) = yaprox(i) - teta(j)*y_val(i-j);
        elseif i==j
        yaprox(i) = yaprox(i) + teta(j)*y_val(1);
        end
    end
end

for i=1:length(u_val)
    for j =1:nb
        if i>j
        yaprox(i) = yaprox(i) + teta(na+j)*u_val(i-j);
        elseif i==j
        yaprox(i) = yaprox(i) + teta(na+j)*u_val(1);
        end
    end
end
e = y_val - yaprox;
MSE = 1/length(e) * sum(e.^2);
figure;
plot(y_val);
hold on
plot(yaprox,'Color','red');
legend('yval','yaprox');
title(['MSE', num2str(MSE)]); %eroare mare din cauza la date de pe motor

