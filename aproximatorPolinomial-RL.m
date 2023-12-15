load('proj_fit_08');
M = 15;
MIN_MSE_val = inf;
MIN_MSE_id = inf;
%date identificare
x1_id = id.X{1};
x2_id = id.X{2};
y_id = id.Y;

%date validare
x1_val = val.X{1};
x2_val = val.X{2};
y_val = val.Y;

figure;
mesh(x1_id,x2_id,y_id);
title('Date identificare');
figure;
mesh(x1_val,x2_val,y_val);
title('Date validare');

%initializare MSE si y_final ca fiind goale inafara loop-ului de training pentru a nu
%pierde valorile
MSE_val = [];
MSE_id = [];
y_final_val = [];
y_final_id = [];
for z = 1:M

%initializare matrici necesare ca fiind goale
phi_id = [];
phi_val = [];
theta = [];

%phi identificare
for i = 1:length(x1_id)
    for j = 1:length(x2_id)
        polinom = genereazaPolinom(x1_id(i), x2_id(j), z);
        phi_id = [phi_id; polinom];
    end
end

%phi validare
for i = 1:length(x1_val)
    for j = 1:length(x2_val)
        polinom = genereazaPolinom(x1_val(i), x2_val(j), z);
        phi_val = [phi_val; polinom];
    end
end

%phi_3d = reshape(phi,[length(x1_id),length(x2_id),width(phi)]);
y_concatenat = reshape(y_id,[],1);

%generare iesirii aproximate
theta = phi_id\y_concatenat;
y_aprox_id = phi_id * theta;
y_aprox_val = phi_val * theta;
y_aprox2d_val = reshape(y_aprox_val,[length(y_val),width(y_val)]);
y_aprox2d_id = reshape(y_aprox_id,[length(y_id),width(y_id)]);

%calcul MSE
e_val = y_val-y_aprox2d_val;
e_val = reshape(e_val,[length(e_val)*length(e_val),1]);
MSE_val(z) = 1/length(e_val) * sum(e_val.^2);
e_id = y_id - y_aprox2d_id;
e_id = reshape(e_id,[length(e_id)*length(e_id),1]);
MSE_id(z) = 1/length(e_id) * sum(e_id.^2);

if MSE_val(z) < MIN_MSE_val
y_final_val = y_aprox2d_val;
MIN_MSE_val = MSE_val(z);
grad_minim_val = z;
end

if MSE_id(z) < MIN_MSE_id
y_final_id = y_aprox2d_id;
MIN_MSE_id = MSE_id(z);
grad_minim_id = z;
end
end

%afisare rezultate finale
figure;
mesh(x1_val,x2_val,y_final_val);
title(['yaprox val optim, MSE = ',num2str(MIN_MSE_val),',grad = ',num2str(grad_minim_val)]);


figure;
mesh(x1_id,x2_id,y_final_id);
title(['yaprox id optim, MSE = ',num2str(MIN_MSE_id),',grad = ',num2str(grad_minim_id)]);

figure;
plot(MSE_val);
hold on
plot(MSE_id,'Color','Red')
legend('MSE val','MSE id')
title('MSE');

function polinom = genereazaPolinom(x1, x2, m)
    polinom = [];
    
    for grad_x1 = 0:m
        for grad_x2 = 0:(m - grad_x1)
            termen = x1.^grad_x1.* x2.^grad_x2;
            polinom = [polinom, termen];
        end
    end
end
