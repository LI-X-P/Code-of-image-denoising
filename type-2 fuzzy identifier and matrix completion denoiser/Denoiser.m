function [Y] = Denoiser(M, M_Omega, maxiteration, array_Omega)
%%
% INPUT: M is the original image
%        M_Omega is the incomplet image
%        maxiteration is the maximum number of iteration
%        array_Omega is the subset
%%
p = 2;
[m ,n] = size(M_Omega);
X_j = zeros(m, n);
u_j = randn(m, 1);
v_j = randn(n, 1);
MSE = [];
U = [];
V = [];
w = 1;
in = 0;

for iter = 1 : maxiteration
    R_j = M_Omega - X_j.*array_Omega; 
    %% calculate u_j
    in2 = 0;
    while (in2 < 15)
        c_j = v_j;
        u_j0 = u_j;
        a = norm(R_j - u_j * v_j.','fro')^2/norm(R_j,'fro')^2;
        for i = 1 : m
            g_l = R_j(i,:).';
            index = array_Omega(i,:).';

            u_j(i,1) = IRLS(c_j, g_l,p,10,index);
        end
        u_j = w*u_j + (1-w)*u_j0;
        
        %% calculate v_j
        %     d_j = s_j * u_j;
        d_j = u_j;
        v_j0 = v_j;
        for i = 1 : n
            h_q = R_j(:,i);
            index = array_Omega(:,i);
            v_j(i,1) = IRLS(d_j, h_q, p,50,index);
        end
        v_j = w*v_j + (1-w)*v_j0;
        b = norm(R_j - u_j * v_j.','fro')^2/norm(R_j,'fro')^2;
        if a-b < 0.0001
            in2 = in2 + 1;
        end
    end
    U = [U u_j];
    V = [V v_j];
    X_j = X_j + u_j * v_j.';
    %% judgement
    MSE = [MSE, norm(M-X_j,'fro')^2/norm(M,'fro')^2];
    Y = X_j;
    if iter~=1
        step_MSE = MSE(iter-1) - MSE(iter);
        if step_MSE < 0.0001
            in = in  + 1;
        end
        if in > 10
            break;
        end
    end

end
end