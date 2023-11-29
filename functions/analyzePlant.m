function lengthBernoulli = analyzePlant(A,B,C,p_max,n_x,n_u,n_y,epsilon_b,folderData)
    lengthBernoulli = ceil(log(epsilon_b)/log(1-p_max));
    fprintf('Max fading int under which the system remains controllable with probability 1 ')
    fprintf('lasts %u time steps.\n\n',lengthBernoulli)
    cd(folderData)
    save('pendulum.mat', 'A','B','C','p_max','n_x','n_u','n_y','lengthBernoulli')
end