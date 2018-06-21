clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum element size
Hmin = 1e-3;
% Maximum elements size
Hmax = 0.05;
% Step theta^{n+1} = theta^{n} - d * Dtheta^{n}
d = 1e-3;
% Maximum number of iterations
maxIter = 2000000;
% Dirichlet boundary condition on the outer boundary
ue = 1;
% Dirichlet boundary condition on the inner boundary
uc = 0;
% Smoothing of sensitivity field
smoothing = true;
% Reshaping of elements whose edges are lying on the controlled boundary
reshaping = true;
% Mesh motion solver method:
% - Laplace Equation ('laplace')
% - Spring Analogy ('spring')
% - Free Form Deformation ('FFD')
% - Radial Basis Functions ('RBF')
meshMotionSolver = 'RBF';
% Parameters for Laplace method
% Diffusivity (diff)
% - Inverse distance ('invdist')
invdist = true;
m = 2;
% Parameters for FFD method:
% Mi - number of control points along X axis
% Nj - number of control points along Y axis
% Xmin, Xmax - coordinates defining the control rectangle
Mi = 9;
Nj = 9;
Xmin = [-0.7, -0.6];
Xmax = -Xmin;
% Parameters for RBF method
% h - Dimensionless radial distance (r/h)
% RBFtype - Radial Basis Functions
% - Spline type, n odd ('Rn')
% - Thin Plate Spline, n even ('TPSn')
% - Multiquadratic ('MQ')
% - Inverse multiquadratic ('IMQ')
% - Inverse quadratic ('IQ')
% - Gaussian ('GS')
% n - Exponent
h = 2.5*0.3;
RBFtype = 'R';
n = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create geometry and mesh using MATLAB PDE TOOLBOX
[tr, Ie, Ic] = createMesh(Hmin, Hmax, uc, ue);
nodes = tr.Points;
elem = tr.ConnectivityList;

aux = ones(size(nodes, 1), 1)';
aux(Ic) = 0;
Inc = find(aux);

% Mesh edges
fb = freeBoundary(tr);

% Sort the inner circle nodes
Cc = sortNodes(fb, Ic);
% Sort the outer circle nodes
Ce = sortNodes(fb, Ie);

% Cost function difference between iterations
dJ = -1;
J = zeros(maxIter, 1);
J(1) = 1e308;
count = 0;
Ncount = 1;

V = zeros(maxIter, 1);
fsh_min = zeros(maxIter, 1);
fsh_mean = zeros(maxIter, 1);
fsz_min = zeros(maxIter, 1);
fsz_mean = zeros(maxIter, 1);
fss = 1;

iter = 1;

% Assemble FE matrices
[M, K, Cx, Cy, Ae, An, fshape] = assembleMatrices(tr);
Ae0 = Ae;

video = VideoWriter('shapeOptimization');
%video.FrameRate = 25;
video.FrameRate = 50;
video.Quality=100;
open(video);

fighandle = figure(20);
triplot(tr)
axis equal

filename = 'anim.gif';
% Capture the plot as an image
set(fighandle, 'color', 'white');
frame = getframe(fighandle); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 

while (iter <= maxIter && count < Ncount && fss > eps)

    % Assemble FE matrices
    [M, K, Cx, Cy, Ae, An, fshape] = assembleMatrices(tr);
    
    % Mesh quality metrics
    fshape(fshape < 0) = 0;
    fsh_min(iter) = min(fshape);
    fsh_mean(iter) = Ae'*fshape/sum(Ae);
    tau = Ae./Ae0;
    fsize = min(tau, 1./tau);
    fsize(fsize < 0) = 0;
    fsz_min(iter) = min(fsize);
    fsz_mean(iter) = Ae'*fsize/sum(Ae);
    fss = min(sqrt(fsize).*fshape);
    %fprintf('Min area: %d - Min fshape: %d - Max fshape: %d - Mean fshape: %d\n', min(Ae), min(fse), max(fse), Ae'*fse/sum(Ae))
    % Assemble source term
    F = assembleVector(@f, tr);

    Ko = K;
    Fo = F;
    
    % Apply Dirichlet boundary conditions
    [K, F] = dirichletBC(K, F, Ic, Ie, uc, ue);

    % Primal Solution
    u = K\F;
    % Gradients
    ux = M\(Cx*u);
    uy = M\(Cy*u);

    % Assemble vector of target values
    Ud = assembleVector(@target, tr);
    ud = M\Ud;
    
%     figure(300)
%     trisurf(tr.ConnectivityList, tr.Points(:, 1), tr.Points(:, 2), u)
    
    % Compute cost functional
    J(iter) = 0.5*(ud - u)'*M*(ud - u);
    if iter == 1
        dJ = -1;
    else
        dJ = J(iter) - J(iter-1);
    end
    
    if dJ >= eps || min(Ae) <= eps
        count = count + 1;
    else
        count = 0;
    end
    
    % Force term in the adjoint problem
    G = Ud - M*u;
    Go = G;
    % Apply zero boundary conditions on all boundaries
    G(Ic) = 0;
    G(Ie) = 0;

    % Dual Solution
    p = K\G;
    % Gradients
    px = M\(Cx*p);
    py = M\(Cy*p);
    
    % Compute domain volume
    V(iter) = sum(sum(full(M)));
    fprintf("Iteration %i - Volume %f - Cost %f - fss %f\n", iter, V(iter), J(iter), fss);
    fprintf("             - fsize_min %f - fshape_min %f\n", fsz_min(iter), fsh_min(iter));
    
    % Normal vectors on boundary edges and boundary nodes
    [nbe, xbe] = normalBoundaryEdge(tr, fb);
    [nbv, xbv] = normalBoundaryVertex(tr, fb, nbe);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute sensitivity field on the boundary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Nc = size(Ic, 1);
    Ne = size(Ie, 1);
    un = zeros(Nc+Ne, 1);
    pn = zeros(Nc+Ne, 1);
    sens = zeros(Nc+Ne, 1);
    I = [Ic; Ie];
    for i = 1:(Nc + Ne)
        k = I(i);
        nc = nbv(k, :);
        un(k) = nc*[ux(k); uy(k)];
        pn(k) = nc*[px(k); py(k)];
        sens(k) = 0.5*(u(k) - ud(k))^2 - un(k)*pn(k);
    end
    
    if smoothing
        %weight = [3/8 2/8 3/16];
        weight = [7/16 5/32 3/32 1/32];
        %weight = [1/2 6/32 1/32 1/32];
        sens = sensitivitySmoothing(sens, weight, Cc);
    end
    
    r = 0;
    lc = 0;
    for i = 1:Nc
        tc = nodes(Cc(i+1), :) - nodes(Cc(i), :);
        lc = lc + norm(tc); 
        r = r + 0.5*(sens(Cc(i+1)) + sens(Cc(i)))*norm(tc);
    end
    
    sens(Ic) = sens(Ic) - r/lc;
    Sc = max(sens(Ic));
    sens(Ic) = sens(Ic)/Sc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update coordinates of moving boundary nodes
    nci = nodes(Ic, :);
    nodes2 = nodes;
    nodes2(Ic, :) = nodes2(Ic, :) - d*sens(Ic).*nbv(Ic, :);
    %nodes2(Ic, :) = nodes2(Ic, :) - d*[1, 1];
    
    % Tangent reshaping
    if reshaping
        [Cc2, wtan] = tangential(nodes, Cc);
        nodes2(Cc2, :) = nodes2(Cc2, :) + wtan;
    end
    
    ncf = nodes2(Ic, :);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(20)
    axis equal
    
    switch meshMotionSolver
        case 'laplace'   
            v = zeros(size(nodes));
            v(Ic, :) = nodes2(Ic, :) - nodes(Ic, :);
            if invdist
                dist = distance(nodes, Ic);
                distmin = min(dist(Inc));
                dist(Ic) = distmin;
                
                K = assembleKdist(tr, 1./(dist.^m));      
                %K = assembleKdist(tr, An);
                K([Ic; Ie], :) = 0;
                K([Ic; Ie], [Ic; Ie]) = eye(size([Ic; Ie], 1));
            end
            w = K\v;           
            Nout = nodes + w;
        case 'spring'
            Id = [Ic; Ie];
            wd = nodes2(Id, :) - nodes(Id, :);
            [Is, K] = springAnalogy(tr, Id);
            Kss = K(Is, Is);
            Ksd = K(Is, Id);
            ws = -Kss\(Ksd*wd);
            w = zeros(size(nodes));
            w(Id, :) = wd;
            w(Is, :) = ws;
            Nout = nodes + w;
        case 'FFD'
            [nodesIr, Ir] = inRectangle(nodes, Xmin, Xmax);
            [NoutIr, P] = FFD(nodesIr, nci, ncf, Xmin, Xmax, Mi, Nj, 'extraRows');
            Nout = nodes;
            Nout(Ir, :) = NoutIr;
        case 'RBF'
            Nout = RBF(nodes, nodes2, [Ic; Ie], h, RBFtype, n);
        case 'none'
            Nout = nodes;
    end
    
    % Update mesh data
    tr = triangulation(elem, Nout);
    nodes = tr.Points;
    
    % Plot mesh
    figure(20)
    triplot(tr)
    hold on
    axis equal
    if strcmp(meshMotionSolver, 'FFD')
        plot(P(:, 1), P(:, 2), 'ro')
    end
    th = linspace(0, 2*pi, 100);
    plot(0.3*cos(th), 0.3*sin(th), 'black')
    hold off
    
    % Write to video/gif
    if rem(iter-1, 1) == 0
        writeVideo(video, getframe(gcf));

        % Capture the plot as an image 
        set(fighandle, 'color', 'white');
        frame = getframe(fighandle); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    iter = iter + 1;
    
end

writeVideo(video, getframe(gcf));
close(video);

iter = iter - 1;

J = J(1:iter);
figure(71)
semilogy(1:length(J), J, 'r')

V = V(1:iter);
figure(72)
plot(1:length(V), V, 'b')

fsh_min = fsh_min(1:iter);
fsh_mean = fsh_mean(1:iter);
fsz_min = fsz_min(1:iter);
fsz_mean = fsz_mean(1:iter);

dlmwrite('results.csv', [(1:iter)', V, J, fsz_min, fsz_mean, fsh_min, fsh_mean], 'precision', '%.16f')

% figure(51)
% trisurf(elem, nodes(:, 1), nodes(:, 2), u)
% figure(52)
% trisurf(elem, nodes(:, 1), nodes(:, 2), ux)
% figure(53)
% trisurf(elem, nodes(:, 1), nodes(:, 2), uy)
% figure(54)
% trisurf(elem, nodes(:, 1), nodes(:, 2), p)
% figure(55)
% trisurf(elem, nodes(:, 1), nodes(:, 2), px)
% figure(56)
% trisurf(elem, nodes(:, 1), nodes(:, 2), py)

