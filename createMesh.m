function [tr, Ie, Ic] = createMesh(Hmin, Hmax, uc, ue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE GEOMETRY AND MESH USING MATLAB PDE TOOLBOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = createpde;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Geometry
% - Circle hole inside circle ('circle')
% - Square hole inside circle ('square')
% - Square hole inside square ('2square')
geom = 'square';
% Hole center coordinates
chx = 0.2;
chy = 0.2;
% Hole radius
rh = 0.3;
% Circle center coordinates
ccx = 0;
ccy = 0;
% Circle radius
rc = 1;

switch geom
    case 'circle'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Circle with circular hole
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        C1 = [1,chx,chy,rh]';       
        C2 = [1,ccx,ccy,rc]';
        %C1 = [C1;zeros(length(C2) - length(C1),1)];
        C2 = [C2;zeros(length(C1) - length(C2),1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'square'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Circle with rectangular hole
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        l = sqrt(pi)*rh;
        C1 = [3,4,chx-l/2,chx+l/2,chx+l/2,chx-l/2,chy-l/2,chy-l/2,chy+l/2,chy+l/2]';     
        C2 = [1,ccx,ccy,rc]';
        C2 = [C2;zeros(length(C1) - length(C2),1)];  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    case '2square'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Square with rectangular hole
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        l = 1/2;
        d = l/sqrt(2);
        ang = pi/16;
        chx = rc/2;
        chy = rc/2;     
        C1 = [3,4,chx+d*cos(pi/4+ang),chx+d*cos(3*pi/4+ang),chx+d*cos(5*pi/4+ang),...
            chx+d*cos(7*pi/4+ang),chy+d*sin(pi/4+ang),chy+d*sin(3*pi/4+ang),...
            chy+d*sin(5*pi/4+ang),chy+d*sin(7*pi/4+ang)]';
        C2 = [3,4,0,rc,rc,0,0,0,rc,rc]';
        C2 = [C2;zeros(length(C1) - length(C2),1)];
end

gm = [C2,C1];
sf = 'C2-C1';
ns = char('C2','C1');
ns = ns';
g = decsg(gm,sf,ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geometryFromEdges(model,g);

% Plot Geometry
figure(1) 
pdegplot(model,'EdgeLabels','on')
axis equal
xlim([0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

generateMesh(model, 'Hmin', Hmin, 'Hmax', Hmax, 'GeometricOrder', 'linear');

% Plot Mesh
figure(2)
pdeplot(model,'NodeLabels','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set boundary conditions in order to identifiy boundaries with different
% Dirichlet conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch geom
    case 'circle'
        applyBoundaryCondition(model,'dirichlet','Edge',1:4,'u',0);
        applyBoundaryCondition(model,'dirichlet','Edge',5:model.Geometry.NumEdges,'u',1);
    case 'square'
        applyBoundaryCondition(model,'dirichlet','Edge',1:4,'u',1);
        applyBoundaryCondition(model,'dirichlet','Edge',5:model.Geometry.NumEdges,'u',0);
    case '2square'
        applyBoundaryCondition(model,'dirichlet','Edge',[1,2,6,7],'u',0);
        applyBoundaryCondition(model,'dirichlet','Edge',[3,4,5,8],'u',1);
end

f = @(region,state) 32*pi^2*sin(4*pi*region.x).*sin(4*pi*region.y);
% Assemble matrices
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',f);
FEM = assembleFEMatrices(model);

% Indices for boundary nodes
% Ic - Inner circle
% Ie - Outer square
Ic = find(FEM.R);
Ie = find(FEM.R - 1);

% Apply boundary conditions for the actual problem
switch geom
    case 'circle'
        applyBoundaryCondition(model,'dirichlet','Edge',1:4,'u',ue);
        applyBoundaryCondition(model,'dirichlet','Edge',5:model.Geometry.NumEdges,'u',uc);
    case 'square'
        applyBoundaryCondition(model,'dirichlet','Edge',1:4,'u',uc);
        applyBoundaryCondition(model,'dirichlet','Edge',5:model.Geometry.NumEdges,'u',ue);
    case '2square'
        applyBoundaryCondition(model,'dirichlet','Edge',[1,2,6,7],'u',ue);
        applyBoundaryCondition(model,'dirichlet','Edge',[3,4,5,8],'u',uc);
end

% applyBoundaryCondition(model,'dirichlet','Edge',1:4,'u',ue);
% applyBoundaryCondition(model,'dirichlet','Edge',5:model.Geometry.NumEdges,'u',uc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve PDE
results = solvepde(model);
% figure(3)
% pdeplot(model,'XYData',results.NodalSolution,'ZData',results.NodalSolution)
% figure(4)
% pdeplot(model,'XYData',results.XGradients,'ZData',results.XGradients)
% figure(5)
% pdeplot(model,'XYData',results.YGradients,'ZData',results.YGradients)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create triangulation from the PDE TOOLBOX mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tr = triangulation(model.Mesh.Elements', model.Mesh.Nodes');

% Plot results
elem = model.Mesh.Elements';
nodes = model.Mesh.Nodes';
u = results.NodalSolution;
ux = results.XGradients;
uy = results.YGradients;

% figure(11)
% trisurf(elem, nodes(:, 1), nodes(:, 2), u)
% figure(12)
% trisurf(elem, nodes(:, 1), nodes(:, 2), ux)
% figure(13)
% trisurf(elem, nodes(:, 1), nodes(:, 2), uy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end