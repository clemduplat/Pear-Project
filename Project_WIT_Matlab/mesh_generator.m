clear all
model = createpde;
g = geometryFromEdges(model,@circlefunction);
fmesh = generateMesh(model, Hmax = 0.00005,GeometricOrder = 'linear');
%disp(fmesh);
%pdemesh(model);
elements = fmesh.Elements;
nodes = fmesh.Nodes;

%indLeft = find(nodes(2,:) < 0.001);
%indCircle = find((nodes(1,:).^2+nodes(2,:).^2) > 0.99);
%indEqual = union(indCircle,indLeft);
[p,e,t] = meshToPet(model.Mesh);
%disp(e)
%figure()
%plot(nodes(1,indEqual),nodes(2,indEqual),'ro')
%hold on
%for i = 1:max(indEqual)
    %text(nodes(1,i),nodes(2,i), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    %hold on
%end

% Determine the number of nodes
numNodes = size(nodes, 2);

% Initialize a matrix to store x, y, and boundary information
nodeMatrix = zeros(numNodes, 3);

% Copy x and y coordinates to the new matrix
nodeMatrix(:, 1:2) = nodes';

boundarynodes = [];

for k = e(1,:)
    boundarynodes = [boundarynodes, k];
end

% Mark nodes on the boundary with a 1
nodeMatrix(boundarynodes, 3) = 1;

[A,AE] = area(fmesh);

numElements = size(elements, 2);

elementsMatrix = zeros(numElements, 4);

elementsMatrix(:,1:3) = elements';

elementsMatrix(:,4) = AE;

boundary = zeros(size(e,2),3);
boundary(:,1:2) = e(1:2,:)';
lengths = zeros(size(e,2),1);
for i = 1: size(e,2)
    lengths(i) = sqrt((nodeMatrix(boundary(i,1),1) - nodeMatrix(boundary(i,2),1))^2 + (nodeMatrix(boundary(i,1),2)-nodeMatrix(boundary(i,2),2))^2);
end
boundary(:, 3) = lengths;

writematrix(elementsMatrix,'C:\Users\tdieg\team_01\semi_circle\readElements.txt')
writematrix(nodeMatrix,'C:\Users\tdieg\team_01\semi_circle\readNodes.txt')
writematrix(boundary,'C:\Users\tdieg\team_01\semi_circle\readBoundaries.txt')
%% Implementation

nbTriangles = size(elementsMatrix,1);
nbNodes = size(nodeMatrix,1);
nbBoundaryNodes = size(boundary,1);
% Parameters depend on storage conditions:
% storage = 1 if Orchard
% storage = 2 if Shelf life
% storage = 3 if Refrigerator
% storage = 4 if Precooling
% storage = 5 if Disorder inducing
% storage = 6 if Optimal CA
 storage = 6;

 % Diffusitivities
sig_ur = 2.8e-10;
sig_uz = 1.1e-9;
sig_vr = 2.32e-9;
sig_vz = 6.97e-9;

% Respiration kinetic parameters
T_ref = 293.15;
R_g = 8.314;
% Oxygen consumption
V_mu_ref = 2.39e-4;
Ea_vmu_ref = 80200;
% Fermentative carbon dioxide production
V_mfv_ref = 1.61e-4;
Ea_vmfv_ref = 56700;

% Michaelis-Menten constants
K_mu = 0.4103;
K_mv = 27.2438;
K_mfu = 0.1149;

% Respiration constant
r_q = 0.97;

% Convective mass transfer coefficients
rho_u = 7e-7;
rho_v = 7.5e-7;

% Ambient conditions
p_atm = 101300;

list_params = setParameters(storage, p_atm, R_g, Ea_vmfv_ref, Ea_vmu_ref, T_ref, V_mu_ref, V_mfv_ref);
T_cel = list_params(1);
eta_u = list_params(2);
eta_v = list_params(3);
T = list_params(4);
Cu_amb = list_params(5);
Cv_amb = list_params(6);
V_mu = list_params(7);
V_mfv = list_params(8);

Ku = zeros(nbNodes, nbNodes);
Kv = zeros(nbNodes, nbNodes);
Fu = zeros(nbNodes, nbNodes);
Fv1 = zeros(nbNodes, nbNodes);
Fv2 = zeros(nbNodes,1);

for k = 1: nbTriangles
    triangle = elementsMatrix(k, 1:3);
    triangleAir = elementsMatrix(k, 4)./1;
    coord = nodeMatrix(triangle, 1:2);
    r1 = coord(1,1)/1;
    r2 = coord(2,1)/1;
    r3 = coord(3,1)/1;
    z1 = coord(1,2)/1;
    z2 = coord(2,2)/1;
    z3 = coord(3,2)/1;
    Ku(triangle, triangle) = Ku(triangle, triangle) + KuAddition(triangleAir, sig_ur, sig_uz, r1, r2, r3, z1, z2, z3 );
    Kv(triangle, triangle) = Kv(triangle, triangle) + KvAddition(triangleAir, sig_vr, sig_vz, r1, r2, r3, z1, z2, z3 );
    Fu(triangle, triangle) = Fu(triangle, triangle) - partFuAddition(triangleAir, r1, r2, r3);
    Fv1(triangle, triangle) = Fv1(triangle, triangle) + partFuAddition(triangleAir, r1, r2, r3);
    Fv2(triangle) = Fv2(triangle) + partFvAddition(triangleAir, r1, r2, r3);
end

Hu1 = zeros(nbNodes, nbNodes);
Hu2 = zeros(nbNodes, 1);
Hv1 = zeros(nbNodes, nbNodes);
Hv2 = zeros(nbNodes, 1);

for i = 1:size(boundary,1)
    bound = boundary(i,:);
    pt1 = bound(1);
    pt2 = bound(2);

    length = bound(3)/1;
    r1 = nodeMatrix(pt1,1)/1;
    r2 = nodeMatrix(pt2,1)/1;

    Hu1(bound(1:2),bound(1:2)) = Hu1(bound(1:2),bound(1:2)) + rho_u*H1(length, r1, r2);
    Hu2(bound(1:2)) = Hu2(bound(1:2)) + rho_u*Cu_amb*H2(length, r1, r2);
    Hv1(bound(1:2),bound(1:2)) = Hv1(bound(1:2),bound(1:2)) + rho_v*H1(length, r1, r2);
    Hv2(bound(1:2)) = Hv2(bound(1:2)) + rho_v*Cv_amb*H2(length, r1, r2);
end

Cu0 = (Ku - Fu*V_mu/K_mu + Hu1)\Hu2;
Cv0 = (Kv + Hv1)\(V_mfv*Fv2 + Hv2 + r_q*Fv1*Cu0);
%disp(Cu0)
%disp(Cv0)
tolerance = 1e-12;
max_iterations = 10;

[Cu, Cv] = newton_raphson(Cu0, Cv0, tolerance, max_iterations, Ku, Kv, Fu, V_mu, K_mu, K_mv, Hu1, Hu2, Fv1, K_mfu, Hv1, Hv2, V_mfv, r_q);

make_contour_figure(nodeMatrix, Cu, Cv)
