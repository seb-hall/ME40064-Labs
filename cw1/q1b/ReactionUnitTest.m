%% Test 1: test symmetry of the matrix
% % Test that this matrix is symmetric

tol = 1e-14; % test tolerance 
lambda = 2; % reaction coefficient
eID = 1; % element ID

xmin = 0;
xmax = 1;
Ne = 10;
msh = OneDimLinearMeshGen(xmin, xmax, Ne);

elemat = ReactionElemMatrix(lambda, eID, msh);

assert(abs(elemat(1,2) - elemat(2,1)) <= tol)

%% Test 2: test 2 different elements of the same size produce same matrix
% % Test that for two elements of an equispaced mesh, the element matrices 
% % calculated are the same.

tol = 1e-14; % test tolerance
lambda = 5; % reaction coefficient
eID = 1; % element ID

xmin = 0;
xmax = 1;
Ne = 10;
msh = OneDimLinearMeshGen(xmin, xmax, Ne);

elemat1 = ReactionElemMatrix(lambda, eID, msh);

eID = 2; %element ID
elemat2 = ReactionElemMatrix(lambda, eID, msh);

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 3: test that one matrix is evaluated correctly
% % Test that element 1 of the (equispaced) three element mesh problem is evaluated correctly

tol = 1e-14; % test tolerance 
lambda = 2.5; % reaction coefficient
eID = 1; % element ID

xmin = 0;
xmax = 1;
Ne = 3;
msh = OneDimLinearMeshGen(xmin, xmax, Ne);

elemat1 = ReactionElemMatrix(lambda, eID, msh);

matrix = [2, 1; 1, 2];
elemSize = (xmax - xmin) / Ne;
elemat2 = matrix * (lambda * elemSize / 6);

diff = elemat1 - elemat2; % calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); % calculate the total squared error between the matrices
assert(abs(diffnorm) <= tol)

%% Test 4: test that different sized elements in a mesh are evaluated correctly - element 1
% % Test that elements in a non-equally spaced mesh are evaluated correctly

tol = 1e-14; % test tolerance
lambda = 1; % reaction coefficient
eID = 1; % element ID

xmin = 0;
xmax = 1;
Ne = 5;
msh = OneDimSimpleRefinedMeshGen(xmin, xmax, Ne);

elemat1 = ReactionElemMatrix(lambda, eID, msh);

elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1); % get element size from mesh
elemat2 = [2, 1; 1, 2] * (lambda * elemSize / 6); % calculate resultant matrix

diff = elemat1 - elemat2; % calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); % calculate the total squared error between the matrices
assert(abs(diffnorm) <= tol)

%% Test 5: test that different sized elements in a mesh are evaluated correctly - element 4
% % Test that elements in a non-equally spaced mesh are evaluated correctly

tol = 1e-14; % test tolerance
lambda = 1; % reaction coefficient
eID = 4; % element ID

xmin = 0;
xmax = 1;
Ne = 5;
msh = OneDimSimpleRefinedMeshGen(xmin, xmax, Ne);

elemat1 = ReactionElemMatrix(lambda, eID, msh);

elemSize = msh.elem(eID).x(2) - msh.elem(eID).x(1); % get element size from mesh
elemat2 = [2, 1; 1, 2] * (lambda * elemSize / 6); % calculate resultant matrix

diff = elemat1 - elemat2; % calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); % calculate the total squared error between the matrices
assert(abs(diffnorm) <= tol)