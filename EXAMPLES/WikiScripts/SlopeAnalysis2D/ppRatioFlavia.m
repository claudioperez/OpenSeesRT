function ppRatioFlavia(fid,material)

%--------------------------PORE PRESSURE-------------------------------
% node pointer
nodePtr = load('nodeInfo.dat');

% pore pressure data files
pwp = load('porePressure.out');
time = pwp(:,1);
pwp(:,1) = [];
[nStep,nNode] = size(pwp);

% excess pore pressure
for k = 1:nStep
	exPwp(k,:) = abs(pwp(k,:) - pwp(1,:));
end

%-----------------------------STRESS-----------------------------------

% determine material type
if (strcmpi(material,'elastic')==1 || strcmpi(material,'J2')==1)
    nDOF = 3;
elseif (strcmpi(material,'PDMY')==1 || strcmpi(material,'PIMY')==1)
    nDOF = 5;
end

for i = 1:4
    mLoad = sprintf('stress{i} = load(''stress%i.out'');',i);
    eval(mLoad)
    stress{i}(:,1) = [];
end

[nStep,nStress] = size(stress{1});
nElem = nStress/nDOF;

% initial state of stress
for i = 1:4
    iniStress{i} = reshape(stress{i}(1,:), nDOF, nElem);
end

% initial vertical stress
for i = 1:4
    sigV(i,:) = iniStress{i}(2,:);
end

%-----------------------STRESS RECOVERY-------------------------------
P = [4,0,0;0,2.4,0;0,0,2.4];
Q = [1,1,1,1;...
     -sqrt(3/5),sqrt(3/5),sqrt(3/5),-sqrt(3/5);...
     -sqrt(3/5),-sqrt(3/5),sqrt(3/5),sqrt(3/5)];
S = P\Q;
E = [1, -1, -1; 1, 1, -1; 1, 1, 1; 1, -1, 1];
% transformation matrix
TR = E*S;

nStress = zeros(4,nElem);
for k = 1:nElem
    nStress(:,k) = TR*[sigV(1,k);sigV(2,k);sigV(3,k);sigV(4,k)];
end

% load element connectivity information
eleList = load('elementInfo.dat');

% array of nodal information
nodes = zeros(nNode,3);
nodes(:,1) = nodePtr(:,1);

% loop over elements to assign local nodal values to global list
for k = 1:nElem
    ida = nodes(:,1)==eleList(k,2);
    nodes(ida,2) = nodes(ida,2) + nStress(1,k);
    nodes(ida,3) = nodes(ida,3) + 1;
    idb = nodes(:,1)==eleList(k,3);
    nodes(idb,2) = nodes(idb,2) + nStress(2,k);
    nodes(idb,3) = nodes(idb,3) + 1;
	idc = nodes(:,1)==eleList(k,4);
    nodes(idc,2) = nodes(idc,2) + nStress(3,k);
    nodes(idc,3) = nodes(idc,3) + 1;
	idd = nodes(:,1)==eleList(k,5);
    nodes(idd,2) = nodes(idd,2) + nStress(3,k);
    nodes(idd,3) = nodes(idd,3) + 1;
end

% compute average values at nodes which have multiple elements
nodes(:,2) = nodes(:,2)./nodes(:,3);

% compute pore pressure ratio
fullStress = abs(nodes(:,2)*ones(1,nStep));
exPwp = exPwp';
ru = exPwp./fullStress;

% transformation to GiD format
for k = 1:nStep
    fprintf(fid,'Result "a.  PorePressureRatio" "Loading_Analysis"\t%12.5g Scalar OnNodes\n', time(k));
		fprintf(fid,'ComponentNames "Pore Pressure Ratio"\n');
		fprintf(fid,'Values\n');

    for j = 1:nNode
        fprintf(fid, '%d \t %-12.8e\n', nodePtr(j), ru(j,k));
    end
    fprintf(fid,'End Values \n');
	fprintf(fid,' \n');
end

return
