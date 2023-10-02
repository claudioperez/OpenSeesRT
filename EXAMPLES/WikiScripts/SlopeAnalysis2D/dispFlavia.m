function dispFlavia(fid,gtime)

% node pointer
nodePtr = load('nodeInfo.dat');

% displacement data files
gdisp = load('Gdisplacement.out');
pdisp = load('displacement.out');

% define time vector
gdisp(:,1) = [];
ptime = pdisp(:,1);
pdisp(:,1) = [];
time = [gtime;ptime];

% combine into a single array
disp = [gdisp;pdisp];
clear gdisp pdisp

[nStep,nDisp] = size(disp);
nNode = nDisp/2;

for k = 1:nStep
    fprintf(fid,'Result "a.  Nodal Displacements" "Loading_Analysis"\t%12.5g Vector OnNodes\n', time(k));
		fprintf(fid,'ComponentNames "X-Displacement"  "Y-Displacement"\n');
		fprintf(fid,'Values\n');

    u = reshape(disp(k,:), 2, nNode);
    for j = 1:nNode
        fprintf(fid, '%d \t %-12.8e %-12.8e\n', nodePtr(j), u(:,j));
    end
    fprintf(fid,'End Values \n');
	fprintf(fid,' \n');
end

clear disp

return
