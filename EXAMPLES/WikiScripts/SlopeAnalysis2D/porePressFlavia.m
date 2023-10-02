function porePressFlavia(fid,gtime)

% node pointer
nodePtr = load('nodeInfo.dat');

% pore pressure data files
gpwp = load('GporePressure.out');
ppwp = load('porePressure.out');

% define time vector
gpwp(:,1) = [];
ptime = ppwp(:,1);
ppwp(:,1) = [];
time = [gtime;ptime];

% combine into single array
pwp = [gpwp;ppwp];
clear gpwp ppwp

% transformation to GiD format
[nStep,nNode] = size(pwp);

for k = 1:nStep
    fprintf(fid,'Result "a.  Nodal PorePressures" "Loading_Analysis"\t%12.5g Scalar OnNodes\n', time(k));
		fprintf(fid,'ComponentNames "Pore Pressure"\n');
		fprintf(fid,'Values\n');

    for j = 1:nNode
        fprintf(fid, '%d \t %-12.8e\n', nodePtr(j), pwp(k,j));
    end
    fprintf(fid,'End Values \n');
	fprintf(fid,' \n');
end

clear pwp

return
