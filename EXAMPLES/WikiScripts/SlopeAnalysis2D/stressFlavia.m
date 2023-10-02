function stressFlavia(fid,material,gtime)

% determine material type
if (strcmpi(material,'elastic')==1 || strcmpi(material,'J2')==1)
    nDOF = 3;
elseif (strcmpi(material,'PDMY')==1 || strcmpi(material,'PIMY')==1)
    nDOF = 5;
end

% load and combine data
for i = 1:4
    mLoad = sprintf('gstress{i} = load(''Gstress%i.out'');',i);
    eval(mLoad)
    gstress{i}(:,1) = [];
    mLoad = sprintf('pstress{i} = load(''stress%i.out'');',i);
    eval(mLoad)
    if i == 1
        ptime = pstress{i}(:,1);
    end
    pstress{i}(:,1) = [];
    stress{i} = [gstress{i};pstress{i}];
end
clear gstress pstress

time = [gtime;ptime];

[nStep,nStress] = size(stress{1});
nElem = nStress/nDOF;

if (nDOF==5)

for k = 1:nStep
    fprintf(fid,'GaussPoints "stress" ElemType Quadrilateral\n');
	fprintf(fid,'Number of Gauss Points: 4\n');
	fprintf(fid,'Natural Coordinate: Internal\n');
	fprintf(fid,'End Gausspoints\n\n');
	fprintf(fid,'Result "Gauss Point Stress" "Loading_Analysis"\t%12.5g', time(k));
    fprintf(fid,'\tPlainDeformationMatrix OnGaussPoints "stress"\n');
	fprintf(fid,'Values\n');

    for i = 1:4
        gp{i} = reshape(stress{i}(k,:), nDOF, nElem);
    end

    for j = 1:nElem
        fprintf(fid,'%6.0f  ', j);
        for i = 1:4
            fprintf(fid,'%12.6g %12.6g %12.6g %12.6g\n', gp{i}(1,j), gp{i}(2,j), gp{i}(4,j), gp{i}(3,j));
	    end
	end
	fprintf(fid,'End Values \n');
	fprintf(fid,'\n');
end

end

if (nDOF==3)

for k = 1:nStep
    fprintf(fid,'GaussPoints "stress" ElemType Quadrilateral\n');
	fprintf(fid,'Number of Gauss Points: 4\n');
	fprintf(fid,'Natural Coordinate: Internal\n');
	fprintf(fid,'End Gausspoints\n\n');
	fprintf(fid,'Result "Gauss Point Stress" "Loading_Analysis"\t%12.5g', time(k));
    fprintf(fid,'\tPlainDeformationMatrix OnGaussPoints "stress"\n');
	fprintf(fid,'Values\n');

    for i = 1:4
        gp{i} = reshape(stress{i}(k,:), nDOF, nElem);
    end

    for j = 1:nElem
        fprintf(fid,'%6.0f  ', j);
        for i = 1:4
            fprintf(fid,'%12.6g %12.6g %12.6g %12.6g\n', gp{i}(1,j), gp{i}(2,j), gp{i}(3,j), gp{i}(1,j));
	    end
	end
	fprintf(fid,'End Values \n');
	fprintf(fid,'\n');
end

end

clear stress gp

return
