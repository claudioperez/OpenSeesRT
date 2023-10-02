function flaviaWriter(matl,dT,numGravStep)
% creates GiD post-processing results file from OpenSees data

%-------------------------OPEN FILE-----------------------------------

fid = fopen('renameMe.flavia.res','w');

fprintf(2,'* %s\n','Creating flavia.res from FEA. THIS MAY TAKE A FEW MINUTES ...')

fprintf(fid,'GiD Post Results File 1.0 \n\n');

%---------------DEFINE TIME VECTOR FOR GRAVITY ANALYSIS---------------
dg = dT/(numGravStep+1);
gtime = linspace(dg,dT-dg,numGravStep)';

%----------------------DISPLACEMENT-----------------------------------

dispFlavia(fid,gtime)
fprintf(2,'* %s\n','Done with displacements ...') 

%----------------------PORE PRESSURE----------------------------------

%porePressFlavia(fid,gtime)
%fprintf(2,'* %s\n','Done with porePressures ...')

%----------------------PORE PRESSURE RATIO----------------------------

ppRatioFlavia(fid,matl)
fprintf(2,'* %s\n','Done with porePressureRatio ...')

%------------------------STRESS---------------------------------------

%stressFlavia(fid,matl,gtime)
%fprintf(2,'* %s\n','Done with stress ...')

%------------------------CLOSE FILE-----------------------------------
fclose(fid);

return
