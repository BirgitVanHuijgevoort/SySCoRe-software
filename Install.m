%% Install

% Find folder of install file
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
disp('Please add either tensortoolbox or tensor lab to your path')
disp(' You can download the tensortoolbox for matlab from  https://www.tensortoolbox.org')

disp('MPT toolbox is not automatically installed')
try % Add mosek to path
    addpath ~/mosek/9.3/toolbox/r2015a
    addpath ~/tensor_toolbox-v3.2.1/
end
disp('Yalmip and Mosek are not automatically installed')
cd(folder)

