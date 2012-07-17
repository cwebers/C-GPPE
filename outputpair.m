% Copyright (c) 2012, National ICT Australia
% All rights reserved.
%
% The contents of this file are subject to the Mozilla Public License
% Version 2.0 (the "License"); you may not use this file except in
% compliance with the License. You may obtain a copy of the License at
% http://www.mozilla.org/MPL/
%
% Software distributed under the License is distributed on an "AS IS"
% basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
% License for the specific language governing rights and limitations
% under the License.
% function outputmat(matrix,filename)
function outputpair(all_pairs,filename)
nbpairs=size(all_pairs,1);

fid=fopen([filename '/informations.txt'],'wt');
   fprintf(fid,['%d'],size(all_pairs,1));
fclose(fid);



for i=1 : size(all_pairs,1)
    file=['/pair' num2str(i-1) '.txt'];
   
    outputmat(all_pairs{i},[filename file])
end