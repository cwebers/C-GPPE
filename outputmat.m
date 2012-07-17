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

dimline=size(matrix,1);
dimcol=size(matrix,2);
fmt = [repmat(' %d',1,dimcol)];

fid=fopen(filename,'wt');
   fprintf(fid,['%d %d \n'],size(matrix,1),size(matrix,2));
   fprintf(fid,[fmt '\n'],matrix.');
fclose(fid);



