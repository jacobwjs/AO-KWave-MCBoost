% DESCRIPTION:
%       subscript to set literals for saving HDF5 files
%
% ABOUT:
%       author      - Bradley Treeby and Jiri Jaros
%       date        - 9th August 2012
%       last update - 9th August 2012
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2012 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% data type
DATA_TYPE_ATT_NAME          = 'data_type';
MATRIX_DATA_TYPE_MATLAB     = 'single';
MATRIX_DATA_TYPE_C          = 'float';
INTEGER_DATA_TYPE_MATLAB    = 'uint64';
INTEGER_DATA_TYPE_C         = 'long';

% real / complex
DOMAIN_TYPE_ATT_NAME        = 'domain_type';
DOMAIN_TYPE_REAL            = 'real';
DOMAIN_TYPE_COMPLEX         = 'complex';

% file descriptors
FILE_MAJOR_VER_ATT_NAME     = 'major_version';
FILE_MINOR_VER_ATT_NAME     = 'minor_version';
FILE_DESCR_ATT_NAME         = 'file_description';
FILE_CREATION_DATE_ATT_NAME = 'creation_date';
CREATED_BY_ATT_NAME         = 'created_by';

% file type
FILE_TYPE_ATT_NAME          = 'file_type';
HDF_INPUT_FILE              = 'input';
HDF_OUTPUT_FILE             = 'output';
HDF_CHECKPOINT_FILE         = 'checkpoint';

% file version information
HDF_FILE_MAJOR_VERSION      = '1';
HDF_FILE_MINOR_VERSION      = '0';

% compression level
COMPRESSION_LEVEL           = 9;