function [ data ] = load_netsim_data( filepath , params )
% *NETSIM*
%
%	LOAD NETSIM DATA
%
%	Loads data from a NETSIM binary file.
%
%	INPUT
%	filepath: /path/to/file
%
%	OUTPUT
%	data - data values (double)
%

assert( ~isempty( regexp( filepath, '\.bin?$', 'match', 'once' ) ) )

% open and read file
fid = fopen( filepath, 'rb' );
data = fread( fid, 'double' );
data = reshape( data, params.N, [] );

% close file
fclose( fid );

end