function [ si, st ] = load_netsim_spikes( timespath, idspath )
% *NETSTIM*
%
%	LOAD NETSIM SPIKES
%
%	Loads spikes from a NETSIM binary file.
%
%	INPUT
%	filepath: /path/to/file
%
%	OUTPUT
%	ids - spike ids (int32) (note: 1 was added to this vector to adhere to Matlab's 1-based indexing)
%	times - spike times (double)
%

assert( ~isempty( regexp( timespath, 'spk_times\.bin?$', 'match', 'once' ) ) )
assert( ~isempty( regexp( idspath, 'spk_ids\.bin?$', 'match', 'once' ) ) )

fid = fopen(timespath, 'rb');
st = fread(fid, 'double' ) + 1e-4;
fclose(fid);

fid = fopen(idspath, 'rb');
si = fread(fid, 'uint32' ) + 1;
fclose(fid);