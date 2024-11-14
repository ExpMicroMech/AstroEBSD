
% Copyright 2021 EDAX INC
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, this list
% of conditions and the following disclaimer. Redistributions in binary form must
% reproduce the above copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials provided with the
% distribution. Neither the name of the copyright holder nor the names of its
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
%

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% 	LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.[8]

function [pats] = readEdaxUp(filename)
%readEdaxUp read patterns from a up1/up2 file

% determine bitdepth
if ~isfile(filename)
    throw(MException('readEdaxUp:invalid_argument', 'file does not exist'));
end
pixbytes = 0;
if '1' == filename(end)
    pixbytes = 1;
elseif '2' == filename(end)
    pixbytes = 2;
else
    throw(MException('readEdaxUp:invalid_argument', 'file must be *.up1 or *.up2'));
end

% open file and read header
fbytes = dir(filename).bytes;         
f = fopen(filename, 'r');
header = fread(f,4,'int32'); % version, patten width, pattern height, offset to patterns
vers = header(1);
w = header(2); % pattern width
h = header(3); % pattern width
offset = header(4); % offset to raw pattern data in bytes

% determine number of patterns and read
patbytes = w * h * pixbytes;
numpats = floor( (fbytes - offset) / patbytes );
fseek(f, offset, -1);
if 1 == pixbytes
    pats = fread(f, [w*h numpats], 'uint8') / 255;
else if 2 == pixbytes
    pats = fread(f, [w*h numpats], 'uint16') / 65535;
end

% reshape and rescale to [0,1]
pats = reshape(pats, w, h, numpats);
fclose(f);
end

