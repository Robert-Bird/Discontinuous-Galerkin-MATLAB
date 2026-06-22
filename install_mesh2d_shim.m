function install_mesh2d_shim(meshDir)
% install_mesh2d_shim: deploy a size-compatibility shim into a Mesh2D folder
% so refine2 works on recent MATLAB, without modifying any Mesh2D file.
%
% Mesh2D (by D. Engwirda) is third-party, under its own licence, and is NOT
% shipped with this project. When you add a Mesh2D folder into the main
% directory, this writes <meshDir>/private/size.m - a FOLDER-SCOPED shim that
% makes one-argument size(X) return the first dimension for Mesh2D code only.
% That lets refine2's legacy "1:size(X)" work (recent MATLAB rejects it with
% "Colon operands must be real scalars"). It affects Mesh2D only and changes
% nothing elsewhere; no Mesh2D source file is touched.
%
% Called automatically by path_add when a Mesh2d folder is present. Safe to
% run repeatedly - it does nothing if the shim is already installed.
%
% Input(s):
% meshDir - the Mesh2D folder (default 'Mesh2d')

%  Copyright (C) 2026 Robert Bird

if nargin < 1 || isempty(meshDir), meshDir = 'Mesh2d'; end
if ~isfolder(meshDir), return, end                                         % no Mesh2D present -> nothing to do

privDir = fullfile(meshDir,'private');
target  = fullfile(privDir,'size.m');
if isfile(target), return, end                                             % already installed
if ~isfolder(privDir), mkdir(privDir); end

shim = {
'function varargout = size(varargin)'
'% Folder-scoped size shim for Mesh2D on recent MATLAB (installed by'
'% install_mesh2d_shim - part of the host project, NOT Mesh2D). One-argument'
'% size(X) returns the first dimension so refine2''s legacy 1:size(X) works;'
'% every other call passes through. Affects only code in this Mesh2d folder.'
'if nargin == 1 && nargout <= 1'
'    varargout{1} = builtin(''size'', varargin{1}, 1);'
'else'
'    [varargout{1:max(nargout,1)}] = builtin(''size'', varargin{:});'
'end'
'end'
};

fid = fopen(target,'w');
if fid < 0, warning('install_mesh2d_shim:write','Could not write %s',target); return, end
fprintf(fid,'%s\n',shim{:});
fclose(fid);
rehash                                                                     % make the new private function visible now
fprintf('install_mesh2d_shim: installed %s\n',target);
end
