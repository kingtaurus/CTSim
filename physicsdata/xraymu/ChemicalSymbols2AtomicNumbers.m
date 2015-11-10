function [zs,wgts] = ChemicalSymbols2AtomicNumbers(elems)
%  function [zs,wgts] = ChemicalSymbols2AtomicNumbers(elems)
%  returns atomic numbers and weights corresponding to symbols in cell array of strings
%  Error if symbol not found
%  inputs:
%    elems: cell array of strings (case sensitive)
%  outputs:
%    zs: an array of atomic number integers 
%    wgts:atomic weights
% for more information see http://aprendtech.com/wordpress/?p=45
%  REA 5/28/07

zs = zeros(length(elems),1);
wgts = zeros(length(elems),1);

[weights,symbols] = XrayData;
for kelem = 1:length(elems)
  z = strmatch( elems{kelem},symbols,'exact'); 
  if isempty(z)
    error('no match to symbol %s',char(elems{kelem}));
  end
  zs(kelem) = z;
  wgts(kelem) = weights(z);
end % for kelem = 1:length(elems)
