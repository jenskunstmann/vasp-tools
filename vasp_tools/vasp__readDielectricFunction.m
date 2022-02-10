function [real, imag] = vasp__readDielectricFunction(xmlfile)
% read the real and imaginary part of the dielectric function from a given
% vasprun.xml file 
% for xml files the parsing takes forever -> use the python tool:
% ~/prog/python/vasp-scripts/vgetdielectric.py 
% 
% imag(npnts, comp) : imaginary part of the dielectric function
% real(npnts, comp) : real part of the dielectric function, with 'npnts'
%                        energy points, and comp = 1..7
%                           1 = energy
%                           2 = xx
%                           3 = yy
%                           4 = zz
%                           5 = xy
%                           6 = yz
%                           7 = zx                 
% 
% xml help from:
% http://blogs.mathworks.com/community/2010/11/01/xml-and-matlab-navigating-a-tree/
try
   tree = xmlread(xmlfile);
catch
   error('Failed to read XML file %s.',xmlfile);
end

% get the xpath mechanism into the workspace
import javax.xml.xpath.*
factory = XPathFactory.newInstance;
xpath = factory.newXPath;

% compile and evaluate the XPath Expression to read the real and imaginary
% parts of the dielectric function
imag_expr = xpath.compile('modeling/calculation/dielectricfunction/imag/array/set');
imag_node = imag_expr.evaluate(tree, XPathConstants.NODE);
imag = str2num(imag_node.getTextContent);

real_expr = xpath.compile('modeling/calculation/dielectricfunction/real/array/set');
real_node = real_expr.evaluate(tree, XPathConstants.NODE);
real = str2num(real_node.getTextContent);
