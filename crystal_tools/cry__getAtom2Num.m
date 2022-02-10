function [atom2num] = cry__getAtom2Num()
% get a generic cell array that related atomic symbols and atomic numbers
% to be extended to all atomic numbers ...
%
% TO DO: extend the data base with atomic radii and colors

atom2num = {...
 'H',  1
 'He', 2
 'Li', 3
 'Be', 4
 'B',  5
 'C',  6
 'N',  7
 'O',  8
 'F',  9
 'Ne', 10
 'Na', 11
 'Mg', 12
 'Al', 13 
 'Si', 14
 'P',  15
 'S',  16
 'Cl', 17
 'Ar', 18
 'K',  19
 'Ca', 20
 'Sc', 21
 'Ti', 22 
 'Mn', 25
 'Fe', 26
 'Ga', 31
 'As', 33
 'Se', 34
 'Mo', 42
 'Ag', 47
 'Sn', 50
 'W', 74
 'Re', 75
 'Pb', 82
 };