clear all; close all;

X=[[1 2 4 4 12 2 4]' [2 3 4 2 4 1 5]'];
Y=[1 1 1 2 4 2 5]';

deg=2;
gam=3;

[co,A,b,res,monbas]=ridgepoly(X,Y,deg,gam);

[co monbas]

