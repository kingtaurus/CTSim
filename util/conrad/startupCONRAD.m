% initalized conrad interface for matlab
% first created by Martin Berger
% modified by Meng Wu, 2014-8-26

clear all;
close all;
%% Adding path and initialized conrad 


disp('Adding CONRAD_Source files to the Matlab path ...')

addpath(genpath('E:\konrad\CONRAD'),'-end'); % change to your CONRAD installing path

import java.*

c=CONRADmatlab('E:\konrad'); % change to your CONRAD installing path too

c.initialize;

return;

%% test if successfully installed

c.RecoPipeline

return;
