function [mainax, insetax] = vasp__plotDOS(dos)
% make a general DOS plot as specified by the 'dos' structure;
% we can plot an arbitrary number of curves in both single (non
% spin-polarized) and split-panel (spin-polarized) mode
%
% dos.e(pnts)           : energy axis
% dos.ddos(pnts,curve)  : holds the data of the dos curves to be plottet 
% dos.color(curve,comp) : rgb color of each curve, comp=1(r),2(g),3(b)
% dos.lstyle{curve}     : line style for each curve
% dos.legend{curve}     : legend of each curve
% dos.sort{panel,curve} : defines which curves are plottet in which panel
%                         and in what order; the order is important for the
%                         order in the legend box; if panel=2 is empty,
%                         only one panel is plotted 
%                         panel=1 : upper panel, usually spin-up DOS
%                         panel=2 : lower panel, usually spin-down DOS 
% dos.plotinset = 1=yes, 0=no, plot inset? 
% dos.emin = plotting ranges of the energy axis
% dos.emax =  
% dos.dosmin =  only significant for spin-pol. calcs
% dos.dosmax = 
% dos.inset_emin = plotting ranges of the INSET
% dos.inset_emax = 
% dos.inset_dosmin =  only significant for spin-pol. calcs
% dos.inset_dosmax = 
% dos.insetshift = shifts the inset and main plot in y direction by insetshift*dosmax 


% single or double panel plot?
if (size(dos.sort,2) == 1)
    singlepanel = 1;
else
    singlepanel = 0;
end

% look if the main plot already exists, if not, create it
% this is necessary if we want to run the program several times to add
% cuves to the figure
% HINT: This feature has been disabled for vasp_PlotProjDOS_ManyAtoms()
%mainax = findobj(gcf, 'UserData', 'Main Axes');
%if (isempty(mainax))  % if it does not exist
mainax = gca;  
%    set(mainax, 'UserData', 'Main Axes');
%end
hold on     % hold main plot
box on      % draw boxed graph
%xlabel(mainax,'E-E_F (eV)'); 

% look if the inset plot already exists, if not, create it
if(dos.plotinset)    
    % This axis has its origin at RELATIVE location (.15 , .57)
    % The X- and Y-axes lengths are both 0.34 (i.e. 34% of main figure)      
    %insetax = findobj(gcf, 'UserData', 'Inset Axes');
    %if (isempty(insetax))
    insetax = axes('pos',[.15 .57 .34 .34]);
    %    set(insetax, 'UserData', 'Inset Axes');
    %end     
    hold on     % hold inset 
    box on      % draw boxed graph
else
    insetax = 0;
end;  
       
if(singlepanel) %  non-spin-polarized DOS    
    
    % plot several lines; we have to use a loop because otherwise
    % we cannot properly define the color of the individual lines
    for i=dos.sort{1}      
        plot(mainax, dos.e, dos.ddos(:,i), dos.lstyle{i},'color', dos.color(i,:)); 
    end
    plot(mainax, [0 0], [0 dos.insetshift*dos.dosmax], '--k');      % plot vertical line for Fermi level        
    axis(mainax, [dos.emin dos.emax 0 dos.dosmax]);         % range     
    
    % plot inset
    if(dos.plotinset)
        for i=dos.sort{1} 
            plot(insetax, dos.e, dos.ddos(:,i), dos.lstyle{i},'color', dos.color(i,:)); 
        end
        plot(insetax,[0 0], [0 dos.inset_dosmax], '--k');      % plot vertical line for Fermi level
        axis(insetax,[dos.inset_emin dos.inset_emax 0 dos.inset_dosmax]);         % range   
        % allow for some space for the inset on the main axes
        axis(mainax, [dos.emin dos.emax 0 dos.insetshift*dos.dosmax]);
        set(insetax,'xtick',[dos.inset_emin:1:dos.inset_emax]);                 
    end;     
                         
else % double panel, spin-polarized DOS    
    
    for i=dos.sort{1}      % spin-up panel
        plot(mainax, dos.e, dos.ddos(:,i), dos.lstyle{i},'color', dos.color(i,:)); 
    end    
    for i=dos.sort{2}      % spin-down panel
        plot(mainax, dos.e, -dos.ddos(:,i), dos.lstyle{i},'color', dos.color(i,:)); 
    end    
    plot(mainax, [0 0], [dos.dosmin dos.dosmax], '--k');  % plot vertical line for Fermi level
    plot(mainax, [dos.emin dos.emax], [0 0], '-k');      % plot zero line 
    axis(mainax, [dos.emin dos.emax dos.dosmin dos.dosmax]); 
                
    % plot inset
    if(dos.plotinset)
        for i=dos.sort{1}   % spin-up panel   
            plot(insetax, dos.e, dos.ddos(:,i), dos.lstyle{i},'color', dos.color(i,:)); 
        end      
        for i=dos.sort{2}   % spin-down panel
            plot(insetax, dos.e, -dos.ddos(:,i), dos.lstyle{i},'color', dos.color(i,:)); 
        end             
        plot(insetax, [dos.inset_emin dos.inset_emax], [0 0], '-k');      % plot horizontal line
        plot(insetax, [0 0], [dos.inset_dosmin dos.inset_dosmax], '--k');  % plot vertical line for Fermi level
        axis(insetax, [dos.inset_emin dos.inset_emax dos.inset_dosmin dos.inset_dosmax]);         % range
        % allow for some space for the inset on the main axes
        axis(mainax, [dos.emin dos.emax dos.dosmin dos.insetshift*dos.dosmax]); 
        set(insetax,'xtick',[dos.inset_emin:1:dos.inset_emax]);
    end;
end;       
    
% after plotting, remove the ticks of the axes of the inset
if(dos.plotinset)
    set(insetax,'ytick',[ ]);                        
end;  

% the legends have to be resorted in the order they are plotted, 
% i.e as defined in 'dos.sort',  the order of 'curves' is usually different
if (~isempty(dos.legend))
    doslegend = {dos.legend{dos.sort{1}}};
    if (~singlepanel)
        doslegend = {doslegend{:}, dos.legend{dos.sort{2}}};
    end

    % plot the legend
    %legend(mainax, doslegend, 'Location', 'EastOutside'); 
    legend(mainax, doslegend, 'Location', 'Best'); 
end


        




