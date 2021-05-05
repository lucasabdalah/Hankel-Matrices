function hfig = plotpoles(r,large_size,received_symbol,write_bits);
    % function hfig = plotpoles(M,r);
    % Computa a codificao de Gray a partir de uma sequencia de bits de tamanho K.
    %
    % SYNTAX: hfig = scatter_MQAM(r,large_size,received_symbol);
    %
    % INPUTS: 
    %       r : complex number symbol
    %       large_size : to enlarge the circle in the plot
    %       received_symbol : to plot a received symbol
    % 
    % OUTPUTS:
    %       hfig : saved figure properties
    %
    %HISTORY:
    % 2021/04/28: - Lucas Abdalah.
    %
    
    %% Color definitions for graphics
    yellow = [0.9290 0.6940 0.1250];
    black  = [0 0 0];
    M = size(r,2);
    hold on
    
    % Draw a circle radius = d
    % center = [0,0];
    % radii = abs(max(r));
    % viscircles(center,radii,'Color','b',...
    % 'LineStyle',':',...
    % 'LineWidth', 1.0);
    
    %% Scatter plot
    if received_symbol == false
        hfig = scatter(real(r),imag(r), 'filled',...
        'MarkerEdgeColor',yellow,... 
        'MarkerFaceColor',yellow);
        xlabel('Real');
        ylabel('Im');
        title(['Constelacao ',num2str(M),'-PSK']);
        hold on;
        ax = 1.2*abs(max(r));
        line([0,0],[-ax,ax], 'Color', black,...
        'HandleVisibility','off');
        line([-ax,ax],[0,0],'Color', black, ...
        'HandleVisibility','off');
    end
    
    if received_symbol == true
        symbol_color=[204 51 0]./255;
        hfig = scatter(real(r),imag(r), 'filled',...
        'Marker','x',...
        'MarkerEdgeColor',symbol_color,... 
        'MarkerFaceColor',symbol_color);
        hfig.SizeData = 20;
    end
    
    if large_size == true
        hfig.SizeData = 200;
    end
    
    if write_bits == true
        % Generating Gray alphabet and constellation's plot.
        for jj = 1:M
            gray_alfabeto = gray_const(M,false);    
            str = {strjoin(string(gray_alfabeto(jj,:)))};
            text(real(r(jj)),imag(r(jj))+0.1*sqrt(2),...
            str,...
            'FontSize', 6,...
            'HorizontalAlignment', 'center');
        end
    
    end
    
    axis square
end