function VisualizePolicy(sysAbs,pol,DFA,varargin)
    %VISUALIZE Summary of this function goes here
    %   Detailed explanation goes here
    if sysAbs.dim~=2 
        warning('!! Visualization is not possible ')
    end
        
    if ~isempty(varargin)
        try
            sys = varargin{1};
        for P = sys.regions
            plot(P1, 'LineWidth',2)
            axis equal
            hold on
            plot(P2,'LineWidth',2)
        end
        catch
        end
    end
    
    
    if ~iscell(sysAbs.hx)  
        l = [length(sysAbs.hx),length(sysAbs.hx)];
         hx1 = sysAbs.hx;
         hx2 = sysAbs.hx;
         %% Note due to order of computations in gridspace_2d 
         [X2hat, X1hat] = ndgrid(hx2,hx1);

    else
         l = [length(sysAbs.hx{1}),length(sysAbs.hx{2})];
         hx1 = sysAbs.hx{1};
         hx2 = sysAbs.hx{2};
                  [X1hat, X2hat] = ndgrid(hx1,hx2);
    end
    if length(size(pol))>2
    Xn = sysAbs.f_det(sysAbs.states, pol(:,:,DFA.S0))-sysAbs.states;
    else
    Xn = sysAbs.f_det(sysAbs.states, pol(:,:))-sysAbs.states;
    end
    quiver(X1hat,X2hat,reshape(Xn(1,:),l(1),l(2)),reshape(Xn(2,:),l(1),l(2)))
    ylabel('x_2')   
    xlabel('x_1')
    title('Policy at DFA state S0')
end