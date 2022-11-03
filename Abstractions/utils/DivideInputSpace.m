function InputSpace = DivideInputSpace(InputSpace,Act_part,Fb_part);
% Divide input space InputSpace into part for actuation (Act_part)
% and part for feedback (Fb_part), corresponding to the interface function
% u = uhat+K(x-xhat). That is, uhat \in Act_part, K(x-xhat) \in Fb_part and
% we should have Act_part+Fb_part <=1, such that uhat \in InputSpace

%% Initialisation
% Check if division is correct
if Act_part + Fb_part > 1
    warning(['Your finite-state input space is larger than the original input space. Is this what you want? ...' ...
        'Please verify if Act_part and Fb_part are specified correctly.'])
end

% Check if both inputs are given
if isempty(Act_part) & isempty(Fb_part)
    Act_part = 1;
    Fb_part = 0;
elseif isempty(Act_part)
    Act_part = 1-Fb_part;
elseif isempty(Fb_part)
    Fb_part = 1-Act_part;
end

%% Divide input space    
InputSpace{2} = Act_part*InputSpace{1};
InputSpace{3} = Fb_part*InputSpace{1};

end


