% Define GR vector
GR = [0.01:0.001125:0.1, ...
      0.11125:0.01125:1, ...
      1.1125:0.1125:10, ...
      11.125:1.125:100];

GR = unique(GR, 'stable')';   % column vector

% Preallocate p_Ga2O3
p_Ga2O3 = zeros(size(GR));

% Threshold value
threshold = 0.00025;

% Flag to indicate threshold reached
threshold_reached = false;

for i = 1:length(GR)
    if ~threshold_reached
        p_Ga2O3(i) = GR(i)/4000;
        if p_Ga2O3(i) >= threshold
            p_Ga2O3(i) = threshold;   % cap at threshold
            threshold_reached = true;  % keep flat from now on
        end
    else
        % Once threshold reached, keep value constant
        p_Ga2O3(i) = threshold;
    end
end

% Display first few results
disp([GR, p_Ga2O3]);