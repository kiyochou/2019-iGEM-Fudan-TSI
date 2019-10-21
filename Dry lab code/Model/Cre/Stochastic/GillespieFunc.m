function [tvec, xmat] = GillespieFunc(x, Pre, Post, c, iter)
    tt = 0;
    S = (Post-Pre)';                                              % Reaction matrix
    [u, v] = size(S);                                               % u -> #Species, v -> #Reactions
    tvec = zeros(iter+1,1);                                    % Record the time series
    xmat = zeros(iter+1, u);                                  % Record the substance number series
    xmat(1,:) = x;                                                    % Initial substance number
    
    for i = 2:iter
        RateConst = RateConstFunc(x, tt, c);                                       %  Compute the instantaneous rate constant.
        tt = tt + exppdf(rand()*10000, sum(RateConst));                    % Time of the next occurring event.
        Prob = RateConst./sum(RateConst);                                          % Probability distribution calculation
        j = randsrc(1,1,[1:v; Prob]);                                                         % Simulate the type of the occurring event at moment tt.
        x = x+S (:,j)';                                                                              % Change the substance number
        tvec(i) = tt;                                                                                % Record the event occurring moment
        xmat(i+1, :) = x;                                                                           % Record change of substance number
    end
    
end