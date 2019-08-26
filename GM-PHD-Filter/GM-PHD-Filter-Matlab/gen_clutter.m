function clutter = gen_clutter(model)

% This generates clutter that will be appended to the actual observations
% of targets.
nClutter = model.lambda_t;
xrange = model.xrange;
yrange = model.yrange;
clutter = zeros(2, nClutter);
for i = 1:nClutter
    clutterX = rand * (xrange(2) - xrange(1)) + xrange(1); % Random number between xrange(1) and xrange(2), uniformly distributed.
    clutterY = rand * (yrange(2) - yrange(1)) + yrange(1); % Random number between yrange(1) and yrange(2), uniformly distributed.
    
    clutter(1,i) = clutterX;
    clutter(2,i) = clutterY;
end