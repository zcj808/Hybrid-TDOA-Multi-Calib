% transform g.x into the unknown parameter vector
function [low] = high2low(g)
    if g.label=="final"
        high=g.x(2:g.M,:)';
        low=high(:);
        high=g.x(g.M+1:end,1:3)';
        low=[low;high(:)];
        low=[g.x(1,8);low];
    elseif g.label=="init"
        low=[g.x(1,8)];
        for i=2:g.M
            low=[low;g.x(i,1:3)';g.x(i,7:8)'];
        end
        s_locs=g.x(g.M+1:end,1:3)';
        low=[low;s_locs(:)];
    end
end