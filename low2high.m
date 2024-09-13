% Update g.x using parameter vector of non-linear least square in TDOA-GN
function [g] = low2high(low,g)
if g.label=="final"
    g.x(1,8)=low(1);
    g.x(2:g.M,1:end)=reshape(low(2:g.m*(g.M-1)+1),g.m,g.M-1)';
    g.x(g.M+1:end,1:3)=reshape(low(g.m*(g.M-1)+2:end),g.n,g.K)';
elseif g.label=="init"
    g.x(1,8)=low(1);
    for i=2:g.M
        g.x(i,1:3)=low(1+(g.m-3)*(i-2)+1:1+(g.m-3)*(i-2)+3)';
        g.x(i,7:8)=low(1+(g.m-3)*(i-2)+4:1+(g.m-3)*(i-2)+5)';
    end
    g.x(g.M+1:end,1:3)=reshape(low((g.m-3)*(g.M-1)+2:end),g.n,g.K)';
end
