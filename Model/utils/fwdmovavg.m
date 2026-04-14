function y = fwdmovavg(x,window)

y = x;

for t=1:length(x)
   
    if t<length(x)-window
    y(t) = nanmean(x(t:t+window-1));
    else
        y(t) = nanmean(x(t:end));
    end
    
end

