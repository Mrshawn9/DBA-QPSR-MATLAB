function [curve, ptr] = RecordCheckpoint(curve, ptr, nfe, MaxFEs, best)

while ptr <= 10 && nfe >= round(ptr*MaxFEs/10)
    curve(ptr) = best;
    ptr = ptr+1;
end
end
