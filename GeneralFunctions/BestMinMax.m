function [CurrentBestMinMax] = BestMinMax(Data,CompareTo)

[minValue,closestIndex] = min(abs(Data-CompareTo));
CurrentBestMinMax = [Data(closestIndex), (abs(min(Data)-CompareTo)+abs(max(Data)-CompareTo))/2]*100;

end

