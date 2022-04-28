function [result] = preFillCharIQM(value2prefill,lengthString,fillChar)
% [result] = preFillCharIQM(value2prefill,lengthString,fillChar)

if isnumeric(value2prefill),
    result = [char(double(fillChar)*ones(1,lengthString-length(num2str(value2prefill)))) num2str(value2prefill)];
elseif ischar(value2prefill),
    result = [char(double(fillChar)*ones(1,lengthString-length(value2prefill))) value2prefill];
else
    error('Unknown type to prefill.');
end