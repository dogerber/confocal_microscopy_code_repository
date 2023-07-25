function out_logical = input_YN(prompt)

if nargin == 0
prompt = '[Y]es / [N]o ?\n';
end

answer = input(prompt,'s');

if answer == 'Y' || answer == 'y'
    out_logical = true;
else
    out_logical = false;
end



end