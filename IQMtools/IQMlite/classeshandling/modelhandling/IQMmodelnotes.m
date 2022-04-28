function notes = IQMmodelnotes(model)
% IQMmodelnotes: displays the notes in model.
%
% USAGE:
% ======
% [notes] = IQMmodelnotes(model)
%
% model: IQMmodel description of model
%
% Output Arguments:
% =================
% notes: notes, stored in the model

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

ms = struct(model);
notes = ms.notes;


