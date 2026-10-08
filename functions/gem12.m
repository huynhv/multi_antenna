function rgb = gem12(name)
%GEM12 Return RGB value for a named MATLAB gem12 color.
%
%   rgb = gem12(name)
%
%   NAME can be:
%        "blue"
%        "orange"
%        "yellow"
%        "purple"
%        "green"
%        "lightblue"
%        "magenta"
%        "cyan"
%        "lightorange"
%        "lightgreen"
%        "lightpurple"
%        "lightred"
%
%   RGB is returned as a 1x3 vector with values in [0,1].

    C = orderedcolors("gem12");

    names = [
        "blue"
        "orange"
        "yellow"
        "purple"
        "green"
        "lightblue"
        "magenta"
        "cyan"
        "lightorange"
        "lightgreen"
        "lightpurple"
        "lightred"
    ];

    idx = find(strcmpi(string(name), names), 1);

    if isempty(idx)
        error('Unknown gem12 color: "%s"', name);
    end

    rgb = C(idx,:);
end