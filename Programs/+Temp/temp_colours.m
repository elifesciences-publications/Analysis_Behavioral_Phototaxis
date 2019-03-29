function[colour] = temp_colours(varargin)

c1 = colour_palette(0,1);
c2 = colour_palette(0,2);
c3 = colour_palette(0,3);
c4 = colour_palette(0,4);
colour = [c1(4,:); c2(4,:) ; c3(4,:); c4(4,:)];

if nargin > 0
    strcmp(varargin{1}, 'dark')
    colour = [c1(3,:); c2(3,:) ; c3(3,:); c4(3,:)];
end