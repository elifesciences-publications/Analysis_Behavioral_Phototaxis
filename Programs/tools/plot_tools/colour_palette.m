function[colour] = colour_palette(show, colourID)

% color palette
% paletton.com

% *** Complementary color:
complementary_color(1,:) =[0,0.439,0.408];
complementary_color(2,:) = [0.008,0.153,0.141];
complementary_color(3,:) = [0.024,0.341,0.318];
complementary_color(4,:) = [0.106,0.518,0.486];
complementary_color(5,:) = [0.247,0.62,0.592];

% *** Secondary color (1):

secondary_color1(1,:) = [0.545,0.686,0];
secondary_color1(2,:) =[0.192,0.239,0.012];
secondary_color1(3,:) = [0.431,0.533,0.039];
secondary_color1(4,:) = [0.675,0.808,0.165];
secondary_color1(5,:) = [0.843,0.961,0.38];

%*** Secondary color (2):

secondary_color2(1,:) = [0.659,0,0.141];
secondary_color2(2,:) =[0.227,0.012,0.059];
secondary_color2(3,:) =[0.514,0.035,0.141];
secondary_color2(4,:) =[0.776,0.157,0.29];
secondary_color2(5,:) =[0.925,0.369,0.486];

%*** Complement color:

primary_color(1,:) = [0.714,0.463,0];
primary_color(2,:) = [0.247,0.165,0.012];
primary_color(3,:) = [0.557,0.376,0.039];
primary_color(4,:) = [0.839,0.604,0.173];
primary_color(5,:) = [1,0.788,0.396];


n = size(primary_color,1);
if show
    figure;
    hold on
    for i = 1:n
        plot(i, 1, 'o', 'MarkerFaceColor', primary_color(i,:),...
                        'MarkerEdgeColor', primary_color(i,:),...
                        'MarkerSize', 12)
    end
    for i = 1:n
        plot(i, 2, 'o', 'MarkerFaceColor', secondary_color1(i,:),...
                        'MarkerEdgeColor', secondary_color1(i,:),...
                        'MarkerSize', 12)
    end
    for i = 1:n
        plot(i, 3, 'o', 'MarkerFaceColor', secondary_color2(i,:),...
                        'MarkerEdgeColor', secondary_color2(i,:),...
                        'MarkerSize', 12)
    end
    for i = 1:n
        plot(i, 4, 'o', 'MarkerFaceColor', complementary_color(i,:),...
                        'MarkerEdgeColor', complementary_color(i,:),...
                        'MarkerSize', 12)
    end
end

if colourID == 1
    colour = primary_color;
elseif colourID == 2
    colour = secondary_color1;
elseif colourID == 3
    colour = secondary_color2;
else
    colour = complementary_color;
end