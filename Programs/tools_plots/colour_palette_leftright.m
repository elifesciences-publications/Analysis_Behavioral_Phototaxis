function[left_color, right_color] = colour_palette_leftright(show, comp)

% color palette
% paletton.com

% *** Primary color:
left_color(1,:) = [0.271,0.098,0.188];
left_color(2,:) = [0.671,0,0.353];
left_color(3,:) = [0.49,0,0.259];
left_color(4,:) = [0.071,0,0.035];
left_color(5,:) = [0.016,0,0.008];

% *** Secondary color:
if comp    
    right_color(1,:) = [0.176,0.267,0.459];
    right_color(2,:) = [0.227,0.047,0.561];
    right_color(3,:) = [0.169,0.039,0.408];
    right_color(4,:) = [0.024,0.004,0.059];
    right_color(5,:) = [0.004,0,0.012];
else
    right_color(1,:) = [0.333,0.204,0.122];
    right_color(2,:) = [0.82,0.333,0];
    right_color(3,:) = [0.6,0.243,0];
    right_color(4,:) = [0.086,0.035,0];
    right_color(5,:) = [0.02,0.008,0];
end

n = size(left_color,1);
if show
    figure;
    hold on
    for i = 1:n
        plot(i, 1, 'o', 'MarkerFaceColor', left_color(i,:),...
                        'MarkerEdgeColor', left_color(i,:),...
                        'MarkerSize', 12)
    end
    for i = 1:n
        plot(i, 2, 'o', 'MarkerFaceColor', right_color(i,:),...
                        'MarkerEdgeColor', right_color(i,:),...
                        'MarkerSize', 12)
    end
end