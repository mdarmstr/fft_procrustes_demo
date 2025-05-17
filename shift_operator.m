function img_shifted2 = shift_operator(matrix,shift)

% Read the image
img = matrix;
img = sum(img,3);

% Convert to grayscale if it's RGB
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Convert to double for accurate computation (optional, if needed)
img = im2double(img);

% Define first shift (e.g., vertical shift down)
img_shifted1 = circshift(img, [shift(1), 0]);  % shift 40 rows down

% Define second shift (e.g., horizontal shift right)
img_shifted2 = circshift(img_shifted1, [0, shift(2)]);  % shift 40 columns right

end
