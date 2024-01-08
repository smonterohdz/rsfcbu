function Hb_brain = GlobalRegression(Hb_brain)
A = mean(Hb_brain,1);
A_inv = A'/(A*A');
for i = 1:size(Hb_brain,1)
    y = Hb_brain(i,:);
    b = y*A_inv;
    Hb_brain(i,:) = y - A*b;
end
end
