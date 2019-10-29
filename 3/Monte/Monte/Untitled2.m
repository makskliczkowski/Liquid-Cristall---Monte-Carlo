f = fopen('graphT1.0.txt'); 
c = textscan(f,'%s','Delimiter','\n');fclose(f);
c2 = c{1}(~cellfun(@isempty,c{1}))
imshow(c2, []);