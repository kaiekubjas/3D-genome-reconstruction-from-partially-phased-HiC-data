function vis(A,title_text)

imagesc(log10(A));
set(gca,'YDir','normal')
colormap(flipud(gray))
colorbar
set(gca,'DataAspectRatio',[1 1 1])

if nargin == 1
    title(inputname(1))
else
    title(title_text)
end

end