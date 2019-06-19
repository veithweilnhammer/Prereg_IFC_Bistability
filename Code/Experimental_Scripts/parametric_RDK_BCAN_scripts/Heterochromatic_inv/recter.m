function rect = recter(siz, coord, Rect, shift)

if nargin == 3
    shift = [0,0];
elseif nargin == 2
    Rect = get(0, 'ScreenSize');
    shift = [0,0];
elseif nargin == 1
    coord = [0.5,0.5];
    Rect = get(0, 'ScreenSize');
    shift = [0,0];
end

if size(siz,2) == 1
    siz = [siz,siz];
end

if coord(1)>=1
    rectx1 = -siz(:,1)/2+shift(:,1)+coord(:,1); 
    recty1 = -siz(:,2)/2+shift(:,2)+coord(:,2);
    rectx2 = siz(:,1)/2+shift(:,1)+coord(:,1);
    recty2 = siz(:,2)/2+shift(:,2)+coord(:,2);
else
    rectx1 = -siz(:,1)/2+shift(:,1)+Rect(3)*coord(:,1); 
    recty1 = -siz(:,2)/2+shift(:,2)+Rect(4)*coord(:,2);
    rectx2 = siz(:,1)/2+shift(:,1)+Rect(3)*coord(:,1);
    recty2 = siz(:,2)/2+shift(:,2)+Rect(4)*coord(:,2);
end
rect = [rectx1,recty1,rectx2,recty2];

if size(rect,1) > 1
    rect = rect';
end