% Sequence Task Demo

plotnum = 1:20;
count = 1;
center = 5;
overlap = [3,3];
clr = ['rgbc'];
for i = 1:4;
    for ii = 1:5;
        subplot(4,5,plotnum(count));
        if ii == 1;
            plot(center,center,'+w','markersize',5)
        elseif ii == 3
            plot(overlap(1),overlap(2),[clr(i) 's'],'markersize',5);
        else
            r = rectangle('Position',[randi(9) randi(9) 0.5 0.5]);
            set(r,'edgecolor',clr(i));
            set(r,'Facecolor',clr(i));
        end
        set(subplot(4,5,plotnum(count)),'Color','k');
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis([0 10 0 10])
        count = count+1;
    end
end
subtitle('Sequence Task')