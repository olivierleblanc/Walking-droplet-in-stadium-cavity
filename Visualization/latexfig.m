function ax = latexfig(gca, pixel2mm)
%
    set(gca,'TickLabelInterpreter','Latex', 'FontSize', 14.0);
    ax = gca;
    xTick = get(ax,'XTick');
    set(ax,'XTick',downsample(xTick,2));
    yTick = get(ax,'YTick');
    set(ax,'YTick',downsample(yTick,2));
    ax.XTickLabel = round(ax.XTick.*pixel2mm.*1e3)./1e3;
    ax.YTickLabel = round(ax.YTick.*pixel2mm.*1e3)./1e3;
    arr = str2num([ax.XTickLabel]);
    arr = arr-arr(1);
    ax.XTickLabel = num2cell(arr);
    arr = str2num([ax.YTickLabel]);
    arr = arr-arr(1);
    ax.YTickLabel = num2cell(arr);
end

