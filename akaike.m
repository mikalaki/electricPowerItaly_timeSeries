function [aicS] = akaike(x, Tmax)
aicS = zeros(6,6);
for p=1:7
    for q=1:7
        [~,~,~,~,aicS(p, q),~]=fitARMA(x, p-1, q-1, Tmax);
    end
end
figure()
for i=1:7
    plot(aicS(:,i))
    hold on
end
end

