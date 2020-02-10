function [aicS] = akaike(x, Tmax)
aicS = zeros(6,6);
for p=1:7
    for q=1:7
        [~,~,~,~,aicS(p, q),~]=fitARMA(x, p-1, q-1, Tmax);
    end
end
figure()
for i=1:7
    plot([0:6],aicS(:,i),'-o');
    ylabel("AIC");
    xlabel("q");
    hold on
end
hold off
legend("p=0","p=1","p=2","p=3","p=4","p=5","p=6")
end

