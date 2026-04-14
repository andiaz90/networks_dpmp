yyy = ylim ;

zzz = yyy;

if abs(zzz(1))<0.1
    zzz(1)=-0.2;
end
if abs(zzz(2))<0.1
    zzz(2)=0.2;
end


if zzz(2)-zzz(1)<0.01
    zzz(2)=zzz(2)+0.01
    zzz(1)=zzz(1)-0.01
end

ylim([zzz(1) zzz(2)]);
