function [] = parsave2(rprop, home, name, name2)

save(strcat(home,'/',name2,name,'.mat'),'rprop','-v7.3')

end