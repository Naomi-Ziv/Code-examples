function [] = parsave(frames, home, name, name2)

save(strcat(home,'/',name2,name,'.mat'),'frames','-v7.3')

end