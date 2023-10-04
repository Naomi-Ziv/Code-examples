#!/usr/bin/env python

#
for name in ['F1','F6','F12','F12_clone1','F12_clone2','F12_clone3','C1S2','C1S4','C1S6','C3S2','C3S4','C3S6','C2S2','C2S4','C2S6','C4S1','C4S4','C4S6','C5S0','C5S4','C5S7','C6S0','C6S4','C6S6','F2SegPool','Oak','Vineyard']:
     with open('/Users/Naomi/snp/{name}.sort.snp'.format(name=name)) as f:
          x=f.readlines()
          f2=open('/Users/Naomi/snp/{name}.snp'.format(name=name),'w')
          for i in range(27,(len(x)-1)):
               y=x[i].split('\t')
               y2=y[7].split('DP4=')[1].split(';')
               y3=y2[0].split(',')
               f2.write('{a}\t{b}\t{c}\t{d}\t{e}\t{g}\t{h}\t{i}\t{j}\n'.format(a=y[0],b=y[1],c=y[3],d=y[4],e=y[5],g=y3[0],h=y3[1],i=y3[2],j=y3[3]))
          f2.close()


   
