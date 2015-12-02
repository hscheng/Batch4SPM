contrast.txt

在每个contrast之后加6个0【6个头动参数
比如有4个条件那么就是
con1 1 0 0 0 0 0 0 0 0 0 [con1 是这个contrast的名字]


Onset.xls
#如果每个被试的Onset设置是一样的
一个run一个sheet,每个条件就是一列
注意onset必须是从小到大

#如果每个被试的onset不一样
第一列就是被试编号，第二列就是一个条件的onset

*********************************************************
contrast.txt

In each contrast there should add 6 zeros(headmotion parameters)
so when you have 4 conditions,the contrast should be like: 
con1 1 0 0 0 0 0 0 0 0 0 [con1 is the name of the contrast]


Onset.xls

#if all subjects share the same onset setup,then
1 sheet 1 run(DONOT leave black sheet in xls file)
in each sheet every condition is the same order and 1 condition 1 colunm

# if onset differ one subject to another
data in each sheet should organize in the format
1st column is the sub_id,2nd column is onset,3rd col is sub_id,4th col is onset....
every 2 column is 1 condtion