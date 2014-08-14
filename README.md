Batch4SPM
=========
Author: hongshengcheng.math@gmail.com
If you have any question,please feel free to contact me

--------
Matlab script to run SPM in batch mode,Only for Task 
this script can convert the dcm file to nii file(3d);also can do preprocess and 1st level analysis

Note:
1.data should be organized like 'sub**/run*', i.e. each subjects' data should store in 1 folders named witd prefix "sub",and each run should like run1 run2.... 

2.Before run,you should add script in your matlab search path:Set path --> add folder

3.This script only tested on windows,linux or mac may also work well

4.Parallel mode may only support to convert the img format

5.Res Dir

RealignParameter headmotion parameters(which will be stored in folders named with SubID)

rp_Report 3 cut-off level report,txt file contains the bad scan order number


6.Brief guide to run:

How many workers to use?(Default:1) : 0 will initialize the maximum parallel worker number
                                      1 or just ENTER to skip will not trigger the parallel


Select Data Dir choose your data locationg,e.g.D:\DataProc\sub01，D:\DataProc\sub02，D:\DataProc\sub03..., then the DataDir should be D:\DataProc\

Other settings ref the manual of SPM
