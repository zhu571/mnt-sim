====================================================================
该程序结合经验的势垒分布计算熔合截面以及HIVAP程序计算复合核存活几率

请参阅文献 Phys. Rev. C 77, (2008) 014603； Nucl.Phys.A768 (2006) 80

Ning Wang,  2008-08-21，

wangning@gxnu.edu.cn;     ningmany@gmail.com
===================================================================

该程序已经将最大的分波修改为 300. 原来是100.  By C.W.Shen

Originally the L-max is 100.

to change the L-max, you should change the following.

hivap.f          Change the JDIM to that you want.   XXX
common.f         SPECJ(XXX),FLAND(XXX),YFISJ(XXX),YPOPJ(XXX),YERJ(XXX),
                 YFISS(XXX),PER(XXX),PFIS(XXX)
change All of the EJAY to EJAY(5*XXX) in all common block /YR/
change All of the YRSMO to YRSMO(5*XXX) in all common block /DNS/
change ARRAY, YJK to (XXX) in densty.f

－－－－－－－－－－－－－－－－－－－－－－－－－－
2007年5月23日：Ning Wang and Kai Zhao
(1) 这组程序的最大分波为L=300.

(2) 液滴模型的裂变位垒的计算由hivapein.dat中的 FISROT PARAMETER 控制。新
    增加的一个是经验拟合的液滴裂变位垒, 即将 fisrot parameter设为3. 该公式
    取自：Phys. Rev. C 77, 014603 (2008)

(3) 由于Esh=Eexp-Eld, 已经包含了对能，因此在计算中，需要将pairing设置为0.

(4) 熔合（俘获）截面由经验的位垒分布计算，M.Liu, N Wang, et al., Nucl.Phys.A768(2006)80.

(5) 能级密度参数采用了标准的HIVAP（1981）参数 (在seqin.f中).

上述 (2)-(4)的修改，几乎都在sepenf.f中进行。由于读取fisrot parameter的需要，
修改了hivap.f及其调用sepenf.f的参数表（增加了一项）。

fusion部分没有作任何修改。只是在计算Sigma_{fus}时，使用另外一组新的程序linkfus。
同时和2002年从日本拷贝的hivap源文件进行了比较。

这个程序没有考虑准裂变效应，所以对于能够形成超重核的反应体系需要额外计算复合核形成几率PCN. 

------------ Ning Wang, 2008-08-21 --------------
为了方便使用，对输入输出文件做了修改：

1) 输入文件：input.dat  （从HIVAP.F中读入）

      需要输入反应体系弹核与靶核的质量和电荷：A1, Z1, A2, Z2 

      另外还需要输入质心入射能量 Ecm， 最后输入0，代表能量输入结束   

2) 输出文件: SIGMA.DAT 
       从左到右各列数据依次表示：质心系能量、俘获（熔合）截面、总的蒸发残余截面、裂变截面、俘获截面上限、俘获截面下限
             SIGXPN.DAT
       从左到右各列数据依次表示：蒸发质子数、蒸发中子数、质心系能量、俘获（熔合）截面、该反应道的蒸发残余截面、裂变截面
 
--------------------------------------------------------------------------

Mainly changed program: SEPENF.F
another subroutine "linkfus.f" is added to calculate the fusion cross section  
(1) The upper limit of the maximal angular momentum  L=300

(2) For the parameter "FISROT PARAMETER", more option, 3, could be set. If 
    this parameter is set to 3, the empirical formular for macroscopic fission
    barrier will be included. Phys. Rev. C 77, 014603 (2008)

(3) Because in fission barrier calculation, Esh=Eexp-Eld, the pairing energy is already
    included, pairing calculation option should be set to zero.

(4) The fusion cross section is calculated by the method in reference NPA768(2006)80

(5) The level density parameter proposed by Reisdorf in 1981 is adopted.

-------------------------------------------------------------------------
some other changes: (By C.shen, K.Zhao and N.Wang)

(1) hivap.f, we add "CALL LINKFUS()" to incorporate our formula for \sigma_{fus}.

(2) sepenf.f, (IFISRT.NE.2) ---> (IFISRT.EQ.0) for Liquid-drop barrier
                                 (IFISRT.EQ.2) for Sierk's barrier
                                 (IFISRT.EQ.3) for barrier defined by user
    some other models for fission barrier are added and tested.

(3) rowout.f, add "write(1236, )" to output the cross sections.

(4) rot.f, "DM0=BFCPS" is added to avoid the overflow caused by DM0=0 for super-heavy nuclei.
  
(5) msben.f, "CMASS=SMASS+SHLL" is added to some null return value for some nuclei, such as (nuclues with Z=113,A=285)
(6) grids.f, "M=0", "N=0", "Exmax=0" are added by Kai Zhao.

(7) fusio.f, the error function ERFC(x) is added.

(8) file.f, "linkfus.f" is added.

(9) alloc.f, some additional files are created by N.Wang.

(10) to search for "!", one can find the main changes by N. Wang. 

-------------------------------------------------------------------------

------------ Ning Wang, 2008-08-21 --------------
To do our calculations more convenient, we add additional input and output files 

1) input file：input.dat  (use in HIVAP.F)

      input the mass and charge of the reaction system：A1, Z1, A2, Z2 

      input the center-of-mass energy: Ecm  ， '0' to finish Ecm input;     

2) output files: SIGMA.DAT 
       From left to right：center-of-mass energy、capture (fusion) cross section, total evaporation residue cross section, fission cross section, upper limit of the capture cross section,  lower limit of the capure cross section;

             SIGXPN.DAT
       From left to right：number of evaporated protons, number of evaporated neutrons, center-of-mass energy, capture cross section, evaporation cross section at this channel, fission cross section

----------------------------------------------------






