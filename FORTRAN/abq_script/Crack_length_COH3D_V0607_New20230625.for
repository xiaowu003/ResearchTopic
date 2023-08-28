      subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
C constitutive response of cohesive elements using a traction-separation
C     更新日期：2021年2月1日
c    2019年9月12日-----更新能量计算（拆分计算）
C    2019年11月14日----------修正Quads 二次名义应力 损伤判断推迟一步bug      ---完成日期-------2019年11月14日      
C    2019年11月14日----------更换压剪本构算法（罚刚度不变）                  ---完成日期-------2019年11月21日
C    2020年4月13日-----------单元不删除，完全破坏后模拟摩擦行为，调整整体结构---完成日期-------2020年4月30日
C    2020年5月4日 -----------修改压力变化过程的切应力变化(有位移才有应力-无过程)方案1---完成日期------2020年5月22日
C    2020年12月4日 ----------- 更正sigm_shear  = SQRT(sigm_s**2 + sigm_t**2 )不影响主体程序---完成日期------2020年12月4日  
C    2021年1月6日 ----------- 简单将切向位移或者法向位移大于0.5时的 ，所有应力归零 
C    2021年2月1日 ----------- 添加用于实体单元的 VUSDFLD ,用于删除飞溅的实体单元 ---完成日期------2021年2月6日  
C    2021年2月24日 ----------- 当材料名称包含“VUMAT-friction”只执行摩擦本构 ---完成日期------2021年2月24日
C    2021年3月12日 ----------- 增加裂纹模式判断,-1代表未断裂，0-90代表混合角度，90剪切，0拉伸，---完成日期------2021年3月12日
C    2021年6月7日 ----------- 增加初始时刻裂纹模式判断,-1代表未断裂，0-90代表混合角度，90剪切，0拉伸，---完成日期------2021年6月12日  
C    2021年6月18日 ----------- 增加任意时刻裂纹模式判断,-1代表未断裂，0-90代表混合角度，90剪切，0拉伸，---完成日期------2021年6月18日        
C    2022年10月25日 ----------- 增加状态变量STATE(*,14),用于SDEG>=0.99时的裂纹模式判断,-1代表未断裂，0-90代表混合角度，90剪切，0拉伸，---完成日期------2022年10月25日
C    2022年10月26日 ----------- 增加状态变量STATE(*,15),用于统计SDEG>=0.99时的裂纹面积，---完成日期------2022年10月26日
C    2023年6月25日 -----------  增加输入参数Sdeg_L，提高程序适用性。Sdeg_L可以是任意SGED值，用于判断和计算内聚力单元 >=Sdeg_L 时的裂纹类型和裂纹面积，---完成日期------2023年6月25日
C      
C*******************************************    状态变量数：13   **************************
C*******************************************    删除标志号：10   **************************
C ******************************************  材料参数数量：10   **************************
C delt_位移   sigm_应力   epsil_应变   
C The state variables are stored as:
C      STATE(*,1) = deletion flag 删除标志  =0删除  代表内聚力本构结束
C      STATE(*,2) = strain component 33
C      STATE(*,3) = strain component 23
C      STATE(*,4) = strain component 31
C      STATE(*,5) = DMICRT 初始损伤准则标志 >=1开始损伤
C      STATE(*,6) = SIGM_TMAX 界面抗剪强度
c      STATE(*,7) = SDEG 损伤程度0-1 =1完全损伤
c      STATE(*,10) = delete_flag  =0删除  代表单元结束 不在计算摩擦力
c      STATE(*,11) = fracture_mode计算断裂时的断裂类型，仅在完全破坏后的下一增量步计算，更新一次后不再计算 =90.d0为纯剪切，0-90为混合角度比
c      STATE(*,12) = fracture_mode_initial计算初始断裂时的断裂类型，仅在初始破坏后的下一增量步计算，更新一次后不再计算 =90.d0为纯剪切，0-90为混合角度比
c      STATE(*,13) = fracture_b 计算任意时刻的断裂类型， =90.d0为纯剪切，0-90为混合角度比，每一个增量步都更新
c      STATE(*,14) = fracture_mode_SDEG 计算SDEG=0.99的断裂类型， =90.d0为纯剪切，0-90为混合角度比，每一个增量步都更新
c      STATE(*,15) = delete_flag099 用于计算SDEG大于0.99的内聚力单元的裂纹面积
    
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(nblock), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),
     4  stretchOld(*), defgradOld(*),
     5  fieldOld(*), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(*), fieldNew(*),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
      parameter(pi = 3.14159265358979323D0,crit_lens = 1.18D0
     1   ,crit_lenn = 1.1D0)
      INTEGER D_MODE
         
c
C ----------------------------提取材料参数
C
      EN      = props(1)      !杨氏模量
      ES      = props(2)      !切变模量
      SIGM_NMAX = props(3)    !界面抗拉强度
      C       = props(4)      !内聚力(当PHI=0时为界面抗剪强度）
      PHI     = props(5)      !内摩擦角/deg
      GIC     = props(6)      !I型断裂能
      GIIC    = props(7)      !II型断裂能
      ETA     = props(8)      !混合模式系数
      D_MODE  = props(9)      !初始破坏准则标志 M-C=0  Quads=2 
      FMU_c   = props(10)     !完全破坏后的 摩擦系数
      Sdeg_L  = props(11)     !输入任意SDEG值，可用于判断输入SDEG值下的裂纹类型和统计裂纹面积。
C
      ET      = ES                !切变模量
      FMU_i   = tan(PHI/180*PI)   !摩擦角等效摩擦系数--初始破坏处
      FMU_r   = FMU_i             !摩擦角等效摩擦系数--完全破坏处
C-------------------------------------------------------------
      do 100 i = 1,nblock
C-----------------------------进入循环
        pk_n   =EN/charLength(i)                       !罚刚度
        pk_s   =ES/charLength(i)
        pk_t   =ET/charLength(i)
C        
        !计算各项应变以及位移
        epsil_n = stateOld( i,2) + strainInc( i,1)          ! strain component 33
        epsil_s = stateOld( i,3) + strainInc( i,2)          ! strain component 32
        epsil_t = stateOld( i,4) + strainInc( i,3)          ! strain component 31
c        
C
        delt_n = epsil_n * charLength(i)               !位移
        delt_s = 2*epsil_s * charLength(i)
        delt_t = 2*epsil_t * charLength(i)
        delt_shear = SQRT(delt_s**2+delt_t**2)        
C       
        !计算初始断裂位移和应变、完全破坏位移
        SIGM_fr  = pk_n*DIM(-delt_n,0.D0)*FMU_i                !残余摩擦应力--初始破坏处
        SIGM_fi  = pk_n*DIM(-delt_n,0.D0)*FMU_i    !!!初步修改,需进一步调试!!!!
        SIGM_TMAX= C + SIGM_fi                                !仅在压剪和纯剪模式下有意义   
        SIGM_SMAX= SIGM_TMAX                                  !剪切强度
C       
        delt_n0 = SIGM_NMAX/(pk_n)                            !初始断裂位移，不适用于拉剪混合模式
        delt_s0 = SIGM_SMAX/(pk_s)
        delt_t0 = SIGM_TMAX/(pk_t)
        delt_sc = SIGM_fr  /(pk_s)
C          
        epsil_n0 = SIGM_NMAX/EN                               !初始断裂应变
        epsil_s0 = SIGM_SMAX/ES
        epsil_t0 = SIGM_TMAX/ET
C        
        delt_nf= 2.d0*GIC/SIGM_NMAX                           !完全破坏位移
        delt_sf= 2.d0*GIIC/C+delt_sc
        delt_tf= delt_sf
C        
        DMICRT    = stateOld(i,5)                             !继承上一步变量
        SDEG      = stateOld(i,7)                             
        fracture_mode = -1.d0                                 !赋初值  -1表示未完全损伤
        fracture_mode_initial = -1.d0
        fracture_b = -1.d0                                    !fracture_b 表示任意时刻损伤的类型
        fracture_mode_SDEG = -1.d0                            !fracture_mode_SDEG表示SDEG=0.99时候的损伤类型。
        !delt_mf   =0.d0              ！调试
        !delt_m0   =0.d0
        !delt_m    =0.d0
        !SDEG_NEW  =0.d0
C
C-------------------------------------------------------------         
C-------------------------罚刚度混合--------------------------
C
        if (delt_s==0.d0 .and. delt_t==0.d0) then           !避免纯压缩时 分母为零 该值为NaN
          penaltyK_all = pk_n
        else
          penaltyK_all =(pk_n*(DIM(delt_n,0.d0))**2 +         !混合模式罚刚度
     1    pk_s*delt_s**2  +  pk_t*delt_t**2) /
     2    ((DIM(delt_n,0.d0))**2+delt_s**2+delt_t**2)
        end if
C        
        penaltyK_shear = pk_s                               !先不混合，避免当纯压缩时分母为0
C
C--------------------------------------------------------------------------
        !当材料名称包含“VUMAT-FRICTION”只执行摩擦本构 ,注意abaqus会将材料名称cmname全部转为大写字母
        
        if ( index(cmname,'VUMAT-FRICTION') >= 1 ) then
            SDEG =1D0
        end if
C--------------------------------------------------------------------------       
        
        if( SDEG>=1D0 ) then
          SDEG=1D0
C-----------------------------初始损伤判断
C
        else if ( DMICRT < 1.d0 .and.SDEG < 1D0) then                 !判断单元是否已经损伤
            !write(6,*)"no damage"
C              
C       ----------------------若未损伤，根据计算上一增量步的DMICRT
          select case (D_MODE)                                !选着初始损伤判断准则
          case (0)                                            !M-C 纯摩尔库仑
              SIGM_TMAX    = C - pk_n*delt_n*FMU_i             !判断损伤时，coh必处于弹性阶段
              DMICRT_n     = pk_n*DIM(delt_n,0.d0)/SIGM_NMAX  !未损伤时，按弹性计算。
              DMICRT_shear = (penaltyK_shear*delt_shear)/
     1                       SIGM_TMAX 
              DMICRT=MAX(DMICRT_n,DMICRT_shear)
          case (1)                                            !Maxe  最大名义应变 ||未使用||
               DMICRT=0.5d0
      !         DMICRT=MAX(DIM(stateOld( i,2),0.d0)/epsil_n0,  stateOld
      !1         (i,3)/epsil_s0, stateOld( i,4)/epsil_t0)
          case (2)                                            !Quads 二次名义应力
              DMICRT=(pk_n*DIM(delt_n,0.d0)/SIGM_NMAX)**2 +
     1    (pk_s*delt_s/SIGM_SMAX)**2+(pk_t*delt_t/SIGM_TMAX)**2   
          case (3)                                            !Maxs  最大名义应力 ||未使用||
              DMICRT=MAX(DIM(stressOld(i,1),0.d0)/SIGM_NMAX,
     1    stressOld(i,2)/SIGM_SMAX, stressOld(i,3)/SIGM_TMAX)
          end select
C
C         -------------------判断上一增量步是否损伤  
          if (DMICRT >= 1.d0) then                              !若损伤，根据混合模式计算sdeg值
            DMICRT = 1.d0                                       !将dmicrt置为1，保证其在[0,1]范围内

            !write(6,*)"----------damage"
            select case (D_MODE)                                !选择破坏模式
            case (0)
              call  BK_M_C(delt_n0,delt_s0,delt_t0,           
     1        delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc,    
     2        delt_sf, beta, pk_n,penaltyK_shear,SIGM_NMAX,FMU_i,C, 
c     4        delt_mf,delt_m0,delt_m,SDEG_NEW , i,                     !调试用 
     3        SDEG)  
            case (2)
              call  BK_Quads(delt_n0,delt_s0,delt_t0,          
     1        delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,    
     2        delt_sf,delt_tf, delt_nf,beta,pk_n,penaltyK_shear,
     3        SDEG)
            end select
C           
          else                                                !若未损伤，SDEG=0
            SDEG=0.d0   
          end if    
C       ---------------------若损伤，根据混合模式计算SDEG值              
        else if( DMICRT >= 1.d0 .and.SDEG < 1D0) then
          DMICRT = 1.d0
          select case (D_MODE) 
            case (0)
              call  BK_M_C(delt_n0,delt_s0,delt_t0,            
     1        delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc,    
     2        delt_sf, beta, pk_n,penaltyK_shear,SIGM_NMAX,FMU_i,C,
c     4        delt_mf,delt_m0,delt_m, SDEG_NEW ,i,                     !调试用               
     3        SDEG)  
            case (2)
              call  BK_Quads(delt_n0,delt_s0,delt_t0,            
     1        delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,    
     2        delt_sf,delt_tf, delt_nf,beta,pk_n,penaltyK_shear,
     3        SDEG)
            end select
          !write(6,*)"DMICRT",DMICRT
        end if
C------------------------------------------------------------------------------------------
        !计算初始断裂时的断裂类型，仅在初始破坏后的下一增量步计算，更新一次后不再计算
        if (DMICRT == 1) then
           if ( stateOld(i,12) == -1.d0) then
             if ( delt_n <= 0  ) then
                 fracture_mode_initial = 90.d0
             else if( delt_n > 0  ) then
                 fracture_mode_initial = atan(delt_shear/delt_n)/pi*180
             end if
           else
               fracture_mode_initial = stateOld(i,12)
           end if
      end if  
C-------------------------------------------------------------
        if (SDEG >= 1D0 )then  !如果完全损伤，则执行摩擦
          
          stateNew(i,1)= 0     !SDV1，即删除标志
          delete_flag = 1       !SDV10
C-------------------------------------------------------------          
          !计算断裂时的断裂类型，仅在完全破坏后的下一增量步计算，更新一次后不再计算
          if ( stateOld(i,11) == -1.d0) then
              if    ( delt_n <= 0  )then 
                  fracture_mode = 90.d0
              else if( delt_n > 0  )then    
                  fracture_mode = atan(delt_shear/delt_n)/pi*180
              end if
          else
              fracture_mode = stateOld(i,11)
          end if
C          write(6,*),'stateOld',stateOld(i,11)
C---------------------------------------------------------------          
          if(delt_n < 0 )then    
              if (delt_shear >= crit_lens ) then 
                  sigm_n = 0.d0        
                  sigm_s = 0.d0
                  sigm_t = 0.d0 
                  delete_flag = 0
              else
                  call Friction_COH (delt_n, strainInc(i,2), 
     1            strainInc(i,3),
     2            pk_n, pk_s, FMU_c, charLength(i),
     3            stressOld( i,2), stressOld( i,3),
     4            sigm_n, sigm_s, sigm_t )
     5            !sigm_shear,delt_E_shear,delt_E_t,sign,delt_E_s,test1)
              end if
          else!delt_n >= 0
              if (delt_n >= crit_lenn .OR. delt_shear >= crit_lens) then
                  delete_flag = 0
              end if
              sigm_n = 0.d0        
              sigm_s = 0.d0
              sigm_t = 0.d0
          end if
          
        else if ( SDEG<1.0d0) then

          stateNew(i,1)= 1
          delete_flag = 1
C--------------------------------------------------------------------        
      !计算SDEG>0.99时的断裂类型。
          if (SDEG > Sdeg_L) then
           if ( stateOld(i,14) == -1.d0) then
             if ( delt_n <= 0  ) then
                 fracture_mode_SDEG = 90.d0
             else if( delt_n > 0  ) then
                 fracture_mode_SDEG = atan(delt_shear/delt_n)/pi*180
             end if
           else
               fracture_mode_SDEG = stateOld(i,14)
           end if
          end if
          
         if (SDEG > Sdeg_L) then
            delete_flag099 = 0
         else if (SDEG <= Sdeg_L) then
            delete_flag099 = 1
         end if   
C      write(6,*),'delete_flag099',delete_flag099   !调试用
C--------------------------------------------------------------------
       if (delt_n  >= 0.d0 ) then                  !更新应力 对损伤后的卸载有效
            pk_nd = (1.d0-SDEG) * pk_n
            pc    = 0.d0
          else if (delt_n  < 0.d0 ) then 
            pk_nd = pk_n
            delt_shear_D = (delt_s0-delt_sc) * (delt_sf-delt_sc) / 
     1                   ( (delt_sf-delt_sc) - SDEG*(delt_sf-delt_s0) )
     2                   + delt_sc
            pc    = SDEG * delt_sc / delt_shear_D
          end if
C       
          pk_sd = (1.d0-SDEG+pc)*pk_s
          !---------------------------------------
          !调整压力变化阶段切应力的变化
          delt_shearInc   = delt_shear - stateOld(i,8)                !切向总位移增量
          if (delt_shearInc<=0 .and. delt_n<0 .and.
     1        strainInc(i,1)<=0) then
              pk_sd = stateOld(i,6)
          endif

c           !---------------------------------------
          sigm_n      = pk_nd * delt_n
          sigm_s      = pk_sd * delt_s
          sigm_t      = pk_sd * delt_t
          ! write(6,*) "delt_sc",delt_sc
c           ----------------------------------------
          
        end if
        if ( delt_n <= 0  ) then
                fracture_b = 90.d0
            elseif( delt_n > 0  ) then
                fracture_b = atan(delt_shear/delt_n)/pi*180
        endif 
C       
C------------------------------------------------------------
C         
        sigm_shear  = SQRT(sigm_s**2 + sigm_t**2) !断裂后等于摩擦力
        stateNew(i,2) = epsil_n 
        stateNew(i,3) = epsil_s
        stateNew(i,4) = epsil_t
        stateNew(i,5) = DMICRT
        stateNew(i,6) = pk_sd
        stateNew(i,7) = SDEG
        stateNew(i,8) = delt_shear 
        stateNew(i,9) = sigm_shear     
        stateNew(i,10) = delete_flag    
        stateNew(i,11) = fracture_mode  !  ↑↑↑↑↑↑  该行及以上不能修改!            
        stateNew(i,12) = fracture_mode_initial  ! 初始断裂损伤类型判断
        stateNew(i,13) = fracture_b             ! 任意时刻损伤类型判断
        stateNew(i,14) = fracture_mode_SDEG     !SDEG大于0.99的内聚力单元的损伤类型
        stateNew(i,15) = delete_flag099         !用于计算SDEG大于0.99的裂纹面积，0—表示SDEG大于0.99，1—表示SDEG小于等于0.99
C     
        stressNew(i,1) = sigm_n
        stressNew(i,2) = sigm_s
        stressNew(i,3) = sigm_t
C        
C------------------------------------------------------------
C
C--------------------更新内能和非弹性耗散能-------------------
C--更新内能---
        stressPower = 0.5d0*(
     1        ( stressNew(i,1)+stressOld(i,1) )*strainInc( i,1) +
     2   2.D0*( stressNew(i,2)+stressOld(i,2) )*strainInc( i,2) +
     3   2.D0*( stressNew(i,3)+stressOld(i,3) )*strainInc( i,3) ) 
C         
        enerInternNew(i) = enerInternOld(i) + stressPower / density(i)
c        
  100 continue
      return
      end 
C------------------------------------------------------------

C------------------------------------------------------------      
           !线性BK_M_C拉伸截断混合模式子程序（计算sdeg值）
C
      subroutine BK_M_C(delt_n0,delt_s0,delt_t0,              !三个方向的初始断裂位移
     1 delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc,  !位移、罚刚度、断裂能及系数
     2 delt_sf, beta, pk_n,penaltyK_shear,SIGM_NMAX,FMU_i,C,
c     4 delt_mf,delt_m0,delt_m,SDEG_NEW,i,                           !调试用
C    WRITE 
     3 SDEG)                                                  !最大分离量、sdeg
C
        real(kind=8) ::  delt_n0,delt_s0,delt_t0,              !定义变量为双精度浮点数   
     1 delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,delt_sc,      
     2 delt_sf,beta, pk_n,penaltyK_shear,SIGM_NMAX,FMU_i,C,
     3 SDEG,SDEG_NEW,delt_m,delt_m0, mix_rate,delt_mf,beta_0
C
        delt_m     = SQRT((DIM(delt_n,0.d0))**2+delt_shear**2)  !当前分离量 
C      
      if (delt_n>0.d0) then                                  !初始断裂分离量
        beta    = ABS(delt_shear/delt_n )                    !混合比 beta>=0
        beta_0  = pk_n*(C-SIGM_NMAX*FMU_i)/                    !计算delt_m0的分界点
     1  (penaltyK_shear*SIGM_NMAX)
C--------------------------        
        if ( beta <=  beta_0 ) then
          delt_m0 = SIGM_NMAX* SQRT(1+beta**2)/pk_n
        else
          delt_m0 = C*SQRT(1+beta**2)/(penaltyK_shear*beta+
     1               pk_n*FMU_i)
        end if
C--------------------------
        mix_rate= penaltyK_shear*beta**2/
     1             (pk_n+penaltyK_shear*beta**2)
        delt_mf = 2.d0/penaltyK_all/delt_m0*
     1             (GIC+(GIIC-GIC)*(mix_rate)**ETA)
      else                                                  !不然则为压剪模式
        delt_m0 = delt_s0
        delt_mf = delt_sf
       !write(6,*),"delt_mf",delt_mf,"delt_m0",delt_m0
      end if
C
      if (delt_m > delt_m0) then
          SDEG_NEW =( (delt_mf-delt_sc )*(delt_m-delt_m0) )
     1              /( (delt_mf-delt_m0)*(delt_m-delt_sc) )    
      else
          SDEG_NEW =0.d0                                    !取零，不影响实际sdeg取值
      end if   
c          
      !  if (i==75) then
      !     write(6,*),"delt_mf",delt_mf,"delt_m0",delt_m0 
      !2     "SDEG_NEW",SDEG_NEW
       ! else
       !end if
        SDEG  = MAX( SDEG ,SDEG_NEW ) 
      return 
      end subroutine BK_M_C     
C
C------------------------------------------------------------
C      
           !线性BK_Quads能量混合模式子程序        
      subroutine BK_Quads(delt_n0,delt_s0,delt_t0,       !三个方向的初始断裂位移
     1 delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,      !位移、罚刚度、断裂能及系数
     2 delt_sf,delt_tf, delt_nf, beta, pk_n,penaltyK_shear,   !各向完全破坏位移
C    WRITE 
     3 SDEG)                                                  !最大分离量、sdeg
C
       real(kind=8) ::  delt_n0,delt_s0,delt_t0,          !定义变量为双精度浮点数   
     1 delt_n, delt_shear,penaltyK_all,GIC,GIIC,ETA,      
     2 delt_sf,delt_tf, delt_nf, beta, pk_n,penaltyK_shear,
     3 SDEG,SDEG_NEW,delt_m,delt_m0, mix_rate,delt_mf
C       
        delt_m     = SQRT((DIM(delt_n,0.d0))**2+delt_shear**2) !当前分离量
C      
        if (delt_n>0.d0) then                                  !初始断裂分离量
          beta    = delt_shear/delt_n                          !混合比
          delt_m0 = delt_n0*delt_s0*sqrt((1+beta**2)/(delt_s0**2+
     1             (beta*delt_n0)**2))
          delt_mf = 2/penaltyK_all/delt_m0*
     1             (GIC+(GIIC-GIC)*(penaltyK_shear*beta**2/
     2             (pk_n+penaltyK_shear*beta**2))**ETA)
        else 
          delt_m0 = delt_s0
          delt_mf = delt_sf
        end if
C
          SDEG_NEW =(delt_mf*(delt_m-delt_m0))
     1    /((delt_mf-delt_m0)*delt_m)
          SDEG  = MAX(SDEG ,SDEG_NEW ) 
       !write(6,*) "sdeg",sdeg
      return
      end subroutine BK_Quads
C
      !friction model
      !参考 ls_dyna 理论手册
      subroutine Friction_COH (delt_n, strainInc_s, strainInc_t, 
     2           pk_n, pk_s, FMU_c, charLength,
     3           sigm_s_old,sigm_t_old,
C    WRITE 
     4           sigm_n, sigm_s, sigm_t
     5                )
          include 'vaba_param.inc'
          dimension delt_E(2),friction_trial(2),          !定义向量
     1              friction_old(2),friction_new(2)
          
          sigm_n          = pk_n  * delt_n
          friction_max    = abs( FMU_c * sigm_n )
          !write(*,*) "------------"
          !write(*,*) "sigm_n" ,sigm_n,friction_max
      
          friction_old = (/ sigm_s_old , sigm_t_old /)    !f_n
          
          delt_E_s = 2 * strainInc_s * charLength         !切向位移增量
          delt_E_t = 2 * strainInc_t * charLength
          delt_E   = (/ delt_E_s , delt_E_t /)
          
          friction_trial =  friction_old + pk_s * delt_E  !f^*
          friction_trial_norm = norm2 ( friction_trial )
          
          if ( friction_trial_norm <= friction_max )then !norm2--2范数
              friction_new = friction_trial
          else
              friction_new = friction_max * 
     1                       friction_trial / friction_trial_norm
          end if
          sigm_s = friction_new(1)
          sigm_t = friction_new(2)
          return
      end subroutine Friction_COH
      
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
c               %%%%%      VUSDFLD for solid element    %%%%%
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
C*******************************************    状态变量数：5   **************************
C*******************************************    删除标志号：1   **************************
      subroutine vusdfld(
c Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElem, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
c Write only -
     *   stateNew, field )
c
      include 'vaba_param.inc'
c
      dimension jElem(nblock), coordMp(nblock,*), 
     *          direct(nblock,3,3), T(nblock,3,3), 
     *          charLength(nblock), props(nprops), 
     *          stateOld(nblock,nstatev), 
     *          stateNew(nblock,nstatev),
     *          field(nblock,nfieldv)
      character*80 cmname
      logical flag_c
c
c     Local arrays from vgetvrm are dimensioned to 
c     maximum block size (maxblk)
c      parameter( nrData=6 )
c      character*3 cData(maxblk*nrData)
c      dimension rData(maxblk*nrData), jData(maxblk*nrData)
c      
      parameter( crit_disp = 5.D0  )
      
      do 100 k = 1, nblock
          
c      没有使用 field ，将其设为固定值   
          field(k,1) = 1.0
          
C      获取第一步的积分点坐标，作为初始坐标
          X = coordMp(k,1)
          Y = coordMp(k,2)
          Z = coordMp(k,3)
          
        flag_c = totalTime <= dt
        
        !write (6,*) flag_c,stepTime,dt
        !write (6,*) "x",coordMp(nblock,1)
        
        if  (flag_c) then
          stateNew(k,2) = X 
          stateNew(k,3) = Y
          stateNew(k,4) = Z
        else
          stateNew(k,2) = stateOld(k,2)
          stateNew(k,3) = stateOld(k,3)
          stateNew(k,4) = stateOld(k,4)  
        end if
        
c      计算积分点相对初始位置的位移  
        d_x = ( X - stateNew(k,2) )**2
        d_y = ( Y - stateNew(k,3) )**2
        d_z = ( Z - stateNew(k,4) )**2
        distance = SQRT ( d_x + d_y + d_z )
        stateNew(k,5) = distance 
        
        if (distance >= crit_disp) then
          stateNew(k,1) = 0
        else
          stateNew(k,1) = 1
        end if
      
c        CALL XPLB_EXIT
  100 continue
c
      return
      end subroutine vusdfld 
