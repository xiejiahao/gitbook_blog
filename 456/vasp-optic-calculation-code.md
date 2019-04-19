# VASP OPTICS

* `PROGRAM VAMP`
    * `IF (IO%LOPTICS) CALL SET_NABIJ_AUG(P,T_INFO%NTYP)`

    * `IF (IO%LOPTICS) THEN CALL LR_OPTIC`(linear_optics.F)
        * LR_OPTIC开头这一段很迷，设置了EMAX和EMAX_CON   
        
          ```fortran
          EMAX=MAX_ENERGY_UNOCCUPIED(WDES,W)*1.2
          IF (EMAX<=0) THEN
          RETURN
          ENDIF
          IF (OMEGAMAX_OPTIC/=-1) THEN
          EMAX=OMEGAMAX_OPTIC
          ENDIF
         
          EMAX_COND=5
          IF (OMEGAMAX_OPTIC/=-1) THEN
          EMAX_COND=OMEGAMAX_OPTIC
          ENDIF
          ```
        * `NBANDS_CDER=MIN(LAST_FILLED_OPTICS(W)*2, WDES%NB_TOT)`设置NBANDS_CDER，即此子程序中考虑的最大能带数量
        * `CALL LRF_EPSILON`计算介电常数虚部的函数
            * 设置NB_TOT和NBANDS
            ```fortran
            WDES1=WDES
            WDES1%NB_TOT=NBANDS_CDER
            WDES1%NBANDS=NBANDS_CDER/WDES1%NB_PA
            ```
            * 调用真正的计算函数
            
            ```fortran
    * IF (NPAR ==1 .AND. KPAR==1 .AND. LPAW) THEN CALL CALC_NABIJ(optics.F)