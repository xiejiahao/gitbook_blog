# VASP OPTICS

* `PROGRAM VAMP`

  * `IF (IO%LOPTICS) CALL SET_NABIJ_AUG(P,T_INFO%NTYP)`

  * `IF (IO%LOPTICS) THEN CALL LR_OPTIC`\(linear\_optics.F\)

    * LR\_OPTIC开头这一段很迷，设置了EMAX和EMAX\_CON

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

    * `NBANDS_CDER=MIN(LAST_FILLED_OPTICS(W)*2, WDES%NB_TOT)`设置NBANDS\_CDER，即此子程序中考虑的最大能带数量

    * `CALL LRF_EPSILON`计算介电常数虚部的函数

      * 设置NB\_TOT和NBANDS
        ```fortran
        WDES1=WDES
        WDES1%NB_TOT=NBANDS_CDER
        WDES1%NBANDS=NBANDS_CDER/WDES1%NB_PA
        ```
      * 调用真正的计算函数

        ```fortran
        DO IDIR=1,3
        DO JDIR=1,3
          IF (KPOINTS%ISMEAR >= -1) THEN
            CALL EPSILON_IMAG( WDES, W0, CHAM(:,:,:,:,IDIR), CHAM(:,:,:,:,JDIR), & 
            ENERGY_DER(:,:,:,IDIR), ENERGY_DER(:,:,:,JDIR), EFERMI, &
            NEDOS, EPS_IMAG(:,IDIR,JDIR), DELTAE, KPOINTS%ISMEAR, KPOINTS%SIGMA, & 
            LATT_CUR%OMEGA, WPLASMON(IDIR, JDIR), CON(IDIR, JDIR), RTIME)

          ELSEIF (KPOINTS%ISMEAR <= -4) THEN
            CALL EPSILON_IMAG_TET( WDES, W0, CHAM(:,:,:,:,IDIR), CHAM(:,:,:,:,JDIR), EMAX, &
            NEDOS, EPS_IMAG(:,IDIR,JDIR), DELTAE, LATT_CUR%OMEGA, IO, INFO, KPOINTS)
          ENDIF
        ENDDO
        ENDD
        ```

        其中值得注意的是根据ISMEAR的不同有两种计算方法，-1,0是一种，-4,-5是另一种。  
        虽然在注释中提到这两种算法在理论上是等价的，但是实际上结果确实不尽相同。

        * ISMEAR = -1,0的方法

          ```fortran
              SUBROUTINE EPSILON_IMAG( WDES, W0, CHAM1, CHAM2, DER1, DER2, EFERMI, & 
                NEDOS, DOS, DELTAE, ISMEAR, SIGMA, OMEGA, WPLASMON, CON, TAU)
                DOS=0
                WPLASMON=0
                CON=0

                DO ISP=1,WDES%ISPIN
                DO NK=1,WDES%NKPTS

                   IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

          ! N1 could be restricted to occupied bands
                   DO N1=1,WDES%NB_TOT
          ! N1 must be smaller tan NBANDS_CDER
                   IF (N1>NBANDS_CDER) CYCLE

                   ENERGY=W0%CELTOT(N1,NK,ISP)-EFERMI
                   IF (ISMEAR==-1) THEN
                      SFUN=F(-ENERGY,ABS(SIGMA))
                      DFUN=G(ENERGY,ABS(SIGMA))
                   ELSE
                      CALL DELSTP(ISMEAR, ENERGY/ABS(SIGMA), DFUN, SFUN)
                      DFUN=DFUN/ABS(SIGMA)
                   ENDIF
          ! DFUN d f / d epsilon (= delta)
          ! conversion factor 4 pi e^2/ volume to obtain polarization
                   WPLASMON=WPLASMON+DFUN*DER1(N1,NK,ISP)*DER2(N1,NK,ISP)*WDES%RSPIN*WDES%WTKPT(NK)*(4*PI*FELECT/OMEGA)
                   CON=CON          +DFUN*DER1(N1,NK,ISP)*DER2(N1,NK,ISP)*WDES%RSPIN*WDES%WTKPT(NK)*(4*PI*FELECT/OMEGA)*TAU*CONTCON 

          ! N2 could be restricted to empty bands
                   DO N2=N1+1,WDES%NB_TOT
                      WEIGHT=(W0%FERTOT(N1,NK,ISP)-W0%FERTOT(N2,NK,ISP))*WDES%RSPIN*WDES%WTKPT(NK)
                      DECEL =(W0%CELTOT(N2,NK,ISP)-W0%CELTOT(N1,NK,ISP))
          ! contribution <u_m| - i d/dk_i  |u_n> <u_m| i d/dk_j  |u_n>
          !            = <u_m| d/dk_i |u_n> <u_m| d/dk_j |u_n>
          ! note that we use the transposed and conjugated elements since the
          ! second index is constrained to occupied bands
                      A     =CONJG(CHAM1(N2,N1,NK,ISP))*CHAM2(N2,N1,NK,ISP)
          ! conversion factor 4 pi^2 e^2/ volume to obtain polarization
          ! add contribution for negative frequencies as well
                      CALL SLOT( DECEL, ISMEAR, SIGMA, NEDOS, DELTAE,  WEIGHT*A*4*PI*PI*FELECT/OMEGA, DOS)
                      CALL SLOT(-DECEL, ISMEAR, SIGMA, NEDOS, DELTAE, -WEIGHT*A*4*PI*PI*FELECT/OMEGA, DOS)
                   ENDDO
                   ENDDO
                ENDDO
                ENDDO

              END SUBROUTINE EPSILON_IMA
          ```

          * IF \(NPAR ==1 .AND. KPAR==1 .AND. LPAW\) THEN CALL CALC\_NABIJ\(optics.F\)



