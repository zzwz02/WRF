MODULE module_cellpert

! AUTHOR: Domingo Munoz-Esparza (NCAR/RAL)
! Cell perturbation method for turbulence generation in nested WRF-LES

USE module_domain_type
USE module_configure
#ifdef DM_PARALLEL
   USE module_dm
   USE module_comm_dm, ONLY : &
                           HALO_EM_PXY_sub   &
                          ,PERIOD_EM_PXY_sub
#endif

CONTAINS

SUBROUTINE cp_perturb(th_pert,curr_secs,phb,ph_2,              &
                      config_flags,grid,mytask,                &
                      its,ite,jts,jte,kts,kte,                 &
                      ids,ide,jds,jde,kds,kde,                 &
                      ips,ipe,jps,jpe,kps,kpe,                 &
                      ims,ime,jms,jme,kms,kme)

   ! Generates 3d-array of random moist potential temperature perturbations
   ! for the cell perturbation method

   IMPLICIT NONE

   TYPE(domain), INTENT(IN) :: grid
   TYPE (grid_config_rec_type), INTENT (IN) :: config_flags
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th_pert
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: phb,ph_2
   INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte, &
                          ids,ide,jds,jde,kds,kde, &
                          ips,ipe,jps,jpe,kps,kpe, &
                          ims,ime,jms,jme,kms,kme 
   INTEGER, INTENT(IN) :: curr_secs
   INTEGER, INTENT(IN) :: mytask

   ! local variables (start)

   CHARACTER*256 :: message
   INTEGER :: cs, nbcx, nbcy, npx, npy, ncxt, ncyt, opp_x, opp_y, &
              i, j, k, how_many, isc, jsc, indf_x, indf_y 
   INTEGER , ALLOCATABLE , DIMENSION(:) :: iseed
   REAL, DIMENSION(1:kte-1) :: z_vec
   REAL :: cbl_z_fact_2, cbl_z_fact, nczr, nbcz_zi, cbl_uowthst, rand_valtemp

   ! local variables (end)

   cs = config_flags%cell_gppc
   nbcx = config_flags%cell_nbcx
   nbcy = config_flags%cell_nbcy
   npx = ide-1
   npy = jde-1
   ncxt = npx/cs
   ncyt = npy/cs
   opp_x = (ncyt-nbcx)*cs
   opp_y = (ncxt-nbcy)*cs   

   ! -------------------------------------------------- !
   ! first time CPM is executed, set random number seed !
   ! -------------------------------------------------- !
   IF ( curr_secs .EQ. 1 ) THEN
     write(message,*) 'Initializing random seed'
     call wrf_message( trim(message) )
     CALL random_seed(size=how_many)
     IF ( ALLOCATED(iseed)) DEALLOCATE(iseed)
     ALLOCATE(iseed(how_many))
     DO i = 1,how_many
       iseed(i)=mytask
     ENDDO
     call random_seed(put=iseed(1:how_many))
     DEALLOCATE(iseed)
   ENDIF

   ! --------------------------------------  !
   ! calculate vertical hegiht (mass points) !
   ! --------------------------------------- !
   DO k = kts, kte-1
    z_vec(k) = 0.5*((phb(its,k,jts) + ph_2(its,k,jts))+(phb(its,k+1,jts) + ph_2(its,k+1,jts)))/9.81
   ENDDO

   ! --------------------- !
   ! fill in th_pert array !
   ! --------------------- !
   DO k = kts, kte-1

     ! ----------------------------------------------------------------------- !
     ! thermal-variance scaling for CBLs (Munoz-Esparza and Kosovic, MWR 2018) !
     ! ----------------------------------------------------------------------- !
     IF ( config_flags%cell_pert_cbl ) THEN
       cbl_uowthst = config_flags%cell_cbl_uowthst
       nbcz_zi = config_flags%cell_cbl_zitop
       nczr = z_vec(k)
       IF (nczr .lt. nbcz_zi) THEN
         cbl_z_fact_2 = 2.*(nczr/nbcz_zi)**(-2./3.)*(1.-nczr/nbcz_zi)**(4./3.) + 0.94*(nczr/nbcz_zi)**(4./3.)*(1.-nczr/nbcz_zi)**(-2./3.)
         cbl_z_fact = sqrt (cbl_z_fact_2) * cbl_uowthst
       ELSE
         cbl_z_fact = 1.0
       ENDIF
     ELSE
       cbl_z_fact = 1.0
     ENDIF

     DO i = its, ite
       IF ( MOD(i-1, cs) .EQ. 0. ) THEN
         DO j = jts, jte
           IF ( MOD(j-1, cs) .EQ. 0. ) THEN
             call random_number (rand_valtemp)
             DO isc = i, min(i+cs-1,ite)
               DO jsc = j, min(j+cs-1,jte)
                 IF ( (i .GE. its) .AND. (i .LE. ite) .AND. (j .GE. jts) .AND. (j .LE. jte)  ) THEN
                   th_pert(isc,k,jsc) = (( rand_valtemp * 2 ) - 1) * cbl_z_fact
                 ENDIF
               ENDDO
             ENDDO
           ENDIF
         ENDDO
       ENDIF
     ENDDO

   ENDDO ! k = kts, kte

#ifdef DM_PARALLEL
#      include "HALO_EM_PXY.inc"
#      include "PERIOD_EM_PXY.inc"
#endif

   DO k = kts, kte-1

     IF ( ( MOD(its-1, cs) .NE. 0. ) .AND. ( MOD(jts-1, cs) .EQ. 0. ) )THEN
       indf_x = cs - MOD(its-1,cs) -1
       DO isc = its, its + indf_x
         DO jsc = jts, jte
           th_pert(isc,k,jsc) = th_pert(its-1,k,jsc)
         ENDDO
       ENDDO
     ELSEIF ( ( MOD(its-1, cs) .EQ. 0. ) .AND. ( MOD(jts-1, cs) .NE. 0. ) )THEN
       indf_y = cs - MOD(jts-1,cs) -1
       DO jsc = jts, jts + indf_y
         DO isc = its, ite
           th_pert(isc,k,jsc) = th_pert(isc,k,jts-1)
         ENDDO
       ENDDO
     ELSEIF ( ( MOD(its-1, cs) .NE. 0. ) .AND. ( MOD(jts-1, cs) .NE. 0. ) )THEN
       indf_x = cs - MOD(its-1,cs) -1
       indf_y = cs - MOD(jts-1,cs) -1
       DO isc = its, its + indf_x
         DO jsc = jts, jts + indf_y
           th_pert(isc,k,jsc) = th_pert(its-1,k,jts-1)
         ENDDO
       ENDDO
       DO isc = its + indf_x +1 , ite
         DO jsc = jts, jts + indf_y
           th_pert(isc,k,jsc) = th_pert(isc,k,jts-1)
         ENDDO
       ENDDO
       DO jsc = jts + indf_y +1 , jte
         DO isc = its, its + indf_x
           th_pert(isc,k,jsc) = th_pert(its-1,k,jsc)
         ENDDO
       ENDDO
     ENDIF

   ENDDO ! k = kts, kte


END SUBROUTINE cp_perturb

SUBROUTINE cp_apply_4bdy(th_pert,t_2,phb,ph_2,                 &
                      config_flags,grid,                       &
                      its,ite,jts,jte,kts,kte,                 &
                      ids,ide,jds,jde,kds,kde,                 &
                      ips,ipe,jps,jpe,kps,kpe,                 &
                      ims,ime,jms,jme,kms,kme)

   ! Apply theta perturbations according to the cell perturbation method

   IMPLICIT NONE

   TYPE(domain), INTENT(IN) :: grid
   TYPE (grid_config_rec_type), INTENT (IN) :: config_flags
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_2
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: th_pert,phb,ph_2
   INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte, &
                          ids,ide,jds,jde,kds,kde, &
                          ips,ipe,jps,jpe,kps,kpe, &
                          ims,ime,jms,jme,kms,kme

   ! local variables (start)

   CHARACTER*256 :: message
   INTEGER :: cs, nbcx, nbcy, npx, npy, ncxt, ncyt, opp_x, opp_y, &
              i, j, k, nbcz_s, nbcz_e, nbczRL_s, nbczRL_e, &
              ibs, ibe, jbs, jbe, b_dist, b_limit, itf, jtf
   REAL :: nbcz_s_m, nbcz_e_m, nbczRL_s_m, nbczRL_e_m, cpert_amp, &
           cpert_ampRL, z_tmp
   REAL, DIMENSION(1:kte-1) :: z_vec

   ! local variables (end)

   cs = config_flags%cell_gppc
   nbcx = config_flags%cell_nbcx
   nbcy = config_flags%cell_nbcy
   npx = ide-1
   npy = jde-1
   ncxt = npx/cs
   ncyt = npy/cs
   opp_x = (ncyt-nbcx)*cs
   opp_y = (ncxt-nbcy)*cs
   ibs = ids
   ibe = ide-1
   jbs = jds
   jbe = jde-1
   itf = min(ite,ibe)
   jtf = min(jte,jbe)

   IF ( config_flags%cell_tvcp ) THEN
     nbcz_e_m = grid%zpert_tvcp
     cpert_amp = grid%dthpert_tvcp
     IF ( config_flags%cell_pert_RL) THEN
       nbczRL_s_m = grid%zpertRLs_tvcp
       nbczRL_e_m = grid%zpertRLe_tvcp
       cpert_ampRL = grid%dthpertRL_tvcp
     ENDIF
   ELSE
     nbcz_e_m = config_flags%cell_ztop
     IF ( config_flags%cell_pert_cbl ) THEN
       cpert_amp = 1.0
     ELSE
       cpert_amp = config_flags%cell_pert_amp
     ENDIF
     IF ( config_flags%cell_pert_RL) THEN
       nbczRL_s_m = config_flags%cellRL_zbottom
       nbczRL_e_m = config_flags%cellRL_ztop
       cpert_ampRL = config_flags%cellRL_pert_amp
     ENDIF
   ENDIF
   nbcz_s_m = config_flags%cell_zbottom

   ! ------------------------------------------------------- !
   ! determine vertical grid point indices for perturbations !
   ! ------------------------------------------------------- !
   DO k = kts, kte-1
     z_vec(k) = 0.5*((phb(its,k,jts) + ph_2(its,k,jts))+(phb(its,k+1,jts) + ph_2(its,k+1,jts)))/9.81
   ENDDO
   DO k = kts, kte-1
     z_tmp = z_vec(k)
     nbcz_s = k
     IF (z_tmp .ge. nbcz_s_m) EXIT
   ENDDO
   nbcz_s = max(max(nbcz_s,config_flags%cell_kbottom),kts)
   DO k = kts, kte-1
     z_tmp = z_vec(k)
     nbcz_e = k
     IF (z_tmp .ge. nbcz_e_m) EXIT
   ENDDO
   nbcz_e = min(nbcz_e,kte-1)

   IF ( config_flags%cell_pert_RL) THEN

     DO k = kts, kte-1
       z_tmp = z_vec(k)
       nbczRL_s = k
       IF (z_tmp .ge. nbczRL_s_m) EXIT
     ENDDO
     nbczRL_s = max(max(nbcz_e,nbczRL_s),kts)
     DO k = kts, kte-1
       z_tmp = z_vec(k)
       PRINT *,'k,z_tmp=',k,',',z_tmp,'m'
       nbczRL_e = k
       IF (z_tmp .ge. nbczRL_e_m) EXIT
     ENDDO
     nbczRL_e = min(max(nbcz_s,nbczRL_e),kte-1)

   ENDIF

   ! ------------------------- !
   ! print general information !
   ! ------------------------- !
   write(message,*) 'DOING CELL PERT 4 BOUNDARIES'
   call wrf_message( trim(message) )
   write(message,*) 'Lowest vertical grid point perturbed BL = ',nbcz_s 
   call wrf_message( trim(message) )
   write(message,*) 'Highest vertical grid point perturbed BL = ',nbcz_e
   call wrf_message( trim(message) )
   IF (.NOT. config_flags%cell_pert_cbl ) THEN
     write(message,*) 'Perturbation amplitude = ',cpert_amp,' K'
     call wrf_message( trim(message) )
   ELSE
     write(message,*) 'Perturbation amplitude based on thermal variance scaling in CBLs, cell_cbl_uowthst =',config_flags%cell_cbl_uowthst
     call wrf_message( trim(message) )
   ENDIF
   IF ( config_flags%cell_pert_RL) THEN
     IF (nbczRL_s_m .gt. nbczRL_e_m ) THEN
       write(message,*) 'cellRL_zbottom > cellRL_ztop -> Not applying RL perturbations'
       call wrf_message( trim(message) )
     ELSE
       write(message,*) 'Lowest vertical grid point perturbed RL = ',nbczRL_s
       call wrf_message( trim(message) )
       write(message,*) 'Highest vertical grid point perturbed RL = ',nbczRL_e
       call wrf_message( trim(message) )
       write(message,*) 'Perturbation amplitude RL = ',cpert_ampRL,' K'
       call wrf_message( trim(message) )
     ENDIF
   ELSE
     write(message,*) 'Not applying RL perturbations'
     call wrf_message( trim(message) )
   ENDIF

   ! ------------------------------------- !
   ! add cell perturbations to theta field !
   ! ------------------------------------- !
   IF (jts - jbs .lt. cs*nbcx) THEN
     DO j = jts, min(jtf,jbs+cs*nbcx-1)
       b_dist = j - jbs
       b_limit = b_dist
       DO k = nbcz_s, nbcz_e
         DO i = max(its,b_limit+ibs), min(itf,ibe-b_limit)
           t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
         ENDDO
       ENDDO
     ENDDO
   ENDIF
   IF (jbe - jtf .lt. cs*nbcx) THEN
     DO j = max(jts,jbe-cs*nbcx+1), jtf
       b_dist = jbe - j
       b_limit = b_dist
       DO k = nbcz_s, nbcz_e
         DO i = max(its,b_limit+ibs), min(itf,ibe-b_limit)
           t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
         ENDDO
       ENDDO
     ENDDO
   ENDIF
   IF (its - ibs .lt. cs*nbcy) THEN
     DO i = its, min(itf,ibs+cs*nbcy-1)
       b_dist = i - ibs
       DO k = nbcz_s, nbcz_e
         DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
           t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
         ENDDO
       ENDDO
     ENDDO
   ENDIF
   IF (ibe - itf .lt. cs*nbcy) THEN
     DO i = max(its,ibe-cs*nbcy+1), itf
       b_dist = ibe - i
       DO k = nbcz_s, nbcz_e
         DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
           t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
         ENDDO
       ENDDO
     ENDDO
   ENDIF

   IF ( config_flags%cell_pert_RL) THEN

     IF (jts - jbs .lt. cs*nbcx) THEN
       DO j = jts, min(jtf,jbs+cs*nbcx-1)
         b_dist = j - jbs
         b_limit = b_dist
         DO k = nbczRL_s, nbczRL_e
           DO i = max(its,b_limit+ibs), min(itf,ibe-b_limit)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_ampRL
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF (jbe - jtf .lt. cs*nbcx) THEN
       DO j = max(jts,jbe-cs*nbcx+1), jtf
         b_dist = jbe - j
         b_limit = b_dist
         DO k = nbczRL_s, nbczRL_e
           DO i = max(its,b_limit+ibs), min(itf,ibe-b_limit)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_ampRL
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF (its - ibs .lt. cs*nbcy) THEN
       DO i = its, min(itf,ibs+cs*nbcy-1)
         b_dist = i - ibs
         DO k = nbczRL_s, nbczRL_e
           DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_ampRL
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF (ibe - itf .lt. cs*nbcy) THEN
       DO i = max(its,ibe-cs*nbcy+1), itf
         b_dist = ibe - i
         DO k = nbczRL_s, nbczRL_e
           DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_ampRL
           ENDDO
         ENDDO
       ENDDO
     ENDIF

   ENDIF

END SUBROUTINE cp_apply_4bdy

SUBROUTINE cp_apply_2bdy(th_pert,t_2,phb,ph_2,                 &
                      config_flags,grid,                       &
                      its,ite,jts,jte,kts,kte,                 &
                      ids,ide,jds,jde,kds,kde,                 &
                      ips,ipe,jps,jpe,kps,kpe,                 &
                      ims,ime,jms,jme,kms,kme)

   ! Apply theta perturbations according to the cell perturbation method

   IMPLICIT NONE

   TYPE(domain), INTENT(IN) :: grid
   TYPE (grid_config_rec_type), INTENT (IN) :: config_flags
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_2
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: th_pert,phb,ph_2
   INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte, &
                          ids,ide,jds,jde,kds,kde, &
                          ips,ipe,jps,jpe,kps,kpe, &
                          ims,ime,jms,jme,kms,kme

   ! local variables (start)

   CHARACTER*256 :: message
   INTEGER :: cs, nbcx, nbcy, npx, npy, ncxt, ncyt, opp_x, opp_y, &
              i, j, k, nbcz_s, nbcz_e, nbczRL_s, nbczRL_e, &
              ibs, ibe, jbs, jbe, b_dist, b_limit, itf, jtf
   REAL :: nbcz_s_m, nbcz_e_m, nbczRL_s_m, nbczRL_e_m, cpert_amp, &
           cpert_ampRL, z_tmp
   REAL, DIMENSION(1:kte-1) :: z_vec

   ! local variables (end)

   cs = config_flags%cell_gppc
   nbcx = config_flags%cell_nbcx
   nbcy = config_flags%cell_nbcy
   npx = ide-1
   npy = jde-1
   ncxt = npx/cs
   ncyt = npy/cs
   opp_x = (ncyt-nbcx)*cs
   opp_y = (ncxt-nbcy)*cs
   ibs = ids
   ibe = ide-1
   jbs = jds
   jbe = jde-1
   itf = min(ite,ibe)
   jtf = min(jte,jbe)

   IF ( config_flags%cell_tvcp ) THEN
     nbcz_e_m = grid%zpert_tvcp
     cpert_amp = grid%dthpert_tvcp
     IF ( config_flags%cell_pert_RL) THEN
       nbczRL_s_m = grid%zpertRLs_tvcp
       nbczRL_e_m = grid%zpertRLe_tvcp
       cpert_ampRL = grid%dthpertRL_tvcp
     ENDIF
   ELSE
     nbcz_e_m = config_flags%cell_ztop
     IF ( config_flags%cell_pert_cbl ) THEN
       cpert_amp = 1.0
     ELSE
       cpert_amp = config_flags%cell_pert_amp
     ENDIF
     IF ( config_flags%cell_pert_RL) THEN
       nbczRL_s_m = config_flags%cellRL_zbottom
       nbczRL_e_m = config_flags%cellRL_ztop
       cpert_ampRL = config_flags%cellRL_pert_amp
     ENDIF
   ENDIF
   nbcz_s_m = config_flags%cell_zbottom

   ! ------------------------------------------------------- !
   ! determine vertical grid point indices for perturbations !
   ! ------------------------------------------------------- !
   DO k = kts, kte-1
     z_vec(k) = 0.5*((phb(its,k,jts) + ph_2(its,k,jts))+(phb(its,k+1,jts) + ph_2(its,k+1,jts)))/9.81
   ENDDO
   DO k = kts, kte-1
     z_tmp = z_vec(k)
     nbcz_s = k
     IF (z_tmp .ge. nbcz_s_m) EXIT
   ENDDO
   nbcz_s = max(max(nbcz_s,config_flags%cell_kbottom),kts)
   DO k = kts, kte-1
     z_tmp = z_vec(k)
     nbcz_e = k
     IF (z_tmp .ge. nbcz_e_m) EXIT
   ENDDO
   nbcz_e = min(nbcz_e,kte-1)

   IF ( config_flags%cell_pert_RL) THEN

     DO k = kts, kte-1
       z_tmp = z_vec(k)
       nbczRL_s = k
       IF (z_tmp .ge. nbczRL_s_m) EXIT
     ENDDO
     nbczRL_s = max(max(nbcz_e,nbczRL_s),kts)
     DO k = kts, kte-1
       z_tmp = z_vec(k)
       nbczRL_e = k
       IF (z_tmp .ge. nbczRL_e_m) EXIT
     ENDDO
     nbczRL_e = min(max(nbcz_s,nbczRL_e),kte)

   ENDIF

   ! ------------------------- !
   ! print general information !
   ! ------------------------- !
   IF (config_flags%cell_pert_2d_opt .EQ. 0) THEN
     write(message,*) 'DOING CELL PERT TWO BOUNDARIES (SOUTH & WEST)'
     call wrf_message( trim(message) )
   ELSE IF (config_flags%cell_pert_2d_opt .EQ. 1) THEN
     write(message,*) 'DOING CELL PERT TWO BOUNDARIES (NORTH & WEST)'
     call wrf_message( trim(message) )
   ELSE
     write(message,*) 'NOT A SUPPORTED OPTION. SET cell_pert_2d_opt TO 0 OR 1'
     call wrf_message( trim(message) )
   ENDIF

   write(message,*) 'Lowest vertical grid point perturbed BL = ',nbcz_s
   call wrf_message( trim(message) )
   write(message,*) 'Highest vertical grid point perturbed BL = ',nbcz_e
   call wrf_message( trim(message) )
   IF (.NOT. config_flags%cell_pert_cbl ) THEN
     write(message,*) 'Perturbation amplitude = ',cpert_amp,' K'
     call wrf_message( trim(message) ) 
   ELSE
     write(message,*) 'Perturbation amplitude based on thermal variance scaling in CBLs, cell_cbl_uowthst =',config_flags%cell_cbl_uowthst
     call wrf_message( trim(message) ) 
   ENDIF
   IF ( config_flags%cell_pert_RL) THEN
     IF (nbczRL_s_m .gt. nbczRL_e_m ) THEN
       write(message,*) 'cellRL_zbottom > cellRL_ztop -> Not applying RL perturbations'
       call wrf_message( trim(message) ) 
     ELSE
       write(message,*) 'Lowest vertical grid point perturbed RL = ',nbczRL_s
       call wrf_message( trim(message) ) 
       write(message,*) 'Highest vertical grid point perturbed RL = ',nbczRL_e
       call wrf_message( trim(message) ) 
       write(message,*) 'Perturbation amplitude RL = ',cpert_ampRL,' K'
       call wrf_message( trim(message) ) 
     ENDIF
   ELSE
     write(message,*) 'Not applying RL perturbations'
     call wrf_message( trim(message) ) 
   ENDIF

   ! ------------------------------------- !
   ! add cell perturbations to theta field !
   ! ------------------------------------- !
   IF (config_flags%cell_pert_2d_opt .EQ. 0 ) THEN ! south and west boundaries

     IF (jts - jbs .lt. cs*nbcx) THEN
       DO j = jts, min(jtf,jbs+cs*nbcx-1)
         b_dist = j - jbs
         b_limit = b_dist
         DO k = nbcz_s, nbcz_e
           DO i = max(its,b_limit+ibs), min(itf,ibe)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF (its - ibs .lt. cs*nbcy) THEN
       DO i = its, min(itf,ibs+cs*nbcy-1)
         b_dist = i - ibs
         DO k = nbcz_s, nbcz_e
           DO j = max(jts,b_dist+jbs+1), min(jtf,jbe)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
           ENDDO
         ENDDO
       ENDDO
     ENDIF     

     IF ( config_flags%cell_pert_RL) THEN
        IF (jts - jbs .lt. cs*nbcx) THEN
          DO j = jts, min(jtf,jbs+cs*nbcx-1)
            b_dist = j - jbs
            b_limit = b_dist
            DO k = nbczRL_s, nbczRL_e
              DO i = max(its,b_limit+ibs), min(itf,ibe)
                t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_ampRL
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (its - ibs .lt. cs*nbcy) THEN
          DO i = its, min(itf,ibs+cs*nbcy-1)
            b_dist = i - ibs
            DO k = nbczRL_s, nbczRL_e
              DO j = max(jts,b_dist+jbs+1), min(jtf,jbe)
                t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_ampRL
              ENDDO
            ENDDO
          ENDDO
        ENDIF
     ENDIF

   ELSEIF (config_flags%cell_pert_2d_opt .EQ. 1 ) THEN ! north and west boundaries

     IF (its - ibs .lt. cs*nbcy) THEN
       DO i = its, min(itf,ibs+cs*nbcy-1)
         b_dist = i - ibs
         DO k = nbcz_s, nbcz_e
           DO j = max(jts,jbs), min(jtf,jbe)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF (jbe - jtf .lt. cs*nbcx) THEN
       DO j = max(jts,jbe-cs*nbcx+1), jtf
         b_dist = jbe - j
         b_limit = b_dist
         DO k = nbcz_s, nbcz_e
           DO i = max(its,cs*nbcx+1), min(itf,ibe)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
           ENDDO
         ENDDO
       ENDDO
     ENDIF

     IF ( config_flags%cell_pert_RL) THEN
       IF (its - ibs .lt. cs*nbcy) THEN
         DO i = its, min(itf,ibs+cs*nbcy-1)
           b_dist = i - ibs
           DO k = nbczRL_s, nbczRL_e
             DO j = max(jts,jbs), min(jtf,jbe)
               t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
             ENDDO
           ENDDO
         ENDDO
       ENDIF
       IF (jbe - jtf .lt. cs*nbcx) THEN
         DO j = max(jts,jbe-cs*nbcx+1), jtf
           b_dist = jbe - j
           b_limit = b_dist
           DO k = nbczRL_s, nbczRL_e
             DO i = max(its,cs*nbcx+1), min(itf,ibe)
               t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
             ENDDO
           ENDDO
         ENDDO
       ENDIF
     ENDIF

   ENDIF


END SUBROUTINE cp_apply_2bdy

SUBROUTINE cp_apply_1bdy(th_pert,t_2,phb,ph_2,                 &
                      config_flags,grid,                       &
                      its,ite,jts,jte,kts,kte,                 &
                      ids,ide,jds,jde,kds,kde,                 &
                      ips,ipe,jps,jpe,kps,kpe,                 &
                      ims,ime,jms,jme,kms,kme)

   ! Apply theta perturbations according to the cell perturbation method

   IMPLICIT NONE

   TYPE(domain), INTENT(IN) :: grid
   TYPE (grid_config_rec_type), INTENT (IN) :: config_flags
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_2
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: th_pert,phb,ph_2
   INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte, &
                          ids,ide,jds,jde,kds,kde, &
                          ips,ipe,jps,jpe,kps,kpe, &
                          ims,ime,jms,jme,kms,kme

   ! local variables (start)

   CHARACTER*256 :: message
   INTEGER :: cs, nbcx, nbcy, npx, npy, ncxt, ncyt, opp_x, opp_y, &
              i, j, k, nbcz_s, nbcz_e, nbczRL_s, nbczRL_e, &
              ibs, ibe, jbs, jbe, b_dist, b_limit, itf, jtf
   REAL :: nbcz_s_m, nbcz_e_m, nbczRL_s_m, nbczRL_e_m, cpert_amp, &
           cpert_ampRL, z_tmp
   REAL, DIMENSION(1:kte-1) :: z_vec

   ! local variables (end)

   cs = config_flags%cell_gppc
   nbcx = config_flags%cell_nbcx
   nbcy = config_flags%cell_nbcy
   npx = ide-1
   npy = jde-1
   ncxt = npx/cs
   ncyt = npy/cs
   opp_x = (ncyt-nbcx)*cs
   opp_y = (ncxt-nbcy)*cs
   ibs = ids
   ibe = ide-1
   jbs = jds
   jbe = jde-1
   itf = min(ite,ibe)
   jtf = min(jte,jbe)

   IF ( config_flags%cell_tvcp ) THEN
     nbcz_e_m = grid%zpert_tvcp
     cpert_amp = grid%dthpert_tvcp
     IF ( config_flags%cell_pert_RL) THEN
       nbczRL_s_m = grid%zpertRLs_tvcp
       nbczRL_e_m = grid%zpertRLe_tvcp
       cpert_ampRL = grid%dthpertRL_tvcp
     ENDIF
   ELSE
     nbcz_e_m = config_flags%cell_ztop
     IF ( config_flags%cell_pert_cbl ) THEN
       cpert_amp = 1.0
     ELSE
       cpert_amp = config_flags%cell_pert_amp
     ENDIF
     IF ( config_flags%cell_pert_RL) THEN
       nbczRL_s_m = config_flags%cellRL_zbottom
       nbczRL_e_m = config_flags%cellRL_ztop
       cpert_ampRL = config_flags%cellRL_pert_amp
     ENDIF
   ENDIF
   nbcz_s_m = config_flags%cell_zbottom

   ! ------------------------------------------------------- !
   ! determine vertical grid point indices for perturbations !
   ! ------------------------------------------------------- !
   DO k = kts, kte-1
     z_vec(k) = 0.5*((phb(its,k,jts) + ph_2(its,k,jts))+(phb(its,k+1,jts) + ph_2(its,k+1,jts)))/9.81
   ENDDO
   DO k = kts, kte-1
     z_tmp = z_vec(k)
     nbcz_s = k
     IF (z_tmp .ge. nbcz_s_m) EXIT
   ENDDO
   nbcz_s = max(max(nbcz_s,config_flags%cell_kbottom),kts)
   DO k = kts, kte-1
     z_tmp = z_vec(k)
     nbcz_e = k
     IF (z_tmp .ge. nbcz_e_m) EXIT
   ENDDO
   nbcz_e = min(nbcz_e,kte-1)

   ! ------------------------- !
   ! print general information !
   ! ------------------------- !
   write(message,*) 'DOING CELL PERT ONE BOUNDARY (WEST)'
   call wrf_message( trim(message) )
   write(message,*) 'Lowest vertical grid point perturbed BL = ',nbcz_s
   call wrf_message( trim(message) )
   write(message,*) 'Highest vertical grid point perturbed BL = ',nbcz_e
   call wrf_message( trim(message) )
   IF (.NOT. config_flags%cell_pert_cbl ) THEN
     write(message,*) 'Perturbation amplitude = ',cpert_amp,' K'
     call wrf_message( trim(message) )
   ELSE
     write(message,*) 'Perturbation amplitude based on thermal variance scaling in CBLs, cell_cbl_uowthst =',config_flags%cell_cbl_uowthst
     call wrf_message( trim(message) )
   ENDIF
   IF ( config_flags%cell_pert_RL) THEN
     IF (nbczRL_s_m .gt. nbczRL_e_m ) THEN
       write(message,*) 'cellRL_zbottom > cellRL_ztop -> Not applying RL perturbations'
       call wrf_message( trim(message) )
     ELSE
       write(message,*) 'Lowest vertical grid point perturbed RL = ',nbczRL_s
       call wrf_message( trim(message) )
       write(message,*) 'Highest vertical grid point perturbed RL = ',nbczRL_e
       call wrf_message( trim(message) )
       write(message,*) 'Perturbation amplitude RL = ',cpert_ampRL,' K'
       call wrf_message( trim(message) )
     ENDIF
   ELSE
     write(message,*) 'Not applying RL perturbations'
     call wrf_message( trim(message) )
   ENDIF

   ! ------------------------------------- !
   ! add cell perturbations to theta field !
   ! ------------------------------------- !
   IF (its - ibs .lt. cs*nbcy) THEN
     DO i = its, min(itf,ibs+cs*nbcy-1)
       b_dist = i - ibs
       DO k = nbcz_s, nbcz_e
         DO j = max(jts,jbs), min(jtf,jbe)
           t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
         ENDDO
       ENDDO
     ENDDO
   ENDIF

   IF ( config_flags%cell_pert_RL) THEN

     IF (its - ibs .lt. cs*nbcy) THEN
       DO i = its, min(itf,ibs+cs*nbcy-1)
         b_dist = i - ibs
         DO k = nbczRL_s, nbczRL_e
           DO j = max(jts,jbs), min(jtf,jbe)
             t_2(i,k,j) = t_2(i,k,j) + th_pert(i,k,j) * cpert_amp
           ENDDO
         ENDDO
       ENDDO
     ENDIF

   ENDIF

END SUBROUTINE cp_apply_1bdy

END MODULE module_cellpert
