
!WRF:MODEL_LAYER:BOUNDARY
!
MODULE module_bc_em

   USE module_bc, ONLY: set_physical_bc2d, set_physical_bc3d, spec_bdytend, &
                        spec_bdytend_perturb, relax_bdytend_tile, relax_bdytend, &
                        spec_bdytend_perturb_chem
   USE module_configure, ONLY: grid_config_rec_type
   USE module_wrf_error
   USE module_model_constants, ONLY: R_d, R_v, T0

! --------- LJMB ---------
! For intra-processor communications:
#ifdef DM_PARALLEL
   USE module_dm
   USE module_comm_dm, ONLY : &
                           HALO_EM_MOM_PERT_U_sub   &
                          ,PERIOD_EM_MOM_PERT_U_sub &
                          ,HALO_EM_MOM_PERT_V_sub   &
                          ,PERIOD_EM_MOM_PERT_V_sub &
                          ,HALO_EM_MOM_PERT_W_sub   &
                          ,PERIOD_EM_MOM_PERT_W_sub
#endif

   USE module_domain_type

   ! In WRFV3, the module_domain_type is defined
   ! in a separaate source file, frame/module_domain_type.F
   ! This enables splitting off the alloc_space_field routine
   ! into a separate file, reducing the size of module_domain

   !  Now that a "domain" TYPE exists, we can use it to store a few pointers
   !  to this type.  These are primarily for use in traversing the linked list.
   !  The "head_grid" is always the pointer to the first domain that is
   !  allocated.  This is available and is not to be changed.  The others are
   !  just temporary pointers.

! ------- END LJMB --------


CONTAINS

!------------------------------------------------------------------------

   SUBROUTINE spec_bdyupdate_ph( ph_save, field,           &
                               field_tend, mu_tend, muts,  &
                               c1, c2, dt,                 &
                               variable_in, config_flags,  &
                               spec_zone,                  &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

!  This subroutine adds the tendencies in the boundary specified region.
!  spec_zone is the width of the outer specified b.c.s that are set here.
!  (JD August 2000)

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      INTEGER,      INTENT(IN   )    :: spec_zone
      CHARACTER,    INTENT(IN   )    :: variable_in
      REAL,         INTENT(IN   )    :: dt


      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(INOUT) :: field
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: field_tend, ph_save
      REAL,  DIMENSION( ims:ime , jms:jme ), INTENT(IN   ) :: mu_tend, muts
      REAL,  DIMENSION( kms:kme ), INTENT(IN   ) :: c1, c2
      TYPE( grid_config_rec_type ) config_flags

      CHARACTER  :: variable
      INTEGER    :: i, j, k, ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: b_dist, b_limit

!     Local array

      REAL,  DIMENSION( its:ite , jts:jte ) :: mu_old
      LOGICAL    :: periodic_x

      periodic_x = config_flags%periodic_x

      variable = variable_in

      IF (variable == 'U') variable = 'u'
      IF (variable == 'V') variable = 'v'
      IF (variable == 'M') variable = 'm'
      IF (variable == 'H') variable = 'h'

      ibs = ids
      ibe = ide-1
      itf = min(ite,ide-1)
      jbs = jds
      jbe = jde-1
      jtf = min(jte,jde-1)
      ktf = kde-1
      IF (variable == 'u') ibe = ide
      IF (variable == 'u') itf = min(ite,ide)
      IF (variable == 'v') jbe = jde
      IF (variable == 'v') jtf = min(jte,jde)
      IF (variable == 'm') ktf = kte
      IF (variable == 'h') ktf = kte

      IF (jts - jbs .lt. spec_zone) THEN
! Y-start boundary
        DO j = jts, min(jtf,jbs+spec_zone-1)
          b_dist = j - jbs
          b_limit = b_dist
          IF(periodic_x)b_limit = 0
          DO k = kts, ktf
            DO i = max(its,b_limit+ibs), min(itf,ibe-b_limit)

              MU_OLD(i,j) = MUTS(i,j) - dt*MU_TEND(i,j)

              field(i,k,j) = field(i,k,j)*(c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) + &
                   dt*field_tend(i,k,j)/(c1(k)*muts(i,j)+c2(k)) +               &
                   ph_save(i,k,j)*((c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) - 1.)

            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF (jbe - jtf .lt. spec_zone) THEN
! Y-end boundary
        DO j = max(jts,jbe-spec_zone+1), jtf
          b_dist = jbe - j
          b_limit = b_dist
          IF(periodic_x)b_limit = 0
          DO k = kts, ktf
            DO i = max(its,b_limit+ibs), min(itf,ibe-b_limit)

              MU_OLD(i,j) = MUTS(i,j) - dt*MU_TEND(i,j)

              field(i,k,j) = field(i,k,j)*(c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) + &
                   dt*field_tend(i,k,j)/(c1(k)*muts(i,j)+c2(k)) +               &
                   ph_save(i,k,j)*((c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) - 1.)

            ENDDO
          ENDDO
        ENDDO
      ENDIF

    IF(.NOT.periodic_x)THEN
      IF (its - ibs .lt. spec_zone) THEN
! X-start boundary
        DO i = its, min(itf,ibs+spec_zone-1)
          b_dist = i - ibs
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)

              MU_OLD(i,j) = MUTS(i,j) - dt*MU_TEND(i,j)

              field(i,k,j) = field(i,k,j)*(c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) + &
                   dt*field_tend(i,k,j)/(c1(k)*muts(i,j)+c2(k)) +               &
                   ph_save(i,k,j)*((c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) - 1.)

            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (ibe - itf .lt. spec_zone) THEN
! X-end boundary
        DO i = max(its,ibe-spec_zone+1), itf
          b_dist = ibe - i
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)

              MU_OLD(i,j) = MUTS(i,j) - dt*MU_TEND(i,j)

              field(i,k,j) = field(i,k,j)*(c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) + &
                   dt*field_tend(i,k,j)/(c1(k)*muts(i,j)+c2(k)) +               &
                   ph_save(i,k,j)*((c1(k)*mu_old(i,j)+c2(k))/(c1(k)*muts(i,j)+c2(k)) - 1.)

            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF

   END SUBROUTINE spec_bdyupdate_ph

!------------------------------------------------------------------------

   SUBROUTINE relax_bdy_dry ( config_flags,                                    &
                              ru_tendf, rv_tendf, ph_tendf, t_tendf,           &
                              rw_tendf, mu_tend, c1h, c2h, c1f, c2f,           &
                              ru, rv, ph, t,                                   &
                              w, mu, mut,                                      &
                              u_bxs,u_bxe,u_bys,u_bye,                         &
                              v_bxs,v_bxe,v_bys,v_bye,                         &
                              ph_bxs,ph_bxe,ph_bys,ph_bye,                     &
                              t_bxs,t_bxe,t_bys,t_bye,                         &
                              w_bxs,w_bxe,w_bys,w_bye,                         &
                              mu_bxs,mu_bxe,mu_bys,mu_bye,                     &
                              u_btxs,u_btxe,u_btys,u_btye,                     &
                              v_btxs,v_btxe,v_btys,v_btye,                     &
                              ph_btxs,ph_btxe,ph_btys,ph_btye,                 &
                              t_btxs,t_btxe,t_btys,t_btye,                     &
                              w_btxs,w_btxe,w_btys,w_btye,                     &
                              mu_btxs,mu_btxe,mu_btys,mu_btye,                 &
                              spec_bdy_width, spec_zone, relax_zone,           &
                              dtbc, fcx, gcx,             &
                              ids,ide, jds,jde, kds,kde,  & ! domain dims
                              ims,ime, jms,jme, kms,kme,  & ! memory dims
                              ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                              its, ite, jts, jte, kts, kte)
   IMPLICIT NONE

   !  Input data.
   TYPE( grid_config_rec_type ) config_flags

   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte
   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, spec_zone, relax_zone

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(IN   ) :: ru,     &
                                                                      rv,     &
                                                                      ph,     &
                                                                      w,      &
                                                                      t
   REAL , DIMENSION( ims:ime , jms:jme  ) , INTENT(IN   )          :: mu  , &
                                                                      mut
   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(INOUT) :: ru_tendf, &
                                                                      rv_tendf, &
                                                                      ph_tendf, &
                                                                      rw_tendf, &
                                                                      t_tendf
   REAL , DIMENSION( ims:ime , jms:jme  ) , INTENT(INOUT)          :: mu_tend

   REAL , DIMENSION( kms:kme ) , INTENT(IN   ) :: c1h, c2h, c1f, c2f

   REAL , DIMENSION( spec_bdy_width) , INTENT(IN   ) :: fcx, gcx

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: u_bxs,u_bxe, &
                                                                               v_bxs,v_bxe, &
                                                                               ph_bxs,ph_bxe, &
                                                                               w_bxs,w_bxe, &
                                                                               t_bxs,t_bxe, &
                                                                               u_btxs,u_btxe, &
                                                                               v_btxs,v_btxe, &
                                                                               ph_btxs,ph_btxe, &
                                                                               w_btxs,w_btxe, &
                                                                               t_btxs,t_btxe

   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: u_bys,u_bye, &
                                                                               v_bys,v_bye, &
                                                                               ph_bys,ph_bye, &
                                                                               w_bys,w_bye, &
                                                                               t_bys,t_bye, &
                                                                               u_btys,u_btye, &
                                                                               v_btys,v_btye, &
                                                                               ph_btys,ph_btye, &
                                                                               w_btys,w_btye, &
                                                                               t_btys,t_btye


   REAL,  DIMENSION( jms:jme , 1:1     , spec_bdy_width    ), INTENT(IN   ) :: mu_bxs,mu_bxe, &
                                                                               mu_btxs,mu_btxe

   REAL,  DIMENSION( ims:ime , 1:1     , spec_bdy_width    ), INTENT(IN   ) :: mu_bys,mu_bye, &
                                                                               mu_btys,mu_btye
   REAL, INTENT(IN   ) :: dtbc

! changed to tile dimensions, 20090923, JM
   REAL , DIMENSION( its-1:ite+1 , kts:kte, jts-1:jte+1  ) :: rfield
   INTEGER :: i_start, i_end, j_start, j_end, i, j, k

           CALL relax_bdytend ( ru, ru_tendf,             &
                               u_bxs,u_bxe,u_bys,u_bye,u_btxs,u_btxe,u_btys,u_btye, &
                               'u'        , config_flags, &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
           CALL relax_bdytend ( rv, rv_tendf,             &
                               v_bxs,v_bxe,v_bys,v_bye,v_btxs,v_btxe,v_btys,v_btye, &
                               'v'        , config_flags, &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

! rfield will be calculated beyond tile limits because relax_bdytend
!   requires a 5-point stencil, and this avoids need for inter-tile/patch
!   communication here
           i_start = max(its-1, ids)
           i_end = min(ite+1, ide-1)
           j_start = max(jts-1, jds)
           j_end = min(jte+1, jde-1)

           CALL mass_weight ( ph , mut , rfield , c1f, c2f, &
                              ids,ide, jds,jde, kds,kde,    &        ! domain dims
                              ims,ime, jms,jme, kms,kme,    &        ! memory dims
                              its-1,ite+1 , jts-1,jte+1 ,   &        ! rfield dims
                              kts,kte,                      &        ! rfield
                              i_start,i_end, j_start,j_end, kts,kte)  ! tile dims

           CALL relax_bdytend_tile ( rfield, ph_tendf,             &
                               ph_bxs,ph_bxe,ph_bys,ph_bye, ph_btxs,ph_btxe,ph_btys,ph_btye,       &
                               'h'        , config_flags, &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte,  &
                               its-1, ite+1, jts-1,jte+1,kts,kte )  ! dims of first argument

           CALL mass_weight ( t,   mut , rfield , c1h, c2h, &
                              ids,ide, jds,jde, kds,kde,    &         ! domain dims
                              ims,ime, jms,jme, kms,kme,    &         ! memory dims
                              its-1,ite+1 , jts-1,jte+1 ,   &        ! rfield dims
                              kts,kte,                      &        ! rfield
                              i_start,i_end, j_start,j_end, kts,kte-1) ! tile dims

           CALL relax_bdytend_tile ( rfield, t_tendf,              &
                               t_bxs,t_bxe,t_bys,t_bye, t_btxs,t_btxe,t_btys,t_btye,       &
                               't'        , config_flags, &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte,  &
                               its-1, ite+1, jts-1,jte+1,kts,kte )  ! dims of first argument

           CALL relax_bdytend ( mu, mu_tend,               &
                               mu_bxs,mu_bxe,mu_bys,mu_bye, mu_btxs,mu_btxe,mu_btys,mu_btye,                &
                               'm'        , config_flags,  &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, 1  ,1  ,  & ! domain dims
                               ims,ime, jms,jme, 1  ,1  ,  & ! memory dims
                               ips,ipe, jps,jpe, 1  ,1  ,  & ! patch  dims
                               its,ite, jts,jte, 1  ,1   )

         IF( config_flags%nested) THEN

           i_start = max(its-1, ids)
           i_end = min(ite+1, ide-1)
           j_start = max(jts-1, jds)
           j_end = min(jte+1, jde-1)

           CALL mass_weight ( w , mut , rfield , c1f, c2f, &
                              ids,ide, jds,jde, kds,kde,   &        ! domain dims
                              ims,ime, jms,jme, kms,kme,   &        ! memory dims
                              its-1,ite+1 , jts-1,jte+1 ,  &        ! rfield dims
                              kts,kte,                     &        ! rfield
                              i_start,i_end, j_start,j_end, kts,kte) ! tile dims

           CALL relax_bdytend_tile ( rfield, rw_tendf,             &
                               w_bxs,w_bxe,w_bys,w_bye, w_btxs,w_btxe,w_btys,w_btye,       &
                               'h'        , config_flags, &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte,  &
                               its-1, ite+1, jts-1,jte+1,kts,kte )  ! dims of first argument

        END IF
   END SUBROUTINE relax_bdy_dry
!------------------------------------------------------------------------
   SUBROUTINE relax_bdy_scalar ( scalar_tend,                &
                                 scalar, mu, c1h, c2h,       &
                                 scalar_bxs,scalar_bxe,scalar_bys,scalar_bye, &
                                 scalar_btxs,scalar_btxe,scalar_btys,scalar_btye, &
                                 spec_bdy_width, spec_zone, relax_zone,       &
                                 dtbc, fcx, gcx,             &
                                 config_flags,               &
                                 ids,ide, jds,jde, kds,kde,  & ! domain dims
                                 ims,ime, jms,jme, kms,kme,  & ! memory dims
                                 ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                 its, ite, jts, jte, kts, kte)
   IMPLICIT NONE

   !  Input data.
   TYPE( grid_config_rec_type ) config_flags

   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte
   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, spec_zone, relax_zone

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(IN   ) :: scalar
   REAL , DIMENSION( ims:ime , jms:jme  ) , INTENT(IN   ) :: mu
   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(INOUT) :: scalar_tend
   REAL , DIMENSION( spec_bdy_width) , INTENT(IN   ) :: fcx, gcx

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: scalar_bxs,scalar_bxe, &
                                                                               scalar_btxs,scalar_btxe
   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: scalar_bys,scalar_bye, &
                                                                               scalar_btys,scalar_btye
   REAL , DIMENSION( kms:kme ) , INTENT(IN   ) :: c1h, c2h
   REAL, INTENT(IN   ) :: dtbc
!Local
   INTEGER :: i,j,k, i_start, i_end, j_start, j_end
   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) :: rscalar

! rscalar will be calculated beyond tile limits because relax_bdytend
!   requires a 5-point stencil, and this avoids need for inter-tile/patch
!   communication here
           i_start = max(its-1, ids)
           i_end = min(ite+1, ide-1)
           j_start = max(jts-1, jds)
           j_end = min(jte+1, jde-1)

           CALL mass_weight ( scalar , mu , rscalar, c1h, c2h, &
                              ids,ide, jds,jde, kds,kde,   &          ! domain dims
                              ims,ime, jms,jme, kms,kme,   &          ! memory dims
                              ims,ime, jms,jme, kms,kme,   &          ! rfield dims
                              i_start,i_end, j_start,j_end, kts,kte-1) ! tile dims

           CALL relax_bdytend (rscalar, scalar_tend,             &
                               scalar_bxs,scalar_bxe,scalar_bys,scalar_bye, scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,       &
                               'q'        , config_flags, &
                               spec_bdy_width, spec_zone, relax_zone, &
                               dtbc, fcx, gcx,             &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )


   END SUBROUTINE relax_bdy_scalar

!------------------------------------------------------------------------
   SUBROUTINE spec_bdy_dry ( config_flags,                        &
                             ru_tend, rv_tend, ph_tend, t_tend,   &
                             rw_tend, mu_tend,                    &
                             u_bxs,u_bxe,u_bys,u_bye,             &
                             v_bxs,v_bxe,v_bys,v_bye,             &
                             ph_bxs,ph_bxe,ph_bys,ph_bye,         &
                             t_bxs,t_bxe,t_bys,t_bye,             &
                             w_bxs,w_bxe,w_bys,w_bye,             &
                             mu_bxs,mu_bxe,mu_bys,mu_bye,         &
                             u_btxs,u_btxe,u_btys,u_btye,         &
                             v_btxs,v_btxe,v_btys,v_btye,         &
                             ph_btxs,ph_btxe,ph_btys,ph_btye,     &
                             t_btxs,t_btxe,t_btys,t_btye,         &
                             w_btxs,w_btxe,w_btys,w_btye,         &
                             mu_btxs,mu_btxe,mu_btys,mu_btye,     &
                             spec_bdy_width, spec_zone,           &
                             ids,ide, jds,jde, kds,kde,  & ! domain dims
                             ims,ime, jms,jme, kms,kme,  & ! memory dims
                             ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                             its, ite, jts, jte, kts, kte)
   IMPLICIT NONE

   !  Input data.
   TYPE( grid_config_rec_type ) config_flags


   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte
   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, spec_zone

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(OUT  ) :: ru_tend, &
                                                                      rv_tend, &
                                                                      ph_tend, &
                                                                      rw_tend, &
                                                                      t_tend
   REAL , DIMENSION( ims:ime , jms:jme  ) , INTENT(OUT  )          :: mu_tend

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: u_bxs,u_bxe,  &
                                                                               v_bxs,v_bxe,  &
                                                                               ph_bxs,ph_bxe, &
                                                                               w_bxs,w_bxe, &
                                                                               t_bxs,t_bxe,  &
                                                                               u_btxs,u_btxe, &
                                                                               v_btxs,v_btxe, &
                                                                               ph_btxs,ph_btxe, &
                                                                               w_btxs,w_btxe, &
                                                                               t_btxs,t_btxe

   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: u_bys,u_bye,  &
                                                                               v_bys,v_bye,  &
                                                                               ph_bys,ph_bye, &
                                                                               w_bys,w_bye, &
                                                                               t_bys,t_bye,  &
                                                                               u_btys,u_btye, &
                                                                               v_btys,v_btye, &
                                                                               ph_btys,ph_btye, &
                                                                               w_btys,w_btye, &
                                                                               t_btys,t_btye

   REAL,  DIMENSION( jms:jme , 1:1 ,     spec_bdy_width    ), INTENT(IN   ) :: mu_bxs,mu_bxe, &
                                                                               mu_btxs,mu_btxe

   REAL,  DIMENSION( ims:ime , 1:1 ,     spec_bdy_width    ), INTENT(IN   ) :: mu_bys,mu_bye, &
                                                                               mu_btys,mu_btye
         CALL spec_bdytend (   ru_tend,                &
                               u_bxs,u_bxe,u_bys,u_bye, u_btxs,u_btxe,u_btys,u_btye,    &
                               'u'     , config_flags, &
                               spec_bdy_width, spec_zone, &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
         CALL spec_bdytend (   rv_tend,                &
                               v_bxs,v_bxe,v_bys,v_bye, v_btxs,v_btxe,v_btys,v_btye,    &
                               'v'     , config_flags, &
                               spec_bdy_width, spec_zone, &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
         CALL spec_bdytend (   ph_tend,                &
                               ph_bxs,ph_bxe,ph_bys,ph_bye, ph_btxs,ph_btxe,ph_btys,ph_btye,    &
                               'h'     , config_flags, &
                               spec_bdy_width, spec_zone, &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
         CALL spec_bdytend (   t_tend,                &
                               t_bxs,t_bxe,t_bys,t_bye, t_btxs,t_btxe,t_btys,t_btye,    &
                               't'     , config_flags, &
                               spec_bdy_width, spec_zone, &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
         CALL spec_bdytend (   mu_tend,                &
                               mu_bxs,mu_bxe,mu_bys,mu_bye, mu_btxs,mu_btxe,mu_btys,mu_btye,       &
                               'm'     , config_flags, &
                               spec_bdy_width, spec_zone, &
                               ids,ide, jds,jde, 1  ,1  ,  & ! domain dims
                               ims,ime, jms,jme, 1  ,1  ,  & ! memory dims
                               ips,ipe, jps,jpe, 1  ,1  ,  & ! patch  dims
                               its,ite, jts,jte, 1  ,1   )

         if(config_flags%nested)                           &
         CALL spec_bdytend (   rw_tend,                    &
                               w_bxs,w_bxe,w_bys,w_bye, w_btxs,w_btxe,w_btys,w_btye,                  &
                               'h'     , config_flags,     &
                               spec_bdy_width, spec_zone,  &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

   END SUBROUTINE spec_bdy_dry


!------------------------------------------------------------------------
! KRS 9/13/2012: New subroutine spec_bdy_dry_perturb.
! Perturbation field passed in: field_(u,v,t)_tend_perturb
! field_(u,v,t)_tend_perturb=r(u,v,t)_tendf_stoch if perturb_bdy=1
! field_(u,v,t)_tend_perturb=User provided patterns if perturb_bdy=2
! This routine calls spec_bdytend_perturb in share/module_bc.F for u,v,t.
!------------------------------------------------------------------------

   SUBROUTINE spec_bdy_dry_perturb ( config_flags,                                   &
                                     ru_tend, rv_tend, t_tend,mu_2, mub, c1, c2,     &
                    msfu, msfv, msft,                               &
                                     field_u_tend_perturb,field_v_tend_perturb,field_t_tend_perturb, &
      spec_bdy_width, spec_zone,  &
                                     kme_stoch,                  & ! stoch dims
                                     ids,ide, jds,jde, kds,kde,  & ! domain dims
      ims,ime, jms,jme, kms,kme,  & ! memory dims
                                     ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                     its, ite, jts, jte, kts, kte)
   IMPLICIT NONE

   !  Input data.
   TYPE( grid_config_rec_type ) config_flags


   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte, &
                                            kme_stoch
   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, spec_zone

   REAL , DIMENSION( ims:ime , kms:kme ,jms:jme  ) , INTENT(INOUT) :: ru_tend, &
                                                                      rv_tend, &
                                                                      t_tend

   REAL,  DIMENSION( ims:ime ,           jms:jme ), INTENT(INOUT  )   :: mu_2
   REAL,  DIMENSION( ims:ime ,           jms:jme ), INTENT(INOUT  )   :: mub
   REAL,  DIMENSION( ims:ime ,           jms:jme ), INTENT(IN     )   :: msfu
   REAL,  DIMENSION( ims:ime ,           jms:jme ), INTENT(IN     )   :: msfv
   REAL,  DIMENSION( ims:ime ,           jms:jme ), INTENT(IN     )   :: msft
   REAL,  DIMENSION(           kms:kme           ), INTENT(IN     )   :: c1, c2

   REAL,  DIMENSION( ims:ime , kms:kme_stoch , jms:jme ), INTENT(IN     )   :: field_u_tend_perturb, &
                                                                         field_v_tend_perturb, &
                                                                         field_t_tend_perturb

         CALL spec_bdytend_perturb (   ru_tend,            &
                               field_u_tend_perturb,       &
                               mu_2,mub, c1, c2,           &
                               'u', msfu, config_flags,    &
                               spec_bdy_width, spec_zone,  &
                               kme_stoch,                  & ! stoch  dims
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )
         CALL spec_bdytend_perturb (   rv_tend,            &
                               field_v_tend_perturb,       &
                               mu_2,mub, c1, c2,           &
                               'v', msfv, config_flags,    &
                               spec_bdy_width, spec_zone,  &
                               kme_stoch,                  & ! stoch  dims
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

         CALL spec_bdytend_perturb (   t_tend,             &
                               field_t_tend_perturb,       &
                               mu_2,mub, c1, c2,           &
                               't', msft, config_flags,    &
                               spec_bdy_width, spec_zone,  &
                               kme_stoch,                  & ! stoch  dims
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )

   END SUBROUTINE spec_bdy_dry_perturb

!------------------------------------------------------------------------
!----------------------------- START LJMB -------------------------------
! 6/19/2017: New subroutines:
! 1) "mk_cells"
!    Creates an array of the size of the domain (grid%mom_pert_u/v).
!    The array contains perturbations all over organized as cells.
! 2) "mom_pert"
!    Modifies u, v and w tendencies (ru_tend, rv_tend, rw_tend).
!    Adds random perturbations to cells of a certain size @ inflow bdries.
!------------------------------------------------------------------------

   SUBROUTINE mk_cells_u (mom_pert_u,mom_pert_ktop,                     &
                        mom_pert_kbot,cell_width,                       &
                        mom_pert_amp,config_flags,grid,                 &
                        its,ite,jts,jte,kts,kte,                        &
                        ids,ide,jds,jde,kds,kde,                        &
                        ips,ipe,jps,jpe,kps,kpe,                        &
                        ims,ime,jms,jme,kms,kme)

      ! __________ DECLARATIONS ___________
      IMPLICIT NONE

      ! Array dimentions:
      INTEGER,INTENT(IN)::ids,ide,jds,jde,kds,kde,                      &
                          its,ite,jts,jte,kts,kte,                      &
                          ips,ipe,jps,jpe,kps,kpe,                      &
                          ims,ime,jms,jme,kms,kme

      ! Cell characteristics:
      REAL,INTENT(IN)::mom_pert_ktop,mom_pert_kbot,mom_pert_amp
      INTEGER,INTENT(IN)::cell_width

      ! Grid pointer for intra-processor communications:
      TYPE(domain),INTENT(IN)::grid

      ! Structure that contains run-time configuration (namelist) data for
      ! domain
      TYPE(grid_config_rec_type),INTENT(IN)::config_flags

      ! Array of Perturbations:
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT)::mom_pert_u

      ! Other variables:
      INTEGER::i,j,k,icell,jcell,kcell
      REAL::rand_amp,cell_pt_x,cell_pt_y,cell_pt_z

      ! ________ END DECLARATIONS _________

      CALL random_seed ! Must be called to randomnize the random # generator (-_-)

      ! Fill up array:
      ! - Begin by finding cell boundaries within tile.
      ! - After finding the boundaries, fill up the rest of that cell 

      !For each perturbed k:
      DO k = mom_pert_kbot,mom_pert_ktop

         ! For all x's in the tile:
         DO i = its,ite-1

            ! If this z falls on the transition between cells: 
            IF (MOD(k-1,cell_width) .EQ. 0) THEN

               ! For all y's in the tile:
               DO j = jts,jte-1

                  ! If this y falls on the transition between cells:
                  IF (MOD(j-1,cell_width) .EQ. 0) THEN

                     ! Define a random number from -1 to 1:
                     CALL random_number(rand_amp)
                     rand_amp = ((rand_amp*2)-1)*mom_pert_amp

                     ! For all the gps in this cell (or rest of tile) in z:
                     DO kcell = k,min(k+cell_width,kte)

                        IF (kcell .LE. mom_pert_ktop) THEN

                           ! For all the gps in this cell (or rest of tile) in y:
                           DO jcell = j,min(j+cell_width,jte)

                              ! Assign value to this cell in the final array!
                              mom_pert_u(i,kcell,jcell) = rand_amp

                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

! Pass information between processors:
#ifdef DM_PARALLEL
#      include "HALO_EM_MOM_PERT_U.inc"
#      include "PERIOD_EM_MOM_PERT_U.inc"
#endif

      ! Now fill the gaps between tiles!

      ! For all perturbed points in x:
      DO i = its,ite

         ! If the first point of the tile falls on a y-cell intersection but not z: 
         IF ((MOD(kts-1,cell_width) .NE. 0) .AND. (MOD(jts-1,cell_width) .EQ. 0)) THEN

            ! its with respect to current cell:
            cell_pt_z = cell_width - MOD(kts-1,cell_width) -1

            ! For all z-points from the start of the tile to the end of the cell:
            DO kcell = kts,kts+cell_pt_z

               IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                  ! For all y-points in the tile:
                  DO jcell = jts,jte

                     ! Assign value to this grid-cell:
                     mom_pert_u(i,kcell,jcell) = mom_pert_u(i,kts-1,jcell)

                  ENDDO
               ENDIF
            ENDDO

         ! If the first point of the tile falls on a z-cell intersection but not y: 
         ELSEIF ((MOD(kts-1,cell_width) .EQ. 0) .AND. (MOD(jts-1,cell_width) .NE. 0)) THEN

            ! jts with respect to current cell:
            cell_pt_y = cell_width - MOD(jts-1,cell_width) - 1

            ! For all y-points from the start of the tile to the end of the cell:
            DO jcell = jts,jts+cell_pt_y

               ! For all k-points in the tile:
               DO kcell = kts,kte
                  IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN
                     ! Assign value to this grid-cell:
                     mom_pert_u(i,kcell,jcell) = mom_pert_u(i,kcell,jts-1)
                  ENDIF
               ENDDO
            ENDDO

         ! If the first point of the tile doesn't fall on an z- or y-cell intersection: 
         ELSEIF ((MOD(kts-1,cell_width) .NE. 0) .AND. (MOD(jts-1,cell_width) .NE. 0)) THEN

            ! kts and jts with respect to current cell:
            cell_pt_z = cell_width - MOD(kts-1,cell_width) - 1
            cell_pt_y = cell_width - MOD(jts-1,cell_width) - 1

            ! For all z-points from the start of the tile to the end of the cell:
            DO kcell = kts,kts+cell_pt_z

               IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                  ! Now for y: 
                  DO jcell = jts,jts+cell_pt_y

                     ! Assign value to this grid-cell:
                     mom_pert_u(i,kcell,jcell) = mom_pert_u(i,kts-1,jts-1)

                  ENDDO

               ENDIF

            ENDDO

            ! For z's after end of cell to end of tile:
            DO kcell = kts+cell_pt_z+1,kte

               IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                  !For lower y-boundary:
                  DO jcell = jts,jts+cell_pt_y

                     ! Assign values! 
                     mom_pert_u(i,kcell,jcell) = mom_pert_u(i,kcell,jts-1)

                  ENDDO

               ENDIF

            ENDDO

            ! For y's after end of cell to end of tile:
            DO jcell = jts+cell_pt_y+1,jte

               !For lower x-boundary:
               DO kcell = kts,kts+cell_pt_z
                  IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN
                     ! Assign values! 
                     mom_pert_u(i,kcell,jcell) = mom_pert_u(i,kts-1,jcell)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE mk_cells_u

!------------------------------------------------------------------------

   SUBROUTINE mk_cells_v (mom_pert_v,mom_pert_ktop,                     &
                        mom_pert_kbot,cell_width,                       &
                        mom_pert_amp,config_flags,grid,                 &
                        its,ite,jts,jte,kts,kte,                        &
                        ids,ide,jds,jde,kds,kde,                        &
                        ips,ipe,jps,jpe,kps,kpe,                        &
                        ims,ime,jms,jme,kms,kme)

      ! __________ DECLARATIONS ___________
      IMPLICIT NONE

      ! Array dimentions:
      INTEGER,INTENT(IN)::ids,ide,jds,jde,kds,kde,                      &
                          its,ite,jts,jte,kts,kte,                      &
                          ips,ipe,jps,jpe,kps,kpe,                      &
                          ims,ime,jms,jme,kms,kme

      ! Cell characteristics:
      REAL,INTENT(IN)::mom_pert_ktop,mom_pert_kbot,mom_pert_amp
      INTEGER,INTENT(IN)::cell_width

      ! Grid pointer for intra-processor communications:
      TYPE(domain),INTENT(IN)::grid

      ! Structure that contains run-time configuration (namelist) data for
      ! domain
      TYPE(grid_config_rec_type),INTENT(IN)::config_flags

      ! Array of Perturbations:
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT)::mom_pert_v

      ! Other variables:
      INTEGER::i,j,k,icell,jcell,kcell
      REAL::rand_amp,cell_pt_x,cell_pt_y,cell_pt_z

      ! ________ END DECLARATIONS _________

      ! We will not call it for y, or our perturbations will look the same as x
      CALL random_seed ! Must be called to randomnize the random # generator (-_-)

      ! Fill up array:
      ! - Begin by finding cell boundaries within tile.
      ! - After finding the boundaries, fill up the rest of that cell 

      !For each perturbed z:
      DO k = mom_pert_kbot,mom_pert_ktop

         ! For all x's in the tile:
         DO i = its,ite-1

            ! If this x falls on the transition between cells: 
            IF (MOD(i-1,cell_width) .EQ. 0) THEN

               ! For all y's in the tile:
               DO j = jts,jte-1

                  ! If this y falls on the transition between cells:
                  IF (MOD(k-1,cell_width) .EQ. 0) THEN

                     ! Define a random number from -1 to 1:
                     CALL random_number(rand_amp)
                     rand_amp = ((rand_amp*2)-1)*mom_pert_amp

                     ! For all the gps in this cell (or rest of tile) in x:
                     DO icell = i,min(i+cell_width,ite)

                        ! For all the gps in this cell (or rest of tile) in y:
                        DO kcell = k,min(k+cell_width,kte)

                           IF (kcell .LE. mom_pert_ktop) THEN

                              ! Assign value to this cell in the final array!
                              mom_pert_v(icell,kcell,j) = rand_amp

                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

! Pass information between processors:
#ifdef DM_PARALLEL
#      include "HALO_EM_MOM_PERT_V.inc"
#      include "PERIOD_EM_MOM_PERT_V.inc"
#endif

      ! Now fill the gaps between tiles!

      ! For all perturbed points in y:
      DO j = jts,jte

         ! If the first point of the tile falls on a y-cell intersection but not x: 
         IF ((MOD(its-1,cell_width) .NE. 0) .AND. (MOD(kts-1,cell_width) .EQ. 0)) THEN

            ! its with respect to current cell:
            cell_pt_x = cell_width - MOD(its-1,cell_width) -1

            ! For all x-points from the start of the tile to the end of the cell:
            DO icell = its,its+cell_pt_x

               ! For all z-points in the tile:
               DO kcell = kts,kte

                  IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                     ! Assign value to this grid-cell:
                     mom_pert_v(icell,kcell,j) = mom_pert_v(its-1,kcell,j)

                  ENDIF
               ENDDO
            ENDDO

         ! If the first point of the tile falls on a x-cell intersection but not y: 
         ELSEIF ((MOD(its-1,cell_width) .EQ. 0) .AND. (MOD(kts-1,cell_width) .NE. 0)) THEN

            ! jts with respect to current cell:
            cell_pt_z = cell_width - MOD(kts-1,cell_width) - 1

            ! For all y-points from the start of the tile to the end of the cell:
            DO kcell = kts,kts+cell_pt_z

               IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                  ! For all x-points in the tile:
                  DO icell = its,ite

                     ! Assign value to this grid-cell:
                     mom_pert_v(icell,kcell,j) = mom_pert_v(icell,kts-1,j)

                  ENDDO
               ENDIF
            ENDDO

         ! If the first point of the tile doesn't fall on an x- or y-cell intersection: 
         ELSEIF ((MOD(its-1,cell_width) .NE. 0) .AND. (MOD(kts-1,cell_width) .NE. 0)) THEN

            ! its and jts with respect to current cell:
            cell_pt_x = cell_width - MOD(its-1,cell_width) - 1
            cell_pt_z = cell_width - MOD(kts-1,cell_width) - 1

            ! For all x-points from the start of the tile to the end of the cell:
            DO icell = its,its+cell_pt_x

               ! Now for y: 
               DO kcell = kts,kts+cell_pt_z
                  IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                     ! Assign value to this grid-cell:
                     mom_pert_v(icell,kcell,j) = mom_pert_v(its-1,kts-1,j)

                  ENDIF
               ENDDO
            ENDDO

            ! For x's after end of cell to end of tile:
            DO icell = its+cell_pt_x+1,ite

               !For lower y-boundary:
               DO kcell = kts,kts+cell_pt_z
                  IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN
                     ! Assign values! 
                     mom_pert_v(icell,kcell,j) = mom_pert_v(icell,kts-1,j)

                  ENDIF
               ENDDO
            ENDDO

            ! For y's after end of cell to end of tile:
            DO kcell = kts+cell_pt_z+1,kte
               IF ((kcell .GE. mom_pert_kbot).AND.(kcell .LE. mom_pert_ktop)) THEN

                  !For lower x-boundary:
                  DO icell = its,its+cell_pt_x

                     ! Assign values! 
                     mom_pert_v(icell,kcell,j) = mom_pert_v(its-1,kcell,j)

                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE mk_cells_v

!------------------------------------------------------------------------

   SUBROUTINE mk_cells_w (mom_pert_w,mom_pert_ktop,                     &
                        mom_pert_kbot,cell_width,                       &
                        mom_pert_amp,config_flags,grid,                 &
                        its,ite,jts,jte,kts,kte,                        &
                        ids,ide,jds,jde,kds,kde,                        &
                        ips,ipe,jps,jpe,kps,kpe,                        &
                        ims,ime,jms,jme,kms,kme)

      ! __________ DECLARATIONS ___________
      IMPLICIT NONE

      ! Array dimentions:
      INTEGER,INTENT(IN)::ids,ide,jds,jde,kds,kde,                      &
                          its,ite,jts,jte,kts,kte,                      &
                          ips,ipe,jps,jpe,kps,kpe,                      &
                          ims,ime,jms,jme,kms,kme

      ! Cell characteristics:
      REAL,INTENT(IN)::mom_pert_ktop,mom_pert_kbot,mom_pert_amp
      INTEGER,INTENT(IN)::cell_width

      ! Grid pointer for intra-processor communications:
      TYPE(domain),INTENT(IN)::grid

      ! Structure that contains run-time configuration (namelist) data for domain
      TYPE(grid_config_rec_type),INTENT(IN)::config_flags

      ! Array of Perturbations:
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT)::mom_pert_w

      ! Other variables:
      INTEGER::i,j,k,icell,jcell
      REAL::rand_amp,cell_pt_x,cell_pt_y

      ! ________ END DECLARATIONS _________

      ! We will not call it for y, or our perturbations will look the same as x
      CALL random_seed ! Must be called to randomnize the random # generator (-_-)

      ! Fill up array:
      ! - Begin by finding cell boundaries within tile.
      ! - After finding the boundaries, fill up the rest of that cell

      !For each perturbed z:
      DO k = mom_pert_kbot,mom_pert_ktop

         ! For all x's in the tile:
         DO i = its,ite-1

            ! If this x falls on the transition between cells:
            IF (MOD(i-1,cell_width) .EQ. 0) THEN

               ! For all y's in the tile:
               DO j = jts,jte-1

                  ! If this y falls on the transition between cells:
                  IF (MOD(j-1,cell_width) .EQ. 0) THEN

                     ! Define a random number from -1 to 1:
                     CALL random_number(rand_amp)
                     rand_amp = ((rand_amp*2)-1)*mom_pert_amp

                     ! For all the gps in this cell (or rest of tile) in x:
                     DO icell = i,min(i+cell_width,ite)

                        ! For all the gps in this cell (or rest of tile) in y:
                        DO jcell = j,min(j+cell_width,jte)

                           ! Assign value to this cell in the final array!
                           mom_pert_w(icell,k,jcell) = rand_amp

                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

! Pass information between processors:
#ifdef DM_PARALLEL
#      include "HALO_EM_MOM_PERT_W.inc"
#      include "PERIOD_EM_MOM_PERT_W.inc"
#endif

      ! Now fill the gaps between tiles!

      ! For all perturbed points in z:
      DO k = mom_pert_kbot,mom_pert_ktop

         ! If the first point of the tile falls on a y-cell intersection but not
         ! x:
         IF ((MOD(its-1,cell_width) .NE. 0) .AND. (MOD(jts-1,cell_width) .EQ. 0)) THEN

            ! its with respect to current cell:
            cell_pt_x = cell_width - MOD(its-1,cell_width) -1

            ! For all x-points from the start of the tile to the end of the
            ! cell:
            DO icell = its,its+cell_pt_x

               ! For all y-points in the tile:
               DO jcell = jts,jte

                  ! Assign value to this grid-cell:
                  mom_pert_w(icell,k,jcell) = mom_pert_w(its-1,k,jcell)

               ENDDO
            ENDDO

         ! If the first point of the tile falls on a x-cell intersection but not
         ! y:
         ELSEIF ((MOD(its-1,cell_width) .EQ. 0) .AND. (MOD(jts-1,cell_width) .NE. 0)) THEN

            ! jts with respect to current cell:
            cell_pt_y = cell_width - MOD(jts-1,cell_width) - 1

            ! For all y-points from the start of the tile to the end of the
            ! cell:
            DO jcell = jts,jts+cell_pt_y

               ! For all x-points in the tile:
               DO icell = its,ite

                  ! Assign value to this grid-cell:
                  mom_pert_w(icell,k,jcell) = mom_pert_w(icell,k,jts-1)

               ENDDO
            ENDDO

         ! If the first point of the tile doesn't fall on an x- or y-cell
         ! intersection:
         ELSEIF ((MOD(its-1,cell_width) .NE. 0) .AND. (MOD(jts-1,cell_width) .NE. 0)) THEN

            ! its and jts with respect to current cell:
            cell_pt_x = cell_width - MOD(its-1,cell_width) - 1
            cell_pt_y = cell_width - MOD(jts-1,cell_width) - 1

            ! For all x-points from the start of the tile to the end of the
            ! cell:
            DO icell = its,its+cell_pt_x

               ! Now for y:
               DO jcell = jts,jts+cell_pt_y

                  ! Assign value to this grid-cell:
                  mom_pert_w(icell,k,jcell) = mom_pert_w(its-1,k,jts-1)

               ENDDO
            ENDDO

            ! For x's after end of cell to end of tile:
            DO icell = its+cell_pt_x+1,ite

               !For lower y-boundary:
               DO jcell = jts,jts+cell_pt_y

                  ! Assign values!
                  mom_pert_w(icell,k,jcell) = mom_pert_w(icell,k,jts-1)

               ENDDO
            ENDDO

            ! For y's after end of cell to end of tile:
            DO jcell = jts+cell_pt_y+1,jte

               !For lower x-boundary:
               DO icell = its,its+cell_pt_x

                  ! Assign values!
                  mom_pert_w(icell,k,jcell) = mom_pert_w(its-1,k,jcell)

               ENDDO
            ENDDO
         ENDIF
      ENDDO

   END SUBROUTINE mk_cells_w

! -----------------------------------------------------------------------

   SUBROUTINE mom_pert (field_tend,                                     &
                        mom_pert_uv,config_flags,                       &
                        its,ite,jts,jte,kts,kte,                        &
                        ids,ide,jds,jde,kds,kde,                        &
                        ims,ime,jms,jme,kms,kme)

      ! __________ DECLARATIONS ___________
      IMPLICIT NONE

      ! Structure that contains run-time configuration (namelist) data for domain
      TYPE(grid_config_rec_type),INTENT(IN)::config_flags

      ! Tile dimentions:
      INTEGER,INTENT(IN)::its,ite,jts,jte,kts,kte

      ! Memory dimentions:
      INTEGER,INTENT(IN)::ims,ime,jms,jme,kms,kme

      ! Domain dimentions:
      INTEGER,INTENT(IN)::ids,ide,jds,jde,kds,kde

      ! Input tendency array to be perturbed:
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT)::field_tend

      ! Output array of perturbations (perts only on boundaries):
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT)::mom_pert_uv

      ! Other variables:
      INTEGER :: i,j,k,mom_pert_gps,m_pert_locx,m_pert_locy
      REAL,DIMENSION(ims:ime,kms:kme,jms:jme)::temp_mom_pert


      ! ________ END DECLARATIONS _________

      ! Store temporary array with previously calculated perturbations
      temp_mom_pert = mom_pert_uv

      ! Reset perturbation arrays to all zeros (before adding perts @ bdries)
      mom_pert_uv(:,:,:) = 0

      ! Extract number of perturbed gps at boundaries:
      mom_pert_gps = config_flags%mom_pert_gps

      ! Extract location of perturbed boundaries:
      m_pert_locx = config_flags%m_pert_locx
      m_pert_locy = config_flags%m_pert_locy

!------ If the current memory size includes the perturbed region in -x:
      IF (its-ids .LT. mom_pert_gps-1 ) THEN

         ! If -x is perturbed:
         IF (m_pert_locx .EQ. -1) THEN

              ! For all z's:
              DO k = kts,kte

                 ! For all perturbed x's:
                 DO i = its,MIN(ite,ids+mom_pert_gps-1)

                    ! For all y's:
                    DO j = jts,jte

                       ! Add perturbation where needed:
                       field_tend(i,k,j) = field_tend(i,k,j) +                &
                                           temp_mom_pert(i,k,j)
                       mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +              &
                                            temp_mom_pert(i,k,j)

                    ENDDO
                 ENDDO

                 ! If the current memory size ALSO includes the perturbed region in -y:
                 IF ((jts-jds .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. -1)) THEN

                    ! For all x's (ignoring the perturbed x-region, since it was
                    ! already perturbed):
                    DO i = ids+mom_pert_gps,ite

                       ! For all perturbed y's:
                       DO j = jts,MIN(jte,jds+mom_pert_gps-1)

                          ! Add perturbation where needed:
                          field_tend(i,k,j) = field_tend(i,k,j) +             &
                                              temp_mom_pert(i,k,j)
                          mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                               temp_mom_pert(i,k,j)

                       ENDDO
                    ENDDO

                 ! If the current memory size ALSO includes the perturbed region in +y:
                 ELSEIF ((jde-jte .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. 1)) THEN

                    ! For all x's (ignoring the perturbed x-region, since it was
                    ! already perturbed):
                    DO i = ids+mom_pert_gps,ite

                       ! For all perturbed y's:
                       DO j = MAX(jts,jde-mom_pert_gps+1),jte

                          ! Add perturbation where needed:
                          field_tend(i,k,j) = field_tend(i,k,j) +             &
                                              temp_mom_pert(i,k,j)
                          mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                               temp_mom_pert(i,k,j)

                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO

         ! If -x is NOT perturbed:
         ELSEIF (m_pert_locx .NE. -1) THEN

              ! If the current memory size includes the perturbed region in -y:
              IF ((jts-jds .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. -1)) THEN

                 ! For all z's:
                 DO k = kts,kte

                    ! For all x's: 
                    DO i = its,ite

                       ! For all perturbed y's:
                       DO j = jts,MIN(jte,jds+mom_pert_gps-1)

                          ! Add perturbation where needed:
                          field_tend(i,k,j) = field_tend(i,k,j) +             &
                                              temp_mom_pert(i,k,j)
                          mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                               temp_mom_pert(i,k,j)

                       ENDDO
                    ENDDO
                 ENDDO

              ! If the current memory size includes the perturbed region in +y:
              ELSEIF ((jde-jte .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. 1)) THEN

                 ! For all z's:
                 DO k = kts,kte

                    ! For all x's:
                    DO i = its,ite

                       ! For all perturbed y's:
                       DO j = MAX(jts,jde-mom_pert_gps+1),jte

                          ! Add perturbation where needed:
                          field_tend(i,k,j) = field_tend(i,k,j) +             &
                                              temp_mom_pert(i,k,j)
                          mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                               temp_mom_pert(i,k,j)

                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
         ENDIF

!---- If the current memory size does not include the perturbed +x pr -x region,
      ! but it includes the perturbed region in -y:
      ELSEIF ((jts-jds .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. -1)) THEN

         ! For all z's:
         DO k = kts,kte

            ! For all x's:
            DO i = its,ite

               ! For all perturbed y's:
               DO j = jts,MIN(jte,jds+mom_pert_gps-1)

                  ! Add perturbation where needed:
                  field_tend(i,k,j) = field_tend(i,k,j) +                &
                                      temp_mom_pert(i,k,j)
                  mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +              &
                                       temp_mom_pert(i,k,j)

               ENDDO
            ENDDO
         ENDDO

      ! If the current memory size does not include the perturbed +x or -x region,
      ! but it includes the perturbed region in +y:
      ELSEIF ((jde-jte .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. 1)) THEN
         ! For all z's:
         DO k = kts,kte

            ! For all x's:
            DO i = its,ite

               ! For all perturbed y's:
               DO j = MAX(jts,jde-mom_pert_gps+1),jte

                  ! Add perturbation where needed:
                  field_tend(i,k,j) = field_tend(i,k,j) +                &
                                      temp_mom_pert(i,k,j)
                  mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +              &
                                       temp_mom_pert(i,k,j)

               ENDDO
            ENDDO
         ENDDO

!------ If the current memory size includes the perturbed region in +x:

      ELSEIF (ide-ite .LT. mom_pert_gps-1 ) THEN

         ! If +x is perturbed:
         IF (m_pert_locx .EQ. 1) THEN

            ! For all z's:
            DO k = kts,kte

               ! For all perturbed x's:
               DO i = MAX(its,ide-mom_pert_gps+1),ite

                  ! For all y's:
                  DO j = jts,jte

                     ! Add perturbation where needed:
                     field_tend(i,k,j) = field_tend(i,k,j) +                &
                                         temp_mom_pert(i,k,j)
                     mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +              &
                                          temp_mom_pert(i,k,j)

                  ENDDO
               ENDDO

               ! If the current memory size ALSO includes the perturbed region in -y:
               IF ((jts-jds .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. -1)) THEN

                  ! For all x's:
                  DO i = its,ite

                     ! For all perturbed y's:
                     DO j = jts,MIN(jte,jds+mom_pert_gps-1)

                        ! Add perturbation where needed:
                        field_tend(i,k,j) = field_tend(i,k,j) +             &
                                            temp_mom_pert(i,k,j)
                        mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                             temp_mom_pert(i,k,j)

                     ENDDO
                  ENDDO
               ENDIF

               ! If the current memory size ALSO includes the perturbed region in +y:
               IF ((jde-jte .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. 1)) THEN

                  ! For all x's:
                  DO i = its,ite

                     ! For all perturbed y's:
                     DO j = MAX(jts,jde-mom_pert_gps+1),jte

                        ! Add perturbation where needed:
                        field_tend(i,k,j) = field_tend(i,k,j) +             &
                                            temp_mom_pert(i,k,j)
                        mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                             temp_mom_pert(i,k,j)

                     ENDDO
                  ENDDO
               ENDIF
            ENDDO

         ! If +x is NOT perturbed:
         ELSEIF (m_pert_locx .NE. 1) THEN

            ! If the current memory size ALSO includes the perturbed region in -y:
            IF ((jts-jds .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. -1)) THEN

               ! For all z's:
               DO k = kts,kte

                  ! For all x's (ignoring the perturbed x-region, since it was
                  ! already perturbed):
                  DO i = its,ide-mom_pert_gps+1

                     ! For all perturbed y's:
                     DO j = jts,MIN(jte,jds+mom_pert_gps-1)

                        ! Add perturbation where needed:
                        field_tend(i,k,j) = field_tend(i,k,j) +             &
                                            temp_mom_pert(i,k,j)
                        mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                             temp_mom_pert(i,k,j)

                     ENDDO
                  ENDDO
               ENDDO

            ! If the current memory size ALSO includes the perturbed region in +y:
            ELSEIF ((jde-jte .LT. mom_pert_gps-1) .AND. (m_pert_locy .EQ. 1)) THEN

              ! For all z's:
               DO k = kts,kte

                  ! For all x's (ignoring the perturbed x-region, since it was
                  ! already perturbed):
                  DO i = its,ide-mom_pert_gps+1

                     ! For all perturbed y's:
                     DO j = MAX(jts,jde-mom_pert_gps+1),jte

                        ! Add perturbation where needed:
                        field_tend(i,k,j) = field_tend(i,k,j) +             &
                                            temp_mom_pert(i,k,j)
                        mom_pert_uv(i,k,j) = mom_pert_uv(i,k,j) +           &
                                             temp_mom_pert(i,k,j)

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF


   END SUBROUTINE mom_pert

!------------------------------- END LJMB -------------------------------
!------------------------------------------------------------------------


!------------------------------------------------------------------------
   SUBROUTINE spec_bdy_chem_perturb (periodic_x,                &
                                     field_bdy_tend_xs, field_bdy_tend_xe, &
                                     field_bdy_tend_ys, field_bdy_tend_ye, &
                                     field_scalar_perturb, &
                                     spec_bdy_width, spec_zone,  &
                                     kme_stoch,                  & ! stoch  dims
                                     ids,ide, jds,jde, kds,kde,  & ! domain dims
      ims,ime, jms,jme, kms,kme,  & ! memory dims
                                     ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                     its, ite, jts, jte, kts, kte)
   IMPLICIT NONE


   LOGICAL ,               INTENT(IN   ) :: periodic_x
   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte, &
                                            kme_stoch
   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, spec_zone

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),        INTENT(IN   ) :: field_scalar_perturb
   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ), INTENT(INOUT) :: field_bdy_tend_xs, field_bdy_tend_xe
   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ), INTENT(INOUT) :: field_bdy_tend_ys, field_bdy_tend_ye


         CALL spec_bdytend_perturb_chem ( field_bdy_tend_xs, field_bdy_tend_xe, &
                                          field_bdy_tend_ys, field_bdy_tend_ye, &
                                          field_scalar_perturb, 'c',       &
                                          periodic_x,                 &
                                          spec_bdy_width, spec_zone,  &
                                          kme_stoch,                  &
                                          ids,ide, jds,jde, kds,kde,  & ! domain dims
                                          ims,ime, jms,jme, kms,kme,  & ! memory dims
                                          ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                                          its,ite, jts,jte, kts,kte )

   END SUBROUTINE spec_bdy_chem_perturb   



!------------------------------------------------------------------------
   SUBROUTINE spec_bdy_scalar ( scalar_tend,    &
                          scalar_bxs,scalar_bxe,scalar_bys,scalar_bye,  &
                          scalar_btxs,scalar_btxe,scalar_btys,scalar_btye, &
                          spec_bdy_width, spec_zone,                   &
                          config_flags,               &
                          ids,ide, jds,jde, kds,kde,  & ! domain dims
                          ims,ime, jms,jme, kms,kme,  & ! memory dims
                          ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                          its, ite, jts, jte, kts, kte)
   IMPLICIT NONE

   !  Input data.
   TYPE( grid_config_rec_type ) config_flags


   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte
   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, spec_zone

   REAL , DIMENSION( ims:ime , kms:kme, jms:jme  ) , INTENT(OUT  ) :: scalar_tend

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: scalar_bxs,scalar_bxe, &
                                                                               scalar_btxs,scalar_btxe

   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width    ), INTENT(IN   ) :: scalar_bys,scalar_bye, &
                                                                               scalar_btys,scalar_btye

!Local
   INTEGER :: i,j,k


         CALL spec_bdytend (   scalar_tend,                &
                               scalar_bxs,scalar_bxe,scalar_bys,scalar_bye, scalar_btxs,scalar_btxe,scalar_btys,scalar_btye,    &
                               'q'     , config_flags, &
                               spec_bdy_width, spec_zone, &
                               ids,ide, jds,jde, kds,kde,  & ! domain dims
                               ims,ime, jms,jme, kms,kme,  & ! memory dims
                               ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                               its,ite, jts,jte, kts,kte )


   END SUBROUTINE spec_bdy_scalar

!------------------------------------------------------------------------

   SUBROUTINE set_phys_bc_dry_1( config_flags, u_1, u_2, v_1, v_2,   &
                                 rw_1, rw_2, w_1, w_2,           &
                                 t_1, t_2, tp_1, tp_2, pp, pip,  &
                                 ids,ide, jds,jde, kds,kde,      &
                                 ims,ime, jms,jme, kms,kme,      &
                                 ips,ipe, jps,jpe, kps,kpe,      &
                                 its,ite, jts,jte, kts,kte      )

!
!  this is just a wrapper to call the boundary condition routines
!  for each variable
!

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      TYPE( grid_config_rec_type ) config_flags

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &
           u_1,u_2, v_1, v_2, rw_1, rw_2, w_1, w_2,                  &
           t_1, t_2, tp_1, tp_2, pp, pip



      CALL set_physical_bc3d( u_1  , 'u', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( u_2  , 'u', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( v_1  , 'v', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( v_2  , 'v', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( rw_1 , 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( rw_2 , 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( w_1  , 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( w_2  , 'w', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( t_1, 'p', config_flags,                 &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( t_2, 'p', config_flags,                 &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( tp_1, 'p', config_flags,                &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( tp_2, 'p', config_flags,                &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( pp , 'p', config_flags,                 &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( pip , 'p', config_flags,                &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )

  END SUBROUTINE set_phys_bc_dry_1

!--------------------------------------------------------------

   SUBROUTINE set_phys_bc_dry_2( config_flags,                     &
                                 u_1, u_2, v_1, v_2, w_1, w_2,     &
                                 t_1, t_2, ph_1, ph_2, mu_1, mu_2, &
                                 ids,ide, jds,jde, kds,kde,        &
                                 ims,ime, jms,jme, kms,kme,        &
                                 ips,ipe, jps,jpe, kps,kpe,        &
                                 its,ite, jts,jte, kts,kte        )

!
!  this is just a wrapper to call the boundary condition routines
!  for each variable
!

      IMPLICIT NONE

      TYPE( grid_config_rec_type ) config_flags

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &
         u_1, u_2, v_1, v_2, w_1, w_2,                       &
         t_1, t_2, ph_1, ph_2

      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: &
                             mu_1, mu_2


      CALL set_physical_bc3d( u_1, 'U', config_flags,           &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( u_2, 'U', config_flags,           &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( v_1 , 'V', config_flags,          &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( v_2 , 'V', config_flags,          &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( w_1, 'w', config_flags,           &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( w_2, 'w', config_flags,           &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( t_1, 'p', config_flags,           &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( t_2, 'p', config_flags,           &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( ph_1 , 'w', config_flags,         &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( ph_2 , 'w', config_flags,         &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc2d( mu_1, 't', config_flags, &
                              ids, ide, jds, jde,  &
                              ims, ime, jms, jme,  &
                              ips, ipe, jps, jpe,  &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( mu_2, 't', config_flags, &
                              ids, ide, jds, jde,  &
                              ims, ime, jms, jme,  &
                              ips, ipe, jps, jpe,  &
                              its, ite, jts, jte  )

   END SUBROUTINE set_phys_bc_dry_2

!------------------------------------------------------------------------

   SUBROUTINE set_phys_bc_smallstep_1( config_flags, ru_1, du, rv_1, dv,   &
                                       ids,ide, jds,jde, kds,kde,      &
                                       ims,ime, jms,jme, kms,kme,      &
                                       ips,ipe, jps,jpe, kps,kpe,      &
                                       its,ite, jts,jte, kts,kte      )

!
!  this is just a wrapper to call the boundary condition routines
!  for each variable
!

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      TYPE( grid_config_rec_type ) config_flags

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &
           ru_1,du, rv_1, dv

      CALL set_physical_bc3d( ru_1  , 'u', config_flags,              &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kde )
      CALL set_physical_bc3d( du , 'u', config_flags,                 &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kde )
      CALL set_physical_bc3d( rv_1  , 'v', config_flags,              &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kde )
      CALL set_physical_bc3d( dv  , 'v', config_flags,                &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kde )

  END SUBROUTINE set_phys_bc_smallstep_1

!-------------------------------------------------------------------

   SUBROUTINE rk_phys_bc_dry_1( config_flags, u, v, rw, w,  &
                                muu, muv, mut, php, alt, p, &
                                ids,ide, jds,jde, kds,kde,  &
                                ims,ime, jms,jme, kms,kme,  &
                                ips,ipe, jps,jpe, kps,kpe,  &
                                its,ite, jts,jte, kts,kte  )

!
!  this is just a wrapper to call the boundary condition routines
!  for each variable
!

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      TYPE( grid_config_rec_type ) config_flags

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                    &
                                INTENT(INOUT) ::  u, v, rw, w, php, alt, p
      REAL, DIMENSION( ims:ime, jms:jme ),                             &
                                INTENT(INOUT) ::    muu, muv, mut

      CALL set_physical_bc3d( u  , 'u', config_flags,             &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( v  , 'v', config_flags,             &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d(rw , 'w', config_flags,              &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( w , 'w', config_flags,              &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( php , 'w', config_flags,            &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( alt, 't', config_flags,             &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc3d( p, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc2d( muu, 'u', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              ips, ipe, jps, jpe,      &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( muv, 'v', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              ips, ipe, jps, jpe,      &
                              its, ite, jts, jte  )

      CALL set_physical_bc2d( mut, 't', config_flags,  &
                              ids, ide, jds, jde,      &
                              ims, ime, jms, jme,      &
                              ips, ipe, jps, jpe,      &
                              its, ite, jts, jte  )

  END SUBROUTINE rk_phys_bc_dry_1

!------------------------------------------------------------------------

  SUBROUTINE rk_phys_bc_dry_2( config_flags, u, v, w,      &
                               t, ph, mu,                  &
                               ids,ide, jds,jde, kds,kde,  &
                               ims,ime, jms,jme, kms,kme,  &
                               ips,ipe, jps,jpe, kps,kpe,  &
                               its,ite, jts,jte, kts,kte  )

!
!  this is just a wrapper to call the boundary condition routines
!  for each variable
!

      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

      TYPE( grid_config_rec_type ) config_flags

      REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &
                             u, v, w, t, ph

      REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: &
                             mu

      CALL set_physical_bc3d( u   , 'U', config_flags,            &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( v   , 'V', config_flags,            &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( w  , 'w', config_flags,             &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( t, 'p', config_flags,               &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )
      CALL set_physical_bc3d( ph  , 'w', config_flags,            &
                              ids, ide, jds, jde, kds, kde,       &
                              ims, ime, jms, jme, kms, kme,       &
                              ips, ipe, jps, jpe, kps, kpe,       &
                              its, ite, jts, jte, kts, kte )

      CALL set_physical_bc2d( mu, 't', config_flags, &
                              ids, ide, jds, jde,    &
                              ims, ime, jms, jme,    &
                              ips, ipe, jps, jpe,    &
                              its, ite, jts, jte    )

  END SUBROUTINE rk_phys_bc_dry_2

!---------------------------------------------------------------------

   SUBROUTINE zero_bdytend  (                                                  &
                              u_btxs,u_btxe,u_btys,u_btye,                     &
                              v_btxs,v_btxe,v_btys,v_btye,                     &
                              ph_btxs,ph_btxe,ph_btys,ph_btye,                 &
                              t_btxs,t_btxe,t_btys,t_btye,                     &
                              w_btxs,w_btxe,w_btys,w_btye,                     &
                              mu_btxs,mu_btxe,mu_btys,mu_btye,                 &
                              moist_btxs,moist_btxe,   &
                              moist_btys,moist_btye,   &
                              scalar_btxs,scalar_btxe,   &
                              scalar_btys,scalar_btye,   &
                              spec_bdy_width,n_moist,n_scalar,                 &
                              ids,ide, jds,jde, kds,kde,  & ! domain dims
                              ims,ime, jms,jme, kms,kme,  & ! memory dims
                              ips,ipe, jps,jpe, kps,kpe,  & ! patch  dims
                              its,ite, jts,jte, kts,kte   )
   IMPLICIT NONE

   !  Input data.

   INTEGER ,               INTENT(IN   ) :: spec_bdy_width, n_moist,n_scalar

   INTEGER ,               INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            ips, ipe, jps, jpe, kps, kpe, &
                                            its, ite, jts, jte, kts, kte

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width    ), INTENT(INOUT) :: u_btxs,u_btxe, &
                                                                               v_btxs,v_btxe, &
                                                                               ph_btxs,ph_btxe, &
                                                                               w_btxs,w_btxe, &
                                                                               t_btxs,t_btxe

   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width    ), INTENT(INOUT) :: u_btys,u_btye, &
                                                                               v_btys,v_btye, &
                                                                               ph_btys,ph_btye, &
                                                                               w_btys,w_btye, &
                                                                               t_btys,t_btye

   REAL,  DIMENSION( jms:jme , 1:1     , spec_bdy_width    ), INTENT(INOUT) :: mu_btxs,mu_btxe
   REAL,  DIMENSION( ims:ime , 1:1     , spec_bdy_width    ), INTENT(INOUT) :: mu_btys,mu_btye

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width , n_moist ), INTENT(INOUT) :: moist_btxs,moist_btxe
   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width , n_moist ), INTENT(INOUT) :: moist_btys,moist_btye

   REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width , n_scalar ), INTENT(INOUT) :: scalar_btxs,scalar_btxe
   REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width , n_scalar ), INTENT(INOUT) :: scalar_btys,scalar_btye


! setting bdy tendencies to zero during DFI

       CALL wrf_debug( 10, 'In zero_bdytend, setting bdy tendencies to 0 during DFI' )
       u_btxs = 0.
       u_btxe = 0.
       u_btys = 0.
       u_btye = 0.
       v_btxs = 0.
       v_btxe = 0.
       v_btys = 0.
       v_btye = 0.
       t_btxs = 0.
       t_btxe = 0.
       t_btys = 0.
       t_btye = 0.
       ph_btxs = 0.
       ph_btxe = 0.
       ph_btys = 0.
       ph_btye = 0.
       mu_btxs = 0.
       mu_btxe = 0.
       mu_btys = 0.
       mu_btye = 0.
       moist_btxs = 0.
       moist_btxe = 0.
       moist_btys = 0.
       moist_btye = 0.
       scalar_btxs = 0.
       scalar_btxe = 0.
       scalar_btys = 0.
       scalar_btye = 0.


!  ENDIF

   END SUBROUTINE zero_bdytend

!---------------------------------------------------------------------

  SUBROUTINE set_w_surface( config_flags, znw, fill_w_flag,              &
                            w, ht, u, v, cf1, cf2, cf3, rdx, rdy,        &
                            msftx, msfty,                                &
                            ids, ide, jds, jde, kds, kde,                &
                            ims, ime, jms, jme, kms, kme,                &
                            its, ite, jts, jte, kts, kte                 )
     implicit none

     TYPE( grid_config_rec_type ) config_flags

     INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte

     REAL :: rdx, rdy, cf1, cf2, cf3

     REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                                 INTENT(IN   ) ::  u,       &
                                                                   v

     REAL , DIMENSION(  ims:ime , kms:kme, jms:jme ) ,                      &
                                              INTENT(INOUT) ::  w

     REAL , DIMENSION(  ims:ime , jms:jme ) , INTENT(IN   ) ::  ht,         &
                                                                msftx,      &
                                                                msfty
     REAL , DIMENSION( kms:kme ) , INTENT(IN   ) ::  znw
  
     LOGICAL, INTENT(IN   ) :: fill_w_flag


     INTEGER :: i,j,k
     INTEGER :: ip1,im1,jp1,jm1
     INTEGER :: ip1_limit,im1_limit,jp1_limit,jm1_limit

!  set kinematic lower boundary condition on W

!  Comments on directional map scale factors:
!  Chain rule: if Z=Z(X,Y) [true at the surface] then
!  dZ/dt = dZ/dX * dX/dt + dZ/dY * dY/dt, U=dX/dt, V=dY/dt
!  using capitals to denote actual values
!  in mapped values, u=U, v=V, z=Z, 1/dX=mx/dx, 1/dY=my/dy
!     => w = dz/dt = mx u dz/dx + my v dz/dy
!  [where dz/dx is just the surface height change between x
!   gridpoints, and dz/dy is the change between y gridpoints]
!  [NB - cf1, cf2 and cf3 do vertical weighting of u or v values
!   nearest the surface]

! get indices for points next to edge of domain

     jm1_limit = jds        ! No periodic BC's
     jp1_limit = jde-1
     im1_limit = ids
     ip1_limit = ide-1

     IF ( config_flags%periodic_x ) THEN
       im1_limit = ids-1
       ip1_limit = ide
     ENDIF

     IF ( config_flags%periodic_y ) THEN
       jm1_limit = jds-1
       jp1_limit = jde
     ENDIF

     DO j = jts,min(jte,jde-1)
       jm1 = max(j-1, jm1_limit)
       jp1 = min(j+1, jp1_limit)
     DO i = its,min(ite,ide-1)
       im1 = max(i-1, im1_limit)
       ip1 = min(i+1, ip1_limit)

         w(i,1,j)=  msfty(i,j)*                              &
                  .5*rdy*(                                   &
                           (ht(i,jp1)-ht(i,j  ))             &
          *(cf1*v(i,1,j+1)+cf2*v(i,2,j+1)+cf3*v(i,3,j+1))    &
                          +(ht(i,j  )-ht(i,jm1))             &
          *(cf1*v(i,1,j  )+cf2*v(i,2,j  )+cf3*v(i,3,j  ))  ) &
                    +msftx(i,j)*                             &
                  .5*rdx*(                                   &
                           (ht(ip1,j)-ht(i,j  ))             &
          *(cf1*u(i+1,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j))    &
                          +(ht(i  ,j)-ht(im1,j))             &
          *(cf1*u(i  ,1,j)+cf2*u(i  ,2,j)+cf3*u(i  ,3,j))  )

      ENDDO
      ENDDO

! Fill the atmospheric w field with smoothly decaying values
      IF (fill_w_flag) THEN
        DO j = jts,min(jte,jde-1)
        DO k = kts+1,kte
        DO i = its,min(ite,ide-1)
          w(i,k,j) = w(i,1,j)*znw(k)*znw(k)
        ENDDO
        ENDDO
        ENDDO
      ENDIF

  END SUBROUTINE set_w_surface

  SUBROUTINE lbc_fcx_gcx ( fcx , gcx , spec_bdy_width , &
                           spec_zone , relax_zone , dt , spec_exp , &
                           specified , nested )
  
     IMPLICIT NONE
  
     INTEGER , INTENT(IN) :: spec_bdy_width , spec_zone , relax_zone
     REAL , INTENT(IN) :: dt , spec_exp
     LOGICAL , INTENT(IN) :: specified , nested
     REAL , DIMENSION(spec_bdy_width) :: fcx , gcx
  
     ! Local variables.
  
     INTEGER :: loop
     REAL :: spongeweight
  
     IF (specified) THEN
       
       ! Arrays for specified boundary conditions
       
       DO loop = spec_zone + 1, spec_zone + relax_zone
         fcx(loop) = 0.1 / dt * (spec_zone + relax_zone - loop) / (relax_zone - 1)
         gcx(loop) = 1.0 / dt / 50. * (spec_zone + relax_zone - loop) / (relax_zone - 1)
         spongeweight=exp(-(loop-(spec_zone + 1))*spec_exp)
         fcx(loop) = fcx(loop)*spongeweight
         gcx(loop) = gcx(loop)*spongeweight
       ENDDO   
       
     ELSE IF (nested) THEN
       
       ! Arrays for specified boundary conditions
       
       DO loop = spec_zone + 1, spec_zone + relax_zone
         fcx(loop) = 0.1 / dt * (spec_zone + relax_zone - loop) / (relax_zone - 1)
         gcx(loop) = 1.0 / dt / 50. * (spec_zone + relax_zone - loop) / (relax_zone - 1)
!        spongeweight=EXP(-(loop-2)/3.)
!        fcx(loop) = fcx(loop)*spongeweight
!        gcx(loop) = gcx(loop)*spongeweight
!        fcx(loop) = 0.
!        gcx(loop) = 0.
       ENDDO
       
     ENDIF
  
  END SUBROUTINE lbc_fcx_gcx

!------------------------------------------------------------------------

   SUBROUTINE theta_and_thetam_lbc_only (                              &
                                theta_to_thetam,                       &
                                mub,                                   &
                                mu_bdy_xs, mu_bdy_xe,                  &
                                mu_bdy_ys, mu_bdy_ye,                  &
                                mu_bdy_tend_xs, mu_bdy_tend_xe,        &
                                mu_bdy_tend_ys, mu_bdy_tend_ye,        &
                                orig_t_bdy_xs, orig_t_bdy_xe,          &
                                orig_t_bdy_ys, orig_t_bdy_ye,          &
                                orig_t_bdy_tend_xs, orig_t_bdy_tend_xe,&
                                orig_t_bdy_tend_ys, orig_t_bdy_tend_ye,&
                                moist_bdy_xs, moist_bdy_xe,            &
                                moist_bdy_ys, moist_bdy_ye,            &
                                moist_bdy_tend_xs, moist_bdy_tend_xe,  &
                                moist_bdy_tend_ys, moist_bdy_tend_ye,  &
                                spec_bdy_width,                        &
                                dt_interval,                           &
                                ids,ide, jds,jde, kds,kde,             &
                                ims,ime, jms,jme, kms,kme,             &
                                ips,ipe, jps,jpe, kps,kpe,             &
                                its,ite, jts,jte, kts,kte              )

      IMPLICIT NONE

      !  This routine is called from the solve_em routine.  The purpose is to
      !  convert the thermal lateral boundary conditions between dry potential
      !  temperature and moist potential temperature.  The first argument is a
      !  flag telling us the direction of the conversion:
      !     True  = convert dry to moist potential temp
      !     False = convert moist to dry potential temp

      LOGICAL,      INTENT(IN   )    :: theta_to_thetam

      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      INTEGER,      INTENT(IN   )    :: spec_bdy_width
      REAL   ,      INTENT(IN   )    :: dt_interval

      REAL,  DIMENSION( ims:ime ,           jms:jme        ), INTENT(IN   ) :: mub

      REAL,  DIMENSION( jms:jme , 1       , spec_bdy_width ), INTENT(IN   ) :: mu_bdy_xs, mu_bdy_xe
      REAL,  DIMENSION( ims:ime , 1       , spec_bdy_width ), INTENT(IN   ) :: mu_bdy_ys, mu_bdy_ye
      REAL,  DIMENSION( jms:jme , 1       , spec_bdy_width ), INTENT(IN   ) :: mu_bdy_tend_xs, mu_bdy_tend_xe
      REAL,  DIMENSION( ims:ime , 1       , spec_bdy_width ), INTENT(IN   ) :: mu_bdy_tend_ys, mu_bdy_tend_ye

      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ), INTENT(INOUT) :: orig_t_bdy_xs, orig_t_bdy_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ), INTENT(INOUT) :: orig_t_bdy_ys, orig_t_bdy_ye
      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ), INTENT(INOUT) :: orig_t_bdy_tend_xs, orig_t_bdy_tend_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ), INTENT(INOUT) :: orig_t_bdy_tend_ys, orig_t_bdy_tend_ye

      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ), INTENT(IN   ) :: moist_bdy_xs, moist_bdy_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ), INTENT(IN   ) :: moist_bdy_ys, moist_bdy_ye
      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ), INTENT(IN   ) :: moist_bdy_tend_xs, moist_bdy_tend_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ), INTENT(IN   ) :: moist_bdy_tend_ys, moist_bdy_tend_ye

      !  Local variables
#ifdef _OPENMP
      INTEGER, EXTERNAL :: omp_get_thread_num
#endif

      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ) :: t_bdy_xs, t_bdy_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ) :: t_bdy_ys, t_bdy_ye
      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ) :: t_bdy_tend_xs, t_bdy_tend_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ) :: t_bdy_tend_ys, t_bdy_tend_ye

      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ) :: new_t_bdy_xs, new_t_bdy_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ) :: new_t_bdy_ys, new_t_bdy_ye
      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width ) :: new_t_bdy_tend_xs, new_t_bdy_tend_xe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width ) :: new_t_bdy_tend_ys, new_t_bdy_tend_ye

      INTEGER    :: i, j, k, ii, jj

      REAL :: dt_time_until_next_lbc

      REAL :: mu_old_bdy_xs         , mu_old_bdy_xe
      REAL :: mu_old_bdy_ys         , mu_old_bdy_ye
      REAL :: mu_new_bdy_xs         , mu_new_bdy_xe
      REAL :: mu_new_bdy_ys         , mu_new_bdy_ye

      REAL :: t_old_bdy_xs          , t_old_bdy_xe
      REAL :: t_old_bdy_ys          , t_old_bdy_ye
      REAL :: t_new_bdy_xs          , t_new_bdy_xe
      REAL :: t_new_bdy_ys          , t_new_bdy_ye
      REAL :: t_old_bdy_tend_xs     , t_old_bdy_tend_xe
      REAL :: t_old_bdy_tend_ys     , t_old_bdy_tend_ye

      REAL :: moist_old_bdy_xs      , moist_old_bdy_xe
      REAL :: moist_old_bdy_ys      , moist_old_bdy_ye
      REAL :: moist_new_bdy_xs      , moist_new_bdy_xe
      REAL :: moist_new_bdy_ys      , moist_new_bdy_ye
      REAL :: moist_old_bdy_tend_xs , moist_old_bdy_tend_xe
      REAL :: moist_old_bdy_tend_ys , moist_old_bdy_tend_ye

      INTEGER :: i_min, i_max, j_min, j_max

      !  IF      (       theta_to_thetam ) THEN
      !     Convert dry potential temperature to theta_m
      !     Defined as: theta_m = ( theta   + T0 ) * ( 1. + (R_v/R_d) Qv ) - T0
      !  ELSE
      !     Convert dry potential temperature to theta_m
      !     Defined as: theta   = ( theta_m + T0 ) / ( 1. + (R_v/R_d) Qv ) - T0
      !  END IF

      !  We want the current value and the tendency, using information mostly
      !  from the lateral boundary file.  In that file, the thermal variable
      !  is a potential temperature with the T0 offset removed (theta-300).  Both
      !  the moisture variable and the potential temperature are coupled
      !  (multiplied by total dry column pressure).  And to add one more complication,
      !  the MU variable in the lateral boundary array is perturbation only.

      !  Since we need to end up with lateral boundary values that are coupled,
      !  we need to first DECOUPLE T and Qv, compute Tm, and then couple that.  As
      !  there is a need for the lateral tendency also, we compute the T and Qv
      !  values at the two boundary times (previous/current and next).  These two
      !  times are adequate to get us a tendency.  For the tendency, we need to have
      !  coupled values for the T (or Tm) at both times, which gives us a coupled
      !  tendency.  We cannot have an uncoupled tendency and somehow multiply that
      !  by some intermediate/average column pressure.

      !  This routine's purpose is to manufacture a lateral boundary set of arrays
      !  (all eight of them) for the thermal field.  Depending on the logical flag
      !  passed in, this will either be dry potential temperature or moist potential
      !  temperature.

      !  The i_min, i_max for the south and north boundaries depends on if we are doing
      !  serial, OpenMP, or MPI.  For OpenMP, we do not want any overlap between tiles that
      !  are on the same task (either OpenMP only, or OpenMP+MPI).

      IF      ( its .EQ. ids ) THEN
         i_min = its
      ELSE IF ( its .EQ. ips ) THEN
         i_min = ims
      ELSE
         i_min = its
      END IF
      i_min = MAX(ids,i_min)

      IF      ( ite .EQ. ide ) THEN
         i_max = ite
      ELSE IF ( ite .EQ. ipe ) THEN
         i_max = ime
      ELSE
         i_max = ite
      END IF
      i_max = MIN(i_max,ide-1)

      !  South and north lateral boundaries.  This is the i-extent of its through ite, but j only
      !  goes to within spec_bdy_width of the top and bottom (north and south) boundaries.

      !  South boundary: i,k,j
      !  jj increasing

      DO jj = MAX(jts,1) , MIN(jte,jde-1,spec_bdy_width)
         j = jj
         DO k = kts , kte-1
            DO i = i_min , i_max
              t_bdy_ys     (i,k,j) = orig_t_bdy_ys     (i,k,j)
              t_bdy_tend_ys(i,k,j) = orig_t_bdy_tend_ys(i,k,j)
            END DO
         END DO
      END DO
      DO jj = MAX(jts,1) , MIN(jte,jde-1,spec_bdy_width)
         j = jj
         DO k = kts , kte-1
            DO i = i_min , i_max
               mu_old_bdy_ys         =   mu_bdy_ys(i,1,j)    + mub(i,jj)
               t_old_bdy_ys          = ( t_bdy_ys(i,k,j)                                            ) / mu_old_bdy_ys
               moist_old_bdy_ys      = ( moist_bdy_ys(i,k,j)                                        ) / mu_old_bdy_ys
               mu_new_bdy_ys         =   mu_old_bdy_ys       + mu_bdy_tend_ys(i,1,j)   *dt_interval
               t_new_bdy_ys          = ( t_bdy_ys(i,k,j)     + t_bdy_tend_ys(i,k,j)    *dt_interval ) / mu_new_bdy_ys
               moist_new_bdy_ys      = ( moist_bdy_ys(i,k,j) + moist_bdy_tend_ys(i,k,j)*dt_interval ) / mu_new_bdy_ys
               t_old_bdy_tend_ys     = ( t_new_bdy_ys        - t_old_bdy_ys                         ) / dt_interval
               moist_old_bdy_tend_ys = ( moist_new_bdy_ys    - moist_old_bdy_ys                     ) / dt_interval
               IF ( theta_to_thetam ) THEN
                  new_t_bdy_ys(i,k,j) = ( ( ( t_old_bdy_ys + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_ys ) ) - T0 ) * mu_old_bdy_ys
                  new_t_bdy_tend_ys(i,k,j) = ( ( mu_new_bdy_ys * ( ( t_new_bdy_ys + T0 ) * ( 1. + (R_v/R_d) * moist_new_bdy_ys ) - T0 ) ) - &
                                           ( mu_old_bdy_ys * ( ( t_old_bdy_ys + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_ys ) - T0 ) ) ) / dt_interval
               ELSE
                  new_t_bdy_ys(i,k,j) = ( ( ( t_old_bdy_ys + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_ys ) ) - T0 ) * mu_old_bdy_ys
                  new_t_bdy_tend_ys(i,k,j) = ( ( mu_new_bdy_ys * ( ( t_new_bdy_ys + T0 ) / ( 1. + (R_v/R_d) * moist_new_bdy_ys ) - T0 ) ) - &
                                           ( mu_old_bdy_ys * ( ( t_old_bdy_ys + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_ys ) - T0 ) ) ) / dt_interval
               END IF
            END DO
         END DO
      END DO

      !  North boundary: i,k,j
      !  jj decreasing

      DO jj = MIN(jde-1,jte) , MAX(jde-spec_bdy_width,jts) , -1
         j = jde-jj
         DO k = kts , kte-1
            DO i = i_min , i_max
               t_bdy_ye     (i,k,j) = orig_t_bdy_ye     (i,k,j)
               t_bdy_tend_ye(i,k,j) = orig_t_bdy_tend_ye(i,k,j)
            END DO
         END DO
      END DO
      DO jj = MIN(jde-1,jte) , MAX(jde-spec_bdy_width,jts) , -1
         j = jde-jj
         DO k = kts , kte-1
            DO i = i_min , i_max
               mu_old_bdy_ye         =   mu_bdy_ye(i,1,j)    + mub(i,jj)
               t_old_bdy_ye          = ( t_bdy_ye(i,k,j)                                            ) / mu_old_bdy_ye
               moist_old_bdy_ye      = ( moist_bdy_ye(i,k,j)                                        ) / mu_old_bdy_ye
               mu_new_bdy_ye         =   mu_old_bdy_ye       + mu_bdy_tend_ye(i,1,j)   *dt_interval
               t_new_bdy_ye          = ( t_bdy_ye(i,k,j)     + t_bdy_tend_ye(i,k,j)    *dt_interval ) / mu_new_bdy_ye
               moist_new_bdy_ye      = ( moist_bdy_ye(i,k,j) + moist_bdy_tend_ye(i,k,j)*dt_interval ) / mu_new_bdy_ye
               t_old_bdy_tend_ye     = ( t_new_bdy_ye        - t_old_bdy_ye                         ) / dt_interval
               moist_old_bdy_tend_ye = ( moist_new_bdy_ye    - moist_old_bdy_ye                     ) / dt_interval
               IF ( theta_to_thetam ) THEN
                  new_t_bdy_ye(i,k,j) = ( ( ( t_old_bdy_ye + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_ye ) ) - T0 ) * mu_old_bdy_ye
                  new_t_bdy_tend_ye(i,k,j) = ( ( mu_new_bdy_ye * ( ( t_new_bdy_ye + T0 ) * ( 1. + (R_v/R_d) * moist_new_bdy_ye ) - T0 ) ) - &
                                           ( mu_old_bdy_ye * ( ( t_old_bdy_ye + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_ye ) - T0 ) ) ) / dt_interval
               ELSE
                  new_t_bdy_ye(i,k,j) = ( ( ( t_old_bdy_ye + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_ye ) ) - T0 ) * mu_old_bdy_ye
                  new_t_bdy_tend_ye(i,k,j) = ( ( mu_new_bdy_ye * ( ( t_new_bdy_ye + T0 ) / ( 1. + (R_v/R_d) * moist_new_bdy_ye ) - T0 ) ) - &
                                           ( mu_old_bdy_ye * ( ( t_old_bdy_ye + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_ye ) - T0 ) ) ) / dt_interval
               END IF
            END DO
         END DO
      END DO

      !  The j_min, j_max for the west and east boundaries depends on if we are doing
      !  serial, OpenMP, or MPI.  For OpenMP, we do not want any overlap between tiles that
      !  are on the same task (either OpenMP only, or OpenMP+MPI).

      IF      ( jts .EQ. jds ) THEN
         j_min = jts
      ELSE IF ( jts .EQ. jps ) THEN
         j_min = jms
      ELSE
         j_min = jts
      END IF
      j_min = MAX(jds,j_min)

      IF      ( jte .EQ. jde ) THEN
         j_max = jte
      ELSE IF ( jte .EQ. jpe ) THEN
         j_max = jme
      ELSE
         j_max = jte
      END IF
      j_max = MIN(j_max,jde-1)

      !  West and east lateral boundaries.  This is the j-extent of jts through jte, but i only
      !  goes to within spec_bdy_width of the left and right (west and east) boundaries.

      !  West boundary: j,k,i
      !  ii increasing

      DO ii = MAX(its,1) , MIN(ite,ide-1,spec_bdy_width)
         i = ii
         DO k = kts , kte-1
            DO j = j_min , j_max
               t_bdy_xs     (j,k,i) = orig_t_bdy_xs     (j,k,i)
               t_bdy_tend_xs(j,k,i) = orig_t_bdy_tend_xs(j,k,i)
            END DO
         END DO
      END DO
      DO ii = MAX(its,1) , MIN(ite,ide-1,spec_bdy_width)
         i = ii
         DO k = kts , kte-1
            DO j = j_min , j_max
               mu_old_bdy_xs         =   mu_bdy_xs(j,1,i)    + mub(ii,j)
               t_old_bdy_xs          = ( t_bdy_xs(j,k,i)                                            ) / mu_old_bdy_xs
               moist_old_bdy_xs      = ( moist_bdy_xs(j,k,i)                                        ) / mu_old_bdy_xs
               mu_new_bdy_xs         =   mu_old_bdy_xs       + mu_bdy_tend_xs(j,1,i)   *dt_interval
               t_new_bdy_xs          = ( t_bdy_xs(j,k,i)     + t_bdy_tend_xs(j,k,i)    *dt_interval ) / mu_new_bdy_xs
               moist_new_bdy_xs      = ( moist_bdy_xs(j,k,i) + moist_bdy_tend_xs(j,k,i)*dt_interval ) / mu_new_bdy_xs
               t_old_bdy_tend_xs     = ( t_new_bdy_xs        - t_old_bdy_xs                         ) / dt_interval
               moist_old_bdy_tend_xs = ( moist_new_bdy_xs    - moist_old_bdy_xs                     ) / dt_interval
               IF ( theta_to_thetam ) THEN
                  new_t_bdy_xs(j,k,i) = ( ( ( t_old_bdy_xs + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_xs ) ) - T0 ) * mu_old_bdy_xs
                  new_t_bdy_tend_xs(j,k,i) = ( ( mu_new_bdy_xs * ( ( t_new_bdy_xs + T0 ) * ( 1. + (R_v/R_d) * moist_new_bdy_xs ) - T0 ) ) - &
                                           ( mu_old_bdy_xs * ( ( t_old_bdy_xs + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_xs ) - T0 ) ) ) / dt_interval
               ELSE
                  new_t_bdy_xs(j,k,i) = ( ( ( t_old_bdy_xs + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_xs ) ) - T0 ) * mu_old_bdy_xs
                  new_t_bdy_tend_xs(j,k,i) = ( ( mu_new_bdy_xs * ( ( t_new_bdy_xs + T0 ) / ( 1. + (R_v/R_d) * moist_new_bdy_xs ) - T0 ) ) - &
                                           ( mu_old_bdy_xs * ( ( t_old_bdy_xs + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_xs ) - T0 ) ) ) / dt_interval
               END IF
            END DO
         END DO
      END DO

      !  East boundary: j,k,i
      !  ii decreasing

      DO ii = MIN(ide-1,ite) , MAX(ide-spec_bdy_width,its) , -1
         i = ide-ii
         DO k = kts , kte-1
            DO j = j_min , j_max
               t_bdy_xe     (j,k,i) = orig_t_bdy_xe     (j,k,i)
               t_bdy_tend_xe(j,k,i) = orig_t_bdy_tend_xe(j,k,i)
            END DO
         END DO
      END DO
      DO ii = MIN(ide-1,ite) , MAX(ide-spec_bdy_width,its) , -1
         i = ide-ii
         DO k = kts , kte-1
            DO j = j_min , j_max
               mu_old_bdy_xe         =   mu_bdy_xe(j,1,i)    + mub(ii,j)
               t_old_bdy_xe          = ( t_bdy_xe(j,k,i)                                            ) / mu_old_bdy_xe
               moist_old_bdy_xe      = ( moist_bdy_xe(j,k,i)                                        ) / mu_old_bdy_xe
               mu_new_bdy_xe         =   mu_old_bdy_xe       + mu_bdy_tend_xe(j,1,i)   *dt_interval
               t_new_bdy_xe          = ( t_bdy_xe(j,k,i)     + t_bdy_tend_xe(j,k,i)    *dt_interval ) / mu_new_bdy_xe
               moist_new_bdy_xe      = ( moist_bdy_xe(j,k,i) + moist_bdy_tend_xe(j,k,i)*dt_interval ) / mu_new_bdy_xe
               t_old_bdy_tend_xe     = ( t_new_bdy_xe        - t_old_bdy_xe                         ) / dt_interval
               moist_old_bdy_tend_xe = ( moist_new_bdy_xe    - moist_old_bdy_xe                     ) / dt_interval
               IF ( theta_to_thetam ) THEN
                  new_t_bdy_xe(j,k,i) = ( ( ( t_old_bdy_xe + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_xe ) ) - T0 ) * mu_old_bdy_xe
                  new_t_bdy_tend_xe(j,k,i) = ( ( mu_new_bdy_xe * ( ( t_new_bdy_xe + T0 ) * ( 1. + (R_v/R_d) * moist_new_bdy_xe ) - T0 ) ) - &
                                           ( mu_old_bdy_xe * ( ( t_old_bdy_xe + T0 ) * ( 1. + (R_v/R_d) * moist_old_bdy_xe ) - T0 ) ) ) / dt_interval
               ELSE
                  new_t_bdy_xe(j,k,i) = ( ( ( t_old_bdy_xe + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_xe ) ) - T0 ) * mu_old_bdy_xe
                  new_t_bdy_tend_xe(j,k,i) = ( ( mu_new_bdy_xe * ( ( t_new_bdy_xe + T0 ) / ( 1. + (R_v/R_d) * moist_new_bdy_xe ) - T0 ) ) - &
                                           ( mu_old_bdy_xe * ( ( t_old_bdy_xe + T0 ) / ( 1. + (R_v/R_d) * moist_old_bdy_xe ) - T0 ) ) ) / dt_interval
               END IF
            END DO
         END DO
      END DO

      !  Put the final values for the tendencies into the arrays that get passed
      !  back out to the calling routine.

      DO jj = MAX(jts,1) , MIN(jte,jde-1,spec_bdy_width)
         j = jj
         DO k = kts , kte-1
            DO i = i_min , i_max
              orig_t_bdy_ys     (i,k,j) = new_t_bdy_ys     (i,k,j)
              orig_t_bdy_tend_ys(i,k,j) = new_t_bdy_tend_ys(i,k,j)
            END DO
         END DO
      END DO

      DO jj = MIN(jde-1,jte) , MAX(jde-spec_bdy_width,jts) , -1
         j = jde-jj
         DO k = kts , kte-1
            DO i = i_min , i_max
               orig_t_bdy_ye     (i,k,j) = new_t_bdy_ye     (i,k,j)
               orig_t_bdy_tend_ye(i,k,j) = new_t_bdy_tend_ye(i,k,j)
            END DO
         END DO
      END DO

      DO ii = its , MIN(ite,ide-1,spec_bdy_width)
         i = ii
         DO k = kts , kte-1
            DO j = j_min , j_max
               orig_t_bdy_xs     (j,k,i) = new_t_bdy_xs     (j,k,i)
               orig_t_bdy_tend_xs(j,k,i) = new_t_bdy_tend_xs(j,k,i)
            END DO
         END DO
      END DO

      DO ii = MIN(ide-1,ite) , MAX(ide-spec_bdy_width,its) , -1
         i = ide-ii
         DO k = kts , kte-1
            DO j = j_min , j_max
               orig_t_bdy_xe     (j,k,i) = new_t_bdy_xe     (j,k,i)
               orig_t_bdy_tend_xe(j,k,i) = new_t_bdy_tend_xe(j,k,i)
            END DO
         END DO
      END DO

   END SUBROUTINE theta_and_thetam_lbc_only

!------------------------------------------------------------------------

   SUBROUTINE mass_weight  ( field , mut, rfield , c1 , c2 , &
                             ids,ide, jds,jde, kds,kde,      &       ! domain dims
                             ims,ime, jms,jme, kms,kme,      &       ! memory dims
                             irs,ire, jrs,jre, krs,kre,      &       ! rfield dims
                             its,ite, jts,jte, kts,kte )             ! tile dims

      IMPLICIT NONE

      INTEGER , INTENT(IN   ) :: ids,ide, jds,jde, kds,kde, &
                                 ims,ime, jms,jme, kms,kme, &
                                 irs,ire, jrs,jre, krs,kre, &
                                 its,ite, jts,jte, kts,kte
      REAL , DIMENSION(ims:ime, kms:kme, jms:jme) , INTENT(IN   ) :: field
      REAL , DIMENSION(ims:ime,          jms:jme) , INTENT(IN   ) :: mut
      REAL , DIMENSION(         kms:kme         ) , INTENT(IN   ) :: c1, c2
      REAL , DIMENSION(irs:ire, krs:kre, jrs:jre) , INTENT(  OUT) :: rfield

      !  Local loop counters

      INTEGER :: i , j , k

      DO j = jts , jte
         DO k = kts , kte
            DO i = its , ite
               rfield(i,k,j) = field(i,k,j) * (c1(k)*mut(i,j)+c2(k))
            END DO
         END DO
      END DO

   END SUBROUTINE mass_weight
  
END MODULE module_bc_em
