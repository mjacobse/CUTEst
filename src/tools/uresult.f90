!-*-*-*-*-*-*-  C U T E S T    U R E S U L T    S U B R O U T I N E  -*-*-*-*-*-

      SUBROUTINE CUTEST_uresult( status, output, n, X, X_l, X_u )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( OUT ) :: status
      INTEGER, INTENT( IN )  :: output, n
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X, X_l, X_u

!  local variables

      INTEGER ( KIND = 8 ) :: sysclock_count, sysclock_count_rate
      REAL :: time_now, time_user
      REAL ( KIND = wp ) :: time_real
      REAL ( KIND = wp ) :: f, box_complementarity, box_dual_feasibility
      REAL ( KIND = wp ), DIMENSION( n ) :: Dl, Z
      REAL ( KIND = wp ) :: nc2of, nc2og, nc2oh, nhvpr
      CHARACTER ( LEN=21 ), PARAMETER :: fmt = "(A, SS, ES25.17E3, A)"

      CALL CPU_TIME( time_now )
      CALL SYSTEM_CLOCK( sysclock_count, sysclock_count_rate )
      time_user = time_now - CUTEST_data_global%sttime
      time_real = (sysclock_count - CUTEST_data_global%sttime_sysclock_count) /&
                  (1.0_wp * sysclock_count_rate)
      WRITE( output, "(A)" ) "{"
      WRITE( output, fmt ) '"time_user":', time_user, ","
      WRITE( output, fmt ) '"time_real":', time_real, ","

      nc2of = 0.0_wp
      nc2og = 0.0_wp
      nc2oh = 0.0_wp
      nhvpr = 0.0_wp
      DO i = 1, SIZE( CUTEST_work_global )
        nc2of = nc2of + CUTEST_work_global( i )%nc2of
        nc2og = nc2og + CUTEST_work_global( i )%nc2og
        nc2oh = nc2oh + CUTEST_work_global( i )%nc2oh
        nhvpr = nhvpr + CUTEST_work_global( i )%nhvpr
      END DO
      WRITE( output, fmt ) '"obj_function_evals":', nc2of, ","
      WRITE( output, fmt ) '"obj_gradient_evals":', nc2og, ","
      WRITE( output, fmt ) '"obj_hessian_evals":', nc2oh, ","
      WRITE( output, fmt ) '"hessian_vector_products":', nhvpr, ","
      WRITE( output, fmt ) '"constr_function_evals":', 0.0_wp, ","
      WRITE( output, fmt ) '"constr_gradient_evals":', 0.0_wp, ","
      WRITE( output, fmt ) '"constr_hessian_evals":', 0.0_wp, ","

      CALL CUTEST_ufn( status, n, X, f )
      WRITE( output, fmt ) '"objective_value":', f, ","
      WRITE( output, fmt ) '"primal_feasibility":', MAX( 0.0_wp,               &
        MAXVAL( X_l - X ), MAXVAL( X - X_u ) ), ","

      CALL CUTEST_ugr( status, n, X, Dl )
      Z = CUTEST_box_multipliers( X, X_l, X_u, Dl )
      Dl = Dl + Z
      WRITE( output, fmt ) '"stationarity":', MAXVAL( ABS( Dl ) ), ","

      box_dual_feasibility = 0.0_wp
      box_complementarity = 0.0_wp
      DO i = 1, n
        IF ( X_u( i ) == X_l( i ) ) THEN
          CYCLE
        END IF
        IF ( ABS( X_u( i ) - X( i ) ) < ABS( X( i ) - X_l( i ) ) ) THEN
          box_dual_feasibility = MAX( box_dual_feasibility, -Z( i ) )
          box_complementarity = MAX( box_complementarity,                      &
                                     ( X( i ) - X_u( i ) ) * Z( i ) )
        ElSE
          box_dual_feasibility = MAX( box_dual_feasibility, Z( i ) )
          box_complementarity = MAX( box_complementarity,                      &
                                     ( X_l( i ) - X( i ) ) * Z( i ) )
        END IF
      END DO
      WRITE( output, fmt ) '"box_dual_feasibility":', box_dual_feasibility, ","
      WRITE( output, fmt ) '"box_complementarity":', box_complementarity, ","

      WRITE( output, fmt ) '"general_dual_feasibility":', 0.0_wp, ","
      WRITE( output, fmt ) '"general_complementarity":', 0.0_wp, ","
      IF ( n > 1 ) THEN
        WRITE( output, fmt ) '"X":[', X( 1 ), ","
        DO i = 2, (n - 1)
          WRITE( output, "(SS, ES25.17E3, A)" ) X( i ), ","
        END DO
        WRITE( output, "(SS, ES25.17E3, A)" ) X( n ), "]"
      ELSE IF ( n > 0 ) THEN
        WRITE( output, fmt ) '"X":[', X( 1 ), "]"
      END IF

      WRITE( output, "(A)" ) "}"
      status = 0
      RETURN

!  End of subroutine CUTEST_uresult

      END SUBROUTINE CUTEST_uresult
