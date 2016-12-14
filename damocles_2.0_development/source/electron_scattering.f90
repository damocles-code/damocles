MODULE electron_scattering

    use input

    IMPLICIT NONE

contains

    SUBROUTINE N_e_const(ES_const)
        REAL, INTENT(OUT)    :: ES_const
        REAL                 :: R_max_ES
        REAL                 :: lum_q
        REAL, PARAMETER      :: Q_Halpha_5=6.71E-25     !at 5000K
        REAL, PARAMETER      :: Q_Halpha_10=3.56E-25    !at 10000K
        REAL, PARAMETER      :: Q_Halpha_20=1.83E-24    !at 20000K

        IF (ES_temp==5000) THEN
            lum_q=L_Halpha/Q_Halpha_5
        ELSE IF (ES_temp==10000) THEN
            lum_q=L_Halpha/Q_Halpha_10
        ELSE IF (ES_temp==20000) THEN
            lum_q=L_Halpha/Q_Halpha_20
        END IF


        IF (lg_ES) THEN
        IF ((3-2*gas_geometry%rho_power==0) .OR. (1-gas_geometry%rho_power==0)) THEN
            PRINT*,'You have selected a density profile with exponent 1.5 or 1.0 - need an alternative calculation in this case'
            STOP
        END IF
        END IF

        !At some point, include coding for separate line emitting region and gas region
        R_max_ES=gas_geometry%R_max
        
        ES_const=(lum_q**0.5)*(((3-2*gas_geometry%rho_power)/(4*pi))/(((R_max_ES*1E15)**(3-2*gas_geometry%rho_power)-(gas_geometry%R_min*1E15)**(3-2*gas_geometry%rho_power))))**0.5

        PRINT*,'av e- density',(ES_const*1E20*((1E15)**(-gas_geometry%rho_power))*((R_max_ES)**(3-gas_geometry%rho_power)-(gas_geometry%R_min)**(3-gas_geometry%rho_power)))/((3-dust_geometry%rho_power)*((R_max_ES)**3-(gas_geometry%R_min)**3))
        PRINT*,''
        PRINT*,'e- optical depth',ES_const*6.6E-5*((1E15*R_max_ES)**(1-gas_geometry%rho_power)-(1E15*gas_geometry%R_min)**(1-gas_geometry%rho_power))/(1-gas_geometry%rho_power)

    END SUBROUTINE



END MODULE
