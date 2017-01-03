MODULE class_packet

USE globals

    IMPLICIT NONE

TYPE packet_obj
    REAL :: r  !current radius of packet location in cm
    REAL :: v  !current velocity of packet in km/s
    INTEGER :: iG   !current cell ID of packet

END TYPE

TYPE(packet_obj) :: packet

END MODULE class_packet
