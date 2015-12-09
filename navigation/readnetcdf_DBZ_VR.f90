
    PROGRAM readnetcdf 
!
! Program to read netcdf aircraft sweepfiles.
! Huaqing Cai April 2010
! Niles Oien  March 2010
!   This program reads in each netcdf file and output a text file
!   Text files are number 1 through nfile

    IMPLICIT none 
    INCLUDE 'netcdf.inc'

!  Variable Declaration 
    INTEGER ntimes    ! Number of times rays were collected --- Number of rays
    INTEGER nranges   ! Number of range gates per time  --- Number of gates for each ray
    INTEGER nsweeps   ! Number of sweep in each netcdf file; usually it is one

    INTEGER*4 counter,T_RAYS  ! Ray number and total rays after each sweep
    INTEGER RETVAL          ! Return value for netcdf calls.
    INTEGER nid             ! NETCDF ID
    INTEGER i,j,k,ii,jj,kk  ! General purpose
    INTEGER NUM_MISSING     ! Number of gates that had missing data

    LOGICAL FIRST     ! Used to get min, max

!  Scaler variable for coccrection factors
 
    REAL azimuth_correction   
    REAL elevation_correction
    REAL range_correction
    REAL longitude_correction
    REAL latitude_correction
    REAL pressure_altitude_correction
! Change radar_altitude_correction to altitude_correction
!   REAL radar_altitude_correction
    REAL altitude_correction
! Change ew_gound_speed_correction to eastward_velocity_correction
!   REAL ew_gound_speed_correction
    REAL eastward_velocity_correction
! Change ns_ground_speed_correction to northward_velocity_correction
!   REAL ns_ground_speed_correction
    REAL northward_velocity_correction
    REAL vertical_velocity_correction
    REAL heading_correction
    REAL roll_correction
    REAL pitch_correction
    REAL drift_correction
    REAL rotation_correction
    REAL tilt_correction

! Miscellaneous variables

    REAL MAX,MIN      ! max, min values

    REAL BAD_DATA_VALUE
    PARAMETER(BAD_DATA_VALUE=-999.0) ! Value we use for bad data

! One dimensional dynamical array of azimuths, one per time

    INTEGER,ALLOCATABLE,DIMENSION(:) :: sweep_number
    
    REAL*8,ALLOCATABLE,DIMENSION(:) :: time   
    REAL,ALLOCATABLE,DIMENSION(:) :: range 
    REAL,ALLOCATABLE,DIMENSION(:) :: azimuth      
    REAL,ALLOCATABLE,DIMENSION(:) :: elevation   
    REAL*8,ALLOCATABLE,DIMENSION(:) :: latitude
    REAL*8,ALLOCATABLE,DIMENSION(:) :: longitude
    REAL*8,ALLOCATABLE,DIMENSION(:) :: altitude
    REAL,ALLOCATABLE,DIMENSION(:) :: altitude_agl
    REAL,ALLOCATABLE,DIMENSION(:) :: heading 
    REAL,ALLOCATABLE,DIMENSION(:) :: roll
    REAL,ALLOCATABLE,DIMENSION(:) :: pitch 
    REAL,ALLOCATABLE,DIMENSION(:) :: drift
    REAL,ALLOCATABLE,DIMENSION(:) :: rotation 
    REAL,ALLOCATABLE,DIMENSION(:) :: tilt 
! Change ew_velocity to eastward_velocity
! Change ns_velocity to northward_velocity
!   REAL,ALLOCATABLE,DIMENSION(:) :: ew_velocity 
!   REAL,ALLOCATABLE,DIMENSION(:) :: ns_velocity
    REAL,ALLOCATABLE,DIMENSION(:) :: eastward_velocity 
    REAL,ALLOCATABLE,DIMENSION(:) :: northward_velocity 

    REAL,ALLOCATABLE,DIMENSION(:) :: vertical_velocity 
! Change ew_wind to eastward_wind, ns_wind to northward_wind
!   REAL,ALLOCATABLE,DIMENSION(:) :: ew_wind
!   REAL,ALLOCATABLE,DIMENSION(:) :: ns_wind 
    REAL,ALLOCATABLE,DIMENSION(:) :: eastward_wind 
    REAL,ALLOCATABLE,DIMENSION(:) :: northward_wind 

    REAL,ALLOCATABLE,DIMENSION(:) :: vertical_wind 

! Two dimensional dynamical array of DBZ, VR, SW, NCP, etc 
   
    REAL,ALLOCATABLE,DIMENSION(:,:) :: DBZ  
    REAL,ALLOCATABLE,DIMENSION(:,:) :: NCP 
    REAL,ALLOCATABLE,DIMENSION(:,:) :: VR 
    REAL,ALLOCATABLE,DIMENSION(:,:) :: SW 
 
! Two dimensional dynamical array of 2 byte ints as stored in file
    INTEGER*2,ALLOCATABLE,DIMENSION(:,:) :: I2DATA    
    REAL SCALE  ! Scale used to convert ints to reals
    REAL OFFSET ! Offset used to convert ints to reals
    INTEGER*2 MISSING ! Value used to indicate missing data in the netcdf file

! Variables for start and end time for one netcdf file

! Change start_time to time_coverage_start
! Change end_time to time_coverage_end

!   CHARACTER(len=64) start_time,end_time
    CHARACTER(len=64) time_coverage_start,time_coverage_end 

! Variables for input file list
    CHARACTER(len=100) filename,outfilename , argu
    INTEGER  nfile     ! number of netcdf files

!   Initialization
    counter = 0
    T_RAYS = 0

!====================================================================
!   Read the text file contains all the file names need to be read  
!   The file name is: input_file_list.txt
!====================================================================
! get the input file name from command line
       CALL GETARG(1, argu)

       open(10, file=argu, status='old')
       read(10,*) nfile
       do ii=1,nfile
         read(10,*)filename 
         print *, 'Total Files =',nfile,'file=',ii,'  ',filename
     
! Open the output file for storing ray information
       write(outfilename,'(i10)') ii
       outfilename =  TRIM(ADJUSTL(outfilename)) // '.txt'

       open(20, file=outfilename, status='new')
!
! Open the netcdf file.
!
	RETVAL = NF_OPEN(TRIM(filename), NF_NOWRITE, nid)
	CALL CHECK_NETCDF(RETVAL, 'Failed to open input file')
!
! Read the number of sweeps, times and ranges. These are
! dimensions rather than variables and so are stored differently.
!
        RETVAL = NF_INQ_DIMID(nid, 'sweep', I)
        CALL CHECK_NETCDF(RETVAL, 'Failed to find a dimension named sweep')
        RETVAL = NF_INQ_DIMLEN(NID, I, NSWEEPS)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read number of sweeps')

	RETVAL = NF_INQ_DIMID(nid, 'time', I)
	CALL CHECK_NETCDF(RETVAL, 'Failed to find a dimension named time')
	RETVAL = NF_INQ_DIMLEN(NID, I, NTIMES)
	CALL CHECK_NETCDF(RETVAL, 'Failed to read number of times')

	RETVAL = NF_INQ_DIMID(NID, 'range', I)
	CALL CHECK_NETCDF(RETVAL, 'Failed to find a dimension named range')
	RETVAL = NF_INQ_DIMLEN(NID, I, NRANGES)
	CALL CHECK_NETCDF(RETVAL, 'Failed to read number of ranges')

	PRINT *, NSWEEPS,' sweeps ', NTIMES, ' beams found, each with ', NRANGES, ' gates.'
!
! Read scalar variables - all the correction factors.
!
        RETVAL = NF_INQ_VARID(NID, 'azimuth_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named elevation_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, azimuth_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read azimuth correction.')

        RETVAL = NF_INQ_VARID(NID, 'elevation_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named elevation_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, elevation_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read elevation correction.')

        RETVAL = NF_INQ_VARID(NID, 'range_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named range_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, range_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read range correction.')

! Change range_correction unit from meters to km:
        range_correction = range_correction/1000.0

        RETVAL = NF_INQ_VARID(NID, 'longitude_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named longitude_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, longitude_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read longitude correction.')

        RETVAL = NF_INQ_VARID(NID, 'latitude_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named latitude_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, latitude_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read  latitude correction.')

        RETVAL = NF_INQ_VARID(NID, 'pressure_altitude_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named pressure_altitude_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, pressure_altitude_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read pressure altitude correction.')

! change  pressure_altitude_correction unit from meters to km:
        pressure_altitude_correction = pressure_altitude_correction/1000.0

        RETVAL = NF_INQ_VARID(NID, 'altitude_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named altitude_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, altitude_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read altitude correction.')

! change  altitude_correction unit from meters to km:
        altitude_correction = altitude_correction/1000.0

        RETVAL = NF_INQ_VARID(NID, 'eastward_velocity_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named eastward_velocity_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, eastward_velocity_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read eastward_velocity_correction.')

        RETVAL = NF_INQ_VARID(NID, 'northward_velocity_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named northward_velocity_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, northward_velocity_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read northward_velocity_correction.')

        RETVAL = NF_INQ_VARID(NID, 'vertical_velocity_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named vertical_velocity_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, vertical_velocity_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read vertical_velocity correction.')

        RETVAL = NF_INQ_VARID(NID, 'heading_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named heading_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, heading_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read heading correction.')

        RETVAL = NF_INQ_VARID(NID, 'roll_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named roll_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, roll_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read roll correction.')

        RETVAL = NF_INQ_VARID(NID, 'pitch_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named pitch_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, pitch_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read pitch correction.')

        RETVAL = NF_INQ_VARID(NID, 'drift_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named drift_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, drift_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read drift correction.')

        RETVAL = NF_INQ_VARID(NID, 'rotation_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named rotation_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, rotation_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read rotation correction.')

        RETVAL = NF_INQ_VARID(NID, 'tilt_correction', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named tilt_correction found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, tilt_correction)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read tilt correction.')

!  Get the start_time  of the sweep
!
        RETVAL = NF_INQ_VARID(NID, 'time_coverage_start', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named time_coverage_start found.')
        RETVAL = NF_GET_VAR_TEXT(NID, I, time_coverage_start)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read time_coverage_start.')

        PRINT *,'start_time of the sweep: ', TRIM(time_coverage_start(1:20)) 

        PRINT *, 'Correction Factors : '
        PRINT *, 'azimuth_correction =',azimuth_correction
        PRINT *, 'elevation_correction =',elevation_correction
        PRINT *, 'range_correction =',range_correction
        PRINT *, 'longitude_correction =',longitude_correction
        PRINT *, 'latitude_correction =',latitude_correction
        PRINT *, 'pressure_altitude_correction =',pressure_altitude_correction
        PRINT *, 'altitude_correction =',altitude_correction
        PRINT *, 'eastward_velocity_correction =',eastward_velocity_correction
        PRINT *, 'northward_velocity_correction =',northward_velocity_correction
        PRINT *, 'vertical_velocity_correction =',vertical_velocity_correction
        PRINT *, 'heading_correction =',heading_correction
        PRINT *, 'roll_correction =',roll_correction
        PRINT *, 'pitch_correction =',pitch_correction
        PRINT *, 'drift_correction =',drift_correction
        PRINT *, 'rotation_correction =',rotation_correction
        PRINT *, 'tilt_correction =',tilt_correction


!
! Read one dimensional array - the time, range, azimuths, elevation, etc
!
!  ========== Sweep ============
        ALLOCATE(sweep_number(nsweeps)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'sweep_number', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named sweep_number found.')
        RETVAL = NF_GET_VAR_INT(NID, I, sweep_number)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read sweep_number.')

        DO I=1, NSWEEPS
!        PRINT *, I, ' sweep_number is ', sweep_number(I)
        END DO

!  ========== Time ============
        ALLOCATE(time(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'time', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named time found.')
        RETVAL = NF_GET_VAR_DOUBLE(NID, I, time)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read time.')

        DO I=1, NTIMES
!        PRINT *, I, ' time is ', time(I)
        END DO
 
!  ========== range ============
        ALLOCATE(range(nranges)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'range', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named range found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, range)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read range.')


! Change range from meters to km, since in the netcdf file, range is in 
! meters, not km anymore
        DO I=1, NRANGES
          range(I) = range(I)/1000.0
        END DO

        DO I=1, NRANGES
!        PRINT *, I, ' range is ', range(I)
        END DO

!  ========== azimuth ============
        ALLOCATE(azimuth(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'azimuth', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named azimuth found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, azimuth)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read azimuth.')

        DO I=1, NTIMES
!        PRINT *, I, ' azimuth is ', azimuth(I)
        END DO

!  ========== elevation ============
        ALLOCATE(elevation(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'elevation', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named elevation found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, elevation)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read elevation.')

        DO I=1, NTIMES
!        PRINT *, I, ' elevation is ', elevation(I)
        END DO

!  ========== latitude ============
        ALLOCATE(latitude(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'latitude', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named latitude found.')
        RETVAL = NF_GET_VAR_DOUBLE(NID, I, latitude)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read latitude.')

        DO I=1, NTIMES
!        PRINT *, I, ' latitude is ', latitude(I)
        END DO

!  ========== longitude ============
        ALLOCATE(longitude(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'longitude', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named longitude found.')
        RETVAL = NF_GET_VAR_DOUBLE(NID, I, longitude)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read longitude.')

        DO I=1, NTIMES
!        PRINT *, I, ' longitude is ', longitude(I)
        END DO

!  ========== altitude ============
        ALLOCATE(altitude(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'altitude', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named altitude found.')
        RETVAL = NF_GET_VAR_DOUBLE(NID, I, altitude)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read altitude.')

        DO I=1, NTIMES
!        PRINT *, I, ' altitude is ', altitude(I)
        END DO
! change unit of altitude from meter to km:
        DO I=1, NTIMES
          altitude(I) = altitude(I)/1000.0
        END DO




!  ========== altitude_agl ============
        ALLOCATE(altitude_agl(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'altitude_agl', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named altitude_agl found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, altitude_agl)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read altitude_agl.')

        DO I=1, NTIMES
!        PRINT *, I, ' altitude_agl is ', altitude_agl(I)
        END DO

! change unit of altitude_agl from meter to km:

        DO I=1, NTIMES
         altitude_agl(I) = altitude_agl(I)/1000.0
        END DO


!  ========== heading ============
        ALLOCATE(heading(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'heading', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named heading found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, heading)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read heading.')

        DO I=1, NTIMES
!        PRINT *, I, ' heading is ', heading(I)
        END DO

!  ========== roll ============
        ALLOCATE(roll(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'roll', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named roll found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, roll)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read roll.')

        DO I=1, NTIMES
!        PRINT *, I, ' roll is ', roll(I)
        END DO

!  ========== pitch ============
        ALLOCATE(pitch(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'pitch', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named pitch found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, pitch)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read pitch.')

        DO I=1, NTIMES
!        PRINT *, I, ' pitch is ', pitch(I)
        END DO

!  ========== drift ============
        ALLOCATE(drift(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'drift', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named drift found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, drift)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read drift.')

        DO I=1, NTIMES
!        PRINT *, I, ' drift is ', drift(I)
        END DO

!  ========== rotation ============
        ALLOCATE(rotation(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'rotation', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named rotation found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, rotation)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read rotation.')

        DO I=1, NTIMES
!        PRINT *, I, ' rotation is ', rotation(I)
        END DO

!  ========== tilt ============
        ALLOCATE(tilt(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'tilt', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named tilt found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, tilt)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read tilt.')

        DO I=1, NTIMES
!        PRINT *, I, ' tilt is ', tilt(I)
        END DO

!  ========== eastward_velocity ============
        ALLOCATE(eastward_velocity(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'eastward_velocity', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named eastward_velocity found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, eastward_velocity)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read eastward_velocity.')

        DO I=1, NTIMES
!        PRINT *, I, ' eastward_velocity is ', eastward_velocity(I)
        END DO

!  ========== northward_velocity ============
        ALLOCATE(northward_velocity(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'northward_velocity', I)
        CALL CHECK_NETCDF(RETVAL,'No variable named northward_velocity found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, northward_velocity)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read northward_velocity.')

        DO I=1, NTIMES
!        PRINT *, I, ' northward_velocity is ', northward_velocity(I)
        END DO

!  ========== vertical_velocity ============
        ALLOCATE(vertical_velocity(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'vertical_velocity', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named vertical_velocity found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, vertical_velocity)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read vertical_velocity.')

        DO I=1, NTIMES
!        PRINT *, I, ' vertical_velocity is ', vertical_velocity(I)
        END DO

!  ========== eastward_wind ============
        ALLOCATE(eastward_wind(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'eastward_wind', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named eastward_wind found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, eastward_wind)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read eastward_wind.')

        DO I=1, NTIMES
!        PRINT *, I, ' eastward_wind is ', eastward_wind(I)
        END DO

!  ========== northward_wind ============
        ALLOCATE(northward_wind(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'northward_wind', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named northward_wind found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, northward_wind)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read northward_wind.')

        DO I=1, NTIMES
!        PRINT *, I, ' northward_wind is ', northward_wind(I)
        END DO

!  ========== vertical_wind ============
        ALLOCATE(vertical_wind(ntimes)) ! Dynamic memory allocation
        RETVAL = NF_INQ_VARID(NID, 'vertical_wind', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named vertical_wind found.')
        RETVAL = NF_GET_VAR_REAL(NID, I, vertical_wind)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read vertical_wind.')

        DO I=1, NTIMES
!        PRINT *, I, ' vertical_wind is ', vertical_wind(I)
         IF(abs(vertical_wind(I)).gt.100.0) THEN
            vertical_wind(I) = -999.0
         ENDIF
        END DO

!
! Read two-dimensional arrays such as DBZ. Read DBZ integers and scale them into real dbz values.
!
! ========== DBZ ======================
        ALLOCATE(DBZ(NRANGES,NTIMES))
        ALLOCATE(I2DATA(NRANGES,NTIMES))

        RETVAL = NF_INQ_VARID(NID, 'DBZ', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named DBZ found.')
        RETVAL = NF_GET_VAR_INT2(NID, I, I2DATA)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read dbz.')
!
! Also need to get scale, offset which are attributes for this variable.
!

! Change attibutes: from missing_value to: _FillValue

        RETVAL = NF_GET_ATT_REAL(NID,I,'scale_factor', SCALE)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get scale attribute for DBZ variable.')
        RETVAL = NF_GET_ATT_REAL(NID,I,'add_offset', OFFSET)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get offset attribute for DBZ variable.')
!       RETVAL = NF_GET_ATT_INT2(NID,I,'missing_value', MISSING)
        RETVAL = NF_GET_ATT_INT2(NID,I,'_FillValue', MISSING)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get missing attribute for DBZ variable.')
!
! Apply scale, offset to get dbz values. Count how many were missing as well.
!
        NUM_MISSING = 0
        FIRST = .TRUE.
        DO J=1,NTIMES
         DO K=1,NRANGES
          IF (I2DATA(K,J) .EQ. MISSING) THEN
           DBZ(K,J) = BAD_DATA_VALUE
           NUM_MISSING = NUM_MISSING + 1
          ELSE
           DBZ(K,J) = SCALE*I2DATA(K,J) + OFFSET
           IF (FIRST) THEN
                MIN = DBZ(K,J)
                MAX = MIN
                FIRST = .FALSE.
           ELSE 
                IF (DBZ(K,J) .LT. MIN) MIN = DBZ(K,J)
                IF (DBZ(K,J) .GT. MAX) MAX = DBZ(K,J)
           END IF
          END IF
         END DO
        END DO

        DEALLOCATE(I2DATA) ! We are done with the interger array

        PRINT *, NUM_MISSING, ' of ', NTIMES*NRANGES, ' gates had missing dbz values.'
        PRINT *, 'DBZ Min, Max ', MIN, ' to ', MAX
        

! ========== NCP ======================
        ALLOCATE(NCP(NRANGES,NTIMES))
        ALLOCATE(I2DATA(NRANGES,NTIMES))

        RETVAL = NF_INQ_VARID(NID, 'NCP', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named NCP found.')
        RETVAL = NF_GET_VAR_INT2(NID, I, I2DATA)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read ncp.')
!
! Also need to get scale, offset which are attributes for this variable.
!
        RETVAL = NF_GET_ATT_REAL(NID,I,'scale_factor', SCALE)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get scale attribute for NCP variable.')
        RETVAL = NF_GET_ATT_REAL(NID,I,'add_offset', OFFSET)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get offset attribute for NCP variable.')
        RETVAL = NF_GET_ATT_INT2(NID,I,'_FillValue', MISSING)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get missing attribute for NCP variable.')
!
! Apply scale, offset to get ncp values. Count how many were missing as well.
!
        NUM_MISSING = 0
        FIRST = .TRUE.
        DO J=1,NTIMES
         DO K=1,NRANGES
          IF (I2DATA(K,J) .EQ. MISSING) THEN
           NCP(K,J) = BAD_DATA_VALUE
           NUM_MISSING = NUM_MISSING + 1
          ELSE
           NCP(K,J) = SCALE*I2DATA(K,J) + OFFSET
           IF (FIRST) THEN
                MIN = NCP(K,J)
                MAX = MIN
                FIRST = .FALSE.
           ELSE
                IF (NCP(K,J) .LT. MIN) MIN = NCP(K,J)
                IF (NCP(K,J) .GT. MAX) MAX = NCP(K,J)
           END IF
          END IF
         END DO
        END DO

        DEALLOCATE(I2DATA) ! We are done with the interger array

        PRINT *, NUM_MISSING, ' of ', NTIMES*NRANGES, ' gates had missing ncp values.'
        PRINT *, 'NCP Min, Max ', MIN, ' to ', MAX


! ========== VR ======================
        ALLOCATE(VR(NRANGES,NTIMES))
        ALLOCATE(I2DATA(NRANGES,NTIMES))

        RETVAL = NF_INQ_VARID(NID, 'VR', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named VR found.')
        RETVAL = NF_GET_VAR_INT2(NID, I, I2DATA)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read VR.')
!
! Also need to get scale, offset which are attributes for this variable.
!
        RETVAL = NF_GET_ATT_REAL(NID,I,'scale_factor', SCALE)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get scale attribute for VR variable.')
        RETVAL = NF_GET_ATT_REAL(NID,I,'add_offset', OFFSET)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get offset attribute for VR variable.')
        RETVAL = NF_GET_ATT_INT2(NID,I,'_FillValue', MISSING)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get missing attribute for VR variable.')
!
! Apply scale, offset to get dbz values. Count how many were missing as well.
!
        NUM_MISSING = 0
        FIRST = .TRUE.
        DO J=1,NTIMES
         DO K=1,NRANGES
          IF (I2DATA(K,J) .EQ. MISSING) THEN
           VR(K,J) = BAD_DATA_VALUE
           NUM_MISSING = NUM_MISSING + 1
          ELSE
           VR(K,J) = SCALE*I2DATA(K,J) + OFFSET
           IF (FIRST) THEN
                MIN = VR(K,J)
                MAX = MIN
                FIRST = .FALSE.
           ELSE
                IF (VR(K,J) .LT. MIN) MIN = VR(K,J)
                IF (VR(K,J) .GT. MAX) MAX = VR(K,J)
           END IF
          END IF
         END DO
        END DO

        DEALLOCATE(I2DATA) ! We are done with the interger array

        PRINT *, NUM_MISSING, ' of ', NTIMES*NRANGES, ' gates had missing vr values.'
        PRINT *, 'VR Min, Max ', MIN, ' to ', MAX


! ========== SW ======================
        ALLOCATE(SW(NRANGES,NTIMES))
        ALLOCATE(I2DATA(NRANGES,NTIMES))

        RETVAL = NF_INQ_VARID(NID, 'SW', I)
        CALL CHECK_NETCDF(RETVAL, 'No variable named SW found.')
        RETVAL = NF_GET_VAR_INT2(NID, I, I2DATA)
        CALL CHECK_NETCDF(RETVAL, 'Failed to read sw.')
!
! Also need to get scale, offset which are attributes for this variable.
!
        RETVAL = NF_GET_ATT_REAL(NID,I,'scale_factor', SCALE)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get scale attribute for SW variable.')
        RETVAL = NF_GET_ATT_REAL(NID,I,'add_offset', OFFSET)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get offset attribute for SW variable.')
        RETVAL = NF_GET_ATT_INT2(NID,I,'_FillValue', MISSING)
        CALL CHECK_NETCDF(RETVAL, 'Failed to get missing attribute for SW variable.')
!
! Apply scale, offset to get dbz values. Count how many were missing as well.
!
        NUM_MISSING = 0
        FIRST = .TRUE.
        DO J=1,NTIMES
         DO K=1,NRANGES
          IF (I2DATA(K,J) .EQ. MISSING) THEN
           SW(K,J) = BAD_DATA_VALUE
           NUM_MISSING = NUM_MISSING + 1
          ELSE
           SW(K,J) = SCALE*I2DATA(K,J) + OFFSET
           IF (FIRST) THEN
                MIN = SW(K,J)
                MAX = MIN
                FIRST = .FALSE.
           ELSE
                IF (SW(K,J) .LT. MIN) MIN = SW(K,J)
                IF (SW(K,J) .GT. MAX) MAX = SW(K,J)
           END IF
          END IF
         END DO
        END DO

        DEALLOCATE(I2DATA) ! We are done with the interger array

        PRINT *, NUM_MISSING, ' of ', NTIMES*NRANGES, ' gates had missing sw values.'
        PRINT *, 'SW Min, Max ', MIN, ' to ', MAX


!
! Print three beams for testing
!
        DO I=1,3
!        PRINT *, 'Beam ', I, ' at azimuth ', azimuth(I)
         DO J=1,NRANGES
!           PRINT *, 'Gate ', J, ' is ', DBZ(J,I)
         END DO
        END DO

!
! Write the ray informatin to the output file
!        
        DO I=1,NTIMES
           counter = I+T_RAYS
           write(20,101)counter,filename,sweep_number(nsweeps) &
              ,NTIMES,NRANGES &
              ,time_coverage_start(1:4),time_coverage_start(6:7) &
              ,time_coverage_start(9:10) &
              ,time_coverage_start(12:13),time_coverage_start(15:16) &
              ,time_coverage_start(18:19),time(I) & 
              ,azimuth(I),elevation(I),latitude(I),longitude(I),altitude(I) & 
              ,altitude_agl(I),heading(I),roll(I),pitch(I),drift(I) & 
              ,rotation(I),tilt(I),eastward_velocity(I) & 
              ,northward_velocity(I),vertical_velocity(I),eastward_wind(I) &
              ,northward_wind(I) & 
              ,vertical_wind(I),azimuth_correction,elevation_correction & 
              ,range_correction,longitude_correction,latitude_correction & 
              ,pressure_altitude_correction,altitude_correction & 
              ,eastward_velocity_correction,northward_velocity_correction & 
              ,vertical_velocity_correction,heading_correction & 
              ,roll_correction,pitch_correction,drift_correction & 
              ,rotation_correction,tilt_correction 
           write(20,102)counter,(range(J), J=1,nranges)     
           write(20,102)counter,(DBZ(J,I),J=1,nranges)
           write(20,102)counter,(NCP(J,I),J=1,nranges)
           write(20,102)counter,(VR(J,I),J=1,nranges)
           write(20,102)counter,(SW(J,I),J=1,nranges)
        ENDDO
 101   format(I10,2x,a50,3I10,a5,5a3,d20.8,2f10.4,3d20.8,29f10.4)
 102   format(I10,800f10.4) 

! Increase the ray counter by NTIMES; this is done after each sweep
        T_RAYS = T_RAYS + NTIMES          

! Close the file just read
	RETVAL = NF_CLOSE(NID)
	IF (RETVAL .NE. 0) THEN
	 PRINT *, 'INFO : Unable to close ', TRIM(FILENAME) ! Really just info, not abnormal termination
	END IF


! ====== DEALLOCATE all dynamical arrays  ============
!  One dimensional arrays
        DEALLOCATE(sweep_number)

        DEALLOCATE(time)
        DEALLOCATE(range)
        DEALLOCATE(azimuth)
        DEALLOCATE(elevation)
        DEALLOCATE(latitude)
        DEALLOCATE(longitude)
        DEALLOCATE(altitude)
        DEALLOCATE(altitude_agl)
        DEALLOCATE(heading)
        DEALLOCATE(roll)
        DEALLOCATE(pitch)
        DEALLOCATE(drift)
        DEALLOCATE(rotation)
        DEALLOCATE(tilt)
        DEALLOCATE(eastward_velocity)
        DEALLOCATE(northward_velocity)
        DEALLOCATE(vertical_velocity)
        DEALLOCATE(eastward_wind)
        DEALLOCATE(northward_wind)
        DEALLOCATE(vertical_wind)

! two dimensional arrays
        DEALLOCATE(DBZ)
        DEALLOCATE(NCP)
        DEALLOCATE(VR)
        DEALLOCATE(SW)
! close the output file
        close(20)

    enddo       !   end of netcdf file loop, go to next file


    STOP
    END

!---------------------------------------------------
!
! Routine to check if a netcdf call has gone well.
!

	SUBROUTINE CHECK_NETCDF(RETVAL, ERRSTR)

	INTEGER RETVAL
	CHARACTER(*) ERRSTR

	IF (RETVAL .NE. 0) THEN
	 PRINT *, TRIM(ERRSTR)
	 PRINT *, 'Netcdf return value was ', RETVAL
	 PRINT *, 'Abnormal termination.'
	STOP
	ENDIF

	RETURN
	END

