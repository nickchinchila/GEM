      ! This program converts a .txt file into a .dat file for GEM 
      ! observations at 2300 MHz in Brazil. The read statement in line
	! 172 is for files in may 1999 when ch14 substituted ch11 in the
	! DAS box for acquiring the azimuth encoder data. See also line 207
	! Temperature statistics to control the gain susceptability are 
	! based on a 2*TA_min detection level of 0.016K (Tsys = 57K, tau =
	! 0.56002s and BW = 100MHz), Gain = 47K/V, etagain = -0.022058/deg C, 
	! S_cal = 1.32276V at T2_cal = 307.5991 in 99_05_2207 file, such that
	! 
	!			Delta_T2_max = [etagain*(S_cal*Gain/0.016 - 1)]^(-1)
	!						 = 0.0117 deg C
	!
	! 03/07/2003:	the time check between filenames and 1st frame 
	!				timestamps is logged for convenient inspection.
	! 30/08/2003:	Bishop data processing requires adaptation. Arrays
	!				resized from 16500 to 60000. Use V_max = 9.982.
	!				Change read statement, but don't use Odetics time.
	!				Write statement uses one more digit for elevation.
	! 31/08/2003:	create filegroup 6 to speed up Bishop data iteration 
	!				for Vmax.
	! 27/12/2004:	new approach to file reading by keeping track of 
	!				frame number differences.
	! 31/08/2005:	adapted for 2.3 GHz data, including Cesar's time bug.
	!				Offset to timeini inplemented to synchronize with filename.
	!				T1 measures the temperature at the HEMT, but T2 at the diode
	!				is more stable and is used for corrections later.
	! 07/10/2005:	found bug in timestamp for filename when timestamp
	!				increases beyond 24 hr: the day is not increased. Fixed.
      !****************************************************************
      !$DEBUG 
      INTEGER*4	len_name(50),jfini(50),jffin(50),jcount(50)

      REAL*8      encoff,azoff,time, elevat, signal, azimut, T1, T2, T3,
     &            T4, HI, Tns, Vns, timeini, ml, timefin,mdtaux,mdtime,
     &			Vmaxaz,timefra, timeaux, timediff, tioff,dtime(50),
     &			temp(60000,4),tempavg(4,2),delta(4),tem(60000),sdtime,
     &			sT1(50),dT1(50),sT2(50),dT2(50),sT3(50),dT3(50),
     &			sT4(50),dT4(50)
     
      CHARACTER   filenum*8,filelist*9,tag*17,shh*2,smm*2,sms*4,
     &			fullname(50)*100,newlist*9,syear*3,smonth*2,sday*2,
     &			s2year*2
	
      PRINT *,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     & %%%%%%%%%%%%%%%%%%%%'
      PRINT *,'               GEM 0 : CONVERSION OF DIGITAL TO PHYSICAL 
     &UNITS'
      PRINT *,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     &%%%%%%%%%%%%%%%%%%%%'

	PRINT *,' Single file    : [ -> 1]'
	PRINT *,' Multiple files : [ -> 2]'
	READ(*,*) jlist
	PRINT *,' '
	PRINT *,' SELECT the observational site:'
	PRINT *,'	1	->	Bishop'
	PRINT *,'	2	->	Cachoeira Paulista'
	READ(*,*) jsite
	PRINT *,' '
	
	IF (jsite.GT.1) THEN

		PRINT *,' WERE the observations taken in 1998 or after 03/04/05?'
		PRINT *,'	1	:	YES'
		PRINT *,'	2	:	NO'
		READ(*,*) j199
		PRINT *,' '

		IF (j199.EQ.1) THEN
			jobs	=	2
		ELSE
			jobs	=	3
		END IF

	ELSE

		jobs	=	1

	END IF
	

      !::::::::::::::::::::: encoders section ::::::::::::::::::::::::::

	IF (jsite.EQ.1) THEN

		V_180   = 4.989     ! effective half scale for elevation			! Bishop
		Vmaxel  = 9.980     ! maximum voltage reading of elevation encoder

	ELSE

		V_180   = 4.974     ! effective half scale for elevation			! Cachoeira
	    Vmaxel  = 9.995     ! maximum voltage reading of elevation encoder

	END IF

	PRINT *,' SET Vmax?'! for the azimuth encoder
	PRINT *,'    1 -> YES'
     	PRINT *,'    2 ->  NO'
	READ(*,*) jset
	PRINT *,' '

	IF (jset.EQ.1) THEN

		PRINT *,' Vmax:'
		READ(*,*) Vmaxaz
		PRINT *,' '

	ELSE IF (jsite.EQ.1) THEN

		Vmaxaz	= 9.9211 ! +/- 0.0010 V as given by Vmax1 from in 
						! Vmax.for using all the original files in 
						! fileli500.

	ELSE

C		Vmaxaz  = 9.89485   ! Weighted mean using iterative process with 
                          ! azimuth.for in 15 files of filelit2.txt. The
                          ! precise value gives 9.8947 +/- 0.0015 Volts, 
                          ! whereas the encoder step gives 0.0024 Volts.

		Vmaxaz  = 9.8993 ! Weighted mean using iterative process with
                          ! azimuth.for in 102 files in the 1st group of
						! filelists (12,18,13,19,20,21). The 1-sigma
                          ! variance is +/- 0.0004 Volts.

C		Vmaxaz	= 9.7373	! Weighted mean using in 2004 iterative process as
							! above with +/- 0.0047 Volts with fileli601 

C		Vmaxaz	= 9.7887	! Weighted mean in 2005 using iterative process as
							! above with +/- 0.0052 Volts with fileli638 

	END IF

      ! offset in degrees. Must be zero to use GEM 1::::::::::::::::::::

      azoff   =   0.D0                ! Uncalibrated encoder    
      encoff  = 360.0*(1.0-1.0/4096)  !angle of encoder step at V_max
	
      !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	l	=	0
	IF (jlist.EQ.1) THEN

		l	=	l	+	1

		PRINT *,' ENTER the .txt file name:'
		READ(*,'(A)') fullname(l)

		len_name(l)	= LEN_TRIM(fullname(l))

		OPEN(3,	file	='c:\GEM\2300mhz\Temps\Tem',
     &                    //fullname(l)(1:(len_name(l)-4))//'.dat',			
     &			status  = 'unknown',
     &			err     = 200)

		WRITE(3,33)

	ELSE

		PRINT *,' ENTER the file list name:'
		READ(*,*) filelist

	END IF
	PRINT *,' '

	PRINT *,' Evaluate Temperature statistics alone?'
	PRINT *,'     1  ->  YES'
	PRINT *,'     2  ->   NO'
	READ(*,*) jval
	PRINT *,' '

	IF (jval.EQ.1) GOTO 185

	IF (jlist.EQ.2) THEN

		OPEN(1, file	= 'c:\GEM\2300mhz\txt\'//filelist//'.txt', 
     &			status  = 'old',
     &			err     = 200)

		DO WHILE (.NOT.EOF(1))

			l	=	l	+	1
			
			READ(1,'(A)') fullname(l)

		END DO

		CLOSE(1)

		OPEN(6,file = 'c:\GEM\2300mhz\dat\time_check_'//filelist//
     &	'.dat',
     &		status  = 'unknown',
     &		err     = 200)

		WRITE(6,5)
    5 FORMAT('        filename       filetime  filespan frametime fram',
     &'espan N_frames    dtime')

		OPEN(9,file = 'c:\GEM\2300mhz\dat\time_offset_'//filelist//
     &	'.dat',
     &		status  = 'unknown',
     &		err     = 200)

		WRITE(9,6)
    6 FORMAT('   offset   mdtime  sdtime  iteration')

		OPEN(7,file	= 'c:\GEM\2300mhz\dat\log_'//filelist//'.dat',
     &		status  = 'unknown',
     &		err     = 200)

		WRITE(7,71)
   71 FORMAT('file    init    last    diff   count')

		OPEN(3,	file	='c:\GEM\2300mhz\Temps\Tem'//filelist//'.dat',			
     &			status  = 'unknown',
     &			err     = 200);

		WRITE(3,33)
   33 FORMAT('    T1      sT      dT1     T2      sT2     dT2     T3    
     &  sT3     dT3     T4      sT4     dT4     rate')

		DO ll = 1,l

			len_name(ll)	= LEN_TRIM(fullname(ll))

			OPEN(1,	file	= 'c:\GEM\2300mhz\txt\'
     &						//fullname(ll)(1:len_name(ll)), 
     &				status  = 'old',
     &				err     = 200)

			jk	=	0

			IF (jobs.EQ.1) THEN

				DO i = 1,3

					READ(1,*)

				END DO

			END IF

			DO WHILE (.NOT.EOF(1))

				READ(1,*,err=111) jframe
			
				jk	=	jk	+	1

				IF (jk.EQ.1) jfini(ll) = jframe

			END DO
  111			jffin(ll)	= jframe
			jcount(ll)	= jk

			CLOSE(1)

			WRITE(7,77) ll,jfini(ll),jffin(ll),
     &					(jffin(ll)-jfini(ll)+1),jk
   77 FORMAT(1x,I3,4(2x,I6))
			PRINT *,' file: ',ll,'   recorded: ',jk,'   sampled: ',
     &(jffin(ll)-jfini(ll)+1)

		END DO

		PRINT *,' '
		PRINT *,' ENTER the filelist name of the converted .dat files:'
		READ(*,'(A)') newlist

	END IF

	! Start the conversion :::::::::::::::::::::::::::::::::::::::::::::::::

	print *,fullname(l)(1:len_name(l))

	mit	=	0
	mdtime	=	0.D0
	mdtaux	=	0.D0
	tioff	=	0.D0

	DO WHILE ((DABS(mdtime).LT.DABS(mdtaux)).OR.(mit.LE.1))						! begin of offset loop

		mit	=	mit + 1

		IF (jlist.EQ.2) THEN

			OPEN(8, file	= 'c:\GEM\2300mhz\txt\'//newlist//'.txt',
     &				status  = 'unknown',
     &				err     = 200)

		END IF

	DO ll = 1,l																	! begin of loop over files in list

		OPEN(1,	file	= 'c:\GEM\2300mhz\txt\',
     &						fullname(ll)(1:len_name(ll)), 
     &			status  = 'old',
     &			err     = 200);			

	! Sample file name: 104_10_2120034594.txt ::::::::::::::::::::::::::::::

		PRINT *,' '

		IF (((jffin(ll)-jfini(ll)+1).EQ.jcount(ll)).OR.(jlist.EQ.1)) 
     &	THEN

			ju	=	INDEX(fullname(ll),'_')

			READ(fullname(ll)(1:ju-1),*) jyear
			jyear	= 1900 + jyear

			READ(fullname(ll)(ju+1:ju+2),*) month

			READ(fullname(ll)(ju+4:ju+5),*) jday

			READ(fullname(ll)(ju+6:len_name(ll)-4),*) jfilenum

			jhh = jfilenum/10**6 
			jmm = (jfilenum - jhh*10**6)/10**4
			jss = (jfilenum - jhh*10**6 - jmm*10**4)/10**2
			ml  = (jfilenum - jhh*10**6 - jmm*10**4 -jss*10**2)/100.D0 

			IF (ml.GE.0.50D0) THEN
			
				jss = jss - 1												! Cesar's bug
				
				ms	=	(jss + ml)*100							 

				WRITE(sms,'(I4)') ms

				IF (ms.EQ.0) THEN
					sms = '0000'
				ELSE IF (ms/10.EQ.0) THEN
					sms = '000'//sms(4:4)
				ELSE IF (ms/100.EQ.0) THEN
					sms	= '00'//sms(3:4)
				ELSE IF (ms/1000.EQ.0) THEN
					sms	= '0'//sms(2:4)
				END IF

				fullname(ll)	= 
     &			fullname(ll)(1:len_name(ll)-8)//sms//'.txt'

			END IF
      
			timeini = jhh*3600.D0 + jmm*60.D0 + jss*1.D0 + ml
			
			IF (jlist.EQ.2) THEN
			
				IF (jyear.LT.2000) THEN

					WRITE(8,221) fullname(ll)(1:ju+5),
     &				fullname(ll)(ju+6:len_name(ll)-4),jday,month,jyear

  221 FORMAT(A8,'	',A8,2('	',I2),'	',I4)									! filelist of good file names in usual format

				ELSE
				
					WRITE(8,222) fullname(ll)(1:ju+5),
     &				fullname(ll)(ju+6:len_name(ll)-4),jday,month,jyear

  222 FORMAT(A9,'	',A8,2('	',I2),'	',I4)									! filelist of good file names in usual format

					
				END IF

				OPEN(2,	file	= 'c:\GEM\2300mhz\dat\',
     &						fullname(ll)(1:len_name(ll)-4)//'.dat', 
     &					status  = 'unknown',
     &					err     = 200);

				tag	= fullname(ll)(1:len_name(ll)-4)
				
			END IF
			
			j99a	=	0													! skip mo frames							                      
			j99b	=	1

		ELSE


			DO i = 1,99

				READ(1,*) jframe											! exclude syndrome99 ...

			END DO

			j99a	=	99	 												! 99 after last frame in previous file

			j99b	= jffin(ll) - (jcount(ll)-100) - jffin(ll-1)			! 99 + lost frames after last frame in previous file.

			timeini	=	timefin + j99b*0.56002D0							! timestamp of the 100th frame in the syndrome99 affected file
      
			jhh = INT(timeini/3600)
			jmm = INT((timeini - jhh*3600)/60)
			ms	= (timeini - jhh*3600 - jmm*60)*100

			IF ((timeini - timeaux).LT.0.D0) jhh = jhh + 24					! syndrome99 names haven't got a datetag yet ... need update

			IF (jhh.GE.24) THEN

				jhh	= jhh - 24

				jday	= jday + 1

				IF (jday.GT.31) THEN
					
					month = month + 1

					IF (month.GT.12) THEN

						month	= month - 12

						jyear	= jyear + 1

					END IF

				END IF

			END IF

			IF (jyear.LT.2000) THEN

				WRITE(s2year,'(I2)') (jyear-1900)
			ELSE
				WRITE(syear,'(I3)') (jyear-1900)

			END IF

			WRITE(smonth,'(I2)') month
			IF (month/10.EQ.0) smonth = '0'//smonth(2:2)

			WRITE(sday,'(I2)') jday
			IF (jday/10.EQ.0) sday = '0'//sday(2:2)

			WRITE(shh,'(I2)') jhh
			IF (jhh.EQ.0) shh = '00'
			IF (jhh/10.EQ.0) shh = '0'//shh(2:2)

			WRITE(smm,'(I2)') jmm
			IF (jmm.EQ.0) smm = '00'
			IF (jmm/10.EQ.0) smm = '0'//smm(2:2)

			WRITE(sms,'(I4)') ms
			IF (ms.EQ.0) THEN
				sms = '0000'
			ELSE IF (ms/10.EQ.0) THEN
				sms = '000'//sms(4:4)
			ELSE IF (ms/100.EQ.0) THEN
				sms	= '00'//sms(3:4)
			ELSE IF (ms/1000.EQ.0) THEN
				sms	= '0'//sms(2:4)
			END IF

			jfilenum	= jhh*10**6 + jmm*10**4 + ms

			filenum		= shh//smm//sms 

			IF (jyear.LT.2000) THEN

				WRITE(8,223) s2year,smonth,sday,filenum,jday,month,
     &						jyear											! corrected new name of affected file
  223 FORMAT(A2,2('_',A2),'	',A8,2('	',I2),'	',I4)
  
				OPEN(2,	file	= 'c:\GEM\2300mhz\dat\'//s2year
     &				//'_'//smonth//'_'//sday//filenum//'.dat.', 
     &				status  = 'unknown',
     &				err     = 200);

				tag	= s2year//'_'//smonth//'_'//sday//filenum   
			
			ELSE
			
				WRITE(8,224) syear,smonth,sday,filenum,jday,month,
     &						jyear											! corrected new name of affected file
  224 FORMAT(A3,2('_',A2),'	',A8,2('	',I2),'	',I4)
  
				OPEN(2,	file	= 'c:\GEM\2300mhz\dat\'//syear
     &				//'_'//smonth//'_'//sday//filenum//'.dat.', 
     &				status  = 'unknown',
     &				err     = 200);

				tag	= syear//'_'//smonth//'_'//sday//filenum 
				
			END IF				                     						                   

		END IF				

		PRINT *, ' Doing file : ',tag 
		PRINT *,' '
		                     

		IF ((ll.EQ.1).OR.((jfini(ll)-jffin(ll-1)).NE.1)) THEN				! files that start a sequence of contiguous files
		
			timefin		=	timeini - 0.56002D0
			timefra		=	(jffin(ll)-jfini(ll)+1)*0.56002D0
			timediff	=	0.D0                        

C			tioff	=	0.D0
C
C			PRINT *,' ADD offset to 1st timestamp?'
C			PRINT *,'	1	->	Yes'
C			PRINT *,'	2	->	No'
C			READ(*,*) jof
C			PRINT *,' '
C			IF (jof.EQ.1) THEN
C
C				PRINT *,' ENTER offset (in s):'
C				READ(*,*) tioff
C
				timeini	= timeini + tioff
C	
C		END IF

			mdtaux	=	mdtime

			mdtime	=	tioff

			lm	=	1

			dtime(lm)	=	tioff

C			WRITE(*,2) timeini
			WRITE(6,4)	ll,tag,timeini,timediff,(timefin + .56002D0),
     &					timefra,(NINT(timefra/0.56002D0)),dtime(lm)

		ELSE

			timediff = timeini - timeaux									! timeaux : timeini of previous file

			IF (timediff.LT.0.D0) timediff = 86400.D0 + timediff

C			IF (ll.EQ.2) timefra = timefra + 0.56002						! timefra is referenced to last frame in the previous file

			frametime	= timefin + j99b*.56002D0							! time of 1st frame in file

			jframespan	= jffin(ll) - jfini(ll) + j99a + 1					! difference between 1st frames

			timefra		= jframespan*0.56002D0

			IF (DABS(timeini - frametime).LT.0.56002D0) THEN

				lm		= lm + 1

				dtime(lm)	= timeini - frametime							! GPS synchronized files should have dtime = 0
			
				mdtime	= mdtime + dtime(lm)

			END IF

C			WRITE(*,3)	timediff,frametime,timefra,lm,dtime,mdtime			! timefin is the time of the last frame in the previous file
			WRITE(6,4)	ll,tag,timeini,timediff,frametime,
     &					timefra,jframespan,dtime(lm)

		END IF

    2 FORMAT(	'  Timestamp from file name  : ',F8.2,' sec')
    3 FORMAT(	'  Time stamp difference : ',F8.2,' sec'/
     &		'  Time of 1st frame in file : ',F8.2,' sec'/
     &        '  Difference between 1st frames : ',F8.2,' sec'/
     &		'  offset index: ',I2,' increment: ',F8.5,' sec'/
     &		'  accumulated offset: ',F8.5,' sec')
    4 FORMAT(1x,I2,1x,A17,4(2x,F8.2),3x,I5,2x,F9.6)
          
		k = 0
	
		DO WHILE (.NOT.EOF(1))												! begin of loop over individual file

			k = k+1

			SELECT CASE (jobs)

				CASE(1)
		
			read(1,*,err=100) jframe, x, jelevat, jsignal, jazimut,			! 1994 Bishop data with Odetics
     &			jT1, jT2, jT3, jT4, j, jHI, j, jTns, jVns, jHII				! does not apply to 2.3 GHz (maybe in Colombia)
     
				CASE(2)     		 
          
			read(1,*,err=100) jframe, jelevat, jsignal, jazimut, 
     &				jT1, jT2, jT3, jT4, j, jHI, j, jTns, jVns, j			! 1998 observations

				CASE(3)

			read(1,*,err=100) jframe, jelevat, jsignal, j, 
     &				jT1, jT2, jT3, jT4, j, jHI, j, jTns, jVns, jazimut		! 1999 observations

			END SELECT

 
      ! pointing direction :::::::::::::::::::::::::::::::::::::::::
          
			elevat  = (du_V(jelevat) - V_180)/Vmaxel*encoff					! in degrees from zenith

			azimut  = du_V(jazimut)/Vmaxaz*encoff + azoff					! in degrees

      ! radiometer signal ::::::::::::::::::::::::::::::::::::::::::
          
			signal  = du_V(jsignal)											! in Volts    
          
      ! temperature sensors readouts in Celsius ::::::::::::::::::::
                     
			T1      = du_V(jT1)*10.0	! @ HEMT amplifier::::::::::::		                                         
			T2      = du_V(jT2)*10.0    ! @ detector diode & filter:::		! displays the best stability 
			T3      = du_V(jT3)*10.0    ! @ ambient at feed ::::::::::
			T4      = du_V(jT4)*10.0    ! @ electronics box ::::::::::
			
			Temp(k,1)	=	T1												! For housekeeping statistics.
			Temp(k,2)	=	T2
			Temp(k,3)	=	T3
			Temp(k,4)	=	T4
          	
      ! noise source :::::::::::::::::::::::::::::::::::::::::::::::
          
		    Tns	= 333.0-du_V(jTns)*46.29									! temperature in Celsius  
          
			Vns	= du_V(jVns)												! broadcast signal in V 
          
      ! impulse proportional counter :::::::::::::::::::::::::::::::
C                                                           
C			RFI	= du_V(jRFI)												! only 408 and 1465 ::::::::::
C		    
      ! heaters ::::::::::::::::::::::::::::::::::::::::::::::::::::
                  
			HI	= du_V(jHI)													! heater voltage in V ::::::::
		 
      ! time stamp assignment (exact value depends on filelist):::::
	! filelist4 : 98_08_10 to 98_08_12 uses 0.560020 sec :::::::::
                                                  
			IF ((ll.EQ.1).OR.((jfini(ll)-jffin(ll-1)).NE.1)) THEN			! files that start a sequence of contiguous files
		
				timefra	=	(k - 1)*.560020D0

				time	=	timeini + timefra
			
			ELSE 
		    
				timefra	=	(jframe - jfin)*.560020D0

				time    =	timefin + timefra
			                                    
			END IF
					    
			IF (time.GE.86400.) time = time - 86400.
			
	! Data output option :::::::::::::::::::::::::::::::::::::::::
      ! The 1st write statement delivers the full data conversion
      ! The 2nd is sufficient for calibration files.
	  ! The 3rd is used for calibrating the azimuth
	  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

C			SELECT CASE (jobs)
C
C				CASE (1,2)
C
C					WRITE(2,60) jframe, time, elevat, signal, azimut, 		! full data output
C     &							T1, T2, T3, T4, HI, Tns, Vns     
C				CASE (3)
C     
					WRITE(2,70) jframe, time, elevat, signal, azimut,		
     &							T1, T2, T3, T4, HI, Tns, Vns             
C                
C          write(2,80) jframe, time, signal, T1, T2, T3, T4, T5, HI,
C          write(2,90) jframe, time, elevat, azimut
C
C			END SELECT

		END DO																! end of individual loop

  100		CLOSE(1)
		CLOSE(2)
	
		timefin	=	time
		jfin	=	jframe
		timeaux	=	timeini
      
   60 format(1x,I6,1x,F8.2,1x,F7.2,1x,F9.6,1x,F10.6,4(1x,F7.4),
     &	   1x,F6.3,1x,F8.4,1x,F6.3)

   70 format(1x,I6,1x,F8.2,1x,F7.2,1x,F9.6,1x,F10.6,4(1x,F7.4),
     &	   1x,F6.3,1x,F8.4,1x,F6.3)

   80 format(1x,I5,3x,F8.2,7(3x,F10.6))
        
   90 format(1x,I6,3x,F8.2,3(F12.6))
        
C		WRITE(*,150) k         
C  150 FORMAT('  Data points in this file : ', I5)
C		PRINT *,' '

		IF (mit.EQ.1) THEN													! Temp stats are written only once.
	
		DO jj = 1,4
			tempavg(jj,1)	=	0.D0
			tempavg(jj,2)	=	0.D0
		END DO              

		DO jj = 1,4

			DO i = 1,k
				tempavg(jj,1)	=	tempavg(jj,1) + temp(i,jj)
				tem(i)			=	temp(i,jj)
			END DO
			
			CALL hpsortdb(k,tem)
			delta(jj)	=	(tem(k) - tem(1))/2

		END DO

		DO jj = 1,4	
			tempavg(jj,1)	=	tempavg(jj,1)/k 	
		END DO

		DO jj = 1,4
			DO i = 1,k
				tempavg(jj,2)	=	tempavg(jj,2) 
     &							+	(temp(i,jj) - tempavg(jj,1))**2
			END DO
		END DO

		DO jj = 1,4	

			tempavg(jj,2)	=	SQRT(tempavg(jj,2)/(k-1))
     &							*3600/(k*0.56002D0)	! 1-sigma/hour 	
			WRITE(*,170) tempavg(jj,1),tempavg(jj,2),delta(jj)
  170	FORMAT(' Mean T : ',F8.4,' +/- ',F7.4,'   Delta : ',F7.4)		

		END DO

		WRITE(3,180)(tempavg(jj,1),tempavg(jj,2),delta(jj), jj = 1,4),
     &				(3600/(k*0.56002D0)),tag
  180	FORMAT(4(F8.4,2(F8.4)),1x,F4.2,A16)

		END IF																! End of temp stats

C		PAUSE ' Press ENTER to process next file'

	END DO																	! End of files loop

	IF (jlist.EQ.1) THEN

		mdtime	=	1.D0

	ELSE

		mdtime	=	mdtime/(lm*1.D0)

		sdtime	=	0.D0

		DO i = 1,lm

			sdtime	= sdtime + (dtime(i) - mdtime)**2.D0

		END DO

		sdtime	= DSQRT(sdtime/(lm - 1))

		OPEN(1, file	= 'c:\GEM\2300mhz\txt\'//filelist//'.txt',			! reset fullnames for new time offset 
     &			status  = 'old',
     &			err     = 200)

		l	=	0

		DO WHILE (.NOT.EOF(1))

			l	=	l	+	1
			
			READ(1,'(A)') fullname(l)

		END DO

		CLOSE(1)

	END IF

	PRINT *,' '

	WRITE(*,182) mit, tioff, mdtime, sdtime
  182 FORMAT('           iteration : ',I2,'  =>  offset : ',F8.5,' sec'/
     &       ' < time difference > : ',F8.5,' sec ',' +/- ',F7.5,' sec')
	WRITE(9,183) tioff, mdtime, sdtime, mit 
  183 FORMAT(2(1x,F8.5),1x,F7.5,1x,I2)

		tioff	=	tioff + mdtime											! new time offset

	END DO																	! end of offset loop 

	CLOSE(3)
	CLOSE(7)
	CLOSE(8)
	CLOSE(9)

	lj	=	l

  185	IF (lj.GT.1) THEN

		OPEN(4,	file = 'c:\GEM\2300mhz\Temps\Tem'//filelist//'.dat',		
     &			status  = 'old',
     &			err		= 200);

		READ(4,*) 															! header line

		DO l = 1,lj

			READ(4,*) x,sT1(l),dT1(l),x,sT2(l),dT2(l),x,sT3(l),dT3(l),
     &				  x,sT4(l),dT4(l)

		END DO

		CALL hpsortdb(lj,sT1)
		CALL hpsortdb(lj,sT2)
		CALL hpsortdb(lj,sT3)
		CALL hpsortdb(lj,sT4)
		CALL hpsortdb(lj,dT1)
		CALL hpsortdb(lj,dT2)
		CALL hpsortdb(lj,dT3)
		CALL hpsortdb(lj,dT4)

		PRINT *,' '
		PRINT *,'                                1-sigma & Half Range'
		PRINT *,'                                ===================='
		PRINT *,' '

		PRINT *,'             T1            T2            T3        
     &             T4'
		PRINT *,' '


		WRITE(*,190) sT1(1),dT1(1),sT2(1),dT2(1),sT3(1),dT3(1),
     &				 sT4(1),dT4(1),
     &				 sT1(lj),dT1(lj),sT2(lj),
     &				 dT2(lj),sT3(lj),dT3(lj),
     &				 sT4(lj),dT4(lj)
    
    	PRINT *,' '
  190 FORMAT(' MIN : ',8(F7.4),//,' MAX : ',8(F7.4))
	WRITE(*,191) filelist
  191 FORMAT('    sigma cut in T2      N        in      =-> ',A9,' <-=')
  	PRINT *,' '

		k	=	0
		m	=	0
		DO WHILE (k.NE.l)

			m	=	m	+	1

			k	=	k	+	1

			DO WHILE ((sT2(k).LT.(m*0.0117)).AND.(k.LE.l))
			
				k	=	k	+	1
				
			END DO

			k	=	k	-	1
			
			IF (k.GT.kaux) THEN 

				WRITE(*,197) m,(m*0.0117),k
  197 FORMAT('      <',I2,' : ',F6.4,'       ',I2)

			END IF

			kaux	=	k

		END DO

		PRINT *,' '

		CLOSE(4)

C	PAUSE ' continue ...'
	
      END IF

      CLOSE(6)
      
  200 stop ' Conversion done'
      end
      
      function du_V(jdu)        
          jbits   = 2**15
          du_V    = 10.0*float(jdu - jbits)/float(jbits)
      END

	INCLUDE 'hpsortdb.for'