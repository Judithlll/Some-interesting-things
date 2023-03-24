# 天体测量学入门  实践二(NOVAS)





## 【题目】

![image-20220428150730162](C:\Users\LZX\AppData\Roaming\Typora\typora-user-images\image-20220428150730162.png)

![image-20220428150802368](C:\Users\LZX\AppData\Roaming\Typora\typora-user-images\image-20220428150802368.png)



## 【实践原理】（使用linux系统）

### 1.NOVAS

#### 利用NOVAS的程序包（下载地址：https://ascl.net/）仿写程序进行天体视位置和站心位置的计算。

#### 在下载的程序包中本实验主要用到了：

#### NOVAS_F3.1_Guide

#### NOVAS_F3.1_solsys2.f

#### NOVAS_F3.1.f

#### example.f

#### 按照NOVAS_F3.1)Guide中对各种子程序的描述以及example.f中的实例程序，本实验用到以下几个子程序：

```fortran
!     获得高精度位置
      HIACC 
!     获得某颗星的视位置
	  APSTAR ( TJD, N, RAI, DECI, PMRA, PMDEC, PARLAX, RADVEL, RA, DEC )
!	  获得某颗星的站心位置
	  TPSTAR ( UJD, GLON, GLAT, HT, RA, DEC)
!	  获得太阳系内某天体的视位置
	  APPLAN ( TJD, L, N, RA, DEC, DIS )
!	  获得太阳系内某天体的站心位置
	  TPPLAN ( UJD, GLON, GLAT, HT, RA, DEC, DIS )
!	  获得儒略日时间
	  JULDAT ( IYEAR, MONTH, IDAY, HOUR,   UTCJD )
!	  将赤道坐标转换为地平坐标
	  ZDAZ ( UJD, XP, YP, GLON, GLAT, HT, RA, DEC, IREFR, ZD, AZ, RAR, DECR )
```

#### 另外，NOVAS_F3.1_solsys2.f需要利用JPL历表，所以方便起见，把上次实验中生成的jplsub.f复制到NOVAS程序包下。

### 2.时间生成

#### 由于上述各子程序所需的时间不在同一时间系统下，所以在计算前先将儒略日转换成本实验所需的各时间系统下的时间，代码实现如下：

```fortran
!     ESTABLISH TIME ARGUMENTS
      CALL JULDAT  ( IYEAR, MONTH, IDAY, HOUR,   UTCJD )
      TTJD = UTCJD + ( LEAPS + 32.184D0 ) / 86400.D0
      UT1JD = UTCJD + UT1UTC / 86400.D0
      DELTAT = 32.184D0 + LEAPS - UT1UTC
      CALL SETDT ( DELTAT )
      WRITE ( *, * )
      WRITE ( *, * ) 'TT and UT1 Julian Dates and Delta-T:'
      WRITE ( *, 9040 ) TTJD, UT1JD, DELTAT
```

### 3.查询极移和UT1-UTC的值

#### 在IERS网站上进行查询（https://datacenter.iers.org/data/latestVersion/6_BULLETIN_A_V2013_016.txt)，查询得到以下数据：

```
                    MJD      x(arcsec)   y(arcsec)   UT1-UTC(sec)     
       2022  5  1  59700       0.0893      0.4713     -0.09766
```



## 【实践过程】

### 1.将所需的程序包放到同一文件夹下。

### 2.题目一：

#### 1）将题目数据赋给相关变量

#### 2）利用子程序得到视位置和站心位置

#### 3）将得到的数据转换为时分秒或度分秒的形式

#### 4）将站心位置转换为地平坐标

### 总代码如下：

```fortran
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      IMPLICIT INTEGER ( I-N )
      
      DIMENSION STAR(6), OBSERV(6), SKYPOS(7),
     .     POS(3), VEL(3), POSE(3), VTER(3), VCEL(3),
     .     SSS(3), VALUES(500) 
     
      CHARACTER NAMES(500)*6, TXTEPH*158, FILNAM*80

      INTEGER DENUM

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           ) 

      DATA STAR, OBSERV / 12 * 0.D0 /

      DATA IYEAR, MONTH, IDAY, HOUR, LEAPS, UT1UTC, XP, YP /
     .     2022, 5, 1, 23.2619D0, 33, 0.0976630D0, 
     .     +0.089287D0, +0.47129D0/
     
      DATA GLON, GLAT, HT / -70.7365D0, -30.2408D0, 2378.D0 /
      
      DATA LU, FILNAM / 20, 'SS_EPHEM.TXT' /


!     FORMAT STATEMENTS FOR OUTPUT

 9010 FORMAT ( 1X, 'JPL ephemeris DE', I3, ' open. Start JD = ', F10.2, 
     .         '  End JD = ', F10.2 ) 
 9020 FORMAT ( A158 )
 9030 FORMAT ( 1X, 3 ( F15.10, 8X ) )
 9040 FORMAT ( 1X, 2 ( F15.6, 8X ), 1F16.11 )
 9050 FORMAT ( 1X, 2 ( F15.10, 8X ), 1F15.12 )
 9060 FORMAT ( 1X, 2 ( F16.11, 8X ), 1F15.10 )
 9070 FORMAT ( 2X, I3,A1,1X,I2,A1,1X,F6.3,A1)
 
!     CHECK WHICH JPL EPHEMERIS IS AVAILABLE
!     REMOVE THIS BLOCK FOR USE WITH SOLSYS 1
      CALL CONST ( NAMES, VALUES, SSS, N )
      DENUM = VALUES ( 1 )
      EPHBEG = SSS ( 1 )
      EPHEND = SSS ( 2 )
      WRITE ( *, 9010 )  DENUM, EPHBEG, EPHEND  
      
!     VERIFY THAT TEXT EPHEMERIS FILE IS AVAILABLE
!     UNCOMMENT THIS BLOCK FOR USE WITH SOLSYS 1
*     OPEN ( UNIT=LU, FILE=FILNAM, STATUS='UNKNOWN' )
*     READ ( LU, 9020 ) TXTEPH
*     WRITE ( *, 9020 ) TXTEPH
*     CLOSE ( LU )

  
!     SETUP CALLS
!     HIGH ACCURACY AND EQUINOX MODE ARE NOVAS DEFAULTS.
      CALL HIACC
      CALL EQINOX
      IEARTH = IDSS ( 'EARTH' )
    
!     WRITE OUT ASSUMED LONGITUDE, LATITUDE, HEIGHT (ITRS = WGS-84)
      WRITE ( *, * )
      WRITE ( *, * ) 'Geodetic location:'
      WRITE ( *, 9030 ) GLON, GLAT, HT
    
!     ESTABLISH TIME ARGUMENTS
      CALL JULDAT  ( IYEAR, MONTH, IDAY, HOUR,   UTCJD )
      TTJD = UTCJD + ( LEAPS + 32.184D0 ) / 86400.D0
      UT1JD = UTCJD + UT1UTC / 86400.D0
      DELTAT = 32.184D0 + LEAPS - UT1UTC
      CALL SETDT ( DELTAT )
      WRITE ( *, * )
      WRITE ( *, * ) 'TT and UT1 Julian Dates and Delta-T:'
      WRITE ( *, 9040 ) TTJD, UT1JD, DELTAT
      
!     APPARENT AND TOPOCENTRIC PLACE OF STAR FK6 1307 = GRB 1830
      RA2000 = 14.5713366194D0
      DC2000 = -12.5195545833D0
      PMRA   = -354.45D0
      PMDEC  = 595.35D0
      PARX   =  164.99D0
      RV     =  0D0
      CALL APSTAR ( TTJD, IEARTH, RA2000, DC2000, PMRA, PMDEC, PARX, RV,
     .     RA, DEC )
      CALL TPSTAR ( UT1JD, GLON, GLAT, HT,   RAT, DECT )

!     CONVERT RA TO hms, DEC to dfs
      ih=int(RA)
      im=int((RA-ih)*60)
      as=((RA-ih)*60-im)*60
      
      idd=int(DEC)
      idm=abs(int((DEC-idd)*60))
      ads=abs(((DEC-idd)*60-int((DEC-idd)*60))*60)
      WRITE(*,*)
      WRITE(*,*)"geocentric positions:"
      WRITE(*,9070,advance="no")ih,"h",im,"m",as,"s"
      WRITE(*,9070)idd,"d",idm,"m",ads,"s"

      ih=int(RAT)
      im=int((RAT-ih)*60)
      as=((RAT-ih)*60-im)*60
      
      idd=int(DECT)
      idm=abs(int((DECT-idd)*60))
      ads=abs(((DECT-idd)*60-int((DECT-idd)*60))*60)
      WRITE(*,*)
      WRITE(*,*)"topocentric positions:"
      WRITE(*,9070,advance="no")ih,"h",im,"m",as,"s"
      WRITE(*,9070)idd,"d",idm,"m",ads,"s"
      CALL ZDAZ ( UT1JD, XP, YP, GLON, GLAT, HT, RAT, DECT, 1,
     .              ZD, AZ, RAR, DECR )
      WRITE ( *, * )
      WRITE ( *, * ) 'zenith distance and azimuth:'
      WRITE ( *, 9030 ) ZD, AZ
      end
```

### 输出：

```fortran
 JPL ephemeris DE405 open. Start JD = 2451536.50  End JD = 2525008.50

 Geodetic location:
  -70.7365000000         -30.2408000000        2378.0000000000

 TT and UT1 Julian Dates and Delta-T:
  2459701.470000         2459701.469247          65.08633700000

 geocentric positions:
   14h 35m 29.388s  -12d 36m 51.449s

 topocentric positions:
   14h 35m 29.391s  -12d 36m 51.389s

 zenith distance and azimuth:
   75.7999596585          96.5295969608
```

#### 可以看到地心视位置（ geocentric positions）为：

$$
14h 35m 29.388s    ，    -12d 36m 51.449s
$$

#### 站心位置（ topocentric positions）为：

$$
   14h 35m 29.391s，  -12d 36m 51.389s
$$

#### 天顶距和方位角（zenith distance and azimuth）为：

$$
 75.7999596585   ，       96.5295969608
$$

### 3.题目二

#### 1）将观测站位置改为兴隆观测站：（官网查询）

```fortran
      DATA GLON, GLAT, HT / 117.577222D0, 40.395833D0, 900.D0 /
```

#### 2）利用循环得到月球和其他七大行星的视位置和站心坐标

#### 3）将得到的数据转换为时分秒或度分秒的形式

#### 4）将赤道坐标转换为地平坐标

#### 5）利用天顶距判断能否观测到这颗行星

### 代码如下：

```fortran
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      IMPLICIT INTEGER ( I-N )
      
      DIMENSION STAR(6), OBSERV(6), SKYPOS(7),
     .     POS(3), VEL(3), POSE(3), VTER(3), VCEL(3),
     .     SSS(3), VALUES(500),RA(8),DEC(8),RAT(8), 
     .     DECT(8),ZD(8),AZ(8)
     
      CHARACTER NAMES(500)*6, TXTEPH*158, FILNAM*80
      CHARACTER*3 ID(8)

      INTEGER DENUM

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           ) 

      DATA STAR, OBSERV / 12 * 0.D0 /

      DATA IYEAR, MONTH, IDAY, HOUR, LEAPS, UT1UTC, XP, YP /
     .     2022, 5, 1, 23.2619D0, 33, 0.0976630D0, 
     .     +0.089287D0, +0.47129D0/
     
      DATA GLON, GLAT, HT / 117.577222D0, 40.395833D0, 900.D0 /
      
      DATA LU, FILNAM / 20, 'SS_EPHEM.TXT' /


!     FORMAT STATEMENTS FOR OUTPUT

 9010 FORMAT ( 1X, 'JPL ephemeris DE', I3, ' open. Start JD = ', F10.2, 
     .         '  End JD = ', F10.2 ) 
 9020 FORMAT ( A158 )
 9030 FORMAT ( 1X, 3 ( F15.10, 8X ) )
 9040 FORMAT ( 1X, 2 ( F15.6, 8X ), 1F16.11 )
 9050 FORMAT ( 1X, 2 ( F15.10, 8X ), 1F15.12 )
 9060 FORMAT ( 1X, 2 ( F16.11, 8X ), 1F15.10 )
 9070 FORMAT ( 2X, I3,A1,1X,I2,A1,1X,F6.3,A1)
 
!     CHECK WHICH JPL EPHEMERIS IS AVAILABLE
!     REMOVE THIS BLOCK FOR USE WITH SOLSYS 1
      CALL CONST ( NAMES, VALUES, SSS, N )
      DENUM = VALUES ( 1 )
      EPHBEG = SSS ( 1 )
      EPHEND = SSS ( 2 )
      WRITE ( *, 9010 )  DENUM, EPHBEG, EPHEND  
      
!     VERIFY THAT TEXT EPHEMER
!     IS FILE IS AVAILABLE
!     UNCOMMENT THIS BLOCK FOR USE WITH SOLSYS 1
*     OPEN ( UNIT=LU, FILE=FILNAM, STATUS='UNKNOWN' )
*     READ ( LU, 9020 ) TXTEPH
*     WRITE ( *, 9020 ) TXTEPH
*     CLOSE ( LU )

!     SETUP CALLS
!     HIGH ACCURACY AND EQUINOX MODE ARE NOVAS DEFAULTS.
      CALL HIACC
      CALL EQINOX
      IEARTH = IDSS ( 'EARTH' )
    
!     WRITE OUT ASSUMED LONGITUDE, LATITUDE, HEIGHT (ITRS = WGS-84)
      WRITE ( *, * )
      WRITE ( *, * ) 'Geodetic location:'
      WRITE ( *, 9030 ) GLON, GLAT, HT
    
!     ESTABLISH TIME ARGUMENTS
      CALL JULDAT  ( IYEAR, MONTH, IDAY, HOUR,   UTCJD )
      TTJD = UTCJD + ( LEAPS + 32.184D0 ) / 86400.D0
      UT1JD = UTCJD + UT1UTC / 86400.D0
      DELTAT = 32.184D0 + LEAPS - UT1UTC
      CALL SETDT ( DELTAT )
      WRITE ( *, * )
      WRITE ( *, * ) 'TT and UT1 Julian Dates and Delta-T:'
      WRITE ( *, 9040 ) TTJD, UT1JD, DELTAT

!     APPARENT AND TOPOCENTRIC PLACE OF THE MOON
      

      ID=(/'MOO','MER','VEN','MAR','JUP','SAT','URA','NEP'/)
      l=1
      DO 40 WHILE(l.LE.8)
      NAME=IDSS(ID(l))
      WRITE(*,*)
      WRITE(*,*)"********************",ID(l),"********************"
      CALL APPLAN ( TTJD, NAME, IEARTH,   RA(l), DEC(l), DIS )
      CALL TPPLAN ( UT1JD, GLON, GLAT, HT,   RAT(l), DECT(l), DIST )
      WRITE ( *, * ) 'The geocentric positions are:'
      ih=int(RA(l))
      im=abs(int((RA(l)-ih)*60))
      as=abs(((RA(l)-ih)*60-int((RA(l)-ih)*60))*60)
   
      idd=int(DEC(l))
      idm=abs(int((DEC(l)-idd)*60))
      ads=abs(((DEC(l)-idd)*60-int((DEC(l)-idd)*60))*60)

      WRITE(*,9070,advance="no")ih,"h",im,"m",as,"s"
      WRITE(*,9070)idd,"d",idm,"m",ads,"s"

      WRITE ( *, * )
      WRITE ( *, * ) 'The topocentric positions are:'
      ih=int(RAT(l))
      im=abs(int((RAT(l)-ih)*60))
      as=abs(((RAT(l)-ih)*60-int((RAT(l)-ih)*60))*60)

      idd=int(DECT(l))
      idm=abs(int((DECT(l)-idd)*60))
      ads=abs(((DECT(l)-idd)*60-int((DECT(l)-idd)*60))*60)
      WRITE(*,9070,advance="no")ih,"h",im,"m",as,"s"
      WRITE(*,9070)idd,"d",idm,"m",ads,"s"

      WRITE ( *, * )
      WRITE(*,*)"the zenith distances and azimuths are:"
      CALL ZDAZ ( UT1JD, XP, YP, GLON, GLAT, HT, RAT(l), DECT(l), 1,
     .              ZD(l), AZ(l), RAR, DECR )
      WRITE ( *, 9030 ) ZD(l), AZ(l)
      l=l+1
40    CONTINUE
      WRITE(*,*)
      i=1
      DO 10 WHILE(i.LE.8)
      IF(ZD(i).GT.90)THEN
      WRITE(*,*)"      ",ID(i)," is not visible"
      ELSE
      WRITE(*,*)"      ",ID(i)," is visible"
      END IF
      i=i+1
10    CONTINUE

      end
```

### 输出：

```fortran
 JPL ephemeris DE405 open. Start JD = 2451536.50  End JD = 2525008.50

 Geodetic location:
  117.5772220000          40.3958330000         900.0000000000

 TT and UT1 Julian Dates and Delta-T:
  2459701.470000         2459701.469247          65.08633700000

 ********************MOO********************
 The geocentric positions are:
    3h 27m  4.027s   18d 58m  6.709s

 The topocentric positions are:
    3h 30m  0.583s   18d 25m 32.217s

 the zenith distances and azimuths are:
   75.2586970548          78.1104728384

 ********************MER********************
 The geocentric positions are:
    3h 56m  0.248s   23d  7m 59.029s

 The topocentric positions are:
    3h 56m  0.859s   23d  7m 52.215s

 the zenith distances and azimuths are:
   77.0894766407          70.4276625427

 ********************VEN********************
 The geocentric positions are:
   23h 59m 12.365s   -1d 33m 33.797s

 The topocentric positions are:
   23h 59m 12.623s   -1d 33m 39.577s

 the zenith distances and azimuths are:
   51.8364350890         135.5424440196

 ********************MAR********************
 The geocentric positions are:
   22h 58m 41.787s   -8d 12m 13.490s

 The topocentric positions are:
   22h 58m 41.890s   -8d 12m 17.496s

 the zenith distances and azimuths are:
   51.4351128480         156.5747436424

 ********************JUP********************
 The geocentric positions are:
   23h 55m  3.671s   -1d 42m 15.049s

 The topocentric positions are:
   23h 55m  3.727s   -1d 42m 16.077s

 the zenith distances and azimuths are:
   51.4114371073         136.7618801425

 ********************SAT********************
 The geocentric positions are:
   21h 47m 54.846s  -14d 22m 42.980s

 The topocentric positions are:
   21h 47m 54.863s  -14d 22m 43.690s

 the zenith distances and azimuths are:
   54.7562648254         179.2728759145

 ********************URA********************
 The geocentric positions are:
    2h 49m  0.121s   15d 52m 25.769s

 The topocentric positions are:
    2h 49m  0.147s   15d 52m 25.462s

 the zenith distances and azimuths are:
   69.1125158338          86.5465027745

 ********************NEP********************
 The geocentric positions are:
   23h 41m 56.007s   -3d 11m 41.456s

 The topocentric positions are:
   23h 41m 56.028s   -3d 11m 41.645s

 the zenith distances and azimuths are:
   51.0651671079         141.3589250536

       MOO is visible
       MER is visible
       VEN is visible
       MAR is visible
       JUP is visible
       SAT is visible
       URA is visible
       NEP is visible
```

#### 以上分别是各星的位置，且在这一天题述各星均可见。



## 【思考】

### 就一点，应该把转换为时分秒、度分秒的程序写成一个子程序以便调用，不然程序就像上面一样有点冗长。
