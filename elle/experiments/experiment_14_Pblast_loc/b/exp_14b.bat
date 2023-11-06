:: script generated by shelle v1.24
@echo OFF
setlocal enableextensions enabledelayedexpansion
if NOT ERRORLEVEL 0 echo extensions not enabled, script may fail
cd %ELLEPATH%\..\experiments\experiment_14_Pblast_loc\b
:: for debugging
:: set DEBUG=echo
:: for running
set DEBUG=
set DELETE=del /Q
set COPY=copy /y
set MOVE=rename
set LOCAL=%HOMEDRIVE%
set BIN=%ELLEPATH%\binwx
:: set BIN=..\..\..\elle\binwx
set STARTFILE=exp14b.elle
set T=0
set LASTTIME=30
set TIMESTEP=1
set GGSTAGES=003
set GSTAGES=000
set PBSTAGES=015
set USERNAME=myprocess
set USERSTAGES=000
set U1=
set U2=
set U3=
set U4=
set U5=
set U6=
set GBDIFFSTAGES=000
set GBDIFFKAPPA=0.0001
set LATTDIFFSTAGES=000
set EXCHSTAGES=000
set EXCHKAPPA=0.0001
set EXCHTEMP=0.0
set EXCHDIFFONLY=0
set EXPANDAREA=0
set EXPANDSPEED=10
set MANUELSPEED=0.00
set SAVEINTVL=1
set SAVEALL=15
set BASINFILE=exp14b.in
set ROOTNAME=exp_14b
set GV_SAVE=0
set PLOTAX=0
set LIST=1 2 3 13
set R=%ROOTNAME%
:: set the input for the first process
if %T% EQU 0 (
   if NOT EXIST %STARTFILE% GOTO END
   echo copy %STARTFILE% to tmp.elle
   %DEBUG% %COPY% %STARTFILE% tmp.elle
   echo copy %STARTFILE% to %R%.000.00.elle
   %DEBUG% %COPY% %STARTFILE% %R%.000.00.elle
   :: randomise quartz orientations
   for %%i in ( %LIST% ) do (
   if %%i EQU 7 (
   echo %BIN%\tidy -i tmp.elle -u 0 1 -s 1 -f 1 -n
   %DEBUG% start /w %BIN%\tidy -i tmp.elle -u 0 1 -s 1 -f 1 -n
   if NOT ERRORLEVEL 0 GOTO END
   echo %COPY% tidy001.elle tmp.elle
   %DEBUG% %COPY% tidy001.elle tmp.elle
   if NOT ERRORLEVEL 0 GOTO END
   )
   )
if NOT EXIST FD.sols mkdir FD.sols
if NOT EXIST FD.out mkdir FD.out
   if NOT ERRORLEVEL 0 GOTO END
)

:: loop through until time is up

:LOOP
:: update time variable
set /a T=%T%+%TIMESTEP%
set /a CNT=%T% %% %SAVEINTVL%
set /a LEFT=%LASTTIME%-%T%
set TSTR=%T%
if %T% LSS 100 set TSTR=0%T%
if %T% LSS 10 set TSTR=00%T%
:: save file for ghostview display
if %GV_SAVE% == 1 (
echo %BIN%\showelleps -i tmp.elle
%DEBUG% start /w %BIN%\showelleps -i tmp.elle
)

echo TIME=%TSTR%

:: run the loop
for %%i in ( %LIST% ) do (

:: create the polygon input from the elle file
if %%i EQU 1 (
echo %BIN%\elle2poly tmp.elle tmp.poly pb
%DEBUG% start /w %BIN%\elle2poly tmp.elle tmp.poly pb
)
:: run basil
if %%i EQU 2 (
echo %BIN%\basil %BASINFILE%
%DEBUG% start /w /MIN %BIN%\basil %BASINFILE%
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %COPY% FD.sols\basil FD.sols\%R%.%TSTR%
   %DEBUG% %COPY% FD.sols\basil FD.sols\%R%.%TSTR%
   echo %COPY% FD.out\basil.out FD.out\%R%.out.%TSTR%
   %DEBUG% %COPY% FD.out\basil.out FD.out\%R%.out.%TSTR%
   if NOT ERRORLEVEL 0 GOTO END
   )
)
:: convert to elle format
if %%i EQU 3 (
echo %BIN%\basil2elle -i FD.sols\basil -r 2 -o tmp2.elle
%DEBUG% start /w %BIN%\basil2elle -i FD.sols\basil -r 2 -o tmp2.elle
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% tmp2.elle tmp.elle
%DEBUG% %COPY% tmp2.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% tmp2.elle %R%.%TSTR%.03.elle
   %DEBUG% %MOVE% tmp2.elle %R%.%TSTR%.03.elle
   if NOT ERRORLEVEL 0 GOTO END
   )
)
:: let tbh happen
if %%i EQU 4 (
echo %BIN%\elle_tbh -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_tbh -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% tbh.calc001.elle tmp.elle
%DEBUG% %COPY% tbh.calc001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% tbh.calc001.elle %R%.%TSTR%.04.elle
   %DEBUG% %MOVE% tbh.calc001.elle %R%.%TSTR%.04.elle
   if %PLOTAX% EQU 1 (
   echo %BIN%\plotaxes -i %R%.%TSTR%.04.elle -s 1 -f 1 -n
   %DEBUG% start /w %BIN%\plotaxes -i %R%.%TSTR%.04.elle -s 1 -f 1 -n
   )
) ELSE (
   if %SAVEALL% EQU 4 %DEBUG% %MOVE% tbh.calc001.elle %R%.%TSTR%.04.elle
)
)
:: let grains split into subgrains
if %%i EQU 5 (
echo %BIN%\elle_split -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_split -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% split001.elle tmp.elle
%DEBUG% %COPY% split001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% split001.elle %R%.%TSTR%.05.elle
   %DEBUG% %MOVE% split001.elle %R%.%TSTR%.05.elle
) ELSE (
   if %SAVEALL% EQU 5 %DEBUG% %MOVE% split001.elle %R%.%TSTR%.05.elle
   )
)
:: let grains boundaries migrate by everything
if %%i EQU 6 (
set S=%GSTAGES%
echo %BIN%\elle_gbm -i tmp.elle -s !S! -f !S! -n
%DEBUG% start /w %BIN%\elle_gbm -i tmp.elle -s !S! -f !S! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% gbm!S!.elle tmp.elle
%DEBUG% %COPY% gbm!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% gbm!S!.elle %R%.%TSTR%.06.elle
   %DEBUG% %MOVE% gbm!S!.elle %R%.%TSTR%.06.elle
) ELSE (
   if %SAVEALL% EQU 6 %DEBUG% %MOVE% gbm!S!.elle %R%.%TSTR%.06.elle
   )
)
:: grow pblast
if %%i EQU 8 (
set S=%PBSTAGES%
echo %BIN%\elle_pblast -i tmp.elle -s !S! -f !S! -n
%DEBUG% start /w %BIN%\elle_pblast -i tmp.elle -s !S! -f !S! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% pblast!S!.elle tmp.elle
%DEBUG% %COPY% pblast!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% pblast!S!.elle %R%.%TSTR%.08.elle
   %DEBUG% %MOVE% pblast!S!.elle %R%.%TSTR%.08.elle
) ELSE (
   if %SAVEALL% EQU 8 %DEBUG% %MOVE% pblast!S!.elle %R%.%TSTR%.08.elle
   )
)
:: do oof calculation
if %%i EQU 9 (
echo %BIN%\elle2poly tmp.elle tmp.poly pb
%DEBUG% start /w %BIN%\elle2poly tmp.elle tmp.poly pb
echo %BIN%\triangle -pzAen tmp.elle.poly
%DEBUG% start /w %BIN%\triangle -pzAen tmp.elle.poly
echo %BIN%\triangle -pzAen tmp.elle.1.poly
%DEBUG% start /w %BIN%\triangle -pzAen tmp.elle.1.poly
echo %BIN%\triangle -pzrAPna0.000625 tmp.elle.2.poly
%DEBUG% start /w %BIN%\triangle -pzrAPna0.000625 tmp.elle.2.poly
%DEBUG% %DELETE% tmp.elle.*.edge tmp.elle.*.neigh tmp.elle.1* tmp.elle.2*
echo %BIN%\poly2goof -i tmp.elle -n
%DEBUG start /w %BIN%\poly2goof -i tmp.elle -n
if NOT ERRORLEVEL 0 GOTO END
echo %LOCAL%\oof\oof -grid elle.goof -text -quit
%DEBUG% %LOCAL%\oof\oof -grid elle.goof -text -quit
if NOT ERRORLEVEL 0 GOTO END
echo %BIN%\goof2elle -i tmp.elle -s 1 -f 1 -n
%DEBUG% %BIN%\goof2elle -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% goof2elle001.elle tmp.elle
%DEBUG% %COPY% goof2elle001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% goof2elle001.elle %R%.%TSTR%.09.elle
   %DEBUG% %MOVE% goof2elle001.elle %R%.%TSTR%.09.elle
) ELSE (
   if %SAVEALL% EQU 9 %DEBUG% %MOVE% goof2elle001.elle %R%.%TSTR%.09.elle
   )
)
:: let let nucleation/recrystallisation happen if dislocden is high enough
if %%i EQU 10 (
echo %BIN%\elle_disloc_rx -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_disloc_rx -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% disloc_rx001.elle tmp.elle
%DEBUG% %COPY% disloc_rx001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% disloc_rx001.elle %R%.%TSTR%.10.elle
   %DEBUG% %MOVE% disloc_rx001.elle %R%.%TSTR%.10.elle
) ELSE (
   if %SAVEALL% EQU 10 %DEBUG% %MOVE% disloc_rx001.elle %R%.%TSTR%.10.elle
   )
)
:: let rotation recrystallisation on only split flynns happen
if %%i EQU 11 (
echo %BIN%\elle_angle_rx -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_angle_rx -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% angle_rx001.elle tmp.elle
%DEBUG% %COPY% angle_rx001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% angle_rx001.elle %R%.%TSTR%.11.elle
   %DEBUG% %MOVE% angle_rx001.elle %R%.%TSTR%.11.elle
) ELSE (
   if %SAVEALL% EQU 11 %DEBUG% %MOVE% angle_rx001.elle %R%.%TSTR%.11.elle
   )
)
:: let grains recover
if %%i EQU 12 (
echo %BIN%\elle_recovery -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_recovery -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% recover001.elle tmp.elle
%DEBUG% %COPY% recover001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% recover001.elle %R%.%TSTR%.12.elle
   %DEBUG% %MOVE% recover001.elle %R%.%TSTR%.12.elle
) ELSE (
   if %SAVEALL% EQU 12 %DEBUG% %MOVE% recover001.elle %R%.%TSTR%.12.elle
   )
)
:: let repositioning happen
if %%i EQU 13 (
echo %BIN%\reposition -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\reposition -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% repos.elle tmp.elle
%DEBUG% %COPY% repos.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% repos.elle %R%.%TSTR%.13.elle
   %DEBUG% %MOVE% repos.elle %R%.%TSTR%.13.elle
) ELSE (
   if %SAVEALL% EQU 13 %DEBUG% %MOVE% repos.elle %R%.%TSTR%.13.elle
   )
)
:: do statistics
if %%i EQU 14 (
echo %BIN%\elle_stats -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_stats -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% stats001.elle tmp.elle
%DEBUG% %COPY% stats001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% stats001.elle %R%.%TSTR%.14.elle
   %DEBUG% %MOVE% stats001.elle %R%.%TSTR%.14.elle
) ELSE (
   if %SAVEALL% EQU 14 %DEBUG% %MOVE% stats001.elle %R%.%TSTR%.14.elle
   )
)
:: check for small angle vertices
if %%i EQU 15 (
echo %BIN%\checkangle -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\checkangle -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% checkangle001.elle tmp.elle
%DEBUG% %COPY% checkangle001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% checkangle001.elle %R%.%TSTR%.15.elle
   %DEBUG% %MOVE% checkangle001.elle %R%.%TSTR%.15.elle
) ELSE (
   if %SAVEALL% EQU 15 %DEBUG% %MOVE% checkangle001.elle %R%.%TSTR%.15.elle
   )
)
:: recalculate viscosity
if %%i EQU 16 (
echo %BIN%\elle_viscosity -i tmp.elle -s 1 -f 1 -u 1 -n
%DEBUG% start /w %BIN%\elle_viscosity -i tmp.elle -s 1 -f 1 -u 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% viscosity001.elle tmp.elle
%DEBUG% %COPY% viscosity001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% viscosity001.elle %R%.%TSTR%.16.elle
   %DEBUG% %MOVE% viscosity001.elle %R%.%TSTR%.16.elle
) ELSE (
   if %SAVEALL% EQU 16 %DEBUG% %MOVE% viscosity001.elle %R%.%TSTR%.16.elle
   )
)
:: let grains merge
if %%i EQU 17 (
echo %BIN%\elle_merge -i tmp.elle -s 1 -f 1 -n
%DEBUG% start /w %BIN%\elle_merge -i tmp.elle -s 1 -f 1 -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% merge001.elle tmp.elle
%DEBUG% %COPY% merge001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% merge001.elle %R%.%TSTR%.17.elle
   %DEBUG% %MOVE% merge001.elle %R%.%TSTR%.17.elle
) ELSE (
   if %SAVEALL% EQU 17 %DEBUG% %MOVE% merge001.elle %R%.%TSTR%.17.elle
   )
)
:: let grains boundaries migrate by surface energy
if %%i EQU 18 (
set S=%GGSTAGES%
echo %BIN%\elle_gg -i tmp.elle -s !S! -f !S! -n
%DEBUG% start /w %BIN%\elle_gg -i tmp.elle -s !S! -f !S! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% growth!S!.elle tmp.elle
%DEBUG% %COPY% growth!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% growth!S!.elle %R%.%TSTR%.18.elle
   %DEBUG% %MOVE% growth!S!.elle %R%.%TSTR%.18.elle
) ELSE (
   if %SAVEALL% EQU 18 %DEBUG% %MOVE% growth!S!.elle %R%.%TSTR%.18.elle
   )
)
:: grow pblast by expand algorithm
if %%i EQU 19 (
set S=%PBSTAGES%
set P=-u %EXPANDSPEED% 0 0 0 0 %EXPANDAREA%
echo %BIN%\elle_expand -i tmp.elle -s !S! -f !S! !P! -n
%DEBUG% start /w %BIN%\elle_expand -i tmp.elle -s !S! -f !S! !P! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% expand!S!.elle tmp.elle
%DEBUG% %COPY% expand!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% expand!S!.elle %R%.%TSTR%.19.elle
   %DEBUG% %MOVE% expand!S!.elle %R%.%TSTR%.19.elle
) ELSE (
   if %SAVEALL% EQU 19 %DEBUG% %MOVE% expand!S!.elle %R%.%TSTR%.19.elle
   )
)
:: gb diffussion
if %%i EQU 20 (
set S=%GBDIFFSTAGES%
set P=-u %GBDIFFKAPPA%
echo %BIN%\elle_gbdiff -i tmp.elle -s !S! -f !S! !P! -n
%DEBUG% start /w %BIN%\elle_gbdiff -i tmp.elle -s !S! -f !S! !P! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% gbdiff!S!.elle tmp.elle
%DEBUG% %COPY% gbdiff!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% gbdiff!S!.elle %R%.%TSTR%.20.elle
   %DEBUG% %MOVE% gbdiff!S!.elle %R%.%TSTR%.20.elle
) ELSE (
   if %SAVEALL% EQU 20 %DEBUG% %MOVE% gbdiff!S!.elle %R%.%TSTR%.20.elle
   )
)
:: lattice diffusion
if %%i EQU 21 (
set S=%LATTDIFFSTAGES%
echo %BIN%\elle_diff -i tmp.elle -s !S! -f !S! -n
%DEBUG% start /w %BIN%\elle_diff -i tmp.elle -s !S! -f !S! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% diff!S!.elle tmp.elle
%DEBUG% %COPY% diff!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% diff!S!.elle %R%.%TSTR%.21.elle
   %DEBUG% %MOVE% diff!S!.elle %R%.%TSTR%.21.elle
) ELSE (
   if %SAVEALL% EQU 21 %DEBUG% %MOVE% diff!S!.elle %R%.%TSTR%.21.elle
   )
)
:: let homogenoeus simple shear happen
if %%i EQU 22 (
set S=%MANUELSPEED%
echo %BIN%\elle_manuel -i tmp.elle -s 1 -f 1 -v !S! -n
%DEBUG% start /w %BIN%\elle_manuel -i tmp.elle -s 1 -f 1 -v !S! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% manuel001.elle tmp.elle
%DEBUG% %COPY% manuel001.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% manuel001.elle %R%.%TSTR%.22.elle
   %DEBUG% %MOVE% manuel001.elle %R%.%TSTR%.22.elle
) ELSE (
   if %SAVEALL% EQU 22 %DEBUG% %MOVE% manuel001.elle %R%.%TSTR%.22.elle
   )
)
:: exchange reaction
if %%i EQU 23 (
set S=%EXCHSTAGES%
set P=-u %EXCHKAPPA% 0 0 %EXCHDIFFONLY% %EXCHTEMP%
echo %BIN%\elle_exchange -i tmp.elle -s !S! -f !S! !P! -n
%DEBUG% start /w %BIN%\elle_exchange -i tmp.elle -s !S! -f !S! !P! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% exchange!S!.elle tmp.elle
%DEBUG% %COPY% exchange!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% exchange!S!.elle %R%.%TSTR%.23.elle
   %DEBUG% %MOVE% exchange!S!.elle %R%.%TSTR%.23.elle
) ELSE (
   if %SAVEALL% EQU 23 %DEBUG% %MOVE% exchange!S!.elle %R%.%TSTR%.23.elle
   )
)
:: userdefined algorithm
if %%i EQU 99 (
set S=%USERSTAGES%
set P=-u %U1% %U2% %U3% %U4% %U5% %U6%
echo %BIN%\myprocess -i tmp.elle -s !S! -f !S! !P! -n
%DEBUG% start /w %BIN%\myprocess -i tmp.elle -s !S! -f !S! !P! -n
if NOT ERRORLEVEL 0 GOTO END
echo %COPY% myprocess!S!.elle tmp.elle
%DEBUG% %COPY% myprocess!S!.elle tmp.elle
if NOT ERRORLEVEL 0 GOTO END
if %CNT% EQU 0 (
   echo %MOVE% myprocess!S!.elle %R%.%TSTR%.99.elle
   %DEBUG% %MOVE% myprocess!S!.elle %R%.%TSTR%.99.elle
) ELSE (
   if %SAVEALL% EQU 99 %DEBUG% %MOVE% myprocess!S!.elle %R%.%TSTR%.99.elle
   )
)
if NOT ERRORLEVEL 0 GOTO END
)
if NOT ERRORLEVEL 0 GOTO END
if %LEFT% NEQ 0 GOTO LOOP
:END
endlocal
