:: script generated by shelle v1.24
@echo OFF
:: for debugging
:: set DEBUG=echo
:: for running
set DEBUG=
setlocal enableextensions enabledelayedexpansion
if NOT ERRORLEVEL 0 echo extensions not enabled, script may fail
if NOT DEFINED ELLEPATH (
echo Please set your environment variable ELLEPATH
echo so that ELLEPATH\binwx is the location of the
echo Elle executables
GOTO END
)
cd %ELLEPATH%\..\experiments\experiment_12_Grooves\a
set BIN=%ELLEPATH%\binwx
set STARTFILE=res50.elle
%DEBUG% start %BIN%\elle_latte -i %STARTFILE% -u 0
