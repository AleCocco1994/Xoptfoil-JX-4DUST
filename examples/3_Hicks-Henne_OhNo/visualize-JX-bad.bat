@echo off
set Airfoil=JX-bad

rem ... is xoptfoil_visualizer-jx in the default directory?  
set LocalPath=..\..\windows\bin\
if not exist %LocalPath%xoptfoil_visualizer-jx.exe set LocalPath=

%LocalPath%xoptfoil_visualizer-jx -c %Airfoil% -o 3
