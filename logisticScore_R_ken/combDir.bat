del E:\Temp\Validations /Q
del E:\Temp\RLogisticScore\DtaOut /Q
del E:\Temp\RLogisticScore\DtaOut1 /Q
del E:\Temp\RLogisticScore\DtaOut2 /Q
del E:\Temp\RLogisticScore\DtaOut3 /Q

cd c:\Program Files (x86)\R\R-2.4.1\bin
echo source("E:/Temp/RLogisticScore/runComb.r") | Rterm.exe --vanilla
exits