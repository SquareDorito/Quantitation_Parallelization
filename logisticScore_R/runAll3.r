source('E:/Temp/RLogisticScore/getDMods.r')
source('E:/Temp/RLogisticScore/getSMods.r')
source('E:/Temp/RLogisticScore/pseq.r')
source('E:/Temp/RLogisticScore/outmodsFile.r')
source('E:/Temp/RLogisticScore/getmass.r')
source('E:/Temp/RLogisticScore/getavgmass.r')
source('E:/Temp/RLogisticScore/vectorFind.r')
source('E:/Temp/RLogisticScore/getNumAfter.r')
source('E:/Temp/RLogisticScore/isolate.r')
source('E:/Temp/RLogisticScore/pepcalc.r')
source('E:/Temp/RLogisticScore/findTop.r')
source('E:/Temp/RLogisticScore/uscore.r')
source('E:/Temp/RLogisticScore/filterFile.r')
source('E:/Temp/RLogisticScore/validate.r')
source('E:/Temp/RLogisticScore/makeValidArray.r')
source('E:/Temp/RLogisticScore/parse_dt_out2.r')
source('E:/Temp/RLogisticScore/matComp.r')
source('E:/Temp/RLogisticScore/predictMain.r')
source('E:/Temp/RLogisticScore/PredictGui.r')
source('E:/Temp/RLogisticScore/logreg2.R')
source('E:/Temp/RLogisticScore/commandMain3.R')
source('E:/Temp/RLogisticScore/combineDtaOut.r')
dtaDir='E:/Temp/dta'
outDir='E:/Temp/out'
commandMain3('E:/Temp/RLogisticScore/Reduced.Rdata',dtaDir,outDir)    ##need to clear folder \\proteome\Data\LogisticTempFiles\DtaOut and ...\Validations first!
#q(save="no")
