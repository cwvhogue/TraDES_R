#------------------------------
# General TRADES data analysis, standard Histogram plotting, 
# Energy calibration, Boltzmann ensemble extraction and Ramachandran plotting 
#
# (c) CWV Hogue June 2012. National University of Singapore
# Licence for the following R code is granted under the Creative Commons licence.
# Fine print is here: http://creativecommons.org/licenses/by-nc-sa/2.5/legalcode.
#
#
#  Instructions/Comments come before each function
#  Use source("TRADES_R_Analysis_Package.r") to load these functions.
#   
#
# BEFORE YOU START:
# You must first run a simulation with the TRADES package and create a *.log file of results.
# Minimum Run: trades - 300,000 structures per run  
#
# Also works with concatenated runs on the same *.trj file e.g.  30000 for 10 runs
# Set up a *.bat or shell script to run trades with structures numbered consecutively:
# 1-30000, 30001-60000, 60001-70000, ..., 270001-300000
# Concatenate your log files in numerical order to make one large log file.
#
# In the Windows command line this can be done with:
# copy /b file1.log + file2.log + file3.log bigfile.log

# MAKE Sure your log file has exactly a multiple of 300000 structures!
# Why?  Because it is the minimum to approach saturation of conformational space of 14-22 amino acids
# AND because it is easily factored for the Boltzmann fitting.



#---Start Here---------------------------------------------------------
#
# For starting ensembles of about 300,000 structures:
# Use 
#>TRADES.AnalyzeLog()->Expt_Results  - this will return graphs in your current working directory - use getwd() to see where - and returns summary stats
#
# For SMALL TraDES samples 
# Use 
#>TRADES.AnalyzeLog(Boltzmann=FALSE)
#
# TRADES.AnalyzeLog
# Needs no parameters, will prompt you for location of log file
# -parses the log file (wait...) 
# -Look in your getwd() for the output files:
# Graphs and files output: 
# xxxxxx_Histograms.png
# xxxxxx_EnergyScale.png
# xxxxxx_EnergyPlot.png
# xxxxxx_nnn_SpBoltzmann_State_T25.png
# xxxxxx_BoltzmannScanEnergy_T25.png
# xxxxxx_BoltzmannScanHist_T25.png
# xxxxxx_BoltzmannScanData_T25.csv
# xxxxxx_BoltzmanScanIter_T25.png
# xxxxxx_Ensemble_{nnn{_SpBS_T25_{Oblate/Prolate}.txt
# xxxxxx is the base name of your TraDES run embedded in the logfile

# Make sure your Structs_per_Boltz_State value is a number that divides your ensemble into even piecies
# i.e. (Sample Size) / (Structus_per_Boltz_State) has no remainder!
# you may need to remove lines from the end of your logfile to make it a round number.

# if you have a new trades log, use as-is, no paraeters needed.
# if you have a foldtraj log + a solvateL log, set zfile= the filename of the solvateL log
# so that TRADES.readEnergyZfile can patch the two together to mimic a TraDES-2 logfile.


TRADES.AnalyzeLog<-function(file, Boltzmann = TRUE, zfile, ExperimentalRgyr=0, ExperimentalNCDist=0, T=25){
   library(MASS)
   # main parser for foldtraj or trades format log files
   logF<-TRADES.readlog(logfile=file)
   if (logF$logformat == "foldtraj") {
       cat("foldtraj logs require *.csv file made by solvateL to complete analysis:")
       if (!missing(zfile)) logF<-TRADES.readEnergyZfile(logF,zfile)
   }
  # correct log structures for energy signs, Zhang score conversions, delete unused objects, compact others into stat values...
   logF<-TRADES.AdjustLog(logF) 
   cat("\nComputing the maximal length..\n")
   cat("Computing the size of conformational space and saturation peptide length..\n")
   Length_space<-TRADES.LengthAndSpace(logF$protein,length(logF$values$Structure))
   if (logF$logformat == "TraDES") {
     cat("Computing scaled energy..\n")
     cat("Plotting energy terms..\n")
     filename<-paste(logF$BaseName,"_EnergyScale.png",sep="")
     logF<-TRADES.EnergyScale(logF,file=filename)
     filename<-paste(logF$BaseName,"_EnergyPlot_T25.png",sep="")
     TRADES.EnergyPlot(logF,file=filename) 
     # plot basic histograms 
     Stats<-TRADES.LogAnalyze(logF)
     # plot Ellipsoid data
     Ellipsoid<- TRADES.Ellipsoid(logF, plot=TRUE)
   }
   else {
      cat("Cannot scale energy from old style foldtraj log.\n Option 1: Use trades with your old *.trj file\n Option 2: Reprocess old foldtraj *.val data with solvateL and pass in .csv filename\n")
         Boltzmann=FALSE
   }
   if (Boltzmann == TRUE) {
#     fileout=paste(logF$BaseName,"_AnalyzeLogSavedB4Boltzmann.Rdata",sep="")
#     cat("Saving R object:",fileout,"..\n",sep=" ")
#     values<-list(Length_Space = Length_space, Stats = Stats, Ellipsoid = Ellipsoid, Data=logF )
#     save(values,file=fileout)
     Boltzmann_Ensembles<-NULL
     Boltzmann_Ensembles<-TRADES.BoltzmannSearch(logF, T=T, ExperimentalRgyr=ExperimentalRgyr, ExperimentalNCDist=ExperimentalNCDist)
     cat("Completed.\n")
     values<-list(Boltzmann_Ensembles= Boltzmann_Ensembles,  Length_Space = Length_space, Stats = Stats, Ellipsoid = Ellipsoid, Data=logF )
   }
   else {
    values<-list(Length_Space = Length_space, Stats = Stats, Ellipsoid = Ellipsoid, Data=logF )
   }
fileout=paste(logF$BaseName,"_AnalyzeLogSaved.Rdata",sep="")
cat("Saving R object:",fileout,"\n",sep=" ")
save(values,file=fileout)
values
}

#####
# for Cn3D concatenation - Cn3D limit is about 50 models of 400 aa
# so 20000 amino acids total in N models
# To get M the number of models, 20000 / sequence length
# Sample M structures from optimal ensemble filename list in R
# copy sampled structures to holding directory
# rename structures to sequential numbering scheme, saving scheme translation
# use bin2prt to convert to ascii
# use concatmodels to concatenate to single *.prt
# use prt2bin to convert to *.val
# look at in Cn3D



#----------------------------
#
#  Note - most of the functions below are already called by TRADES.AnalyzeLog
#  Use them with the log R structure you already parsed in - if you want to re-run them to change parameters
#
#----------------------------



#-----------------------------------------------------------------------------------------------------------------
# TRADES.LengthAndSpace is a tool to estimate Proline dependent parameters of protein length and conformational space
# Main outputs
# A - the maximal extended length of the protein from alpha carbon of N to alpha carbon of C terminus
# B - the size of conformational space anticipated for sampling at 0.5 A and 6.0 A RMSD coverage    
# C - the length of sequence that will have on average, its conformational space saturated at 0.5 A and 0.6 A RMSD
# D - base parameters for reference including the %Gly - if this is grossly above 7.5% (e.g. keratin tails) expect underestimation of conformaitonal space.
# Peptides < 5 - this will not work.
# IF your protein is small, or your region of interest is small (like a binding motif) you can adjust the 
# sample size to fit the size of conformational space value. 
# 
# Proline rich proteins have less conformational space, require less sampling, and are shorter in length
# 
#Inputs are a character string (no breaks) representing the sequence, and size of ensemble
#Outputs are: 
# Input (seq)
# Input sample ensemble size (size)
# maximum extended length accounting for Pro in Angstroms (Extended_length)
# Size of conformational space at 2.5A RMSD (Conf_Space)
# Size of conformational space at 6.0A RMSD (Conf_Space6A)
# Size of conformational space as a random coil with PPII bias (Conf_Space_Coil)
# Average peptide length for 2.5 A saturation of conformational space given size of ensemble (Sat_Len) in residues
# Average peptide length for 6.0 A saturation of conformational space given size of ensemble (Sat_Len6A) in residues
# Average peptide length for saturation of space with coil propensity typical of an intrinsic disordered protein (Sat_Len_Coil)
# Sequence Length (Seq_Len)
# # of Pro (nPro)
# Number of NON Pro PPII constrained residues (nPrePro)
# # of Gly (nGly)
# # of Free Gly not constrained by Pro (nFreeGly)
# Average Omega value associated with 2.5A RMSD calculation (can range from 1.6 for poly-pro to 7.0224 for poly-gly) 
# For a sequence without Pro/Gly average Omega will be 2.65 

#Omega values 2.65 for 2.5A RMSD and 1.35 for 6A RMSD from Feldman & Hogue New Hope For Brute Force paper (2002)
#Extended lengths of 3.31A/aa for Beta, 3.1A/aa for polyproline II determined by PDB search over 5 residue lengths. 
#Omega value for Gly estimated at 2.65*2.65 (7.0 ish) owing to Ramachandran symetry 
#PMID:11746699

#Proline Omega value 1.6 from Hamburger, Ferreon, Whitten and Hilser (2004) Biochemistry 43:9790 
#Other residue Omega cited in this reference for non-pro is 7.9 which is outdated, too high 
# - but close to the Gly estimate used above. 
#PMID:15274633

# Omega for Coil (1.95) is based on the average Omega for a long sequence 
# "AAPAAPAAPAAPAAP..." which is constrained to PPII for 2/3 of residues
# This represents the reduced conformational space explored by an intrinsic disordered protein
# with a PPII bias in conformations as (shown by several NMR experiments)




TRADES.LengthAndSpace <- function(seq, size) {

seq_len <- nchar(seq)
seq_use <- substr(seq,2,seq_len-1)
gregexpr("G",seq_use)[[1]] -> glylocs
gregexpr("P",seq_use)[[1]] -> prolocs
npro <- length(prolocs)
ngly <- length(glylocs)
if (prolocs[[1]][1] < 0) npro = 0
if (glylocs[[1]][1] < 0) ngly = 0
npreproGly <- 0
nprepro <- 0
if (npro > 0) {
  strsplit(seq_use,"P")[[1]]->pieces
  piecelen <- nchar(pieces)
#strsplit makes only one piece if P on end of seq_use, so check for that condition...
  if (length(piecelen) > 1)  {
    prepro<-substr(pieces, piecelen, piecelen)
    npreproGly <- sum(na.omit(charmatch(prepro,"G",nomatch=0)))-sum(na.omit(charmatch(prepro,"",nomatch=0)))
    nprepro <- length(pieces[nchar(pieces) > 0])-1
  } 
}
if (nprepro < 0) nprepro <- 0
#see whether ignored C-terminal residue is proline, adjust nprepro if the penultimate redisue is also not a P
if ( charmatch( substr(seq,seq_len,seq_len),"P",nomatch=0) ) {
if ( !charmatch( substr(seq_use,nchar(seq_use),nchar(seq_use)),"P",nomatch=0) )  { nprepro <- nprepro + 1 }
if (charmatch(substr(seq_use,nchar(seq_use),nchar(seq_use)),"G",nomatch=0) ) {npreproGly <- npreproGly + 1 }
}
freegly <- ngly - npreproGly
if (freegly < 0) freegly = 0
#this is the number of Gly not in front of Pro, unconstrained.

# length calculations corrected to alpha-N to alpha-C distance, which TraDES reports
remainder <- seq_len-npro-nprepro
extended <- (npro *3.1) + (nprepro * 3.1) + (remainder * 3.31)
extended <- extended - (extended/seq_len)

#Note that the extended conformation reported is a C-alpha to C-alpha distance
#Conformational space is not counted for N- and C- terminal residues - they do not have dihedral angle constraints
# need to match Pro/Gly on ends, it will affect the estimates - esp for short sequences.

sequ_len <- nchar(seq_use)
normal <- sequ_len - freegly - npro - nprepro
if (normal < 0) normal = 0

conf_space = 1.6^(npro + nprepro) * 2.65^(normal) * (2.65*2.65)^(freegly)
conf_omega_ave = (1.6 * ((npro + nprepro)/(sequ_len))   
               +  2.65 * ((normal)/(sequ_len))
               +  2.65 * 2.65 * ((freegly)/(sequ_len)) )
conf_space6A = 1.35^(sequ_len)
conf_space_PPII = 1.95^(sequ_len)
sat_len <- log(size)/log(conf_omega_ave)
sat_len6A <- log(size)/log(1.35)
sat_lenPPII <- log(size)/log(1.95)


# recompute numbers of residues for entire sequence
gregexpr("G",seq)[[1]] -> glylocs
gregexpr("P",seq)[[1]] -> prolocs
npro <- length(prolocs)
ngly <- length(glylocs)
if (prolocs[[1]][1] < 0) npro = 0
if (glylocs[[1]][1] < 0) ngly = 0

values <-list(Seq = seq, Extended_length_A = extended,  SampleSize = size,
Conf_Space = conf_space, Conf_Space_6A = conf_space6A, Conf_Space_Coil = conf_space_PPII,
Sat_Len = as.integer(round(sat_len)), Sat_Len6A = as.integer(round(sat_len6A)), 
Sat_Len_Coil = as.integer(round(sat_lenPPII)), Seq_Len = seq_len,
nPro = npro, NprePro = nprepro, nGly = ngly, nFreeGly = freegly, Omega_ave = conf_omega_ave
)
values
}
##################################################3
#TRADES.EnergyScale applies the appropriate scaling terms for Bryant, Crease, Zhang, VSCORE potentials


TRADES.EnergyScale<-function(logF, file){

if (mode(logF$Data) == "list")
   logF<-logF$Data

if (logF$logformat != "TraDES") {
 
 cat("TRADES.EnergyScale requies a TraDES-2 format logfile.")
 return(NULL)
}   

T=25

if (!missing(file)) {
            png(filename=file,width=1500, height=500)
            par(mfrow=c(1,3),cex=1.8)
}
else par(mfrow=c(1,3))


R<- 1.9858775
kt = (R * (T + 273.15))/1000
len<- logF$values$N[1]

#This is empirically determined by testing to match ensemble selected Rgyr
#or rNC with known distribution from experiment.
#The LenRatio corrects for shorter/longer sequences as a baseline backbone configurational entropy terms
LenRatio = logF$seqlength/306
Factor = 1.0 + (3.5 * LenRatio)
CreaseScale = Factor

ZhangEnergy <- logF$values$Zhang1/21
BryantEnergy <- as.vector (-kt * logF$values$Bryant3)
CreaseEntropy <- CreaseScale/kt * as.vector(logF$values$Crease3 / logF$seqlength)
BTotalEnergy <- BryantEnergy + as.vector(kt * CreaseEntropy)


BryantStats<-TRADES.histogram(BryantEnergy, plot=FALSE)
VoronoiStats<-TRADES.histogram(logF$values$VSCORE1, plot=FALSE)
# note that we don't use the peak or mean or std dev 
# in case the distributions are skewed or multimodal - we use the extreme 
# left and right half-peak maximum values, and their average as the center
# also we are not scaling to match Zhang... 
VoronoiCenter<-(VoronoiStats$data_Left_HM + VoronoiStats$data_Right_HM) / 2
BryantCenter<-(BryantStats$data_Left_HM + BryantStats$data_Right_HM) / 2
VZero<-logF$values$VSCORE1 - VoronoiCenter
cat("Voronoi Center:", VoronoiCenter, "Bryant Center:" , BryantCenter, " ")
Vscale = BryantStats$data_FWHM / VoronoiStats$data_FWHM 
cat("Vscale:",Vscale,"\n")
Vfinal = Vscale * VZero + BryantCenter


logF$values$EnergyZhang<-ZhangEnergy
logF$values$EnergyBryant<-BryantEnergy
logF$values$EntropyCrease<-CreaseEntropy
logF$values$EnergyFreeBC<-BTotalEnergy
logF$values$EnergyVSCORE <- Vfinal
logF$values$EnergyFreeVC<- Vfinal + as.vector(kt * CreaseEntropy)



d=density(na.omit(logF$values$EntropyCrease))
plot(d, main=" ", ylab="Density", xlab="Crease Entropy (kcal/mol)" )
v=density(na.omit(logF$values$EnergyVSCORE))
b=density(na.omit(logF$values$EnergyBryant))
plot(b, ylab="Density", main = "Bryant-blue  VSCORE-black", xlab= "(kcal/mol)",col="blue")
lines(v,col="black")
e=density(na.omit(logF$values$EnergyFreeVC))
f=density(na.omit(logF$values$EnergyFreeBC))
plot(f, ylab="Density", main ="^G BC-blue, ^G VC-black",xlab="kcal/mol)",col="blue")
lines(e, col="black")

par(col="black")
if (!missing(file)) dev.off()

logF
}




TRADES.EnergyPlot<-function(logF, T=25, file) {
if (mode(logF$Data) == "list")
   log<-logF$Data


if (logF$logformat != "TraDES") {
 
 cat("TRADES.EnergyPlot requies a TraDES-2 format logfile.")
 return(NULL)
}   
   
   
Lab.palette <-colorRampPalette(c("white","lightyellow", "lightcyan","cyan", "lightskyblue", "lightseagreen", "yellowgreen" ,"yellow", "goldenrod", "orange", "orange4", "firebrick", "darkred", "red", "darkmagenta", "magenta", "hotpink", "pink","lightpink","white"), space = "Lab")
if (!missing(file)) {
            png(filename=file,width=4000, height=3000)
            par(mfrow=c(3,4),pty="s",cex=3.6,lwd=4)
}
else par(mfrow=c(3,4))


#calculate ylimits


yVSCOREMin<-min(logF$values$EnergyVSCORE)-5
yVSCOREMax<-max(logF$values$EnergyVSCORE)+5
yZhangMin<-min(logF$values$EnergyZhang)-5
yZhangMax<-max(logF$values$EnergyZhang)+5
yBryantMin<-min(logF$values$EnergyBryant)-5
yBryantMax<-max(logF$values$EnergyBryant)+5
yCreaseMin<-min(logF$values$EntropyCrease)-5
yCreaseMax<-max(logF$values$EntropyCrease)+5
yFreeBCMax<-max(logF$values$EnergyFreeBC)+5
yFreeBCMin<-min(logF$values$EnergyFreeBC)-5
yFreeVCMax<-max(logF$values$EnergyFreeVC)+5
yFreeVCMin<-min(logF$values$EnergyFreeVC)-5

mins <-c(yVSCOREMin,yBryantMin,yCreaseMin,yFreeBCMin,yFreeVCMin)
maxs <-c(yVSCOREMin,yBryantMax,yFreeBCMax,yFreeVCMax)
yTotalMin <-min(mins)
yTotalMax <-max(maxs)
TotalLims <- c(yTotalMin,yTotalMax)



#Free Energy Plots
smoothScatter(logF$values$EnergyFreeBC ~ logF$values$NC, ylab="Free Energy BC (kcal/mol)",ylim=TotalLims, xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
abline(h=0)
smoothScatter(logF$values$EnergyFreeBC ~ logF$values$Rgyr, ylab="Free Energy BC (kcal/mol)",ylim=TotalLims, xlab="Rgyr (Angstroms)", colramp=Lab.palette)
abline(h=0)
smoothScatter(logF$values$EnergyFreeVC ~ logF$values$NC, ylab="Free Energy VC (kcal/mol)", ylim=TotalLims, xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
abline(h=0)
smoothScatter(logF$values$EnergyFreeVC ~ logF$values$Rgyr, ylab="Free Energy VC (kcal/mol)",ylim=TotalLims, xlab="Rgyr (Angstroms)", colramp=Lab.palette)
abline(h=0)



#Correlation Plots
smoothScatter( logF$values$EnergyZhang ~ logF$values$EnergyVSCORE, ylab="Zhang Energy (kcal/mol)", ylim=c(yZhangMin,yZhangMax) ,xlim=TotalLims, xlab="VSCORE Energy (kcal/mol)", colramp=Lab.palette)
abline(h=0)
abline(v=0)
smoothScatter( logF$values$EnergyBryant ~ logF$values$EnergyVSCORE, ylab="Bryant Energy (kcal/mol)",ylim=TotalLims,xlim=TotalLims, xlab="VSCORE Energy (kcal/mol)", colramp=Lab.palette)
abline(h=0)
abline(v=0)
smoothScatter( logF$values$EntropyCrease ~ logF$values$EnergyVSCORE, ylab="Crease Entropy (kcal/mol)", ylim=c(yTotalMin,yCreaseMax),xlim=TotalLims, xlab="VSCORE Energy (kcal/mol)", colramp=Lab.palette)
abline(h=0)
abline(v=0)
smoothScatter(logF$values$EnergyFreeVC ~ logF$values$EnergyFreeBC, ylab="Free Energy VC (kcal/mol)",ylim=TotalLims, xlim=TotalLims, xlab="Free Energy BC", colramp=Lab.palette)
abline(h=0)
abline(v=0)

#Rgyr Plots
smoothScatter(logF$values$EnergyZhang ~ logF$values$Rgyr, ylab="Zhang Energy (kcal/mol)",ylim=c(yZhangMin,yZhangMax), xlab="Rgyr (Angstroms)", colramp=Lab.palette)
abline(h=0)
smoothScatter(logF$values$EnergyVSCORE ~ logF$values$Rgyr, ylab="VSCORE Energy (kcal/mol)", ylim=TotalLims, xlab="Rgyr (Angstroms)", colramp=Lab.palette)
abline(h=0)
smoothScatter(logF$values$EnergyBryant ~ logF$values$Rgyr, ylab="Bryant Energy (kcal/mol)",ylim=TotalLims, xlab="Rgyr (Angstroms)", colramp=Lab.palette)
abline(h=0)
smoothScatter(logF$values$EntropyCrease ~ logF$values$Rgyr, ylab="Crease Entropy (kcal/mol)",ylim=c(yCreaseMin,yCreaseMax), xlab="Rgyr (Angstroms)", colramp=Lab.palette)
abline(h=0)


#Rnc Plots
#smoothScatter(logF$values$EnergyZhang ~ logF$values$NC, ylab="Zhang Energy (kcal/mol)",ylim=TotalLims, xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
#abline(h=0)
#smoothScatter(logF$values$EnergyVSCORE ~ logF$values$NC, ylab="VSCORE Energy (kcal/mol)", ylim=TotalLims, xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
#abline(h=0)
#smoothScatter( logF$values$EnergyBryant ~ logF$values$NC, ylab="Bryant Energy (kcal/mol)",ylim=TotalLims, xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
#abline(h=0)
#smoothScatter(logF$values$EntropyCrease ~ logF$values$NC, ylab="Crease Entropy (kcal/mol)",ylim=TotalLims, xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
#abline(h=0)


if (!missing(file)) dev.off()

}





#----------------------------------------------------------------------------------------
# 
# TRADES.LogAnalyze  - Generates histograms and statistics about sampled structures
# >TRADES.LogAnalyze(Expt_Log)->Expt_Results  - 
# this will make a graph "xxxxxxHistograms.png" in your current working directory and return summary stats
# The parameter "pdf=TRUE" will cause it to write "Histograms.pdf" instead of ".png"
# The Expt_Results will look like this, first stats for the 4 histograms, 
# FWHM is full width at half max for skewed peaks
# The number of the MinStruc and MaxStruc is the ordinal number of the structure with minimum and maximum values, so you can retrieve it.
# SolventStructs is the number of analyzed structures that were not NA marked when zwipe=T
# The last 3 numbers are the percentages of DSSP Helix, DSSP Beta sheet, and Elongated unpaired Beta strands (Ecaca) found in the ensemble 
#--------------Example Run of a log with 300000 structues ----------------- 
#> TRADES.LogAnalyze(CasSD_Beta)
#$HistStats
#$HistStats$Rgyr
#  data_Mean data_Std data_Median data_Peak data_FWHM data_Max data_MaxStruc
#1  69.90528 16.82931    67.91227  64.70227  38.71002 168.6062        194546
#  data_Min data_MinStruc
#1 27.44937        167040
#
#$HistStats$NCDist
#  data_Mean data_Std data_Median data_Peak data_FWHM data_Max data_MaxStruc
#1   162.674 67.70578    157.3337  147.1926  166.2128 511.7277        278694
#  data_Min data_MinStruc
#1 4.625165        206013
#
#$HistStats$Voronoi
#  data_Mean  data_Std data_Median data_Peak data_FWHM data_Max data_MaxStruc
#1  12.22159 0.5772207    12.28773   12.3873  1.222791 14.02713        184298
#  data_Min data_MinStruc
#1 7.835475         52872
#
#$HistStats$Zhang
#  data_Mean data_Std data_Median data_Peak data_FWHM      data_Max
#1 -5.669293 4.473612   -4.660024 -1.841243  7.896344 -4.761905e-05
#  data_MaxStruc  data_Min data_MinStruc
#1         45372 -49.16138        158388
#
#
#$SolventStructs
#[1] 249178
#
#$PctHelix
#[1] 0.001222222
#
#$PctEdssp
#[1] 0.07449129
#
#$PctEcaca
#[1] 44.11987


TRADES.LogAnalyze<-function(logF, pdf=FALSE){
if (mode(logF$Data) == "list")
   logF<-logF$Data
if (pdf == TRUE) TRADES.HistogramsPDF(logF)->HistStats
else TRADES.HistogramsPNG(logF)->HistStats
length(logF$values$Structure) - sum(is.na(logF$values$Rgyr))->NStructs
TotalResidues<-logF$seqlength * length(logF$values$Structure)
PctHelix<-100*sum(na.omit(logF$values$Helix))/TotalResidues
PctEdssp<-100*sum(na.omit(logF$values$Edssp))/TotalResidues
PctEcaca<-100*sum(na.omit(logF$values$Ecaca))/TotalResidues
values=list(HistStats=HistStats, Structs=NStructs, PctHelix = PctHelix, PctEdssp=PctEdssp, PctEcaca=PctEcaca)
values
}
 
##############################################################################################
# TRADES.AdjustLog
# This function corrects the energy terms to proper sign and kcal/mol.
# the scaling factor 21 from Zhang correction for variation in atomic coordination number - averaged over atom types
# it is already RT/(qr*2) where qr is the average coordination number for 18 atom types.
# Zhang potentials are Counts * log odds score - summed
# Zhang correction factor is /21 kcal/mol as per paper
# Bryant and Crease need an inversion of the sign.


TRADES.AdjustLog<-function(logF, clean=TRUE) {
if (mode(logF$Data) == "list")
   logF<-logF$Data
   
   
len<- logF$values$N[1]

#Reduce some common log file terms to summaries
if (clean==TRUE) {
logF$seqlength<-len

logF$Time_mean <- mean(logF$values$Time)
logF$Time_std <- sd(logF$values$Time)

logF$Tries_mean <- mean(logF$values$Tries)
logF$Tries_std <- sd(logF$values$Tries)

logF$Crashes_mean <- mean(logF$values$Crashes)
logF$Crashes_std <- sd(logF$values$Crashes)

logF$BadBB_mean <- mean(logF$values$BadBB)
logF$BadBB_std <- sd(logF$values$BadBB)

#Remove some of the redundant terms and TRADES debugging terms in the table
logF$values$Time<-NULL
logF$values$Tries<-NULL
logF$values$ViolatedConstr<-NULL
logF$values$BadBB<-NULL
logF$values$Crashes<-NULL
logF$values$N<-NULL
}

#Make some real file names for structures

StrucNames<-NULL
len_data<-length(logF$values$Structure)
no_names<-length(logF$BaseNamesUnique)
pieces<-sum(logF$BaseNameCounts)
logsize<-0
#30000
linect<- 1
for (i in 1:no_names) {
logsize <- logsize + ((len_data/pieces) * logF$BaseNameCounts[i])
StrucNames<-c(StrucNames,paste(logF$BaseNamesUnique[i],"_",sprintf("%07d",logF$values$Structure[linect:logsize]),".val",sep=""))
linect<-linect+((len_data/pieces) * logF$BaseNameCounts[i])
}

logF$values$Filenames<-StrucNames
cat("Structure Filenames should look like this:\n")
cat(head(logF$values$Filenames,5),sep="\n")
cat(tail(logF$values$Filenames,5), sep="\n")

#Adjust Energy Terms According to Old/New foldtraj/TraDES log format

   
if (length(logF$logformat) == 0) {
  logF$values$Zhang <- logF$values$Zhang/21
	logF$values$Crease <- -logF$values$Crease
# Old behavior on object parsed with old code
# Return the logF handed in
	logF
}
else {
	if (logF$logformat == "foldtraj") {
#foldtraj log object, parsed with new code
		logF$values$Zhang <- logF$values$Zhang/21
		logF$values$Crease <- -logF$values$Crease
	} 
logF
}
}


 
 


TRADES.HistogramsPDF<-function(logF, filename="Histograms.pdf"){
if (mode(logF$Data) == "list")
   logF<-logF$Data
paste(logF$BaseName,"_",filename,sep="")->filename
pdf(filename, colormodel="cmyk" )
par(mfrow=c(2,2))
TRADES.histogram(logF$values$Rgyr, name="Radius of Gyration",  xlabel="Angstroms", peak=T, rug=T)->Rgyr
TRADES.histogram(logF$values$NCdist, name="N-C Distance", xlabel="Angstroms", peak=T, rug=T)->NCDist
TRADES.histogram(logF$values$EnergyFreeBC, name="Free Energy BC", xlabel="(kcal/mol)", peak=T, rug=T)->FreeEnergyVC
TRADES.histogram(logF$values$EnergyFreeVC, name="Free Energy VC",xlabel="(kcal/mol)", peak=T, rug=T)->FreeEnergyBC
dev.off()
values=list(Rgyr = Rgyr, NCDist=NCDist, FreeEnergyVC=FreeEnergyVC, FreeEnergyBC=FreeEnergyBC)
values
}


TRADES.HistogramsPNG<-function(logF, filename="Histograms.png"){
if (mode(logF$Data) == "list")
   log<-logF$Data
paste(logF$BaseName,"_",filename,sep="")->filename
png(filename,width=3000, height=3000)
par(mfrow=c(2,2), cex=3.6,lwd=4)
TRADES.histogram(logF$values$Rgyr, name="Radius of Gyration",  xlabel="Angstroms", peak=T, rug=T)->Rgyr
TRADES.histogram(logF$values$NCdist, name="N-C Distance", xlabel="Angstroms", peak=T, rug=T)->NCDist
TRADES.histogram(logF$values$EnergyFreeBC, name="Free Energy BC", xlabel="(kcal/mol)", peak=T, rug=T)->FreeEnergyVC
TRADES.histogram(logF$values$EnergyFreeVC, name="Free Energy VC",xlabel="(kcal/mol)", peak=T, rug=T)->FreeEnergyBC
dev.off()
values=list(Rgyr = Rgyr, NCDist=NCDist, FreeEnergyVC=FreeEnergyVC, FreeEnergyBC=FreeEnergyBC)
values
}

TRADES.histogram<-function(data, name="Histogram", xlabel=" ", peak=FALSE, xlimits=xlimits, rug=FALSE, plot=TRUE, breaks=100){

analysis<-data.frame(
  data_N=numeric(1),
  data_Mean=numeric(1),
  data_Std=numeric(1),
  data_Median=numeric(1),
  data_Peak=numeric(1),
  data_FWHM=numeric(1),
  data_Left_HM=numeric(1),
  data_Right_HM=numeric(1),
  data_Max=numeric(1),
  data_MaxStruc=numeric(1),
  data_Min=numeric(1),
  data_MinStruc=numeric(1))
analysis$data_N=length(na.omit(data))
analysis$data_Mean<-mean(na.omit(data))
analysis$data_Std<-sd(na.omit(data))
analysis$data_Median<-median(na.omit(data))
d=density(na.omit(data))
analysis$data_Peak<-d$x[which.max(d$y)]
hm<-d$y[which.max(d$y)] / 2
d.y <- d$y
d.x <- d$x
above <- d.x[which(d.y>hm)]
analysis$data_FWHM = max(above) - min(above)
analysis$data_Left_HM = min(above)
analysis$data_Right_HM = max(above)

analysis$data_Max=max(na.omit(data))
as.logical(data==analysis$data_Max)-> Max_vals
analysis$data_MaxStruc=which(Max_vals)[1]
analysis$data_Min=min(na.omit(data))
as.logical(data==analysis$data_Min)-> Min_vals
analysis$data_MinStruc=which(Min_vals)[1]

#paste(name,".pdf",sep="")->filename
#pdf(filename,colormodel="cmyk")
#paste("Histogram of", name) ->Title
if(plot == TRUE) {
if (!missing(xlimits))
hist(na.omit(data), breaks=breaks, main=name, prob=TRUE, xlab=xlabel, xlim=xlimits, col="gray")
else
hist(na.omit(data), breaks=breaks, main=name, prob=TRUE, xlab=xlabel, col="gray")

par(col="red")
abline(v=analysis$data_Mean, lwd=2)
abline(v=analysis$data_Mean+analysis$data_Std,lty="dotted")
abline(v=analysis$data_Mean-analysis$data_Std,lty="dotted")
abline(v=analysis$data_Median,lty="dashed")
if (peak) {
par(col="blue")
abline(v=analysis$data_Peak, lwd=2)
abline(v=max(above), lty="dotted")
abline(v=min(above), lty="dotted")
lines(d, lwd=2)
}
if (rug)
{
rug(data)
}
par(col="black")
#dev.off()
}
analysis
}





#------------------------------------
# Ramachandran Distribution Plotting
# First you must create a phi/psi dataset from your TRADES output - requires you save in *.val format (foldtraj parameter -a T)

# C:\Users\chogue\Dropbox\TraDES>ramangles -

# RamAngles:
# Report Phi, Psi from TraDES *.val file
#  Creates up to 20 *.csv files for R plotting.
#   arguments:
#
#  -f  Input VAL File Name (NO EXTENSION). [File In]
#  -s  Foldtraj Range Start Number (optinal) [Integer]  Optional
#    default = 0
#    range from 1 to 9999999
#  -r  Foldtraj Range (optional) [Integer]  Optional
#    default = 0
#    range from 1 to 50000

# ramangles.exe reads a range of *.val files (up to 50,000) and creates 20 *.CSV file, one for each amino acid
# Suggestion is that more than 30,000 sample *.val files will hurt performance of R, and not change the graph.
# Each ramangle.exe generated .CSV file has 4 elements per row: structure #, residue #, phi, and psi
# Example output:
# [.]              [..]             A_Cas_3ST_.csv   *///no Cys!///
# D_Cas_3ST_.csv   E_Cas_3ST_.csv   F_Cas_3ST_.csv   G_Cas_3ST_.csv
# H_Cas_3ST_.csv   I_Cas_3ST_.csv   K_Cas_3ST_.csv   L_Cas_3ST_.csv
# M_Cas_3ST_.csv   N_Cas_3ST_.csv   P_Cas_3ST_.csv   Q_Cas_3ST_.csv
# R_Cas_3ST_.csv   S_Cas_3ST_.csv   T_Cas_3ST_.csv   V_Cas_3ST_.csv
# W_Cas_3ST_.csv   Y_Cas_3ST_.csv

# If you are missing amino acids in your sequence use dummy *.csv files that looks like this, made in Excel:
# ----file C_BLANK_.csv is used to fill in for a missing Cys in the above example
# 0,0,0,0
# 0,0,0,0
# 0,0,0,0
# 0,0,0,0

# Now with all 20 CSV files you can plot one or all 20 Ramachandran Plots

 
# A single Ramachandran plot can be made with TRADES.Ramaplot()
# it prompts you to choose the file with the GUI.

# The function for making arrays of 20 Ramachandran plots is:
# TRADES.Plot20RamaPNG()  or TRADES.Plot20RamaPDF()
# These are calibrated to make reproduction quality PNG (rgb) or PDF (cmyk) files.  
# It prompts for a multiple selection of 20 files representing all 20 amino acids.  
# Each file begins with the uppercase letter of the amino acid represented. You need to multiple selecta all 20 of these
# with (Shift- or control-click) from the GUI as input.
#
#
#----------------------------------------
#Ramachandran Summary Table
#
# This function uses the same picking method as the Ramachandran Graph Plotter
# It creates a table of values for assignment of % alpha, % beta, % ppII, etc to each residue type
# To use - assign to a variable
# coil_space<-TRADES.RamaSpaceSummary()
# write.csv(coil_space$Table, file="Coil_Space.csv)
# 
# - load into Excel and enjoy.
#
#
# Picking Individual Residues - Serine 29 = FQSPP   203 = PPSVS  190 = LSSSH
# If you want to look at Ramachandran space of a specific residue 
# like these to check Pro conformational influence on left- and right- side of Ser 
# Go try this... Load in one single *.csv file from ramangles output (described above)
# with the residue type. Then you can select the array of phi-psi angles with just that
# residue number - here are the three examples:
#
# serine<-read.csv(file.choose(),header=F)
# which(serine$V2 == 29)->ser_29_idx
# serine[ser_29_idx,]->ser_29
# which(serine$V2 == 203)->ser_203_idx
# serine[ser_203_idx,]->ser_203
# which(serine$V2 == 190)->ser_190_idx
# serine[ser_190_idx,]->ser_190
# str(ser_29)
# 'data.frame':   30000 obs. of  4 variables:
# $ V1: int  1 2 3 4 5 6 7 8 9 10 ...
# $ V2: int  29 29 29 29 29 29 29 29 29 29 ...
# $ V3: num  -68.9 -105.7 -62.5 -79.6 -130.1 ...
# $ V4: num  -26.5 127.4 -29.2 -19.4 -23 ...
# NOW this can go into either of the two basic functions:
# TRADES.RamaSum("S 29", ser_29) 
# TRADES.RamaPlot("S 29", ser_29)
# Can concatenate the columns together like this to 
# to make a table to see how TraDES sampling chooses pre-proline residues...
#
# TRADES.RamaSum("S 29", ser_29)$Values->Single_Serine
# Single_Serine
#                   S 29
# Structures 30000.000000
# Positions      1.000000
# % Alpha-R     16.606667
# % Beta        37.500000
# % PPII        43.020000
# % Alpha-L      1.680000
# % epsilon      1.193333
#
# Single_Serine <- cbind(Single_Serine,TRADES.RamaSum("S 190", ser_190)$Values)
# Single_Serine
#                   S 29        S 190
# Structures 30000.000000 30000.000000
# Positions      1.000000     1.000000
# % Alpha-R     16.606667    29.593333
# % Beta        37.500000    31.590000
# % PPII        43.020000    35.816667
# % Alpha-L      1.680000     1.893333
# % epsilon      1.193333     1.106667
# Single_Serine <- cbind(Single_Serine,TRADES.RamaSum("S 203", ser_203)$Values)
# round(Single_Serine, digits=5)
#                  S 29       S 190       S 203
# Structures 30000.00000 30000.00000 30000.00000
# Positions      1.00000     1.00000     1.00000
# % Alpha-R     16.60667    29.59333    29.07000
# % Beta        37.50000    31.59000    32.38333
# % PPII        43.02000    35.81667    36.03000
# % Alpha-L      1.68000     1.89333     1.54000
# % epsilon      1.19333     1.10667     0.97667

# then write it out to a csv file
# write.csv(Single_Serine, file="Proline_Effect_on_Serine.csv")
# and format it in Excel for publication.  
#


TRADES.Plot20RamaPNG<-function(filename="RamaPlot.png"){
png(filename,width=6000, height=6000)
choose.files(caption="Select 20 CSV files from TRADES Ramangle - In alphabetical order!")->filelist
par(mfrow=c(4,5),pty="s", cex=3.6,lwd=4)
TRADES.RamaPlot("A", aa=read.csv(filelist[1],header=F))
TRADES.RamaPlot("C", aa=read.csv(filelist[2],header=F))
TRADES.RamaPlot("D", aa=read.csv(filelist[3],header=F))
TRADES.RamaPlot("E", aa=read.csv(filelist[4],header=F))
TRADES.RamaPlot("F", aa=read.csv(filelist[5],header=F))
TRADES.RamaPlot("G", aa=read.csv(filelist[6],header=F))
TRADES.RamaPlot("H", aa=read.csv(filelist[7],header=F))
TRADES.RamaPlot("I", aa=read.csv(filelist[8],header=F))
TRADES.RamaPlot("K", aa=read.csv(filelist[9],header=F))
TRADES.RamaPlot("L", aa=read.csv(filelist[10],header=F))
TRADES.RamaPlot("M", aa=read.csv(filelist[11],header=F))
TRADES.RamaPlot("N", aa=read.csv(filelist[12],header=F))
TRADES.RamaPlot("P", aa=read.csv(filelist[13],header=F))
TRADES.RamaPlot("Q", aa=read.csv(filelist[14],header=F))
TRADES.RamaPlot("R", aa=read.csv(filelist[15],header=F))
TRADES.RamaPlot("S", aa=read.csv(filelist[16],header=F))
TRADES.RamaPlot("T", aa=read.csv(filelist[17],header=F))
TRADES.RamaPlot("V", aa=read.csv(filelist[18],header=F))
TRADES.RamaPlot("W", aa=read.csv(filelist[19],header=F))
TRADES.RamaPlot("Y", aa=read.csv(filelist[20],header=F))
dev.off()
}


TRADES.Plot20RamaPDF<-function(filename="RamaPlot.pdf"){
pdf(filename, colormodel="cmyk", useDingbats=TRUE, width=16, height=16 )
choose.files(caption="Select 20 CSV files from TRADES Ramangle - In alphabetical order!")->filelist
par(mfrow=c(4,5),pty="s")
TRADES.RamaPlot("A", aa=read.csv(filelist[1],header=F))
TRADES.RamaPlot("C", aa=read.csv(filelist[2],header=F))
TRADES.RamaPlot("D", aa=read.csv(filelist[3],header=F))
TRADES.RamaPlot("E", aa=read.csv(filelist[4],header=F))
TRADES.RamaPlot("F", aa=read.csv(filelist[5],header=F))
TRADES.RamaPlot("G", aa=read.csv(filelist[6],header=F))
TRADES.RamaPlot("H", aa=read.csv(filelist[7],header=F))
TRADES.RamaPlot("I", aa=read.csv(filelist[8],header=F))
TRADES.RamaPlot("K", aa=read.csv(filelist[9],header=F))
TRADES.RamaPlot("L", aa=read.csv(filelist[10],header=F))
TRADES.RamaPlot("M", aa=read.csv(filelist[11],header=F))
TRADES.RamaPlot("N", aa=read.csv(filelist[12],header=F))
TRADES.RamaPlot("P", aa=read.csv(filelist[13],header=F))
TRADES.RamaPlot("Q", aa=read.csv(filelist[14],header=F))
TRADES.RamaPlot("R", aa=read.csv(filelist[15],header=F))
TRADES.RamaPlot("S", aa=read.csv(filelist[16],header=F))
TRADES.RamaPlot("T", aa=read.csv(filelist[17],header=F))
TRADES.RamaPlot("V", aa=read.csv(filelist[18],header=F))
TRADES.RamaPlot("W", aa=read.csv(filelist[19],header=F))
TRADES.RamaPlot("Y", aa=read.csv(filelist[20],header=F))
dev.off()
}


TRADES.RamaPlot <- function(Title, aa) {
if (missing(aa)) aa <- read.csv(file.choose(), header=F)
if (missing(Title)) Title<- "Ramachandran Plot"

Lab.palette <-
       colorRampPalette(c("white","lightyellow", "lightcyan","cyan", "lightskyblue", "lightseagreen", "yellowgreen" ,"yellow", "goldenrod", "orange", "orange4", "firebrick", "darkred", "red", "darkmagenta", "magenta", "hotpink", "pink","lightpink","white"), space = "Lab")
par(las=1,pty="s")
smoothScatter(aa$V3,aa$V4, xlab=expression(phi), ylab=expression(psi), asp=1, ylim=c(-180,180), xlim=c(-180,180), main=Title, xaxs="i", yaxs="i", axes=F, frame.plot=F, colramp=Lab.palette)
prange<- seq(-180, 180, by=40)
axis(1, at = prange)
axis(2, at = prange)
abline(v=0)
abline(h=0)
}


TRADES.RamaSpaceSummary<-function(){
choose.files(caption="Select 20 CSV files from TRADES Ramangle - In alphabetical order!")->filelist


data_A<-TRADES.RamaSum("A", aa=read.csv(filelist[1],header=F))
data_C<-TRADES.RamaSum("C", aa=read.csv(filelist[2],header=F))
data_D<-TRADES.RamaSum("D", aa=read.csv(filelist[3],header=F))
data_E<-TRADES.RamaSum("E", aa=read.csv(filelist[4],header=F))
data_F<-TRADES.RamaSum("F", aa=read.csv(filelist[5],header=F))
data_G<-TRADES.RamaSum("G", aa=read.csv(filelist[6],header=F))
data_H<-TRADES.RamaSum("H", aa=read.csv(filelist[7],header=F))
data_I<-TRADES.RamaSum("I", aa=read.csv(filelist[8],header=F))
data_K<-TRADES.RamaSum("K", aa=read.csv(filelist[9],header=F))
data_L<-TRADES.RamaSum("L", aa=read.csv(filelist[10],header=F))
data_M<-TRADES.RamaSum("M", aa=read.csv(filelist[11],header=F))
data_N<-TRADES.RamaSum("N", aa=read.csv(filelist[12],header=F))
data_P<-TRADES.RamaSum("P", aa=read.csv(filelist[13],header=F))
data_Q<-TRADES.RamaSum("Q", aa=read.csv(filelist[14],header=F))
data_R<-TRADES.RamaSum("R", aa=read.csv(filelist[15],header=F))
data_S<-TRADES.RamaSum("S", aa=read.csv(filelist[16],header=F))
data_T<-TRADES.RamaSum("T", aa=read.csv(filelist[17],header=F))
data_V<-TRADES.RamaSum("V", aa=read.csv(filelist[18],header=F))
data_W<-TRADES.RamaSum("W", aa=read.csv(filelist[19],header=F))
data_Y<-TRADES.RamaSum("Y", aa=read.csv(filelist[20],header=F))


rnames=c("Structures", "Positions", "% Alpha-R", "% Beta", "% PPII", "% Alpha-L", "% epsilon")
cnames=c("% Total", "Average", "Max", "Min")
values=matrix(nrow=7, ncol=4, dimnames=list(rnames,cnames))
count=0
if (! is.na(data_A$Values[1])) {values<-cbind(values, data_A$Values[,1])
cnames<-c(cnames,"A") 
count <- count +1 }
if (! is.na(data_C$Values[1])) {values<-cbind(values, data_C$Values[,1])
cnames<-c(cnames,"C") 
count <- count +1}
if (! is.na(data_D$Values[1])) {values<-cbind(values, data_D$Values[,1])
cnames<-c(cnames,"D") 
count <- count +1}
if (! is.na(data_E$Values[1])) {values<-cbind(values, data_E$Values[,1])
cnames<-c(cnames,"E") 
count <- count +1}
if (! is.na(data_F$Values[1])) {values<-cbind(values, data_F$Values[,1]) 
cnames<-c(cnames,"F") 
count <- count +1}
if (! is.na(data_G$Values[1])) {values<-cbind(values, data_G$Values[,1])
cnames<-c(cnames,"G") 
count <- count +1}
if (! is.na(data_H$Values[1])) {values<-cbind(values, data_H$Values[,1])
cnames<-c(cnames,"H") 
count <- count +1}
if (! is.na(data_I$Values[1])) {values<-cbind(values, data_I$Values[,1])
cnames<-c(cnames,"I") 
count <- count +1}
if (! is.na(data_K$Values[1])) {values<-cbind(values, data_K$Values[,1])
cnames<-c(cnames,"K") 
count <- count +1}
if (! is.na(data_L$Values[1])) {values<-cbind(values, data_L$Values[,1])
cnames<-c(cnames,"L") 
count <- count +1}
if (! is.na(data_M$Values[1])) {values<-cbind(values, data_M$Values[,1])
cnames<-c(cnames,"M") 
count <- count +1}
if (! is.na(data_N$Values[1])) {values<-cbind(values, data_N$Values[,1])
cnames<-c(cnames,"N") 
count <- count +1}
if (! is.na(data_P$Values[1])) {values<-cbind(values, data_P$Values[,1])
cnames<-c(cnames,"P") 
count <- count +1}
if (! is.na(data_Q$Values[1])) {values<-cbind(values, data_Q$Values[,1])
cnames<-c(cnames,"Q") 
count <- count +1}
if (! is.na(data_R$Values[1])) {values<-cbind(values, data_R$Values[,1])
cnames<-c(cnames,"R") 
count <- count +1}
if (! is.na(data_S$Values[1])) {values<-cbind(values, data_S$Values[,1])
cnames<-c(cnames,"S") 
count <- count +1}
if (! is.na(data_T$Values[1])) {values<-cbind(values, data_T$Values[,1])
cnames<-c(cnames,"T") 
count <- count +1}
if (! is.na(data_V$Values[1])) {values<-cbind(values, data_V$Values[,1])
cnames<-c(cnames,"V") 
count <- count +1}
if (! is.na(data_W$Values[1])) {values<-cbind(values, data_W$Values[,1])
cnames<-c(cnames,"W") 
count <- count +1}
if (! is.na(data_Y$Values[1])) {values<-cbind(values, data_Y$Values[,1])
cnames<-c(cnames,"Y") 
count <- count +1}
dimnames(values)<-list(rnames,cnames)
idx=count+4


#total fraction of conform (%/100) occupied by num of residue i *   / total residue count in population 
values[1,1]<-values[1,5]
values[2,1]<-sum(values[2,5:idx])
values[3,1]<-100*(sum( (values[3,5:idx]/100) * (values[2,5:idx]*values[1,1]) ) /(values[1,1]*values[2,1]))
values[4,1]<-100*(sum( (values[4,5:idx]/100) * (values[2,5:idx]*values[1,1]) ) /(values[1,1]*values[2,1]))
values[5,1]<-100*(sum( (values[5,5:idx]/100) * (values[2,5:idx]*values[1,1]) ) /(values[1,1]*values[2,1]))
values[6,1]<-100*(sum( (values[6,5:idx]/100) * (values[2,5:idx]*values[1,1]) ) /(values[1,1]*values[2,1]))
values[7,1]<-100*(sum( (values[7,5:idx]/100) * (values[2,5:idx]*values[1,1]) ) /(values[1,1]*values[2,1]))

#average
values[1,2]<-values[1,5]
values[2,2]<-mean(values[2,5:idx])
values[3,2]<-mean(values[3,5:idx])
values[4,2]<-mean(values[4,5:idx])
values[5,2]<-mean(values[5,5:idx])
values[6,2]<-mean(values[6,5:idx])
values[7,2]<-mean(values[7,5:idx])

#max
values[1,3]<-values[1,5]
values[2,3]<-max(values[2,5:idx])
values[3,3]<-max(values[3,5:idx])
values[4,3]<-max(values[4,5:idx])
values[5,3]<-max(values[5,5:idx])
values[6,3]<-max(values[6,5:idx])
values[7,3]<-max(values[7,5:idx])

#min
values[1,4]<-values[1,5]
values[2,4]<-min(values[2,5:idx])
values[3,4]<-min(values[3,5:idx])
values[4,4]<-min(values[4,5:idx])
values[5,4]<-min(values[5,5:idx])
values[6,4]<-min(values[6,5:idx])
values[7,4]<-min(values[7,5:idx])


output=list(A=data_A,C=data_C,D=data_D,E=data_E,F=data_F,G=data_G,H=data_H,I=data_I,K=data_K,L=data_L,M=data_M,N=data_N,P=data_P,Q=data_Q,R=data_R,S=data_S,T=data_T,V=data_V,W=data_W,Y=data_Y, Table=values)
output
}


TRADES.RamaSum <- function(Title, aa) {
if (missing(aa)) aa <- { file.choose()-> filename
                         read.csv(filename, header=F) }
if (missing(Title)) Title<- "x"
phi<-aa$V3
psi<-aa$V4
beta<-length(which((psi > 50.00000000) & (phi <= -100.00000000)))
beta<- beta + length(which( (psi <= -100.00000000) &   (phi <= -100.00000000)))                    
ppII<-length(which( (psi > 50.00000000) & (phi > -100.00000000) &  (phi <= 0.00000000) ) )
ppII<- ppII + length(which( (psi <= -100.00000000) &  (phi > -100.00000000) & (phi <= 0.00000000) ) )                    
epsilon<-length(which( (psi > 100.00000000) &  (phi > 0.00000000) ) )
epsilon <- epsilon + length(which( (psi <= -50.00000000) & (phi > 0.00000000) ) )                    
alphaR<-length(which( (psi > -100.00000000) & (psi <= 50.00000000)  &   (phi <= 0.00000000) ) )  
alphaL<-length(which( (psi > -50.00000000)  & (psi <= 100.00000000) &   (phi > 0.00000000) ) )
Length<-length(psi)
Res<-unique(aa$V2)
ResNo<-length(Res)
Nstruc<-length(which(aa$V2 == Res[1]))
added=beta+ppII+epsilon+alphaR+alphaL
beta<-100*beta/added
ppII<-100*ppII/added
epsilon<-100*epsilon/added
alphaR<-100*alphaR/added
alphaL<-100*alphaL/added
added=beta+ppII+epsilon+alphaR+alphaL 
rnames=c("Structures", "Positions", "% Alpha-R", "% Beta", "% PPII", "% Alpha-L", "% epsilon")
cnames=(c(Title))
values=matrix(nrow=7, ncol=1, dimnames=list(rnames,cnames))
values[1,1]<-Nstruc
values[2,1]<-ResNo
values[3,1]<-alphaR
values[4,1]<-beta
values[5,1]<-ppII
values[6,1]<-alphaL
values[7,1]<-epsilon

output=list(Values=values,Positions=Res,Filename=filename )
if (length(phi) > 10) output
else list(Values=NA, Positions=NA, Filename=filename)
}

 

#-------------------------------------------------------------------------
# These are the functions that parse foldtraj logs and create the R objects

# TRADES.readEnergy reads a separate Z*.csv file (or concatenated ones) to 
# paste on additional energy terms to a log for testing
# The additional energy terms are created by ramangles with the -p switch
# The useful energy terms will be incorporated into foldtraj and TRADES.readlog
# after which this function will become deprecated.
# THIS IS an intermediate function to use if old foldtraj logs are combined with solvateL logs ... For testing purposes only.
# but it should work if you need to reanalyze old foldtraj data with solvateL and get at the Boltzmann refinement.
# the input logfile should be freshly parsed from TRADES.readlog()

TRADES.readEnergyZfile=function(logF, Zfile) {
if (missing(Zfile)) Zfile <- file.choose()

if (mode(logF$Data) == "list")
   logF<-logF$Data

zdata<-read.csv(Zfile,stringsAsFactors=FALSE)
dimnames(zdata)[[1]] -> FileNameList
dimnames(zdata)[[2]] -> headings

version_headings = c("Filename","Zhang1","Bryant4","Bryant3","Crease4","VoronoiF","VTotTerm","VSolvTerm","VSSTerm","VPDBatomN","VASATotal","SolvAtoms","SdX","SdY","SdZ","Ea","Eb","Ec")
if(all(headings == version_headings)) {


zdata[,1]->Zhang1
Zhang1 <- as.numeric(Zhang1[-(which(Zhang1 == headings[2]))])
#Zhang1<-Zhang1/21
#zdata[,2]->Bryant4
#Bryant4 <- -as.numeric(Bryant4[-(which(Bryant4 == headings[3]))])
zdata[,3]->Bryant3
Bryant3 <- -as.numeric(Bryant3[-(which(Bryant3 == headings[4]))])
#zdata[,4]->Crease4
#Crease4 <- - as.numeric(Crease4[-(which(Crease4 == headings[5]))])
#zdata[,5]->VoronoiF
#VoronoiF <- as.numeric(VoronoiF[-(which(VoronoiF == headings[6]))])
#zdata[,6]->VTotTerm
#VTotTerm <- as.numeric(VTotTerm[-(which(VTotTerm == headings[7]))])
#zdata[,7]->VSolvTerm
#VSolvTerm <- as.numeric(VSolvTerm[-(which(VSolvTerm == headings[8]))])
#zdata[,8]->VSSTerm
#VSSTerm <- as.numeric(VSSTerm[-(which(VSSTerm == headings[9]))])
#as.numeric(zdata[1,9])->VPDBatomN
#zdata[,10]->VASATotal
#VASATotal <- as.numeric(VASATotal[-(which(VASATotal == headings[11]))])
zdata[,11]->SolvAtoms
SolvAtoms <- as.numeric(SolvAtoms[-(which(SolvAtoms == headings[12]))])
#zdata[,12]->SdX
#SdX <- as.numeric(SdX[-(which(SdX == headings[13]))])
#zdata[,13]->SdY
#SdY <- as.numeric(SdY[-(which(SdY == headings[14]))])
#zdata[,14]->SdZ
#SdZ <- as.numeric(SdZ[-(which(SdZ == headings[15]))])
zdata[,15]->Ea
Ea <- as.numeric(Ea[-(which(Ea == headings[16]))])
zdata[,16]->Eb
Eb <- as.numeric(Eb[-(which(Eb == headings[17]))])
zdata[,16]->Ec
Ec <- as.numeric(Ec[-(which(Ec == headings[18]))])

FileNameList <- FileNameList[-which(FileNameList == headings[1])]

if(all(FileNameList == logF$values$Filenames)) {
# everything is in order, attach and finish
logF$values$Zhang1 <- Zhang1
#logF$values$Bryant4 <- Bryant4
logF$values$Bryant3 <- Bryant3
#logF$values$Crease4 <- Crease4
logF$values$Crease3 <- -logF$values$Crease
logF$values$VSCORE1 <- logF$values$Voronoi
logF$values$Voronoi <- NULL
logF$values$Zhang <- NULL
logF$values$Crease <- NULL
#logF$values$VoronoiF <- VoronoiF
#logF$values$VTotTerm <- VTotTerm
#$logF$values$VSolvTerm <- VSolvTerm
#logF$values$VSSTerm <- VSSTerm
#logF$values$VASATotal <- VASATotal
logF$values$SolvAtoms <- SolvAtoms
#logF$values$SdX <- SdX
#logF$values$SdY <- SdY
#logF$values$SdZ <- SdZ
logF$values$Ea <- Ea
logF$values$Eb <- Eb
logF$values$Ec <- Ec
#$logF$VPDBatomN <- VPDBatomN

# NOW equivalent to a TraDES-2 log
logF$logformat="TraDES-2"
logF
}
else
{
cat("vector of filenames do not match - looking to resolve\n")
cat("Length of CSV FileNameList:",length(FileNameList),"\n")
cat("Length of log Filenames:",length(logF$values$Filenames),"\n")

NULL
}

}
else
{
cat("Error: Input file",Zfile,"data not recognized\n")
cat("Specified file missing the expected header line:",version_headings,"\n")
cat("Instead these headings were read:",headings,"\n")
NULL
}
}



TRADES.readlog=function(logfile) {
if (missing(logfile)) logfile <- file.choose()

# File must start with a Foldtraj log on line 1 or 2
# DEAL WITH FIRST LINE SKIPPING (Tolerant of start at line 1 or 2)
# match the expected string, keep track of starting line number in skip variable.

startline<-0
FoldtrajVer <- scan(logfile, what=character(0), skip=startline, nlines=1)
if (length(FoldtrajVer != 0))  {cat("Loading Log Line 1:",FoldtrajVer,"\n",sep=" ") }
if (length(FoldtrajVer) == 0) {
    startline<-1
    FoldtrajVer <- scan(logfile, what=character(0), skip=startline, nlines=1)
    if (length(FoldtrajVer != 0))  {
	cat("Loading Log Line 2:",FoldtrajVer,"\n",sep=" ") 
   }
   else {
    return(1)
   }
   }
   
if ((FoldtrajVer[1] != "Foldtraj") && FoldtrajVer[1] != "TraDES"){ 
   cat("Foldtraj/TraDES-2 Header Not Found","\n",sep=" ")
   return(1)}
TraDES_LogFormat <- FoldtrajVer[1]   

# Read in the rest of the TOPMOST logfile information
Version=paste(FoldtrajVer,collapse=" ")
Line2 <- scan(logfile, what=character(0),sep='\t',skip=startline+1,nlines=1)
TrajFile<-paste(substring(Line2[1],18,100),".trj",sep="",collapse=NULL)

Protein <- scan(logfile, what=character(0), skip=startline+2, nlines=1)
Protein <- Protein[2]
Conditions <- scan(logfile, what=character(0), sep='\t',skip=startline+3, nlines=1)

# !!!! TO DO Need to break out the individual conditions too!  
# AND Above - they change as well depending on Unfoldtraj mode
# Cannot assume Conditions[12] is the method!!! - There are two diff tokens : 
# Trajectory Distribution: 3-State Secondary Structure Prediction
# and NA - based on Unfoldtraj - will say "Compared to Native Structure: 1YU5" on Line2, 
# but no "Trajectory Distribution.." entry on Conditions line
 


TrajMethod<-Conditions[12]
Conditions[-12]
Conditions<-c(Conditions,Line2)
Conditions[2]<-paste(substr(Conditions[2],9,100),": ",substr(Conditions[2],1,7),sep="",collapse=NULL)
Conditions[3]<-paste(substr(Conditions[3],9,100),": ",substr(Conditions[3],1,7),sep="",collapse=NULL)
System <- paste(scan(logfile, what=character(0), sep=NULL, skip=startline+4,nlines=1),sep="",collapse=" ")
StartDate <-scan(logfile, what=character(0),skip=startline+5,nlines=1)
StartDate <-paste(StartDate[3:7],sep="",collapse=" ")
Headers<-scan(logfile, what=character(0), sep='\t',skip=startline+7,nlines=1)

width<-length(Headers)
cat("Logfile Header Count:",width,"\n",sep=" ")

# Read in the tab delimited headers - indicates width of the table

# NOW READ IN THE REST OF THE FILE as Tab Delimited Strings
TempValues<-scan(logfile, what=character(0), sep='\t', skip = startline+8, quiet=TRUE)
cat("Values Scanned as Strings:",length(TempValues),"\n",sep=" ")


if (TraDES_LogFormat == "Foldtraj") {
# REMOVE INCOMPLETE ROWS
# incomplete rows are marked at the end "INCOMPLETE"
	incomplete <- which(TempValues == "INCOMPLETE")
# mark the rest of the row backwards by adding on 
# the entire set of list indexes minus 1 backwards to the table width
	for (i in 0:width) {
		if (i==0) expanded <-incomplete
		else expanded <- c(expanded, incomplete - i)
		i <- i+1 
	}
# now the object expanded has the list of all the values to remove for each INCOMPLETE row
	Complete<-TempValues[-expanded]
	cat("INCOMPLETE Rows removed:", length(incomplete),"Raw Strings:",length(Complete),"\n",sep=" ")
}
else {
	Complete<-TempValues
}

# next look for any end markers left in by Foldtraj - remove them in case concatenation left them in
incomplete<-which(Complete == "----------------------------------------")
Complete<-Complete[-incomplete]
cat("End marker Rows removed:", length(incomplete),"Raw Strings:",length(Complete),"\n",sep=" ")



# remove the semicolon end of line markers 
incomplete<-which(Complete == ";")
Complete<-Complete[-incomplete]
cat("Semicolons removed:", length(incomplete),"Raw Strings:",length(Complete),"\n",sep=" ")


# Find any other contatenated HEADERS buried in the file and add them to the top level objects

#foldtraj log
FTLine<- which(substr(Complete,1,8) == "Foldtraj")
if (length(FTLine) == 0) {
#TraDES-2 log
     FTLine<- which(substr(Complete,1,6) == "TraDES")

}


if (length(FTLine) !=0) {
         TLine<-which(substr(Complete,1,16) == "Trajectory File:")
         GLine<-which(substr(Complete,1,12) == "# Generated:")
         StLine<-which(substr(Complete,1,19) == "Start Numbering at:")
         CoLine<-which(substr(Complete,1,29) == "Compared to Native Structure:")
         BaLine<-which(substr(Complete,1,25) == "Structure File Base Name:")
         RsLine<-which(substr(Complete,1,12) == "Random Seed:")
         SqLine<-which(substr(Complete,1,9) == "Sequence:")

# !!!! Need to break out the individual conditions too!  AND Above - they change as well depending on Unfoldtraj mode

	 CLine<-which(substr(Complete,1,7) == "Folding")
	 CLine<-c(CLine, CLine+1, CLine+2, CLine+3, CLine+4,CLine+5,CLine+6,CLine+7,CLine+8,CLine+9,CLine+10)
	 CLine<-sort(CLine) 
	 MLine<-which(substr(Complete,1,24) == "Trajectory Distribution:")
	 SyLine<-which(substr(Complete,1,6) == "System")
	 JLine<-which(substr(Complete,1,3) == "Job") 
	 HLine<-which(Complete == "Structure")
	 Hremoved<-length(HLine)
	 for (i in 1:(width-1)) {
	   tHLine<-HLine + i
	   if (i == 1) aHLine<-HLine
	   aHLine<-c(aHLine,tHLine ) }
	 HLine<-aHLine

	# organize all the headers into one set of addresses, then extract from Complete
	LLines<-sort(c(FTLine, TLine, CoLine, StLine, GLine, BaLine, RsLine, SqLine, CLine, MLine, SyLine, JLine, HLine))

	HeadExtra<-Complete[LLines]
	
	SqLines<-Complete[SqLine]
	Seqs<-c(Protein,(substring(SqLines,11,6000)))
	SeqsUniq<-unique(Seqs)
	if(length(SeqsUniq) != 1) {
	cat("WTF? - Different protein sequences detected! Are you kidding me?\n")
	cat(SeqsUniq,sep="\n")
	cat("ERROR: Violation of space-time continuum - please correct your log file!\n")
	stop
	}
	TLines<-Complete[TLine] 
	TrajFiles_t<-c(substring(TLines,18,100))
	TrajFiles_t<-paste(TrajFiles_t,".trj",sep="",collapse=NULL)
	TrajFiles<-c(TrajFile,TrajFiles_t)
	TrajFilesUniq<-unique(TrajFiles)
	BaLines<-Complete[BaLine]
	BaseName<-substring(Conditions[14],27,100)
	BaseNames<-c(BaseName,(substring(BaLines,27,100)))
	BaseNamesUnique<-unique(BaseNames)
	BaseNameCounts<-NULL
	for (i in 1:length(BaseNamesUnique)) {
	 BaseNameCounts<-c(BaseNameCounts,length(BaseNames[BaseNames==BaseNamesUnique[i]]))
	 cat("Base Ensemble:",BaseNamesUnique[i], " from ", BaseNameCounts[i], "log file pieces\n") 
	}
	# now all the additional headers are removed
	Complete<-Complete[-LLines]
                    BaseName<-paste(BaseNamesUnique,sep="",collapse="")
	cat("Full Headers extracted:", Hremoved,"Raw Strings:",length(Complete),"\n",sep=" ")
} else { 
# Only one set of trj file listings in log, set all to that one.
	HeadExtra<- "0" 
	TrajFilesUniq<-TrajFile
	BaseName<-substring(Conditions[14],27,100)
	BaseNames <- BaseName
	BaseNamesUnique<-substring(Conditions[14],27,100)
	BaseNameCounts<-1
}

cat("Final Table Size:", length(Complete),"Columns:", width, '\n', sep=" ")

# Match columns with titles and push into data.frame by column.
Values<-TRADES.logframe(width,Complete)
# Values<-Complete

# Stitch the objects together and return the data frame with the parsed information and columns
# This object is ready for plotting
logF<-list(logformat=TraDES_LogFormat, startdate=StartDate, trajfile=TrajFilesUniq, BaseName = BaseName, BaseNamesUnique = BaseNamesUnique,
BaseNames = BaseNames, BaseNameCounts = BaseNameCounts,
protein=Protein, original=logfile, trajmethod=TrajMethod, versioN=Version, conditions=Conditions, 
cpu=System, logextra=HeadExtra,  headings=Headers, values=Values)
logF
}



TRADES.logframe = function( width, Complete) {
columnsz<-length(Complete)/width
if (width == 20) {
# foldtraj logs (width = 20,21,22,23)  
# NEW TraDES-2 logs (width = 29,30)
# no RMSD or CHARMM options, default  log format
# make a data frame with the shortened version of the Headings
	Values<-data.frame(
          Structure=integer(columnsz), 
          Time=integer(columnsz),
          Tries=integer(columnsz),
          BadBB=integer(columnsz),
          Crashes=integer(columnsz),
          ViolatedConstr=character(columnsz),
          N=integer(columnsz),
          Rgyr=numeric(columnsz),
          HRgyr=numeric(columnsz),
          NCdist=numeric(columnsz),
          Rn=numeric(columnsz),
          Cn=numeric(columnsz),
          ASA=numeric(columnsz),
          HASA=numeric(columnsz),
          Helix=integer(columnsz),
          Edssp=integer(columnsz),
          Ecaca=integer(columnsz),
          Zhang=numeric(columnsz),
          Voronoi=numeric(columnsz),
          Crease=numeric(columnsz))
          colnames(Values)<-c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca","Zhang","Voronoi","Crease")
# Now convert the array of strings into a column vector of the appropriate type
          Values$Structure<-as.integer(Complete[ seq(from = 1, to = length(Complete), by = width ) ]) 
          Values$Time<-as.integer(Complete[ seq(from = 2, to = length(Complete), by = width ) ])
          Values$Tries<-as.integer(Complete[ seq(from = 3, to = length(Complete), by = width ) ])
          Values$BadBB<-as.integer(Complete[ seq(from = 4, to = length(Complete), by = width ) ])
          Values$Crashes<-as.integer(Complete[ seq(from = 5, to = length(Complete), by = width ) ])
          Values$ViolatedConstr<-Complete[ seq(from = 6, to = length(Complete), by = width ) ]
          Values$N<-as.integer(Complete[ seq(from = 7, to = length(Complete), by = width ) ])
          Values$Rgyr<-as.numeric(Complete[ seq(from = 8, to = length(Complete), by = width ) ])
          Values$HRgyr<-as.numeric(Complete[ seq(from = 9, to = length(Complete), by = width ) ])
          Values$NCdist<-as.numeric(Complete[ seq(from = 10, to = length(Complete), by = width ) ])
          Values$Rn<-as.numeric(Complete[ seq(from = 11, to = length(Complete), by = width ) ])
          Values$Cn<-as.numeric(Complete[ seq(from = 12, to = length(Complete), by = width ) ])
          Values$ASA<-as.numeric(Complete[ seq(from = 13, to = length(Complete), by = width ) ])
          Values$HASA<-as.numeric(Complete[ seq(from = 14, to = length(Complete), by = width ) ])
          Values$Helix<-as.integer(Complete[ seq(from = 15, to = length(Complete), by = width ) ])
          Values$Edssp<-as.integer(Complete[ seq(from = 16, to = length(Complete), by = width ) ])
          Values$Ecaca<-as.integer(Complete[ seq(from = 17, to = length(Complete), by = width ) ])
          Values$Zhang<-as.numeric(Complete[ seq(from = 18, to = length(Complete), by = width ) ])
          Values$Voronoi<-as.numeric(Complete[ seq(from = 19, to = length(Complete), by = width ) ])
          Values$Crease<-(as.numeric(Complete[ seq(from = 20, to = length(Complete), by = width ) ]))
           } else if(width == 21) {
# includes RMSD column in the log
     Values<-data.frame(
          Structure=integer(columnsz), 
          Time=integer(columnsz),
          Tries=integer(columnsz),
          BadBB=integer(columnsz),
          Crashes=integer(columnsz),
          ViolatedConstr=integer(columnsz),
          N=integer(columnsz),
          Rgyr=numeric(columnsz),
          HRgyr=numeric(columnsz),
          NCdist=numeric(columnsz),
          Rn=numeric(columnsz),
          Cn=numeric(columnsz),
          ASA=numeric(columnsz),
          HASA=numeric(columnsz),
          Helix=integer(columnsz),
          Edssp=integer(columnsz),
          Ecaca=integer(columnsz),
          RMSD=numeric(columnsz),
          Zhang=numeric(columnsz),
          Voronoi=numeric(columnsz),
          Crease=numeric(columnsz))
          colnames(Values)<-c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca","RMSD", "Zhang","Voronoi","Crease")
 # Now convert the array of strings into a column vector of the appropriate type
          Values$Structure<-as.integer(Complete[ seq(from = 1, to = length(Complete), by = width ) ]) 
          Values$Time<-as.integer(Complete[ seq(from = 2, to = length(Complete), by = width ) ])
          Values$Tries<-as.integer(Complete[ seq(from = 3, to = length(Complete), by = width ) ])
          Values$BadBB<-as.integer(Complete[ seq(from = 4, to = length(Complete), by = width ) ])
          Values$Crashes<-as.integer(Complete[ seq(from = 5, to = length(Complete), by = width ) ])
          Values$ViolatedConstr<-Complete[ seq(from = 6, to = length(Complete), by = width ) ]
          Values$N<-as.integer(Complete[ seq(from = 7, to = length(Complete), by = width ) ])
          Values$Rgyr<-as.numeric(Complete[ seq(from = 8, to = length(Complete), by = width ) ])
          Values$HRgyr<-as.numeric(Complete[ seq(from = 9, to = length(Complete), by = width ) ])
          Values$NCdist<-as.numeric(Complete[ seq(from = 10, to = length(Complete), by = width ) ])
          Values$Rn<-as.numeric(Complete[ seq(from = 11, to = length(Complete), by = width ) ])
          Values$Cn<-as.numeric(Complete[ seq(from = 12, to = length(Complete), by = width ) ])
          Values$ASA<-as.numeric(Complete[ seq(from = 13, to = length(Complete), by = width ) ])
          Values$HASA<-as.numeric(Complete[ seq(from = 14, to = length(Complete), by = width ) ])
          Values$Helix<-as.integer(Complete[ seq(from = 15, to = length(Complete), by = width ) ])
          Values$Edssp<-as.integer(Complete[ seq(from = 16, to = length(Complete), by = width ) ])
          Values$Ecaca<-as.integer(Complete[ seq(from = 17, to = length(Complete), by = width ) ])
          Values$RMSD<-as.integer(Complete[ seq(from = 18, to= length(Complete), by = width ) ])
          Values$Zhang<-as.numeric(Complete[ seq(from = 19, to = length(Complete), by = width ) ])
          Values$Voronoi<-as.numeric(Complete[ seq(from = 20, to = length(Complete), by = width ) ])
          Values$Crease<-(as.numeric(Complete[ seq(from = 21, to = length(Complete), by = width ) ]))
           } else  if(width== 22) { 
# includes CHARMM energy (with extra unused column)
     Values<-data.frame(
          Structure=integer(columnsz), 
          Time=integer(columnsz),
          Tries=integer(columnsz),
          BadBB=integer(columnsz),
          Crashes=integer(columnsz),
          ViolatedConstr=integer(columnsz),
          N=integer(columnsz),
          Rgyr=numeric(columnsz),
          HRgyr=numeric(columnsz),
          NCdist=numeric(columnsz),
          Rn=numeric(columnsz),
          Cn=numeric(columnsz),
          ASA=numeric(columnsz),
          HASA=numeric(columnsz),
          Helix=integer(columnsz),
          Edssp=integer(columnsz),
          Ecaca=integer(columnsz),
          Zhang=numeric(columnsz),
          Voronoi=numeric(columnsz),
          Crease=numeric(columnsz),
          CHARMM=numeric(columnsz))
          colnames(Values)<-c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca","Zhang","Voronoi","Crease","CHARMM")
# Now convert the array of strings into a column vector of the appropriate type
          Values$Structure<-as.integer(Complete[ seq(from = 1, to = length(Complete), by = width ) ]) 
          Values$Time<-as.integer(Complete[ seq(from = 2, to = length(Complete), by = width ) ])
          Values$Tries<-as.integer(Complete[ seq(from = 3, to = length(Complete), by = width ) ])
          Values$BadBB<-as.integer(Complete[ seq(from = 4, to = length(Complete), by = width ) ])
          Values$Crashes<-as.integer(Complete[ seq(from = 5, to = length(Complete), by = width ) ])
          Values$ViolatedConstr<-Complete[ seq(from = 6, to = length(Complete), by = width ) ]
          Values$N<-as.integer(Complete[ seq(from = 7, to = length(Complete), by = width ) ])
          Values$Rgyr<-as.numeric(Complete[ seq(from = 8, to = length(Complete), by = width ) ])
          Values$HRgyr<-as.numeric(Complete[ seq(from = 9, to = length(Complete), by = width ) ])
          Values$NCdist<-as.numeric(Complete[ seq(from = 10, to = length(Complete), by = width ) ])
          Values$Rn<-as.numeric(Complete[ seq(from = 11, to = length(Complete), by = width ) ])
          Values$Cn<-as.numeric(Complete[ seq(from = 12, to = length(Complete), by = width ) ])
          Values$ASA<-as.numeric(Complete[ seq(from = 13, to = length(Complete), by = width ) ])
          Values$HASA<-as.numeric(Complete[ seq(from = 14, to = length(Complete), by = width ) ])
          Values$Helix<-as.integer(Complete[ seq(from = 15, to = length(Complete), by = width ) ])
          Values$Edssp<-as.integer(Complete[ seq(from = 16, to = length(Complete), by = width ) ])
          Values$Ecaca<-as.integer(Complete[ seq(from = 17, to = length(Complete), by = width ) ])
          Values$Zhang<-as.numeric(Complete[ seq(from = 18, to = length(Complete), by = width ) ])
          Values$Voronoi<-as.numeric(Complete[ seq(from = 19, to = length(Complete), by = width ) ])
          Values$Crease<-(as.numeric(Complete[ seq(from = 20, to = length(Complete), by = width ) ]))
          Values$CHARMM<-as.numeric(Complete[ seq(from = 22, to = length(Complete), by = width ) ])
           } else if(width == 23) {
# includes RMSD and CHARMM energy (with extra unused colunm)
     Values<-data.frame(
          Structure=integer(columnsz), 
          Time=integer(columnsz),
          Tries=integer(columnsz),
          BadBB=integer(columnsz),
          Crashes=integer(columnsz),
          ViolatedConstr=integer(columnsz),
          N=integer(columnsz),
          Rgyr=numeric(columnsz),
          HRgyr=numeric(columnsz),
          NCdist=numeric(columnsz),
          Rn=numeric(columnsz),
          Cn=numeric(columnsz),
          ASA=numeric(columnsz),
          HASA=numeric(columnsz),
          Helix=integer(columnsz),
          Edssp=integer(columnsz),
          Ecaca=integer(columnsz),
          RMSD=numeric(columnsz),
          Zhang=numeric(columnsz),
          Voronoi=numeric(columnsz),
          Crease=numeric(columnsz),
          CHARMM=numeric(columnsz))
          colnames(Values)<- c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca","RMSD","Zhang","Voronoi","Crease","CHARMM")
# Now convert the array of strings into a column vector of the appropriate type
          Values$Structure<-as.integer(Complete[ seq(from = 1, to = length(Complete), by = width ) ]) 
          Values$Time<-as.integer(Complete[ seq(from = 2, to = length(Complete), by = width ) ])
          Values$Tries<-as.integer(Complete[ seq(from = 3, to = length(Complete), by = width ) ])
          Values$BadBB<-as.integer(Complete[ seq(from = 4, to = length(Complete), by = width ) ])
          Values$Crashes<-as.integer(Complete[ seq(from = 5, to = length(Complete), by = width ) ])
          Values$ViolatedConstr<-Complete[ seq(from = 6, to = length(Complete), by = width ) ]
          Values$N<-as.integer(Complete[ seq(from = 7, to = length(Complete), by = width ) ])
          Values$Rgyr<-as.numeric(Complete[ seq(from = 8, to = length(Complete), by = width ) ])
          Values$HRgyr<-as.numeric(Complete[ seq(from = 9, to = length(Complete), by = width ) ])
          Values$NCdist<-as.numeric(Complete[ seq(from = 10, to = length(Complete), by = width ) ])
          Values$Rn<-as.numeric(Complete[ seq(from = 11, to = length(Complete), by = width ) ])
          Values$Cn<-as.numeric(Complete[ seq(from = 12, to = length(Complete), by = width ) ])
          Values$ASA<-as.numeric(Complete[ seq(from = 13, to = length(Complete), by = width ) ])
          Values$HASA<-as.numeric(Complete[ seq(from = 14, to = length(Complete), by = width ) ])
          Values$Helix<-as.integer(Complete[ seq(from = 15, to = length(Complete), by = width ) ])
          Values$Edssp<-as.integer(Complete[ seq(from = 16, to = length(Complete), by = width ) ])
          Values$Ecaca<-as.integer(Complete[ seq(from = 17, to = length(Complete), by = width ) ])
          Values$RMSD<-as.integer(Complete[ seq(from = 18, to= length(Complete), by = width ) ])
          Values$Zhang<-as.numeric(Complete[ seq(from = 19, to = length(Complete), by = width ) ])
          Values$Voronoi<-as.numeric(Complete[ seq(from = 20, to = length(Complete), by = width ) ])
          Values$Crease<-(as.numeric(Complete[ seq(from = 21, to = length(Complete), by = width ) ]))
          Values$CHARMM<-as.numeric(Complete[ seq(from = 23, to = length(Complete), by = width ) ])
           }
     else if(width == 29) {
# NEW TRADES-2 LOG FORMAT 
     Values<-data.frame(
          Structure=integer(columnsz), 
          Time=integer(columnsz),
          Tries=integer(columnsz),
          BadBB=integer(columnsz),
          Crashes=integer(columnsz),
          ViolatedConstr=integer(columnsz),
          N=integer(columnsz),
          Rgyr=numeric(columnsz),
          HRgyr=numeric(columnsz),
          NCdist=numeric(columnsz),
          Rn=numeric(columnsz),
          Cn=numeric(columnsz),
          ASA=numeric(columnsz),
          HASA=numeric(columnsz),
          Helix=integer(columnsz),
          Edssp=integer(columnsz),
          Ecaca=integer(columnsz),
#         RMSD=numeric(columnsz),
          Zhang1=numeric(columnsz),
          VSCORE1=numeric(columnsz),
          Bryant3=numeric(columnsz),
          Crease3=numeric(columnsz),
 #         VPDBatomN=integer(columnsz),
          SolvAtoms=integer(columnsz),
 #         SdX=numeric(columnsz),
 #         SdY=numeric(columnsz),
 #         SdZ=numeric(columnsz),
          Ea=numeric(columnsz),
          Eb=numeric(columnsz),
          Ec=numeric(columnsz))
#          colnames(Values)<- c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca",
#          "Zhang1","VSCORE1","Bryant3","Crease3","VPDBatomN","SolvAtoms1","SdX","SdY","SdZ","Ea","Eb","Ec")
          colnames(Values)<- c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca",
          "Zhang1","VSCORE1","Bryant3","Crease3","SolvAtoms1","Ea","Eb","Ec")
# Now convert the array of strings into a column vector of the appropriate type
          Values$Structure<-as.integer(Complete[ seq(from = 1, to = length(Complete), by = width ) ]) 
          Values$Time<-as.integer(Complete[ seq(from = 2, to = length(Complete), by = width ) ])
          Values$Tries<-as.integer(Complete[ seq(from = 3, to = length(Complete), by = width ) ])
          Values$BadBB<-as.integer(Complete[ seq(from = 4, to = length(Complete), by = width ) ])
          Values$Crashes<-as.integer(Complete[ seq(from = 5, to = length(Complete), by = width ) ])
          Values$ViolatedConstr<-Complete[ seq(from = 6, to = length(Complete), by = width ) ]
          Values$N<-as.integer(Complete[ seq(from = 7, to = length(Complete), by = width ) ])
          Values$Rgyr<-as.numeric(Complete[ seq(from = 8, to = length(Complete), by = width ) ])
          Values$HRgyr<-as.numeric(Complete[ seq(from = 9, to = length(Complete), by = width ) ])
          Values$NCdist<-as.numeric(Complete[ seq(from = 10, to = length(Complete), by = width ) ])
          Values$Rn<-as.numeric(Complete[ seq(from = 11, to = length(Complete), by = width ) ])
          Values$Cn<-as.numeric(Complete[ seq(from = 12, to = length(Complete), by = width ) ])
          Values$ASA<-as.numeric(Complete[ seq(from = 13, to = length(Complete), by = width ) ])
          Values$HASA<-as.numeric(Complete[ seq(from = 14, to = length(Complete), by = width ) ])
          Values$Helix<-as.integer(Complete[ seq(from = 15, to = length(Complete), by = width ) ])
          Values$Edssp<-as.integer(Complete[ seq(from = 16, to = length(Complete), by = width ) ])
          Values$Ecaca<-as.integer(Complete[ seq(from = 17, to = length(Complete), by = width ) ])
#        Values$RMSD<-as.integer(Complete[ seq(from = 18, to= length(Complete), by = width ) ])
          Values$Zhang1<-as.numeric(Complete[ seq(from = 18, to = length(Complete), by = width ) ])
          Values$VSCORE1<-as.numeric(Complete[ seq(from = 19, to = length(Complete), by = width ) ])
          Values$Bryant3<-as.numeric(Complete[ seq(from = 20, to = length(Complete), by = width ) ])
          Values$Crease3<-as.numeric(Complete[ seq(from = 21, to = length(Complete), by = width ) ])
#        Values$VPDBatomN<-as.integer(Complete[ seq(from = 22, to = length(Complete), by = width ) ])
          Values$SolvAtoms1<-as.integer(Complete[ seq(from = 23, to = length(Complete), by = width ) ])
#        Values$SdX<-as.numeric(Complete[ seq(from = 24, to = length(Complete), by = width ) ])
#        Values$SdY<-as.numeric(Complete[ seq(from = 25, to = length(Complete), by = width ) ])
#        Values$SdZ<-as.numeric(Complete[ seq(from = 26, to = length(Complete), by = width ) ])
          Values$Ea<-as.numeric(Complete[ seq(from = 27, to = length(Complete), by = width ) ])
          Values$Eb<-as.numeric(Complete[ seq(from = 28, to = length(Complete), by = width ) ])
          Values$Ec<-as.numeric(Complete[ seq(from = 29, to = length(Complete), by = width ) ])
           }            
     else if(width == 30) {
# NEW TRADES-2 LOG FORMAT WITH RMSD 
     Values<-data.frame(
          Structure=integer(columnsz), 
          Time=integer(columnsz),
          Tries=integer(columnsz),
          BadBB=integer(columnsz),
          Crashes=integer(columnsz),
          ViolatedConstr=integer(columnsz),
          N=integer(columnsz),
          Rgyr=numeric(columnsz),
          HRgyr=numeric(columnsz),
          NCdist=numeric(columnsz),
          Rn=numeric(columnsz),
          Cn=numeric(columnsz),
          ASA=numeric(columnsz),
          HASA=numeric(columnsz),
          Helix=integer(columnsz),
          Edssp=integer(columnsz),
          Ecaca=integer(columnsz),
          RMSD=numeric(columnsz),
          Zhang1=numeric(columnsz),
          VSCORE1=numeric(columnsz),
          Bryant3=numeric(columnsz),
          Crease3=numeric(columnsz),
# Single Number not used in computations
# VPDBatomN=integer(columnsz),
          SolvAtoms=integer(columnsz),
# These are not used in computations
#         SdX=numeric(columnsz),
#         SdY=numeric(columnsz),
#         SdZ=numeric(columnsz),
          Ea=numeric(columnsz),
          Eb=numeric(columnsz),
          Ec=numeric(columnsz))
#          colnames(Values)<- c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca",
#          "RMSD","Zhang1","VSCORE1","Bryant3","Crease3","VPDBatomN","SolvAtoms1","SdX","SdY","SdZ","Ea","Eb","Ec")
          colnames(Values)<- c("Structure","Time","Tries","BadBB","Crashes","ViolatedConstr","N","Rgyr","HRgyr","NCdist","Rn","Cn","ASA","HASA","Helix","Edssp","Ecaca"
          ,"RMSD","Zhang1","VSCORE1","Bryant3","Crease3","SolvAtoms1","Ea","Eb","Ec")
# Now convert the array of strings into a column vector of the appropriate type
          Values$Structure<-as.integer(Complete[ seq(from = 1, to = length(Complete), by = width ) ]) 
          Values$Time<-as.integer(Complete[ seq(from = 2, to = length(Complete), by = width ) ])
          Values$Tries<-as.integer(Complete[ seq(from = 3, to = length(Complete), by = width ) ])
          Values$BadBB<-as.integer(Complete[ seq(from = 4, to = length(Complete), by = width ) ])
          Values$Crashes<-as.integer(Complete[ seq(from = 5, to = length(Complete), by = width ) ])
          Values$ViolatedConstr<-Complete[ seq(from = 6, to = length(Complete), by = width ) ]
          Values$N<-as.integer(Complete[ seq(from = 7, to = length(Complete), by = width ) ])
          Values$Rgyr<-as.numeric(Complete[ seq(from = 8, to = length(Complete), by = width ) ])
          Values$HRgyr<-as.numeric(Complete[ seq(from = 9, to = length(Complete), by = width ) ])
          Values$NCdist<-as.numeric(Complete[ seq(from = 10, to = length(Complete), by = width ) ])
          Values$Rn<-as.numeric(Complete[ seq(from = 11, to = length(Complete), by = width ) ])
          Values$Cn<-as.numeric(Complete[ seq(from = 12, to = length(Complete), by = width ) ])
          Values$ASA<-as.numeric(Complete[ seq(from = 13, to = length(Complete), by = width ) ])
          Values$HASA<-as.numeric(Complete[ seq(from = 14, to = length(Complete), by = width ) ])
          Values$Helix<-as.integer(Complete[ seq(from = 15, to = length(Complete), by = width ) ])
          Values$Edssp<-as.integer(Complete[ seq(from = 16, to = length(Complete), by = width ) ])
          Values$Ecaca<-as.integer(Complete[ seq(from = 17, to = length(Complete), by = width ) ])
          Values$RMSD<-as.integer(Complete[ seq(from = 18, to= length(Complete), by = width ) ])
          Values$Zhang1<-as.numeric(Complete[ seq(from = 19, to = length(Complete), by = width ) ])
          Values$VSCORE1<-as.numeric(Complete[ seq(from = 20, to = length(Complete), by = width ) ])
          Values$Bryant3<-as.numeric(Complete[ seq(from = 21, to = length(Complete), by = width ) ])
          Values$Crease3<-as.numeric(Complete[ seq(from = 22, to = length(Complete), by = width ) ])
 #        Values$VPDBatomN<-as.integer(Complete[ seq(from = 23, to = length(Complete), by = width ) ])
          Values$SolvAtoms1<-as.integer(Complete[ seq(from = 24, to = length(Complete), by = width ) ])
 #        Values$SdX<-as.numeric(Complete[ seq(from = 25, to = length(Complete), by = width ) ])
 #        Values$SdY<-as.numeric(Complete[ seq(from = 26, to = length(Complete), by = width ) ])
 #        Values$SdZ<-as.numeric(Complete[ seq(from = 27, to = length(Complete), by = width ) ])
          Values$Ea<-as.numeric(Complete[ seq(from = 28, to = length(Complete), by = width ) ])
          Values$Eb<-as.numeric(Complete[ seq(from = 29, to = length(Complete), by = width ) ])
          Values$Ec<-as.numeric(Complete[ seq(from = 30, to = length(Complete), by = width ) ])
           }           
Values
}

##################################################################
# TRADES.Ellipsoid returns some basic information about how 
# Oblate and Prolate forms are partitioned by Rgyr and NCdist
# Plots a graph if called with plot=TRUE

TRADES.Ellipsoid<-function(logF, plot=TRUE) {

if (mode(logF$Data) == "list")
   logF<-logF$Data
   

if(length(logF$BaseNamesUnique) > 1) bn <- paste(substring(logF$BaseName,1,4),"_xxx",sep="")
else bn<-logF$BaseNamesUnique[1]


cat("Plotting Ellipsoid Parameters on Entire Ensemble..\n")
filename<-paste(bn,"_Ellipsoid.png",sep="")
   
if (plot == TRUE)
{
 png(filename, width=2000, height=2000)
 par(mfrow=c(2,2),cex=2.3)
}

EVolumes <- logF$values$Ea * logF$values$Eb * logF$values$Ec * 4 * 3.1415926 / 3
BProlate <-  logF$values$Ea > (logF$values$Eb + logF$values$Ec)
BOblate <- !BProlate
#EForm takes two logical vectors and converts them to -1 and +1 
EForm <- as.integer(BProlate) - as.integer(BOblate)
SignedVolumes <- EVolumes * EForm    
#Prolates are Positive, Oblates are negative volumes
Oblates<- SignedVolumes[SignedVolumes < 0]
Prolates<- SignedVolumes[SignedVolumes >= 0]
OVCEnergy <- logF$values$EnergyVC[SignedVolumes < 0]
PVCEnergy <- logF$values$EnergyVC[SignedVolumes >= 0]
ORgyr <- logF$values$Rgyr[SignedVolumes < 0]
PRgyr <- logF$values$Rgyr[SignedVolumes >= 0]
ONCdist <- logF$values$NCdist[SignedVolumes < 0]
PNCdist <- logF$values$NCdist[SignedVolumes >= 0]

#Plot two distribution curves of Oblate, Prolate, scaled 

#Plot Energy, Rgyr distributions of Oblate, Prolate
# Free Energy distributions split into Oblate, Prolate
# NCDist?
#Return Histogram Data
EllipsoidStats<-NULL

titlestring<-paste("Prolate:",length(Prolates),sep=" ")
EllipsoidStats$ProlateRgyr <-TRADES.histogram(PRgyr,rug=TRUE,peak=TRUE,
name=titlestring,xlab= "Rgyr",breaks=300, plot=FALSE)

EllipsoidStats$ProlateNCdist <-TRADES.histogram(PNCdist,rug=TRUE,peak=TRUE,
name=titlestring,xlab= "N - C Distance (Angstroms)",breaks=300,plot=FALSE)

xlimitsRgyr=c(0,EllipsoidStats$ProlateRgyr$data_Max)
xlimitsNCdist=c(0,EllipsoidStats$ProlateNCdist$data_Max)

TRADES.histogram(PRgyr,rug=TRUE,peak=TRUE,
name=titlestring,xlab= "Rgyr",breaks=300, xlimits=xlimitsRgyr, plot=plot)

TRADES.histogram(PNCdist,rug=TRUE,peak=TRUE,
name=titlestring,xlab= "N - C Distance (Angstroms)",breaks=300, xlim=xlimitsNCdist, plot=plot)

titlestring<-paste("Oblate: ",length(Oblates),sep=" ")
EllipsoidStats$OblateRgyr <- TRADES.histogram(ORgyr,rug=TRUE,peak=TRUE,
name=titlestring,xlab= "Rgyr",breaks=300, xlimits=xlimitsRgyr, plot=plot)

EllipsoidStats$OblateNCdist <-TRADES.histogram(ONCdist,rug=TRUE,peak=TRUE,
name=titlestring,xlab= "N - C Distance (Angstroms)",breaks=300, xlim=xlimitsNCdist, plot=plot)

NOblate <- sum(BOblate)
NProlate <- sum(BProlate)
PtoOratio <- NProlate / NOblate

if (plot == TRUE)
{
 dev.off() 
}

Ellipsoid_Summary<-list(NOblate = NOblate, NProlate=NProlate, PtoOratio = PtoOratio, EllipsoidStats=EllipsoidStats)
Ellipsoid_Summary
}





#-------------------------------------------------------
# TRADES.BoltzmannSearch tries to find the optimum
# bandwidth for Boltzmann state partitioning calculations
# and in turn, find the best, smallest ensemble posible
# This dumps a lot of graphs which should be inspected by eye 



########################################################################################################
# Main Function for Hydrodynamic Boltzmann Ensemble Refinement
# Requires a sample size that is divisible by 300000 for simplifying the factoring of groups
#


TRADES.BoltzmannSearch<-function(logF, T=25, plot=TRUE, ExperimentalRgyr=0, ExperimentalNCDist=0) {

if (mode(logF$Data) == "list")
   logF<-logF$Data

#calib_size <- 300000
sample_size<-length(logF$values$Structure)


#frac_size <- sample_size / calib_size

if (sample_size %% 300000 != 0) {
 cat("Boltzman Search cannot start\n - Adjust sample size to make it exactly divisible by 300000 please.\n")
 return(0)
}
 
interval_list = c(300, 400, 500, 600, 1000, 1200, 1500, 3000, 3750, 5000, 7500, 10000, 15000, 20000, 30000, 50000, 60000, 100000, 150000, 300000)
stru_per_Boltz_state = sample_size / interval_list

#maintain at least 10 structures per Boltzmann state for averaging. 
#This can go as low as 5 but averaging code fails if lower than that
trim <- stru_per_Boltz_state >= 10
interval_list<-interval_list[trim]
stru_per_Boltz_state<-stru_per_Boltz_state[trim]

#R object to collect all the results
Ensemble_List<-list()


#Start of iterations  - Selects a group size, finds the Boltzmann ensemble, then lowers the group size, repeats...
for (i in 1:length(interval_list)) { 
  cat("\n---------------------------\nIteration: ")
  cat(i)
  cat(" of ")
  cat(length(interval_list))
  cat("\nStructures Per Boltzmann State: ")
  cat(stru_per_Boltz_state[i])
  cat("\n")

#Call to the main Boltzmann procedure, which dumps out group-size dependent graphs , returns data 
  Ensemble <- TRADES.Boltzmann(logF, plot=plot, T = T, intervals=interval_list[i], file = TRUE)
  
#Returned Values in Ensemble:
#Ensemble$T
#Ensemble$SpBS
#Ensemble$intervals
#Ensemble$Boltzmann_States_Selected_Ensemble
# [,1] Struc Number, 
# [,2] Free VC Energy at input T 
# [,3] is Rgyr
# [,4] is NCDist 
# [,5] is the signed Boltzmann State value
#Ensemble$Boltzmann_States_Selected_Files
#Ensemble$Ensemble_Stats 
#                       $Energy - overall Energy Distribution, mean, peak, FWHM...
#                       $Rgyr
#                       $NCdist
#   - if >10 Prolate:   $ProlateEnergy  
#                       $ProlateRgyr    
#                       $ProlateNCdist  
#   - if >10 Oblate:    $OblateEnergy   
#                       $OblateRgyr
#                       $OblateNCdist
  
  Ensemble_List[[i]] <- Ensemble
  
} # End of Iterations
cat("Iterations Completed.\n")
list_len = i


##################################################################################################
# Loop over Ensemble List and 
# Collect information about Joint/Prolate/Oblate scores/fwhm

#scoring matrix has up to 3 rows per Boltzmann run, for all, prolate, oblate 
cnames=c("T","SpBS","States","Form","Chosen_Str","Energy_Mean","Energy_FWHM","Rgyr_Mean","Rgyr_FWHM","NCDist_Mean","NCDist_FWHM","Index")
e_rgy_values<-matrix(ncol=12, nrow=(list_len + 2) * 3, dimnames=list(NULL,cnames))



##################################################################
#Fill in the scoring structure, Draw the Histograms
j = 1
for (i in 1:list_len) {
  e_rgy_values[j,1]=Ensemble_List[[i]]$T
  e_rgy_values[j,2]=Ensemble_List[[i]]$SpBS
  e_rgy_values[j,3]=Ensemble_List[[i]]$States
  e_rgy_values[j,4]=0 # Prolate + Oblate
  e_rgy_values[j,5]=length(Ensemble_List[[i]]$Boltzmann_States_Selected_Ensemble[,1])
  e_rgy_values[j,6]=Ensemble_List[[i]]$Ensemble_Stats$Energy$data_Mean
  e_rgy_values[j,7]=Ensemble_List[[i]]$Ensemble_Stats$Energy$data_FWHM
  e_rgy_values[j,8]=Ensemble_List[[i]]$Ensemble_Stats$Rgyr$data_Mean
  e_rgy_values[j,9]=Ensemble_List[[i]]$Ensemble_Stats$Rgyr$data_FWHM
  e_rgy_values[j,10]=Ensemble_List[[i]]$Ensemble_Stats$NCdist$data_Mean
  e_rgy_values[j,11]=Ensemble_List[[i]]$Ensemble_Stats$NCdist$data_FWHM
  e_rgy_values[j,12]=i
  j = j + 1
  
  test<-Ensemble_List[[i]]$Ensemble_Stats$ProlateEnergy
  if (mode(test) != "NULL")
  { #Prolate Forms
    e_rgy_values[j,1]=Ensemble_List[[i]]$T
    e_rgy_values[j,2]=Ensemble_List[[i]]$SpBS
    e_rgy_values[j,3]=Ensemble_List[[i]]$States
    e_rgy_values[j,4]=1
    e_rgy_values[j,5]=Ensemble_List[[i]]$Ensemble_Stats$ProlateEnergy$data_N
    e_rgy_values[j,6]=Ensemble_List[[i]]$Ensemble_Stats$ProlateEnergy$data_Mean
    e_rgy_values[j,7]=Ensemble_List[[i]]$Ensemble_Stats$ProlateEnergy$data_FWHM
    e_rgy_values[j,8]=Ensemble_List[[i]]$Ensemble_Stats$ProlateRgyr$data_Mean
    e_rgy_values[j,9]=Ensemble_List[[i]]$Ensemble_Stats$ProlateRgyr$data_FWHM
    e_rgy_values[j,10]=Ensemble_List[[i]]$Ensemble_Stats$ProlateNCdist$data_Mean
    e_rgy_values[j,11]=Ensemble_List[[i]]$Ensemble_Stats$ProlateNCdist$data_FWHM
    e_rgy_values[j,12]=i
    j = j + 1
  }
  test<-Ensemble_List[[i]]$Ensemble_Stats$OblateEnergy
  if (mode(test) != "NULL")
  { #Oblate Forms
    e_rgy_values[j,1]=Ensemble_List[[i]]$T
    e_rgy_values[j,2]=Ensemble_List[[i]]$SpBS
    e_rgy_values[j,3]=Ensemble_List[[i]]$States
    e_rgy_values[j,4]=-1
    e_rgy_values[j,5]=Ensemble_List[[i]]$Ensemble_Stats$OblateEnergy$data_N
    e_rgy_values[j,6]=Ensemble_List[[i]]$Ensemble_Stats$OblateEnergy$data_Mean
    e_rgy_values[j,7]=Ensemble_List[[i]]$Ensemble_Stats$OblateEnergy$data_FWHM
    e_rgy_values[j,8]=Ensemble_List[[i]]$Ensemble_Stats$OblateRgyr$data_Mean
    e_rgy_values[j,9]=Ensemble_List[[i]]$Ensemble_Stats$OblateRgyr$data_FWHM
    e_rgy_values[j,10]=Ensemble_List[[i]]$Ensemble_Stats$OblateNCdist$data_Mean
    e_rgy_values[j,11]=Ensemble_List[[i]]$Ensemble_Stats$OblateNCdist$data_FWHM
    e_rgy_values[j,12]=i
    j = j + 1
  }
}

j = j - 1
# j is the number of rows to report 

if(length(logF$BaseNamesUnique) > 1) bn <- paste(substring(logF$BaseName,1,4),"_xxx",sep="")
else bn<-logF$BaseNamesUnique[1]

#truncate the n/a 
e_rgy_values<-e_rgy_values[1:j,]

cat("Writing out table of Boltzmann results..\n")
filename<-paste(bn,"_BoltzmannScanData_T",T,".csv",sep="")
write.csv(e_rgy_values, file=filename)

###########################################################################################
# Plots
# Draw the long plot of histograms to check quality of fits... 
#  Histograms Showing Best Hydrodynamic Ensembles chosen by Energy (Prolate or Oblate)


if (plot == TRUE)
{
 cat("Plotting Histograms of Selected Boltzmann Ensembles..\n")
 filename<-paste(bn,"_BoltzmannScanHist_T",T,".png",sep="")
 png(filename, width=1500, height=(500 * j ))
 par(mfrow=c(j,3),cex=1.3)

 for (i in 1:j) {
   fEnsemble_length <- e_rgy_values[i,5]
   fspbs<-e_rgy_values[i,2]
   fTemp <- e_rgy_values[i,1]
   fEnergy <- e_rgy_values[i,6]
   fNCdist <- e_rgy_values[i,10]
   fRgyr <- e_rgy_values[i,8]
   Index <- e_rgy_values[i,12]
   if  (e_rgy_values[i,4] == 0) { #Joint Row
     plot_values_Energy <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,2])
     plot_values_Rgyr <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,3])
     plot_values_NCdist <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,4])
     form = "All"
   } else
   if  (e_rgy_values[i,4] == 1) { #Prolate Row
     Prolates <- Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,5] > 0
     plot_values_Energy <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Prolates  ,2])
     plot_values_Rgyr <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Prolates  ,3])
     plot_values_NCdist <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Prolates  ,4])   
     form = "Prolate"
   } else
   if  (e_rgy_values[i,4] == -1) { #Oblate Row
     Oblates <- Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,5] <= 0
     plot_values_Energy <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Oblates ,2])
     plot_values_Rgyr <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Oblates ,3])
     plot_values_NCdist <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Oblates ,4])
     form = "Oblate"
   }
   name<-paste(bn," N:",fEnsemble_length," T:",fTemp," SpBS:",fspbs,sep="")
   TRADES.histogram(plot_values_Energy, plot=plot, name=name, xlabel="Energy (kcal/mol)", peak=TRUE, rug=TRUE)
   name<-paste("Energy:",round(fEnergy,1) ," Rgyr:",round(fRgyr,1) , sep="")
   TRADES.histogram(plot_values_Rgyr, plot=plot,name=name, xlabel="Rgyr (Angstroms)", peak=TRUE, rug=TRUE)
   name<-paste("NCdist:",round(fNCdist,1) ," Form:",form, sep="")
   TRADES.histogram(plot_values_NCdist, plot=plot, name=name, xlabel="N - C Distance (Angstroms)", peak=TRUE, rug=TRUE)
  }
  dev.off()
 }
 
 
 
#############################################################3 
#Draw the long Energy Plots

if(plot==TRUE) { #plot the summary graphs 
 filename<-paste(bn,"_BoltzmannScanEnergy_T",T,".png",sep="")
  png(filename, width=1500, height=(500 * (j + 1)))
  par(mfrow=c(j+1,3),cex=1.3)
  
# plot the big smoothScatter image for each run to compare to the starting ensemble
 Lab.palette <-colorRampPalette(c("white","lightyellow", "lightcyan","cyan", "lightskyblue", "lightseagreen", "yellowgreen" ,"yellow", "goldenrod", "orange", "orange4", "firebrick", "darkred", "red", "darkmagenta", "magenta", "hotpink", "pink","lightpink","white"), space = "Lab")
 ylimits<-c(min(logF$values$EnergyFreeVC), max(logF$values$EnergyFreeVC))
 Rxlimits<-c(min(logF$values$Rgyr), max(logF$values$Rgyr))
 NCxlimits<-c(min(logF$values$NCdist), max(logF$values$NCdist))
 #Plot the starting energy surfaces
 Starting_RgyrStats<-TRADES.histogram(as.vector(logF$values$Rgyr),plot=FALSE)
 Starting_EnergyStats<-TRADES.histogram(as.vector(logF$values$EnergyFreeVC),plot=FALSE)
 name<-paste("N:",sample_size," ",bn," T:",T,"C",sep="")
 smoothScatter(logF$values$NCdist ~ logF$values$Rgyr, main=name,ylim=NCxlimits, xlim=Rxlimits, xlab="Rgyr (Angstroms)", ylab="N - C Distance (Angstroms)", colramp=Lab.palette)
 abline(h=0)
 name<-paste("Rg Mean:", round(Starting_RgyrStats$data_Mean,2)," FWHM:", round(Starting_RgyrStats$data_FWHM,1), sep="")
 smoothScatter(logF$values$EnergyFreeVC ~ logF$values$Rgyr, ylab="Energy (kcal/mol)",main=name,xlim=Rxlimits, ylim=ylimits, xlab="Rgyr (Angstroms)", colramp=Lab.palette)
 abline(h=0)
 name<-paste("Energy Mean:",round(Starting_EnergyStats$data_Mean,2)," FWHM:",round(Starting_EnergyStats$data_FWHM,1),sep="")
 smoothScatter(logF$values$EnergyFreeVC ~ logF$values$NCdist, main=name,xlim=NCxlimits, ylim=ylimits, ylab="Energy (kcal/mol)", xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
 abline(h=0)

for (i in 1:j) {
   fEnsemble_length <- e_rgy_values[i,5]
   fspbs<-e_rgy_values[i,2]
   fTemp <- e_rgy_values[i,1]
   fEnergy <- e_rgy_values[i,6]
   fNCdist <- e_rgy_values[i,10]
   fRgyr <- e_rgy_values[i,8]
   fRgyr_FWHM<-e_rgy_values[i,9]
   fEnergy_FWHM<-e_rgy_values[i,7]
   Index <- e_rgy_values[i,12]
   if  (e_rgy_values[i,4] == 0) { #Joint Row
     plot_values_Energy <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,2])
     plot_values_Rgyr <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,3])
     plot_values_NCdist <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,4])
     form = "All"
   } else
   if  (e_rgy_values[i,4] == 1) { #Prolate Row
     Prolates <- Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,5] > 0
     plot_values_Energy <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Prolates  ,2])
     plot_values_Rgyr <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Prolates  ,3])
     plot_values_NCdist <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Prolates  ,4])   
     form = "Prolate"
   } else
   if  (e_rgy_values[i,4] == -1) { #Oblate Row
     Oblates <- Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[,5] <= 0
     plot_values_Energy <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Oblates ,2])
     plot_values_Rgyr <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Oblates ,3])
     plot_values_NCdist <- as.vector(Ensemble_List[[Index]]$Boltzmann_States_Selected_Ensemble[ Oblates ,4])
     form = "Oblate"
   }
 #Plot the energy surface of each increment
 name<-paste("Form: ",form," N:",fEnsemble_length," SpBS:",fspbs,sep="")
 smoothScatter(plot_values_NCdist  ~ plot_values_Rgyr , main=name, ylim=NCxlimits, xlim=Rxlimits, xlab="Rgyr (Angstroms)", ylab="N - C Distance (Angstroms)", colramp=Lab.palette)
 abline(h=0)
 name<-paste("Rg: ", round(fRgyr,2)," FWHM: ", round(fRgyr_FWHM,1), sep="")
 smoothScatter(plot_values_Energy ~ plot_values_Rgyr , main=name, ylab="Energy (kcal/mol)",xlim=Rxlimits, ylim=ylimits, xlab="Rgyr (Angstroms)", colramp=Lab.palette)
 abline(h=0)
 name<-paste("Energy Mean:", round(fEnergy,2)," FWHM: ", round(fEnergy_FWHM,1), sep="")
 smoothScatter(plot_values_Energy  ~ plot_values_NCdist , main=name, xlim=NCxlimits, ylim=ylimits, ylab="Energy (kcal/mol)", xlab="N - C Distance (Angstroms)", colramp=Lab.palette)
 abline(h=0)

 }
 dev.off()
} # Long Energy Plot done...

 
 
#######################################################
# Select Best Ensemble With Prolate / Oblate Metrics
# Optimize for Smallest (Rgyr FWHM) and lowest Energy.
#


# Exclude N < 100 structures in ensemble because FWHM are unreliable as are E values
remove_list <- e_rgy_values[,5] < 100
if (sum(as.numeric(remove_list)) != 0) {
  sel_matrix <- e_rgy_values[-which(remove_list),]
} else {
 sel_matrix <- e_rgy_values
}
# Exclude "All", keep only Prolate or Oblate data sets in contention
remove_list <- sel_matrix[,4] == 0
if (sum(as.numeric(remove_list)) !=0) { #keeps it from returning 0 later on for iForm.
  sel_matrix <- sel_matrix[-which(remove_list),]
}
# take the minima of the Rgyr FWHM and the Energy
fwhm_min<-which.min(sel_matrix[,9])
e_min<-which.min(sel_matrix[,6])

# Compute each iteration's difference from the minima
delta_fwhm<-sel_matrix[fwhm_min,9] - sel_matrix[,9]
delta_e<-sel_matrix[e_min,6] - sel_matrix[,6]
# Multiply the values
pickfn <- (delta_fwhm + delta_e)/2
picked <- which.max(pickfn)
iPicked <- sel_matrix[picked,12]
iForm <- sel_matrix[picked,4]
pForm = "_"
if (iForm == 1) { pForm = "Prolate" }
if (iForm == -1) { pForm = "Oblate" }


#####################################################
# Go back and check the other data in this iteration
# to see if the prolate/oblate separation is too close
# to call by energy difference. 
# - if two significant fractions with <1 Kcal/mol
# we report the mixture of the two

all_picked <- 0
pro_picked <- 0
ob_picked <- 0
Iteration <- e_rgy_values[which(e_rgy_values[,12]==iPicked),] 
for (k in 1:length(Iteration[,1])) {
  if (Iteration[k,4] == 0) { all_picked <-k }
  if (Iteration[k,4] == 1) { pro_picked <-k }
  if (Iteration[k,4] == -1) { ob_picked <-k }
}
# only use if there are both oblate and prolate and number greater than 100 
if (pro_picked!=0 && ob_picked!=0) {
  if (Iteration[ob_picked,5] >= 100 && Iteration[pro_picked,5] >=100) {
   # if the Energy diff between Ob, Prolate is less than 1 Kcal/mol, USE the ALL dataset
   if (abs(Iteration[pro_picked,6] - Iteration[ob_picked,6]) < 1.0) {
    #overwrite the sel_matrix with the all data
    sel_matrix[picked,] <- Iteration[all_picked,]
    iForm <- 0
    pForm<- "Prolate+Oblate"
    }
  }
}

#######################################################
# Show the selection curves 
#

if (plot == TRUE) { 
 #plot the summary graphs 
 filename<-paste(bn,"_BoltzmanScanIter_T",T,".png",sep="")
 png(filename, width=1500, height=1000 )
 par(mfrow=c(2,3),cex=2)
 plot(sel_matrix[,8] ~ sel_matrix[,12],type="s", xlab="Iteration", ylab="Mean Rgyr  (Angstroms)")
 plot(sel_matrix[,6] ~ sel_matrix[,12],type="s", xlab="Iteration", ylab="Mean Energy (kcal/mol)")
 plot(sel_matrix[,5] ~ sel_matrix[,12],type="s", xlab="Iteration", ylab="Ensemble Size")
 plot(sel_matrix[,9] ~ sel_matrix[,12],type="s", xlab="Iteration", ylab="Rgyr FWHM")
 plot(pickfn ~ sel_matrix[,12],type="s", xlab="Iteration", ylab="Optimal Rgyr FWHM & Energy")
 dev.off()
}
 
 

cat("Optimal Ensemble is ",pForm, " Iteration: ",iPicked," with SpBS=",sel_matrix[picked,2]," and ",round(sel_matrix[picked,5])," Structures\n",sep="")
cat("Writing File of these Structure Names\n")

#OK export a file with the list of structures...
# Extract the filenames that belong to the prolate or oblate sub ensemble...

if (iForm == 1) { #Prolate Row
     Prolates <- Ensemble_List[[iPicked]]$Boltzmann_States_Selected_Ensemble[,5] > 0
     file_list <- Ensemble_List[[iPicked]]$Boltzmann_States_Selected_Files[ Prolates ]
   } else
   if  (iForm == -1) { #Oblate Row
     Oblates <- Ensemble_List[[iPicked]]$Boltzmann_States_Selected_Ensemble[,5] <= 0
     file_list <- Ensemble_List[[iPicked]]$Boltzmann_States_Selected_Files[ Oblates ]
   } else
    if  (iForm == 0) { #All row
     file_list <- Ensemble_List[[iPicked]]$Boltzmann_States_Selected_Files
    }   


cat(" Head of filenames\n")
cat(head(file_list),sep="\n")
cat(" Tail of filenames\n")
cat(tail(file_list),sep="\n")
cat(" length of filenames:\n")
cat(length(file_list))
cat("\n")
filename<-paste(bn,"Ensemble_",sel_matrix[picked,2],"_SpBS_T",T,"_",pForm,".txt",sep="")
fileCon<-file(filename)
writeLines(file_list, fileCon)
flush(fileCon)
close(fileCon)


cat("Optimal Ensemble Parameters:\n")
cat(sel_matrix[picked,])
cat("\n")

values <- list(Ensemble_List=Ensemble_List, Ensemble_Summary=e_rgy_values, Optimal_Ensemble=sel_matrix[picked,], Optimal_Files=file_list)
values
} #End of TRADES.BoltzmannSearch


  



############################################################################################################################
############################################################################################################################
# Boltzmann function using Hydrodynamic Ellipsoid state partitioning into groups and mean VC energy per group
# After computing Botlzmann probability of each group, retrieves representative low-energy 3D structures from the sample set
# passed in logF as a selected ensemble.  
# This function is driven by TRADES.BoltzmannSearch but can be called
# independently, with variations in Temperature once the best "intervals" rage is known
# "intervals" is the number of groups to split the sample into.
# sample size is retrieved from logF
# Boltzmann state group size (spBS) is sample size / intervals 
# So if you see a spBS - divide it by sample size to get intervals, then pass that value into the function.

TRADES.Boltzmann<-function(logF, T=25, intervals, plot=TRUE, file=TRUE) {

if (mode(logF$Data) == "list")
   logF<-logF$Data

   
#############################################################################
# Set up the Boltzmann State Partition Function - Groups of Ellipsoidal Volume + NCdist vector
# Oblate elliposids have their sign flipped
# 
# Ea*Eb*Ec*(4/3)(pi) + NCdist is chosen as the state function
# Sensitive to volumetric changes with small increments in NCdist
# And the function scales crudely like the states of gas particles in a baloon of radius Rgyr

cat("Setting up Hydrodynamic Boltzmann States..\n")
EVolumes <- logF$values$Ea * logF$values$Eb * logF$values$Ec * 4 * 3.1415926 / 3
BProlate <-  logF$values$Ea > (logF$values$Eb + logF$values$Ec)
BOblate <- !BProlate
#EForm takes two logical vectors and converts them to -1 and +1 
EForm <- as.integer(BProlate) - as.integer(BOblate)
States <- (EVolumes + logF$values$NCdist) * EForm
#This adds on the NC distance and flips sign for Oblates    
#Prolates are Positive Boltzmann states, Oblates are negative Boltzmann states
# States object remains unsorted and in numerical order of logF

#Find the boundary values for Boltzmann States
EStatesP <- States[BProlate]
EStatesO <- States[BOblate]
minStateP <- min(EStatesP)
maxStateP <- max(EStatesP)
minStateO <- min(EStatesO)
maxStateO <- max(EStatesO)

#cat("Omin,Omax,Pmin,Pmax\n")
#cat(minStateO)
#cat("\n")
#cat(maxStateO)
#cat("\n")
#cat(minStateP)
#cat("\n")
#cat(maxStateP)
#cat("\n")


cat("Computing Individual and Group Energies..\n")
####################################################################################
# Recalculate VC Energy according to input T
# delta S is scaled according to T by T/298
#
R<- 1.9858775
T=T+273.15
kt = (R * T)/1000
CorrectedEnergy <- (logF$values$EnergyVSCORE) + (kt * logF$values$EntropyCrease)
#recompute the energy in case a change in T is in the input parameters

#####################################################################################
# Set up a structure to sort given Corrected Energy & Boltzman State values
# and group width "iwidth"

structures<-length(logF$values$Structure)
if (missing(intervals)) intervals <- structures/1000

#Default is 1000 structures per Boltzman State
#Cut the States into intervals - fewer = better averaging.
iwidth<-structures / intervals

SortSet <- matrix(States, nrow=structures, ncol=1)
#sets up a copy of the Boltzmann States object for sorting and grouping grouping
SortSet<-cbind(SortSet, CorrectedEnergy)
SortSet<-cbind(SortSet,logF$values$Structure)
#connect States, Energy and Structure number for sorting by State value
# then sort the matrix by State value... 
Sorted<-SortSet[order(SortSet[,1], decreasing=FALSE),]
#Sorted now is sorted by Boltzmann state value, oblates are -, prolates are +

#Select intervals of iWidth
cuts<-seq.int(from = 1, to = structures, by=iwidth)

#Ranges of States  in ascending order
Ncuts<-Sorted[cuts,1]
Zsplit<-split(Sorted[,2],f=findInterval(Sorted[,1],Ncuts))
#Zsplit now holds the groups of size iwidth, sorted by Boltzmann stae 
lapply(Zsplit,mean)->Energies
#This computes the mean energy of each group
Energies<-sapply(Energies,as.numeric)
#This removes the statistical R fluff and just returns the mean energy value for each 

cat("Computing Boltzmann Group Probabilities..\n")
#########################################################################################
# Do the Partition function to compute Z, then Probabilities with the Boltzmann equation 
# Energies in kcal/mol
# Energies are mean values per Hydrodynamic Boltzmann state group
#   intervals is the number of Boltzmann states for Z
#   iwidth is the number of structures per Boltzmann state (spBS) 

Z <- sum((exp (-Energies/kt)))
Probabilities = (exp (-Energies/kt)) / Z

########################################################################################
# Now we have the probabilities for each group
# Next - these need to be converted back to individual representative structures 
# We will compute some group boundaries, and normalize the Probabilities according to the
# x-width of each group
# We will build an object States_Ensemble to hold information for retrieving
# a probability-proportionate number of low-energy group structures out of logF
# (and intentionally try to skip the high-energy group structures!)

Groups<-as.matrix(Probabilities)
Groups<-cbind(Groups,Ncuts)
#Ncuts has the maximum state (upper limit) for each interval
dimnames(Groups)<-list(NULL,c("Probabilities","States"))
Boundaries<-cbind(Groups,c(Groups[2:intervals,2],maxStateP))
#cat(Boundaries[,2])
#cat("\n")
#cat(Boundaries[,3])
#cat("\n")

#Boundaries is used assemble a matrix with the max value of the state interval at Boundaries[,2]
# and minimum values of the state intervals at the topmost value Boundaries[,3]
# as each group spans a different width on the Boltzmann state axis, we convert probabilities to 
# Heights associated with the midpoint of each group on the Boltzmann state axis
Heights<-Boundaries[,1]*intervals/(Boundaries[,3]-Boundaries[,2])
#Heights are the probability values (i.e. in Boundaries[,1])  adjusted for the graph area of the group bin,
# to plot it at midway between upper and lower bounds of group bin
States_Mid<-Boundaries[,2]+ (Boundaries[,3]-Boundaries[,2])/2
# now GroupBoxes[,5] has midpoint values for X axis ...
norm_freq <- Heights/max(Heights)
# norm_freq is the normalized height value (group area probability)

States_Ensemble<-Boundaries[,3]
States_Ensemble<-cbind(States_Ensemble,Boundaries[,2])
States_Ensemble<-cbind(States_Ensemble,Heights)
States_Ensemble<-cbind(States_Ensemble,States_Mid)
States_Ensemble <-cbind(States_Ensemble, norm_freq)

################################################################################################################
#NOW States_Ensemble has organized the Ensemble group information for retrieval from logF with:
#[,1] Upper X value, [,2] Lower X value, [,3] corrected Prob Y at group midpoint, [,4] midpoint X value, [,5], normalized Y height 
#and again X is Boltzmann state space which is a signed hydrodynamic quantity (ellipsoid volume + NCdist) + for prolate, - for oblate

cat("Selecting Representative Low-Energy Structures..\n")
################################################################################################
#Next to select low-energy structure representatives from ensemble group
#Using the States_Ensemble information to pull matching rows of data out of logF

Boltzmann_States_Selected_Ensemble<-NULL
Boltzmann_States_Selected_Files<-NULL
#cat(head(States_Ensemble,1))
for (i in 1 : length(States_Ensemble[,1]))  {  #loop over every state in the Boltzmann state space 
#cat("\nProbability Value:",States_Ensemble[i,5]," at i=",i,"\n")
   if (States_Ensemble[i,5] > 0.01) {   #If significant Boltzmann probability, take representative low-energy structures
        # we use the normalized Y height to choose structure groups having group probability more than 1%. 
        as.logical(States >=States_Ensemble[i,2])-> left
        as.logical(States <=States_Ensemble[i,1])-> right
        States_mid<-(right & left)
        # this grabs just the structures out of the original States structure that fall into the group
        logF$values$Structure[States_mid] -> States_Interval_Structures
        # now we have the Structure numbers matching the group
        Num_Interval <- length(States_Interval_Structures)
        # this is the number of structures in the interval
        #changed to corrected energy values from T=25...
        States_Interval_Structures <- cbind(States_Interval_Structures,CorrectedEnergy[States_mid])
        # Now we have [,1] Struc Number, [,2] Free VC Energy at input T
        States_Interval_Structures <- cbind(States_Interval_Structures,logF$values$Rgyr[States_mid])
        # [,3] is Rgyr
        States_Interval_Structures <- cbind(States_Interval_Structures,logF$values$NCdist[States_mid])
        # [,4] is NCDist
        States_Interval_Structures <- cbind(States_Interval_Structures,States[States_mid])
        # [,5] is the signed Boltzmann State value
        States_Interval_Filenames <- logF$values$Filenames[States_mid]
        # We keep the filename strings out of the matrix for sorting purposes...
        Sort_Group<-as.matrix(States_Interval_Structures,ncol=6,nrow=Num_Interval)
        dimnames(Sort_Group)<-list(NULL,c("Structure","VC Energy","Rgyr","NCDist","Boltzmann State"))
        #We are going to sort the group by  VC energy
        # and take by priority the lowest energy structures 
        # and accumulate them in numbering according to the frequency 0.01 = 5 structures
        extract <- iwidth/2
        States_Int_Sorted<-Sort_Group[order(Sort_Group[,2], decreasing=FALSE),]
        States_Filenames_Sorted<-States_Interval_Filenames[order(Sort_Group[,2],decreasing=FALSE)] #had to be sorted separately outside of matrix
        #Now that they are sorted by VC energy, we can extract the structures that are low energy and append
        # to the final ensemble as best representatives of the group
        for (j in 1 : as.integer(round(States_Ensemble[i,5] * extract))) {
                   #this for loop counts up to the normalized y height * 1/2 group width
                   # - on each j, it selects the next low-energy structures from the group
                   # The Upper bound is frequency * 1/2 the group size 
                   # so that this never chooses from the "high-energy" half of the group 
                   # and leaves them out for the final ensemble
                   # So for a single group with Boltzmann probability 1, this will only choose the lower 1/2 of structures
                   # -NOTE -some structures can still be selected that have VC Energy > 0. 
           Boltzmann_States_Selected_Ensemble<-rbind(Boltzmann_States_Selected_Ensemble, States_Int_Sorted[j,])
           Boltzmann_States_Selected_Files<-c(Boltzmann_States_Selected_Files, States_Filenames_Sorted[j])
        }   
        #This makes up an ensemble by choosing low-energy structures from the best probability groups
        #The final ensemble can be made up of a combination of prolate/oblate structures, which are recorded.
   }
}
###################################################################################
#Now we have collected the representative structures and filenames into the objects:
#Boltzmann_States_Selected_Ensemble
#Boltzmann_States_Selected_Files

cat("Plotting..\n");
###################################################################################
# Next Plot Selected Ensemble Values and accumulate Total/Oblate/Prolate statistics
#

if (plot == TRUE) {

  if (file==TRUE) {
    filename<-paste(logF$BaseName,"_",iwidth,"_SpBoltzmann_State_T",round(T-273,0),".png",sep="")
    png(filename=filename,width=4000, height=4000)
    par(mfrow=c(2,2),cex=5.6, lwd=2)
   }
  else par(mfrow=c(2,2))


if(length(logF$BaseNamesUnique) > 1) bn <- paste(substring(logF$BaseName,1,4),"_xxx",sep="")
else bn<-logF$BaseNamesUnique[1]
name<-paste(bn," SpBS:",iwidth," Ensemble Size:",length(Boltzmann_States_Selected_Ensemble[,1]),sep="")

##############################################################################################################################
#These are graphical gynmastics to try to get the energy curve scaled and plotted overtop of the Boltzmann state distribution.
#We want the y axis zero line to be in the right place on the curve and scale the whole thing top-bottom to fit the graph.
#There is no scale on the energy axis, but it is "to scale". 
#move to zero
Energies_Norm<- (Energies - max(Energies))
#scale to 1
Energies_Norm<-  (Energies_Norm) / min(Energies_Norm)
#move back to zero
pos<-which(Energies > 0)
neg<- which(Energies <= 0)
Energies_lt0 <- Energies_Norm[neg]
Energies_gt0 <- Energies_Norm[pos]
if (length(Energies_gt0) == 0)
 offset = 0
if (length(Energies_lt0) == 0)
 offset = - 1
else offset = - (Energies_lt0[1] + Energies_lt0[length(Energies_lt0)])/2
Energies_Norm <- Energies_Norm + offset


########################################################################################################
#Now plot the Boltzman states, overlaid with Energy curve in blue and selected ensemble in red
#Zee "moustache" plot with oblate/prolate separated on left/right of Y axis
t<-hist(States,breaks=1000,plot=FALSE)
t_Peak<-t$density[which.max(t$density)]
t$density <- t$density/t_Peak
par(col="black")
ylimits<-c(min(-Energies_Norm),1)
xlimits<-c(min(States_Ensemble[,4]),max(States_Ensemble[,4]))
plot(t$density ~ t$mids,type="l",ylim=ylimits, xlim=xlimits, main=name, ylab="Normalized Values", xlab="Boltzmann State Groups (Ellipsoid Vol. + N:C Dist.)")
leg_x <- (xlimits[2]-xlimits[1])/2
leg_y <- 0.9 
legend(leg_x,leg_y, c("Initial States","Selected States","Mean State Energy"), lwd=2, cex=0.7, bty="n", pt.bg="white", lty=1, col = c("black","red","blue"))
abline(h=0)

#THIS curve could use some modifications for the part near the Y axis so the curve is split
#into oblate/prolate parts.  
par(col="blue")
lines(-Energies_Norm  ~ States_Ensemble[,4])
par(col="red")
lines(States_Ensemble[,5] ~ States_Ensemble[,4])
par(col="black")


#Boltzmann_States_Selected_Ensemble:
# [,1] Struc Number, 
# [,2] Free VC Energy at input T 
# [,3] is Rgyr
# [,4] is NCDist 
# [,5] is the signed Boltzmann State value


# Now do Energy Before/After Histogram Split out Prolate / Oblate Components

Ensemble_Stats<-list()
Ensemble_Stats$Energy<-TRADES.histogram(as.vector(Boltzmann_States_Selected_Ensemble[,2]), plot=FALSE)
name_e<-paste("Energy Mean:", round(Ensemble_Stats$Energy$data_Mean,1)," FWHM:", round(Ensemble_Stats$Energy$data_FWHM,1), sep="")

#Energy Distribution of Starting Ensemble - Black Line
e<-density(CorrectedEnergy)
e_Peak<-e$y[which.max(e$y)]
e$y <- e$y/e_Peak
plot(e,xlab="Energy (kcal/mol)",ylab="Normalized Frequency", main=name_e)

#Energy Distribution of Selected Ensemble - Red Line (Combined prolate + oblate)
f<-density(Boltzmann_States_Selected_Ensemble[,2])
f_Peak<-f$y[which.max(f$y)]
f$y <- f$y/f_Peak
par(col="red")
lines(f)

#Split out contributions:
prolate_rows <- Boltzmann_States_Selected_Ensemble[,5] > 0
oblate_rows <- Boltzmann_States_Selected_Ensemble[,5] <= 0
e_prolate <- Boltzmann_States_Selected_Ensemble[ prolate_rows , 2]
e_oblate <- Boltzmann_States_Selected_Ensemble[ oblate_rows, 2]
 
#Energy Distribution of Selected Low-Energy Prolate Ensemble Structures
#Chosen from Boltzmann Hydrodynamic Group with probability > 1%
if (length(e_prolate) > 10) {
cat("Prolate Ensemble Members Found:",length(e_prolate),"\n")
  f_prolate<-density(e_prolate)
  f_prolate_Peak<-f_prolate$y[which.max(f_prolate$y)]
  f_prolate$y <- f_prolate$y/f_prolate_Peak
  par(col="purple")
  lines(f_prolate)
#log the energy stats for Prolate
Ensemble_Stats$ProlateEnergy<-TRADES.histogram(as.vector(e_prolate), plot=FALSE)
}

#Energy Distribution of Selected Low-Energy Oblate Ensemble Structures
#Chosen from Boltzmann Hydrodynamic Group with probability > 1%
if (length(e_oblate) > 10) {
cat("Oblate Ensemble Members Found:",length(e_oblate),"\n")
  f_oblate<-density(e_oblate)
  f_oblate_Peak<-f_oblate$y[which.max(f_oblate$y)]
  f_oblate$y <- f_oblate$y/f_oblate_Peak
  par(col="orange")
  lines(f_oblate)
#log the energy stats for Oblate
Ensemble_Stats$OblateEnergy<-TRADES.histogram(as.vector(e_oblate), plot=FALSE)
}

par(col="black")



#Now the Rgyr Before/After Histogram

Ensemble_Stats$Rgyr<-TRADES.histogram(as.vector(Boltzmann_States_Selected_Ensemble[,3]), plot=FALSE)
name_r<-paste("Rgyr Mean:", round(Ensemble_Stats$Rgyr$data_Mean,2)," FWHM:", round(Ensemble_Stats$Rgyr$data_FWHM,1), sep="")


#Plot the cumulative Rgyr Histogram:
r<-density(logF$values$Rgyr)
r_Peak<-r$y[which.max(r$y)]
r$y <- r$y/r_Peak
plot(r,xlab="Rgyr (Angstroms)",ylab="Normalized Frequency", main=name_r)
s<-density(Boltzmann_States_Selected_Ensemble[,3])
s_Peak<-s$y[which.max(s$y)]
s$y <- s$y/s_Peak
par(col="red")
lines(s)


#Split out contributions:
r_prolate <- Boltzmann_States_Selected_Ensemble[ prolate_rows , 3]
r_oblate <- Boltzmann_States_Selected_Ensemble[ oblate_rows, 3]
 
#Rgyr Distribution of Selected Low-Energy Prolate Ensemble Structures
#Chosen from Boltzmann Hydrodynamic Group with probability > 1%
if (length(r_prolate) > 10) {
  s_prolate<-density(r_prolate)
  s_prolate_Peak<-s_prolate$y[which.max(s_prolate$y)]
  s_prolate$y <- s_prolate$y/s_prolate_Peak
  par(col="purple")
  lines(s_prolate)
#log the energy stats for Prolate
Ensemble_Stats$ProlateRgyr<-TRADES.histogram(as.vector(r_prolate), plot=FALSE)
}

#Rgyr Distribution of Selected Low-Energy Oblate Ensemble Structures
#Chosen from Boltzmann Hydrodynamic Group with probability > 1%
if (length(r_oblate) > 10) {
  s_oblate<-density(r_oblate)
  s_oblate_Peak<-s_oblate$y[which.max(s_oblate$y)]
  s_oblate$y <- s_oblate$y/s_oblate_Peak
  par(col="orange")
  lines(s_oblate)
#log the energy stats for Oblate
Ensemble_Stats$OblateRgyr<-TRADES.histogram(as.vector(r_oblate), plot=FALSE)
}


par(col="black")




#Now the NCDist Before/After Histogram:

Ensemble_Stats$NCdist<-TRADES.histogram(as.vector(Boltzmann_States_Selected_Ensemble[,4]), plot=FALSE)
name_nc<-paste("N-C dist. Mean:", round(Ensemble_Stats$NCdist$data_Mean,1)," FWHM:", round(Ensemble_Stats$NCdist$data_FWHM,1), sep="")


n<-density(logF$values$NCdist)
n_Peak<-n$y[which.max(n$y)]
n$y <- n$y/n_Peak
plot(n,xlab="N - C Distance (Angstroms)",ylab="Normalized Frequency", main=name_nc)
o<-density(Boltzmann_States_Selected_Ensemble[,4])
o_Peak<-o$y[which.max(o$y)]
o$y <- o$y/o_Peak
par(col="red")
lines(o)

#Split out contributions:
n_prolate <- Boltzmann_States_Selected_Ensemble[ prolate_rows , 4]
n_oblate <- Boltzmann_States_Selected_Ensemble[ oblate_rows, 4]
 
#NCDist Distribution of Selected Low-Energy Prolate Ensemble Structures
#Chosen from Boltzmann Hydrodynamic Group with probability > 1%
if (length(n_prolate) > 10) {
  o_prolate<-density(n_prolate)
  o_prolate_Peak<-o_prolate$y[which.max(o_prolate$y)]
  o_prolate$y <- o_prolate$y/o_prolate_Peak
  par(col="purple")
  lines(o_prolate)
#log the energy stats for Prolate
Ensemble_Stats$ProlateNCdist<-TRADES.histogram(as.vector(n_prolate), plot=FALSE)
}

#NCDist Distribution of Selected Low-Energy Oblate Ensemble Structures
#Chosen from Boltzmann Hydrodynamic Group with probability > 1%
if (length(n_oblate) > 10) {
  o_oblate<-density(n_oblate)
  o_oblate_Peak<-o_oblate$y[which.max(o_oblate$y)]
  o_oblate$y <- o_oblate$y/o_oblate_Peak
  par(col="orange")
  lines(o_oblate)
#log the energy stats for Oblate
Ensemble_Stats$OblateNCdist<-TRADES.histogram(as.vector(n_oblate), plot=FALSE)
}


par(col="black")


if (file == TRUE) dev.off()
}
cat("TRADES.Boltzmann: Ensemble Completed.\n")
##############################################################################################
#Plotting is done
#Return the selected ensemble, structure file list and stats accrued
#

Boltzmann_All<-list(T = T, SpBS = iwidth, States = intervals, Boltzmann_States_Selected_Ensemble = Boltzmann_States_Selected_Ensemble, 
Boltzmann_States_Selected_Files = Boltzmann_States_Selected_Files, 
Ensemble_Stats = Ensemble_Stats)


Boltzmann_All
} #End of TRADES.Boltzmann

 
