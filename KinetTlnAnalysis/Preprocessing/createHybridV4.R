#Create a hybrid dataset from a normalized steady state and the direct lig PAL-seq data
#V2 fixes a bug in version 1 and allows for comand line args. 
#V3 tries a local smoothing algorithm.
#V4 changes the previous versions so that the datasets are not rescaled to the PAL-seq data.
# Also for V4: bug fixed, really 30-70 matching.

args <- commandArgs(trailingOnly = TRUE)

#For debugging:
# args <- c("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Data/Run3AllData/TailsUridylation20190613_ShortTailMod_NormModHMM/miR155_uninduced_R0/tail_lengths_annot_end_only_CLEANED20190615_reformat.txt", "/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7.txt", "/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190610.txt")

# args <- c("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Data/Run2AllDataReProcess/CombinedData/miR-155_uninduced/tail_lengths_annot_end_only_reformat_UPDATE20190305.txt","/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v6.txt","test.txt")

library(tidyverse)

ShortTailNorm <- function(rowX, WindowSize){
	adjustment = sum(rowX[(252 + CUTOFF - WindowSize + 1):(252 + CUTOFF + WindowSize + 1)])/sum(rowX[(CUTOFF - WindowSize + 1):(CUTOFF + WindowSize + 1)])
	rowX[1:(CUTOFF)] = rowX[253:(253+CUTOFF-1)]
	rowX[1:(CUTOFF)] = rowX[1:(CUTOFF)]/adjustment
	# rowX[1:251] = rowX[1:251]/sum(rowX[1:251])*rowX['expr'] ##Expression normalization.
	return(rowX)
}
####### UNCOMMENT TO RUN ##########
DIR <- read_tsv(args[1]) #DIR direct 3p Lig dataset
SPL <- read_tsv(args[2],col_names = FALSE)
####### UNCOMMENT TO RUN ##########

## TROUBLESHOOTING CODE, UNCOMMENT TO TEST ##
# args <- c("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Data/Run2AllDataReProcess/CombinedData/miR-155_uninduced/tail_lengths_annot_end_only_reformat.txt","/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v6.txt","test.txt")
# RowDIR = c(-(-20:20)^2  + 400,c(0:209))
# names(RowDIR) = c(paste0("tl",0:250))
# DIR <- bind_rows(RowDIR) %>% add_column(accession = "TestRNA", .before = "tl0")
# RowSPL = c(0:250)
# names(RowSPL) = c(paste0("tl",0:250))
# SPL <- bind_rows(RowSPL) %>% add_column(accession = "TestRNA", .after = "tl250")
## TROUBLESHOOTING CODE, UNCOMMENT TO TEST ##

colnames(SPL) <- c(paste0("SPL_",0:250),"accession")
SPL <- SPL[,c(252,1:251)] #make accession the first column
DIR <- DIR[which(rowSums(DIR[,-1]) >= 50),]

print(dim(DIR))
print(dim(SPL))

CUTOFF = 50
WindowSize = 20

#Normalization
SPL$expr           = rowSums(SPL[,-1])

SPLDIR             = inner_join(SPL,DIR,by="accession")
SPLDIR[,-1]        = t(apply(SPLDIR[,-1],1,function(x){ShortTailNorm(x, WindowSize)}))
print(dim(SPLDIR))

SPLDIR_write <- SPLDIR[,c(2:252,1)] #switch the accession code back to the last column
 
write_delim(SPLDIR_write,path=args[3],delim = "\t",col_names = FALSE)