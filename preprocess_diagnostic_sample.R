rm(list=ls())

library(minfi)
library(limma)
library(openxlsx)
library(optparse)
library(Rtsne)
library(RSpectra)

source(file.path("R","MNPprocessIDAT_functions.R"))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="IDAT file location with basename", metavar="character"),
  make_option(c("-b", "--batch_correct"), type="character", default="results/ba.coef.RData",
              help="Location of batch correction coefficient RData file [default= %default]", metavar="character"),
  make_option(c("-t", "--sample_type"), type="character", default="FFPE", 
              help="Type of Sample (FFPE/Frozen) [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_idat_base = opt$file
sample_type = opt$sample_type
ba_corr = opt$batch_correct

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Input file must be supplied!", call.=FALSE)
}

message("reading in IDAT file ...",Sys.time())
RGset <- read.metharray(input_idat_base,verbose=TRUE)

message("running normalization ...",Sys.time())
Mset <- MNPpreprocessIllumina(RGset)

message("probe filtering ...",Sys.time())
amb.filter <- read.table(file.path("filter","amb_3965probes.vh20151030.txt"),header=F)
epic.filter <- read.table(file.path("filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
snp.filter <- read.table(file.path("filter","snp_7998probes.vh20151030.txt"),header=F)
xy.filter <- read.table(file.path("filter","xy_11551probes.vh20151030.txt"),header=F)
rs.filter <- grep("rs",rownames(Mset))
ch.filter <- grep("ch",rownames(Mset))

# filter CpG probes
remove <- unique(c(match(amb.filter[,1], rownames(Mset)),
                   match(epic.filter[,1], rownames(Mset)),
                   match(snp.filter[,1], rownames(Mset)),
                   match(xy.filter[,1], rownames(Mset)),
                   rs.filter,
                   ch.filter))
Mset_filtered <- Mset[-remove,]

message("performing batchadjustment ...",Sys.time())
methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)

load(ba_corr)
#methy.ba = 2^log2((methy + 1) + get(sample_type, methy.coef))
#unmethy.ba = 2^log2((unmethy + 1) + get(sample_type, unmethy.coef))
methy.ba = 2^(log2(methy + 1) + get(sample_type, methy.coef))
unmethy.ba = 2^(log2(unmethy + 1) + get(sample_type, unmethy.coef))

# recalculate betas, illumina like
diag_betas <- methy.ba / (methy.ba +unmethy.ba +100)
diag_betas <- as.data.frame(t(diag_betas))

message("preprocessing finished ...",Sys.time())

########################
message("tSNE starting ...",Sys.time())
source(file.path("R","RSpectra_pca.R"))

message("loading preprocessed data for tSNE ...",Sys.time())
load(file.path("results","betas.ba.RData"))

# methylation classes
y = append(anno$`methylation class:ch1`, "DIAGNOSTIC")
y <- as.factor(y)

# Implement a way to add the "diag_betas" to this betas object
betas = rbind(betas, diag_betas)

# sd filtering to 32k probes
betas <- betas[,order(-apply(betas,2,sd))[1:32000]]


# calculate first 94 PCs
pca <- prcomp_svds(betas,k=94)

# calculate tSNE
res <- Rtsne(pca$x,pca=F,max_iter=2500,theta=0,verbose=T)

tsne_df = data.frame(tsne1 = res$Y[,1], tsne2 = res$Y[,2], grp = y)
write.csv(tsne_df, file.path("results","tsne_diag.csv"))

# scatterplot tSNE
pdf(file.path("results","tsne_diag.pdf"))
plot(res$Y,pch=19,col=y)
dev.off()
