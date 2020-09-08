### This is to analyze some scRNA-seq data for shiyan in CRC          
### Fetch supplementary files from GEO
### Do differential expression analysis for Y1 high vs low groups
###
library(GEOquery);
setwd('/Users/sujunchen/Document/labs/2020/scRNA/crc');
source('/Users/sujunchen/Dropbox/me/singleCell/code/enrich_gprofiler.R');
getGEOSuppFiles(GEO, makeDirectory = TRUE, baseDir = getwd(),
  fetch_files = TRUE, filter_regex = NULL);
samp <- read.table('rnaseq_sample.txt');
for(s in samp$V1){
	a <- getGEOSuppFiles(s, fetch_files = TRUE)
};
flist <- list.files('./', 'gz', full.names = TRUE, recursive = TRUE);
mdat <- do.call(cbind, sapply(flist, function(x) read.table(gzfile(x))[, 2, drop = FALSE]));
tmp <- read.table(gzfile(flist[1]));
rownames(mdat) <- tmp$V1;
saveRDS(mdat, file = paste0(Sys.Date(), 'GSE97693_CRC.rds'));
####
colnames(mdat) <- gsub('.*Rna_', '', colnames(mdat))
dat.t <- mdat[, grep('TPM', colnames(mdat))];
colnames(dat.t) <- gsub('\\..*', '', colnames(dat.t));
ref.t <- data.frame(name = colnames(dat.t), samp = gsub('_.*', '', colnames(dat.t)), 
	type = gsub('_.*', '', gsub('CRC[0-9]+_', '', colnames(dat.t))));
dat.t['YTHDF1', dat.t['YTHDF1', ]>0]
###
dat.f <- mdat[, grep('FPKM', colnames(mdat))];
colnames(dat.f) <- gsub('\\..*', '', colnames(dat.f));
ref.f <- data.frame(name = colnames(dat.f), samp = gsub('_.*', '', colnames(dat.f)), 
	type = gsub('_.*', '', gsub('CRC[0-9]+_', '', colnames(dat.f))));

dat.f <- cbind(dat.f[, dat.f['YTHDF1', ]>0.45], dat.f[, dat.f['YTHDF1', ]<0.45]);
groups <- rep(c('high', 'low'), times = c(ncol(dat.f[, dat.f['YTHDF1', ]>0.45]), ncol(dat.f[, dat.f['YTHDF1', ]<0.45])));
pval <- data.frame(pval = apply(dat.f, 1, function(x) wilcox.test(x~groups)$p.value));
pval$lfc <- apply(dat.f+1, 1, function(x) log2(mean(x[groups=='high'])/mean(x[groups=='low'])));
pval$fdr <- p.adjust(pval$pval);
pval <- na.omit(pval);
saveRDS(pval, file = paste0(Sys.Date(), 'diff_fpkm.rds'));
###
up <- rownames(pval[pval$fdr<0.05&pval$lfc>0.58, ]);
ego.up <- enrich_gprofiler(up, rownames(pval));
#save(pval, ego.up, file = paste0(Sys.Date(), 'diff_enrich_fpkm_fc1.5.rda'));
###
up <- rownames(pval[pval$fdr<0.05&pval$lfc>1, ]);
ego.up2 <- enrich_gprofiler(up, rownames(pval));

up <- rownames(pval[pval$fdr<0.05&pval$lfc>0, ]);
ego.up0 <- enrich_gprofiler(up, rownames(pval));
save(ego.up, ego.up2, ego.up0, file = paste0(Sys.Date(), '_diff_enrich_fpkm_0.45.rda'))
###
up <- rownames(pval[pval$fdr<0.05&pval$lfc>1, ]);
to.plot <- dat.f[up, ];
to.plot <- data.frame(t(scale(t(to.plot))));
dbreaks <- c(seq(-18, -1), seq(-0.6, 0.6, 0.05), seq(1, 18));
annotation_col <- data.frame(Y1 = groups);
rownames(annotation_col) <- colnames(to.plot);
pdf(paste0(Sys.Date(), '_diff_fc2.pdf'), height = 5)
pheatmap(to.plot, breaks = dbreaks, cluster_cols = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(length(dbreaks)),
	annotation_col = annotation_col, show_colnames = FALSE, fontsize_row = 2.5);
dev.off();
#### high vs non
dat.f <- cbind(dat.f[, dat.f['YTHDF1', ]>20], dat.f[, dat.f['YTHDF1', ]==0]);
groups <- rep(c('high', 'low'), times = c(ncol(dat.f[, dat.f['YTHDF1', ]>20]), ncol(dat.f[, dat.f['YTHDF1', ]==0])));
pval <- data.frame(pval = apply(dat.f, 1, function(x) wilcox.test(x~groups)$p.value));
pval$lfc <- apply(dat.f+1, 1, function(x) log2(mean(x[groups=='high'])/mean(x[groups=='low'])));
pval$fdr <- p.adjust(pval$pval);
pval <- na.omit(pval);
saveRDS(pval, file = paste0(Sys.Date(), 'diff_fpkm_20.rds'));

up <- rownames(pval[pval$fdr<0.05&pval$lfc>0.58, ]);
ego.up <- enrich_gprofiler(up, rownames(pval));
#save(pval, ego.up, file = paste0(Sys.Date(), 'diff_enrich_fpkm_fc1.5.rda'));
###
up <- rownames(pval[pval$fdr<0.05&pval$lfc>1, ]);
ego.up2 <- enrich_gprofiler(up, rownames(pval));
save(ego.up, ego.up2, file = paste0(Sys.Date(), '_diff_enrich_fpkm_20.rda'));

up <- rownames(pval[pval$fdr<0.05&pval$lfc>1, ]);
to.plot <- dat.f[up, ];
to.plot <- data.frame(t(scale(t(to.plot))));
dbreaks <- c(seq(-18, -1), seq(-0.6, 0.6, 0.05), seq(1, 18));
annotation_col <- data.frame(Y1 = groups);
rownames(annotation_col) <- colnames(to.plot);
pdf(paste0(Sys.Date(), '_diff_fc2_fpkm_20.pdf'), height = 5)
pheatmap(to.plot, breaks = dbreaks, cluster_cols = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(length(dbreaks)),
	annotation_col = annotation_col, show_colnames = FALSE, show_rownames = FALSE);
dev.off();
#### high 103 vs non
dat.f <- cbind(dat.f[, dat.f['YTHDF1', ]>103], dat.f[, dat.f['YTHDF1', ]==0]);
groups <- rep(c('high', 'low'), times = c(ncol(dat.f[, dat.f['YTHDF1', ]>103]), ncol(dat.f[, dat.f['YTHDF1', ]==0])));
pval <- data.frame(pval = apply(dat.f, 1, function(x) wilcox.test(x~groups)$p.value));
pval$lfc <- apply(dat.f+1, 1, function(x) log2(mean(x[groups=='high'])/mean(x[groups=='low'])));
pval$fdr <- p.adjust(pval$pval);
pval <- na.omit(pval);
saveRDS(pval, file = paste0(Sys.Date(), 'diff_fpkm_103.rds'));

up <- rownames(pval[pval$fdr<0.05&pval$lfc>0.58, ]);
ego.up <- enrich_gprofiler(up, rownames(pval));
#save(pval, ego.up, file = paste0(Sys.Date(), 'diff_enrich_fpkm_fc1.5.rda'));
###
up <- rownames(pval[pval$fdr<0.05&pval$lfc>0, ]);
ego.up0 <- enrich_gprofiler(up, rownames(pval));
save(ego.up, ego.up0, file = paste0(Sys.Date(), '_diff_enrich_fpkm_103.rda'));

up <- rownames(pval[pval$fdr<0.05&pval$lfc>0, ]);
to.plot <- dat.f[up, ];
to.plot <- data.frame(t(scale(t(to.plot))));
dbreaks <- c(seq(-18, -1), seq(-0.6, 0.6, 0.05), seq(1, 18));
annotation_col <- data.frame(Y1 = groups);
rownames(annotation_col) <- colnames(to.plot);
pdf(paste0(Sys.Date(), '_diff_fc0_fpkm_103.pdf'), height = 5)
pheatmap(to.plot, breaks = dbreaks, cluster_cols = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(length(dbreaks)),
	annotation_col = annotation_col, show_colnames = FALSE, show_rownames = TRUE);
dev.off();
###
dat.f <- mdat[, grep('FPKM', colnames(mdat))];
colnames(dat.f) <- gsub('\\..*', '', colnames(dat.f));
dat.f <- dat.f[, colnames(dat.f)%in%ref.f[grep('PT', ref.f$type), ]$name];
###
dat.f <- cbind(dat.f[, dat.f['YTHDF1', ]>50], dat.f[, dat.f['YTHDF1', ]==0]);
groups <- rep(c('high', 'low'), times = c(ncol(dat.f[, dat.f['YTHDF1', ]>50]), ncol(dat.f[, dat.f['YTHDF1', ]==0])));
pval <- data.frame(pval = apply(dat.f, 1, function(x) wilcox.test(x~groups)$p.value));
pval$lfc <- apply(dat.f+1, 1, function(x) log2(mean(x[groups=='high'])/mean(x[groups=='low'])));
pval$fdr <- p.adjust(pval$pval);
pval <- na.omit(pval);
saveRDS(pval, file = paste0(Sys.Date(), 'diff_fpkm_50_PT.rds'));

up <- rownames(pval[pval$pval<0.001&pval$lfc>1, ]);
ego.up2 <- enrich_gprofiler(up, rownames(pval));
save(ego.up2, file = paste0(Sys.Date(), '_diff_enrich_fpkm_50_PT.rda'));

up <- rownames(pval[pval$pval<0.001&pval$lfc>1, ]);
to.plot <- dat.f[up, ];
to.plot <- data.frame(t(scale(t(to.plot))));
dbreaks <- c(seq(-18, -1), seq(-0.6, 0.6, 0.05), seq(1, 18));
annotation_col <- data.frame(Y1 = groups);
rownames(annotation_col) <- colnames(to.plot);
pdf(paste0(Sys.Date(), '_diff_fc2_fpkm_50_PT.pdf'), height = 5)
pheatmap(to.plot, breaks = dbreaks, cluster_cols = FALSE, color = colorRampPalette(c('blue', 'white', 'red'))(length(dbreaks)),
	annotation_col = annotation_col, show_colnames = FALSE, show_rownames = FALSE);
dev.off();
