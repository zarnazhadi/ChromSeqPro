# This script loads the input data for the GBMPROM analysis and performs an 
# initial analysis
# Alastair Droop, 2021-04-16

# Load the necessary libraries:
library(GenomicRanges)
library(jsonlite)
library(ggplot2)

# A function to write a given message to stderr:
log.message <- function(..., verbose=NA){
    if(is.na(verbose)){
        verb <- get0('.verbose', ifnotfound=TRUE)
    } else {
        verb <- verbose
    }
    if(identical(verb, TRUE)){
        message(sprintf(...))
        flush(stderr())
    }
}

# A function to quietly load a vector of libraries from character strings:
loadLibrary <- function(x, verbose=NA){
    log.message('loading library "%s"', x, verbose=verbose)
    if(identical(interactive(), TRUE)){
        res <- require(x, character.only=TRUE, quietly=TRUE)
    } else {
        res <- suppressWarnings(suppressPackageStartupMessages(require(x, character.only=TRUE, quietly=TRUE)))
    }
    if(!identical(res, TRUE)) error('failed to load package "%s"', x)
}

# Set the analysis arguments:
if(!identical(interactive(), TRUE)){
    loadLibrary('argparse', verbose=FALSE)
    parser <- ArgumentParser(description='Perform initial analysis of ChIP/CUT&RUN-seq input data')
    parser$add_argument('-c', '--comparisons', dest='comparisons', default=NULL, help='provide analysis.json file')
    parser$add_argument('-i', '--input', dest='input.dir', default=NULL, help='provide input directory')
    parser$add_argument('-e', '--expression', dest='expression', default=NULL, help='promoter-expression.txt')
    parser$add_argument('-o', '--output', dest='result.dir', default=NULL, help='provide results directory')
    parser$add_argument('-g', '--genome', dest='genome', default=NULL, help='provide genome')
    parser$add_argument('-t', '--threshold', dest='threshold', default=NULL, help='threshold')
    args <- parser$parse_args()    
} else {
args <- list(
    'comparisons' = file.path('.', 'analysis.json'),
    'input.dir' = file.path('..', '..', 'results'),
    'expression' = file.path('..', '..', '..', '..', 'input', 'expression', 'promoter-expression.txt'),
    'genome' = file.path('..', '..', '..', '..', 'metadata', 'genome', 'GRCh38-genome.txt'),
    'result.dir' = file.path('..', 'results'),
    'threshold' = 0.00001
)
}

# Load the genome:
message(sprintf('loading genome data from "%s"...', args$genome))
genome <- read.table(args$genome, sep='\t', colClasses=c('character', 'numeric', 'logical', 'character'), col.names=c('chr', 'length', 'isCircular', 'genome'))
genome <- Seqinfo(genome$chr, seqlengths=genome$length, isCircular=genome$isCircular, genome=genome$genome)

# Load the expression data:
message(sprintf('loading expression data from "%s"...', args$expression))
expression <- read.table(args$expression, sep='\t', header=TRUE, col.names=c('chr', 'pos', 'logFC', 'P', 'R'))
expression$start <- expression$pos
expression$end = expression$start
expression$pos <- NULL
expression <- expression[expression$chr %in% seqlevels(genome), ]
expression <- makeGRangesFromDataFrame(expression, keep.extra.columns=TRUE, seqinfo=genome)
expression <- sort(expression)

# Load the sample file:
message(sprintf('loading sample data from "%s"...', args$comparisons))
samples <- jsonlite::fromJSON(args$comparisons)
sample.labels <- names(samples)
names(sample.labels) <- sapply(samples, '[[', 'label')
samples <- as.data.frame(sapply(sort(setdiff(unique(unlist(lapply(samples, names))), 'label')), function(ctype){
    sapply(samples, '[[', ctype)
}))
rownames(samples) <- names(sample.labels)

# Define the valid datasets:
message('building input file list...')
valid.datasets <- list()
for(sample.id in rownames(samples)){
    sample.label <- sample.labels[sample.id]
    for(ctype in colnames(samples)[which(samples[sample.id,]==TRUE)]){
        dataset.label <- sprintf('%s%s', sample.id, ctype)
        dataset.filename <- file.path(args$input.dir, sprintf('%s_%s', sample.label, ctype), sprintf('%s_%s-signal.txt', sample.label, ctype))
        if(!file.exists(dataset.filename)){
            stop(sprintf('%s signal file "%s" missing', dataset.label, dataset.filename))
        }
        valid.datasets[[dataset.label]] <- dataset.filename
    }
}

# Load the raw data:
message('loading target data...')
d.raw <- sapply(valid.datasets, function(filename){
    output <- read.table(filename, sep='\t', header=TRUE, colClasses=c('character', rep('numeric', 10), 'factor'))
    return(makeGRangesFromDataFrame(output, keep.extra.columns=TRUE, seqinfo=genome))
}, simplify=FALSE)

# Build the data:
d <- GRanges(
    seqnames = seqnames(d.raw[[1]]),
    ranges = ranges(d.raw[[1]]),
    seqinfo=genome
)

# Merge with expression data:
message('calculating expression data overlaps...')
o <- findOverlaps(query=resize(d, width=1, fix='center'), subject=expression)
d$expression_logFC <- NA
mcols(d)[queryHits(o), 'expression_logFC'] <- expression[subjectHits(o)]$logFC
d$expression_P <- NA
mcols(d)[queryHits(o), 'expression_P'] <- expression[subjectHits(o)]$P
d$expression_R <- NA
mcols(d)[queryHits(o), 'expression_R'] <- expression[subjectHits(o)]$R

# Merge the binding data:
message('merging sample data...')
for(dataset in names(d.raw)){
    # Add both value and p-value to output here!
    mcols(d)[[sprintf('%s_lambda_mod', dataset)]] <- d.raw[[dataset]]$lambda_mod
    mcols(d)[[sprintf('%s_pvalue', dataset)]] <- d.raw[[dataset]]$p_value
    mcols(d)[[sprintf('%s_padj', dataset)]] <- d.raw[[dataset]]$p_adj
    mcols(d)[[dataset]] <- as.numeric(d.raw[[dataset]]$p_adj <= args$threshold)
}

# Calculate all possible promoter statuses:
status.label <- paste(rownames(samples), collapse=':')

# Get the ctype status factors:
for(ctype in colnames(samples)){
    status.levels <- expand.grid(lapply(1:nrow(samples), function(i){c('0', '1')}))
    for(i in which(!samples[[ctype]])){
        status.levels[, i] <- rep('0', nrow(status.levels))
    }
    status.levels <- unique(apply(status.levels, 1, paste, collapse=''))
    res <- rep('', length(d))
    for(sample.id in rownames(samples)){
        col <- sprintf('%s%s', sample.id, ctype)
        if(col %in% names(valid.datasets)){
            v <- mcols(d)[, col]
        } else {
            v <- rep('0', length(d))
        }
        res <- sprintf('%s%s', res, as.character(v))
    }
    mcols(d)[[ctype]] <- factor(res, levels=status.levels)
}

# Write the data to the resuts file:
res.filename <- file.path(args$result.dir, 'original-results.txt')
message(sprintf('writing combined data to "%s"...', res.filename))
write.table(d, file=res.filename, sep='\t', quote=FALSE, row.names=FALSE)



# A function to plot p-value histograms:
plot.pHist <- function(x, col, filename){
    g <- ggplot(x, aes_string(col, fill='label'))
    g <- g + geom_histogram(bins=100, show.legend=FALSE)
    g <- g + scale_y_continuous(trans=scales::pseudo_log_trans())
    g <- g + theme_classic()
    g <- g + facet_grid(rows=vars(label))
    g <- g + labs(x='p-value', y='Count')
    ggsave(g, filename=filename, width=12, height=16, units='in')
}

# Extract the p-values for each dataset:
d.pvalues <- data.frame(
    'label' = factor(rep(names(d.raw), each=length(d.raw[[1]]))),
    'pvalue' = unlist(sapply(d.raw, function(i){mcols(i)$p_value}, simplify=FALSE)),
    'padj' = unlist(sapply(d.raw, function(i){mcols(i)$p_adj}, simplify=FALSE))
)

# Plot the p-value histograms:
plot.pHist(d.pvalues, 'pvalue', file.path(args$result.dir, 'pvalue-density.pdf'))
plot.pHist(d.pvalues, 'padj', file.path(args$result.dir, 'padj-density.pdf'))




# Extract the expression by status:
d.expression <- rbind.data.frame(
    data.frame(expression=d$expression_P, status=d$P, ctype=factor('P', levels=c('P', 'R'))),
    data.frame(expression=d$expression_R, status=d$R, ctype=factor('R', levels=c('P', 'R')))
)
d.expression <- d.expression[complete.cases(d.expression),]
d.expression <- d.expression[d.expression$expression > 0,]

# A function to plot the expression data boxplots:
plot.expression <- function(x, filename, notch=TRUE){
    g <- ggplot(x, aes(x=status, y=expression, fill=ctype))
    g <- g + geom_boxplot(outlier.shape = NA, notch=FALSE)
    g <- g + theme_classic()
    g <- g + coord_cartesian(ylim=c(0, quantile(x$expression, 0.9)))
    g <- g + theme(axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_blank())
    g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=0, vjust=0.5, size=10, family='mono'), plot.subtitle=element_text(colour='grey25'), legend.title=element_blank())
    g <- g + labs(y='Expression')
    ggsave(g, filename=filename, width=16, height=8, units='in')
}

# Plot the expression data:
plot.expression(d.expression[d.expression$ctype=='R',], file.path(args$result.dir, 'expression-R.pdf'), notch=TRUE)
plot.expression(d.expression, file.path(args$result.dir, 'expression-PR.pdf'), notch=TRUE)




# A function to plot logFC by status change:
plot.dStatus.logFC <- function(x, filename){
    g <- ggplot(x, aes(x=status, y=logFC))
    g <- g + geom_boxplot(outlier.shape = NA, fill='grey95')
    g <- g + theme_classic()
    # g <- g + coord_cartesian(ylim=c(0, quantile(x$logFC, 0.9)))
    g <- g + theme(axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_blank())
    g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, size=4, family='mono'), plot.subtitle=element_text(colour='grey25'), legend.title=element_blank())
    g <- g + labs(x=sprintf('Change in Promoter Status P -> R', status.label), y='LogFC')
    ggsave(g, filename=filename, width=16, height=8, units='in')
}

# # Plot the expression by promoter status change:
d.status <- data.frame(
    'status' = factor(sprintf('%s->%s', d$P, d$R)),
    'logFC' = d$expression_logFC
)
d.status <- d.status[complete.cases(d.status), ]
plot.dStatus.logFC(d.status, file.path(args$result.dir, 'logFC-P_R.pdf'))
