for (pkg in c('Rsamtools', 'magrittr', 'foreach', 'data.table')) {
    suppressMessages(library(pkg, character.only = TRUE))
}


allowed.chromosome.names <- c(
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "X",
    "Y"
)


# TODO: Add documentation
measureBAMStats <- function(bam.file) {
    bam <- Rsamtools::BamFile(bam.file)
    bam.info <- Rsamtools::seqinfo(bam)
    sequence.length = GenomeInfoDb::seqlengths(bam.info)
    
    # Filtered chromosomes
    chrom.list <- GenomeInfoDb::seqnames(bam.info) %>%
                .[. %in% allowed.chromosome.names]
    
    read.length <- calcMedianReadLength(bam.file)

    output <- foreach( # Run extractChromStats on all chromosomes in parallel
        c = 1:length(chrom.list),
        .options.multicore = list(preschedule = FALSE)
        ) %dopar% {
            selected.chrom = chrom.list[c]
            selected.chrom.length = sequence.length[selected.chrom]
            return(extractChromStats(bam.file, selected.chrom, selected.chrom.length))
        }
    
    MedianInsertSize    <- sapply(output, function(x){x[[1]]}) %>% unlist() %>% median(., na.rm=T)
    TotalReadCount      <- sapply(output, function(x){x[[2]]}) %>% sum()
    TotalSplitReadCount <- sapply(output, function(x){x[[3]]}) %>% sum()

    ReferenceLength     <- sequence.length[names(sequence.length) %in% chrom.list] %>% sum()

    Stats <- data.table(
        Chromosome = paste(chrom.list, collapse=';'),
        ReferenceLength,
        TotalReadNumber = TotalReadCount,
        SplitReadNumber = TotalSplitReadCount,
        MedianInsertSize,
        ReadLength = read.length
    )
    return(Stats)
}




#' Reads a BAM file and calculates the median read length
#' 
#' @return integer Median read length value
#' @example
#' calcMedianReadLength('/path/to/sample.bam')
#' 100
calcMedianReadLength <- function(bam.file) {
    scanBam(
        file = BamFile(bam.file, yieldSize=1000), # TODO: Why only read 1000 reads when there are much more?
        param = ScanBamParam(what='seq')
    ) %>%
    data.frame() %>%
    .[,1] %>%
    nchar() %>%
    median()
}

# TODO: Finish documentation
#' Extracts insert_size, number of reads, number of splitting reads
#' from a single chromosome in a BAM file. 
#' @param bam.file Path to BAM file
#' @param chrom string of chromosome name
#' @param chrom.length integer of chromosome length
extractChromStats <- function(bam.file, chrom, chrom.length, mapqFilter = 0) {
    # TODO: Remove the chrom.length parameter by implicitly finding out the end position of the chromosome!
    scan.range = GRanges(
        seqnames = chrom,
        ranges = IRanges(1, chrom.length)
    )

    scan.parameters <- ScanBamParam(
        which = scan.range,
        what = c('cigar', 'isize'),
        mapqFilter = mapqFilter
    )
    scan.result <- scanBam(bam.file, param = scan.parameters) %>% data.table::rbindlist()
    
    isize <- abs(scan.result$isize)
    ReadCount <- nrow(scan.result)
    # Splitting reads are reads with at least 1 soft-clipped base
    splitReadCount <- scan.result$cigar %>% .[grepl('S', .)] %>% length()
    
    list(isize, ReadCount, splitReadCount)
}
