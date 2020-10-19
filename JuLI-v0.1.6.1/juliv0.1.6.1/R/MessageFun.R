#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

MessageFun=function(step){
  Messages=c('A path of CaseBam file is required.',
             'Paths of reference files (Refgene and Gap) are required.',
             'Row number of VCF is zero.',
             'Paths of reference files (Refgene,Cosmic,Pfam,and Uniprot) are required.',
             'A path of vcf file is required.',
             'A reference fasta is required.',
             'Paths of CaseBam files are required.')
  writeLines(paste0(paste(rep('#',75),collapse=""),'\n',Messages[step],'\n',paste(paste(rep('#',75),collapse=""))))
}