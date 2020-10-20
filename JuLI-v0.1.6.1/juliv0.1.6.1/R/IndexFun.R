#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()


#' Indexes the given Case and Control BAM file if not already indexed.
#'
#' @param CaseBam Path to Case BAM file
#' @param ControlBam Path to Control BAM file (Can be NULL)
#' @returns NULL
IndexFun <- function(CaseBam, ControlBam) {
  # If Case BAM file not indexed, index the BAM file
  caseBamIndexMissing <- Rsamtools::BamFile(CaseBam)$index %>% isEmpty()
  if (caseBamIndexMissing) {
    indexBam(CaseBam) %>% invisible() # Index the BAM file
  }
  
  # If Control BAM file exists but not indexed, index the Control BAM file
  ControlBamExists <- !isEmpty(ControlBam)
  if (ControlBamExists) {
    ControlBamIndexMissing <- BamFile(ControlBam)$index %>% isEmpty()
    if (ControlBamIndexMissing) {
      indexBam(ControlBam) %>% invisible()
    }
  }
}

