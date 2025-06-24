#' Read and Convert Sparse Matrix to Full Data Table
#'
#' This function reads a sparse Matrix Market file along with associated
#' feature and barcode files, then converts the matrix to a full matrix
#' and returns it as a `data.table`.
#'
#' @param matrix_file Path to the `.mtx` or `.mtx.gz` sparse matrix file.
#' @param features_file Path to the `.tsv` or `.tsv.gz` features file.
#' @param barcodes_file Path to the `.tsv` or `.tsv.gz` barcodes file.
#'
#' @return A `data.table` with genes as row names and cells as columns.
#' @export
#'
#' @importFrom Matrix readMM
#' @importFrom data.table fread as.data.table
#'
#' @examples
#' \dontrun{
#' dt <- read_full_matrix("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
#' }
read_full_matrix <- function(matrix_file, features_file, barcodes_file) {
  # Read matrix, features, and barcodes
  sparse_mat <- Matrix::readMM(matrix_file)
  features <- data.table::fread(features_file, header = FALSE)
  barcodes <- data.table::fread(barcodes_file, header = FALSE)

  # Assign row and column names
  rownames(sparse_mat) <- features$V1
  colnames(sparse_mat) <- barcodes$V1

  # Convert to dense matrix, then to data.table
  full_mat <- as.matrix(sparse_mat)
  dt_matrix <- data.table::as.data.table(full_mat, keep.rownames = "Gene")

  return(dt_matrix)
}

