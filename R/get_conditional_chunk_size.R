get_conditional_chunk_size <- function(local_list) {
  chunk_size <- local_list[["chunk_size"]]
  if (is.null(chunk_size)) chunk_size <- 1000L
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size)) {
    stop("local$chunk_size must be a positive integer.", call. = FALSE)
  }
  chunk_size
}
