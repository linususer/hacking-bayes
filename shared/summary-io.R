# Writing and reading the tab-separated summaries in <project>/summaries/.

# Write `summary_df` to `summary_path` (nothing happens if it is NULL) and report it.
write_summary <- function(summary_df, summary_path, label) {
  if (is.null(summary_path)) return(invisible(summary_df))
  dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
  write.table(summary_df, summary_path, sep = "\t", quote = TRUE, row.names = FALSE)
  cat("=> Saved", label, "summary to:", summary_path, "\n")
  invisible(summary_df)
}

# A summary is either a data.frame (or a list of columns) or the path to a summary file.
read_summary <- function(summary) {
  if (is.character(summary)) {
    read.table(summary, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(summary)
  }
}
