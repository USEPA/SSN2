# Copy the mf04p .ssn data to a local directory and read it into R
# When modeling with your .ssn object, you will load it using the relevant
# path to the .ssn data on your machine
copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(
  temp_path,
  predpts = c("pred1km", "CapeHorn"),
  overwrite = TRUE
)
ssn_create_distmat(
  mf04p,
  predpts = c("pred1km", "CapeHorn"),
  overwrite = TRUE,
  among_predpts = TRUE
)
