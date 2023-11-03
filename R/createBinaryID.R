createBinaryID <- function(ssn, overwrite) {
  ## If binaryID.db exists
  if (file.exists("binaryID.db") == TRUE) {
    if (overwrite == TRUE) {
      unlink("binaryID.db")
    } else {
      message("binaryID.db already exists - no changes were made to binaryID.db table\n")
      mm <- TRUE
    }
  }

  ## If binaryID.db does not exist
  options(show.error.messages = FALSE)
  m <- try(if (file.exists("binaryID.db") == FALSE) {
    mm <- FALSE
    ## Define database driver
    driver <- RSQLite::SQLite()

    ## Define connection
    db.name <- "binaryID.db"
    connect <- dbConnect(SQLite(), db.name)

    ## get number of networks from observed sites attribute table...
    net.no <- unique(ssn$edges$netID)

    ## read data into SQLite directly from file
    for (i in seq_len(length(net.no))) {
      network <- paste("net", net.no[i], sep = "")
      file.name <- paste("netID", net.no[i], ".dat", sep = "")

      if (dbExistsTable(connect, network)) {
        dbRemoveTable(connect, network)
      }

      dbWriteTable(connect, network, read.table(
        file = file.name, header = TRUE, sep = ",",
        colClasses = c("numeric", "character")
      ), overwrite = TRUE, row.names = FALSE)

      ## Check to ensure binary files were imported to SQLite database
      if (i == length(net.no)) {
        if (length(dbListTables(connect)) != length(net.no)) {
          dbDisconnect(connect)
          stop("ERROR: binary tables did not import to SQLite database properly")
        }
      }
    }

    ## close the connection and driver
    dbDisconnect(connect)
  }, silent = TRUE)
  options(show.error.messages = TRUE)

  if (mm != TRUE) {
    if (m != TRUE) {
      dbDisconnect(connect)
      stop("ERROR: binary tables did not import to SQLite database properly")
    }
  }
}
