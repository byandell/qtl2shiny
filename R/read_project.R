#' @export
#' @importFrom tools file_ext file_path_sans_ext
read_project <- function(project_info, dataname, columns, filetype) {
  # Read data frame or matrix in some file format.
  
  if(!nrow(project_info))
    return(NULL)
  
  # Taxa and project paths.
  taxapath <- file.path(project_info$directory,
                        project_info$taxa)
  projectpath <- file.path(taxapath,
                           project_info$project)
  
  # Compare file roots in project path to dataname.
  match_filename <- function(dataname, filepath) {
    filenames <- list.files(filepath)
    m <- grep(dataname, filenames)
    if(!length(m)) {
      fileroots <- tools::file_path_sans_ext(filenames)
      m <- grep(tools::file_path_sans_ext(dataname), fileroots)
      if(length(m) > 1) {
        cat(paste("multiple", dataname, "matches"),
            file = stderr())
      }
    }
    if(length(m))
      file.path(filepath, filenames[m])
    else
      NULL
  }
  
  datapath <- match_filename(dataname, projectpath)
  if(is.null(datapath))
    datapath <- match_filename(dataname, taxapath)
  if(is.null(datapath))
    return(NULL)
  
  # File type in order of preference. First use filetype if supplied.
  # Then use extensions of datapath.
  # Watch out for CSV, as may need to preserve column characteristics.
  filetypes <- c("fst","feather","rds","csv")
  datatypes <- tools::file_ext(datapath)
  if(missing(filetype)) {
    # pick in order of filetypes
    m <- match(filetypes, datatypes, nomatch = 0)
    if(!any(m > 0)) {
      cat(paste("file type", paste(datatypes, collapse = ","), "not in approved list"),
          file = stderr())
      return(NULL)
    }
    m <- m[m>0][1]
    filetype <- datatypes[m]
  } else {
    filetype <- match.arg(filetype, filetypes)
    if(!(filetype %in% datatypes)) {
      cat(paste("file type", filetype, "not in directory"),
          file = stderr())
      return(NULL)
    }
  }
  
  # If more that one in datapath, pick by order of filetypes.
  m <- match(filetype, datatypes)
  datapath <- datapath[m]
  
  out <- switch(filetype,
         fst     = fst::read_fst(datapath, columns),
         feather = feather::read_feather(datapath, columns),
         rds     = readRDS(datapath),
         csv     = read.csv(datapath, stringsAsFactors = FALSE))
  
  if(filetype %in% c("rds","csv")) {
    # Pick columns post hoc.
    if(is.data.frame(out) | is.matrix(out))
      out <- out[, columns, drop = FALSE]
  }
  out
}
