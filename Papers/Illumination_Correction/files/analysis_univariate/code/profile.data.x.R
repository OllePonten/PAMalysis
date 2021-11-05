

factors               <- function(obj, ...) UseMethod("factors", obj)
feats                 <- function(obj, ...) UseMethod("feats", obj)
prune.feats           <- function(obj, ...) UseMethod("prune.feats", obj)
rename.feats          <- function(obj, ...) UseMethod("rename.feats", obj)

factors.profile.data.x    <- function(obj, ...) obj$data[obj$factor_cols]
feats.profile.data.x      <- function(obj, ...) obj$data[obj$feat_cols]


profile.data.x <- function(cf) {
  'Construct a profile.data.x object by reading a yml config file'
  
  # ------------ Setup --------------
  
  # Load the config file
  cfg <- yaml.load_file(cf)
  cfg$cwd <- dirname(cf)
  if (!is.null(cfg$profile_dir)) 
    cfg$cwd <- file.path(cfg$cwd, cfg$profile_dir)
  
  
  # Get name of data file (both csv and rda)
  cf_name <- sapply(strsplit(basename(cf),"\\."), 
                    function(x) paste(x[1:(length(x)-1)], collapse=".")) 
  if(is.null(cfg$profile_file))
    cfg$profile_file <- paste(cf_name, "csv", sep=".")
  
  if(is.null(cfg$profile_file_binary))
    cfg$profile_file_binary <- paste(cf_name, "rda", sep=".")
  
  frda <- file.path(cfg$cwd, cfg$profile_file_binary)
  fcsv <- file.path(cfg$cwd, cfg$profile_file)
  
  
  # Read binary file, or if it doesn't exist, then csv. 
  # Save binary file if it doesn't exist
  
  if (file.exists(frda)) {
    data <- readRDS(frda)
  } else {
    data <- data.frame(read.csv(fcsv, header=TRUE))
    saveRDS(data, file=frda)
  }
  
  # Create profile.data.x object
  obj <- list(data=data, cfg=cfg)
  class(obj) <- "profile.data.x"
  
  
  # ------------ Get features --------------
  allfeat_cols <- c(obj$cfg$feat_start:length(names(obj$data)))
  exclude_cols <- NULL
  if (!is.null(obj$cfg$exclude_cols)) {
    for (f in obj$cfg$exclude_cols)
      exclude_cols <- c(exclude_cols, grep(f, names(obj$data)))  
  }
  obj$feat_cols <- names(obj$data)[setdiff(allfeat_cols, exclude_cols)]
  
  
  obj$cfg$treatment_tag <-  "Treatment_"
  
  
  #------- Rename meta-data column names
  dbname <- obj$cfg$mapping$dbname
  metadata_tag <- obj$cfg$mapping$metadata_tag
  
  #  browser() 
  
  metadata_names <- gsub(paste(dbname, '.', sep=''), '', 
                         gsub(metadata_tag, '',
                              names(obj$data)[grep(metadata_tag,
                                                   names(obj$data))]))
  metadata_names_processed <- c()
  for (s2 in metadata_names) {
    expm <- obj$cfg$mapping$explicit_mappings
    s1 <- if (!is.null(expm[[s2]])) expm[[s2]] else s2
    s3 <- paste(dbname, paste(metadata_tag, s2, sep=''), sep=".")
    names(obj$data)[names(obj$data)==s3] <- s1
    s4 <- s1
    
    lastchar <- str_sub(s1,-1,-1)
    # convert to factor unless there is * or | in the name
    cast_as_factor <- (!is.null(obj$cfg$mapping$cast_as_factor)
                       && obj$cfg$mapping$cast_as_factor)
    if (!(lastchar %in% c('*', '|')) && cast_as_factor) {
      obj$data[,s1] <- factor(obj$data[,s1])
    }
    
    # convert to numeric if there is | in the name
    if (lastchar == '|') {
      obj$data[,s1] <- as.numeric(obj$data[,s1])  
      # Trim the name
      s4 <- str_sub(s1, 1,nchar(s1)-1)
      names(obj$data)[names(obj$data)==s1] <- s4	
    }
    
    # convert to character if there is * in the name
    if (lastchar == '*') {
      obj$data[,s1] <- as.character(obj$data[,s1])				
      # Trim the name
      s4 <- str_sub(s1, 1,nchar(s1)-1)
      names(obj$data)[names(obj$data)==s1] <- s4	
    }
    
    # Append X to header names starting with a digit
    if (str_sub(s4, 1,1) %in% as.character(seq(10)-1)) {      
      names(obj$data)[names(obj$data)==s4] <- paste('X', s4, sep='')
      s4 <- paste('X', s4, sep='')
    }
    
    metadata_names_processed <- c(metadata_names_processed, s4)
    
  }
  
  obj$factor_cols <- metadata_names_processed
  
  #------- Format plate names
  if (!is.null(obj$cfg$plate_names)) {
    obj$data$Plate                <- factor(obj$data$Plate,
                                            levels=obj$cfg$plate_names)
    obj$data$Plate_               <- obj$data$Plate
    obj$factor_cols               <- c(obj$factor_cols, 'Plate_')
    levels(obj$data$Plate)        <- obj$cfg$plate_names_abbrev
  } else {
    obj$data$Plate_               <- obj$data$Plate    
  }
  
  #------- Remove No-treatment
  if (!is.null(obj$cfg$exclude_nt) && obj$cfg$exclude_nt) {
    obj$data <- obj$data[obj$data$Treatment != obj$cfg$no_treat_name,]      
    obj$data <- droplevels(obj$data)
  }
  
  
  #------- Create a treatment column by concatenating two columns
  if (!is.null(obj$cfg$create_treatment_by_combining)) {
    cols <- str_split(obj$cfg$create_treatment_by_combining, ',')[[1]]
    obj$data$Treatment <- paste(obj$data[,cols[1]], obj$data[,cols[2]], sep='_')
    obj$factor_cols <- c(obj$factor_cols, 'Treatment')
    
  }
  
  
  #------- Format treatment names
  if (!is.null(obj$cfg$auto_treatment_synonym)) {
    syn <- obj$cfg$auto_treatment_synonym
    # TODO: Handle the case that Type does not exist
    df <- unique(obj$data[,c(syn, 'Treatment', 'Type')])
    names(df)[names(df)==syn] <- 'syn'
    df <- df[with(df, order(syn, Treatment)), ]
    df0 <- ddply(df, .(syn), index_factor_levels)
    df <- df0
    df <- df[with(df, order(Type, syn, Treatment)), ]	  
    names(df)[names(df)=='Treatment'] <- 'name_Treatment'
    names(df)[names(df)=='abbrev'] <- 'name_TreatmentAbbrev'
    names(df)[names(df)=='syn'] <- sprintf('name_%s',
                                           obj$cfg$auto_treatment_synonym)
    df$Type <- NULL
    obj$cfg$tab_treatment_synonym <- df	  
  }
  
  if (!is.null(obj$cfg$treatment_synonyms)) {
    df <- convert_str_to_dataframe(obj$cfg$treatment_synonym)
    obj$cfg$tab_treatment_synonym <- df
  }
  
  if (!is.null(obj$cfg$auto_treatment_synonym) ||
        !is.null(obj$cfg$treatment_synonyms)) {
    obj$data$Treatment <- factor(obj$data$Treatment, levels = df$name_Treatment)
    df <- obj$cfg$tab_treatment_synonym
    sfnames <- names(df)[!(names(df) %in% c('name_Treatment'))]
    for (sfname in sfnames) 
      obj$data[,gsub('name_', '', sfname)] <-
      create_synonym_factor(obj$data$Treatment, df[,sfname])
    
    # do this because the syn is already present in the table
    sfnames <- setdiff(sfnames, c(sprintf('name_%s',
                                          obj$cfg$auto_treatment_synonym)))
    
    obj$factor_cols <- c(obj$factor_cols, gsub('name_', '', sfnames))
    
  }
  
  obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)]
  obj
}

prune.feats.profile.data.x <- function(obj, feat_cols, keep = FALSE, ...) { 
  'Remove features'
  if (!keep)
    obj$feat_cols <- setdiff(obj$feat_cols, feat_cols)
  else
    obj$feat_cols <- feat_cols
  obj$data <- obj$data[,c(obj$factor_cols, obj$feat_cols)] 
  obj 
}

rename.feats.profile.data.x <- function(obj, new_names,  ...) { 
  'Rename features'
  
  X <- feats(obj)
  names(X) <- new_names
  obj$feat_cols <- names(X)
  obj$data <- cbind(X, factors(obj))
  obj 
  
}