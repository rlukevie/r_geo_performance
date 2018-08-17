#library(sf)
#library(rgdal)
# library(microbenchmark)
#library(ggplot2)
library(benchmarkme)
library(readr)
library(rlist)
library(stringr)




# test_performance <- function(FUN, ..., fun_string = "", times = 1, testtype, package, path) {

test_performance <- function(FUN, parstr, fun_string, times = 1, remove_output_files = FALSE, testtype, package, path) {
  # get system information
  sysname <- Sys.info()[["sysname"]]
  
  # remove all output files
  if (remove_output_files) {
    do.call(file.remove, list(list.files("data_output", full.names = TRUE)))
  }
  
  memory_diffs <- vector(mode = "numeric", length = times)
  durations <- vector(mode = "numeric", length = times)
  for (i in 1:times) {
    gc(reset = TRUE) 
    Sys.sleep(0.5)
    if (sysname == "Linux") {
      gc_res <- gc()
      memory_before <- gc_res[11] + gc_res[12]
    } 
    if (sysname == "Windows") {
      memory_before <- memory.size()
    }
    
    
    # print(parstr)
    # print(fun_string)
    # if (parstr != "NULL = ''") {
    #   start_time <- Sys.time()
    #   eval(parse(text = paste0(fun_string, "(", parstr, ")")))
    # }
    # else {
    #   start_time <- Sys.time()
    #   eval(parse(text = fun_string))
    # }
    # durations[i] <- (Sys.time() - start_time) #* 1000 # in milliseconds
    
    if (parstr != "FALSE = ''") {
      eval(parse(text = paste0("time <- system.time(", fun_string, "(", parstr, ")", ")")))
    } else {
      eval(parse(text = paste0("time <- system.time(", fun_string, ")")))
    }
    
    if (sysname == "Windows") {
      time <- time[["elapsed"]]
    }
    durations[i] <- time
    
    Sys.sleep(0.5)
    if (sysname == "Linux") {
      gc_res <- gc()
      memory_after <- gc_res[11] + gc_res[12]
      memory_diffs[i] <- memory_after - memory_before
    } 
    if (sysname == "Windows") {
      memory_after <- memory.size()
      memory_diffs[i] <- memory_after - memory_before
    }
  }
  memory_diff_mean <- mean(memory_diffs)
  memory_diff_median <- median(memory_diffs)
  execution_time_mean <- mean(durations)
  execution_time_median <- median(durations)
  
  # create result data.frame
  # parameters <- paste(paste0(names(list(...)), " = ", list(...)), collapse = " ; ")
  if ((length(grep("_s", fun_string)) > 0) | (length(grep("_s", parstr)) > 0)) {
    size = "small"
  } else if ((length(grep("_m", fun_string)) > 0) | (length(grep("_m", parstr)) > 0)) {
    size = "medium"
  } else if ((length(grep("_l", fun_string)) > 0) | (length(grep("_l", parstr)) > 0)) {
    size = "large"
  } else {
    size = NA
  }
  
  datatype <- c(0,0,0,0) # point, line, poly, raster
  if ((length(grep("point", fun_string)) > 0) | (length(grep("point", parstr)) > 0)) {
    datatype[1] <- 1
  }
  if ((length(grep("line", fun_string)) > 0) | (length(grep("line", parstr)) > 0)) {
    datatype[2] <- 1
  }
  if ((length(grep("poly", fun_string)) > 0) | (length(grep("poly", parstr)) > 0)) {
    datatype[3] <- 1
  }
  if ((length(grep("raster", fun_string)) > 0) | (length(grep("raster", parstr)) > 0)) {
    datatype[4] <- 1
  }
  # print(as.character(fun_string))
  result <- data.frame(datetime = Sys.time(),
                       testtype = testtype,
                       package = package,
                       FUN = fun_string,
                       size = size,
                       point = datatype[1],
                       line = datatype[2],
                       poly = datatype[3],
                       raster = datatype[4],
                       parameters = parstr,
                       times = times,
                       memory_diff_mean = memory_diff_mean,
                       memory_diff_median = memory_diff_median,
                       execution_time_mean = execution_time_mean,
                       execution_time_median = execution_time_median,
                       # benchmark = summary(benchmark),
                       sysname = sysname,
                       cpu = get_cpu()$model_name,
                       stringsAsFactors = FALSE)
  # result$benchmark.expr <- NULL
  
  # write to csv
  if (file.exists(path)) {
    append = TRUE
  } else {
    append = FALSE
  }
  write_csv(result,
            path = path,
            append = append)
  
  # return
  return(result)
}

test_performance_grid <- function(parameter_grid) {
  gc(reset = TRUE)
  Sys.sleep(1)
  # check if structure of config.yml is correct
  if (parameter_grid$testtype == "") {
    stop("Testtype [testtype] not defined in config.yml. Insert something like 'testtype: \"Read vector data: shapefiles\"'")
  } else {
    testtype <- parameter_grid$testtype
    parameter_grid$testtype <- NULL
  }
  
  if ("times" %in% names(parameter_grid)) {
    times <- parameter_grid$times
    parameter_grid$times <- NULL
  } else {
    stop("Number of performance evaluations [times] not defined in config.yml. Insert something like 'times: 5'")
  }
  
  if ("path" %in% names(parameter_grid)) {
    path <- parameter_grid$path
    parameter_grid$path <- NULL
  } else {
    stop("Path to csv file [path] not defined in config.yml. Insert something like 'path: \"test_results.csv\"'")
  }
  
  # create a list of all possible function-parameter-combinations
  sections <- names(parameter_grid)
  functions <- list()
  for (section in sections) {
    functions <- c(functions, parameter_grid[[section]][["function"]])
  }
  gridlist <- list()
  for (s in sections) {
    gridlist[[s]] <- expand.grid(parameter_grid[[s]][["param"]], stringsAsFactors = FALSE)
  }
  # loop through these combinations
  for (section in sections) {
    for (param_combination in 1:nrow(gridlist[[section]])) {
      
      # read package name
      package <- parameter_grid[[section]]$package
      
      # read function name
      f <- parameter_grid[[section]][["function"]]
      
      # construct parameter string
      if (length(gridlist[[section]][param_combination, ]) == 1) {
        names_list <- list(colnames(gridlist[[section]]))
      } else {
        names_list <- names(gridlist[[section]][param_combination, ])
      }
      # print(names_list)
      # strange bug with parameters called 'y'
      names_list[names_list == "TRUE"] <- "y"
      
      values <- as.list(gridlist[[section]][param_combination, ])
      values <- lapply(values, 
                       function(x) {if ((typeof(x) == "character") & (substring(x, 1, 1) != "!")) {paste0("'", x, "'")} 
                         else if (substring(x, 1, 1) == "!") {substring(x, 2)}
                         else {x}})
      
      parvec <- c()
      for (i in 1:length(names_list)) {
        parvec <- c(parvec, paste0(names_list[i], " = ", values[i]))
      }
      parstr <- paste(parvec, collapse = ", ")
      # print(parstr)
      # print(f)
      # print(parstr)
      # fun_string = paste0(f, "(", parstr, ")")
      fun_string = f
      print(paste0("fun_string: ", fun_string))
      # start performance test and print some results on console
      print(paste0("+++ ", f, "{", package, "}"))
      print(paste0("    Test iterations:     ", times))
      print(paste0("    Parameters:          ", parstr))
      # test_performance_string = paste0("result <- test_performance(", 
      #                                  "FUN = ", f,
      #                                  ", ", parstr,
      #                                  ", fun_string = ", fun_string,
      #                                  ", times = ", times,
      #                                  ", testtype = '", testtype, "'",
      #                                  ", package = '", package, "'",
      #                                  ", path = '", path, "'",
      #                                  ")")
      # print(test_performance_string)
      eval(parse(text = paste0("result <- test_performance(", 
                               "FUN = '", f, "'",
                               ", parstr = \"", parstr, "\"",
                               ", fun_string = \"", fun_string, "\"",
                               ", times = ", times,
                               ", testtype = '", testtype, "'",
                               ", package = '", package, "'",
                               ", path = '", path, "'",
                               ")")))
      # print(paste0("    Function call:       ", fun_string))
      print(paste0("    Mean execution time: ", round(result[1, "execution_time_mean"], digits = 2), " s"))
      print(paste0("    Mean memory used:    ", round(result[1, "memory_diff_mean"], digits = 2), " Mb"))
      print("---")
    }
  }
}
