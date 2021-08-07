# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#'  Diagnostics for messy joins
#'
#'  given a data frame with two ways to group rows,
#'  summarize and give examples of situations where the mapping is not 1-1
#'
#' @param x_cols tidyselect specification of a set of columns defining objects
#' @param y_cols tidyselect specification of a set of columns defining objects
#'
#' data <- data.frame(
#' 	x=c(1,2,NA,3,4,5,6,6,6,7,7),
#' 	y=c("a",NA,"c","d","d","d","e","f","g","h","h"))
#'
#' data %>% summarize_map(
#'   x_cols = x),
#'   y_cols = y))
#' X<-[x]:
#'   |X|: 7                          # number of groups
#'   |is.na.X|: 1                    # number of groups with NA in atleaset 1 col
#'   range(|x|:X): 1, 3              # size range of groups
#' Y<-[y]:
#'   |Y|: 7
#'   |is.na.Y|: 1
#'   range(|y|:Y): 1, 3
#' [X U Y]:                          # grouping by the union of xcols and ycols
#'   |X U Y|: 8
#'   |is.na.XUY|: 2
#'   range(|z|:X U Y): 1, 2
#' [X @ Y]:
#'   |X ~ Y|: 5
#'   |X:X < Y|, |Y:Y < X|: 1, 1
#'   |X:X > Y|, |Y:Y < X|: 3, 3
#' $is.na.X
#'    x y
#' 1 NA c
#'
#' $is.na.Y
#'   x    y
#' 1 2 <NA>
#'
#' $dup.XUY
#'   x y
#' 1 7 h
#' 2 7 h
#'
#' $dup.X
#'   x y
#' 1 6 e
#' 2 6 f
#' 3 6 g
#'
#' $dup.Y
#'   y x
#' 1 d 3
#' 2 d 4
#' 3 d 5
#' @export
summarize_map <- function(
	data,
  x_cols,
  y_cols,
  n_examples = 4,
  verbose = FALSE) {

	# convert column selections named vectors of column indices into data
	x_cols <- tidyselect::eval_select(enquo(x_cols), data)
	y_cols <- tidyselect::eval_select(enquo(y_cols), data)
	xUy_cols <- union(x_cols, y_cols)
	names(xUy_cols) <- names(data[xUy_cols])

	if(verbose) {
		cat("The following is a report of the relationship between two different ways of identifying instances\n")
	}

	# example rows
	problems <- list()

	count_xUy <- data %>%
		dplyr::count(dplyr::across(xUy_cols)) %>%
    dplyr::ungroup()
	count_x <- count_xUy %>%
    dplyr::count(dplyr::across(names(x_cols)), name = "size") %>%
    dplyr::ungroup()
	count_y <- count_xUy %>%
		dplyr::count(dplyr::across(names(y_cols)), name = "size") %>%
    dplyr::ungroup()

	if (verbose) {
		cat("\nProperties of X identifiers:\n")
	}
	cat("X<-[", paste(names(x_cols), collapse = ", "), "]:\n", sep = "")
	cat("  |X|: ", count_x %>% na.omit(method = "r") %>% nrow, sep = "")

	na_count <- data %>%
		dplyr::select(x_cols) %>%
		complete.cases() %>%
		magrittr::not() %>%
		sum()
	cat(ifelse(na_count == 0, "", paste0(" (", na_count, " NA)")), "\n", sep = "")

	size_dist <- count_x %>%
		stats::na.omit(method="r") %>%
		dplyr::count(size) %>%
		dplyr::ungroup()
	if (nrow(size_dist) < 12) {
		cat("  count*size: ",
			paste(size_dist$n, size_dist$size, sep = "*", collapse = ", "),
			"\n", sep = "")
	} else {
		top <- 1:6
		bottom <- (nrow(size_dist) - 6+1):nrow(size_dist)
		cat("  count*size: ",
			paste(
				size_dist$n[top],
				size_dist$size[top], sep = "*", collapse = ", "),
			", ... ",
			paste(
				size_dist$n[bottom],
				size_dist$size[bottom], sep = "*", collapse = ", "),
			"\n", sep="")
	}

	if (verbose) {
		cat("\nProperties of the Y identifiers:\n")
	}
	cat("Y<-[", paste(names(y_cols), collapse = ", "), "]:\n", sep = "")
	cat("  |Y|: ", count_y %>% na.omit(method = "r") %>% nrow, sep = "")
	na_count <- data %>%
		dplyr:::select(y_cols) %>%
		complete.cases() %>%
		magrittr::not() %>%
		sum()
	cat(ifelse(na_count == 0, "", paste0(" (", na_count, " NA)")), "\n", sep = "")

	size_dist <- count_y %>%
		na.omit(method = "r") %>%
		dplyr::count(size) %>%
		dplyr::ungroup()
	if (nrow(size_dist) < 12) {
		cat("  count*size: ",
			paste(size_dist$n, size_dist$size, sep = "*", collapse = ", "),
			"\n", sep = "")
	} else {
		top <- 1:6
		bottom <- (nrow(size_dist) - 6+1):nrow(size_dist)
		cat("  count*size: ",
			paste(
				size_dist$n[top],
				size_dist$size[top], sep = "*", collapse = ", "),
			", ... ",
			paste(
				size_dist$n[bottom],
				size_dist$size[bottom], sep = "*", collapse = ", "),
			"\n", sep="")
	}

	if (verbose) {
		cat("\nProperties of the intersection of union of the X and Y identifiers:\n")
	}
	cat("[X U Y]:\n")
	cat("  |X U Y|: ", count_xUy %>% na.omit(method = "r") %>% nrow, sep = "")
	na_count <- data %>%
    dplyr:::select(!!!xUy_cols) %>%
    complete.cases %>%
    magrittr::not() %>%
    sum()
	cat(ifelse(na_count == 0, "", paste0(" (", na_count, " NA)")), "\n", sep = "")

	size_dist <- count_xUy %>%
		na.omit(method = "r") %>%
		dplyr::rename(size = n) %>%
		dplyr::count(size) %>%
		dplyr::ungroup()
	if (nrow(size_dist) < 12) {
		cat("  count*size: ",
			paste(size_dist$n, size_dist$size, sep = "*", collapse = ", "),
			"\n", sep="")
	} else {
		top <- 1:6
		bottom <- (nrow(size_dist) - 6+1):nrow(size_dist)
		cat("  count*size: ",
			paste(
				size_dist$n[top],
				size_dist$size[top], sep = "*", collapse = ", "),
			", ... ",
			paste(
				size_dist$n[bottom],
				size_dist$size[bottom], sep = "*", collapse = ", "),
			"\n", sep = "")
	}


	count_xUy <- count_xUy %>% na.omit(method = "r")

	if (verbose) {
		cat("Properties of the intersection of the X and Y identifiers:\n")
	}
	cat("[X @ Y]:\n")
	if (verbose) {
		cat("  Number of X and Y identifiers that are 1 to 1:\n")
	}
	cat("  |X ~ Y|: ",
		count_xUy %>%
  		dplyr::semi_join(
				count_x %>% filter(size == 1),
        by = names(x_cols)) %>%
			dplyr::semi_join(
        count_y %>% filter(size == 1),
        by = names(y_cols)) %>%
			nrow,
		"\n", sep = "")

	if (verbose) {
		cat("  Number of X and Y identifiers where an X identifier maps to multiple Y identifiers:\n")
	}
	cat(
    "  |X:X < Y|, |Y:X < Y|: ",
		count_xUy %>%
      dplyr::semi_join(
        count_x %>% dplyr::filter(size > 1),
        by = names(x_cols)) %>%
      nrow,
    ", ",
		count_xUy %>%
		  dplyr::count(dplyr::across(names(x_cols)), name = "size") %>%
		  dplyr::filter(size > 1) %>%
      nrow,
		"\n", sep = "")

	if (verbose) {
		cat("  Number of X and Y identifiers where a Y identifier maps to multiple X identifiers:\n")
	}
	cat(
    "  |X:X > Y|, |Y:X > Y|: ",
		count_xUy %>%
      dplyr::semi_join(
        count_y %>%
        dplyr::filter(size > 1),
      by = names(y_cols)) %>%
    nrow,
    ", ",
		count_xUy %>%
      dplyr::count(dplyr::across(names(y_cols)), name = "size") %>%
      dplyr::filter(size > 1) %>%
      nrow,
		"\n", sep="")

	#is.na.X
	ex_rows <- data %>%
    dplyr:::select(x_cols) %>%
    complete.cases() %>%
    magrittr::not() %>%
    which()
	if (length(ex_rows)) {
		if (!is.null(n_examples) && (n_examples < length(ex_rows))) {
			ex_rows <- ex_rows %>% sample(n_examples, replace = FALSE)
		}
		problems$is.na.X <- data %>%
				dplyr::slice(ex_rows) %>%
				dplyr::arrange(dplyr::across(x_cols))
	}

	#is.na.Y
	ex_rows <- data %>%
		dplyr:::select(y_cols) %>%
		complete.cases() %>%
		magrittr::not() %>%
    which()
	if (length(ex_rows)) {
		if (!is.null(n_examples) && (n_examples < length(ex_rows))) {
			ex_rows <- ex_rows %>% sample(n_examples, replace = FALSE)
		}
		problems$is.na.Y <- data %>%
				dplyr::slice(ex_rows) %>%
				dplyr::arrange(dplyr::across(y_cols))
	}

	#dup.X
	dup.X <- count_xUy %>%
		dplyr::filter(n == 1) %>%
		dplyr::count(dplyr::across(names(x_cols)), name = "size") %>%
		dplyr::filter(size > 1) %>%
		dplyr::ungroup() %>%
		dplyr:::select(-size)
	if (nrow(dup.X) > 1) {
		if (!is.null(n_examples) && (n_examples < nrow(dup.X))) {
			dup.X <- dup.X %>% dplyr::sample_n(n_examples, replace = FALSE)
		}
		problems$dup.X <- dup.X %>%
			dplyr::left_join(data, by = names(x_cols)) %>%
			dplyr::arrange(dplyr::across(names(x_cols)))
	}

	#dup.Y
	dup.Y <- count_xUy %>%
		dplyr::filter(n == 1) %>%
		dplyr::count(dplyr::across(names(y_cols)), name = "size") %>%
		dplyr::filter(size > 1) %>%
		dplyr::ungroup() %>%
		dplyr:::select(-size)
	if (nrow(dup.Y) > 1) {
		if (!is.null(n_examples) && (n_examples < nrow(dup.Y))) {
			dup.Y <- dup.Y %>% dplyr::sample_n(n_examples, replace = FALSE)
		}
		problems$dup.Y <- dup.Y %>%
			dplyr::left_join(data, by = names(ycols)) %>%
			dplyr::arrange(dplyr::across(names(y_cols)))
	}

	#dup.XUY
	dup.XUY <- count_xUy %>%
		dplyr::filter(n > 1) %>%
		dplyr:::select(-n)
	if (nrow(dup.XUY) > 1) {
		if (!is.null(n_examples) && (n_examples < nrow(dup.XUY))) {
			dup.XUY <- dup.XUY %>% dplyr::sample_n(n_examples, replace = FALSE)
		}
		problems$dup.XUY <- dup.XUY %>%
			dplyr::left_join(data, by = names(xUy_cols)) %>%
			dplyr::arrange(dplyr::across(names(xUy_cols)))
	}
	if (verbose) {
		cat("Returned instances where:\n")
		cat("\tis.na.X: The X identifier is NA\n")
		cat("\tis.na.Y: The Y identifier is NA\n")
		cat("\tdup.X:   The X identifier is not unique\n")
		cat("\tdup.Y:   The Y identifier is not unique\n")
		cat("\tdup.XUY: The X and Y identifiers together are not unique\n")
	}
	problems
}

