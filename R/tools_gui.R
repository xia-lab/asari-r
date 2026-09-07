# Corresponds to Python asari/tools/gui.py: An injectable R front-end adaptation layer for the experimental desktop GUI.

# Corresponds to run_program.
run_program <- function(params) {
  runner <- get0("run_asari", mode = "function", inherits = TRUE)
  if (!is.function(runner)) stop("run_asari is unavailable.")
  runner(params)
}

# Corresponds to select_directory: the GUI adapter is used first, and interactive R sessions can fall back to tcltk.
select_directory <- function(params) {
  chooser <- getOption("asariR.select_directory")
  directory <- if (is.function(chooser)) {
    chooser(title = "Select Input Directory")
  } else if (interactive() && requireNamespace("tcltk", quietly = TRUE)) {
    as.character(tcltk::tk_choose.dir(caption = "Select Input Directory"))
  } else {
    stop("Directory selection requires an interactive UI adapter.")
  }
  params$input <- directory
  params
}

# Create a mutable object corresponding to TextRedirector.
TextRedirector <- function(text_widget) TextRedirector__init__(text_widget)

# Corresponds to TextRedirector.__init__.
TextRedirector__init__ <- function(text_widget) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("TextRedirector", "environment")
  self$text_widget <- text_widget
  self$write <- function(s) TextRedirector_write(self, s)
  self$flush <- function() TextRedirector_flush(self)
  self
}

# Corresponds to TextRedirector.write.
TextRedirector_write <- function(self, s) {
  widget <- self$text_widget
  if (is.environment(widget) && is.function(widget$insert)) {
    widget$insert("end", s)
    if (is.function(widget$see)) widget$see("end")
  } else if (is.function(widget)) widget(s)
  invisible(NULL)
}

# Corresponds to empty TextRedirector.flush.
TextRedirector_flush <- function(self) { invisible(self); invisible(NULL) }

# Set button state, compatible with Python config or R callback.
.gui_button_state <- function(button, state) {
  if (is.environment(button) && is.function(button$config)) button$config(state = state)
  else if (is.function(button)) button(state)
  invisible(NULL)
}

# Corresponds to run_program_thread: captures the output and writes it to the widget and restores the button.
run_program_thread <- function(params, text_widget, continue_button) {
  redirector <- TextRedirector(text_widget)
  .gui_button_state(continue_button, "disabled")
  on.exit(.gui_button_state(continue_button, "normal"), add = TRUE)
  output <- capture.output(run_program(params), type = "output")
  if (length(output)) redirector$write(paste0(paste(output, collapse = "\n"), "\n"))
  invisible(NULL)
}

# Corresponds to start_program_gui: the front-end adapter can be displayed asynchronously and executed synchronously by default.
start_program_gui <- function(params) {
  adapter <- getOption("asariR.start_program_gui")
  if (is.function(adapter)) adapter(params, run_program_thread) else run_program(params)
  params
}

# The following two functions correspond to the internal callbacks of show_disclaimer.
.gui_on_accept <- function(state) { state$accepted <- TRUE; invisible(NULL) }
.gui_on_decline <- function(state) { state$accepted <- FALSE; invisible(NULL) }

# Corresponds to show_disclaimer.
show_disclaimer <- function() {
  handler <- getOption("asariR.show_disclaimer")
  if (is.function(handler)) return(isTRUE(handler()))
  if (!interactive()) return(FALSE)
  answer <- utils::menu(c("I Accept", "I Decline"), title = "Asari GUI is Experimental")
  identical(answer, 1L)
}

# Corresponds to create_ui internal run_callback: convert the edited value according to the original type.
.gui_run_callback <- function(values, types) {
  result <- list()
  for (key in names(values)) {
    value <- values[[key]]
    type <- types[[key]]
    result[[key]] <- if (identical(type, "logical")) {
      as.logical(value)
    } else if (identical(type, "integer")) {
      parsed <- suppressWarnings(as.integer(value)); if (is.na(parsed)) stop("Invalid integer for ", key); parsed
    } else if (identical(type, "double")) {
      parsed <- suppressWarnings(as.numeric(value)); if (is.na(parsed)) stop("Invalid float for ", key); parsed
    } else as.character(value)
  }
  result
}

# Corresponds to create_ui: only exposes bool/int/float/str/None fields supported by Python GUI.
create_ui <- function(data) {
  supported <- vapply(data, function(value) is.logical(value) || is.integer(value) ||
    is.numeric(value) || is.character(value) || is.null(value), FALSE)
  values <- lapply(data[supported], function(value) if (is.null(value)) "NONE" else value)
  types <- lapply(data[supported], function(value) if (is.null(value)) "character" else typeof(value))
  editor <- getOption("asariR.edit_parameters")
  if (is.function(editor)) values <- editor(values, types)
  result <- .gui_run_callback(values, types)
  lapply(result, function(value) if (identical(value, "NONE")) NULL else value)
}

# Corresponds to create_selection_ui internal set_value.
.gui_set_value <- function(value, data, key) {
  data[[key]] <- value
  data
}

# Corresponds to create_selection_ui.
create_selection_ui <- function(options, data, key) {
  selector <- getOption("asariR.select_option")
  value <- if (is.function(selector)) selector(options, key) else if (interactive()) {
    choice <- utils::menu(options, title = "Select an Option")
    if (choice == 0L) return(list()) else options[[choice]]
  } else stop("Option selection requires an interactive UI adapter.")
  .gui_set_value(value, data, key)
}

# Corresponds to the main_gui overall process.
main_gui <- function() {
  if (!show_disclaimer()) return(invisible(NULL))
  params <- select_directory(PARAMETERS)
  if (is.null(params$input) || !nzchar(params$input)) stop("No directory selected. Exiting.")
  params <- create_ui(params)
  if (is.null(params$autoheight)) stop("Parameters issue (code 1). Exiting.")
  params <- tryCatch({
    updater <- get0("update_peak_detection_params", mode = "function", inherits = TRUE)
    if (!is.function(updater)) stop("update_peak_detection_params unavailable")
    updated <- updater(params)
    create_selection_ui(SUBCOMMANDS, updated, "run_gui")
  }, error = function(error) stop("Parameters issue (code 2). Exiting."))
  tryCatch(start_program_gui(params), error = function(error) {
    stop("Thanks for trying Asari GUI: ", conditionMessage(error))
  })
  invisible(NULL)
}
