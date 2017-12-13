select_range <- function(scan_window, rng) {
  if(is.null(selected <- scan_window))
    selected <- rng
  # Make sure selected is within range
  selected[1] <- max(selected[1], rng[1])
  selected[2] <- min(selected[2], rng[2])
  if(selected[1] >= selected[2])
    selected <- rng
  selected
}