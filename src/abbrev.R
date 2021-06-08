dap <- function() {
  dlba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)
}
drp <- function() {
  dlba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1)
}
dar <- function() {
  dlba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)
}
drr <- function() {
  dlba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)
}
pap <- function() {
  plba_norm(rt, A, b_acc, t0, drifts$AccPrice, 1)
}
prp <- function() {
  plba_norm(rt, A, b_rej, t0, drifts$RejPrice, 1) 
}
par <- function() {
  plba_norm(rt, A, b_acc, t0, drifts$AccRating, 1)
}
prr <- function() {
  plba_norm(rt, A, b_rej, t0, drifts$RejRating, 1)
}

# A macro
*jll"py$jggn
# B macro
df)"pPn
