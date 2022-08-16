here::here()


# gamma distribution
f_distribution <- distcrete::distcrete(
  "gamma",
  shape = 4, 
  scale = 1.25,
  interval = 1,
  w = 1
)
f = f_distribution$d(1:100)


# gamma distribution
w_distribution <- distcrete::distcrete(
  "gamma",
  shape = 56.25,
  scale = 0.26,
  interval = 1,
  w = 1
)
w = w_distribution$d(1:100)

devtools::load_all(here::here())

data <- outbreaker_data(dates = 1:1000, # dates of onset
                w_dens = w, #w.dens, # generation time distribution
                f_dens = f #f.dens # incubation period distribution
)

data[["log_w_dens"]][,999]

res <- outbreaker(data = data)



# 
# ## Replace any 0s with minimum value in w_dens
# if(any(!is.finite(log(w)))){
#   is_positive <- is.finite(log(w))
#   to_replace <- !is_positive
#   val_replace <- min(w[is_positive])
#   w[to_replace] <- val_replace
# }
# 
# 
# ## add an exponential tail summing to 1e-4 to 'w'
# ## to cover the span of the outbreak
# ## (avoids starting with -Inf temporal loglike)
# if (length(w) <1000) {
#   length_to_add <- (1000-length(w)) + 10 # +10 to be on the safe side
#   val_to_add <- stats::dexp(seq_len(length_to_add), 1)
#   val_to_add <- 1e-4*(val_to_add/sum(val_to_add))
#   w <- c(w, val_to_add)
# }


