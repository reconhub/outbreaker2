## Unexported data prep (formerly add_convolutions); used by several tests.
add_convolutions <- function(data, config) {
  outbreaker2:::update_data_with_config(data = data, config = config)
}
