# Suppress R CMD check NOTE about undefined global variables used in NSE
utils::globalVariables(c(
  "drug", "protein", "dose", "x", "y", "y_pred"
))

# Utility functions (if needed in future)

# custom step interpolation (left-continuous)
step_interpolate = function(x, y, new_x) {
  order_idx = order(x)
  x = x[order_idx]
  y = y[order_idx]
  sapply(new_x, function(xi) {
    idx = max(which(x <= xi), na.rm = TRUE)
    return(y[idx])
  })
}


# F-test utility function for isotonic regression model significance
f_test_isotonic = function(x, y, y_pred) {
  n = length(y)
  null_pred = rep(mean(y), n)

  sse_full = sum((y - y_pred)^2)
  sse_null = sum((y - null_pred)^2)

  df_full = n - length(unique(x))
  df_null = n - 1

  f_stat = ((sse_null - sse_full) / (df_null - df_full)) / (sse_full / df_full)
  p_value =  pf(f_stat, df_null - df_full, df_full)

  return(list(
    SSE_Full = sse_full,
    SSE_Null = sse_null,
    F_statistic = f_stat,
    P_value = p_value
  ))
}
