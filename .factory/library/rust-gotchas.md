# Rust Code Gotchas for dlnm

## Unicode in Rust Doc Comments

Do NOT use Unicode characters (e.g., × multiplication sign, ² superscript) in Rust doc comments or string literals that are exposed to R via extendr. R's default encoding (ASCII) cannot handle them, causing `pkgload::load_all()` to fail when parsing the auto-generated `R/extendr-wrappers.R`.

**Fix**: Use ASCII equivalents: `x` instead of `×`, `^2` instead of `²`.

## Group Boundary Pattern

When computing group boundaries from a factor/integer vector for the Rust fused_crossbasis kernel, use:

```r
grp_rle <- rle(as.integer(group))
group_ends <- cumsum(grp_rle$lengths)
group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)
```

This produces 1-indexed start/end pairs for each group. The Rust function expects these as integer vectors.

## Testing R Fallback (with_disabled_rust pattern)

To test that R fallback works when Rust is unavailable, use namespace binding manipulation:

```r
with_disabled_rust <- function(fn_name, code) {
  ns <- asNamespace("dlnm")
  original <- get(fn_name, envir = ns)
  unlockBinding(fn_name, ns)
  assign(fn_name, function(...) stop("Rust unavailable"), envir = ns)
  lockBinding(fn_name, ns)
  on.exit({
    unlockBinding(fn_name, ns)
    assign(fn_name, original, envir = ns)
    lockBinding(fn_name, ns)
  })
  force(code)
}
```

This temporarily replaces the Rust wrapper with a throwing stub, triggering the tryCatch fallback in crossbasis()/crosspred().
