# return a plot with txt in the center
plot.text <- function(txt) {
  ggplot(data.frame(x=0,y=0,label=txt)) + geom_text(aes(x=x,y=y,label=label)) + theme_void()
}

# adjust to nearest 100. This is to avoid flickering when size of region changes by small amounts in browser
img.size.round <- function(x) {
  (x %/% 200) * 200
}

# this needs to default to TRUE even if it is missing
use.cache <- reactive({
  if (is.null(input$opt.use.cache)) {
    TRUE
  } else {
    isTruthy(input$opt.use.cache)
  }
})

# create a plot from plot.func, save it as a PNG using the size of the output region, cache it using key,
# and return image info, setting the class to the output.id and its id to key.
renderCacheImage <- function(plot.func, key, width, height=width, opt.use.cache=use.cache(), progress=NULL) {

  write.log(glue("WxH = {width}x{height}"))
  
  if (is.null(width) | is.null(height)) {
    stop("Missing width or height}")
  }
  fn <- glue("{cache.dir}/{key}_{width}_{height}.png")
  
  if (!file.exists(fn) || !opt.use.cache) {
    write.log(glue("Generating plot {fn}"))
    if (!is.null(progress)) progress$set(value=0.1)
    tryCatch({
      a.plot <- plot.func(progress)
      png(fn, width=width, height=height)
      print(a.plot)
      dev.off()
    }, error = function(e) { dev.off(); unlink(fn) })
    if (!is.null(progress)) progress$set(value=1)
  } else {
    if (!is.null(progress)) progress$set(value=0.9, detail="Retrieved cached plot")
    write.log(glue("Retrieving cached {fn}"))
  }
  
  list(src=fn, width=width, height=height, id=key)
}
