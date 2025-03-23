## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
  collapse=TRUE,
  warning=FALSE,
  message=FALSE,
  fig.height=8, fig.width=8,
  comment="#>",
  fig.path="jambaVignette-"
);
library(jamba);

## ----mixedSort, echo=FALSE----------------------------------------------------
x <- sort(c(
  "miR-12","miR-1","miR-122","miR-1b",
  "miR-1a","miR-2", "miR-22",
  "ABCA2", "ABCA12"));
df1 <- data.frame(
  miRNA=x,
  sort_rank=seq_along(x),
  mixedSort_rank=order(jamba::mixedOrder(x)),
  check.names=FALSE,
  stringsAsFactors=FALSE);
df2 <- jamba::mixedSortDF(df1);
df2;

## ----plotSmoothScatter, fig.show='hold', fig.alt="Example of smooth scatterplots using the default R function, and enhanced function in the jamba package."----
plotSmoothScatter(doTest=TRUE);

## ----imageByColors1, fig.show='hold', fig.alt="Example color matrix as plotted using the image by colors function()"----
a1 <- c("red","blue")[c(1,1,2)];
b1 <- c("yellow","orange")[c(1,2,2)];
c1 <- c("purple","orange")[c(1,2,2)];
d1 <- c("purple","green")[c(1,2,2)];
df1 <- data.frame(a=a1, b=b1, c=c1, d=d1);
imageByColors(df1, cellnote=df1);

## ----imageByColors2, fig.show='hold', fig.alt="Example color matrix with rotated and grouped labels as plotted using the image by colors function."----
imageByColors(df1,
   cellnote=df1, 
   useRaster=TRUE,
   #adjBy="column",
   cexCellnote=list(c(1.5,1.5,1),
      c(1,1.5), 
      c(1.6,1.2), 
      c(1.6,1.5)),
   srtCellnote=list(c(90,0,0), 
      c(0,45), 
      c(0,0,0), 
      c(0,90,0)));

## ----log_ticks_1, fig.show='hold', fig.alt="Example showing log tick axis labels with positive and negative values, and P-values."----
# example showing volcano plot features
set.seed(123);
n <- 1000;
vdf <- data.frame(lfc=rnorm(n) * 2)
vdf$`-log10 (padj)` <- abs(vdf$lfc) * abs(rnorm(n))
plotSmoothScatter(vdf, xaxt="n", yaxt="n", xlab="Fold change",
   main="Volcano plot\ndisplayBase=2")
logFoldAxis(1)
pvalueAxis(2)


## ----plot_dens_log, fig.show='hold', fig.alt="Example polygon density plot."----
withr::with_par(list("mar"=c(6, 4, 4, 2), "mfrow"=c(2, 2)), {
   withr::local_options(list("scipen"=7));
set.seed(123);
plotPolygonDensity(10^(3+rnorm(2000)),
   breaks=50,
   cex.axis=1,
   main="normal-scaled x-axis");
plotPolygonDensity(10^(3+rnorm(2000)),
   log="x",
   breaks=50,
   main="log-scaled x-axis");
plotPolygonDensity((3+rnorm(2000))^2,
   cex.axis=1,
   breaks=50,
   main="normal-scaled x-axis");
plotPolygonDensity((3+rnorm(2000))^2,
   cex.axis=1,
   xScale="sqrt",
   breaks=50,
   main="");
drawLabels(preset="topright",
   txt="sqrt-scaled x-axis",
   labelCex=1.5)
})

## ----draw_labels, fig.show='hold', fig.alt="Example showing text labels, and shadow text labels on a colored background."----
par("mfrow"=c(1,1))
plotPolygonDensity((3+rnorm(2000))^2,
   cex.axis=1,
   xScale="sqrt",
   breaks=50,
   main="");
drawLabels(preset="bottom",
   txt="sqrt-scaled x-axis",
   text_fn=jamba::shadowText,
   labelCex=1.5)

## ----printdebug, results='asis'-----------------------------------------------
printDebugHtml("printDebugHtml(): ",
  "Output is colorized: ",
  head(LETTERS, 8))

withr::with_options(list(jam.htmlOut=TRUE, jam.comment=FALSE), {
  printDebugHtml(c("printDebug() using withr::with_options(): "),
    c("Output should be colorized: "),
    head(LETTERS, 8));
})

## ----get_color_ramp, fig.show='hold', fig.alt="Example of color gradients using the Brewer package."----
showColors(list(
   Reds=getColorRamp("Reds"),
   RdBu=getColorRamp("RdBu"),
   RdBu_r=getColorRamp("RdBu_r"),
   `RdBu_r, lens=5`=warpRamp(getColorRamp("RdBu_r"), lens=5),
   `RdBu_r, lens=-5`=warpRamp(getColorRamp("RdBu_r"), lens=-5),
   `rainbow2(15)`=rainbow2(15)
   ));

## ----setPrompt, echo=TRUE-----------------------------------------------------
setPrompt("jambaVignette");
# {jambaVignette}-R-3.6.0_10789>

## ----jargs, message=TRUE------------------------------------------------------
# all args
jargs(plotSmoothScatter)

# args with "y" in the name
jargs(plotSmoothScatter, "^y")


## ----sdim---------------------------------------------------------------------
L <- list(LETTERS=LETTERS,
   letters=letters,
   lettersDF=data.frame(LETTERS, letters));
sdim(L);

L2 <- list(List1=L,
   List2=L);

sdim(L2);
ssdim(L2)

## ----paste_order--------------------------------------------------------------
a1 <- factor(c("mutant", "control")[c(1,1,2)],
   levels=c("control", "mutant"));
b1 <- factor(c("vehicle", "treated")[c(2,1,1)],
   levels=c("vehicle", "treated"));
d1 <- c("purple","green")[c(1,2,2)];
df2 <- data.frame(a=a1, b=b1, d=d1);
df2;
pasteByRow(df2);
pasteByRowOrdered(df2);

df3 <- data.frame(df2,
   pasteByRowOrdered=pasteByRowOrdered(df2));
mixedSortDF(df3, byCols="pasteByRowOrdered")


## ----make_names---------------------------------------------------------------
x <- rep(head(letters, 4), c(2,4,1,5));
x;

makeNames(x);

nameVector(x);

y <- nameVector(x);
nameVectorN(y);

lapply(nameVectorN(head(y)), function(i){
   i
})

## ----cpaste-------------------------------------------------------------------
set.seed(123);
x <- lapply(seq_len(6), function(i){
   paste0("Gene", 
      sample(LETTERS,
         sample(c(1,1,2,5,9), 1),
         replace=TRUE));
});
cPaste(x);
cPasteU(x);
cPasteSU(x);

data.frame(cPaste=cPaste(x),
   cPasteU=cPasteU(x),
   cPasteSU=cPasteSU(x))

## ----kable_coloring, results='asis'-------------------------------------------
expt_df <- data.frame(
  Sample_ID="",
  Treatment=rep(c("Vehicle", "Dex"), each=6),
  Genotype=rep(c("Wildtype", "Knockout"), each=3),
  Rep=paste0("rep", c(1:3)))
expt_df$Sample_ID <- pasteByRow(expt_df[, 2:4])

# define colors
colorSub <- c(Vehicle="palegoldenrod",
  Dex="navy",
  Wildtype="gold",
  Knockout="firebrick",
  nameVector(color2gradient("grey48", n=3, dex=10), rep("rep", 3), suffix=""),
  nameVector(
    color2gradient(n=3,
      c("goldenrod1", "indianred3", "royalblue3", "darkorchid4")),
    expt_df$Sample_ID))
if (requireNamespace("kableExtra", quietly=FALSE)) {
   kbl <- kable_coloring(
     expt_df,
     caption="Experiment design table showing categorical color assignment.",
     colorSub)
}

