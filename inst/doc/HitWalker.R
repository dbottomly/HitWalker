### R code from vignette source 'HitWalker.Rnw'

###################################################
### code chunk number 1: HitWalker.Rnw:17-31
###################################################
stopifnot(require(HitWalker))
stopifnot(require(RSQLite))
stopifnot(require(HitWalkerData))

data(params)

db.con <- dbConnect("SQLite", hitwalker.db.path())

graph.obj <- loadGraph(graph.file.path(), examp.prior.param)
    
test.out <- run.prioritization.patient(db.con, "08-00102", graph.obj, examp.prior.param)

stopifnot(dbDisconnect(db.con))



###################################################
### code chunk number 2: HitWalker.Rnw:39-42
###################################################
    
plot(test.out, examp.graph.param)



###################################################
### code chunk number 3: HitWalker.Rnw:50-57
###################################################

sum.dta <- summary(test.out)

disp.cols <- c("symbol", "protein", "chr", "pos", "ref", "alt", "type.rank") 

head(sum.dta[,disp.cols])



###################################################
### code chunk number 4: HitWalker.Rnw:119-124
###################################################

maxPlotVars(examp.graph.param) <- 5
    
plot(test.out, examp.graph.param)



###################################################
### code chunk number 5: HitWalker.Rnw:140-148
###################################################

plot.new()
HitWalker:::draw.triangle(.2,.75,.25,.2,.2,"white","black",2,1,"up_std")
HitWalker:::draw.triangle(.2,.25,.25,.2,.2,"white","black",2,1,"down_std")

HitWalker:::draw.triangle(.75,.75,.25,.25,.25,"white","black",2,1,"up_half")
HitWalker:::draw.triangle(.75,.25,.25,.25,.25,"white","black",2,1,"down_half")



###################################################
### code chunk number 6: HitWalker.Rnw:160-171
###################################################

plot.new()
use.color <- c("white", "black")
HitWalker:::draw.triangle(.2,.75,.25,.2,.2,use.color,"black",2,1,"up_std")
use.color <- append(use.color, "blue")
HitWalker:::draw.triangle(.2,.25,.25,.2,.2,use.color,"black",2,1,"down_std")
use.color <- append(use.color, "yellow")
HitWalker:::draw.triangle(.75,.75,.25,.25,.25,use.color,"black",2,1,"up_half")
use.color <- append(use.color, "red")
HitWalker:::draw.triangle(.75,.25,.25,.25,.25,use.color,"black",2,1,"down_half")



