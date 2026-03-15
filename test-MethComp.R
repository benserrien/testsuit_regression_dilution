# exploratie met het R package "MethComp"
# => dit laat wel toe om met herhaalde metingen te werken


library(MethComp)


data(fat)
str(fat)

vis <- Meth(fat, 2, 1, 3, 5)

pw <- to.wide( vis )
par( mar=c(3,3,1,1) )
with(pw, plot( SL ~ KL, pch=16, xlim=range(vis$y), ylim=range(vis$y) ) )
abline( 0,1 )


plot( vis )
plot( perm.repl( vis ) )

par( mar=c(3,3,3,3), mgp=c(3,1,0)/1.6 )
BA.plot( vis, conn.repl=FALSE )


DA.reg( vis )
