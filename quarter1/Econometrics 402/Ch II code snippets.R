#
# Ch II Code snippets
#
library(DataAnalytics)
data(multi)
out=with(multi,
     lm(Sales~p1+p2)
)
lmSumm(out)

#
# demonstrate rgl package by making 3d scatterplot and plotting the plane
#
library(rgl)
data(multi)
with(multi,
     plot3d(p1,p2,Sales,type="s",size=1.5,col="green")
)
out=with(multi,
         lm(Sales~p1+p2)
)
b0=out$coef[1]
b1=out$coef[2]
b2=out$coef[3]
planes3d(b1,b2,-1,b0,alpha=.5,col="light grey")
#  plots a plane with equation  ax + by + cz +d=0  x is p1, y is p2, z is "Sales"
#
# let's draw lines on the plane to illustrate the slope of plane wrt to p1, holding p2 constant
p2star=8
z0=b0+b1*min(p1)+b2*p2star
z1=b0+b1*max(p1)+b2*p2star
x=c(min(p1),max(p1))
y=c(p2star,p2star)
z=c(z0,z1)
lines3d(x,y,z,lwd=2,col="blue")
#
# now slope of plane wrt to p2, holding p1 constant
p1star=5
z0=b0+b1*p1star+b2*min(p2)
z1=b0+b1*p1star+b2*max(p2)
y=c(min(p2),max(p2))
x=c(p1star,p1star)
z=c(z0,z1)
lines3d(x,y,z,lwd=2,col="red")

#
# multifactor models
#
library(reshape2)
data(Vanguard)
Van=Vanguard[,c(1,2,5)]   # grab relevant cols
V_reshaped=dcast(Van,date~ticker,value.var="mret")
data(riskFactors)
Van_risk=merge(V_reshaped,riskFactors,by="date")
outsl=lm(VWNFX ~ RmRf,data=Van_risk)
outml=lm(VWNFX ~ RmRf + SMB + HML,data=Van_risk)

