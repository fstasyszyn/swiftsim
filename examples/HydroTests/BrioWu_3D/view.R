library("fields")
library("akima")
library("rhdf5")

do_tube <- function(file){
hea=h5readAttributes(file,"/Header")
xx=h5read(file,"/PartType0/Coordinates")
rr=h5read(file,"/PartType0/Densities")
vv=h5read(file,"/PartType0/Velocities")
bb=h5read(file,"/PartType0/Bfield")
divb=h5read(file,"/PartType0/divB")
z=xx[3,]
jj=which(z<0.1)
x=xx[1,]
y=xx[2,]
print(c("Tiempo: ", hea$Time))

BB=bb[1,]*bb[1,]+bb[2,]*bb[2,]+bb[3,]*bb[3,]
VV=vv[1,]*vv[1,]+vv[2,]*vv[2,]+vv[3,]*vv[3,]

par(mfrow=c(2,3))
plot(x,rr,main="Density",pch='.',xlim=c(1,3))
plot(x,vv[1,],main="Velocity x",pch='.',xlim=c(1,3))
plot(x,vv[2,],main="Velocity y",pch='.',xlim=c(1,3))
plot(x,vv[3,],main="Velocity z",pch='.',xlim=c(1,3))
plot(x,bb[1,],main="Bfld x",pch='.',xlim=c(1,3))
plot(x,bb[2,],main="Bfld y",pch='.',xlim=c(1,3))
}

do_heat <- function(file){
#file="OrszagTangVortex_0010.hdf5"
#file="OrszagTangVortex_0005.hdf5"
hea=h5readAttributes(file,"/Header")
xx=h5read(file,"/PartType0/Coordinates")
rr=h5read(file,"/PartType0/Densities")
vv=h5read(file,"/PartType0/Velocities")
bb=h5read(file,"/PartType0/Bfield")
divb=h5read(file,"/PartType0/divB")
z=xx[3,]
jj=which(z<0.1)
x=xx[1,]
y=xx[2,]
print(c("Tiempo: ", hea$Time))

BB=bb[1,]*bb[1,]+bb[2,]*bb[2,]+bb[3,]*bb[3,]
VV=vv[1,]*vv[1,]+vv[2,]*vv[2,]+vv[3,]*vv[3,]

adata0=interp(x[jj],y[jj],(divb[jj]),nx=128,ny=128)
adata1=interp(x[jj],y[jj],log10(rr[jj]),nx=128,ny=128)
adata2=interp(x[jj],y[jj],log10(BB[jj]),nx=128,ny=128)
adata3=interp(x[jj],y[jj],log10(VV[jj]),nx=128,ny=128)

par(mfrow=c(2,2))
image.plot(adata0,main=c("DivB at ",hea$Time))
image.plot(adata1,main="Density")
image.plot(adata2,main="B2")
image.plot(adata3,main="V2")

return(0)
}
