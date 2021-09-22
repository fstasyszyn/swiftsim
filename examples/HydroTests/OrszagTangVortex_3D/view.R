library("fields")
library("akima")
library("rhdf5")

do_heat <- function(file){
#file="OrszagTangVortex_0010.hdf5"
#file="OrszagTangVortex_0005.hdf5"
hea=h5readAttributes(file,"/Header")
xx=h5read(file,"/PartType0/Coordinates")
rr=h5read(file,"/PartType0/Densities")
bb=h5read(file,"/PartType0/Bfield")
divb=h5read(file,"/PartType0/divB")
z=xx[3,]
jj=which(z<0.2)
x=xx[1,]
y=xx[2,]
print(c("Tiempo: ", hea$Time))

BB=bb[1,]*bb[1,]+bb[2,]*bb[2,]+bb[3,]*bb[3,]

adata=interp(x[jj],y[jj],(divb[jj]),nx=128,ny=128)
#adata=interp(x[jj],y[jj],log10(rr[jj]),nx=128,ny=128)
#adata=interp(x[jj],y[jj],log10(BB[jj]),nx=128,ny=128)
image.plot(adata)
return(0)
}
