# Makefile for navigation corrections software

FC = gfortran
FFLAGS = -g
LDFLAGS = 

# Libraries
NETCDF = /usr/local
LIBS = -L${NETCDF}/lib -lnetcdf -lnetcdff -lcurl -hdf5 -lhdf5_hl
INCLUDES = -I${NETCDF}/include

all: readnetcdf_DBZ_VR cns_eldo_cai

readnetcdf_DBZ_VR.o: readnetcdf_DBZ_VR.f90
	${FC} ${FFLAGS} ${LDFLAGS} -c $? ${INCLUDES}

*.o:
	${FC} ${FFLAGS} ${LDFLAGS} -c $? ${INCLUDES}

cns_eldo_cai: cns_eldo_cai.o chol_inv.o
	${FC} ${FFLAGS} ${LDFLAGS} -o $@ $? ${LIBS} ${INCLUDES}

readnetcdf_DBZ_VR: readnetcdf_DBZ_VR.o
	${FC} ${FFLAGS} ${LDFLAGS} -o $@ $? ${LIBS} ${INCLUDES}

clean:
	rm -f core *.o *~ 

