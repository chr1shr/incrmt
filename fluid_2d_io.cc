#include <cstring>
#include <cmath>
#include <vector>

#include "mat.hh"
#include "common.hh"
#include "sim_type.hh"
#include "fluid_2d.hh"

/** Saves the header file.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void fluid_2d::save_header(double duration,int frames) {
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,"%s/header",filename);
    FILE *outf=safe_fopen(bufc,f_num==0?"w":"a");
    fprintf(outf,"%g %g %d\n",time,time+duration,frames);
    fclose(outf);
}

/** Outputs the tracer positions in a binary format that can be read by
 * Gnuplot.
 * \param[in] sn the current frame number to append to the filename. */
void fluid_2d::output_tracers(const int sn) {

    // Assemble the output filename and open the output file
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,"%s/trace.%d",filename,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the tracer positions in batches of 128 floats
    int j,tbatch=(ntrace<<1)/128,tres=(ntrace<<1)%128;
    float *fp,*fe=buf+128;
    double *tp=tm,*te=tm+(ntrace<<2);
    for(j=0;j<tbatch;++j) {
        fp=buf;
        while(fp<fe) {
            *(fp++)=*(tp++);
            *(fp++)=*tp;tp+=3;
        }
        fwrite(buf,sizeof(float),128,outf);
    }

    // Output the remaining tracer positions, if any
    if(tres>0) {
        fp=buf;
        do {
            *(fp++)=*(tp++);
            *(fp++)=*tp;tp+=3;
        } while(tp<te);
        fwrite(buf,sizeof(float),tres,outf);
    }

    // Close the file
    fclose(outf);
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void fluid_2d::write_files(int sn) {

    // Output the velocity components
    if(fflags&1) output("u",0,sn);
    if(fflags&2) output("v",1,sn);
    if(fflags&4) output("spd",2,sn);

    // Output the pressure, vorticity, density, and divergence
    if(fflags&8) output("p",3,sn);
    if(fflags&16) output("w",4,sn);
    if(!const_rho&&(fflags&32)) output("rho",5,sn);
    if(fflags&4096) output("div",12,sn);

    // Output the solid fields
    if(n_obj>0&&(fflags&8640)) {
        if(fflags&512) {
            int_box ib(ml,0,nl,0);
            for(obj_field **op=olist;op<oe;op++) (*op)->bound(ib);
            ib.extend(2);
            ib.trim(0,ml,0,nl);
            if(fflags&64) output("X",6,sn,ib);
            if(fflags&128) output("Y",7,sn,ib);
            if(fflags&256) output("phi",8,sn,ib);
        } else {
            if(fflags&64) output("X",6,sn);
            if(fflags&128) output("Y",7,sn);
            if(fflags&256) output("phi",8,sn);
        }
        if(fflags&8192) output("J",9,sn);
    }

    // Output the diagnostic fields
    if(fflags&2048) {
        output("ds",10,sn);
        output("dp",11,sn);
    }

    // Output tracer positions if there are any
    if(ntrace>0) output_tracers(sn);

    // Do any problem-related output
    for(obj_field **op=olist;op<oe;op++) (*op)->obj->write(sn,time);
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void fluid_2d::output(const char *prefix,const int mode,const int sn) {
    if(fflags&1024) {
        int_box ib(mode==4?1:0,ml,mode==4?1:0,nl);
        output(prefix,mode,sn,ib);
    } else {
        int dis=(1<<mode)&(8|16|512|2048)?1:2;
        int_box ib(2,ml-dis,2,nl-dis);
        output(prefix,mode,sn,ib);
    }
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void fluid_2d::output(const char *prefix,const int mode,const int sn,int_box &ib) {

    // Set constants used in grid output
    int i,j,si=ib.isize(),sj=ib.jsize();
    if(si<=0||sj<=0) return;

    // Assemble the output filename and open the output file
    char *bufc=reinterpret_cast<char*>(buf);
    sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
    FILE *outf=safe_fopen(bufc,"wb");

    // Output the first line of the file
    const bool cen=!((1<<mode)&(8|16|512|2048));
    double lx=ax+dx*(ib.li-(cen?1.5:2.)),
           ly=ay+dy*(ib.lj-(cen?1.5:2.));
    float *bp=buf+1,*be=bp+si;
    *buf=si;
    for(i=0;i<si;i++) *(bp++)=lx+i*dx;
    fwrite(buf, sizeof(float), si+1, outf);

    // Output the field values to the file
    int ijr=ib.li+ml*ib.lj;
    for(j=0;j<sj;++j,ijr+=ml) {
        bp=buf+1;
        *buf=ly+j*dy;
        if(mode<6||mode>9) {

            // Output a line of one the global fields, such as velocity,
            // vorticity, or density
            field *fp=fbase+ijr;
            switch (mode) {
                case 0: while(bp<be) *(bp++)=(fp++)->u;break;
                case 1: while(bp<be) *(bp++)=(fp++)->v;break;
                case 2: while(bp<be) *(bp++)=(fp++)->spd();break;
                case 3: while(bp<be) *(bp++)=(fp++)->p;break;
                case 4: while(bp<be) *(bp++)=vorticity(fp++);break;
                case 5: while(bp<be) *(bp++)=(fp++)->rho;break;
                case 10: while(bp<be) *(bp++)=(fp++)->c0;break;
                case 11: while(bp<be) *(bp++)=(fp++)->c1;break;
                case 12: while(bp<be) *(bp++)=divergence(fp++);
            }
        } else {

            // Output a line of one of the solid-based fields, such as the
            // reference map or the volume deviation
            int ij=ijr;
            switch (mode) {
                case 6: while(bp<be) {*(bp++)=min_phi_obj(ij)->sbase[ij].X;ij++;} break;
                case 7: while(bp<be) {*(bp++)=min_phi_obj(ij)->sbase[ij].Y;ij++;} break;
                case 8: while(bp<be) *(bp++)=min_phi(ij++);break;
                case 9: while(bp<be) *(bp++)=vol_change(ij++);
            }
        }
        fwrite(buf,sizeof(float),si+1,outf);
    }
    fclose(outf);
}

/** Compares the simulation fields to another simulation with higher
 * resolution.
 * \param[in] f2d the higher-resolution simulation.
 * \param[in] err an array in which to store the l2, l1, and l_infinity
 *                comparisons. */
void fluid_2d::l_comparison(fluid_2d &f2d,double *err) {
    static double ste[5]={1,-1./16.,9./16,9./16.,-1/16.};

    // Check the resolutions agree
    const int &pm=f2d.m,&pn=f2d.n;
    if(pm%m!=0||pn%n!=0)
        fatal_error("Grid dimension is not a multiple of the finer grid",1);

    // Initialize accumulators
    int dm=pm/m,dn=pn/n,pml=f2d.ml,num_t=omp_get_max_threads();
    double *o=new double[24*num_t],*stx,*sty;
    for(double *pp=o;pp<o+24*num_t;pp++) *pp=0.;

    // Compute the offsets needed for matching up the cell-centered
    // fields
    int li,ui,lj,uj;
    if(dm&1) {ui=(li=dm/2)+1;stx=ste;}
    else {ui=(li=dm/2-2)+4;stx=ste+1;}
    if(dn&1) {uj=(lj=dn/2)+1;sty=ste;}
    else {uj=(lj=dn/2-2)+4;sty=ste+1;}

    // Loop over this grid, and compare to the fine grid using bilinear
    // interpolation if necessary
#pragma omp parallel for
    for(int j=0;j<(y_prd?n:n+1);j++) {
        double *op=o+24*omp_get_thread_num(),*oq,
                     tf=y_prd||(j>0&&j<n)?1:0.5,tf2,delu,delv,delsq;
        for(int i=0;i<(x_prd?m:m+1);i++) {
            tf2=x_prd||(i>0&&i<m)?tf:0.5*tf;
            field *fp=fm+ml*j+i,*fp2=f2d.fm+pml*dn*j+dm*i,*fp3;
            int ij=ml*(2+j)+(i+2),ii,jj;
            double &del=fp->c0,&delp=fp->c1,*sx,*sy,mphi;

            if(j<n&&i<m) {

                // Interpolate the velocity
                delu=delv=0.;
                for(sy=sty,jj=lj;jj<uj;jj++,sy++) for(sx=stx,ii=li;ii<ui;ii++,sx++) {
                    fp3=fp2+jj*pml+ii;
                    delu+=*sy*(*sx)*fp3->u;
                    delv+=*sy*(*sx)*fp3->v;
                }

                // Compute the field differences
                delu-=fp->u;delv-=fp->v;

                // Compute the velocity contributions
                mphi=min_phi(ij);
                oq=op+(mphi>eps?0:(mphi<-eps?1:2));
                del=sqrt(delsq=delu*delu+delv*delv);
                *oq+=1.;
                oq[6]+=delsq;
                oq[12]+=del;
                if(oq[18]<del) oq[18]=del;
            }

            // Compute material type at this cell corner, assuming all edges
            // are fluid
            mphi=0;
            for(jj=-1;jj<3;jj++) for(ii=-1;ii<3;ii++)
                mphi+=ste[jj+2]*ste[ii+2]*min_phi(ij+jj*ml+ii);
            oq=op+(mphi>eps?0:(mphi<-eps?1:2));

            // Compute the pressure contributions
            delp=fp2->p-fp->p;
            oq[3]+=1.;
            oq[9]+=delp*delp*tf2;
            oq[15]+=fabs(delp)*tf2;
            if(oq[21]<delp) oq[21]=delp;
        }
    }

    // Aggregate the results from the different threads
    double fac=1./(mn),s,*ep=err;
    int k,l=0,a;
    for(a=0;a<6;a++,ep+=4) {
        for(l=0;l<3;l++) {
            s=o[3*a+l];
            for(k=1;k<num_t;k++) s+=o[3*a+l+24*k];
            ep[l]=s*fac;
        }
        ep[3]=*ep+ep[1]+ep[2];
    }
    for(;a<8;a++,ep+=4) {
        for(l=0;l<3;l++) {
            s=o[3*a+l];
            for(k=1;k<num_t;k++) if(s<o[3*a+l+24*k]) s=o[3*a+l+24*k];
            err[4*a+l]=s;
        }
        ep[3]=*ep>ep[1]?*ep:ep[1];
        ep[3]=ep[2]>ep[3]?ep[2]:ep[3];
    }

    // Take the square root of the L2 norms
    for(ep=err+8;ep<err+16;ep++) *ep=sqrt(*ep);
}

/** Calculates the total momentum of the fluid and solid.
 * \param[out] (momx,momy) the total momentum vector. */
void fluid_2d::momentum(double &momx,double &momy) {

    // Loop over all of the grid points and sum up their momenta
#pragma omp parallel for reduction(+:momx,momy)
    for(int j=0;j<n;j++) {
        for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++) {
            momx+=fp->rho*fp->u;
            momy+=fp->rho*fp->v;
        }
    }

    // Scale the result by the area integration factor
    momx*=dx*dy;
    momy*=dx*dy;
}

/** Computes the volumetric change at a corner grid point.
 * \param[in] ij the index of the corner point to consider. */
double fluid_2d::vol_change(int ij) {
    obj_field *op=min_corner_phi_obj(ij);
    if(corner_phi(op->phi_base+ij)>0) return 0;
    s_field *sp=op->sbase+ij;

    // Compute the Jacobian of the reference map field
    double Xx=0.5*xsp*(sp->X+sp[-ml].X-sp[-1].X-sp[-ml-1].X),
           Xy=0.5*ysp*(sp->X+sp[-1].X-sp[-ml].X-sp[-ml-1].X),
           Yx=0.5*xsp*(sp->Y+sp[-ml].Y-sp[-1].Y-sp[-ml-1].Y),
           Yy=0.5*ysp*(sp->Y+sp[-1].Y-sp[-ml].Y-sp[-ml-1].Y);

    // Return the volumetric change, based on 1/det(ref. map Jacobian)-1
    return 1./(Xx*Yy-Xy*Yx)-1;
}
