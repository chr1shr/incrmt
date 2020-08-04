#include "heap.hh"
#include "levelset.hh"

/** The levelset constructor.
 * \param[in] m_ the number of horizontal gridpoints.
 * \param[in] n_ the number of vertical gridpoints.
 * \param[in] ax_ the minimum horizontal coordinate.
 * \param[in] ay_ the minimum vertical coordinate.
 * \param[in] dx_ the horizontal grid spacing.
 * \param[in] dy_ the vertical grid spacing.
 * \param[in] phi_down_ the negative phi cutoff value for the narrow band.
 * \param[in] phi_up_ the positive phi cutoff value for the narrow band. */
levelset::levelset(int m_,int n_,double ax_,double ay_,double dx_,double dy_,double phi_down_,double phi_up_) :
    m(m_), n(n_), mn(m*n), me(m+1), ne(n+1), mne(me*ne), ax(ax_), ay(ay_),
    bx(ax+m*dx_), by(ay+n*dy_), dx(dx_), dy(dy_), xsp(1/dx), ysp(1/dy),
    dx2(dx*dx), dy2(dy*dy), phi_down(phi_down_), phi_up(phi_up_),
    u_mem(u_init_memory), bpc(mn+1), sgc(1), mem_t1(temp_init_memory),
    mem_p1(temp_init_memory), mem_up(heap_init_memory),
    mem_down(heap_init_memory), s(new int[mn]), ss(new int[mne]),
    bp(new unsigned int[mn]), sg(new unsigned int[mne]), phi(new double[mn]),
    vp(new double[mn]), u(new int[u_mem]), h_up(new int[mem_up]),
    h_down(new int[mem_down]), t1(new int[mem_t1]), p1(new double[mem_p1]),
    sort(phi), failsafe(false) {
    for(unsigned int *p=bp;p<bp+mn;p++) *p=0;
    for(unsigned int *p=sg;p<sg+mne;p++) *p=0;
}

/** The levelset destructor is responsible for clearing the dynamically
 * allocated memory. */
levelset::~levelset() {
    delete [] h_down;delete [] h_up;
    delete [] p1;delete [] t1;
    delete [] u;
    delete [] sg;delete [] bp;
    delete [] ss;delete [] s;
    delete [] vp;delete [] phi;
}

/** Doubles the size of the temporary integer array. */
void levelset::add_memory_t1() {
#ifdef MEMORY_MESSAGES
    fputs(stderr,"Add temporary integer memory\n");
#endif
    mem_t1<<=1;
    int *pt1(new int[mem_t1]);
    for(int i=0;i<q1;i++) pt1[i]=t1[i];
    delete [] t1;
    t1=pt1;
}

/** Doubles the size of the temporary floating point array. */
void levelset::add_memory_p1() {
#ifdef MEMORY_MESSAGES
    fprintf(stderr,"Add temporary floating point memory\n");
#endif
    mem_p1<<=1;
    double *pp1(new double[mem_p1]);
    for(int i=0;i<r1;i++) pp1[i]=p1[i];
    delete [] p1;
    p1=pp1;
}

/** Carries out a bicubic interpolation of the level set function at a point,
 * also calculating the level set gradient at that point.
 * \param[in] (x,y) the position to interpolate at.
 * \param[out] (phix,phiy) the position */
double levelset::bicubic(double x,double y,double& phix,double& phiy) {
    int i,j,ij;
    const double dxf=0.5,dyf=0.5,d2f=0.25;
    double temp0,temp1,temp2,temp3;
    double plbx,prbx,pltx,prtx;
    double plby,prby,plty,prty;
    double plbxy,prbxy,pltxy,prtxy;
    double me,mf;
    double e=(x-ax)*xsp;
    double f=(y-ay)*ysp;
    i=int(e);j=int(f);
#include "built/bicubic.cc"
}

/** Carries out a bilinear interpolation of any function at a given point.
 * \param[in] fe a pointer to the function to interpolate.
 * \param[in] (x,y) the position to interpolate at.
 * \return the interpolated value. */
double levelset::bilinear(double *fe,double x,double y) {
    int i,j,ij;
    double me,e=(x-ax)*xsp,f=(y-ay)*ysp;
    i=int(e);
    if(i<0) i=0;
    if (i>m-2) i=m-2;
    j=int(f);
    if(j<0) j=0;
    if (j>n-2) j=n-2;
    e-=double(i);f-=double(j);me=1-e;ij=i+m*j;
    return (1-f)*(me*fe[ij]+e*fe[ij+1])+f*(me*fe[ij+m]+e*fe[ij+m+1]);
}

/** Carries out the modified Newton--Raphson iteration of Chopp, to find a
 * point where the bicubic interpolation of the level set function is zero,
 * which is closest to a given grid point.
 * \param[in] i the horizontal grid index.
 * \param[in] ij the grid index.
 * \param[out] (x0,y0) the coordinates of the input grid point.
 * \param[out] (x,y) the coordinates of the closest point where the bicubic
 *                   interpolation vanishes.
 * \return True if the result was a success, false otherwise. */
bool levelset::newton(int i,int ij,double &x0,double &y0,double &x,double &y) {
    const double dyy=dy/(double(m)),step_limit=sqrt(dx*dx+dy*dy);
    int iter=0,cutoff=0;
    double phi,phix,phiy,delx,dely,m0,m1;
    x0=x=ax+i*dx;y0=y=ay+(ij-i)*dyy;
    if(failsafe) return newton_failsafe(i,ij,x0,y0,x,y);
    phi=bicubic(x,y,phix,phiy);
    do {
        delx=x0-x;dely=y0-y;
        m0=delx*phix+dely*phiy+phi;
        m1=phix*phix+phiy*phiy;
        iter++;
        m0/=m1;
        if(iter>100||!std::isfinite(m0)||m0>step_limit) return newton_failsafe(i,ij,x0,y0,x,y);
        x+=delx-m0*phix;
        y+=dely-m0*phiy;
        phi=bicubic(x,y,phix,phiy);
        if(abs(phi)<1e-10) cutoff++;else cutoff=0;
    } while (cutoff<8);
    return abs(x0-x)<3*dx&&abs(y0-y)<3*dy?true:newton_failsafe(i,ij,x0,y0,x,y);
}

/** In cases where the Newton iteration routine fails, this routine provides a
 * lower-order backup following the first order interpolation procedure
 * described in chapter 11 of "Level Set Methods and Fast Marching Methods".
 * \param[in] i the horizontal grid index.
 * \param[in] ij the grid index.
 * \param[in] (x0,y0) the coordinates of the input grid point.
 * \param[out] (x,y) the coordinates of the closest point where the
 *                   interpolation of the level set vanishes.
 * \return True if the result was a success, false otherwise. */
bool levelset::newton_failsafe(int i,int ij,double &x0,double &y0,double &x,double &y) {
    int ts=7-s[ij];bool hf,vf;double tdis,cphi=phi[ij];

    // Look at gridpoints to the left and right
    if(i>0&&s[ij-1]==ts) {
        x=cphi/(phi[ij-1]-cphi);
        if(i<m-1&&s[ij+1]==ts) {
            tdis=cphi/(cphi-phi[ij+1]);
            if(tdis<-x) x=tdis;
        }
        hf=true;
    } else {
        if(i<m-1&&s[ij+1]==ts) {
            x=cphi/(cphi-phi[ij+1]);
            hf=true;
        } else hf=false;
    }

    // Look at gridpoints up and down
    if(ij>=m&&s[ij-m]==ts) {
        y=cphi/(phi[ij-m]-cphi);
        if(ij<mn-m&&s[ij+m]==ts) {
            tdis=cphi/(cphi-phi[ij+m]);
            if(tdis<-y) y=tdis;
        }
        vf=true;
    } else {
        if(ij<mn-m&&s[ij+m]==ts) {
            y=cphi/(cphi-phi[ij+m]);
            vf=true;
        } else vf=false;
    }

    if(hf) {
        if(vf) {

            // If there's a valid intersection in both the
            // horizontal and vertical directions, then we do the
            // diagonal calculation
            y*=dy;tdis=x*dx;
            cphi=tdis*y/(tdis*tdis+y*y);
            x=y*cphi+x0;y=tdis*cphi+y0;
            return true;
        } else {

            // There is only an intersection in the horizontal
            // direction
            x=x*dx+x0;y=y0;return true;
        }
    } else {
        if(vf) {

            // The is only an intersection in the vertical
            // direction
            x=x0;y=y*dy+y0;return true;
        } else return false;
    }
}

/** Builds a narrow band level set, by scanning the phi array and looking for
 * the zero contour. */
void levelset::build_band() {
    double x0,y0,x,y;
    int i,j,ij,k,l;

    // Reset the band position counters
    u_end=0;u_start=u_mem;u_mid=0;

    // Reset the temporary buffer counters
    q1=r1=0;

    // Scan the grid and search for all edges where the level set function
    // changes sign. First consider the bottom row.
    s[0]=phi[0]>0?7:0;
    for(ij=1;ij<m;ij++) {
        if(phi[ij]>0) {
            if(s[ij-1]<4) {
                if (s[ij-1]==0) {s[ij-1]=3;add_negative(ij-1);}
                s[ij]=4;add_positive(ij);
            } else s[ij]=7;
        } else {
            if(s[ij-1]<4) s[ij]=0;
            else {
                if (s[ij-1]==7) {s[ij-1]=4;add_positive(ij-1);}
                s[ij]=3;add_negative(ij);
            }
        }
    }

    // Now consider the rest of the grid
    for(j=1;j<n;j++) {
        if(phi[ij]>0) {
            if(s[ij-m]<4) {
                if (s[ij-m]==0) {s[ij-m]=3;add_negative(ij-m);}
                s[ij]=4;add_positive(ij);
            } else s[ij]=7;
        } else {
            if(s[ij-m]<4) s[ij]=0;
            else {
                if (s[ij-m]==7) {s[ij-m]=4;add_positive(ij-m);}
                s[ij]=3;add_negative(ij);
            }
        }
        ij++;
        for(i=1;i<m;i++,ij++) {
            if(phi[ij]>0) {
                if(s[ij-1]<4) {
                    if (s[ij-1]==0) {s[ij-1]=3;add_negative(ij-1);}
                    if (s[ij-m]==0) {s[ij-m]=3;add_negative(ij-m);}
                    s[ij]=4;add_positive(ij);
                } else if (s[ij-m]<4) {
                    if (s[ij-m]==0) {s[ij-m]=3;add_negative(ij-m);}
                    s[ij]=4;add_positive(ij);
                } else s[ij]=7;
            } else {
                if(s[ij-1]>=4) {
                    if (s[ij-1]==7) {s[ij-1]=4;add_positive(ij-1);}
                    if (s[ij-m]==7) {s[ij-m]=4;add_positive(ij-m);}
                    s[ij]=3;add_negative(ij);
                } else if (s[ij-m]>=4) {
                    if (s[ij-m]==7) {s[ij-m]=4;add_positive(ij-m);}
                    s[ij]=3;add_negative(ij);
                } else s[ij]=0;
            }
        }
    }

    // Store initial counts of the number of points either side of the
    // interface
    num_n=u_start;
    num_p=u_end;

    // Compute new values for points either side of the interface by calling
    // the Newton iteration. Store the values in the temporary floating point
    // array
    for(k=0;k<u_end;k++) {
        ij=u[k];i=ij%m;
        if (!newton(i,ij,x0,y0,x,y)) levelpp_fatal_error("Newton iteration failed.",1);
        x-=x0;y-=y0;
        if(r1==mem_p1) add_memory_p1();
        p1[r1++]=sqrt(x*x+y*y);
    }
    for(k=u_mem-1;k>=u_start;k--) {
        ij=u[k];i=ij%m;
        if (!newton(i,ij,x0,y0,x,y)) levelpp_fatal_error("Newton iteration failed.",1);
        x-=x0;y-=y0;
        if(r1==mem_p1) add_memory_p1();
        p1[r1++]=-sqrt(x*x+y*y);
    }

    // Update the phi field with the new values
    for(k=l=0;k<u_end;k++) phi[u[k]]=p1[l++];
    for(k=u_mem-1;k>=u_start;k--) phi[u[k]]=p1[l++];

    // Set the counters of the number and negative and postive boundary points.
    // num_n is initially set to the memory reference of the last negative
    // point. This value is used when the fast marched points are meshed with
    // the boundary points.
    num_p=u_end;
    num_n=u_start;

    // Sort the points in the band
    sort.sort(u,u_end);
    sort.sort(u+u_start,u_mem-u_start);

    // Positive fast march
    heap<dir_positive> up(*this,mem_up,h_up);
    for(k=0;k<u_end;k++) {
        ij=u[k];i=ij%m;
        up.stack_test(i,ij);
    }
    up.build_heap();
    up.march(phi_up);

    // Negative fast march
    heap<dir_negative> down(*this,mem_down,h_down);
    for(k=u_start;k<u_mem;k++) {
        ij=u[k];i=ij%m;
        down.stack_test(i,ij);
    }
    down.build_heap();
    down.march(phi_down);

    // Sort the boundary between the initial points and the fast marched points
    sort_boundary(num_n,u_start,u_mem);
    sort_boundary(num_p,0,u_end);

    // Clean up rest of phi grid
#ifdef GRID_CLEANING
    clean();
#endif

    // Fix the number of negative boundary points
    num_n=u_mem-num_n;

    // Just in case there was no interface, set u_start to zero
    if (u_start==u_mem) u_start=0;
}

void levelset::sort_boundary(int uv,int u_min,int u_max) {
    int k,l,us;

    // Check to make sure that there are actually some values to sort
    if(uv==u_max) return;

    k=us=search_down_simple(uv,u_min,phi[u[uv]]);
    l=uv;

    q1=0;
    while(k<uv) {
        if(q1==mem_t1) add_memory_t1();
        if(phi[u[k]]<phi[u[l]]) {
            t1[q1++]=u[k++];
        } else {
            t1[q1++]=u[l++];
            if(l==u_max) {
                while(k<uv) {
                    if(q1==mem_t1) add_memory_t1();
                    t1[q1++]=u[k++];
                }
                break;
            }
        }
    }

    k=0;
    while(k<q1) u[us++]=t1[k++];
}

void levelset::sort_boundary_wrap() {
    int k,l,us;

    // Check to make sure that there are actually some values to sort
    if(u_end==0||u_start==0||u_start==u_mem) return;

    k=us=search_down_simple(u_mem,u_start,phi[u[0]]);
    l=0;

    q1=0;
    while(k<u_mem) {
        if(q1==mem_t1) add_memory_t1();
        if(phi[u[k]]<phi[u[l]]) {
            t1[q1++]=u[k++];
        } else {
            t1[q1++]=u[l++];
            if(l==u_end) {
                while(k<u_mem) {
                    if(q1==mem_t1) add_memory_t1();
                    t1[q1++]=u[k++];
                }
                break;
            }
        }
    }

    k=0;while(us<u_mem) u[us++]=t1[k++];
    us=0;while(k<q1) u[us++]=t1[k++];
}

int levelset::search_down_simple(int uv,int u_min,double tphi) {
    int k=0,l=1;
    while(uv-l>u_min) {
        if(phi[u[uv-l]]<tphi) return bisect(uv-l,uv-k,tphi);
        k=l;l<<=1;
    }
    return phi[u[u_min]]<tphi?bisect(u_min,uv-k,tphi):u_min;
}

/** This routine will clean up the grid, setting all points that are outside
 * the narrow band to their extremal values. */
void levelset::clean() {
    for(int ij=0;ij<mn;ij++) {
        if(s[ij]<2) {phi[ij]=phi_down;vp[ij]=0;}
        if(s[ij]>5) {phi[ij]=phi_up;vp[ij]=0;}
    }
}

/** This diagnostic routine will scan the u[] array and print error messages
 * for any pairs of entries which are out of order. */
void levelset::check_sorted() {
    int i;
    if (u_end>u_start) {
        for(i=u_start;i<u_end-1;i++) check_two_points(i,i+1);
    } else {
        for(i=u_start;i<u_mem-1;i++) check_two_points(i,i+1);
        if(u_end>0) {
            check_two_points(u_mem-1,0);
            for(i=0;i<u_end-1;i++) check_two_points(i,i+1);
        }
    }
}

/** Checks that the narrow band is correctly divided into a negative half and a
 * positive half. */
void levelset::check_divided() {
    if(u_end!=u_mid) if(phi[u[u_mid]]<-tolerance) {
        printf("Minus number on positive side\n%d %d %d\n",u_start,u_mid,u_end);
        int i=u_mid,j,ij;
        for(j=0;j<4;j++) i=ul_down(i);
        for(j=0;j<9;j++) {ij=u[j];printf("%d %d %g %d %d\n",i,ij,phi[ij],ij%m,ij/m);}
    }

    if(u_start!=u_mid) if(phi[u[ul_down(u_mid)]]>tolerance) {
        printf("+> %d %g %g\nPlus number on negative side\n%d %d %d\n",
               ul_down(u_mid),phi[ul_down(u_mid)],tolerance,u_start,u_mid,u_end);
        int i=u_mid,j,ij;
        for(j=0;j<4;j++) i=ul_down(i);
        for(j=0;j<9;j++) {ij=u[i];printf("%d %d %g %d %d\n",i,ij,phi[ij],ij%m,ij/m);i=ul_up(i);}
    }
}

/** Saves the band structure to a file for debugging and analysis.
 * \param[in] filename the name of the file to write to. */
void levelset::output_band(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) levelpp_fatal_error("Can't save level set band file",1);
    int k=u_start,j=0,ij;
    while(k!=u_end) {
        ij=u[k];
        fprintf(fp,"%d %d %d %d %d %d %g",k,j,ij,ij%m,ij/m,s[ij],phi[ij]);
        k=ul_up(k);j++;
    }
    fclose(fp);
}

/** Checks that the internal counts of the number of negative and positive
 * boundary points agree with the grid. */
void levelset::check_counts() {
    int anum_n=0,anum_p=0,ij=0;
    while(ij<mn) {
        if(s[ij]==3) anum_n++;
        if(s[ij]==4) anum_p++;
        ij++;
    }
    printf("Thinks: num_n=%d, num_p=%d\n"
           "Actual: num_n=%d, num_n=%d\n",num_n,num_p,anum_n,anum_p);
}

/** This small routine checks to see if two points in the u[] array are out of
 * order and prints an error message if they are.
 * \param[in] i the first point to test.
 * \param[in] j the second point to test. */
void levelset::check_two_points(int i,int j) {
    if(phi[u[i]]>phi[u[j]]+tolerance)
        printf("Error: phi[u[%d]]=%g and phi[u[%d]]=%g\n%d %d %d %d %d %d\n",
               i,phi[u[i]],j,phi[u[j]],u[i],u[j],u[i]%m,u[i]/m,u[j]%m,u[j]/m);
}

/** Adds a point to the positive side of the narrow band list, increasing
 * memory if necessary.
 * \param[in] ij the gridpoint index to add. */
void levelset::add_positive(int ij) {
    u[u_end]=ij;u_end=ul_up(u_end);
    if (u_end==u_start) add_u_memory_init();
}

/** Adds a point to the negative side of the narrow band list, increasing
 * memory if necessary.
 * \param[in] ij the gridpoint index to add. */
void levelset::add_negative(int ij) {
    u_start=ul_down(u_start);u[u_start]=ij;
    if (u_end==u_start) add_u_memory_init();
}

/** Increases the memory for the narrow band. This routine can only be called
 * during the initial band-building phase as it also adjusts the num_n pointer,
 * but does not deal with back pointers. */
void levelset::add_u_memory_init() {
#ifdef MEMORY_MESSAGES
    puts("Add list memory");
#endif
    if(u_mem==mn+1) levelpp_fatal_error("Narrow band memory growing beyond the total number of gridpoints.",1);

    int i,u_nmem=(u_mem<<1);
    if(u_nmem>mn+1) u_nmem=mn+1;
    int dis=(u_nmem-u_mem),*pu=new int[u_nmem];
    for(i=0;i<u_end;i++) pu[i]=u[i];
    while(i<u_mem) {
        pu[i+dis]=u[i];i++;
    }
    u_start+=dis;
    num_n+=dis;
    delete [] u;
    u_mem=u_nmem;
    u=pu;
}

int levelset::search_up(int iu,double iphi) {
#ifdef SIMPLE_SEARCH
    while(iu!=u_end) {
        if(phi[u[iu]]>iphi) return iu;
        iu=ul_up(iu);
    }
    return u_end;
#else
    int k=0,l=1;
    if(u_end<iu) {
        while(iu+l<u_mem) {
            if(phi[u[iu+l]]>iphi) return bisect(iu+k,iu+l,iphi);
            k=l;l<<=1;
        }
        if(phi[u[u_mem-1]]>iphi) return bisect(iu+k,u_mem-1,iphi);
        if(u_end==0||phi[u[0]]>iphi) return 0;
        k=0;l=1;iu=0;
    } else if(u_end==iu) return iu;
    while(iu+l<u_end) {
        if(phi[u[iu+l]]>iphi) return bisect(iu+k,iu+l,iphi);
        k=l;l<<=1;
    }
    return bisect(iu+k,u_end,iphi);
#endif
}

int levelset::search_down(int iu,double iphi) {
#ifdef SIMPLE_SEARCH
    while(iu!=u_start) {
        iu=ul_down(iu);
        if(phi[u[iu]]<iphi) return ul_up(iu);
    }
    return u_start;
#else
    int k=0,l=1;
    if(u_start>iu) {
        while(iu>=l) {
            if(phi[u[iu-l]]<iphi) return bisect(iu-l,iu-k,iphi);
            k=l;l<<=1;
        }
        if(phi[u[0]]<iphi) return bisect(0,iu-k,iphi);
        if(phi[u[u_mem-1]]<iphi) return 0;
        k=0;l=1;iu=u_mem-1;
    } else if(u_start==iu) return iu;
    while(iu-l>=u_start) {
        if(phi[u[iu-l]]<iphi) return bisect(iu-l,iu-k,iphi);
        k=l;l<<=1;
    }
    return phi[u[u_start]]<iphi?bisect(u_start,iu-k,iphi):u_start;
#endif
}

int levelset::bisect(int lu,int uu,double iphi) {
    int i;
    while(uu-lu>1) {
        i=(lu+uu)>>1;
        if (phi[u[i]]<iphi) lu=i;else uu=i;
    }
    return uu;
}

void levelset::change_to_neg(int i,int ij) {
    int e;
    if(s[ij]==4) num_p--;
    s[ij]=3;
    if (ij<m||check_neg(i,ij-m)) e=1;else e=0;
    if (i==0||check_neg(i-1,ij-1)) e++;
    if (i==m-1||check_neg(i+1,ij+1)) e++;
    if (ij>=mn-m||check_neg(i,ij+m)) e++;
    if (e==4) s[ij]=2;else num_n++;
}

void levelset::change_to_pos(int i,int ij) {
    int e;
    if(s[ij]==3) num_n--;
    s[ij]=4;
    if (ij<m||check_pos(i,ij-m)) e=1;else e=0;
    if (i==0||check_pos(i-1,ij-1)) e++;
    if (i==m-1||check_pos(i+1,ij+1)) e++;
    if (ij>=mn-m||check_pos(i,ij+m)) e++;
    if (e==4) s[ij]=5;else num_p++;
}

bool levelset::check_up(int i,int ij) {
    if(s[ij]==5) return false;
    s[ij]=(ij>=m&&s[ij-m]==5)||(i>0&&s[ij-1]==5)||(i<m-1&&s[ij+1]==5)||(ij<mn-m&&s[ij+m]==5)?6:7;
    return true;
}

bool levelset::check_down(int i,int ij) {
    if(s[ij]==2) return false;
    s[ij]=(ij>=m&&s[ij-m]==2)||(i>0&&s[ij-1]==2)||(i<m-1&&s[ij+1]==2)||(ij<mn-m&&s[ij+m]==2)?1:0;
    return true;
}

bool levelset::check_pos(int i,int ij) {
    if(s[ij]<4) {
        if(s[ij]==2) {num_n++;s[ij]=3;}
        return false;
    }
    if((ij>=m&&s[ij-m]<4)||(i>0&&s[ij-1]<4)||(i<m-1&&s[ij+1]<4)||(ij<mn-m&&s[ij+m]<4)) {
        if(s[ij]==5) {num_p++;s[ij]=4;}
    } else {
        if(s[ij]==4) {num_p--;s[ij]=5;}
    }
    return true;
}

bool levelset::check_neg(int i,int ij) {
    if(s[ij]>3) {
        if(s[ij]==5) {num_p++;s[ij]=4;}
        return false;
    }
    if((ij>=m&&s[ij-m]>3)||(i>0&&s[ij-1]>3)||(i<m-1&&s[ij+1]>3)||(ij<mn-m&&s[ij+m]>3)) {
        if(s[ij]==2) {num_n++;s[ij]=3;}
    } else {
        if(s[ij]==3) {num_n--;s[ij]=2;}
    }
    return true;
}

/** Computes the modulus of the gradient of the levelset function phi, using
 * centered differences. It is meant to be applied near the interface and
 * assumes that the values of phi exist at all of the gridpoints that are
 * referenced.
 * \param[in] i the horizontial grid point index.
 * \param[in] ij the grid point memory location.
 * \return the modules of the gradient of the levelset function. */
double levelset::centered_deriv(int i,int ij) {
    double phix,phiy;
    phix=i==0?(phi[ij+1]-phi[ij])*xsp:(i==m-1?(phi[ij]-phi[ij-1])*xsp:(phi[ij+1]-phi[ij-1])*xsp*0.5);
    phiy=ij<m?(phi[ij+m]-phi[ij])*ysp:(ij>=mn-m?(phi[ij]-phi[ij-m])*ysp:(phi[ij+m]-phi[ij-m])*ysp*0.5);
    return sqrt(phix*phix+phiy*phiy);
}

/** Counts the frequency of each number in the s array and prints
 * the results to the standard output. */
void levelset::count_s_array() {
    int i,ij,b[8];
    for(i=0;i<8;i++) b[i]=0;
    for(ij=0;ij<mn;ij++) {
        i=s[ij];
        if(i>=0&&i<8) b[i]++;
        else levelpp_fatal_error("s array element out of bounds",1);
    }
    for(i=0;i<8;i++) printf("%d ",b[i]);
    printf("%d %d",num_n,num_p);
    if (num_n!=b[3]||num_p!=b[4]) printf(" error");
    puts("");
}

/** Resets the back pointers. */
void levelset::reset_back_pointers() {
    for(int ij=0;ij<mn;ij++) bp[ij]=0;
    bpc=mn+1;
}

/** Resets the back pointers. */
void levelset::reset_staggered_back_pointers() {
    for(int ije=0;ije<mne;ije++) sg[ije]=0;
    sgc=1;
}

double levelset::laplacian(int i,int ij) {
    return (i!=0&&i!=m-1?xsp*xsp*(phi[ij+1]+phi[ij-1]-2*phi[ij]):0)
          +(ij>=m&&ij<mn-m?ysp*ysp*(phi[ij+m]+phi[ij-m]-2*phi[ij]):0);
}

void levelset::post_move_update(double dt) {
    int c,k,l,i,ij;
    u_start=u_mem;
    u_end=0;

    for(c=0;c<num_p;c++) {
        ij=t1[c];
        phi[ij]+=vp[ij]*dt;
        u[u_end++]=ij;
        if(u_end==u_start) levelpp_fatal_error("no dynamic memory extension implemented yet",1);
    }
    while(c<num_p+num_n) {
        ij=t1[c];
        phi[ij]+=vp[ij]*dt;
        u[--u_start]=ij;
        if(u_end==u_start) levelpp_fatal_error("no dynamic memory extension implemented yet",1);
        c++;
    }

    // Positive fast march
    heap<dir_positive> up(*this,mem_up,h_up);
    for(k=0;k<u_end;k++) {
        ij=u[k];i=ij%m;
        up.stack_test(i,ij);
    }
    up.build_heap();
    up.march(phi_up);

    // Negative fast march
    heap<dir_negative> down(*this,mem_down,h_down);
    for(k=u_start;k<u_mem;k++) {
        ij=u[k];i=ij%m;
        down.stack_test(i,ij);
    }
    down.build_heap();
    down.march(phi_down);

    // Sort the points in the band
    sort.sort(u,num_p);
    sort.sort(u+(u_mem-num_n),num_n);

    // Sort the boundary between the initial points and the fast marched points
    sort_boundary(u_mem-num_n,u_start,u_mem);
    sort_boundary(num_p,0,u_end);
    sort_boundary_wrap();

    check_sorted();

    // Just in case there was no interface, set u_start to zero
    if (u_start==u_mem) u_start=0;

    if(u_end==0||phi[u[0]]>0) u_mid=search_down(0,0.0);
    else u_mid=search_up(0,0.0);

    l=u_mid;
    while(l!=u_end) {
        ij=u[l];
        if(s[ij]<4) {i=ij%m;change_to_pos(i,ij);}
        l=ul_up(l);
    }
    l=u_mid;
    while(l!=u_start) {
        l=ul_down(l);
        ij=u[l];
        if(s[ij]>3) {i=ij%m;change_to_neg(i,ij);}
    }

    check_divided();

    //Clean up rest of phi[] grid
    clean();
}
